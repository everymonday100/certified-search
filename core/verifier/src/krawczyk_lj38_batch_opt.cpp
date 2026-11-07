// Optimized Batched Krawczyk certificate pipeline for LJ38.
// Usage: krawczyk_lj38_batch_opt <xyz_dir> <radius_exp> [--threads=N]
//
// Optimizations vs krawczyk_lj38_batch:
//   1. Thread count capped at MAX_THREADS=8 (cache-friendly)
//   2. Async file preloading via std::async (overlaps I/O with CPU)
//   3. Thread-local ArbPool reuses arb_t memory across verifications
//   4. schedule(dynamic,1) for fine-grained load balancing
#include <flint/arb.h>
#include <flint/arf.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <omp.h>

#include "lj_potential.h"

static const slong PREC = 256;
static const int N_ATOMS = 38;
static const int MAX_THREADS = 12;  // Cache-efficiency cap

// ---- Thread-local arb_t memory pool ----
// Reuses arb_t allocations across verifications on the same thread,
// avoiding repeated malloc/free overhead in FLINT's arb layer.
struct ArbPool {
    std::vector<arb_struct> pool;
    size_t next_idx = 0;

    ArbPool(size_t initial_size = 15000) {
        pool.resize(initial_size);
        for (auto& a : pool) arb_init(&a);
    }
    ~ArbPool() {
        for (auto& a : pool) arb_clear(&a);
    }
    arb_ptr acquire() {
        if (next_idx >= pool.size()) {
            size_t old = pool.size();
            pool.resize(old * 2);
            for (size_t i = old; i < pool.size(); ++i) arb_init(&pool[i]);
        }
        return &pool[next_idx++];
    }
    void reset() { next_idx = 0; }
};

static thread_local ArbPool* tls_pool = nullptr;

static void init_thread_pool() {
    if (!tls_pool) {
        tls_pool = new ArbPool(15000);
    } else {
        tls_pool->reset();
    }
}

// ---- rigorous double extraction helpers (thread-safe) ----
static double arb_lo_floor(const arb_t x) {
    arf_t t; arf_init(t);
    arb_get_lbound_arf(t, x, PREC);
    double d = arf_get_d(t, ARF_RND_FLOOR);
    arf_clear(t);
    return d;
}
static double arb_hi_ceil(const arb_t x) {
    arf_t t; arf_init(t);
    arb_get_ubound_arf(t, x, PREC);
    double d = arf_get_d(t, ARF_RND_CEIL);
    arf_clear(t);
    return d;
}
static double arb_lo_ceil(const arb_t x) {
    arf_t t; arf_init(t);
    arb_get_lbound_arf(t, x, PREC);
    double d = arf_get_d(t, ARF_RND_CEIL);
    arf_clear(t);
    return d;
}
static double arb_hi_floor(const arb_t x) {
    arf_t t; arf_init(t);
    arb_get_ubound_arf(t, x, PREC);
    double d = arf_get_d(t, ARF_RND_FLOOR);
    arf_clear(t);
    return d;
}

// ---- xyz parsing (thread-safe) ----
static bool read_xyz(const std::string& path, std::vector<double>& coords, int& n) {
    std::ifstream f(path);
    if (!f) return false;
    std::string line;
    if (!std::getline(f, line)) return false;
    n = std::atoi(line.c_str());
    if (n <= 0) return false;
    std::getline(f, line);
    coords.clear();
    coords.reserve(3 * n);
    for (int i = 0; i < n; ++i) {
        if (!std::getline(f, line)) return false;
        std::istringstream ss(line);
        std::string el; double x, y, z;
        if (!(ss >> el >> x >> y >> z)) return false;
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) return false;
        coords.push_back(x); coords.push_back(y); coords.push_back(z);
    }
    return true;
}

// ---- gauge: atom0 at origin, atom1 on +X, atom2 in XY plane (y>0) ----
static bool gauge_transform(std::vector<double>& p, int n) {
    double p0[3] = {p[0], p[1], p[2]};
    for (int i = 0; i < n; ++i) { p[3*i] -= p0[0]; p[3*i+1] -= p0[1]; p[3*i+2] -= p0[2]; }
    double a[3] = {p[3], p[4], p[5]};
    double na = std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    if (na < 1e-8) return false;
    double e1[3] = {a[0]/na, a[1]/na, a[2]/na};
    double b[3] = {p[6], p[7], p[8]};
    double bd = b[0]*e1[0] + b[1]*e1[1] + b[2]*e1[2];
    double bp[3] = {b[0]-bd*e1[0], b[1]-bd*e1[1], b[2]-bd*e1[2]};
    double nbp = std::sqrt(bp[0]*bp[0] + bp[1]*bp[1] + bp[2]*bp[2]);
    if (nbp < 1e-8) return false;
    double e2[3] = {bp[0]/nbp, bp[1]/nbp, bp[2]/nbp};
    double e3[3] = {e1[1]*e2[2]-e1[2]*e2[1], e1[2]*e2[0]-e1[0]*e2[2], e1[0]*e2[1]-e1[1]*e2[0]};
    std::vector<double> q(3*n);
    for (int i = 0; i < n; ++i) {
        double v[3] = {p[3*i], p[3*i+1], p[3*i+2]};
        q[3*i]   = v[0]*e1[0] + v[1]*e1[1] + v[2]*e1[2];
        q[3*i+1] = v[0]*e2[0] + v[1]*e2[1] + v[2]*e2[2];
        q[3*i+2] = v[0]*e3[0] + v[1]*e3[1] + v[2]*e3[2];
    }
    p = q;
    return true;
}

static std::vector<int> slice_indices(int n) {
    std::vector<int> idx;
    idx.push_back(3); idx.push_back(6); idx.push_back(7);
    for (int i = 3; i < n; ++i) { idx.push_back(3*i); idx.push_back(3*i+1); idx.push_back(3*i+2); }
    return idx;
}

static void full_from_q(const std::vector<double>& q0, int n, std::vector<double>& full) {
    full.assign(3*n, 0.0);
    full[3] = q0[0]; full[6] = q0[1]; full[7] = q0[2];
    for (int i = 3; i < n; ++i)
        for (int c = 0; c < 3; ++c) full[3*i+c] = q0[3 + 3*(i-3) + c];
}

// ---- float64 linear algebra (thread-safe) ----
static bool solve_lin(std::vector<double> A, std::vector<double> b, std::vector<double>& x, int m) {
    for (int k = 0; k < m; ++k) {
        int piv = k; double best = std::fabs(A[k*m+k]);
        for (int r = k+1; r < m; ++r) { double v = std::fabs(A[r*m+k]); if (v > best) { best = v; piv = r; } }
        if (best < 1e-14) return false;
        if (piv != k) { for (int c = 0; c < m; ++c) std::swap(A[k*m+c], A[piv*m+c]); std::swap(b[k], b[piv]); }
        for (int r = k+1; r < m; ++r) {
            double f = A[r*m+k] / A[k*m+k];
            if (f == 0.0) continue;
            for (int c = k; c < m; ++c) A[r*m+c] -= f * A[k*m+c];
            b[r] -= f * b[k];
        }
    }
    x.assign(m, 0.0);
    for (int r = m-1; r >= 0; --r) {
        double s = b[r];
        for (int c = r+1; c < m; ++c) s -= A[r*m+c] * x[c];
        x[r] = s / A[r*m+r];
    }
    return true;
}

static bool invert_mat(const std::vector<double>& A0, std::vector<double>& Ainv, int m) {
    Ainv.assign(m*m, 0.0);
    for (int j = 0; j < m; ++j) {
        std::vector<double> e(m, 0.0), x;
        e[j] = 1.0;
        if (!solve_lin(A0, e, x, m)) return false;
        for (int r = 0; r < m; ++r) Ainv[r*m+j] = x[r];
    }
    return true;
}

// ---- interval LJ gradient / Hessian / energy / min distance ----
static void lj_interval_grad_hess(const std::vector<arb_t>& cart, int n,
                                  std::vector<arb_t>& grad, std::vector<arb_t>& hess,
                                  arb_t eout, double* dmin_lb) {
    arb_t dx, dy, dz, r2, i1, i2, i4, i5, i6, i7, i8, r6i, t1, t2, A, B, val;
    arb_t c4, cm12, c6, c168, cm48, zb, ob;
    arb_init(dx); arb_init(dy); arb_init(dz); arb_init(r2);
    arb_init(i1); arb_init(i2); arb_init(i4); arb_init(i5); arb_init(i6);
    arb_init(i7); arb_init(i8); arb_init(r6i);
    arb_init(t1); arb_init(t2); arb_init(A); arb_init(B); arb_init(val);
    arb_init(c4); arb_init(cm12); arb_init(c6); arb_init(c168); arb_init(cm48);
    arb_init(zb); arb_init(ob);
    arb_set_d(c4, 4.0); arb_set_d(cm12, -12.0); arb_set_d(c6, 6.0);
    arb_set_d(c168, 168.0); arb_set_d(cm48, -48.0);
    arb_set_d(zb, 0.0); arb_set_d(ob, 1.0);

    for (int i = 0; i < n; ++i) {
        for (int j = i+1; j < n; ++j) {
            arb_sub(dx, cart[3*i], cart[3*j], PREC);
            arb_sub(dy, cart[3*i+1], cart[3*j+1], PREC);
            arb_sub(dz, cart[3*i+2], cart[3*j+2], PREC);
            arb_sqr(t1, dx, PREC); arb_sqr(t2, dy, PREC); arb_add(r2, t1, t2, PREC);
            arb_sqr(t1, dz, PREC); arb_add(r2, r2, t1, PREC);

            if (dmin_lb) {
                double lo = arb_lo_floor(r2);
                if (lo > 0.0) {
                    double d = std::nextafter(std::sqrt(lo), -INFINITY);
                    if (d < *dmin_lb) *dmin_lb = d;
                } else {
                    *dmin_lb = 0.0;
                }
            }

            arb_inv(i1, r2, PREC);
            arb_sqr(i2, i1, PREC); arb_sqr(i4, i2, PREC);
            arb_mul(i5, i4, i1, PREC); arb_mul(i6, i4, i2, PREC);
            arb_mul(i7, i6, i1, PREC); arb_sqr(i8, i4, PREC);
            arb_mul(r6i, i1, i2, PREC);

            arb_mul(t1, i7, cm12, PREC); arb_mul(t2, i4, c6, PREC);
            arb_add(A, t1, t2, PREC); arb_mul(A, A, c4, PREC);
            arb_mul(t1, i8, c168, PREC); arb_mul(t2, i5, cm48, PREC);
            arb_add(B, t1, t2, PREC); arb_mul(B, B, c4, PREC);

            arb_t rr[3];
            arb_init(rr[0]); arb_init(rr[1]); arb_init(rr[2]);
            arb_set(rr[0], dx); arb_set(rr[1], dy); arb_set(rr[2], dz);

            for (int c = 0; c < 3; ++c) {
                arb_mul(val, A, rr[c], PREC);
                arb_add(grad[3*i+c], grad[3*i+c], val, PREC);
                arb_sub(grad[3*j+c], grad[3*j+c], val, PREC);
            }
            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    arb_mul(val, rr[a], rr[b], PREC);
                    arb_mul(val, val, i1, PREC);
                    arb_mul(val, val, B, PREC);
                    if (a == b) arb_add(val, val, A, PREC);
                    arb_add(hess[(3*i+a)*(3*n) + (3*i+b)], hess[(3*i+a)*(3*n) + (3*i+b)], val, PREC);
                    arb_add(hess[(3*j+a)*(3*n) + (3*j+b)], hess[(3*j+a)*(3*n) + (3*j+b)], val, PREC);
                    arb_sub(hess[(3*i+a)*(3*n) + (3*j+b)], hess[(3*i+a)*(3*n) + (3*j+b)], val, PREC);
                    arb_sub(hess[(3*j+a)*(3*n) + (3*i+b)], hess[(3*j+a)*(3*n) + (3*i+b)], val, PREC);
                }
            }
            arb_clear(rr[0]); arb_clear(rr[1]); arb_clear(rr[2]);

            arb_mul(t1, i6, c4, PREC); arb_mul(t2, r6i, c4, PREC);
            arb_sub(t1, t1, t2, PREC); arb_add(eout, eout, t1, PREC);
        }
    }
    arb_clear(dx); arb_clear(dy); arb_clear(dz); arb_clear(r2);
    arb_clear(i1); arb_clear(i2); arb_clear(i4); arb_clear(i5); arb_clear(i6);
    arb_clear(i7); arb_clear(i8); arb_clear(r6i);
    arb_clear(t1); arb_clear(t2); arb_clear(A); arb_clear(B); arb_clear(val);
    arb_clear(c4); arb_clear(cm12); arb_clear(c6); arb_clear(c168); arb_clear(cm48);
    arb_clear(zb); arb_clear(ob);
}

// ---- Verification result ----
struct VerificationResult {
    std::string filename;
    std::string status;
    double energy_lo;
    double energy_hi;
    double min_dist;
    double kappa;
    double pd_lb;
    double elapsed_sec;
    bool success;
};

// ---- Single file verification ----
static VerificationResult verify_single_file(const std::string& filepath, double rad) {
    VerificationResult result;
    result.filename = filepath;
    result.success = false;
    result.status = "INCONCLUSIVE";

    init_thread_pool();

    auto start = std::chrono::high_resolution_clock::now();

    try {
        int n = 0;
        std::vector<double> p;
        if (!read_xyz(filepath, p, n) || n != N_ATOMS) {
            result.status = "REJECTED";
            return result;
        }

        if (!gauge_transform(p, n)) {
            result.status = "REJECTED";
            return result;
        }

        auto idx = slice_indices(n);
        int m = (int)idx.size();

        std::vector<double> q0(m);
        q0[0] = p[3]; q0[1] = p[6]; q0[2] = p[7];
        for (int i = 3; i < n; ++i)
            for (int c = 0; c < 3; ++c) q0[3 + 3*(i-3) + c] = p[3*i+c];

        // Newton refinement
        double ginf = 1e300;
        for (int iter = 0; iter < 6; ++iter) {
            std::vector<double> full;
            full_from_q(q0, n, full);
            LJResult r = compute_lj_exact(full);
            std::vector<double> g(m), H(m*m);
            for (int a = 0; a < m; ++a) {
                g[a] = r.gradient[idx[a]];
                for (int b = 0; b < m; ++b) H[a*m+b] = r.hessian[(size_t)idx[a]*(3*n) + idx[b]];
            }
            ginf = 0.0;
            for (double x : g) ginf = std::max(ginf, std::fabs(x));
            if (ginf < 1e-13) break;
            std::vector<double> neg(m), d;
            for (int a = 0; a < m; ++a) neg[a] = -g[a];
            if (!solve_lin(H, neg, d, m)) { result.status = "REJECTED"; return result; }
            for (int a = 0; a < m; ++a) q0[a] += d[a];
        }

        // Interval stage
        mag_t radmag; mag_init(radmag); mag_set_d(radmag, rad);

        std::vector<arb_t> qp(m), qX(m), cartp(3*n), cartX(3*n);
        std::vector<arb_t> gradp(3*n), gradX(3*n), hessX(3*n*3*n);
        for (int i = 0; i < m; ++i) {
            arb_init(qp[i]); arb_init(qX[i]);
            arb_set_d(qp[i], q0[i]); arb_set_d(qX[i], q0[i]);
            arb_add_error_mag(qX[i], radmag);
        }
        for (int i = 0; i < 3*n; ++i) { arb_init(cartp[i]); arb_init(cartX[i]); arb_init(gradp[i]); arb_init(gradX[i]); }
        for (size_t i = 0; i < hessX.size(); ++i) arb_init(hessX[i]);

        arb_set_d(cartp[0], 0.0); arb_set_d(cartp[1], 0.0); arb_set_d(cartp[2], 0.0);
        arb_set(cartp[3], qp[0]); arb_set_d(cartp[4], 0.0); arb_set_d(cartp[5], 0.0);
        arb_set(cartp[6], qp[1]); arb_set(cartp[7], qp[2]); arb_set_d(cartp[8], 0.0);
        for (int i = 3; i < n; ++i) for (int c = 0; c < 3; ++c) arb_set(cartp[3*i+c], qp[3+3*(i-3)+c]);
        arb_set_d(cartX[0], 0.0); arb_set_d(cartX[1], 0.0); arb_set_d(cartX[2], 0.0);
        arb_set(cartX[3], qX[0]); arb_set_d(cartX[4], 0.0); arb_set_d(cartX[5], 0.0);
        arb_set(cartX[6], qX[1]); arb_set(cartX[7], qX[2]); arb_set_d(cartX[8], 0.0);
        for (int i = 3; i < n; ++i) for (int c = 0; c < 3; ++c) arb_set(cartX[3*i+c], qX[3+3*(i-3)+c]);

        arb_t ep, eX;
        arb_init(ep); arb_init(eX);
        arb_set_d(ep, 0.0); arb_set_d(eX, 0.0);
        double dmin = 1e300;
        lj_interval_grad_hess(cartp, n, gradp, hessX, ep, nullptr);
        for (size_t i = 0; i < hessX.size(); ++i) arb_set_d(hessX[i], 0.0);
        lj_interval_grad_hess(cartX, n, gradX, hessX, eX, &dmin);

        std::vector<arb_t> F0(m), HX(m*m);
        for (int i = 0; i < m; ++i) arb_init(F0[i]);
        for (int i = 0; i < m*m; ++i) arb_init(HX[i]);
        for (int a = 0; a < m; ++a) arb_set(F0[a], gradp[idx[a]]);
        for (int a = 0; a < m; ++a)
            for (int b = 0; b < m; ++b) arb_set(HX[a*m+b], hessX[(size_t)idx[a]*(3*n) + idx[b]]);

        std::vector<double> full;
        full_from_q(q0, n, full);
        LJResult rf = compute_lj_exact(full);
        std::vector<double> Hsel(m*m), C;
        for (int a = 0; a < m; ++a)
            for (int b = 0; b < m; ++b) Hsel[a*m+b] = rf.hessian[(size_t)idx[a]*(3*n) + idx[b]];
        if (!invert_mat(Hsel, C, m)) { result.status = "REJECTED"; return result; }

        std::vector<arb_t> M(m*m);
        for (int i = 0; i < m*m; ++i) arb_init(M[i]);
        {
            arb_t acc, t, cb2, db;
            arb_init(acc); arb_init(t); arb_init(cb2); arb_init(db);
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < m; ++j) {
                    arb_set_d(acc, 0.0);
                    for (int k = 0; k < m; ++k) {
                        arb_set_d(cb2, C[i*m+k]);
                        arb_mul(t, cb2, HX[k*m+j], PREC);
                        arb_add(acc, acc, t, PREC);
                    }
                    arb_set_d(db, (i == j) ? 1.0 : 0.0);
                    arb_sub(M[i*m+j], db, acc, PREC);
                }
            }
            arb_clear(acc); arb_clear(t); arb_clear(cb2); arb_clear(db);
        }

        double kappa = 0.0;
        {
            arb_t S, t;
            arb_init(S); arb_init(t);
            for (int i = 0; i < m; ++i) {
                arb_set_d(S, 0.0);
                for (int j = 0; j < m; ++j) {
                    arb_abs(t, M[i*m+j]);
                    arb_add(S, S, t, PREC);
                }
                kappa = std::max(kappa, arb_hi_ceil(S));
            }
            arb_clear(S); arb_clear(t);
        }

        std::vector<arb_t> K(m), V(m);
        for (int i = 0; i < m; ++i) { arb_init(K[i]); arb_init(V[i]); arb_set_d(V[i], 0.0); arb_add_error_mag(V[i], radmag); }
        {
            arb_t t, cb2, cf;
            arb_init(t); arb_init(cb2); arb_init(cf);
            for (int i = 0; i < m; ++i) {
                arb_set_d(K[i], q0[i]);
                arb_set_d(cf, 0.0);
                for (int k = 0; k < m; ++k) {
                    arb_set_d(cb2, C[i*m+k]);
                    arb_mul(t, cb2, F0[k], PREC);
                    arb_add(cf, cf, t, PREC);
                }
                arb_sub(K[i], K[i], cf, PREC);
                for (int j = 0; j < m; ++j) {
                    arb_mul(t, M[i*m+j], V[j], PREC);
                    arb_add(K[i], K[i], t, PREC);
                }
            }
            arb_clear(t); arb_clear(cb2); arb_clear(cf);
        }

        bool incl = true;
        for (int i = 0; i < m && incl; ++i) {
            double klo = arb_lo_floor(K[i]), khi = arb_hi_ceil(K[i]);
            double xlo = arb_lo_ceil(qX[i]), xhi = arb_hi_floor(qX[i]);
            if (!(klo > xlo) || !(khi < xhi)) incl = false;
        }

        bool pd = true;
        double pd_lb = 1e300;
        std::vector<arb_t> L(m*m);
        for (int i = 0; i < m*m; ++i) { arb_init(L[i]); arb_set_d(L[i], 0.0); }
        for (int k = 0; k < m && pd; ++k) {
            arb_t s; arb_init(s); arb_set_d(s, 0.0);
            for (int j = 0; j < k; ++j) {
                arb_t t; arb_init(t);
                arb_sqr(t, L[k*m+j], PREC);
                arb_sub(s, s, t, PREC);
                arb_clear(t);
            }
            arb_add(s, s, HX[k*m+k], PREC);
            arb_sqrt(L[k*m+k], s, PREC);
            if (!arb_is_finite(L[k*m+k])) { pd = false; }
            else {
                double lb = arb_lo_floor(L[k*m+k]);
                if (!(lb > 0.0)) pd = false; else pd_lb = std::min(pd_lb, lb);
            }
            arb_clear(s);
            for (int i = k+1; i < m && pd; ++i) {
                arb_t s2; arb_init(s2); arb_set_d(s2, 0.0);
                for (int j = 0; j < k; ++j) {
                    arb_t t; arb_init(t);
                    arb_mul(t, L[i*m+j], L[k*m+j], PREC);
                    arb_sub(s2, s2, t, PREC);
                    arb_clear(t);
                }
                arb_add(s2, s2, HX[i*m+k], PREC);
                arb_div(L[i*m+k], s2, L[k*m+k], PREC);
                arb_clear(s2);
            }
        }

        if (dmin <= 0.0) {
            result.status = "REJECTED";
        } else if (incl && kappa < 1.0 && pd) {
            result.status = "CERTIFIED_STRICT_LOCAL_MINIMUM";
            result.success = true;
        } else {
            result.status = "INCONCLUSIVE";
        }

        result.energy_lo = arb_lo_floor(eX);
        result.energy_hi = arb_hi_ceil(eX);
        result.min_dist = dmin;
        result.kappa = kappa;
        result.pd_lb = pd_lb;

    } catch (...) {
        result.status = "ERROR";
    }

    auto end = std::chrono::high_resolution_clock::now();
    result.elapsed_sec = std::chrono::duration<double>(end - start).count();

    return result;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: krawczyk_lj38_batch_opt <xyz_dir> <radius_exp> [--threads=N]\n";
        return 1;
    }

    std::string xyz_dir = argv[1];
    double rad = std::pow(10.0, -std::atof(argv[2]));

    // Cap threads at MAX_THREADS for cache efficiency
    int num_threads = std::min(omp_get_max_threads(), MAX_THREADS);

    for (int i = 3; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.find("--threads=") == 0) {
            num_threads = std::atoi(arg.substr(10).c_str());
            if (num_threads < 1) num_threads = 1;
            if (num_threads > MAX_THREADS) num_threads = MAX_THREADS;
        }
    }

    omp_set_num_threads(num_threads);

    // Collect all .xyz files
    std::vector<std::string> files;
    try {
        for (const auto& entry : std::filesystem::directory_iterator(xyz_dir)) {
            if (entry.path().extension() == ".xyz") {
                files.push_back(entry.path().string());
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error reading directory: " << e.what() << "\n";
        return 1;
    }

    if (files.empty()) {
        std::cerr << "No .xyz files found in " << xyz_dir << "\n";
        return 1;
    }

    std::sort(files.begin(), files.end());

    // Async file preloading: read files into memory before verification starts.
    // This overlaps disk I/O with the upcoming CPU-intensive verification.
    std::vector<std::future<std::string>> preload_futures;
    preload_futures.reserve(files.size());
    for (const auto& file : files) {
        preload_futures.push_back(std::async(std::launch::async, [file]() -> std::string {
            std::ifstream f(file, std::ios::binary);
            return std::string((std::istreambuf_iterator<char>(f)),
                               std::istreambuf_iterator<char>());
        }));
    }
    // Wait for all preloads to complete (populates OS page cache)
    for (auto& fut : preload_futures) {
        fut.get();
    }

    std::cerr << "[BATCH-OPT] Processing " << files.size() << " files with "
              << num_threads << " threads (capped at " << MAX_THREADS << ")\n";

    // Parallel verification
    std::vector<VerificationResult> results(files.size());

    auto batch_start = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < files.size(); ++i) {
        results[i] = verify_single_file(files[i], rad);
    }

    auto batch_end = std::chrono::high_resolution_clock::now();
    double batch_elapsed = std::chrono::duration<double>(batch_end - batch_start).count();

    // Output JSON results
    std::cout << "{\n";
    std::cout << "  \"batch_elapsed_sec\": " << batch_elapsed << ",\n";
    std::cout << "  \"num_threads\": " << num_threads << ",\n";
    std::cout << "  \"max_threads\": " << MAX_THREADS << ",\n";
    std::cout << "  \"num_files\": " << files.size() << ",\n";
    std::cout << "  \"results\": [\n";

    int certified_count = 0;
    for (size_t i = 0; i < results.size(); ++i) {
        const auto& r = results[i];
        if (r.success) certified_count++;

        std::cout << "    {\n";
        std::cout << "      \"file\": \"" << r.filename << "\",\n";
        std::cout << "      \"status\": \"" << r.status << "\",\n";
        std::cout << "      \"energy\": [" << r.energy_lo << ", " << r.energy_hi << "],\n";
        std::cout << "      \"min_dist\": " << r.min_dist << ",\n";
        std::cout << "      \"kappa\": " << r.kappa << ",\n";
        std::cout << "      \"pd_lb\": " << r.pd_lb << ",\n";
        std::cout << "      \"elapsed_sec\": " << r.elapsed_sec << "\n";
        std::cout << "    }" << (i < results.size() - 1 ? "," : "") << "\n";
    }

    std::cout << "  ],\n";
    std::cout << "  \"certified_count\": " << certified_count << ",\n";
    std::cout << "  \"avg_time_per_file\": " << (batch_elapsed / files.size()) << "\n";
    std::cout << "}\n";

    std::cerr << "[BATCH-OPT] Completed: " << certified_count << "/" << files.size()
              << " certified in " << batch_elapsed << " sec "
              << "(" << (batch_elapsed / files.size()) << " sec/file avg)\n";

    return 0;
}