#include "lj_potential.h"
#include <cmath>
#include <stdexcept>

LJResult compute_lj_exact(const std::vector<double>& coords) {
    if (coords.size() % 3 != 0) {
        throw std::invalid_argument("Coordinates size must be a multiple of 3");
    }
    int N = coords.size() / 3;
    LJResult res;
    res.energy = 0.0;
    res.gradient.assign(3 * N, 0.0);
    res.hessian.assign(9 * N * N, 0.0); // 3N x 3N matrix

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double dx = coords[3*i]   - coords[3*j];
            double dy = coords[3*i+1] - coords[3*j+1];
            double dz = coords[3*i+2] - coords[3*j+2];
            
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < 1e-12) {
                throw std::runtime_error("Atoms overlapping or too close");
            }
            
            // Precompute inverse powers of r^2 to avoid pow()
            double r2_inv = 1.0 / r2;
            double r4_inv = r2_inv * r2_inv;
            double r6_inv = r4_inv * r2_inv;
            double r12_inv = r6_inv * r6_inv;
            
            // Energy: 4 * (r^-12 - r^-6)
            res.energy += 4.0 * (r12_inv - r6_inv);
            
            // Gradient factor A = V'(r)/r = 4 * (-12 r^-14 + 6 r^-8)
            double A = 4.0 * (-12.0 * r12_inv * r2_inv + 6.0 * r6_inv * r2_inv);
            
            double gx = A * dx;
            double gy = A * dy;
            double gz = A * dz;
            
            res.gradient[3*i]   += gx;
            res.gradient[3*i+1] += gy;
            res.gradient[3*i+2] += gz;
            res.gradient[3*j]   -= gx;
            res.gradient[3*j+1] -= gy;
            res.gradient[3*j+2] -= gz;
            
            // Hessian factor B = V''(r) - V'(r)/r = 4 * (168 r^-16 - 48 r^-10)
            double B = 4.0 * (168.0 * r12_inv * r4_inv - 48.0 * r6_inv * r4_inv);
            
            double rr[3] = {dx, dy, dz};
            int idx_i[3] = {3*i, 3*i+1, 3*i+2};
            int idx_j[3] = {3*j, 3*j+1, 3*j+2};
            
            // Fill 3x3 blocks for i-i, j-j, i-j, j-i
            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    double val = (a == b ? A : 0.0) + B * rr[a] * rr[b] * r2_inv;
                    res.hessian[(idx_i[a]) * (3*N) + idx_i[b]] += val;
                    res.hessian[(idx_j[a]) * (3*N) + idx_j[b]] += val;
                    res.hessian[(idx_i[a]) * (3*N) + idx_j[b]] -= val;
                    res.hessian[(idx_j[a]) * (3*N) + idx_i[b]] -= val;
                }
            }
        }
    }
    return res;
}