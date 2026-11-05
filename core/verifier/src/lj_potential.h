#pragma once
#include <string>
#include <vector>

struct LJResult {
    double energy;
    std::vector<double> gradient;  // size 3N
    std::vector<double> hessian;   // size 3N x 3N, row-major
};

// Computes exact LJ 12-6 energy, gradient, and Hessian.
// coords: flat array of size 3N (x0, y0, z0, x1, y1, z1, ...)
LJResult compute_lj_exact(const std::vector<double> &coords);