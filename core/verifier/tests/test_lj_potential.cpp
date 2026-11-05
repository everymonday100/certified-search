#include "../src/lj_potential.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

void test_basic_math() {
    // Test 1: r = 1.0 -> Energy = 0
    std::vector<double> coords1 = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
    LJResult res1 = compute_lj_exact(coords1);
    assert(std::abs(res1.energy - 0.0) < 1e-9);
    
    // Test 2: r = 2^(1/6) ~ 1.122462 -> Minimum Energy = -1.0, Gradient = 0
    double r_min = std::pow(2.0, 1.0/6.0);
    std::vector<double> coords2 = {0.0, 0.0, 0.0, r_min, 0.0, 0.0};
    LJResult res2 = compute_lj_exact(coords2);
    assert(std::abs(res2.energy - (-1.0)) < 1e-9);
    assert(std::abs(res2.gradient[0]) < 1e-9);
    
    std::cout << "✅ Basic LJ math tests passed!" << std::endl;
}

int main() {
    test_basic_math();
    return 0;
}