#include <mpfi.h>
#include <iostream>

void test_mpfi_interval_arithmetic() {
    mpfi_t r, r2, r6, r12, term1, term2, energy, one;
    
    mpfi_init2(r, 256); mpfi_init2(r2, 256); mpfi_init2(r6, 256); mpfi_init2(r12, 256);
    mpfi_init2(term1, 256); mpfi_init2(term2, 256); mpfi_init2(energy, 256); mpfi_init2(one, 256);
    
    // r = [1.1, 1.2]
    mpfi_interv_d(r, 1.1, 1.2);
    
    // r^2
    mpfi_sqr(r2, r);
    // r^6 = (r^2)^3
    mpfi_sqr(r6, r2);
    mpfi_mul(r6, r6, r2);
    // r^12 = (r^6)^2
    mpfi_sqr(r12, r6);
    
    mpfi_set_d(one, 1.0);
    // 1/r^12 and 1/r^6
    mpfi_inv(term1, r12);
    mpfi_inv(term2, r6);
    
    // energy = 4 * (1/r^12 - 1/r^6)
    mpfi_sub(energy, term1, term2);
    mpfi_mul_ui(energy, energy, 4);
    
    // Extract bounds safely using MPFI API
    mpfr_t lower_mpfr, upper_mpfr;
    mpfr_init2(lower_mpfr, 256);
    mpfr_init2(upper_mpfr, 256);
    
    mpfi_get_left(lower_mpfr, energy);
    mpfi_get_right(upper_mpfr, energy);
    
    double lower = mpfr_get_d(lower_mpfr, MPFR_RNDD);
    double upper = mpfr_get_d(upper_mpfr, MPFR_RNDU);
    
    std::cout << "MPFI: r in [1.1, 1.2] -> E in [" << lower << ", " << upper << "]" << std::endl;
    
    mpfr_clear(lower_mpfr);
    mpfr_clear(upper_mpfr);
    
    mpfi_clear(r); mpfi_clear(r2); mpfi_clear(r6); mpfi_clear(r12);
    mpfi_clear(term1); mpfi_clear(term2); mpfi_clear(energy); mpfi_clear(one);
}

int main() {
    test_mpfi_interval_arithmetic();
    return 0;
}