#include <flint/arb.h>
#include <flint/mag.h>
#include <iostream>

void test_arb_wide_interval() {
    arb_t r, r2, r6, r12, term1, term2, energy;
    mag_t rad;
    arf_t lo, hi;
    slong prec = 256;

    arb_init(r); arb_init(r2); arb_init(r6); arb_init(r12);
    arb_init(term1); arb_init(term2); arb_init(energy);
    mag_init(rad);
    arf_init(lo); arf_init(hi);

    // Широкий интервал [1.1, 1.2] через ball arithmetic
    // Центр: 1.15, Радиус: 0.05
    arb_set_d(r, 1.15);           // Создаёт ball с покрытием ошибки округления
    mag_set_d(rad, 0.05);         // Радиус 0.05
    arb_add_error_mag(r, rad);    // r = 1.15 ± (eps + 0.05)

    arb_sqr(r2, r, prec);
    arb_sqr(r6, r2, prec);
    arb_mul(r6, r6, r2, prec);
    arb_sqr(r12, r6, prec);

    arb_inv(term1, r12, prec);
    arb_inv(term2, r6, prec);

    arb_sub(energy, term1, term2, prec);
    arb_mul_ui(energy, energy, 4, prec);

    // Извлечение границ
    arb_get_lbound_arf(lo, energy, prec);
    arb_get_ubound_arf(hi, energy, prec);
    
    double lower = arf_get_d(lo, ARF_RND_FLOOR);
    double upper = arf_get_d(hi, ARF_RND_CEIL);

    std::cout << "Arb (wide):   r in [1.1, 1.2]     -> E in [" << lower << ", " << upper << "]" << std::endl;

    arb_clear(r); arb_clear(r2); arb_clear(r6); arb_clear(r12);
    arb_clear(term1); arb_clear(term2); arb_clear(energy);
    mag_clear(rad);
    arf_clear(lo); arf_clear(hi);
}

void test_arb_narrow_interval() {
    arb_t r, r2, r6, r12, term1, term2, energy;
    mag_t rad;
    arf_t lo, hi;
    slong prec = 256;

    arb_init(r); arb_init(r2); arb_init(r6); arb_init(r12);
    arb_init(term1); arb_init(term2); arb_init(energy);
    mag_init(rad);
    arf_init(lo); arf_init(hi);

    // Узкий интервал [1.122, 1.123] через ball arithmetic
    // Центр: 1.1225, Радиус: 0.0005
    arb_set_d(r, 1.1225);
    mag_set_d(rad, 0.0005);
    arb_add_error_mag(r, rad);

    arb_sqr(r2, r, prec);
    arb_sqr(r6, r2, prec);
    arb_mul(r6, r6, r2, prec);
    arb_sqr(r12, r6, prec);

    arb_inv(term1, r12, prec);
    arb_inv(term2, r6, prec);

    arb_sub(energy, term1, term2, prec);
    arb_mul_ui(energy, energy, 4, prec);

    arb_get_lbound_arf(lo, energy, prec);
    arb_get_ubound_arf(hi, energy, prec);
    
    double lower = arf_get_d(lo, ARF_RND_FLOOR);
    double upper = arf_get_d(hi, ARF_RND_CEIL);

    std::cout << "Arb (narrow): r in [1.122, 1.123] -> E in [" << lower << ", " << upper << "]" << std::endl;

    arb_clear(r); arb_clear(r2); arb_clear(r6); arb_clear(r12);
    arb_clear(term1); arb_clear(term2); arb_clear(energy);
    mag_clear(rad);
    arf_clear(lo); arf_clear(hi);
}

int main() {
    test_arb_wide_interval();
    test_arb_narrow_interval();
    return 0;
}