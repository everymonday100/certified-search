#include <flint/arb.h>
#include <flint/arf.h>
#include <iostream>

void test_arb_ball_arithmetic() {
    arb_t r, r2, r6, r12, term1, term2, energy;
    arf_t a, b, lo, hi;
    slong prec = 256;

    arb_init(r); arb_init(r2); arb_init(r6); arb_init(r12);
    arb_init(term1); arb_init(term2); arb_init(energy);
    arf_init(a); arf_init(b); arf_init(lo); arf_init(hi);

    // r = [1.1, 1.2] через точные arf-концы
    arf_set_d(a, 1.1);
    arf_set_d(b, 1.2);
    arb_set_interval_arf(r, a, b, prec);

    arb_sqr(r2, r, prec);
    arb_sqr(r6, r2, prec);
    arb_mul(r6, r6, r2, prec);
    arb_sqr(r12, r6, prec);

    arb_inv(term1, r12, prec);
    arb_inv(term2, r6, prec);

    arb_sub(energy, term1, term2, prec);
    arb_mul_ui(energy, energy, 4, prec);

    // Извлечение гарантированного интервала
    arb_get_interval_arf(lo, hi, energy, prec);
    double lower = arf_get_d(lo, ARF_RND_FLOOR);  // наружу
    double upper = arf_get_d(hi, ARF_RND_CEIL);   // наружу

    std::cout << "Arb: r in [1.1, 1.2] -> E in [" << lower << ", " << upper << "]" << std::endl;

    arb_clear(r); arb_clear(r2); arb_clear(r6); arb_clear(r12);
    arb_clear(term1); arb_clear(term2); arb_clear(energy);
    arf_clear(a); arf_clear(b); arf_clear(lo); arf_clear(hi);
}

int main() {
    test_arb_ball_arithmetic();
    return 0;
}