#include <flint/arb.h>
#include <flint/mag.h>
#include <flint/arf.h>
#include <iostream>

void test_arb_ball_arithmetic() {
    arb_t r, r2, r6, r12, term1, term2, energy;
    mag_t rad;
    arf_t lo, hi;
    slong prec = 256;

    arb_init(r); arb_init(r2); arb_init(r6); arb_init(r12);
    arb_init(term1); arb_init(term2); arb_init(energy);
    mag_init(rad);
    arf_init(lo); arf_init(hi);

    // Создаём ball, строго покрывающий [1.1, 1.2]
    // 1. Центр: 1.15 (arb_set_d гарантирует покрытие ошибки округления double)
    arb_set_d(r, 1.15);
    
    // 2. Радиус: 0.05 (mag_set_d гарантирует outward rounding)
    mag_set_d(rad, 0.05);
    
    // 3. Добавляем радиус к ball: r = 1.15 ± (machine_eps + 0.05)
    arb_add_error_mag(r, rad);

    // Вычисления с явной точностью
    arb_sqr(r2, r, prec);
    arb_sqr(r6, r2, prec);
    arb_mul(r6, r6, r2, prec);
    arb_sqr(r12, r6, prec);

    arb_inv(term1, r12, prec);
    arb_inv(term2, r6, prec);

    arb_sub(energy, term1, term2, prec);
    arb_mul_ui(energy, energy, 4, prec);

    // Извлечение гарантированных границ (outward rounding по умолчанию)
    arb_get_lbound_arf(lo, energy, prec);  // строгая нижняя граница
    arb_get_ubound_arf(hi, energy, prec);  // строгая верхняя граница
    
    double lower = arf_get_d(lo, ARF_RND_FLOOR);  // округляем вниз (наружу)
    double upper = arf_get_d(hi, ARF_RND_CEIL);   // округляем вверх (наружу)

    std::cout << "Arb: r in [1.1, 1.2] -> E in [" << lower << ", " << upper << "]" << std::endl;

    arb_clear(r); arb_clear(r2); arb_clear(r6); arb_clear(r12);
    arb_clear(term1); arb_clear(term2); arb_clear(energy);
    mag_clear(rad);
    arf_clear(lo); arf_clear(hi);
}

int main() {
    test_arb_ball_arithmetic();
    return 0;
}