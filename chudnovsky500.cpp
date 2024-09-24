#include <iostream>
#include <mpfr.h>

void compute_pi_chudnovsky(mpfr_t pi, long digits) {
    // Set precision in bits (log2(10) ≈ 3.32 bits per decimal digit)
    long precision_bits = digits * 3.32193;  // Approximate conversion from digits to bits

    // Initialize necessary variables with the set precision
    mpfr_set_default_prec(precision_bits);

    mpfr_t sum, term, factorial_k, factorial_3k, factorial_k_3, power_640320, temp, constant;
    mpfr_inits(sum, term, factorial_k, factorial_3k, factorial_k_3, power_640320, temp, constant, NULL);

    // Set up constants used in the Chudnovsky algorithm
    mpfr_set_ui(sum, 0, MPFR_RNDN);                 // sum = 0
    mpfr_set_ui(factorial_k, 1, MPFR_RNDN);          // factorial_k = 1
    mpfr_set_ui(factorial_3k, 1, MPFR_RNDN);         // factorial_3k = 1
    mpfr_set_ui(factorial_k_3, 1, MPFR_RNDN);        // factorial_k_3 = 1

    mpfr_ui_pow_ui(power_640320, 640320, 3, MPFR_RNDN);  // 640320^3
    mpfr_sqrt_ui(temp, 10005, MPFR_RNDN);                // sqrt(10005)
    mpfr_mul_ui(constant, temp, 426880, MPFR_RNDN);      // constant = 426880 * sqrt(10005)

    // Loop for Chudnovsky series
    for (long k = 0; k < (digits / 14) + 1; ++k) {
        // Calculate the numerator: (-1)^k * (6k)! * (13591409 + 545140134 * k)
        mpfr_fac_ui(factorial_k, 6 * k, MPFR_RNDN);    // factorial_k = (6k)!
        mpfr_mul_si(temp, factorial_k, 13591409 + 545140134 * k, MPFR_RNDN);  // temp = (13591409 + 545140134 * k)

        // If k is odd, negate the term
        if (k % 2 == 1) mpfr_neg(temp, temp, MPFR_RNDN);  // (-1)^k factor

        // Calculate the denominator: (3k)! * (k!)^3 * (640320^(3k + 3/2))
        mpfr_fac_ui(factorial_3k, 3 * k, MPFR_RNDN);   // factorial_3k = (3k)!
        mpfr_fac_ui(factorial_k, k, MPFR_RNDN);        // factorial_k = k!
        mpfr_pow_ui(factorial_k_3, factorial_k, 3, MPFR_RNDN);  // factorial_k_3 = (k!)^3

        // Calculate power_640320^(3k)
        mpfr_pow_ui(power_640320, power_640320, 3 * k, MPFR_RNDN);

        // Multiply all the denominators
        mpfr_mul(term, factorial_3k, factorial_k_3, MPFR_RNDN);  // term = (3k)! * (k!)^3
        mpfr_mul(term, term, power_640320, MPFR_RNDN);           // term *= 640320^(3k)

        // Compute final term = numerator / denominator
        mpfr_div(term, temp, term, MPFR_RNDN);  // term = numerator / denominator

        // Add term to sum
        mpfr_add(sum, sum, term, MPFR_RNDN);
    }

    // Multiply sum by constant to get 1/pi, then invert to get pi
    mpfr_mul(sum, sum, constant, MPFR_RNDN);  // sum = constant * sum
    mpfr_ui_div(pi, 1, sum, MPFR_RNDN);       // pi = 1 / sum

    // Clear variables
    mpfr_clears(sum, term, factorial_k, factorial_3k, factorial_k_3, power_640320, temp, constant, NULL);
}

int main() {
    // Set the number of digits of π we want
    long digits = 500000;

    // Create an mpfr_t variable to store the result
    mpfr_t pi;
    mpfr_init2(pi, digits * 3.32193);  // Set precision in bits

    // Compute π
    compute_pi_chudnovsky(pi, digits);

    // Output π with high precision
    mpfr_printf("%.500000Rf\n", pi);  // Print pi with specified digits

    // Clear memory
    mpfr_clear(pi);
    return 0;
}