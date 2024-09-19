#include <iostream>
#include <gmp.h>

using namespace std;

/* 
    Bailey-Borwein-Pouffe Algorithm Function.
    Uses GMP for arbitrary precision to calculate pi.
    Cycle represents 'k' in the BBP algorithm equation.
    Term represents final calculation of Pi.
*/
void bbp(mpf_t term, int cycle){
    mpf_t temp1, temp2, temp3, temp4, divisor; 
    mpf_inits(temp1, temp2, temp3, temp4, divisor, NULL);
    
    //Calculate (1/16^k)
    mpf_set_ui(temp1, 16);
    mpf_pow_ui(temp1, temp1, cycle);
    mpf_ui_div(temp1, 1, temp1);
    
    //Calculate BBP terms 4/(8'k'+1), 2(8'k'+2),...
    mpf_set_ui(temp2, 8*cycle+1);
    mpf_ui_div(temp2, 4, temp2); 

    mpf_set_ui(temp3, 8*cycle+4);
    mpf_ui_div(temp3, 2, temp3);

    mpf_set_ui(temp4, 8*cycle+5);
    mpf_ui_div(temp4, 1, temp4);

    mpf_set_ui(divisor, 8*cycle+6);
    mpf_ui_div(divisor, 1, divisor);

    //term = (temp2 - temp3 - temp4 - divisor) * temp1
    mpf_sub(temp2, temp2, temp3);
    mpf_sub(temp2, temp2, temp4);
    mpf_sub(temp2, temp2, divisor);
    mpf_mul(term, temp1, temp2);

    //clear GMP variables
    mpf_clears(temp1, temp2, temp3, temp4, divisor, NULL);
}

/*
    Calculate pi using BBP Algorithm. 
    Set precision using log base_2 (10)
*/
void calculate_pi(mpf_t pi, unsigned long n_digits) {
    unsigned long precision = n_digits * 3.32193;
    mpf_set_default_prec(precision);

    mpf_t term;
    mpf_init(term);
    
    //Iterate through BBP series
    for (unsigned long k = 0; k < n_digits; k++) {
        bbp(term, k);
        mpf_add(pi, pi, term);   //Add term to pi
    }

    mpf_clear(term);
}

int main() {
    //set to 1000 to test first, MAKE SURE PRINTF MATCHES THIS NUMBER
    unsigned long n_digits = 1000;  //500000 for full test
    //multiply log base_2 (10) to add binary precision
    unsigned long precision = n_digits * 3.32193;

    mpf_t pi;
    mpf_init2(pi, precision);
    mpf_set_ui(pi, 0); //let pi be 0 to start

    calculate_pi(pi, n_digits);
    gmp_printf("Pi to 500000 digits: %.1000Ff\n", pi);    //UPDATE TO MATCH N_DIGITS

    //Clear GMP variable
    mpf_clear(pi);

    return 0;
}