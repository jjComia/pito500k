#include <iostream>
#include <gmp.h>

using namespace std;

/* 
    Bailey-Borwein-Pouffe Algorithm Function.
    Uses GMP for arbitrary precision to calculate pi.
    Parameters: Cycle represents 'k' in the BBP 
                algorithm equation. 
                Term represents final calculation of 
                pi for a given iteration.
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
    Chudnovsky Algorithm Function.
    Uses GMP for arbitrary precision to calculate pi.
    Parameters: Cycle represents 'k' in the Chudnovsky 
                algorithm equation. 
                Term represents final calculation of 
                pi for a given iteration.
*/
void chudnovsky(mpf_t term, int cycle){ 
    mpz_t temp1, temp4, temp5, temp6, temp7;        //numerator will be int based
    mpf_t temp2, temp3, temp6_f, temp7_f, temp4_f;  //denominator will be float based
        
    //initialize
    mpz_inits(temp1, temp4, temp5, temp6, temp7, NULL);
    mpf_inits(temp2, temp3, temp4_f, temp6_f, temp7_f, NULL);

    //calculate large value in numerator
    mpz_set_ui(temp1, 545140134);
    mpz_mul_ui(temp1, temp1, cycle);
    mpz_add_ui(temp1, temp1, 13591409);
    
    //calculate 640320^3k+3/2 in denom.
    mpf_set_ui(temp2, 640320);
    mpf_pow_ui(temp2, temp2, 3 * cycle);    //640320^3k
    mpf_set_d(temp3, 640320.0);
    mpf_sqrt(temp3, temp3);                 //sqrt(640320) = 640320^0.5
    mpf_mul_ui(temp3, temp3, 640320);       //640320^0.5 * 640320 = 640320^1.5
    mpf_mul(temp2, temp2, temp3);           //640320^1.5 * 640320^3k = 640320^(3k+3/2)
    
    //calculating -1^k in numerator
    if(cycle % 2 == 0){
        mpz_set_ui(temp4, 1);
    }
    else{
        mpz_set_si(temp4, -1);
    }

    //numerator and denom. factorials
    mpz_fac_ui(temp5, 6*cycle);             //(6k)!
    mpz_fac_ui(temp6, 3*cycle);             //(3k)!
    mpz_fac_ui(temp7, cycle);               //(k!)^3
    mpz_pow_ui(temp7, temp7, 3);
    
    //numerator = temp4 or temp4_f
    mpz_mul(temp4, temp4, temp5);
    mpz_mul(temp4, temp4, temp1); 

    /*
    convert temp 4, 6 and 7 to floats 
    this makes it easier to combine
    floats and integers with gmp
    */
    mpf_set_z(temp6_f, temp6);
    mpf_set_z(temp7_f, temp7);
    mpf_set_z(temp4_f, temp4);              //numerator = temp4_f

    //denominator
    mpf_mul(temp6_f, temp6_f, temp7_f);
    mpf_mul(temp6_f, temp6_f, temp2);       //denominator = temp6_f

    //numerator / denominator
    mpf_div(term, temp4_f, temp6_f);
    mpf_mul_ui(term, term, 12);

    //clear GMPs
    mpf_clears(temp2, temp3, temp4_f, temp6_f, temp7_f, NULL);
    mpz_clears(temp1, temp4, temp5, temp6, temp7, NULL);
}

/*
    Calculate pi function using Chudnovsky Algorithm. 
    Notes:  If you want to run BBP, uncomment the bbp 
            function call, comment out chudnovsky, 
            change "n_digits/14" in the for loop to 
            "k < n_digits", and comment out inversion
*/
void calculate_pi(mpf_t pi, unsigned long n_digits) {
    mpf_t term;
    mpf_init(term);
    
    //Iterate in range 0 - n/14 (chudnovsky only requires n/14'th terms)
    for (unsigned long k = 0; k < n_digits/14; k++) {
        chudnovsky(term, k);
        //bbp(term,k);
        mpf_add(pi, pi, term);   //Add term to pi
    }

    //invert the sum when using chudnovsky
    mpf_ui_div(pi, 1, pi);
    mpf_clear(term);
}

int main() {
    unsigned long n_digits = 1000;                      //500000 for full test, 1000 to test functionality
    unsigned long precision = n_digits * 3.32193;       //globally set bit precision by multiplying log base_2 (10)

    mpf_t pi;
    mpf_init2(pi, precision);                           //initalizing pi with the gloabl precision
    mpf_set_ui(pi, 0);                                  //let pi be 0 to start

    calculate_pi(pi, n_digits);
    gmp_printf("Pi to 500000 digits: %.1000Ff\n", pi);  //This MUST match n_digits

    //Clear GMP variable
    mpf_clear(pi);

    return 0;
}