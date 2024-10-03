#include <iostream>
#include <gmp.h>

void chudnovsky(mpf_t term, int cycle){ 
    mpz_t temp1, temp4, temp5, temp6, temp7;        
    mpf_t temp2, temp3, temp6_f, temp7_f, temp4_f;  
        
    mpz_inits(temp1, temp4, temp5, temp6, temp7, NULL);
    mpf_inits(temp2, temp3, temp4_f, temp6_f, temp7_f, NULL);

    mpz_set_ui(temp1, 545140134);
    mpz_mul_ui(temp1, temp1, cycle);
    mpz_add_ui(temp1, temp1, 13591409);
    
    mpf_set_ui(temp2, 640320);
    mpf_pow_ui(temp2, temp2, 3 * cycle);    
    mpf_set_d(temp3, 640320.0);
    mpf_sqrt(temp3, temp3);              
    mpf_mul_ui(temp3, temp3, 640320);     
    mpf_mul(temp2, temp2, temp3);           
    
    if(cycle % 2 == 0){
        mpz_set_ui(temp4, 1);
    }
    else{
        mpz_set_si(temp4, -1);
    }


    mpz_fac_ui(temp5, 6*cycle);          
    mpz_fac_ui(temp6, 3*cycle);            
    mpz_fac_ui(temp7, cycle);              
    mpz_pow_ui(temp7, temp7, 3);

    mpz_mul(temp4, temp4, temp5);
    mpz_mul(temp4, temp4, temp1); 

    mpf_set_z(temp6_f, temp6);
    mpf_set_z(temp7_f, temp7);
    mpf_set_z(temp4_f, temp4);            

    mpf_mul(temp6_f, temp6_f, temp7_f);
    mpf_mul(temp6_f, temp6_f, temp2);      

    mpf_div(term, temp4_f, temp6_f);
    mpf_mul_ui(term, term, 12);

    mpf_clears(temp2, temp3, temp4_f, temp6_f, temp7_f, NULL);
    mpz_clears(temp1, temp4, temp5, temp6, temp7, NULL);
}

void calculate_pi(mpf_t pi, unsigned long n_digits) {
    mpf_t term;
    mpf_init(term);
    
    for (unsigned long k = 0; k < n_digits/14; k++) {
        chudnovsky(term, k);
        mpf_add(pi, pi, term);  
    }

    mpf_ui_div(pi, 1, pi);
    mpf_clear(term);
}

int main() {
    unsigned long n_digits = 1000;                   
    unsigned long precision = n_digits * 3.32193;    

    mpf_t pi;
    mpf_init2(pi, precision);                          
    mpf_set_ui(pi, 0);                                

    calculate_pi(pi, n_digits);
    gmp_printf("Pi to 500000 digits: %.1000Ff\n", pi);  

    mpf_clear(pi);

    return 0;
}