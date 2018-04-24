#include <iostream>
#include "/Users/sergeypanov/bin/gmp/include/gmp.h"
#include <random>


void display(std::string msg, mpz_t number){
    std::cout << msg;
    mpz_out_str(stdout, 10, number);
    std::cout << std::endl;
}

// Store random bit_size number into r_number
void set_random_odd(int bit_size, mpz_t r_number){

    srand((int)(time(0)));

    int seed = rand();  // Setup seed

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> dist(3, seed);

    seed = dist(mt);

    gmp_randstate_t r_state;
    gmp_randinit_default(r_state);
    gmp_randseed_ui(r_state, seed);

    mpz_urandomb(r_number, r_state, bit_size);
    mpz_setbit(r_number, bit_size - 1);


    mpz_t reminder;
    mpz_init(reminder);

    mpz_t one;
    mpz_init(one);
    mpz_set_ui(one, 1);




    mpz_t two;
    mpz_init(two);
    mpz_set_ui(two, 2);

    mpz_mod(reminder, r_number, two);


    if (mpz_sgn(reminder) == 0){
        mpz_add(r_number, r_number, one);
    }

    mpz_clear(one);
    mpz_clear(two);
    mpz_clear(reminder);
    gmp_randclear(r_state);
}

// Check if number is prime
bool check_is_prime(mpz_t number){
    bool is_prime = true;

    mpz_t divisor;
    mpz_init(divisor);

    mpz_set_ui(divisor, 2);

    mpz_t quotient;
    mpz_init(quotient);

    mpz_div(quotient, number, divisor); // Div by 2

    mpz_t inc;
    mpz_init(inc);
    mpz_set_ui(inc, 1);


    mpz_t reminder;
    mpz_init(reminder);



    while (mpz_cmp(quotient, divisor) > 0){

        mpz_mod(reminder, number, divisor);

        if (mpz_sgn(reminder) == 0){
            is_prime = false;
            break;
        }

        mpz_add(divisor, divisor, inc);
    }


    mpz_clear(inc);
    mpz_clear(divisor);
    mpz_clear(reminder);
    mpz_clear(quotient);

    return is_prime;
}


// Set number to prime with bit_size bits
void set_prime(int bit_size, mpz_t number){
    set_random_odd(bit_size, number);
//    mpz_nextprime(number, number);
}

// Calculate GCD, store it into res
void set_gcd (mpz_t a, mpz_t b, mpz_t res) {

    mpz_t local_b;
    mpz_init(local_b);
    mpz_set(local_b, b);

    mpz_set(res, a);

    while (mpz_sgn(local_b) > 0){
        mpz_mod(res, res, local_b);
        mpz_swap(res, local_b);
    }
    mpz_clear(local_b);

}



void extended_euclid_gmp(mpz_t a, mpz_t b, mpz_t x, mpz_t y, mpz_t d) {
    mpz_t q, r, x1, x2, y1, y2;

    mpz_init(q);
    mpz_init(r);
    mpz_init(x1);
    mpz_init(x2);
    mpz_init(y1);
    mpz_init(y2);
    if (mpz_sgn(b) == 0) {

        mpz_set(d, a);
        mpz_set_ui(x, 1);
        mpz_set_ui(y, 0);
        return;
    }

    mpz_set_ui(x2, 1);
    mpz_set_ui(x1, 0);
    mpz_set_ui(y2, 0);
    mpz_set_ui(y1, 1);



    mpz_t aux;
    mpz_init(aux);

    while (mpz_sgn(b) > 0) {

        mpz_div(q, a, b);

        mpz_mul(aux, q, b);
        mpz_sub(r, a, aux);

        mpz_mul(aux, q, x1);
        mpz_sub(x, x2, aux);


        mpz_mul(aux, q, y1);
        mpz_sub(y, y2, aux);


        mpz_set(a, b);
        mpz_set(b, r);

        mpz_set(x2, x1);
        mpz_set(x1, x);
        mpz_set(y2, y1);
        mpz_set(y1, y);
    }

    mpz_set(d, a);
    mpz_set(x, x2);
    mpz_set(y, y2);

    mpz_clear(aux);
    mpz_clear(q);
    mpz_clear(r);
    mpz_clear(x1);
    mpz_clear(x2);
    mpz_clear(y1);
    mpz_clear(y2);
}


void inverse_gmp(mpz_t a, mpz_t n, mpz_t x) {
    mpz_t d, y;

    mpz_init(d);
    mpz_init(y);


    extended_euclid_gmp(a, n, x, y, d);

    if (mpz_cmp_ui(d, 1) != 0) mpz_set_ui(x, 0);

    mpz_clear(d);
    mpz_clear(y);
}


///* вычисление a * *x + b * *y = gcd(a, b) = *d */
void extended_euclid(long a, long b, long *x, long *y, long *d) {
    long q, r, x1, x2, y1, y2;
    if (b == 0) {
        *d = a, *x = 1, *y = 0;
        return;
    }

    x2 = 1, x1 = 0, y2 = 0, y1 = 1;
    while (b > 0) {
        q = a / b, r = a - q * b;
        *x = x2 - q * x1, *y = y2 - q * y1;
        a = b, b = r;
        x2 = x1, x1 = *x, y2 = y1, y1 = *y;
    }

    *d = a, *x = x2, *y = y2;
}



long inverse(long a, long n) {
    long d, x, y;
    extended_euclid(a, n, &x, &y, &d);
    if (d == 1) return x;
    return 0;
}




bool solovay_strassen_prime_check(mpz_t number){

    bool is_prime = true;

}

void RSA_prepare(mpz_t p, mpz_t q, mpz_t n, mpz_t e, mpz_t d, int bit_size){

    mpz_t phi;
    mpz_init(phi);


    set_prime(bit_size/2, p);   // Get p

    display("p: ", p);


    set_prime(bit_size/2, q);   //Get q
    display("q: ", q);

    set_gcd(p, q, phi);

    display("gcd: ", phi);

}


int main() {


    mpz_t a, m, res;
    mpz_init(a);
    mpz_init(m);
    mpz_init(res);


    mpz_set_ui(a, 10);
    mpz_set_ui(m, 12);

    inverse_gmp(a, m, res);

    display("invers: ", res);

    std::cout << inverse(10, 12) << std::endl;

    mpz_clear(a);
    mpz_clear(m);
    mpz_clear(res);

    mpz_t p, q, n, e, d;
    mpz_init(p);
    mpz_init(q);
    mpz_init(n);
    mpz_init(e);
    mpz_init(d);


//    RSA_prepare(p, q, n, e, d, 96);

    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(n);
    mpz_clear(e);
    mpz_clear(d);


//
//    char input[1024];
//
//    std::cout << "Enter your number" << std:: endl;
//    std::cin >> input;
//
//    mpz_t n;
//    int flag;
//
//    mpz_init(n);
//    mpz_set_ui(n, 0);
//
//    flag = mpz_set_str(n, input, 10);
//    assert(flag == 0);
//
//    std::cout << "n = ";
//    mpz_out_str(stdout, 10, n);
//    std::cout << std::endl;
//
//
//    mpz_mul(n, n, n);
//    std::cout << "n*n = ";
//    mpz_out_str(stdout, 10, n);
//    std::cout << std::endl;
//
//    mpz_clear(n);
    return 0;
}