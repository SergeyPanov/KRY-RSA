#include <iostream>
#include "/Users/sergeypanov/bin/gmp/include/gmp.h"
#include <random>

class mpz{
public:
    mpz_t value;

    explicit mpz(){
        mpz_init(value);
    };

    ~mpz(){
        mpz_clear(value);
    }

};

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


// Extended Euclidian algorithm for finding inverse
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

// Find inverse for a^(-1) mod n; store it into x
void set_inverse(mpz_t a, mpz_t n, mpz_t x) {
    mpz_t d, y, old_n, old_a;

    mpz_init(d);
    mpz_init(y);
    mpz_init(old_n);
    mpz_init(old_a);

    mpz_set(old_n, n);

    extended_euclid_gmp(old_a, old_n, x, y, d);

    if (mpz_cmp_ui(d, 1) != 0) mpz_set_ui(x, 0);

    if (mpz_sgn(x) == -1) {
        mpz_add(x, n, x);
    }

    mpz_clear(d);
    mpz_clear(y);
    mpz_clear(old_n);
    mpz_clear(old_a);
}


void set_random(int min, mpz_t max, mpz_t dst){


    srand((int)(time(0)));

    int seed = rand();  // Setup seed

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> dist(min, seed);

    seed = dist(mt);


    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, seed);

    mpz_urandomm(dst, state, max);


    if (mpz_cmp_ui(dst, min) < 0){
        mpz_set_ui(dst, min);
    }

    gmp_randclear(state);
}

void set_r_d(mpz &number, mpz &r, mpz &d){

    mpz acc, reminder;


    mpz_set(acc.value, number.value);
    mpz_set_ui(r.value, 0);

    mpz_sub_ui(acc.value, acc.value, 1);    // n - 1


    do{

        mpz_div_ui(acc.value, acc.value, 2);

        mpz_add_ui(r.value, r.value, 1);

        mpz_mod_ui(reminder.value, acc.value, 2);


    }while (mpz_sgn(reminder.value) == 0);

    mpz_set(d.value, acc.value);

}


bool miller_rabin_is_prime(mpz &number, int k){
    mpz r, d, right_bound, a;

    set_r_d(number, r, d);

    for (int i = 0; i < k; ++i) {
        bool is_witness_loop = false;

        mpz n_cpy_dec;

        mpz_set(n_cpy_dec.value, number.value);
        mpz_sub_ui(n_cpy_dec.value, n_cpy_dec.value, 1);    // n_copy <- n - 1

        mpz a;
        set_random(2, n_cpy_dec.value, a.value);    // pick random a from [2 .. n-1)

        mpz x;
        mpz_powm(x.value, a.value, d.value, number.value);  // a <- a^d mod n

        if ( (mpz_cmp_ui(x.value, 1) == 0) ||  (mpz_cmp(x.value, n_cpy_dec.value) == 0)   ){
            continue;
        }

        mpz r_dec;
        mpz_set(r_dec.value, r.value);
        mpz_sub_ui(r_dec.value, r_dec.value, 1);    // r - 1

        while (mpz_cmp_ui(r_dec.value, 0) > 0){


            mpz_powm_ui(x.value, x.value, 2, number.value); // x <- x^2 mod n

            if (mpz_cmp_ui(x.value, 1) == 0){   // if x = 1
                return false;
            }

            if (mpz_cmp(x.value, n_cpy_dec.value) == 0){
                is_witness_loop = true;
                break;
            }

            mpz_sub_ui(r_dec.value, r_dec.value, 1);
        }

        if (is_witness_loop){
            continue;
        }
        return false;
    }
    return true;
}


// Set number to prime with bit_size bits
void set_prime(int bit_size, mpz &number){

    set_random_odd(bit_size, number.value);

    while (!miller_rabin_is_prime(number, 100)){
        mpz_add_ui(number.value, number.value,2);
    }
}

void set_e(int min, mpz &e, mpz &phi){

    mpz  r_n, gcd, init_gcd;
    mpz_set(init_gcd.value, gcd.value);

    set_random(min, phi.value, e.value);

    set_gcd(e.value, phi.value, gcd.value);

    mpz left, right;


    if (mpz_cmp_ui(gcd.value, 1) == 0){
        return;
    }

    mpz_set(left.value, e.value);

    mpz_set(right.value, e.value);

    mpz left_e, right_e;

    while (true){

        if (mpz_cmp(right.value, phi.value) < 0){

            mpz_add_ui(right.value, right.value, 1);

            set_gcd(right.value, phi.value, gcd.value);

            if (mpz_cmp_ui(gcd.value, 1) == 0){
                mpz_set(e.value, right.value);
                return;
            }
        }

        if (mpz_cmp_ui(left.value, 1) > 0){

            mpz_sub_ui(left.value, left.value, 1);

            set_gcd(left.value, phi.value, gcd.value);

            if (mpz_cmp_ui(gcd.value, 1) == 0){
                mpz_set(e.value, left.value);
                return;
            }
        }

    }

}


void RSA_prepare(mpz &p, mpz &q, mpz &n, mpz &e, mpz &d, int bit_size){

    set_prime(bit_size/2, p);   // Get p
    display("p: ", p.value);

    set_prime(bit_size/2, q);   // Get q
    display("q: ", q.value);

    mpz_mul(n.value, p.value, q.value); // Get n
    display("n: ", n.value);


    mpz p_dec, q_dec;

    mpz_sub_ui(p_dec.value, p.value, 1);
    mpz_sub_ui(q_dec.value, q.value, 1);

    mpz phi;
    mpz_mul(phi.value, p_dec.value, q_dec.value);   // Get phi

    display("phi: ", phi.value);

    set_e(1, e, phi);   // Get e

    display("e: ", e.value);


    set_inverse(e.value, phi.value, d.value);

    display("d: ", d.value);


}




int main() {

    mpz p, q, n, e, d;


    RSA_prepare(p, q, n, e, d, 96);



    return 0;
}



///* вычисление a * *x + b * *y = gcd(a, b) = *d */
//void extended_euclid(long a, long b, long *x, long *y, long *d) {
//    long q, r, x1, x2, y1, y2;
//    if (b == 0) {
//        *d = a, *x = 1, *y = 0;
//        return;
//    }
//
//    x2 = 1, x1 = 0, y2 = 0, y1 = 1;
//    while (b > 0) {
//        q = a / b, r = a - q * b;
//        *x = x2 - q * x1, *y = y2 - q * y1;
//        a = b, b = r;
//        x2 = x1, x1 = *x, y2 = y1, y1 = *y;
//    }
//
//    *d = a, *x = x2, *y = y2;
//}



//long inverse(long a, long n) {
//    long d, x, y;
//    extended_euclid(a, n, &x, &y, &d);
//    if (d == 1) return x;
//    return 0;
//}

