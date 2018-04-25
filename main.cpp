#include <iostream>
//#include "/Users/sergeypanov/bin/gmp/include/gmp.h"
#include <random>
#include <gmp.h>
#include <unistd.h>
#include <string>


#define ITERATIONS 1000000
#define K 100

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


void show_hex(mpz &n){
    gmp_printf("%#Zx ", n.value);
}

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

    // Set 2 MSB as 1
    mpz_setbit(r_number, bit_size - 1);
    mpz_setbit(r_number, bit_size - 2);


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
    mpz_set(old_a, a);

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

// Set dst to random number [min .. max - 1)
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
// Set r, d. See Miller-Rabin algorithm
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

// Check if number is prime. Used Miller-Rabin algorithm
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

    while (!miller_rabin_is_prime(number, K)){
        mpz_add_ui(number.value, number.value,2);
    }
}
// Calculate e
void set_e(int min, mpz &e, mpz &phi, int bit_size){
    mpz  r_n, gcd;

    set_random(min, phi.value, e.value);

    // Set 2 MSB as 1
    mpz_setbit(e.value, bit_size - 1);
    mpz_setbit(e.value, bit_size - 2);


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

            if (  (mpz_cmp_ui(gcd.value, 1) == 0) ){
                mpz_set(e.value, right.value);
                return;
            }
        }

        if (mpz_cmp_ui(left.value, 1) > 0){

            mpz_sub_ui(left.value, left.value, 1);

            set_gcd(left.value, phi.value, gcd.value);

            if ( (mpz_cmp_ui(gcd.value, 1) == 0)  ){
                mpz_set(e.value, left.value);
                return;
            }
        }
    }
}

//  Setup all used variables
void RSA_prepare(mpz &p, mpz &q, mpz &n, mpz &e, mpz &d, int bit_size){

    set_prime(bit_size/2, p);   // Get p

    set_prime(bit_size/2, q);   // Get q

    mpz_mul(n.value, p.value, q.value); // Get n

    mpz p_dec, q_dec;

    mpz_sub_ui(p_dec.value, p.value, 1);
    mpz_sub_ui(q_dec.value, q.value, 1);

    mpz phi;
    mpz_mul(phi.value, p_dec.value, q_dec.value);   // Get phi

    set_e(1, e, phi, bit_size);   // Get e

    set_inverse(e.value, phi.value, d.value);
}

// Cipher input message
void cipher(mpz &e, mpz &n, mpz &m, mpz &cm){
    mpz_powm(cm.value, m.value, e.value, n.value);
}
// Decipher
void decipher(mpz &d, mpz &n, mpz &c, mpz &m){
    mpz_powm(m.value, c.value, d.value, n.value);
}

// True if successfully False otherwise. p, q contains multipliers
bool simple_factorization(mpz &n, mpz &p, mpz &q){
    mpz reminder;
    unsigned int divider = 2;
    // Iters 1 000 000 times
    for (int i = 0; i < ITERATIONS; ++i) {

        mpz_mod_ui(reminder.value, n.value, divider);

        if (mpz_sgn(reminder.value) == 0) {
            mpz res;
            mpz mpz_divider;

            mpz_set_ui(mpz_divider.value, divider);

            mpz_div_ui(res.value, n.value, divider);

            if ( miller_rabin_is_prime(res, K) && miller_rabin_is_prime(mpz_divider, K) ){

                mpz_set(p.value, mpz_divider.value);
                mpz_set(q.value, res.value);
                return true;
            }
        }

        ++divider;
    }

    mpz_set_ui(p.value, 0);
    mpz_set_ui(q.value, 0);
    return false;
}


void pollar_rho_factorization_star(mpz &n, mpz &p, mpz &q){


    if (mpz_cmp_ui(n.value, 1) == 0){
        std:: cout << "n pho_factor" << std::endl;
        return;
    }

    mpz two_reminder;
    mpz_mod_ui(two_reminder.value, n.value, 2);
    if (mpz_sgn(two_reminder.value) == 0){
        std::cout << "obe factor is 2" << std::endl;
        return;
    }

    mpz x, y, d, c;
    set_random(2, n.value, x.value);    // Set x <- [2 .. N)

    mpz_set(y.value, x.value);  // Set y <- x


    set_random(1, n.value, c.value);    // Set c <- [2 .. N)

    mpz_set_ui(d.value, 1);

    while (mpz_cmp_ui(d.value, 1) == 0){

        mpz_powm_ui(x.value, x.value, 2, n.value);  // x <- x^2 mod n   x <- f(x)
        mpz aux_x;
        mpz_add(aux_x.value, c.value, n.value);   // aux_x <- c + n
        mpz_add(x.value, x.value, aux_x.value);   // x <- x + aux_x ~ (x^2 mod n) + c + n
        mpz_mod(x.value, x.value, n.value); // x <-  ((x^2 mod n) + c + n) mod n


        mpz_powm_ui(x.value, x.value, 2, n.value);
        mpz_powm_ui(x.value, x.value, 2, n.value);  // y <- f(f(y))
        mpz aux_y;
        mpz_add(aux_y.value, c.value, n.value);   // aux_y <- c + n
        mpz_add(y.value, y.value, aux_y.value);   // y <- y + aux_x ~ (y^2 mod n) + c + n
        mpz_mod(y.value, y.value, n.value); // x <-  ((y^2 mod n) + c + n) mod n


        mpz substract;

        mpz_sub(substract.value, x.value, y.value);
        mpz_abs(substract.value, substract.value);  // Calculate |x - y|

        set_gcd(substract.value, n.value, d.value); // Calculate GDC(|x-y|, n)



        if (mpz_cmp(d.value, n.value) == 0){
            std::cout << "recursion" << std::endl;
            pollar_rho_factorization_star(n, p, q);
        }
    }


    mpz_set(p.value, d.value);
    mpz_div(q.value, n.value, d.value);

}


void fermat_factorization_method(mpz &n){
    mpz n_cpy;
    mpz_set(n_cpy.value, n.value);

    mpz a;
    mpz_sqrt(a.value, n_cpy.value);
    mpz_add_ui(a.value, a.value, 1);    // ceil(sqrt(N))


    mpz a2;
    mpz_mul(a2.value, a.value, a.value);    // a2 <- a*a


    mpz b2;
    mpz_sub(b2.value, a2.value, n_cpy.value);   // b2 <- a*a - N

    while (mpz_perfect_square_p(b2.value) == 0){

        mpz_add_ui(a.value, a.value, 1);

        mpz_mul(a2.value, a.value, a.value);

        mpz_sub(b2.value, a2.value, n_cpy.value);
    }

    display("b2: ", b2.value);

    mpz p, b_sqr;

    mpz_sqrt(b_sqr.value, b2.value);

    display("sqrt b: ", b_sqr.value);
    display("a: ", a.value);

    mpz_add_ui(p.value, a.value, 1);

//    display("n: ", n.value);

    display("p: ", a.value);

}

void RSA_cracker(mpz &n){
    mpz n_cpy;
    mpz_set(n_cpy.value,n.value);

    fermat_factorization_method(n);


//    mpz p, q;
//    if (simple_factorization(n, p, q)){
//        display("p factorized by simple: ", p.value);
//        show_hex(p);
//        std::cout << std::endl;
////        display("q factorized by simple: ", q.value);
//
//    } else{
//        pollar_rho_factorization_star(n_cpy, p, q);
//        display("p factorized by rho: ", p.value);
//        show_hex(p);
//        std::cout << std::endl;
////        display("q factorized by rho: ", q.value);
//    }
}

std::vector< std::string > get_opts(int argc, char* argv[], char key){
    std::vector< std::string > params;

    bool should_read = false;
    char k;
    for (int i = 1; i < argc; ++i) {

        if (argv[i][0] == '-'){
            k = argv[i][1];
            should_read = !should_read;
        }

        if (k == key && should_read){
            params.push_back(argv[i]);
        }
    }

    return params;
}

int main(int argc, char* argv[]) {


    std::vector<std::string> generate = get_opts(argc, argv, 'g');

    std::vector<std::string> encrypt = get_opts(argc, argv, 'e');

    std::vector<std::string> decrypt = get_opts(argc, argv, 'd');

    std::vector<std::string> crack = get_opts(argc, argv, 'b');

    if (generate.size() > 0){
        // Generate
        int bit_size;

        mpz p, q, n, e, d;

        bit_size = std::stoi(generate[generate.size() - 1]);

        RSA_prepare(p, q, n, e, d, bit_size);

        show_hex(p);
        show_hex(q);
        show_hex(n);
        show_hex(e);
        show_hex(d);

//        std::cout << "e: " << mpz_sizeinbase(e.value, 2) << std::endl;
        std::cout << std::endl;

        return 0;
    }

    if (encrypt.size() > 0){
        // Encryption
        mpz e, n, m;

        int flag;

        encrypt[1].erase(0, 2); // Remove 0x
        flag = mpz_set_str(e.value, encrypt[1].c_str(), 16);
        assert(flag == 0);


        encrypt[2].erase(0, 2); // Remove 0x
        flag = mpz_set_str(n.value, encrypt[2].c_str(), 16);
        assert(flag == 0);

        encrypt[3].erase(0, 2); // Remove 0x
        flag = mpz_set_str(m.value, encrypt[3].c_str(), 16);
        assert(flag == 0);


        mpz cm;
        cipher(e, n, m, cm);


        show_hex(cm);
        std::cout << std::endl;

        return 0;
    }

    if (decrypt.size() > 0){
        // Decryption
        mpz d, n, c;

        int flag;

        decrypt[1].erase(0, 2);
        flag = mpz_set_str(d.value, decrypt[1].c_str(), 16);
        assert(flag == 0);

        decrypt[2].erase(0, 2);
        flag = mpz_set_str(n.value, decrypt[2].c_str(), 16);
        assert(flag == 0);

        decrypt[3].erase(0, 2);
        flag = mpz_set_str(c.value, decrypt[3].c_str(), 16);
        assert(flag == 0);

        mpz m;
        decipher(d, n, c, m);
        show_hex(m);
        std::cout << std::endl;

        return 0;

    }

    if (crack.size() > 0){
        mpz n;

        int flag;
        crack[1].erase(0, 2);
        flag = mpz_set_str(n.value, crack[1].c_str(), 16);
        assert(flag == 0);

//        std::cout << mpz_sizeinbase(n.value, 2) << std::endl;
        RSA_cracker(n);
    }

    return 0;
}


