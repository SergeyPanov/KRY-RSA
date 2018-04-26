#include <iostream>
#include <random>
#include <gmp.h>
#include <unistd.h>
#include <string>
#include <assert.h>


#define ITERATIONS 1000000
#define K 100

// Wrapper for mpz_t
class mpz{
public:
    mpz_t value;

    explicit mpz(){
        mpz_init(value);
    };

    ~mpz(){
        mpz_clear(value);
    }

    void display(std::string name){
        std::cout << name << ": ";
        mpz_out_str(stdout, 10, value);
        std::cout << std:: endl;
    }

};

// Display
void show_hex(mpz &n){
    gmp_printf("%#Zx ", n.value);
}


// Store random bit_size odd number into r_number
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

    int p_bit_size = ceil((double)bit_size / 2);
    int q_bit_size = floor((double)bit_size / 2);

    set_prime(p_bit_size, p);   // Get p

    set_prime(q_bit_size, q);   // Get q

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
bool simple_factorization(mpz &n, mpz &p){

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
                return true;
            }
        }
        ++divider;
    }
    mpz_set_ui(p.value, 0);
    return false;
}


// m <- nid(a, b)
void min(mpz& a, mpz& b, mpz& m){
    if (mpz_cmp(a.value, b.value) > 0){
        mpz_set(m.value, b.value);
    } else{
        mpz_set(m.value, a.value);
    }
}

// Brent-Pollard pho implementation
// See http://maths-people.anu.edu.au/~brent/pd/rpb051i.pdf
void rho(mpz &n, mpz& res) {

    mpz mod_2;
    mpz_mod_ui(mod_2.value, n.value, 2);

    if (mpz_cmp_ui(mod_2.value, 0) == 0){
        mpz_set_ui(res.value, 2);
        return;
    }

    mpz y, c, m;
    mpz right_bound;

    mpz_sub_ui(right_bound.value, n.value, 1);  // right_bound <- N - 1

    set_random(1, right_bound.value, y.value);  // y <- [1 .. N - 1)


    set_random(1, right_bound.value, c.value);  // c <- [1 .. N - 1)

    set_random(1, right_bound.value, m.value);  // m <- [1 .. N - 1)


    mpz g, r, q;
    mpz_set_ui(g.value, 1); // g <- 1
    mpz_set_ui(r.value, 1); // r <- 1
    mpz_set_ui(q.value, 1); // q <- 1

    mpz ys, x;

    while (mpz_cmp_ui(g.value, 1) == 0){    // while gcd == 1

        mpz_set(x.value, y.value);  // x <- y

        mpz i;
        mpz_set_ui(i.value, 0); // i <- 0

        while (mpz_cmp(i.value, r.value) < 0){  // for i in [0 .. r)

            mpz_mul(y.value, y.value, y.value); // y <- y*y

            mpz_mod(y.value, y.value, n.value); // y <- y mod N

            mpz_add(y.value, y.value, c.value); // y <- y + c

            mpz_mod(y.value, y.value, n.value); // y <- y mod N

            mpz_add_ui(i.value, i.value, 1);
        }

        mpz k;
        mpz_set_ui(k.value, 0); // k <- 0

        while ( (mpz_cmp(k.value, r.value) < 0) && (mpz_cmp_ui(g.value, 1) == 0) ){ // while k < r && gdc == 1


            mpz_set(ys.value, y.value); // ys <- y

            mpz r_sub_k;
            mpz_sub(r_sub_k.value, r.value, k.value);   // r_sub_k <- r - k

            mpz min_acc;
            min(m, r_sub_k, min_acc);   // min_acc <- min(m, r - k)

            mpz ii;
            mpz_set_ui(ii.value, 0);    // ii <- 0

            while ( mpz_cmp(ii.value, min_acc.value) < 0 ){  // while i < min(m, r-k)

                mpz_mul(y.value, y.value, y.value); // y <- y*y
                mpz_mod(y.value, y.value, n.value); // y <- y mod N
                mpz_add(y.value, y.value, c.value); // y <- y + c
                mpz_mod(y.value, y.value, n.value); // y <- y mod N

                mpz x_sub_y;
                mpz_sub(x_sub_y.value, x.value, y.value);
                mpz_abs(x_sub_y.value, x_sub_y.value);  // abs(x-y)

                mpz_mul(q.value, q.value, x_sub_y.value);   // q <- q*abs(x-y)
                mpz_mod(q.value, q.value, n.value); // q <- q % N


                mpz_add_ui(ii.value, ii.value, 1);
            }

            set_gcd(q.value, n.value, g.value); // g <- gcd(N, q)
            mpz_add(k.value, k.value, m.value); // k <- k + m

        }
        mpz_mul_ui(r.value, r.value, 2);
    }


    if (mpz_cmp(g.value, n.value) == 0){    // if gcd == N

        while (true){

            mpz_mul(ys.value, ys.value, ys.value); // y <- y*y
            mpz_mod(ys.value, ys.value, n.value); // y <- y mod N
            mpz_add(ys.value, ys.value, c.value); // y <- y + c
            mpz_mod(ys.value, ys.value, n.value); // y <- y mod N

            mpz x_sub_ys;
            mpz_sub(x_sub_ys.value, x.value, ys.value);
            mpz_abs(x_sub_ys.value, x_sub_ys.value);    // |x - ys|

            set_gcd(x_sub_ys.value, n.value, g.value);

            if (mpz_cmp_ui(g.value, 1) > 0){
                break;
            }
        }
    }

    mpz_set(res.value, g.value);
}



// Execute factorization for n, write result to res
void RSA_cracker(mpz &n, mpz &res){
    mpz n_cpy;
    mpz_set(n_cpy.value,n.value);


    if (!simple_factorization(n, res)){
        rho(n, res);
    }
}


// Parse input parameters
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

    if (generate.size() == 2){
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

        std::cout << std::endl;

        return 0;
    }

    if (encrypt.size() == 4){
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

    if (decrypt.size() == 4){
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

    if (crack.size() == 2){
        mpz n, res;
        int flag;
        crack[1].erase(0, 2);
        flag = mpz_set_str(n.value, crack[1].c_str(), 16);
        assert(flag == 0);

        RSA_cracker(n, res);

        show_hex(res);
        std::cout << std::endl;
    }

    return 0;
}


