#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>

int lucas_lehmer(int p) {
    if (p == 2) return 1; // Special case where 2^2 - 1 = 3, which is prime

    mpz_t s, m;
    mpz_init_set_ui(s, 4); // Start with S_0 = 4
    mpz_init(m);
    mpz_ui_pow_ui(m, 2, p); // Calculate 2^p
    mpz_sub_ui(m, m, 1); // m = 2^p - 1, the Mersenne number

    time_t start_time, current_time;
    time(&start_time);

    for (int i = 1; i < p - 1; i++) {
        mpz_mul(s, s, s); // s = s^2
        mpz_sub_ui(s, s, 2); // s = s - 2
        mpz_mod(s, s, m); // s = s % m

        // Update progress every 30 seconds
        time(&current_time);
        if (difftime(current_time, start_time) >= 30) {
            printf("\r%d/%d iterations completed", i, p - 2);
            fflush(stdout);
            start_time = current_time; // reset the timer
        }
    }

    int result = mpz_cmp_ui(s, 0) == 0;

    mpz_clears(s, m, NULL);
    return result; // 1 if prime, 0 if not
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <exponent>\n", argv[0]);
        exit(1);
    }

    int p = atoi(argv[1]);
    if (p < 2) {
        fprintf(stderr, "Exponent must be greater than 1.\n");
        exit(1);
    }

    int result = lucas_lehmer(p);
    if (result) {
        printf("\n2^%d - 1 is a Mersenne prime.\n", p);
    } else {
        printf("\n2^%d - 1 is not a Mersenne prime.\n", p);
    }

    return 0;
}
