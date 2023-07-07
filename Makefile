a.out: cdump.c
	gcc -Wall -Wextra -O3 -march=native cdump.c -lm

# && time ./a.out dump.gamma100_A2p100_N1000_f2.83_A0.317_mu0_kn1000_P1_I0.0001
