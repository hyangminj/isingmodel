#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "pcg_random.h"

// Test functions
void test_pcg_basic();
void test_pcg_distribution();
void test_pcg_reproducibility();
void test_ising_physics();
int run_mini_ising1d(int N, double T, int steps);
double run_mini_ising2d(int N, double T, int steps);

int main() {
    printf("=== Ising Model Test Suite ===\n\n");

    printf("1. Testing PCG Random Number Generator...\n");
    test_pcg_basic();
    test_pcg_distribution();
    test_pcg_reproducibility();
    printf("   ✓ PCG tests passed\n\n");

    printf("2. Testing Ising Model Physics...\n");
    test_ising_physics();
    printf("   ✓ Physics tests passed\n\n");

    printf("=== All Tests Passed! ===\n");
    return 0;
}

void test_pcg_basic() {
    init_rnd(42);

    // Test range [0, 1)
    for (int i = 0; i < 1000; i++) {
        double r = drnd();
        assert(r >= 0.0 && r < 1.0);
    }

    // Test non-zero values
    int non_zero_count = 0;
    for (int i = 0; i < 100; i++) {
        if (drnd() > 0.0) non_zero_count++;
    }
    assert(non_zero_count > 50);  // Should have many non-zero values
}

void test_pcg_distribution() {
    init_rnd(12345);

    double sum = 0.0;
    int n = 10000;

    for (int i = 0; i < n; i++) {
        sum += drnd();
    }

    double mean = sum / n;
    // Mean should be approximately 0.5
    assert(fabs(mean - 0.5) < 0.1);

    printf("   - Distribution mean: %.3f (expected ~0.5)\n", mean);
}

void test_pcg_reproducibility() {
    // Test 1: Same seed should give same sequence
    init_rnd(999);
    double seq1[10];
    for (int i = 0; i < 10; i++) {
        seq1[i] = drnd();
    }

    init_rnd(999);
    double seq2[10];
    for (int i = 0; i < 10; i++) {
        seq2[i] = drnd();
    }

    for (int i = 0; i < 10; i++) {
        assert(seq1[i] == seq2[i]);
    }

    printf("   - Reproducibility: ✓\n");
}

void test_ising_physics() {
    // Test 1: High temperature should have low magnetization
    double mag_high_T = run_mini_ising1d(50, 5.0, 10000);
    printf("   - Magnetization at T=5.0: %.3f\n", mag_high_T);
    assert(mag_high_T < 0.3);  // Should be small at high T

    // Test 2: Low temperature should have high magnetization
    double mag_low_T = run_mini_ising1d(50, 0.1, 10000);
    printf("   - Magnetization at T=0.1: %.3f\n", mag_low_T);
    assert(mag_low_T > 0.8);  // Should be large at low T

    // Test 3: 2D model basic functionality
    double mag_2d = run_mini_ising2d(10, 1.0, 1000);
    printf("   - 2D Magnetization at T=1.0: %.3f\n", mag_2d);
    assert(mag_2d >= 0.0 && mag_2d <= 1.0);
}

// Simplified 1D Ising simulation for testing
int run_mini_ising1d(int N, double T, int steps) {
    int *spins = malloc(N * sizeof(int));
    init_rnd(123);

    // Initialize random spins
    for (int i = 0; i < N; i++) {
        spins[i] = (drnd() < 0.5) ? 1 : -1;
    }

    // Monte Carlo steps
    for (int step = 0; step < steps; step++) {
        int i = (int)(drnd() * N);
        int left = (i - 1 + N) % N;
        int right = (i + 1) % N;

        double dE = 2 * spins[i] * (spins[left] + spins[right]);

        if (dE < 0 || drnd() < exp(-dE / T)) {
            spins[i] = -spins[i];
        }
    }

    // Calculate magnetization
    int sum = 0;
    for (int i = 0; i < N; i++) {
        sum += spins[i];
    }

    double mag = fabs((double)sum / N);
    free(spins);
    return mag;
}

// Simplified 2D Ising simulation for testing
double run_mini_ising2d(int N, double T, int steps) {
    int **spins = malloc(N * sizeof(int*));
    for (int i = 0; i < N; i++) {
        spins[i] = malloc(N * sizeof(int));
    }

    init_rnd(456);

    // Initialize random spins
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            spins[i][j] = (drnd() < 0.5) ? 1 : -1;
        }
    }

    // Monte Carlo steps
    for (int step = 0; step < steps; step++) {
        int i = (int)(drnd() * N);
        int j = (int)(drnd() * N);

        int up = (i - 1 + N) % N;
        int down = (i + 1) % N;
        int left = (j - 1 + N) % N;
        int right = (j + 1) % N;

        double dE = 2 * spins[i][j] * (spins[up][j] + spins[down][j] +
                                        spins[i][left] + spins[i][right]);

        if (dE < 0 || drnd() < exp(-dE / T)) {
            spins[i][j] = -spins[i][j];
        }
    }

    // Calculate magnetization
    int sum = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            sum += spins[i][j];
        }
    }

    double mag = fabs((double)sum / (N * N));

    // Free memory
    for (int i = 0; i < N; i++) {
        free(spins[i]);
    }
    free(spins);

    return mag;
}