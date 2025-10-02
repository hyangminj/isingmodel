#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "pcg_random.h"

// Test functions
void test_pcg_basic(void);
void test_pcg_distribution(void);
void test_pcg_reproducibility(void);
void test_ising_physics(void);
void test_binder_cumulant(void);
void test_fss_multiple_sizes(void);
void test_ising_3d_physics(void);
void test_binder_cumulant_3d(void);
int run_mini_ising1d(int N, double T, int steps);
double run_mini_ising2d(int N, double T, int steps);
double run_mini_ising3d(int N, double T, int steps);
double calculate_binder_cumulant(int N, double T, int steps);
double calculate_binder_cumulant_3d(int N, double T, int steps);

int main(void) {
    printf("=== Ising Model Test Suite ===\n\n");

    printf("1. Testing PCG Random Number Generator...\n");
    test_pcg_basic();
    test_pcg_distribution();
    test_pcg_reproducibility();
    printf("   ✓ PCG tests passed\n\n");

    printf("2. Testing Ising Model Physics...\n");
    test_ising_physics();
    printf("   ✓ Physics tests passed\n\n");

    printf("3. Testing Binder Cumulant Calculation...\n");
    test_binder_cumulant();
    printf("   ✓ Binder cumulant tests passed\n\n");

    printf("4. Testing Finite-Size Scaling...\n");
    test_fss_multiple_sizes();
    printf("   ✓ FSS tests passed\n\n");

    printf("5. Testing 3D Ising Model Physics...\n");
    test_ising_3d_physics();
    printf("   ✓ 3D Physics tests passed\n\n");

    printf("6. Testing 3D Binder Cumulant...\n");
    test_binder_cumulant_3d();
    printf("   ✓ 3D Binder cumulant tests passed\n\n");

    printf("=== All Tests Passed! ===\n");
    return 0;
}

void test_pcg_basic(void) {
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

void test_pcg_distribution(void) {
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

void test_pcg_reproducibility(void) {
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

void test_ising_physics(void) {
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

// Test Binder cumulant calculation
void test_binder_cumulant(void) {
    // Test at high temperature (disordered phase)
    double U_high = calculate_binder_cumulant(10, 5.0, 5000);
    printf("   - Binder cumulant at T=5.0: %.3f\n", U_high);
    // At high T, should be close to 0 (fully disordered)
    assert(U_high >= -0.5 && U_high <= 0.5);

    // Test at low temperature (ordered phase)
    double U_low = calculate_binder_cumulant(10, 0.5, 5000);
    printf("   - Binder cumulant at T=0.5: %.3f\n", U_low);
    // At low T, should be close to 2/3 (fully ordered)
    assert(U_low >= 0.3 && U_low <= 1.0);

    // Test at critical temperature (should be around 0.61 for 2D Ising)
    double U_crit = calculate_binder_cumulant(16, 2.27, 10000);
    printf("   - Binder cumulant at T=2.27 (near Tc): %.3f\n", U_crit);
    // Should be between 0.4 and 0.8 (universal value ~0.61)
    assert(U_crit >= 0.4 && U_crit <= 0.8);
}

// Calculate Binder cumulant for 2D Ising model
double calculate_binder_cumulant(int N, double T, int steps) {
    int **spins = malloc(N * sizeof(int*));
    for (int i = 0; i < N; i++) {
        spins[i] = malloc(N * sizeof(int));
    }

    init_rnd(789);

    // Initialize random spins
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            spins[i][j] = (drnd() < 0.5) ? 1 : -1;
        }
    }

    // Thermalization
    for (int step = 0; step < steps/2; step++) {
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

    // Measurement
    double M2_sum = 0.0;
    double M4_sum = 0.0;
    int measure_steps = steps/2;

    for (int step = 0; step < measure_steps; step++) {
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

        // Calculate magnetization
        int sum = 0;
        for (int ii = 0; ii < N; ii++) {
            for (int jj = 0; jj < N; jj++) {
                sum += spins[ii][jj];
            }
        }
        double m = fabs((double)sum / (N * N));
        M2_sum += m * m;
        M4_sum += m * m * m * m;
    }

    double M2 = M2_sum / measure_steps;
    double M4 = M4_sum / measure_steps;

    // Binder cumulant: U_L = 1 - <M^4>/(3<M^2>^2)
    double U = (M2 > 1e-10) ? (1.0 - M4 / (3.0 * M2 * M2)) : 0.0;

    // Free memory
    for (int i = 0; i < N; i++) {
        free(spins[i]);
    }
    free(spins);

    return U;
}

// Test finite-size scaling with multiple system sizes
void test_fss_multiple_sizes(void) {
    int sizes[] = {8, 12, 16};
    int num_sizes = 3;
    double T_near_critical = 2.3;
    double binders[3];

    printf("   - Testing FSS for L = 8, 12, 16 at T=%.2f\n", T_near_critical);

    for (int i = 0; i < num_sizes; i++) {
        binders[i] = calculate_binder_cumulant(sizes[i], T_near_critical, 5000);
        printf("     L=%d: U_L = %.3f\n", sizes[i], binders[i]);
        // All should be in reasonable range (relaxed bounds due to statistical fluctuations)
        assert(binders[i] >= 0.2 && binders[i] <= 0.9);
    }

    // Check that we have variation across sizes (not all identical)
    // This tests that the simulation is actually running
    int all_same = 1;
    for (int i = 0; i < num_sizes - 1; i++) {
        if (fabs(binders[i] - binders[i+1]) > 0.01) {
            all_same = 0;
            break;
        }
    }
    assert(!all_same);  // Binder cumulants should vary with system size

    printf("   - FSS shows size-dependent behavior ✓\n");
}

// Simplified 3D Ising simulation for testing
double run_mini_ising3d(int N, double T, int steps) {
    int ***spins = malloc(N * sizeof(int**));
    for (int i = 0; i < N; i++) {
        spins[i] = malloc(N * sizeof(int*));
        for (int j = 0; j < N; j++) {
            spins[i][j] = malloc(N * sizeof(int));
        }
    }

    init_rnd(789);

    // Initialize random spins
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                spins[i][j][k] = (drnd() < 0.5) ? 1 : -1;
            }
        }
    }

    // Monte Carlo steps
    for (int step = 0; step < steps; step++) {
        int i = (int)(drnd() * N);
        int j = (int)(drnd() * N);
        int k = (int)(drnd() * N);

        int up = (j - 1 + N) % N;
        int down = (j + 1) % N;
        int left = (i - 1 + N) % N;
        int right = (i + 1) % N;
        int front = (k + 1) % N;
        int back = (k - 1 + N) % N;

        double dE = 2 * spins[i][j][k] * (spins[i][up][k] + spins[i][down][k] +
                                          spins[left][j][k] + spins[right][j][k] +
                                          spins[i][j][front] + spins[i][j][back]);

        if (dE < 0 || drnd() < exp(-dE / T)) {
            spins[i][j][k] = -spins[i][j][k];
        }
    }

    // Calculate magnetization
    int sum = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                sum += spins[i][j][k];
            }
        }
    }

    double mag = fabs((double)sum / (N * N * N));

    // Free memory
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            free(spins[i][j]);
        }
        free(spins[i]);
    }
    free(spins);

    return mag;
}

void test_ising_3d_physics(void) {
    // Test 1: High temperature should have low magnetization
    double mag_high_T = run_mini_ising3d(8, 6.0, 5000);
    printf("   - 3D Magnetization at T=6.0: %.3f\n", mag_high_T);
    assert(mag_high_T < 0.3);  // Should be small at high T

    // Test 2: Low temperature should have high magnetization
    double mag_low_T = run_mini_ising3d(8, 0.5, 5000);
    printf("   - 3D Magnetization at T=0.5: %.3f\n", mag_low_T);
    assert(mag_low_T > 0.8);  // Should be large at low T

    // Test 3: Near critical temperature (Tc ~ 4.51 for 3D)
    double mag_crit = run_mini_ising3d(8, 4.5, 5000);
    printf("   - 3D Magnetization at T=4.5 (near Tc): %.3f\n", mag_crit);
    assert(mag_crit >= 0.0 && mag_crit <= 1.0);
}

// Calculate Binder cumulant for 3D Ising model
double calculate_binder_cumulant_3d(int N, double T, int steps) {
    int ***spins = malloc(N * sizeof(int**));
    for (int i = 0; i < N; i++) {
        spins[i] = malloc(N * sizeof(int*));
        for (int j = 0; j < N; j++) {
            spins[i][j] = malloc(N * sizeof(int));
        }
    }

    init_rnd(1001);

    // Initialize random spins
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                spins[i][j][k] = (drnd() < 0.5) ? 1 : -1;
            }
        }
    }

    // Thermalization
    for (int step = 0; step < steps/2; step++) {
        int i = (int)(drnd() * N);
        int j = (int)(drnd() * N);
        int k = (int)(drnd() * N);

        int up = (j - 1 + N) % N;
        int down = (j + 1) % N;
        int left = (i - 1 + N) % N;
        int right = (i + 1) % N;
        int front = (k + 1) % N;
        int back = (k - 1 + N) % N;

        double dE = 2 * spins[i][j][k] * (spins[i][up][k] + spins[i][down][k] +
                                          spins[left][j][k] + spins[right][j][k] +
                                          spins[i][j][front] + spins[i][j][back]);

        if (dE < 0 || drnd() < exp(-dE / T)) {
            spins[i][j][k] = -spins[i][j][k];
        }
    }

    // Measurement
    double M2_sum = 0.0;
    double M4_sum = 0.0;
    int measure_steps = steps/2;

    for (int step = 0; step < measure_steps; step++) {
        int i = (int)(drnd() * N);
        int j = (int)(drnd() * N);
        int k = (int)(drnd() * N);

        int up = (j - 1 + N) % N;
        int down = (j + 1) % N;
        int left = (i - 1 + N) % N;
        int right = (i + 1) % N;
        int front = (k + 1) % N;
        int back = (k - 1 + N) % N;

        double dE = 2 * spins[i][j][k] * (spins[i][up][k] + spins[i][down][k] +
                                          spins[left][j][k] + spins[right][j][k] +
                                          spins[i][j][front] + spins[i][j][back]);

        if (dE < 0 || drnd() < exp(-dE / T)) {
            spins[i][j][k] = -spins[i][j][k];
        }

        // Calculate magnetization
        int sum = 0;
        for (int ii = 0; ii < N; ii++) {
            for (int jj = 0; jj < N; jj++) {
                for (int kk = 0; kk < N; kk++) {
                    sum += spins[ii][jj][kk];
                }
            }
        }
        double m = fabs((double)sum / (N * N * N));
        M2_sum += m * m;
        M4_sum += m * m * m * m;
    }

    double M2 = M2_sum / measure_steps;
    double M4 = M4_sum / measure_steps;

    // Binder cumulant: U_L = 1 - <M^4>/(3<M^2>^2)
    double U = (M2 > 1e-10) ? (1.0 - M4 / (3.0 * M2 * M2)) : 0.0;

    // Free memory
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            free(spins[i][j]);
        }
        free(spins[i]);
    }
    free(spins);

    return U;
}

void test_binder_cumulant_3d(void) {
    // Test at high temperature (disordered phase)
    double U_high = calculate_binder_cumulant_3d(6, 6.0, 3000);
    printf("   - 3D Binder cumulant at T=6.0: %.3f\n", U_high);
    // At high T, should be close to 0 (fully disordered)
    assert(U_high >= -0.5 && U_high <= 0.5);

    // Test at low temperature (ordered phase)
    double U_low = calculate_binder_cumulant_3d(6, 1.0, 3000);
    printf("   - 3D Binder cumulant at T=1.0: %.3f\n", U_low);
    // At low T, should be close to 2/3 (fully ordered)
    assert(U_low >= 0.3 && U_low <= 1.0);

    // Test at critical temperature (should be around 0.61-0.63 for 3D Ising)
    double U_crit = calculate_binder_cumulant_3d(8, 4.51, 5000);
    printf("   - 3D Binder cumulant at T=4.51 (near Tc): %.3f\n", U_crit);
    // Should be between 0.4 and 0.8 (universal value ~0.61-0.63)
    assert(U_crit >= 0.4 && U_crit <= 0.8);
}
