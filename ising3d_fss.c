#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pcg_random.h"

typedef struct Ising {
    int spin;
    struct Ising *up, *down, *left, *right, *front, *back;
} Ising;

int N;
Ising ***lattice;

void initialize() {
    int i, j, k;

    // Allocate 3D lattice
    lattice = (Ising ***)malloc(N * sizeof(Ising **));
    for (i = 0; i < N; i++) {
        lattice[i] = (Ising **)malloc(N * sizeof(Ising *));
        for (j = 0; j < N; j++) {
            lattice[i][j] = (Ising *)malloc(N * sizeof(Ising));
        }
    }

    // Initialize random spins
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            for (k = 0; k < N; k++) {
                lattice[i][j][k].spin = (drnd() < 0.5) ? 1 : -1;
            }
        }
    }

    // Set up periodic boundary conditions (6 neighbors)
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            for (k = 0; k < N; k++) {
                lattice[i][j][k].up = &lattice[i][(j + 1) % N][k];
                lattice[i][j][k].down = &lattice[i][(j - 1 + N) % N][k];
                lattice[i][j][k].right = &lattice[(i + 1) % N][j][k];
                lattice[i][j][k].left = &lattice[(i - 1 + N) % N][j][k];
                lattice[i][j][k].front = &lattice[i][j][(k + 1) % N];
                lattice[i][j][k].back = &lattice[i][j][(k - 1 + N) % N];
            }
        }
    }
}

void mcstep(double T) {
    int i, j, k;
    double dE, r;
    Ising *spin;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            for (k = 0; k < N; k++) {
                spin = &lattice[i][j][k];

                // Calculate energy change if spin is flipped
                dE = 2.0 * spin->spin * (spin->up->spin + spin->down->spin +
                                         spin->left->spin + spin->right->spin +
                                         spin->front->spin + spin->back->spin);

                // Metropolis algorithm
                if (dE < 0.0) {
                    spin->spin = -spin->spin;
                } else {
                    r = drnd();
                    if (r < exp(-dE / T)) {
                        spin->spin = -spin->spin;
                    }
                }
            }
        }
    }
}

double magnetization() {
    int i, j, k;
    double m = 0.0;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            for (k = 0; k < N; k++) {
                m += lattice[i][j][k].spin;
            }
        }
    }

    return fabs(m) / (N * N * N);
}

double energy() {
    int i, j, k;
    double E = 0.0;
    Ising *spin;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            for (k = 0; k < N; k++) {
                spin = &lattice[i][j][k];
                // Count each pair only once (right, up, front neighbors)
                E -= spin->spin * (spin->right->spin + spin->up->spin + spin->front->spin);
            }
        }
    }

    return E / (N * N * N);
}

void cleanup() {
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            free(lattice[i][j]);
        }
        free(lattice[i]);
    }
    free(lattice);
}

int main() {
    double T, m, E, m_sum, E_sum, m2_sum, E2_sum, m4_sum;
    double chi, C, U;
    int i;
    int sizes[] = {4, 6, 8, 10, 12, 16};
    int num_sizes = 6;

    init_rnd(gus());

    printf("# N\tT\tM\tE\tchi\tC\tU\n");

    for (int s = 0; s < num_sizes; s++) {
        N = sizes[s];

        // Temperature range focused near Tc ~ 4.51 for 3D Ising
        for (T = 5.5; T > 3.5; T -= 0.02) {
            initialize();

            // Thermalization
            for (i = 0; i < 200000; i++) {
                mcstep(T);
            }

            // Measurement
            m_sum = E_sum = m2_sum = E2_sum = m4_sum = 0.0;
            for (i = 0; i < 200000; i++) {
                mcstep(T);
                m = magnetization();
                E = energy();
                m_sum += m;
                E_sum += E;
                m2_sum += m * m;
                E2_sum += E * E;
                m4_sum += m * m * m * m;
            }

            m_sum /= 200000.0;
            E_sum /= 200000.0;
            m2_sum /= 200000.0;
            E2_sum /= 200000.0;
            m4_sum /= 200000.0;

            chi = (m2_sum - m_sum * m_sum) * N * N * N / T;
            C = (E2_sum - E_sum * E_sum) * N * N * N / (T * T);

            // Binder cumulant: U_L = 1 - <M^4>/(3<M^2>^2)
            U = (m2_sum > 1e-10) ? (1.0 - m4_sum / (3.0 * m2_sum * m2_sum)) : 0.0;

            printf("%d\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                   N, T, m_sum, E_sum, chi, C, U);

            cleanup();
        }
        printf("\n");  // Blank line between different sizes
    }

    return 0;
}
