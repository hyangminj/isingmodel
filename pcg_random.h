#ifndef PCG_RANDOM_H
#define PCG_RANDOM_H

#include <stdint.h>
#include <time.h>

// PCG Random Number Generator
// Based on PCG-XSH-RR variant, 32-bit output
typedef struct {
    uint64_t state;
    uint64_t inc;
} pcg32_random_t;

// Global PCG state
static pcg32_random_t pcg32_global = {0x853c49e6748fea9bULL, 0xda3e39cb94b95bdbULL};

// Generate random 32-bit integer
static inline uint32_t pcg32_random(void) {
    uint64_t oldstate = pcg32_global.state;
    pcg32_global.state = oldstate * 6364136223846793005ULL + pcg32_global.inc;
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

// Initialize PCG with seed
static inline void pcg32_srandom(uint64_t seed) {
    pcg32_global.state = 0U;
    pcg32_global.inc = (seed << 1u) | 1u;
    pcg32_random();
    pcg32_global.state += seed;
    pcg32_random();
}

// Generate random double in [0, 1)
static inline double pcg32_random_double(void) {
    return pcg32_random() * (1.0 / 4294967296.0);
}

// Function aliases to match twist.c interface
static inline void init_rnd(uint64_t seed) {
    pcg32_srandom(seed);
}

static inline double drnd(void) {
    return pcg32_random_double();
}

static inline uint64_t gus(void) {
    return (uint64_t)time(NULL);
}

#endif // PCG_RANDOM_H