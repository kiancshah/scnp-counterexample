/*
 * bigbottom.h - shared by bigbottom.c and bigsplit.c: 128-lane vectors, the
 * relaxation step of sweep.c, and permutations packed 4 bits per letter
 * (n <= 16; for n >= 2 no permutation packs to 0, which marks an empty slot).
 */
#ifndef BIGBOTTOM_H
#define BIGBOTTOM_H

#include <stdint.h>
#include <string.h>

#define BN 17    /* largest n + 1 */
#define BV 128   /* lanes */

typedef uint8_t bvec_t __attribute__((vector_size(BV), aligned(BV)));

static const char BDIGITS[] = "123456789ABCDEFG";

static inline uint64_t pack(const uint8_t *p, int n)
{
    uint64_t c = 0;
    for (int i = 0; i < n; i++)
        c = c << 4 | p[i];
    return c;
}

static inline void unpack(uint64_t c, int n, uint8_t *p)
{
    for (int i = n - 1; i >= 0; i--) {
        p[i] = (uint8_t)(c & 15);
        c >>= 4;
    }
}

static inline uint64_t mix(uint64_t x)
{
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    return x ^ (x >> 33);
}

/* the vectors are handled in chunks of the target's width: GCC compiles comparisons of
 * vectors wider than that one lane at a time */
#if defined(__AVX512BW__)
#define CW 64
#elif defined(__AVX2__)
#define CW 32
#else
#define CW 16
#endif
typedef uint8_t cvec_t __attribute__((vector_size(CW)));

/* same as relax() in sweep.c, generic vector code */
static inline void brelax(bvec_t *h, uint8_t *f, bvec_t cand, uint8_t fc)
{
    cvec_t c[BV / CW], x[BV / CW], lt[BV / CW];
    uint64_t lo = 0, go = 0;
    memcpy(c, &cand, BV);
    memcpy(x, h, BV);
    for (int k = 0; k < BV / CW; k++) {
        cvec_t gt = (cvec_t)(c[k] > x[k]);
        lt[k] = (cvec_t)(c[k] < x[k]);
        uint64_t l[CW / 8], g[CW / 8];
        memcpy(l, &lt[k], CW);
        memcpy(g, &gt, CW);
        for (int i = 0; i < CW / 8; i++) {
            lo |= l[i];
            go |= g[i];
        }
    }
    if (!lo) {
        if (!go)
            *f |= fc;
        return;
    }
    if (!go) {
        *h = cand;
        *f = fc;
        return;
    }
    for (int k = 0; k < BV / CW; k++)
        x[k] = (x[k] & ~lt[k]) | (c[k] & lt[k]);
    memcpy(h, x, BV);
    *f = 0;
}

/* U from its one-line notation (letters 1..9, A..G); returns n, or -1 */
static inline int parse_big(const char *us, uint8_t *u)
{
    int n = 0;
    unsigned seen = 0;
    for (const char *s = us; *s; s++) {
        const char *c = strchr(BDIGITS, *s);
        if (!c || !*c || n >= BN - 1)
            return -1;
        u[n] = (uint8_t)(c - BDIGITS);
        seen |= 1u << u[n++];
    }
    return n >= 2 && seen == (1u << n) - 1 ? n : -1;
}

#endif
