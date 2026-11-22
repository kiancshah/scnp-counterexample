/*
 * bigbottom.c - the fast criterion from one bottom U for n up to 16, without the
 * tables of S_n (which do not fit for n >= 12).
 *
 *   scnp bigbottom U [--prev BADW_S10.bin.gz] [--maxlen L]
 *
 * Same recursion as sweep.c (H[w] = componentwise min of g over chains U -> w,
 * F[w] = 1 iff one chain attains it), but the covers of a permutation are found
 * from its one-line notation, a level is a hash table keyed by the permutation,
 * and the vectors have 128 lanes (n(n-1)/2 <= 120 intervals); a level stores
 * only the HS = n(n-1)/2 (rounded up to 8) lanes that are used.
 *
 * Prints per length the number of pairs (U, w) and of failures.  With --prev
 * (the bad-w bitmap of S_10) every failing top w is tested for avoiding the
 * known patterns: w avoids them iff every 10-letter subsequence of w is good in
 * S_10 (all patterns have at most 10 letters and bad = contains a pattern there).
 * Failing pairs with an avoiding top are printed as "U w l(w)-l(U)".
 */
#include "scnp.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "bigbottom.h"

typedef struct {
    uint64_t key;  /* permutation, 4 bits per letter; 0 = empty slot */
    uint32_t idx;  /* position in the level arrays */
} slot_t;

static size_t HS = BV;  /* bytes of H stored per node */

typedef struct {
    size_t cap, n;          /* hash capacity (power of 2), number of nodes */
    slot_t *tab;
    uint64_t *perm;         /* node -> permutation */
    uint8_t *H;             /* node -> its HS bytes of H */
    uint8_t *F;
    size_t ncap;            /* capacity of the node arrays */
} level_t;

static inline bvec_t getH(const level_t *L, size_t i)
{
    bvec_t h = {0};
    memcpy(&h, L->H + i * HS, HS);
    return h;
}

static inline void putH(level_t *L, size_t i, bvec_t h)
{
    memcpy(L->H + i * HS, &h, HS);
}

static void *xmalloc(size_t sz)
{
    void *p = malloc(sz ? sz : 1);
    if (!p) {
        fprintf(stderr, "out of memory (%zu bytes)\n", sz);
        exit(1);
    }
    return p;
}

static void level_init(level_t *L, size_t expect)
{
    L->cap = 1024;
    while (L->cap < 2 * expect)
        L->cap <<= 1;
    L->tab = calloc(L->cap, sizeof(slot_t));
    L->ncap = expect > 16 ? expect : 16;
    L->perm = xmalloc(L->ncap * sizeof(uint64_t));
    L->H = xmalloc(L->ncap * HS);
    L->F = xmalloc(L->ncap);
    L->n = 0;
    if (!L->tab) {
        fprintf(stderr, "out of memory (hash)\n");
        exit(1);
    }
}

static void level_free(level_t *L)
{
    free(L->tab);
    free(L->perm);
    free(L->H);
    free(L->F);
}

static void level_grow(level_t *L)
{
    if (L->n + 1 > L->ncap) {           /* node arrays: grow by half */
        size_t nc = L->ncap + L->ncap / 2;
        L->perm = realloc(L->perm, nc * sizeof(uint64_t));
        L->H = realloc(L->H, nc * HS);
        L->F = realloc(L->F, nc);
        if (!L->perm || !L->H || !L->F) {
            fprintf(stderr, "out of memory (grow)\n");
            exit(1);
        }
        L->ncap = nc;
    }
    if (2 * (L->n + 1) > L->cap) {      /* hash: double and rehash */
        free(L->tab);
        L->cap = 2 * L->cap;
        L->tab = calloc(L->cap, sizeof(slot_t));
        if (!L->tab) {
            fprintf(stderr, "out of memory (hash)\n");
            exit(1);
        }
        for (size_t i = 0; i < L->n; i++) {
            size_t h = mix(L->perm[i]) & (L->cap - 1);
            while (L->tab[h].key)
                h = (h + 1) & (L->cap - 1);
            L->tab[h].key = L->perm[i];
            L->tab[h].idx = (uint32_t)i;
        }
    }
}

/* the node of permutation key, creating it (with *fresh = 1) if needed */
static size_t level_get(level_t *L, uint64_t key, int *fresh)
{
    if (2 * (L->n + 1) > L->cap || L->n + 1 > L->ncap)
        level_grow(L);
    size_t h = mix(key) & (L->cap - 1);
    while (L->tab[h].key) {
        if (L->tab[h].key == key) {
            *fresh = 0;
            return L->tab[h].idx;
        }
        h = (h + 1) & (L->cap - 1);
    }
    L->tab[h].key = key;
    L->tab[h].idx = (uint32_t)L->n;
    L->perm[L->n] = key;
    *fresh = 1;
    return L->n++;
}

static uint64_t rank10(const uint8_t *p)
{
    return perm_rank(p, 10);
}

/* w (size n >= 10) avoids the patterns iff no 10-letter subsequence is bad in S_10 */
static int avoids(const uint8_t *w, int n, const uint8_t *bad10)
{
    int pos[10];
    for (int i = 0; i < 10; i++)
        pos[i] = i;
    for (;;) {
        uint8_t q[10];
        for (int i = 0; i < 10; i++) {
            int r = 0;
            for (int j = 0; j < 10; j++)
                r += w[pos[j]] < w[pos[i]];
            q[i] = (uint8_t)r;
        }
        uint64_t r = rank10(q);
        if (bad10[r >> 3] >> (r & 7) & 1)
            return 0;
        int i = 9;
        while (i >= 0 && pos[i] == n - 10 + i)
            i--;
        if (i < 0)
            return 1;
        pos[i]++;
        for (int j = i + 1; j < 10; j++)
            pos[j] = pos[j - 1] + 1;
    }
}

int run_bigbottom(const char *us, const char *prevfile, int maxlen)
{
    uint8_t u[BN];
    int n = parse_big(us, u);
    if (n < 0 || n * (n - 1) / 2 > BV) {
        fprintf(stderr, "bad bottom %s\n", us);
        return 1;
    }
    HS = (size_t)(n * (n - 1) / 2 + 7) / 8 * 8;
    uint8_t *bad10 = NULL;
    if (prevfile && n >= 10) {
        size_t nb = 3628800 / 8;
        bad10 = malloc(nb);
        gzFile g = gzopen(prevfile, "rb");
        if (!g || gzread(g, bad10, (unsigned)nb) != (int)nb) {
            fprintf(stderr, "%s: cannot read the S_10 bitmap\n", prevfile);
            return 1;
        }
        gzclose(g);
    }

    /* interval coordinates and the step vectors delta[a][b] */
    static bvec_t delta[BN][BN];
    int ii[BV], ij[BV], nint = 0;
    for (int i = 0; i + 1 < n; i++)
        for (int j = i; j + 1 < n; j++) {
            ii[nint] = i;
            ij[nint] = j;
            nint++;
        }
    for (int a = 0; a < n; a++)
        for (int b = a + 1; b < n; b++) {
            uint8_t d[BV] = {0};
            for (int I = 0; I < nint; I++)
                d[I] = ii[I] <= a && b - 1 <= ij[I];
            memcpy(&delta[a][b], d, BV);
        }
    int N = n * (n - 1) / 2, lu = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            lu += u[i] > u[j];
    if (maxlen < 0 || maxlen > N - lu)
        maxlen = N - lu;

    level_t cur, nxt;
    level_init(&cur, 16);
    int fresh;
    size_t s0 = level_get(&cur, pack(u, n), &fresh);
    putH(&cur, s0, (bvec_t){0});
    cur.F[s0] = 1;
    uint64_t tp = 1, tf = 0, ta = 0;
    fprintf(stderr, "# bigbottom U = %s (n = %d, l(U) = %d): pairs (U, w) with w >= U\n# %4s %14s %14s %12s\n",
            us, n, lu, "len", "pairs", "SCNP fails", "avoiding");
    fprintf(stderr, "  %4d %14d %14d %12d\n", 0, 1, 0, 0);
    for (int d = 0; d < maxlen && cur.n; d++) {
        free(cur.tab);          /* the hash of the current level is no longer needed */
        cur.tab = NULL;
        level_init(&nxt, cur.n + cur.n / 2);
        for (size_t i = 0; i < cur.n; i++) {
            uint8_t p[BN];
            unpack(cur.perm[i], n, p);
            bvec_t hi = getH(&cur, i);
            for (int a = 0; a < n - 1; a++) {
                unsigned lim = (unsigned)n;
                for (int b = a + 1; b < n; b++)
                    if (p[b] > p[a] && p[b] < lim) {
                        lim = p[b];
                        uint8_t t = p[a];
                        p[a] = p[b];
                        p[b] = t;
                        size_t s = level_get(&nxt, pack(p, n), &fresh);
                        p[b] = p[a];
                        p[a] = t;
                        bvec_t cand = hi + delta[a][b];
                        if (fresh) {
                            putH(&nxt, s, cand);
                            nxt.F[s] = cur.F[i];
                        } else {
                            bvec_t h = getH(&nxt, s);
                            brelax(&h, &nxt.F[s], cand, cur.F[i]);
                            putH(&nxt, s, h);
                        }
                    }
            }
        }
        uint64_t fails = 0, avoid = 0;
        for (size_t s = 0; s < nxt.n; s++)
            if (!nxt.F[s]) {
                fails++;
                if (bad10) {
                    uint8_t w[BN];
                    unpack(nxt.perm[s], n, w);
                    if (avoids(w, n, bad10)) {
                        avoid++;
                        char a[BN], b[BN];
                        for (int k = 0; k < n; k++) {
                            a[k] = BDIGITS[u[k]];
                            b[k] = BDIGITS[w[k]];
                        }
                        a[n] = b[n] = 0;
                        printf("%s %s %d\n", a, b, d + 1);
                    }
                }
            }
        fprintf(stderr, "  %4d %14zu %14llu %12llu\n", d + 1, nxt.n, (unsigned long long)fails,
                (unsigned long long)avoid);
        tp += nxt.n;
        tf += fails;
        ta += avoid;
        level_free(&cur);
        cur = nxt;
    }
    level_free(&cur);
    fprintf(stderr, "  %4s %14llu %14llu %12llu\n", "all", (unsigned long long)tp,
            (unsigned long long)tf, (unsigned long long)ta);
    free(bad10);
    return 0;
}
