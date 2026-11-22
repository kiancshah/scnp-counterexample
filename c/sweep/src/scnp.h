/*
 * scnp.h - shared declarations for the SCNP checker.
 *
 * Conventions (they reproduce the n=6 example u=214365, w=541632):
 *   - permutations of {0..n-1} in one-line notation, printed 1-based;
 *   - covers v < v*t_ab (right multiplication = swap of POSITIONS a<b);
 *   - the step v -> v*t_ab has weight (x_a - x_b) = y_a + ... + y_{b-1},
 *     i.e. the y_k are the simple-root (fundamental-weight) coordinates.
 * D_{u,w} = (1/l!) sum over saturated chains C of m_C, m_C = prod of weights.
 */
#ifndef SCNP_H
#define SCNP_H

#include <stdint.h>
#include <stddef.h>

#define MAXN 11
#define MAXT 55   /* transpositions t_ab (a<b) for n = 11              */
#define MAXL 56   /* interval lengths 0..55 for n = 11                 */
#define VL   64   /* vector lanes; must be >= #intervals [i,j] (55)    */

typedef uint8_t vec_t __attribute__((vector_size(VL), aligned(VL)));
typedef int8_t svec_t __attribute__((vector_size(VL), aligned(VL)));

typedef struct {
    int n;               /* permutations of {0..n-1}                       */
    int N;               /* n(n-1)/2 = length of w0                        */
    int ntr;             /* number of transpositions t_ab                  */
    int nint;            /* number of intervals [i,j], 0 <= i <= j <= n-2  */
    uint32_t nperm;      /* n!                                             */
    uint32_t maxlevel;   /* size of the largest rank level                 */
    uint32_t id, w0;     /* indices of the identity and of w0              */
    uint8_t *perm;       /* nperm*n one-line notation (lexicographic rank) */
    uint8_t *len;        /* Coxeter length                                 */
    uint32_t *off;       /* upper covers of v: cov[off[v] .. off[v+1])     */
    uint32_t *cov;       /* index of the cover v*t_ab                      */
    uint8_t *ctr;        /* transposition id t of that cover               */
    uint32_t *lw0;       /* index of w0*v      (values complemented)       */
    uint32_t *rw0;       /* index of v*w0      (positions reversed)        */
    uint32_t *cw0;       /* index of w0*v*w0                               */
    uint8_t ta[MAXT], tb[MAXT];  /* positions a < b of transposition t     */
    uint8_t ii[VL], ij[VL];      /* interval number I is [ii[I], ij[I]]    */
    vec_t delta[MAXT];   /* delta[t][I] = 1 iff [a, b-1] is contained in I */
} tables_t;

/* tables.c */
uint32_t perm_rank(const uint8_t *p, int n);
void perm_unrank(uint32_t r, int n, uint8_t *p);
int perm_length(const uint8_t *p, int n);
int bruhat_leq(const uint8_t *u, const uint8_t *w, int n);
int parse_perm(const char *s, uint8_t *p);          /* returns n or -1 */
void format_perm(const uint8_t *p, int n, char *out); /* n chars + NUL  */
void tables_build(tables_t *T, int n);
void tables_free(tables_t *T);
static inline const uint8_t *PERM(const tables_t *T, uint32_t v)
{
    return T->perm + (size_t)v * T->n;
}

/* sweep.c: the fast criterion.  For a fixed bottom u it visits every
 * w >= u with l(w) - l(u) <= maxd, level by level.  For each w it computes
 *   H[w][I] = min over chains C: u -> w of g_C(I),
 *   g_C(I)  = #{steps t_ab of C with [a, b-1] contained in I},
 *   F[w]    = 1 iff a single chain attains H[w] in every coordinate,
 * and F[w] = 1 iff SCNP holds for (u, w).                                 */
typedef struct {
    int32_t *slot;       /* nperm; -1 unless in the level being built      */
    uint32_t *node[2];
    vec_t *H[2];
    uint8_t *F[2];
} sweep_work_t;

typedef void (*level_fn)(void *ctx, uint32_t u, int d, const uint32_t *node,
                         const uint8_t *F, const vec_t *H, uint32_t cnt);

void sweep_work_init(sweep_work_t *W, const tables_t *T);
void sweep_work_free(sweep_work_t *W);
uint64_t sweep_from(const tables_t *T, sweep_work_t *W, uint32_t u, int maxd,
                    level_fn fn, void *ctx);

/* run.c: checkpointed production sweep */
int run_sweep(const char *dir, int n, int sym, int nchunks, int maxlen,
              int list, int gz);  /* list: 0 none, 1 orbit reps, 2 all */
int run_bottom(const uint8_t *u, int n, int maxlen, int quiet);
int run_summary(const char *dir);
int load_badw(const char *dir, uint8_t *bits, size_t nbytes);

/* search.c: the fast criterion from a list of bottoms (n <= 11), checkpointed */
int run_search(const char *dir, const char *bottoms, const char *patterns,
               const char *badw_prev, double hours, int chunk);
void commit_file(const char *tmp, const char *final);
void mkdir_p(const char *path);
int file_exists(const char *path);

/* bigbottom.c: the fast criterion from one bottom for n <= 16, no tables */
int run_bigbottom(const char *us, const char *prevfile, int maxlen);
/* bigsplit.c: the same with the levels in files under dir, on all threads */
int run_bigsplit(const char *us, const char *prevfile, int maxlen, const char *dir, double memgb);

/* patterns.c: compare the bad-w set of a finished sweep with patterns */
int run_patterns(const char *dir, const char *patterns);
int perm_contains(const uint8_t *w, int n, const uint8_t *p, int k);

/* greedy.c: are An-Tung-Zhang's greedy chains dominant? */
int run_greedy(int n, const char *patterns);

/* brute.c: explicit-expansion cross checks */
int run_check(int n, int maxlen, int snp);
int run_pair(const tables_t *T, const uint8_t *u, const uint8_t *w,
             int verbose, size_t limit);
int run_pairs_stdin(size_t limit);

#endif
