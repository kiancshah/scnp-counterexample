/*
 * brute.c - explicit-expansion cross checks of the fast criterion.
 *
 * Nothing here uses polytope theory.  For a bottom u we propagate up the
 * Bruhat graph, per node v,
 *   S(v) = supp D_{u,v}: the exponent vectors of all products obtained by
 *          choosing one variable y_k (a <= k < b) per step t_ab of a chain,
 *   M(v) = the distinct multisets of steps of the chains u -> v,
 * and declare SCNP(u,v) iff some multiset in M(v), expanded as a product of
 * linear forms, has exactly |S(v)| monomials.
 *
 * Monomials are packed 6 bits per variable (exponents <= 55, <= 10 vars).
 */
#include "scnp.h"

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    uint8_t c[VL]; /* multiplicity of each transposition t */
} ms_t;

typedef struct {
    uint64_t *a;
    size_t n, cap;
} u64v;

typedef struct {
    ms_t *a;
    size_t n, cap;
} msv;

static void *xrealloc(void *p, size_t sz)
{
    p = realloc(p, sz ? sz : 1);
    if (!p) {
        fprintf(stderr, "out of memory (%zu bytes)\n", sz);
        exit(1);
    }
    return p;
}

static inline void u64v_push(u64v *v, uint64_t x)
{
    if (v->n == v->cap) {
        v->cap = v->cap ? 2 * v->cap : 16;
        v->a = xrealloc(v->a, v->cap * sizeof *v->a);
    }
    v->a[v->n++] = x;
}

static inline void msv_push(msv *v, const ms_t *x)
{
    if (v->n == v->cap) {
        v->cap = v->cap ? 2 * v->cap : 16;
        v->a = xrealloc(v->a, v->cap * sizeof *v->a);
    }
    v->a[v->n++] = *x;
}

static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

static int cmp_ms(const void *a, const void *b)
{
    return memcmp(a, b, sizeof(ms_t));
}

static void u64v_uniq(u64v *v)
{
    if (v->n < 2)
        return;
    qsort(v->a, v->n, sizeof *v->a, cmp_u64);
    size_t k = 1;
    for (size_t i = 1; i < v->n; i++)
        if (v->a[i] != v->a[k - 1])
            v->a[k++] = v->a[i];
    v->n = k;
}

static void msv_uniq(msv *v)
{
    if (v->n < 2)
        return;
    qsort(v->a, v->n, sizeof *v->a, cmp_ms);
    size_t k = 1;
    for (size_t i = 1; i < v->n; i++)
        if (memcmp(&v->a[i], &v->a[k - 1], sizeof(ms_t)))
            v->a[k++] = v->a[i];
    v->n = k;
}

static void u64v_free(u64v *v)
{
    free(v->a);
    memset(v, 0, sizeof *v);
}

static void msv_free(msv *v)
{
    free(v->a);
    memset(v, 0, sizeof *v);
}

/* Support of prod_t (y_a + ... + y_{b-1})^{m[t]}, left in *A. */
static size_t expand(const tables_t *T, const ms_t *m, u64v *A, u64v *B)
{
    A->n = 0;
    u64v_push(A, 0);
    for (int t = 0; t < T->ntr; t++)
        for (int r = 0; r < m->c[t]; r++) {
            B->n = 0;
            for (size_t i = 0; i < A->n; i++)
                for (int k = T->ta[t]; k < T->tb[t]; k++)
                    u64v_push(B, A->a[i] + (1ull << (6 * k)));
            u64v_uniq(B);
            u64v tmp = *A;
            *A = *B;
            *B = tmp;
        }
    return A->n;
}

static void format_mono(uint64_t m, int nvar, char *out, size_t sz)
{
    size_t k = 0;
    out[0] = 0;
    for (int i = 0; i < nvar; i++) {
        int e = (int)(m >> (6 * i) & 63);
        if (!e)
            continue;
        k += (size_t)snprintf(out + k, sz - k, "%sy%d", k ? "*" : "", i + 1);
        if (e > 1)
            k += (size_t)snprintf(out + k, sz - k, "^%d", e);
    }
    if (!k)
        snprintf(out, sz, "1");
}

/* ------------------------------------------------------- brute sweep */

typedef struct {
    int32_t *slot;
    uint32_t *node[2];
    u64v *S[2];
    msv *M[2];
} bwork_t;

typedef void (*blevel_fn)(void *ctx, uint32_t u, int d, const uint32_t *node,
                          u64v *S, msv *M, uint32_t cnt);

static void bwork_init(bwork_t *W, const tables_t *T)
{
    W->slot = malloc(sizeof(int32_t) * T->nperm);
    memset(W->slot, 0xff, sizeof(int32_t) * T->nperm);
    for (int k = 0; k < 2; k++) {
        W->node[k] = malloc(sizeof(uint32_t) * T->maxlevel);
        W->S[k] = calloc(T->maxlevel, sizeof(u64v));
        W->M[k] = calloc(T->maxlevel, sizeof(msv));
    }
}

static void bwork_free(bwork_t *W)
{
    free(W->slot);
    for (int k = 0; k < 2; k++) {
        free(W->node[k]);
        free(W->S[k]);
        free(W->M[k]);
    }
}

/* Returns 0, or -1 if the sets exceeded `limit` elements in one level.
 * If wp != NULL only nodes x <= wp are visited (pair mode).             */
static int brute_sweep(const tables_t *T, bwork_t *W, uint32_t u, int maxd,
                       const uint8_t *wp, size_t limit, blevel_fn fn, void *ctx)
{
    int cur = 0;
    uint32_t ncur = 1;
    W->node[0][0] = u;
    u64v_push(&W->S[0][0], 0);
    ms_t zero;
    memset(&zero, 0, sizeof zero);
    msv_push(&W->M[0][0], &zero);
    fn(ctx, u, 0, W->node[0], W->S[0], W->M[0], 1);
    int rc = 0;

    for (int d = 0; d < maxd && ncur; d++) {
        int nx = cur ^ 1;
        uint32_t nn = 0;
        size_t raw = 0;
        for (uint32_t i = 0; i < ncur && rc == 0; i++) {
            uint32_t v = W->node[cur][i];
            u64v *Sv = &W->S[cur][i];
            msv *Mv = &W->M[cur][i];
            for (uint32_t e = T->off[v]; e < T->off[v + 1]; e++) {
                uint32_t x = T->cov[e];
                int t = T->ctr[e];
                if (wp && !bruhat_leq(PERM(T, x), wp, T->n))
                    continue;
                int32_t s = W->slot[x];
                if (s < 0) {
                    s = (int32_t)nn++;
                    W->slot[x] = s;
                    W->node[nx][s] = x;
                }
                u64v *Sx = &W->S[nx][s];
                msv *Mx = &W->M[nx][s];
                for (size_t j = 0; j < Sv->n; j++)
                    for (int k = T->ta[t]; k < T->tb[t]; k++)
                        u64v_push(Sx, Sv->a[j] + (1ull << (6 * k)));
                for (size_t j = 0; j < Mv->n; j++) {
                    ms_t m = Mv->a[j];
                    m.c[t]++;
                    msv_push(Mx, &m);
                }
                raw += Sv->n * (size_t)(T->tb[t] - T->ta[t]) + Mv->n;
            }
            if (raw > 4 * limit)
                rc = -1;
        }
        size_t tot = 0;
        for (uint32_t s = 0; s < nn; s++) {
            W->slot[W->node[nx][s]] = -1;
            if (rc == 0) {
                u64v_uniq(&W->S[nx][s]);
                msv_uniq(&W->M[nx][s]);
                tot += W->S[nx][s].n + W->M[nx][s].n;
            }
        }
        for (uint32_t i = 0; i < ncur; i++) {
            u64v_free(&W->S[cur][i]);
            msv_free(&W->M[cur][i]);
        }
        cur = nx;
        ncur = nn;
        if (rc == 0 && tot > limit)
            rc = -1;
        if (rc)
            break;
        fn(ctx, u, d + 1, W->node[cur], W->S[cur], W->M[cur], ncur);
    }
    for (uint32_t i = 0; i < ncur; i++) {
        u64v_free(&W->S[cur][i]);
        msv_free(&W->M[cur][i]);
    }
    return rc;
}

/* ------------------------------------------------- expansion cache */

typedef struct {
    ms_t *key;
    uint32_t *val;
    uint8_t *used;
    size_t cap, n;
} mscache_t;

static uint64_t ms_hash(const ms_t *m)
{
    uint64_t w[VL / 8], h = 0x9E3779B97F4A7C15ull;
    memcpy(w, m->c, VL);
    for (int i = 0; i < VL / 8; i++) {
        h ^= w[i];
        h *= 0xBF58476D1CE4E5B9ull;
        h ^= h >> 31;
    }
    return h;
}

static void cache_init(mscache_t *C, size_t cap)
{
    C->cap = cap;
    C->n = 0;
    C->key = malloc(cap * sizeof(ms_t));
    C->val = malloc(cap * sizeof(uint32_t));
    C->used = calloc(cap, 1);
}

static void cache_free(mscache_t *C)
{
    free(C->key);
    free(C->val);
    free(C->used);
}

static size_t cached_expand(const tables_t *T, mscache_t *C, const ms_t *m,
                            u64v *A, u64v *B)
{
    size_t i = ms_hash(m) & (C->cap - 1);
    while (C->used[i]) {
        if (!memcmp(&C->key[i], m, sizeof *m))
            return C->val[i];
        i = (i + 1) & (C->cap - 1);
    }
    size_t v = expand(T, m, A, B);
    if (C->n * 2 > C->cap) { /* full: start over */
        memset(C->used, 0, C->cap);
        C->n = 0;
        i = ms_hash(m) & (C->cap - 1);
    }
    C->used[i] = 1;
    C->key[i] = *m;
    C->val[i] = (uint32_t)v;
    C->n++;
    return v;
}

/* --------------------------------------------- SNP certificate (opt.) */

/* 1 if S = { a : |a| = d, a(X) >= min_{s in S} s(X) for all X }, which
 * certifies that S is the set of lattice points of a polytope.          */
static int snp_certified(const tables_t *T, const u64v *S, int d)
{
    int nv = T->n - 1, ns = 1 << nv;
    int z[1 << (MAXN - 1)];
    for (int X = 0; X < ns; X++)
        z[X] = 1 << 30;
    for (size_t i = 0; i < S->n; i++) {
        int sub[1 << (MAXN - 1)];
        sub[0] = 0;
        for (int X = 1; X < ns; X++) {
            int k = __builtin_ctz((unsigned)X);
            sub[X] = sub[X & (X - 1)] + (int)(S->a[i] >> (6 * k) & 63);
            if (sub[X] < z[X])
                z[X] = sub[X];
        }
    }
    /* enumerate compositions of d into nv parts */
    int a[MAXN], k = 0;
    size_t count = 0;
    a[0] = -1;
    int rem[MAXN + 1];
    rem[0] = d;
    while (k >= 0) {
        if (k == nv - 1) {
            a[k] = rem[k];
            int sub[1 << (MAXN - 1)], ok = 1;
            sub[0] = 0;
            for (int X = 1; X < ns && ok; X++) {
                int j = __builtin_ctz((unsigned)X);
                sub[X] = sub[X & (X - 1)] + a[j];
                ok = sub[X] >= z[X];
            }
            count += (size_t)ok;
            if (count > S->n)
                return 0;
            k--;
            continue;
        }
        if (++a[k] > rem[k]) {
            k--;
            continue;
        }
        rem[k + 1] = rem[k] - a[k];
        k++;
        a[k] = -1;
    }
    return count == S->n;
}

/* ------------------------------------------------------------ check */

typedef struct {
    const tables_t *T;
    uint8_t *fastF;
    uint64_t *fast_pairs;
    mscache_t cache;
    u64v A, B;
    int snp;
    uint64_t pairs[MAXL], fails[MAXL], mism[MAXL], fails_id[MAXL],
        nosnp[MAXL], fpairs[MAXL];
} chk_t;

static void chk_fast_level(void *vctx, uint32_t u, int d, const uint32_t *node,
                           const uint8_t *F, const vec_t *H, uint32_t cnt)
{
    (void)u;
    (void)H;
    chk_t *c = vctx;
    c->fpairs[d] += cnt;
    for (uint32_t s = 0; s < cnt; s++)
        c->fastF[node[s]] = F[s];
}

static void chk_brute_level(void *vctx, uint32_t u, int d, const uint32_t *node,
                            u64v *S, msv *M, uint32_t cnt)
{
    chk_t *c = vctx;
    const tables_t *T = c->T;
    for (uint32_t s = 0; s < cnt; s++) {
        size_t best = 0;
        for (size_t j = 0; j < M[s].n && best < S[s].n; j++) {
            size_t z = cached_expand(T, &c->cache, &M[s].a[j], &c->A, &c->B);
            if (z > best)
                best = z;
        }
        int ok = best == S[s].n;
        c->pairs[d]++;
        if (!ok) {
            c->fails[d]++;
            if (u == T->id)
                c->fails_id[d]++;
        }
        if (ok != c->fastF[node[s]]) {
            c->mism[d]++;
            char a[MAXN + 1], b[MAXN + 1];
            format_perm(PERM(T, u), T->n, a);
            format_perm(PERM(T, node[s]), T->n, b);
#pragma omp critical
            fprintf(stderr, "MISMATCH u=%s w=%s brute=%d fast=%d\n", a, b, ok,
                    c->fastF[node[s]]);
        }
        if (c->snp && !snp_certified(T, &S[s], d))
            c->nosnp[d]++;
    }
}

int run_check(int n, int maxlen, int snp)
{
    tables_t T;
    tables_build(&T, n);
    int N = T.N;
    if (maxlen < 0 || maxlen > N)
        maxlen = N;
    uint64_t P[MAXL] = {0}, Fl[MAXL] = {0}, Mi[MAXL] = {0}, Fi[MAXL] = {0},
             Ns[MAXL] = {0}, FP[MAXL] = {0};
    int aborted = 0;

#pragma omp parallel
    {
        chk_t c;
        memset(&c, 0, sizeof c);
        c.T = &T;
        c.snp = snp;
        c.fastF = calloc(T.nperm, 1);
        cache_init(&c.cache, 1u << 18);
        sweep_work_t SW;
        sweep_work_init(&SW, &T);
        bwork_t BW;
        bwork_init(&BW, &T);
#pragma omp for schedule(dynamic, 1)
        for (uint32_t u = 0; u < T.nperm; u++) {
            int top = N - T.len[u];
            if (top > maxlen)
                top = maxlen;
            sweep_from(&T, &SW, u, top, chk_fast_level, &c);
            if (brute_sweep(&T, &BW, u, top, NULL, (size_t)1 << 30,
                            chk_brute_level, &c)) {
#pragma omp atomic
                aborted++;
            }
        }
#pragma omp critical
        for (int d = 0; d <= N; d++) {
            P[d] += c.pairs[d];
            Fl[d] += c.fails[d];
            Mi[d] += c.mism[d];
            Fi[d] += c.fails_id[d];
            Ns[d] += c.nosnp[d];
            FP[d] += c.fpairs[d];
        }
        free(c.fastF);
        cache_free(&c.cache);
        u64v_free(&c.A);
        u64v_free(&c.B);
        sweep_work_free(&SW);
        bwork_free(&BW);
    }

    printf("# brute-force check of S_%d (all pairs u <= w%s)\n", n,
           maxlen < N ? ", restricted length" : "");
    printf("# %4s %12s %12s %12s %12s%s\n", "len", "pairs", "SCNP fails",
           "fails u=id", "mismatches", snp ? "  SNP-uncertified" : "");
    uint64_t tp = 0, tf = 0, tm = 0, ti = 0, tn = 0;
    int bad = aborted != 0;
    for (int d = 0; d <= maxlen; d++) {
        if (!P[d] && !FP[d])
            continue;
        if (P[d] != FP[d])
            bad = 1;
        printf("  %4d %12llu %12llu %12llu %12llu", d, (unsigned long long)P[d],
               (unsigned long long)Fl[d], (unsigned long long)Fi[d],
               (unsigned long long)Mi[d]);
        if (snp)
            printf("  %16llu", (unsigned long long)Ns[d]);
        printf("\n");
        tp += P[d];
        tf += Fl[d];
        tm += Mi[d];
        ti += Fi[d];
        tn += Ns[d];
    }
    printf("  %4s %12llu %12llu %12llu %12llu", "all", (unsigned long long)tp,
           (unsigned long long)tf, (unsigned long long)ti,
           (unsigned long long)tm);
    if (snp)
        printf("  %16llu", (unsigned long long)tn);
    printf("\n");
    if (aborted)
        printf("# WARNING: %d bottoms aborted (set size limit)\n", aborted);
    if (bad && !aborted)
        printf("# WARNING: fast and brute sweeps visited different pair sets\n");
    printf("# result: %s\n", (tm || bad) ? "FAILED" : "fast criterion agrees with brute force on every pair");
    tables_free(&T);
    return (tm || bad) ? 1 : 0;
}

/* ------------------------------------------------------------- pair */

typedef struct {
    uint32_t w;
    int d;
    int found;
    size_t nS, nM;
    u64v S;
    msv M;
} pair_capture_t;

static void pair_brute_level(void *vctx, uint32_t u, int d, const uint32_t *node,
                             u64v *S, msv *M, uint32_t cnt)
{
    (void)u;
    pair_capture_t *pc = vctx;
    if (d != pc->d)
        return;
    for (uint32_t s = 0; s < cnt; s++)
        if (node[s] == pc->w) {
            pc->found = 1;
            pc->S = S[s];
            pc->M = M[s];
            memset(&S[s], 0, sizeof S[s]); /* take ownership */
            memset(&M[s], 0, sizeof M[s]);
        }
}

typedef struct {
    uint32_t w;
    int d;
    int F;
} fast_capture_t;

static void pair_fast_level(void *vctx, uint32_t u, int d, const uint32_t *node,
                            const uint8_t *F, const vec_t *H, uint32_t cnt)
{
    (void)u;
    (void)H;
    fast_capture_t *fc = vctx;
    if (d != fc->d)
        return;
    for (uint32_t s = 0; s < cnt; s++)
        if (node[s] == fc->w)
            fc->F = F[s];
}

static void u128_str(unsigned __int128 x, char *out)
{
    char tmp[48];
    int k = 0;
    do {
        tmp[k++] = (char)('0' + (int)(x % 10));
        x /= 10;
    } while (x);
    for (int i = 0; i < k; i++)
        out[i] = tmp[k - 1 - i];
    out[k] = 0;
}

/* Chain count and interval size by a level-by-level DP inside [u, w]. */
static void count_chains(const tables_t *T, uint32_t u, const uint8_t *wp, int d,
                         char *chains, size_t *isize)
{
    uint32_t *lv[2];
    unsigned __int128 *cn[2];
    int32_t *slot = malloc(sizeof(int32_t) * T->nperm);
    memset(slot, 0xff, sizeof(int32_t) * T->nperm);
    for (int k = 0; k < 2; k++) {
        lv[k] = malloc(sizeof(uint32_t) * T->maxlevel);
        cn[k] = malloc(sizeof(unsigned __int128) * T->maxlevel);
    }
    int cur = 0, ovf = 0;
    uint32_t ncur = 1;
    long double *ap[2];
    ap[0] = malloc(sizeof(long double) * T->maxlevel);
    ap[1] = malloc(sizeof(long double) * T->maxlevel);
    lv[0][0] = u;
    cn[0][0] = 1;
    ap[0][0] = 1;
    *isize = 1;
    for (int l = 0; l < d; l++) {
        int nx = cur ^ 1;
        uint32_t nn = 0;
        for (uint32_t i = 0; i < ncur; i++) {
            uint32_t v = lv[cur][i];
            for (uint32_t e = T->off[v]; e < T->off[v + 1]; e++) {
                uint32_t x = T->cov[e];
                if (!bruhat_leq(PERM(T, x), wp, T->n))
                    continue;
                int32_t s = slot[x];
                if (s < 0) {
                    s = (int32_t)nn++;
                    slot[x] = s;
                    lv[nx][s] = x;
                    cn[nx][s] = 0;
                    ap[nx][s] = 0;
                }
                unsigned __int128 y = cn[nx][s] + cn[cur][i];
                if (y < cn[nx][s])
                    ovf = 1;
                cn[nx][s] = y;
                ap[nx][s] += ap[cur][i];
            }
        }
        for (uint32_t s = 0; s < nn; s++)
            slot[lv[nx][s]] = -1;
        *isize += nn;
        cur = nx;
        ncur = nn;
    }
    if (ovf)
        snprintf(chains, 64, "%.6Le", ap[cur][0]);
    else
        u128_str(cn[cur][0], chains);
    free(slot);
    for (int k = 0; k < 2; k++) {
        free(lv[k]);
        free(cn[k]);
        free(ap[k]);
    }
}

static void gvec_of(const tables_t *T, const ms_t *m, uint8_t g[VL])
{
    memset(g, 0, VL);
    for (int t = 0; t < T->ntr; t++)
        if (m->c[t])
            for (int I = 0; I < T->nint; I++)
                g[I] += (uint8_t)(m->c[t] * T->delta[t][I]);
}

int run_pair(const tables_t *T, const uint8_t *up, const uint8_t *wp,
             int verbose, size_t limit)
{
    int n = T->n;
    char us[MAXN + 1], ws[MAXN + 1];
    format_perm(up, n, us);
    format_perm(wp, n, ws);
    uint32_t u = perm_rank(up, n), w = perm_rank(wp, n);
    if (!bruhat_leq(up, wp, n)) {
        printf("%s %s: u is not <= w in Bruhat order\n", us, ws);
        return 1;
    }
    int d = T->len[w] - T->len[u];

    char chains[64];
    size_t isize;
    count_chains(T, u, wp, d, chains, &isize);

    fast_capture_t fc = {w, d, -1};
    sweep_work_t SW;
    sweep_work_init(&SW, T);
    sweep_from(T, &SW, u, d, pair_fast_level, &fc);
    sweep_work_free(&SW);

    pair_capture_t pc;
    memset(&pc, 0, sizeof pc);
    pc.w = w;
    pc.d = d;
    bwork_t BW;
    bwork_init(&BW, T);
    int rc = brute_sweep(T, &BW, u, d, wp, limit, pair_brute_level, &pc);
    bwork_free(&BW);

    long best = -1, nmax = 0, ndist = 0;
    int brute = -1, consistent = 1;
    ms_t bestm;
    u64v A = {0}, B = {0}, U = {0}, bestS = {0};
    if (rc == 0 && pc.found) {
        /* expand every distinct step multiset; distinct supports <-> g */
        uint8_t (*g)[VL] = malloc(pc.M.n * VL);
        for (size_t j = 0; j < pc.M.n; j++) {
            long z = (long)expand(T, &pc.M.a[j], &A, &B);
            for (size_t i = 0; i < A.n; i++)
                u64v_push(&U, A.a[i]);
            if (z > best) {
                best = z;
                bestm = pc.M.a[j];
                bestS.n = 0;
                for (size_t i = 0; i < A.n; i++)
                    u64v_push(&bestS, A.a[i]);
            }
            gvec_of(T, &pc.M.a[j], g[j]);
        }
        u64v_uniq(&U);
        consistent = U.n == pc.S.n && !memcmp(U.a, pc.S.a, U.n * sizeof *U.a);
        qsort(g, pc.M.n, VL, cmp_ms);
        for (size_t j = 0; j < pc.M.n; j++) {
            if (j && !memcmp(g[j], g[j - 1], VL))
                continue;
            ndist++;
            int minimal = 1; /* support not strictly inside another one */
            for (size_t k = 0; k < pc.M.n && minimal; k++) {
                int le = 1, ne = 0;
                for (int I = 0; I < T->nint; I++) {
                    le &= g[k][I] <= g[j][I];
                    ne |= g[k][I] != g[j][I];
                }
                if (le && ne)
                    minimal = 0;
            }
            nmax += minimal;
        }
        free(g);
        brute = best == (long)pc.S.n;
    }

    if (!verbose) {
        printf("%s %s %d chains=%s", us, ws, d, chains);
        if (brute >= 0)
            printf(" supp=%zu best_chain=%ld", pc.S.n, best);
        else
            printf(" supp=? best_chain=?");
        printf(" scnp=%s", fc.F ? "yes" : "NO");
        if (brute >= 0 && brute != fc.F)
            printf(" BRUTE_DISAGREES");
        if (!consistent)
            printf(" INCONSISTENT");
        printf("\n");
    } else {
        printf("pair             u = %s, w = %s   (n = %d)\n", us, ws, n);
        printf("lengths          l(u) = %d, l(w) = %d, interval length %d\n",
               T->len[u], T->len[w], d);
        printf("interval         |[u,w]| = %zu\n", isize);
        printf("chains           %s saturated chains u -> w\n", chains);
        if (brute < 0) {
            printf("brute force      skipped (set size limit %zu exceeded)\n", limit);
        } else {
            printf("supp D_{u,w}     %zu monomials in y1..y%d (degree %d)\n",
                   pc.S.n, n - 1, d);
            printf("chain supports   %zu distinct step multisets, %ld distinct "
                   "supports, %ld maximal\n", pc.M.n, ndist, nmax);
            printf("best chain       %ld monomials, steps", best);
            for (int t = 0; t < T->ntr; t++)
                for (int r = 0; r < bestm.c[t]; r++)
                    printf(" (%d%s%d)", T->ta[t] + 1, n > 9 ? "," : "", T->tb[t] + 1);
            printf("\n");
            if (best < (long)pc.S.n) {
                printf("not covered      ");
                size_t i = 0, k = 0, shown = 0;
                while (i < pc.S.n) {
                    if (k < bestS.n && bestS.a[k] == pc.S.a[i]) {
                        i++;
                        k++;
                        continue;
                    }
                    if (k < bestS.n && bestS.a[k] < pc.S.a[i]) {
                        k++;
                        continue;
                    }
                    if (shown < 12) {
                        char mono[128];
                        format_mono(pc.S.a[i], n - 1, mono, sizeof mono);
                        printf("%s%s", shown ? ", " : "", mono);
                    }
                    shown++;
                    i++;
                }
                if (shown > 12)
                    printf(", ... (%zu total)", shown);
                printf("\n");
            }
            if (verbose > 1) {
                printf("support          ");
                for (size_t i = 0; i < pc.S.n; i++) {
                    char mono[128];
                    format_mono(pc.S.a[i], n - 1, mono, sizeof mono);
                    printf("%s%s", i ? " + " : "", mono);
                }
                printf("\n");
            }
            if (!consistent)
                printf("WARNING          union of chain supports != DP support\n");
        }
        printf("SCNP             %s", fc.F ? "holds" : "FAILS");
        printf("   (fast criterion: %s", fc.F ? "holds" : "fails");
        if (brute >= 0)
            printf("; brute force: %s", brute ? "holds" : "fails");
        printf(")\n");
        if (brute >= 0 && brute != fc.F)
            printf("WARNING          fast criterion and brute force disagree\n");
    }
    u64v_free(&A);
    u64v_free(&B);
    u64v_free(&U);
    u64v_free(&bestS);
    u64v_free(&pc.S);
    msv_free(&pc.M);
    return (brute >= 0 && brute != fc.F) || !consistent ? 2 : 0;
}

/* Lines "u w [anything]" on stdin, one summary line per pair. */
int run_pairs_stdin(size_t limit)
{
    char line[1024];
    tables_t T;
    int curn = 0, rc = 0;
    while (fgets(line, sizeof line, stdin)) {
        char a[256], b[256];
        if (line[0] == '#' || sscanf(line, "%255s %255s", a, b) != 2)
            continue;
        uint8_t up[MAXN], wp[MAXN];
        int n1 = parse_perm(a, up), n2 = parse_perm(b, wp);
        if (n1 < 0 || n1 != n2) {
            fprintf(stderr, "bad line: %s", line);
            rc = 1;
            continue;
        }
        if (n1 != curn) {
            if (curn)
                tables_free(&T);
            tables_build(&T, n1);
            curn = n1;
        }
        rc |= run_pair(&T, up, wp, 0, limit);
        fflush(stdout);
    }
    if (curn)
        tables_free(&T);
    return rc;
}
