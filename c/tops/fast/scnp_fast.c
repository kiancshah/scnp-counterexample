/*
 * scnp_fast.c - fast SCNP checker for Postnikov-Stanley polynomials D_{u,w},
 * organised by TOP element w (n <= 11).
 *
 * Conventions: exactly those of ../reference/scnp_reference.py (and of the
 * main repository): steps swap POSITIONS a < b, label x_a + ... + x_{b-1}.
 *
 * Method: ../reference/scnp_criterion.py explains it in plain Python. For a
 * fixed top w, one pass down the lower interval [e, w], one length level at
 * a time, computes for every v <= w
 *     H(v)  = entrywise min over chains C from v to w of z_C   (z_C(R) =
 *             number of labels of C whose block lies inside block R),
 *     F(v)  = 1 iff some chain from v to w attains H(v),
 * and (v, w) is bad (SCNP fails) iff F(v) = 0. H is one byte per block;
 * n = 11 has 55 blocks, so H fits in one 64-byte AVX-512 register.
 *
 * Steps: by default a table of all steps into each permutation of S_n is
 * built once (n = 11: 568 356 480 steps, 2.8 GB, ~11 s); option -T instead
 * generates the steps of each node from its one-line notation (same result,
 * kept as an independent code path for cross-checks). The lexicographic rank
 * of each lower neighbour indexes a per-thread array "slot" that merges
 * duplicates within a level: 4 * n! bytes per thread (160 MB for n = 11).
 *
 * Commands (permutations 1-based, A = 10, B = 11):
 *   scnp_fast pair U W            fast test of one pair (with |[u,w]|)
 *   scnp_fast top W [-v]          all bad u below W
 *   scnp_fast list FILE OUT [-p]  every w in FILE (one per line): per-w line
 *                                 to OUT; with -p all bad pairs to OUT.pairs
 *   scnp_fast all N OUT           every w in S_N: per-length totals, and the
 *                                 bad-top bitmap OUT.bits (bit = lex rank)
 *   scnp_fast avoid N PATS OUT    every w in S_N avoiding all patterns in the
 *                                 comma-separated list PATS
 *   scnp_fast bitmap N LIST OUT   bitmap over ranks of S_N of the perms in LIST
 *                                 (with -b BITS, list also counts the pairs
 *                                 (v, w) with v in BITS, and the bad ones)
 *   scnp_fast contain FILE PATS OUT [-c]
 *                                 pattern-containment bitmask per line (or
 *                                 counts per mask with -c)
 *   scnp_fast batch JOBS list [-p] | batch JOBS witness PATS [maxtry]
 *                                 run list/witness on every "IN OUT" line of
 *                                 JOBS whose OUT does not exist yet (resume)
 *   scnp_fast witness FILE PATS OUT [maxtry]
 *                                 for each w in FILE, look for a bad pair
 *                                 (u, w) by planting a known bad pair of a
 *                                 pattern of w; PATS lists "p:u1/u2/..,p:..";
 *                                 falls back to a full pass when none works
 *   scnp_fast witness-all N PATS OUT [maxtry [RLO RHI]]
 *                                 the same for every w in S_N containing a
 *                                 pattern of PATS (summary + exceptions);
 *                                 RLO RHI: only lexicographic ranks in [RLO, RHI)
 *   options: -T (steps on the fly), -b BITS (class of bottoms, list mode)
 */
#include <omp.h>
#include <stdint.h>
#include <sys/mman.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef __AVX512BW__
#include <immintrin.h>
#endif

/* Build with -DMAXN=12 for n = 12 (66 blocks, two 64-byte halves per H). */
#ifndef MAXN
#define MAXN 11
#endif
#define NBLKMAX (MAXN * (MAXN - 1) / 2)
#define VL (NBLKMAX <= 64 ? 64 : 128)
#define MAXD (NBLKMAX + 2)

typedef uint8_t vec_t __attribute__((vector_size(VL), aligned(VL)));
typedef int8_t svec_t __attribute__((vector_size(VL), aligned(VL)));
typedef uint64_t pk_t; /* packed permutation: entry i (0-based value) in bits 4i..4i+3 */

static int n;                   /* permutations of {0..n-1} internally        */
static int nblk;                /* n(n-1)/2 blocks [i..j] of variables 0..n-2  */
static uint32_t fact[MAXN + 1];
static uint32_t nperm;
static vec_t DEL[MAXN][MAXN];   /* DEL[a][b][R] = 1 iff block [a..b-1] inside R */

static void *xmalloc(size_t sz)
{
    void *p = aligned_alloc(64, (sz + 63) / 64 * 64 + 64);
    if (!p) {
        fprintf(stderr, "out of memory (%zu bytes)\n", sz);
        exit(1);
    }
    return p;
}

/* Large arrays that are accessed at random: 2 MB aligned and marked for
 * transparent huge pages, which saves TLB misses (about 8% at n = 11).
 * (Tried and dropped, both slower at n = 11: a per-level hash table instead
 * of the direct slot array, and sorting each level by rank.) */
static void *xmalloc_huge(size_t sz)
{
    size_t al = (size_t)2 << 20;
    sz = (sz + al - 1) / al * al;
    void *p = aligned_alloc(al, sz);
    if (!p) {
        fprintf(stderr, "out of memory (%zu bytes)\n", sz);
        exit(1);
    }
    madvise(p, sz, MADV_HUGEPAGE);
    return p;
}

static void setup(int nn)
{
    if (nn < 2 || nn > MAXN) {
        fprintf(stderr, "n must be in 2..%d\n", MAXN);
        exit(1);
    }
    n = nn;
    nblk = n * (n - 1) / 2;
    fact[0] = 1;
    for (int i = 1; i <= MAXN; i++)
        fact[i] = fact[i - 1] * (uint32_t)i;
    nperm = fact[n];
    /* blocks [i..j], 0 <= i <= j <= n-2, numbered in a fixed order */
    int bi[VL], bj[VL], R = 0;
    for (int i = 0; i + 1 < n; i++)
        for (int j = i; j + 1 < n; j++) {
            bi[R] = i;
            bj[R] = j;
            R++;
        }
    for (int a = 0; a < n; a++)
        for (int b = a + 1; b < n; b++) {
            uint8_t d[VL] = {0};
            for (int r = 0; r < nblk; r++)
                d[r] = bi[r] <= a && b - 1 <= bj[r];
            memcpy(&DEL[a][b], d, VL);
        }
}

/* ---------- permutations ---------- */

static inline pk_t pack(const uint8_t *p)
{
    pk_t x = 0;
    for (int i = 0; i < n; i++)
        x |= (pk_t)p[i] << (4 * i);
    return x;
}

static inline void unpack(pk_t x, uint8_t *p)
{
    for (int i = 0; i < n && i < MAXN; i++)
        p[i] = (uint8_t)((x >> (4 * i)) & 15u);
}

static inline uint32_t rank_perm(const uint8_t *p)
{
    uint32_t r = 0;
    unsigned used = 0;
    for (int i = 0; i < n; i++) {
        unsigned c = p[i] - (unsigned)__builtin_popcount(used & ((1u << p[i]) - 1u));
        r = r * (uint32_t)(n - i) + c;
        used |= 1u << p[i];
    }
    return r;
}

static void unrank_perm(uint32_t r, uint8_t *p)
{
    unsigned c[MAXN];
    for (int i = n - 1; i >= 0; i--) {
        c[i] = r % (uint32_t)(n - i);
        r /= (uint32_t)(n - i);
    }
    unsigned avail = (1u << n) - 1u;
    for (int i = 0; i < n; i++) {
        unsigned m = avail;
        for (unsigned k = 0; k < c[i]; k++)
            m &= m - 1;
        int v = __builtin_ctz(m);
        p[i] = (uint8_t)v;
        avail &= ~(1u << v);
    }
}

static int perm_len(const uint8_t *p, int k)
{
    int l = 0;
    for (int i = 0; i < k; i++)
        for (int j = i + 1; j < k; j++)
            l += p[i] > p[j];
    return l;
}

/* Tableau criterion: u <= v iff for all k, j: #{i<=k: u(i)>=j} <= #{i<=k: v(i)>=j}. */
static int bruhat_leq(const uint8_t *u, const uint8_t *v)
{
    int cu[MAXN] = {0}, cv[MAXN] = {0};
    for (int k = 0; k < n; k++) {
        for (int j = 0; j <= u[k]; j++)
            cu[j]++;
        for (int j = 0; j <= v[k]; j++)
            cv[j]++;
        for (int j = 0; j < n; j++)
            if (cu[j] > cv[j])
                return 0;
    }
    return 1;
}

/* "214365", "54821A9376", 'A'/'0' = 10, 'B' = 11.  Returns length or -1. */
static int parse_perm(const char *s, uint8_t *p)
{
    int k = 0;
    unsigned seen = 0;
    for (; *s && *s != '\n' && *s != ' ' && *s != '\t' && *s != ','; s++) {
        int v;
        if (*s >= '1' && *s <= '9')
            v = *s - '0';
        else if (*s == 'A' || *s == 'a' || *s == '0')
            v = 10;
        else if (*s == 'B' || *s == 'b')
            v = 11;
        else if (*s == 'C' || *s == 'c')
            v = 12;
        else
            return -1;
        if (k >= MAXN)
            return -1;
        p[k++] = (uint8_t)(v - 1);
    }
    for (int i = 0; i < k; i++) {
        if (p[i] >= k || (seen >> p[i] & 1u))
            return -1;
        seen |= 1u << p[i];
    }
    return k;
}

static void fmt_perm(const uint8_t *p, int k, char *out)
{
    for (int i = 0; i < k; i++)
        out[i] = "123456789ABC"[p[i]];
    out[k] = 0;
}

/* ---------- patterns ---------- */

/* Does w (length nw) contain p (length k)?  Backtracking over positions. */
static int contains_rec(const uint8_t *w, int nw, const uint8_t *p, int k,
                        int t, int start, int *pos)
{
    if (t == k)
        return 1;
    for (int i = start; i <= nw - (k - t); i++) {
        int ok = 1;
        for (int s = 0; s < t && ok; s++)
            ok = (w[pos[s]] < w[i]) == (p[s] < p[t]);
        if (ok) {
            pos[t] = i;
            if (contains_rec(w, nw, p, k, t + 1, i + 1, pos))
                return 1;
        }
    }
    return 0;
}

static int contains(const uint8_t *w, int nw, const uint8_t *p, int k)
{
    int pos[MAXN];
    if (k > nw)
        return 0;
    return contains_rec(w, nw, p, k, 0, 0, pos);
}

typedef struct {
    int k;
    uint8_t p[MAXN];
    int nb;                 /* known bad bottoms of p (for witness mode) */
    uint8_t b[64][MAXN];
} pattern_t;

static int parse_patterns(const char *s, pattern_t *P, int maxp)
{
    int np = 0;
    char buf[4096];
    snprintf(buf, sizeof buf, "%s", s);
    for (char *tok = strtok(buf, ","); tok && np < maxp; tok = strtok(NULL, ",")) {
        pattern_t *q = &P[np];
        memset(q, 0, sizeof *q);
        char *colon = strchr(tok, ':');
        if (colon)
            *colon = 0;
        q->k = parse_perm(tok, q->p);
        if (q->k < 0) {
            fprintf(stderr, "bad pattern %s\n", tok);
            exit(1);
        }
        if (colon) {
            for (char *b = strtok_r(colon + 1, "/", &colon); b; b = strtok_r(NULL, "/", &colon)) {
                if (parse_perm(b, q->b[q->nb]) != q->k) {
                    fprintf(stderr, "bad bottom %s\n", b);
                    exit(1);
                }
                q->nb++;
            }
        }
        np++;
    }
    return np;
}

/* ---------- the sweep ---------- */

/* Fold candidate (cand, fc) into (h, f): h = entrywise min, f = OR of fc
 * over candidates equal to the current minimum. */
static inline void relax(vec_t *h, uint8_t *f, vec_t cand, uint8_t fc)
{
#if defined(__AVX512BW__) && VL == 64
    __m512i a = (__m512i)*h, c = (__m512i)cand;
    __mmask64 le = _mm512_cmple_epu8_mask(c, a);
    __mmask64 ge = _mm512_cmpge_epu8_mask(c, a);
    if (ge == ~0ULL) {          /* cand >= h */
        if (le == ~0ULL)        /* cand == h */
            *f |= fc;
        return;
    }
    if (le == ~0ULL) {          /* cand <= h, cand != h */
        *h = cand;
        *f = fc;
        return;
    }
    *h = (vec_t)_mm512_min_epu8(a, c); /* incomparable */
    *f = 0;
#elif defined(__AVX512BW__)
    /* VL = 128 (n = 12): the same test on two 64-byte halves */
    const __m512i *a = (const __m512i *)h, *c = (const __m512i *)&cand;
    int all_le = 1, all_ge = 1;
    for (int k = 0; k < VL / 64; k++) {
        all_le &= _mm512_cmple_epu8_mask(c[k], a[k]) == ~0ULL;
        all_ge &= _mm512_cmpge_epu8_mask(c[k], a[k]) == ~0ULL;
    }
    if (all_ge) {
        if (all_le)
            *f |= fc;
        return;
    }
    if (all_le) {
        *h = cand;
        *f = fc;
        return;
    }
    __m512i *hm = (__m512i *)h;
    for (int k = 0; k < VL / 64; k++)
        hm[k] = _mm512_min_epu8(a[k], c[k]);
    *f = 0;
#else
    int lo = 0, go = 0;
    for (int i = 0; i < VL; i++) {
        lo |= cand[i] < (*h)[i];
        go |= cand[i] > (*h)[i];
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
    for (int i = 0; i < VL; i++)
        if (cand[i] < (*h)[i])
            (*h)[i] = cand[i];
    *f = 0;
#endif
}

typedef struct {
    uint32_t *slot;         /* nperm entries: 1 + index of the node in the level
                             * being built, 0 if not in it                      */
    size_t slotbytes;
    pk_t *P[2];
    uint32_t *R[2];
    vec_t *H[2];
    uint8_t *F[2];
    size_t cap[2];
} work_t;

/* Optional class of bottoms (option -b FILE: bitmap over lexicographic
 * ranks, bit set = v in the class). Nodes in the class are counted apart. */
static uint8_t *CLS = NULL;

typedef struct {
    uint64_t nodes[MAXD];   /* nodes[d] = #{v <= w (>= lo) : l(w) - l(v) = d} */
    uint64_t bad[MAXD];     /* of those, how many give a bad pair (v, w)     */
    uint64_t cnodes[MAXD];  /* the same, restricted to v in the class CLS    */
    uint64_t cbad[MAXD];
    uint64_t edges;
    int record;             /* store bad bottoms in list                     */
    pk_t *list;
    uint8_t *listd;
    size_t nlist, caplist;
} top_t;

static void work_init(work_t *W)
{
    /* zero pages from the kernel: only the parts actually touched get memory */
    W->slotbytes = sizeof(uint32_t) * (size_t)nperm;
    W->slot = mmap(NULL, W->slotbytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (W->slot == MAP_FAILED) {
        fprintf(stderr, "out of memory (slot array)\n");
        exit(1);
    }
    madvise(W->slot, W->slotbytes, MADV_HUGEPAGE);
    for (int k = 0; k < 2; k++) {
        W->cap[k] = 1024;
        W->P[k] = xmalloc(sizeof(pk_t) * W->cap[k]);
        W->R[k] = xmalloc(sizeof(uint32_t) * W->cap[k]);
        W->H[k] = xmalloc(sizeof(vec_t) * W->cap[k]);
        W->F[k] = xmalloc(W->cap[k]);
    }
}

static void work_free(work_t *W)
{
    munmap(W->slot, W->slotbytes);
    for (int k = 0; k < 2; k++) {
        free(W->P[k]);
        free(W->R[k]);
        free(W->H[k]);
        free(W->F[k]);
    }
}

static void grow(work_t *W, int k)
{
    size_t c = W->cap[k] * 2;
    pk_t *P = xmalloc(sizeof(pk_t) * c);
    uint32_t *R = xmalloc(sizeof(uint32_t) * c);
    vec_t *H = xmalloc(sizeof(vec_t) * c);
    uint8_t *F = xmalloc(c);
    memcpy(P, W->P[k], sizeof(pk_t) * W->cap[k]);
    memcpy(R, W->R[k], sizeof(uint32_t) * W->cap[k]);
    memcpy(H, W->H[k], sizeof(vec_t) * W->cap[k]);
    memcpy(F, W->F[k], W->cap[k]);
    free(W->P[k]);
    free(W->R[k]);
    free(W->H[k]);
    free(W->F[k]);
    W->P[k] = P;
    W->R[k] = R;
    W->H[k] = H;
    W->F[k] = F;
    W->cap[k] = c;
}

static void top_record(top_t *T, pk_t v, int d)
{
    if (T->nlist == T->caplist) {
        T->caplist = T->caplist ? 2 * T->caplist : 256;
        T->list = realloc(T->list, sizeof(pk_t) * T->caplist);
        T->listd = realloc(T->listd, T->caplist);
        if (!T->list || !T->listd) {
            fprintf(stderr, "out of memory (bad list)\n");
            exit(1);
        }
    }
    T->list[T->nlist] = v;
    T->listd[T->nlist] = (uint8_t)d;
    T->nlist++;
}

/* ---------- optional global table of steps into each permutation ----------
 * DOFF[r] .. DOFF[r+1] index the steps v -> x into the permutation x of rank
 * r: DCOV[e] = rank of v, DTR[e] = number of the transposition (a, b).
 * Built once per run (n = 11: 568 356 480 steps, 2.8 GB, ~11 s); without it the
 * steps are generated from the one-line notation (the same steps; the
 * second code path is kept for cross-checks, option -T).                  */
static uint32_t *DOFF;
static uint32_t *DCOV;
static uint8_t *DTR;
static uint8_t TA[NBLKMAX], TB[NBLKMAX];
static vec_t DELT[NBLKMAX];
static int use_table = 1;

static void build_tables(void)
{
    int t = 0;
    for (int a = 0; a < n; a++)
        for (int b = a + 1; b < n; b++) {
            TA[t] = (uint8_t)a;
            TB[t] = (uint8_t)b;
            DELT[t] = DEL[a][b];
            t++;
        }
    int tid[MAXN][MAXN];
    t = 0;
    for (int a = 0; a < n; a++)
        for (int b = a + 1; b < n; b++)
            tid[a][b] = t++;
    DOFF = xmalloc_huge(sizeof(uint32_t) * ((size_t)nperm + 1));
#pragma omp parallel for schedule(static)
    for (uint32_t r = 0; r < nperm; r++) {
        uint8_t x[MAXN];
        unrank_perm(r, x);
        uint32_t c = 0;
        for (int a = 0; a + 1 < n; a++) {
            int lim = -1;
            for (int b = a + 1; b < n; b++)
                if (x[b] < x[a] && x[b] > lim) {
                    lim = x[b];
                    c++;
                }
        }
        DOFF[r + 1] = c;
    }
    DOFF[0] = 0;
    for (uint32_t r = 0; r < nperm; r++)
        DOFF[r + 1] += DOFF[r];
    size_t E = DOFF[nperm];
    DCOV = xmalloc_huge(sizeof(uint32_t) * E);
    DTR = xmalloc_huge(E);
#pragma omp parallel for schedule(static)
    for (uint32_t r = 0; r < nperm; r++) {
        uint8_t x[MAXN];
        unrank_perm(r, x);
        uint32_t e = DOFF[r];
        for (int a = 0; a + 1 < n; a++) {
            int lim = -1;
            for (int b = a + 1; b < n; b++)
                if (x[b] < x[a] && x[b] > lim) {
                    lim = x[b];
                    uint8_t v[MAXN];
                    memcpy(v, x, (size_t)n);
                    v[a] = x[b];
                    v[b] = x[a];
                    DCOV[e] = rank_perm(v);
                    DTR[e] = (uint8_t)tid[a][b];
                    e++;
                }
        }
    }
}

/* setup(k) and, unless -T, the step table; repeated calls are free. */
static void prepare(int k)
{
    static int done_n = 0;
    if (done_n == k)
        return;
    if (done_n) {
        fprintf(stderr, "one n per run\n");
        exit(1);
    }
    setup(k);
    if (k >= 12)
        use_table = 0;  /* the S_12 table would need about 36 GB */
    if (use_table) {
        double t0 = omp_get_wtime();
        build_tables();
        fprintf(stderr, "step table for S_%d: %u steps, %.1f s\n", k, DOFF[nperm],
                omp_get_wtime() - t0);
    }
    done_n = k;
}

/* Add the candidate (cand, fc) for the lower neighbour v (rank r, packed pv)
 * to the level being built. */
static inline void push(work_t *W, int nx, size_t *nn, uint32_t r, pk_t pv,
                        vec_t cand, uint8_t fc)
{
    uint32_t s = W->slot[r];
    if (s == 0) {
        if (*nn == W->cap[nx])
            grow(W, nx);
        s = (uint32_t)(*nn)++;
        W->slot[r] = s + 1;
        W->P[nx][s] = pv;
        W->R[nx][s] = r;
        W->H[nx][s] = cand;
        W->F[nx][s] = fc;
    } else {
        relax(&W->H[nx][s - 1], &W->F[nx][s - 1], cand, fc);
    }
}

static inline pk_t swap_nibbles(pk_t px, int a, int b)
{
    pk_t dx = ((px >> (4 * a)) ^ (px >> (4 * b))) & 15u;
    return px ^ (dx << (4 * a)) ^ (dx << (4 * b));
}

/* One pass down from w over all v <= w (and v >= lo if lo != NULL).
 * Fills T (counts per depth; bad bottoms if T->record).                   */
static void sweep_down(work_t *W, const uint8_t *w, const uint8_t *lo, top_t *T)
{
    int cur = 0;
    size_t ncur = 1;
    memset(T->nodes, 0, sizeof T->nodes);
    memset(T->bad, 0, sizeof T->bad);
    memset(T->cnodes, 0, sizeof T->cnodes);
    memset(T->cbad, 0, sizeof T->cbad);
    T->edges = 0;
    T->nlist = 0;
    W->P[0][0] = pack(w);
    W->R[0][0] = rank_perm(w);
    memset(&W->H[0][0], 0, VL);
    W->F[0][0] = 1;
    T->nodes[0] = 1;
    if (CLS && (CLS[W->R[0][0] >> 3] >> (W->R[0][0] & 7) & 1))
        T->cnodes[0] = 1;
    uint32_t *slot = W->slot;
    int tab = use_table && DOFF != NULL;

    for (int d = 0; ncur; d++) {
        int nx = cur ^ 1;
        size_t nn = 0;
        uint64_t edges = 0;
        for (size_t i = 0; i < ncur; i++) {
            pk_t px = W->P[cur][i];
            uint32_t rx = W->R[cur][i];
            vec_t hx = W->H[cur][i];
            uint8_t fx = W->F[cur][i];
            if (tab) {
                uint32_t e1 = DOFF[rx + 1];
                for (uint32_t e = DOFF[rx]; e < e1; e++) {
                    int t = DTR[e];
                    pk_t pv = swap_nibbles(px, TA[t], TB[t]);
                    if (lo) {
                        uint8_t v[MAXN];
                        unpack(pv, v);
                        if (!bruhat_leq(lo, v))
                            continue;
                    }
                    edges++;
                    push(W, nx, &nn, DCOV[e], pv, hx + DELT[t], fx);
                }
                continue;
            }
            uint8_t x[MAXN], v[MAXN] = {0};
            unpack(px, x);
            for (int a = 0; a + 1 < n; a++) {
                int xa = x[a], lim = -1;
                for (int b = a + 1; b < n; b++) {
                    int xb = x[b];
                    if (xb >= xa || xb < lim)
                        continue;
                    /* v = x t_ab is a step v -> x (no c between with xb < x(c) < xa) */
                    lim = xb;
                    if (lo) {
                        memcpy(v, x, (size_t)n);
                        v[a] = (uint8_t)xb;
                        v[b] = (uint8_t)xa;
                        if (!bruhat_leq(lo, v))
                            continue;
                    }
                    edges++;
                    /* Lehmer-code rank of v from that of x: only the digits at
                     * a and b change, by -(1+k) and +k, where k counts the
                     * positions j > b with xb < x(j) < xa (none lie between). */
                    int k = 0;
                    for (int j = b + 1; j < n; j++)
                        k += x[j] > xb && x[j] < xa;
                    uint32_t r = rx - (uint32_t)(1 + k) * fact[n - 1 - a] +
                                 (uint32_t)k * fact[n - 1 - b];
                    pk_t dx = (pk_t)(xa ^ xb);
                    push(W, nx, &nn, r, px ^ (dx << (4 * a)) ^ (dx << (4 * b)),
                         hx + DEL[a][b], fx);
                }
            }
        }
        T->edges += edges;
        uint64_t nb = 0;
        for (size_t s = 0; s < nn; s++) {
            slot[W->R[nx][s]] = 0;
            if (CLS && (CLS[W->R[nx][s] >> 3] >> (W->R[nx][s] & 7) & 1)) {
                T->cnodes[d + 1]++;
                T->cbad[d + 1] += !W->F[nx][s];
            }
            if (!W->F[nx][s]) {
                nb++;
                if (T->record)
                    top_record(T, W->P[nx][s], d + 1);
            }
        }
        if (nn) {
            T->nodes[d + 1] = nn;
            T->bad[d + 1] = nb;
        }
        cur = nx;
        ncur = nn;
    }
}

static uint64_t sum64(const uint64_t *a)
{
    uint64_t s = 0;
    for (int d = 0; d < MAXD; d++)
        s += a[d];
    return s;
}

/* Is (u, w) bad?  Restricted pass over [u, w].  *size = |[u, w]|. */
static int pair_bad(work_t *W, const uint8_t *u, const uint8_t *w, uint64_t *size)
{
    top_t T;
    memset(&T, 0, sizeof T);
    sweep_down(W, w, u, &T);
    if (size)
        *size = sum64(T.nodes);
    int d = perm_len(w, n) - perm_len(u, n);
    return T.nodes[d] == 1 && T.bad[d] == 1;
}

/* ---------- commands ---------- */

static int cmd_pair(const char *us, const char *ws)
{
    uint8_t u[MAXN], w[MAXN];
    int k = parse_perm(us, u);
    if (k < 0 || parse_perm(ws, w) != k) {
        fprintf(stderr, "bad permutations\n");
        return 1;
    }
    setup(k);
    if (!bruhat_leq(u, w)) {
        printf("%s %s not_leq\n", us, ws);
        return 0;
    }
    work_t W;
    work_init(&W);
    uint64_t size;
    int bad = pair_bad(&W, u, w, &size);
    printf("%s %s len=%d interval=%llu scnp=%s\n", us, ws,
           perm_len(w, n) - perm_len(u, n), (unsigned long long)size,
           bad ? "NO" : "yes");
    work_free(&W);
    return 0;
}

static int cmd_top(const char *ws, int verbose)
{
    uint8_t w[MAXN];
    int k = parse_perm(ws, w);
    if (k < 0)
        return 1;
    setup(k);
    work_t W;
    work_init(&W);
    top_t T;
    memset(&T, 0, sizeof T);
    T.record = 1;
    double t0 = omp_get_wtime();
    sweep_down(&W, w, NULL, &T);
    double t1 = omp_get_wtime();
    printf("w = %s: |[e,w]| = %llu, edges = %llu, bad bottoms = %zu (%.3f s)\n",
           ws, (unsigned long long)sum64(T.nodes), (unsigned long long)T.edges,
           T.nlist, t1 - t0);
    if (verbose)
        for (size_t i = 0; i < T.nlist; i++) {
            uint8_t u[MAXN];
            char s[16];
            unpack(T.list[i], u);
            fmt_perm(u, n, s);
            printf("  %s %d\n", s, T.listd[i]);
        }
    free(T.list);
    free(T.listd);
    work_free(&W);
    return 0;
}

/* Read permutations (first token of each line) from a file. */
static uint8_t *read_perm_file(const char *path, size_t *count, int *k)
{
    FILE *f = strcmp(path, "-") ? fopen(path, "r") : stdin;
    if (!f) {
        perror(path);
        exit(1);
    }
    size_t cap = 1024, cnt = 0;
    uint8_t *A = malloc(cap * MAXN);
    char line[256];
    *k = -1;
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#' || line[0] == '\n')
            continue;
        uint8_t p[MAXN];
        int kk = parse_perm(line, p);
        if (kk < 0 || (*k >= 0 && kk != *k)) {
            fprintf(stderr, "bad line: %s", line);
            exit(1);
        }
        *k = kk;
        if (cnt == cap) {
            cap *= 2;
            A = realloc(A, cap * MAXN);
        }
        memcpy(A + cnt * MAXN, p, MAXN);
        cnt++;
    }
    if (f != stdin)
        fclose(f);
    *count = cnt;
    return A;
}

/* list FILE OUT [-p] */
static int cmd_list(const char *in, const char *out, int pairs)
{
    size_t cnt;
    int k;
    uint8_t *A = read_perm_file(in, &cnt, &k);
    if (cnt == 0) {
        FILE *fo = fopen(out, "w");
        fclose(fo);
        return 0;
    }
    prepare(k);
    char **lines = calloc(cnt, sizeof(char *));
    top_t *res = calloc(cnt, sizeof(top_t));
    uint64_t tot_nodes[MAXD] = {0}, tot_bad[MAXD] = {0}, tot_edges = 0;
    uint64_t tot_cn[MAXD] = {0}, tot_cb[MAXD] = {0};
    double t0 = omp_get_wtime();
#pragma omp parallel
    {
        work_t W;
        work_init(&W);
        uint64_t ln[MAXD] = {0}, lb[MAXD] = {0}, le = 0, lcn[MAXD] = {0}, lcb[MAXD] = {0};
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < cnt; i++) {
            const uint8_t *w = A + i * MAXN;
            top_t *T = &res[i];
            T->record = 1;
            sweep_down(&W, w, NULL, T);
            int lw = perm_len(w, n);
            int dmin = -1, dmax = -1;
            char buf[1024], ws[16];
            int pos = 0;
            uint64_t cn = 0, cb = 0;
            for (int d = 0; d < MAXD; d++) {
                ln[d] += T->nodes[d];
                lb[d] += T->bad[d];
                lcn[d] += T->cnodes[d];
                lcb[d] += T->cbad[d];
                cn += T->cnodes[d];
                cb += T->cbad[d];
                if (T->bad[d]) {
                    if (dmin < 0)
                        dmin = d;
                    dmax = d;
                }
            }
            le += T->edges;
            fmt_perm(w, n, ws);
            pos += snprintf(buf + pos, sizeof buf - pos, "%s\t%d\t%llu\t%llu\t%zu\t%d\t%d\t%llu\t%llu\t",
                            ws, lw, (unsigned long long)sum64(T->nodes),
                            (unsigned long long)T->edges, T->nlist, dmin, dmax,
                            (unsigned long long)cn, (unsigned long long)cb);
            int first = 1;
            for (int d = 0; d < MAXD; d++)
                if (T->bad[d]) {
                    pos += snprintf(buf + pos, sizeof buf - pos, "%s%d:%llu", first ? "" : ",",
                                    d, (unsigned long long)T->bad[d]);
                    first = 0;
                }
            if (first)
                pos += snprintf(buf + pos, sizeof buf - pos, "-");
            lines[i] = strdup(buf);
            if (!pairs) {
                free(T->list);
                free(T->listd);
                T->list = NULL;
                T->listd = NULL;
            }
        }
#pragma omp critical
        {
            for (int d = 0; d < MAXD; d++) {
                tot_nodes[d] += ln[d];
                tot_bad[d] += lb[d];
                tot_cn[d] += lcn[d];
                tot_cb[d] += lcb[d];
            }
            tot_edges += le;
        }
        work_free(&W);
    }
    double t1 = omp_get_wtime();
    char tmp[4096];
    snprintf(tmp, sizeof tmp, "%s.tmp", out);
    FILE *fo = fopen(tmp, "w");
    fprintf(fo, "# w\tl(w)\t|[e,w]|\tedges\tnbad\tminbadlen\tmaxbadlen\tclass_nodes\tclass_bad\tbad_by_len\n");
    for (size_t i = 0; i < cnt; i++) {
        fprintf(fo, "%s\n", lines[i]);
        free(lines[i]);
    }
    fprintf(fo, "# totals: #T len pairs bad class_pairs class_bad\n");
    for (int d = 0; d < MAXD; d++)
        if (tot_nodes[d])
            fprintf(fo, "#T %d %llu %llu %llu %llu\n", d, (unsigned long long)tot_nodes[d],
                    (unsigned long long)tot_bad[d], (unsigned long long)tot_cn[d],
                    (unsigned long long)tot_cb[d]);
    fprintf(fo, "#E %llu edges %.2f s %d threads\n", (unsigned long long)tot_edges, t1 - t0,
            omp_get_max_threads());
    fclose(fo);
    if (pairs) {
        char pp[4096], ptmp[4096];
        snprintf(pp, sizeof pp, "%s.pairs", out);
        snprintf(ptmp, sizeof ptmp, "%s.pairs.tmp", out);
        FILE *fp = fopen(ptmp, "w");
        for (size_t i = 0; i < cnt; i++) {
            char ws[16], us[16];
            fmt_perm(A + i * MAXN, n, ws);
            for (size_t j = 0; j < res[i].nlist; j++) {
                uint8_t u[MAXN];
                unpack(res[i].list[j], u);
                fmt_perm(u, n, us);
                fprintf(fp, "%s %s %d\n", us, ws, res[i].listd[j]);
            }
            free(res[i].list);
            free(res[i].listd);
        }
        fclose(fp);
        rename(ptmp, pp);
    }
    rename(tmp, out);
    free(res);
    free(lines);
    free(A);
    return 0;
}

/* all N OUT: every w in S_N */
static int cmd_all(int k, const char *out)
{
    prepare(k);
    size_t nbytes = (nperm + 7) / 8;
    uint8_t *bits = calloc(nbytes, 1);
    uint64_t tot_nodes[MAXD] = {0}, tot_bad[MAXD] = {0}, tot_edges = 0;
    double t0 = omp_get_wtime();
#pragma omp parallel
    {
        work_t W;
        work_init(&W);
        top_t T;
        memset(&T, 0, sizeof T);
        uint64_t ln[MAXD] = {0}, lb[MAXD] = {0}, le = 0;
#pragma omp for schedule(dynamic, 64)
        for (uint32_t r = 0; r < nperm; r++) {
            uint8_t w[MAXN];
            unrank_perm(r, w);
            sweep_down(&W, w, NULL, &T);
            int isbad = 0;
            for (int d = 0; d < MAXD; d++) {
                ln[d] += T.nodes[d];
                lb[d] += T.bad[d];
                isbad |= T.bad[d] != 0;
            }
            le += T.edges;
            if (isbad) {
#pragma omp atomic
                bits[r >> 3] |= (uint8_t)(1u << (r & 7));
            }
        }
#pragma omp critical
        {
            for (int d = 0; d < MAXD; d++) {
                tot_nodes[d] += ln[d];
                tot_bad[d] += lb[d];
            }
            tot_edges += le;
        }
        work_free(&W);
    }
    double t1 = omp_get_wtime();
    uint64_t np = 0, nb = 0, ntop = 0;
    for (uint32_t r = 0; r < nperm; r++)
        ntop += bits[r >> 3] >> (r & 7) & 1;
    char path[4096];
    snprintf(path, sizeof path, "%s.bits", out);
    FILE *fb = fopen(path, "wb");
    fwrite(bits, 1, nbytes, fb);
    fclose(fb);
    FILE *fo = fopen(out, "w");
    fprintf(fo, "# S_%d, every top w, scnp_fast all\n#  len pairs bad\n", n);
    for (int d = 0; d < MAXD; d++)
        if (tot_nodes[d]) {
            fprintf(fo, "%4d %14llu %14llu\n", d, (unsigned long long)tot_nodes[d],
                    (unsigned long long)tot_bad[d]);
            np += tot_nodes[d];
            nb += tot_bad[d];
        }
    fprintf(fo, "all %llu %llu\nbad_tops %llu\nedges %llu\nseconds %.2f\nthreads %d\n",
            (unsigned long long)np, (unsigned long long)nb, (unsigned long long)ntop,
            (unsigned long long)tot_edges, t1 - t0, omp_get_max_threads());
    fclose(fo);
    free(bits);
    return 0;
}

/* avoid N PATS OUT: build avoiders level by level by inserting the largest
 * value (deleting the largest entry of an avoider gives an avoider). */
static int cmd_avoid(int k, const char *pats, const char *out)
{
    pattern_t *P = malloc(sizeof(pattern_t) * 64);
    int np = parse_patterns(pats, P, 64);
    uint8_t *cur = malloc(MAXN);
    size_t ncur = 1;
    cur[0] = 0;
    for (int m = 2; m <= k; m++) {
        size_t cap = ncur * (size_t)m;
        uint8_t *nxt = malloc(cap * MAXN);
        uint8_t *keep = calloc(cap, 1);
#pragma omp parallel for schedule(dynamic, 256)
        for (size_t i = 0; i < ncur; i++) {
            for (int pos = 0; pos < m; pos++) {
                uint8_t *q = nxt + (i * m + pos) * MAXN;
                int t = 0;
                for (int j = 0; j < m; j++)
                    q[j] = (uint8_t)(j == pos ? m - 1 : cur[i * MAXN + t++]);
                int ok = 1;
                for (int s = 0; s < np && ok; s++)
                    ok = !contains(q, m, P[s].p, P[s].k);
                keep[i * m + pos] = (uint8_t)ok;
            }
        }
        size_t nn = 0;
        for (size_t i = 0; i < cap; i++)
            if (keep[i])
                memmove(nxt + nn++ * MAXN, nxt + i * MAXN, MAXN);
        free(cur);
        free(keep);
        cur = nxt;
        ncur = nn;
        fprintf(stderr, "n=%d: %zu avoiders\n", m, ncur);
    }
    FILE *fo = fopen(out, "w");
    for (size_t i = 0; i < ncur; i++) {
        char s[16];
        fmt_perm(cur + i * MAXN, k, s);
        fprintf(fo, "%s\n", s);
    }
    fclose(fo);
    free(cur);
    free(P);
    return 0;
}

/* witness FILE PATS OUT */
typedef struct {
    int dlen;
    uint8_t u[MAXN];
} cand_t;

typedef struct {
    const uint8_t *w;
    const pattern_t *P;
    cand_t *c;
    int nc, cap;
} occ_ctx;

static void occ_rec(occ_ctx *C, int t, int start, int *pos)
{
    const pattern_t *P = C->P;
    int k = P->k;
    if (t == k) {
        uint8_t vals[MAXN];
        for (int s = 0; s < k; s++)
            vals[s] = C->w[pos[s]];
        /* sort values ascending */
        for (int s = 1; s < k; s++)
            for (int r = s; r > 0 && vals[r - 1] > vals[r]; r--) {
                uint8_t tmp = vals[r];
                vals[r] = vals[r - 1];
                vals[r - 1] = tmp;
            }
        for (int bb = 0; bb < P->nb; bb++) {
            if (C->nc == C->cap) {
                C->cap = C->cap ? 2 * C->cap : 64;
                C->c = realloc(C->c, sizeof(cand_t) * C->cap);
            }
            cand_t *q = &C->c[C->nc++];
            memcpy(q->u, C->w, (size_t)n);
            for (int s = 0; s < k; s++)
                q->u[pos[s]] = vals[P->b[bb][s]];
            q->dlen = perm_len(C->w, n) - perm_len(q->u, n);
        }
        return;
    }
    for (int i = start; i <= n - (k - t); i++) {
        int ok = 1;
        for (int s = 0; s < t && ok; s++)
            ok = (C->w[pos[s]] < C->w[i]) == (P->p[s] < P->p[t]);
        if (ok) {
            pos[t] = i;
            occ_rec(C, t + 1, i + 1, pos);
        }
    }
}

static int cand_cmp(const void *a, const void *b)
{
    return ((const cand_t *)a)->dlen - ((const cand_t *)b)->dlen;
}

/* Plant known bad pairs: for each occurrence of a pattern p (in the order
 * given) at positions i_1 < ... < i_k of w, and each listed bad bottom b of
 * p, let u agree with w off those positions and carry the values of w there
 * in the relative order of b. Then u <= w. Candidates are tried by
 * increasing l(w) - l(u), at most maxtry in all; each is decided exactly by
 * the restricted pass over [u, w]. Returns 1 and the first bad pair found. */
static int find_witness(work_t *W, occ_ctx *C, const uint8_t *w, const pattern_t *P,
                        int np, int maxtry, cand_t *best, int *pat_used, uint64_t *size,
                        int *tries)
{
    *tries = 0;
    for (int s = 0; s < np; s++) {
        C->w = w;
        C->P = &P[s];
        C->nc = 0;
        int pos[MAXN];
        occ_rec(C, 0, 0, pos);
        qsort(C->c, (size_t)C->nc, sizeof(cand_t), cand_cmp);
        for (int c = 0; c < C->nc && *tries < maxtry; c++) {
            if (c > 0 && memcmp(C->c[c].u, C->c[c - 1].u, (size_t)n) == 0)
                continue;
            (*tries)++;
            if (pair_bad(W, C->c[c].u, w, size)) {
                *best = C->c[c];
                *pat_used = s;
                return 1;
            }
        }
    }
    return 0;
}

static int cmd_witness(const char *in, const char *pats, const char *out, int maxtry)
{
    size_t cnt;
    int k;
    uint8_t *A = read_perm_file(in, &cnt, &k);
    pattern_t *P = malloc(sizeof(pattern_t) * 64);
    int np = parse_patterns(pats, P, 64);
    if (cnt == 0) {
        FILE *fo = fopen(out, "w");
        fclose(fo);
        return 0;
    }
    prepare(k);
    char **lines = calloc(cnt, sizeof(char *));
    uint64_t nfull = 0, nnone = 0, ntried = 0;
    double t0 = omp_get_wtime();
#pragma omp parallel
    {
        work_t W;
        work_init(&W);
        occ_ctx C;
        memset(&C, 0, sizeof C);
#pragma omp for schedule(dynamic, 16) reduction(+ : nfull, nnone, ntried)
        for (size_t i = 0; i < cnt; i++) {
            const uint8_t *w = A + i * MAXN;
            char ws[16], us[16], buf[256];
            fmt_perm(w, n, ws);
            int tries = 0, pat_used = -1;
            cand_t best;
            uint64_t size = 0;
            int found = find_witness(&W, &C, w, P, np, maxtry, &best, &pat_used, &size, &tries);
            ntried += (uint64_t)tries;
            if (found) {
                char ps[16];
                fmt_perm(best.u, n, us);
                fmt_perm(P[pat_used].p, P[pat_used].k, ps);
                snprintf(buf, sizeof buf, "%s\t%s\t%d\t%llu\tplant:%s\t%d", ws, us, best.dlen,
                         (unsigned long long)size, ps, tries);
            } else {
                /* fall back to the full pass down from w */
                top_t T;
                memset(&T, 0, sizeof T);
                T.record = 1;
                sweep_down(&W, w, NULL, &T);
                nfull++;
                if (T.nlist) {
                    /* report a bad bottom of least depth */
                    size_t jb = 0;
                    for (size_t j = 1; j < T.nlist; j++)
                        if (T.listd[j] < T.listd[jb])
                            jb = j;
                    uint8_t u[MAXN];
                    unpack(T.list[jb], u);
                    fmt_perm(u, n, us);
                    uint64_t sz;
                    pair_bad(&W, u, w, &sz);
                    snprintf(buf, sizeof buf, "%s\t%s\t%d\t%llu\tfull:%zu\t%d", ws, us,
                             T.listd[jb], (unsigned long long)sz, T.nlist, tries);
                } else {
                    nnone++;
                    snprintf(buf, sizeof buf, "%s\tNONE\t-1\t0\tfull:0\t%d", ws, tries);
                }
                free(T.list);
                free(T.listd);
            }
            lines[i] = strdup(buf);
        }
        free(C.c);
        work_free(&W);
    }
    double t1 = omp_get_wtime();
    char tmp[4096];
    snprintf(tmp, sizeof tmp, "%s.tmp", out);
    FILE *fo = fopen(tmp, "w");
    fprintf(fo, "# w\twitness_u\tlen\t|[u,w]|\thow\ttries\n");
    for (size_t i = 0; i < cnt; i++) {
        fprintf(fo, "%s\n", lines[i]);
        free(lines[i]);
    }
    fprintf(fo, "#S %zu tops, %llu full passes, %llu with no bad bottom, %llu planted pairs tried, %.2f s\n",
            cnt, (unsigned long long)nfull, (unsigned long long)nnone,
            (unsigned long long)ntried, t1 - t0);
    fclose(fo);
    rename(tmp, out);
    free(lines);
    free(A);
    free(P);
    return 0;
}

/* witness-all N PATS OUT [maxtry]: every w in S_N that contains at least
 * one pattern of PATS gets a witness (a bad pair (u, w)) by planting, or else
 * a full pass down from w. Writes OUT (summary), OUT.none (w containing a
 * pattern with no bad bottom at all; expected empty), OUT.full (w that needed
 * a full pass) and OUT.sample (every 1000th witness, for re-checking by the
 * Python reference). */
static int cmd_witness_all(int k, const char *pats, const char *out, int maxtry,
                           uint32_t rlo, uint32_t rhi)
{
    pattern_t *P = malloc(sizeof(pattern_t) * 64);
    int np = parse_patterns(pats, P, 64);
    prepare(k);
    uint64_t hist[64][MAXD];          /* [pattern used][l(w) - l(u)] */
    uint64_t ncont = 0, nfull = 0, nnone = 0, ntried = 0, nfullbad = 0;
    memset(hist, 0, sizeof hist);
    char path[4096];
    snprintf(path, sizeof path, "%s.none", out);
    FILE *fnone = fopen(path, "w");
    snprintf(path, sizeof path, "%s.full", out);
    FILE *ffull = fopen(path, "w");
    snprintf(path, sizeof path, "%s.sample", out);
    FILE *fsamp = fopen(path, "w");
    double t0 = omp_get_wtime();
#pragma omp parallel
    {
        work_t W;
        work_init(&W);
        occ_ctx C;
        memset(&C, 0, sizeof C);
        uint64_t lh[64][MAXD];
        memset(lh, 0, sizeof lh);
#pragma omp for schedule(dynamic, 1024) reduction(+ : ncont, nfull, nnone, ntried, nfullbad)
        for (uint32_t r = rlo; r < (rhi < nperm ? rhi : nperm); r++) {
            uint8_t w[MAXN];
            unrank_perm(r, w);
            int any = 0;
            for (int s = 0; s < np && !any; s++)
                any = contains(w, n, P[s].p, P[s].k);
            if (!any)
                continue;
            ncont++;
            int tries = 0, pat_used = -1;
            cand_t best;
            uint64_t size = 0;
            char ws[16], us[16];
            fmt_perm(w, n, ws);
            int found = find_witness(&W, &C, w, P, np, maxtry, &best, &pat_used, &size, &tries);
            ntried += (uint64_t)tries;
            if (found) {
                lh[pat_used][best.dlen]++;
                if (r % 1000 == 0) {
                    fmt_perm(best.u, n, us);
#pragma omp critical(samp)
                    fprintf(fsamp, "%s %s %d\n", us, ws, best.dlen);
                }
                continue;
            }
            top_t T;
            memset(&T, 0, sizeof T);
            T.record = 1;
            sweep_down(&W, w, NULL, &T);
            nfull++;
            if (T.nlist) {
                nfullbad++;
                size_t jb = 0;
                for (size_t j = 1; j < T.nlist; j++)
                    if (T.listd[j] < T.listd[jb])
                        jb = j;
                uint8_t u[MAXN];
                unpack(T.list[jb], u);
                fmt_perm(u, n, us);
#pragma omp critical(full)
                fprintf(ffull, "%s %s %d %zu\n", us, ws, T.listd[jb], T.nlist);
            } else {
                nnone++;
#pragma omp critical(none)
                fprintf(fnone, "%s\n", ws);
            }
            free(T.list);
            free(T.listd);
        }
#pragma omp critical
        for (int s = 0; s < 64; s++)
            for (int d = 0; d < MAXD; d++)
                hist[s][d] += lh[s][d];
        free(C.c);
        work_free(&W);
    }
    double t1 = omp_get_wtime();
    fclose(fnone);
    fclose(ffull);
    fclose(fsamp);
    char tmpout[4096];
    snprintf(tmpout, sizeof tmpout, "%s.tmp", out);
    FILE *fo = fopen(tmpout, "w");
    fprintf(fo, "# witness-all S_%d, patterns %s, maxtry %d, ranks [%u, %u)\n", n, pats, maxtry,
            rlo, rhi < nperm ? rhi : nperm);
    fprintf(fo, "containing a pattern %llu\nplanted witness found %llu\nfull passes %llu\n"
                "full pass found a bad bottom %llu\nno bad bottom at all %llu\n"
                "planted pairs tried %llu\nseconds %.1f\n",
            (unsigned long long)ncont, (unsigned long long)(ncont - nfull),
            (unsigned long long)nfull, (unsigned long long)nfullbad, (unsigned long long)nnone,
            (unsigned long long)ntried, t1 - t0);
    fprintf(fo, "# witnesses by pattern planted and l(w)-l(u)\n");
    for (int s = 0; s < np; s++) {
        char ps[16];
        fmt_perm(P[s].p, P[s].k, ps);
        for (int d = 0; d < MAXD; d++)
            if (hist[s][d])
                fprintf(fo, "H %s %d %llu\n", ps, d, (unsigned long long)hist[s][d]);
    }
    fclose(fo);
    rename(tmpout, out); /* OUT exists only once the range is finished */
    free(P);
    return 0;
}

/* contain FILE PATS OUT: for the first permutation on each line of FILE,
 * the bitmask of the patterns of PATS it contains (bit s = pattern s).
 * OUT gets "line<TAB>mask" for every line, unless -c: then only the number
 * of lines per mask. */
static int cmd_contain(const char *in, const char *pats, const char *out, int counts_only)
{
    pattern_t *P = malloc(sizeof(pattern_t) * 64);
    int np = parse_patterns(pats, P, 64);
    if (counts_only && np > 20) {
        fprintf(stderr, "contain -c: at most 20 patterns\n");
        return 1;
    }
    FILE *f = strcmp(in, "-") ? fopen(in, "r") : stdin;
    if (!f) {
        perror(in);
        return 1;
    }
    FILE *fo = fopen(out, "w");
    size_t cap = 1 << 16;
    char **L = malloc(sizeof(char *) * cap);
    uint64_t *M = malloc(sizeof(uint64_t) * cap);
    uint64_t *hist = calloc((size_t)1 << (np < 20 ? np : 20), sizeof(uint64_t));
    char line[512];
    for (;;) {
        size_t cnt = 0;
        while (cnt < cap && fgets(line, sizeof line, f)) {
            if (line[0] == '#')
                continue;
            line[strcspn(line, "\n")] = 0;
            L[cnt++] = strdup(line);
        }
        if (!cnt)
            break;
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < cnt; i++) {
            uint8_t w[MAXN];
            int k = parse_perm(L[i], w);
            uint64_t m = 0;
            for (int s = 0; s < np; s++)
                if (k > 0 && contains(w, k, P[s].p, P[s].k))
                    m |= 1ull << s;
            M[i] = m;
        }
        for (size_t i = 0; i < cnt; i++) {
            if (counts_only)
                hist[M[i]]++;
            else
                fprintf(fo, "%s\t%llu\n", L[i], (unsigned long long)M[i]);
            free(L[i]);
        }
    }
    if (counts_only)
        for (size_t m = 0; m < ((size_t)1 << (np < 20 ? np : 20)); m++)
            if (hist[m])
                fprintf(fo, "%zu %llu\n", m, (unsigned long long)hist[m]);
    fclose(fo);
    if (f != stdin)
        fclose(f);
    free(L);
    free(M);
    free(hist);
    free(P);
    return 0;
}

/* bitmap N LIST OUT: bitmap over ranks of S_N, bit set for each perm of LIST */
static int cmd_bitmap(int k, const char *in, const char *out)
{
    setup(k);
    size_t cnt;
    int kk;
    uint8_t *A = read_perm_file(in, &cnt, &kk);
    if (cnt && kk != k) {
        fprintf(stderr, "length mismatch\n");
        return 1;
    }
    size_t nbytes = (nperm + 7) / 8;
    uint8_t *bits = calloc(nbytes, 1);
    for (size_t i = 0; i < cnt; i++) {
        uint32_t r = rank_perm(A + i * MAXN);
        bits[r >> 3] |= (uint8_t)(1u << (r & 7));
    }
    FILE *fo = fopen(out, "wb");
    fwrite(bits, 1, nbytes, fo);
    fclose(fo);
    free(bits);
    free(A);
    return 0;
}

static void load_class(const char *path)
{
    FILE *f = fopen(path, "rb");
    if (!f) {
        perror(path);
        exit(1);
    }
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    CLS = malloc((size_t)sz);
    if (fread(CLS, 1, (size_t)sz, f) != (size_t)sz) {
        fprintf(stderr, "short read %s\n", path);
        exit(1);
    }
    fclose(f);
}

/* batch JOBS list [-p] | batch JOBS witness PATS [maxtry]:
 * JOBS has lines "IN OUT"; jobs whose OUT exists are skipped (resume). */
static int cmd_batch(const char *jobs, int argc, char **argv)
{
    FILE *f = fopen(jobs, "r");
    if (!f) {
        perror(jobs);
        return 1;
    }
    char in[2048], out[2048];
    while (fscanf(f, "%2047s %2047s", in, out) == 2) {
        FILE *g = fopen(out, "r");
        if (g) {
            fclose(g);
            continue;
        }
        double t0 = omp_get_wtime();
        if (!strcmp(argv[0], "list"))
            cmd_list(in, out, argc > 1 && !strcmp(argv[1], "-p"));
        else if (!strcmp(argv[0], "witness") && argc >= 2)
            cmd_witness(in, argv[1], out, argc >= 3 ? atoi(argv[2]) : 200);
        else {
            fprintf(stderr, "batch: unknown mode\n");
            return 1;
        }
        fprintf(stderr, "%s done in %.1f s\n", out, omp_get_wtime() - t0);
    }
    fclose(f);
    return 0;
}

int main(int argc, char **argv)
{
    /* option -T anywhere: generate steps on the fly instead of using the table */
    for (int i = 1; i < argc; i++)
        if (!strcmp(argv[i], "-T")) {
            use_table = 0;
            for (int j = i; j + 1 < argc; j++)
                argv[j] = argv[j + 1];
            argc--;
            i--;
        }
    /* option -b FILE anywhere: class bitmap for list mode */
    for (int i = 1; i + 1 < argc; i++)
        if (!strcmp(argv[i], "-b")) {
            load_class(argv[i + 1]);
            for (int j = i; j + 2 < argc; j++)
                argv[j] = argv[j + 2];
            argc -= 2;
            i--;
        }
    if (argc == 5 && !strcmp(argv[1], "bitmap"))
        return cmd_bitmap(atoi(argv[2]), argv[3], argv[4]);
    if ((argc == 5 || argc == 6 || argc == 8) && !strcmp(argv[1], "witness-all"))
        return cmd_witness_all(atoi(argv[2]), argv[3], argv[4], argc >= 6 ? atoi(argv[5]) : 200,
                               argc == 8 ? (uint32_t)strtoul(argv[6], NULL, 10) : 0,
                               argc == 8 ? (uint32_t)strtoul(argv[7], NULL, 10) : UINT32_MAX);
    if ((argc == 5 || argc == 6) && !strcmp(argv[1], "contain"))
        return cmd_contain(argv[2], argv[3], argv[4], argc == 6 && !strcmp(argv[5], "-c"));
    if (argc >= 4 && !strcmp(argv[1], "batch"))
        return cmd_batch(argv[2], argc - 3, argv + 3);
    if (argc == 4 && !strcmp(argv[1], "pair"))
        return cmd_pair(argv[2], argv[3]);
    if ((argc == 3 || argc == 4) && !strcmp(argv[1], "top"))
        return cmd_top(argv[2], argc == 4);
    if ((argc == 4 || argc == 5) && !strcmp(argv[1], "list"))
        return cmd_list(argv[2], argv[3], argc == 5 && !strcmp(argv[4], "-p"));
    if (argc == 4 && !strcmp(argv[1], "all"))
        return cmd_all(atoi(argv[2]), argv[3]);
    if (argc == 5 && !strcmp(argv[1], "avoid"))
        return cmd_avoid(atoi(argv[2]), argv[3], argv[4]);
    if ((argc == 5 || argc == 6) && !strcmp(argv[1], "witness"))
        return cmd_witness(argv[2], argv[3], argv[4], argc == 6 ? atoi(argv[5]) : 200);
    fprintf(stderr, "usage: see the comment at the top of scnp_fast.c\n");
    return 1;
}
