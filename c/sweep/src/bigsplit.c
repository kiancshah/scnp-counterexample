/*
 * bigsplit.c - bigbottom for intervals that do not fit in memory, on all threads.
 *
 *   scnp bigbottom U --tmp DIR [--mem GB] [--prev BADW_S10.bin.gz] [--maxlen L]
 *
 * Same recursion and same output as bigbottom.c (n <= 16).  The differences:
 *  - a level is a file DIR/level_<d>.bin of records (key, H, F), read through mmap;
 *  - the next level is built in R rounds, round r taking the covers whose hash
 *    falls in part r of R, so that one round fits in --mem GB (default 10); R is
 *    chosen from the exact level sizes, counted beforehand by a DP over prefix
 *    sets, and the number of nodes found at each level is checked against them;
 *  - in a round the threads share one hash table: a CAS on the key claims a
 *    slot, and a spin lock per record guards the relaxation;
 *  - each finished level adds a line to DIR/state.txt (and its avoiding pairs to
 *    DIR/avoid.txt); started again with the same DIR, the run resumes after the
 *    last finished level;
 *  - w avoids the patterns iff it does not contain one, and w contains one iff a
 *    one-letter deletion does, down to S_11 (a bitmap made from the S_10 one).
 */
#include "scnp.h"

#include <fcntl.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

#include "bigbottom.h"

static int n, N;        /* permutations of {0..n-1}, N = n(n-1)/2 lanes */
static size_t RS;       /* record: key (8 bytes), H (N), F, lock; rounded up to 8 */
static bvec_t delta[BN][BN];
static uint8_t *bad10;  /* S_10: bad = contains a pattern */
static uint8_t *cont11; /* S_11: contains a pattern */

/* cnt[k] = number of w >= U with l(w) = l(U) + k.  w >= U iff for every i the
 * set S of the first i letters of w has, for every j, at least as many letters
 * > j as the first i letters of U; DP over S and the inversions so far.   */
static void count_levels(const uint8_t *u, uint64_t *cnt)
{
    size_t ns = (size_t)1 << n;
    uint64_t *f = calloc(ns * (size_t)(N + 1), sizeof(uint64_t));
    uint8_t *reach = calloc(ns, 1);
    int need[BN][BN];
    if (!f || !reach) {
        fprintf(stderr, "out of memory (count_levels)\n");
        exit(1);
    }
    for (int i = 0; i <= n; i++)
        for (int j = 0; j < n; j++) {
            need[i][j] = 0;
            for (int k = 0; k < i; k++)
                need[i][j] += u[k] > j;
        }
    f[0] = 1;
    reach[0] = 1;
    for (size_t S = 0; S < ns; S++) {   /* every subset of S comes before S */
        if (!reach[S])
            continue;
        int i = __builtin_popcountll(S);
        for (int v = 0; v < n; v++) {
            if (S >> v & 1)
                continue;
            size_t T = S | (size_t)1 << v;
            int ok = 1;
            for (int j = 0; j < n && ok; j++)
                ok = __builtin_popcountll(T >> (j + 1)) >= need[i + 1][j];
            if (!ok)
                continue;
            int add = __builtin_popcountll(S >> (v + 1));
            uint64_t *fs = f + S * (size_t)(N + 1), *ft = f + T * (size_t)(N + 1);
            for (int k = 0; k + add <= N; k++)
                ft[k + add] += fs[k];
            reach[T] = 1;
        }
    }
    int lu = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            lu += u[i] > u[j];
    const uint64_t *full = f + (ns - 1) * (size_t)(N + 1);
    for (int k = 0; k <= N; k++)
        cnt[k] = k + lu <= N ? full[k + lu] : 0;
    free(f);
    free(reach);
}

/* q = p (in S_k) with position i deleted and the values standardised */
static void delete_letter(const uint8_t *p, int k, int i, uint8_t *q)
{
    int m = 0;
    for (int j = 0; j < k; j++)
        if (j != i)
            q[m++] = (uint8_t)(p[j] - (p[j] > p[i]));
}

static int contains(const uint8_t *w, int m)
{
    if (m == 10) {
        uint32_t r = perm_rank(w, 10);
        return bad10[r >> 3] >> (r & 7) & 1;
    }
    if (m == 11 && cont11) {
        uint32_t r = perm_rank(w, 11);
        return cont11[r >> 3] >> (r & 7) & 1;
    }
    uint8_t q[BN];
    for (int i = 0; i < m; i++) {
        delete_letter(w, m, i, q);
        if (contains(q, m - 1))
            return 1;
    }
    return 0;
}

static void build_cont11(void)
{
    size_t nb = 39916800 / 8;
    uint8_t *c = calloc(nb, 1);
    if (!c) {
        fprintf(stderr, "out of memory (S_11 bitmap)\n");
        exit(1);
    }
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < nb; i++)     /* one byte per iteration: no shared writes */
        for (int t = 0; t < 8; t++) {
            uint8_t p[11];
            perm_unrank((uint32_t)(8 * i + t), 11, p);
            c[i] |= (uint8_t)(contains(p, 11) << t);
        }
    cont11 = c;
}

static inline bvec_t loadH(const uint8_t *r)
{
    bvec_t h = {0};
    memcpy(&h, r + 8, (size_t)N);
    return h;
}

/* the next level, or the part of it built in one round */
typedef struct {
    uint64_t key;
    uint32_t idx1;      /* record + 1; 0 until the record is set up */
    uint32_t pad;
} islot_t;

typedef struct {
    uint8_t *rec;       /* partcap records */
    islot_t *tab;       /* cap slots, linear probing */
    size_t partcap, cap, cnt;
} part_t;

static size_t find_or_add(part_t *P, uint64_t key, uint64_t m)
{
    size_t h = (size_t)(((unsigned __int128)m * P->cap) >> 64);
    for (;;) {
        islot_t *s = &P->tab[h];
        uint64_t k = __atomic_load_n(&s->key, __ATOMIC_ACQUIRE);
        if (!k) {
            if (__atomic_compare_exchange_n(&s->key, &k, key, 0, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE)) {
                size_t i = __atomic_fetch_add(&P->cnt, 1, __ATOMIC_RELAXED);
                if (i >= P->partcap) {
                    fprintf(stderr, "bigbottom: more nodes in a round than allocated (%zu)\n", P->partcap);
                    exit(1);
                }
                uint8_t *r = P->rec + i * RS;
                memset(r, 0, RS);
                memcpy(r, &key, 8);
                memset(r + 8, 255, (size_t)N);   /* H = +infinity, F = 0 */
                __atomic_store_n(&s->idx1, (uint32_t)(i + 1), __ATOMIC_RELEASE);
                return i;
            }
            /* another thread claimed the slot first: k is its key */
        }
        if (k == key) {
            uint32_t i1;
            while (!(i1 = __atomic_load_n(&s->idx1, __ATOMIC_ACQUIRE)))
                ;
            return i1 - 1;
        }
        if (++h == P->cap)
            h = 0;
    }
}

static void relax_rec(uint8_t *r, bvec_t cand, uint8_t fc)
{
    uint8_t *lock = r + 9 + N;
    while (__atomic_exchange_n(lock, 1, __ATOMIC_ACQUIRE))
        ;
    bvec_t h = loadH(r);
    brelax(&h, r + 8 + N, cand, fc);
    memcpy(r + 8, &h, (size_t)N);
    __atomic_store_n(lock, 0, __ATOMIC_RELEASE);
}

/* round r of R: every cover of every node of cur whose hash is in part r */
static void build_part(const uint8_t *cur, size_t ncur, part_t *P, uint32_t r, uint32_t R)
{
    int sh[BN];
    for (int a = 0; a < n; a++)
        sh[a] = 4 * (n - 1 - a);
#pragma omp parallel for schedule(dynamic, 1024)
    for (size_t i = 0; i < ncur; i++) {
        const uint8_t *rc = cur + i * RS;
        uint64_t key;
        memcpy(&key, rc, 8);
        uint8_t p[BN];
        unpack(key, n, p);
        bvec_t hi = loadH(rc);
        uint8_t fi = rc[8 + N];
        uint64_t ck[BV], cm[BV];
        uint8_t ca[BV], cb[BV];
        size_t ri[BV];
        int nc = 0;
        for (int a = 0; a < n - 1; a++) {
            unsigned lim = (unsigned)n;
            for (int b = a + 1; b < n; b++)
                if (p[b] > p[a] && p[b] < lim) {
                    lim = p[b];
                    uint64_t x = (uint64_t)(p[a] ^ p[b]);
                    uint64_t k2 = key ^ x << sh[a] ^ x << sh[b];
                    uint64_t m = mix(k2);
                    if ((uint32_t)(((m & 0xffffffffu) * R) >> 32) != r)
                        continue;
                    ck[nc] = k2;
                    cm[nc] = m;
                    ca[nc] = (uint8_t)a;
                    cb[nc] = (uint8_t)b;
                    __builtin_prefetch(&P->tab[(size_t)(((unsigned __int128)m * P->cap) >> 64)]);
                    nc++;
                }
        }
        for (int c = 0; c < nc; c++) {
            ri[c] = find_or_add(P, ck[c], cm[c]);
            __builtin_prefetch(P->rec + ri[c] * RS);
            __builtin_prefetch(P->rec + ri[c] * RS + 64);
        }
        for (int c = 0; c < nc; c++)
            relax_rec(P->rec + ri[c] * RS, hi + delta[ca[c]][cb[c]], fi);
    }
}

static void fmt_big(uint64_t key, char *s)
{
    uint8_t p[BN];
    unpack(key, n, p);
    for (int k = 0; k < n; k++)
        s[k] = BDIGITS[p[k]];
    s[n] = 0;
}

static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : x > y;
}

static void level_path(char *buf, size_t sz, const char *dir, int d, const char *ext)
{
    snprintf(buf, sz, "%s/level_%d.%s", dir, d, ext);
}

static void write_all(FILE *f, const void *p, size_t sz, const char *what)
{
    if (sz && fwrite(p, 1, sz, f) != sz) {
        fprintf(stderr, "bigbottom: cannot write %s\n", what);
        exit(1);
    }
}

static void close_synced(FILE *f, const char *what)
{
    if (fflush(f) || fsync(fileno(f)) || fclose(f)) {
        fprintf(stderr, "bigbottom: cannot write %s\n", what);
        exit(1);
    }
}

int run_bigsplit(const char *us, const char *prevfile, int maxlen, const char *dir, double memgb)
{
    uint8_t u[BN];
    n = parse_big(us, u);
    if (n < 0 || n * (n - 1) / 2 > BV) {
        fprintf(stderr, "bad bottom %s\n", us);
        return 1;
    }
    N = n * (n - 1) / 2;
    RS = (size_t)(10 + N + 7) / 8 * 8;
    if (prevfile && n >= 10) {
        size_t nb = 3628800 / 8;
        bad10 = malloc(nb);
        gzFile g = gzopen(prevfile, "rb");
        if (!g || gzread(g, bad10, (unsigned)nb) != (int)nb) {
            fprintf(stderr, "%s: cannot read the S_10 bitmap\n", prevfile);
            return 1;
        }
        gzclose(g);
        if (n >= 12)
            build_cont11();
    }
    double t0 = omp_get_wtime();

    /* interval coordinates and the step vectors delta[a][b] */
    int ii[BV], ij[BV], nint = 0;
    for (int i = 0; i + 1 < n; i++)
        for (int j = i; j + 1 < n; j++) {
            ii[nint] = i;
            ij[nint] = j;
            nint++;
        }
    for (int a = 0; a < n; a++)
        for (int b = a + 1; b < n; b++) {
            uint8_t dv[BV] = {0};
            for (int I = 0; I < nint; I++)
                dv[I] = ii[I] <= a && b - 1 <= ij[I];
            memcpy(&delta[a][b], dv, BV);
        }
    int lu = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            lu += u[i] > u[j];
    if (maxlen < 0 || maxlen > N - lu)
        maxlen = N - lu;
    uint64_t lev[BV + 1];
    count_levels(u, lev);
    size_t budget = (size_t)(memgb * (double)(1ull << 30));

    /* start, or resume after the last finished level */
    char path[4096], path2[4096], line[512];
    mkdir_p(dir);
    snprintf(path, sizeof path, "%s/state.txt", dir);
    int D = -1;
    uint64_t tp = 0, tf = 0, ta = 0;
    fprintf(stderr, "# bigbottom U = %s (n = %d, l(U) = %d): pairs (U, w) with w >= U\n# %4s %14s %14s %12s\n",
            us, n, lu, "len", "pairs", "SCNP fails", "avoiding");
    FILE *st = fopen(path, "r");
    if (st) {
        if (!fgets(line, sizeof line, st) || strncmp(line, "# U ", 4) || strncmp(line + 4, us, strlen(us)) ||
            line[4 + strlen(us)] != '\n') {
            fprintf(stderr, "%s is for another bottom\n", path);
            return 1;
        }
        int d;
        unsigned long long np, nf, na;
        while (fgets(line, sizeof line, st))
            if (sscanf(line, "%d %llu %llu %llu", &d, &np, &nf, &na) == 4) {
                if (d != D + 1 || np != lev[d]) {
                    fprintf(stderr, "%s: unexpected line %s", path, line);
                    return 1;
                }
                fprintf(stderr, "  %4d %14llu %14llu %12llu\n", d, np, nf, na);
                D = d;
                tp += np;
                tf += nf;
                ta += na;
            }
        fclose(st);
    }
    if (D < 0) {
        level_path(path2, sizeof path2, dir, 0, "bin");
        FILE *f = fopen(path2, "wb");
        uint8_t *r = calloc(1, RS);
        uint64_t key = pack(u, n);
        if (!f || !r) {
            fprintf(stderr, "bigbottom: cannot create %s\n", path2);
            return 1;
        }
        memcpy(r, &key, 8);
        r[8 + N] = 1;                   /* H = 0, F = 1 */
        write_all(f, r, RS, path2);
        close_synced(f, path2);
        free(r);
        snprintf(path2, sizeof path2, "%s/avoid.txt", dir);
        FILE *a = fopen(path2, "w");
        st = fopen(path, "w");
        if (!a || !st) {
            fprintf(stderr, "bigbottom: cannot create the files in %s\n", dir);
            return 1;
        }
        fclose(a);
        fprintf(st, "# U %s\n0 1 0 0 1 0.0\n", us);
        close_synced(st, path);
        fprintf(stderr, "  %4d %14d %14d %12d\n", 0, 1, 0, 0);
        D = 0;
        tp = 1;
    }
    else {                              /* the pairs found so far, without a level cut short */
        snprintf(path2, sizeof path2, "%s/avoid.txt", dir);
        FILE *a = fopen(path2, "r");
        char *keep = NULL;
        size_t klen = 0;
        FILE *mem = open_memstream(&keep, &klen);
        while (a && fgets(line, sizeof line, a)) {
            char x[64], y[64];
            int d;
            if (sscanf(line, "%63s %63s %d", x, y, &d) == 3 && d <= D) {
                fputs(line, stdout);
                fputs(line, mem);
            }
        }
        if (a)
            fclose(a);
        fclose(mem);
        a = fopen(path2, "w");
        write_all(a, keep, klen, path2);
        close_synced(a, path2);
        free(keep);
    }
    fflush(stdout);

    for (int d = D; d < maxlen; d++) {
        double tl = omp_get_wtime();
        level_path(path2, sizeof path2, dir, d, "bin");
        size_t ncur = lev[d], nnext = lev[d + 1];
        int fd = open(path2, O_RDONLY);
        struct stat sb;
        if (fd < 0 || fstat(fd, &sb) || (size_t)sb.st_size != ncur * RS) {
            fprintf(stderr, "bigbottom: %s missing or not %zu records\n", path2, ncur);
            return 1;
        }
        uint8_t *cur = mmap(NULL, ncur * RS, PROT_READ, MAP_SHARED, fd, 0);
        if (cur == MAP_FAILED) {
            fprintf(stderr, "bigbottom: cannot map %s\n", path2);
            return 1;
        }
        madvise(cur, ncur * RS, MADV_SEQUENTIAL);

        /* rounds: each holds about nnext / R records and 3/2 as many slots */
        size_t per = RS + 3 * sizeof(islot_t) / 2;
        uint32_t R = (uint32_t)((nnext * per + budget - 1) / budget);
        if (R < 1)
            R = 1;
        part_t P;
        size_t part = nnext / R;
        P.partcap = part + part / 100 + 65536;   /* far above the spread of a hash split */
        P.cap = P.partcap + P.partcap / 2;
        P.rec = malloc(P.partcap * RS);
        P.tab = malloc(P.cap * sizeof(islot_t));
        if (!P.rec || !P.tab) {
            fprintf(stderr, "bigbottom: out of memory for a round of %zu records\n", P.partcap);
            return 1;
        }
        level_path(path, sizeof path, dir, d + 1, "tmp");
        FILE *out = fopen(path, "wb");
        if (!out) {
            fprintf(stderr, "bigbottom: cannot create %s\n", path);
            return 1;
        }
        uint64_t found = 0, fails = 0;
        uint64_t *av = NULL;
        size_t nav = 0, avcap = 0;
        for (uint32_t r = 0; r < R; r++) {
            memset(P.tab, 0, P.cap * sizeof(islot_t));
            P.cnt = 0;
            build_part(cur, ncur, &P, r, R);
            uint64_t f = 0;
#pragma omp parallel for schedule(dynamic, 4096) reduction(+ : f)
            for (size_t s = 0; s < P.cnt; s++) {
                const uint8_t *rec = P.rec + s * RS;
                if (rec[8 + N])
                    continue;
                f++;
                if (!bad10)
                    continue;
                uint64_t key;
                uint8_t w[BN];
                memcpy(&key, rec, 8);
                unpack(key, n, w);
                if (!contains(w, n)) {
#pragma omp critical
                    {
                        if (nav == avcap) {
                            avcap = avcap ? 2 * avcap : 1024;
                            av = realloc(av, avcap * sizeof *av);
                            if (!av) {
                                fprintf(stderr, "out of memory (avoiding tops)\n");
                                exit(1);
                            }
                        }
                        av[nav++] = key;
                    }
                }
            }
            fails += f;
            found += P.cnt;
            write_all(out, P.rec, P.cnt * RS, path);
        }
        close_synced(out, path);
        free(P.rec);
        free(P.tab);
        munmap(cur, ncur * RS);
        close(fd);
        if (found != nnext) {
            fprintf(stderr, "bigbottom: level %d has %llu nodes, the count says %zu\n", d + 1,
                    (unsigned long long)found, nnext);
            return 1;
        }
        level_path(path2, sizeof path2, dir, d + 1, "bin");
        if (rename(path, path2)) {
            fprintf(stderr, "bigbottom: cannot rename %s\n", path);
            return 1;
        }

        /* the avoiding pairs of this level, then the state line, then drop level d */
        qsort(av, nav, sizeof *av, cmp_u64);
        snprintf(path, sizeof path, "%s/avoid.txt", dir);
        FILE *a = fopen(path, "a");
        if (!a) {
            fprintf(stderr, "bigbottom: cannot write %s\n", path);
            return 1;
        }
        for (size_t k = 0; k < nav; k++) {
            char b[BN];
            fmt_big(av[k], b);
            printf("%s %s %d\n", us, b, d + 1);
            fprintf(a, "%s %s %d\n", us, b, d + 1);
        }
        fflush(stdout);
        close_synced(a, path);
        free(av);
        snprintf(path, sizeof path, "%s/state.txt", dir);
        st = fopen(path, "a");
        if (!st) {
            fprintf(stderr, "bigbottom: cannot write %s\n", path);
            return 1;
        }
        fprintf(st, "%d %zu %llu %zu %u %.1f\n", d + 1, nnext, (unsigned long long)fails, nav, R,
                omp_get_wtime() - tl);
        close_synced(st, path);
        level_path(path2, sizeof path2, dir, d, "bin");
        unlink(path2);
        fprintf(stderr, "  %4d %14zu %14llu %12zu\n", d + 1, nnext, (unsigned long long)fails, nav);
        tp += nnext;
        tf += fails;
        ta += nav;
    }
    level_path(path2, sizeof path2, dir, maxlen, "bin");
    unlink(path2);
    fprintf(stderr, "  %4s %14llu %14llu %12llu\n", "all", (unsigned long long)tp,
            (unsigned long long)tf, (unsigned long long)ta);
    fprintf(stderr, "# %.1f s in this run, %d threads, at most %.1f GB per round\n", omp_get_wtime() - t0,
            omp_get_max_threads(), memgb);
    free(bad10);
    free(cont11);
    return 0;
}
