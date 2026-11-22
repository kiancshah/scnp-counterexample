/*
 * search.c - the fast criterion from a list of bottoms, n <= 11, checkpointed.
 *
 *   scnp search -o DIR --bottoms FILE [--patterns P,Q,..] [--prev BADW.bin.gz]
 *               [--hours H] [--chunks K]
 *
 * FILE has one bottom U per line, all in the same S_n ('#' lines are skipped).
 * The tables of S_n are built once and the bottoms are dealt to the threads.
 * For every failing pair (U, w) the top w is set in DIR/badw.bin (bit w&7 of
 * byte w>>3, w = lexicographic rank, as in `scnp sweep`), and if w avoids
 * every pattern the line "U w l(w)-l(U)" is written out.
 *
 * Bottoms are processed in chunks of K consecutive lines (default CHUNK = 16).
 * DIR/chunks/NNNNN.txt has, per bottom, "# U pairs fails avoiding seconds"
 * followed by its pairs with an avoiding top; the chunk is done iff
 * DIR/chunks/NNNNN.stat exists.  Killing the program loses only the chunk in
 * flight; rerunning the same command resumes.  --hours H stops starting new
 * chunks after H hours.
 *
 * Pattern containment is computed for all of S_n at once by deleting letters:
 * v contains a pattern iff v is one, or one of its one-letter deletions
 * contains one.  With --prev (the gzipped bad-w bitmap of S_{n-1}), the
 * containment in S_{n-1} is compared with it, and every bad w of S_n that
 * avoids all patterns is tested for minimality: all its one-letter deletions
 * must be good in S_{n-1}.  Such w are listed in DIR/avoiding_bad_w.txt.
 */
#include "scnp.h"

#include <omp.h>
#include <signal.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <zlib.h>

#define CHUNK 16    /* default bottoms per chunk */
#define MAXPAT 16

static volatile sig_atomic_t s_stop = 0;

static void on_stop(int sig)
{
    (void)sig;
    s_stop = 1;
}

static double wall(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static uint32_t factorial(int k)
{
    uint32_t f = 1;
    for (int i = 2; i <= k; i++)
        f *= (uint32_t)i;
    return f;
}

/* q = p (in S_k) with position i deleted and the values standardised */
static void delete_letter(const uint8_t *p, int k, int i, uint8_t *q)
{
    int m = 0;
    for (int j = 0; j < k; j++)
        if (j != i)
            q[m++] = (uint8_t)(p[j] - (p[j] > p[i]));
}

/* cont[v] = 1 iff v in S_n contains one of the patterns; *prev gets the same
 * array for S_{n-1} (NULL if n-1 is below the smallest pattern size).      */
static uint8_t *containment(int n, uint8_t pat[][MAXN], const int *plen,
                            int npat, uint8_t **prev)
{
    int kmin = n;
    for (int t = 0; t < npat; t++)
        if (plen[t] < kmin)
            kmin = plen[t];
    uint8_t *below = NULL, *cur = NULL;
    *prev = NULL;
    for (int k = kmin; k <= n; k++) {
        uint32_t np = factorial(k);
        cur = calloc(np, 1);
        if (!cur) {
            fprintf(stderr, "out of memory (containment)\n");
            exit(1);
        }
#pragma omp parallel for schedule(static)
        for (uint32_t v = 0; v < np; v++) {
            uint8_t p[MAXN], q[MAXN];
            perm_unrank(v, k, p);
            int c = 0;
            for (int t = 0; t < npat && !c; t++)
                c = plen[t] == k && !memcmp(p, pat[t], (size_t)k);
            for (int i = 0; i < k && !c && below; i++) {
                delete_letter(p, k, i, q);
                c = below[perm_rank(q, k - 1)];
            }
            cur[v] = (uint8_t)c;
        }
        if (k == n)
            *prev = below;
        else
            free(below);
        below = cur;
    }
    return cur;
}

typedef struct {
    const tables_t *T;
    const uint8_t *cont;
    uint8_t *badw;              /* shared bitmap of bad tops */
    uint64_t pairs, fails, avoid;
    double secs;
    char *buf;                  /* pairs with an avoiding top, as text */
    size_t len, cap;
} bres_t;

static void appendf(bres_t *r, const char *fmt, ...)
{
    char line[64];
    va_list ap;
    va_start(ap, fmt);
    int k = vsnprintf(line, sizeof line, fmt, ap);
    va_end(ap);
    if (r->len + (size_t)k + 1 > r->cap) {
        r->cap = r->cap ? 2 * r->cap : 4096;
        r->buf = realloc(r->buf, r->cap);
        if (!r->buf) {
            fprintf(stderr, "out of memory (output)\n");
            exit(1);
        }
    }
    memcpy(r->buf + r->len, line, (size_t)k + 1);
    r->len += (size_t)k;
}

static void on_search_level(void *vctx, uint32_t u, int d, const uint32_t *node,
                            const uint8_t *F, const vec_t *H, uint32_t cnt)
{
    (void)H;
    bres_t *r = vctx;
    const tables_t *T = r->T;
    r->pairs += cnt;
    for (uint32_t s = 0; s < cnt; s++) {
        if (F[s])
            continue;
        uint32_t w = node[s];
        r->fails++;
        __atomic_fetch_or(&r->badw[w >> 3], (uint8_t)(1u << (w & 7)), __ATOMIC_RELAXED);
        if (!r->cont[w]) {
            char a[MAXN + 1], b[MAXN + 1];
            format_perm(PERM(T, u), T->n, a);
            format_perm(PERM(T, w), T->n, b);
            r->avoid++;
            appendf(r, "%s %s %d\n", a, b, d);
        }
    }
}

static void write_file(const char *path, const char *text, size_t len)
{
    char tmp[4200];
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "w");
    if (!f || fwrite(text, 1, len, f) != len || fclose(f)) {
        perror(tmp);
        exit(1);
    }
    commit_file(tmp, path);
}

static int read_bottoms(const char *path, uint8_t (**out)[MAXN], int *n)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        perror(path);
        exit(1);
    }
    char line[256];
    int cnt = 0, cap = 0;
    uint8_t(*b)[MAXN] = NULL;
    *n = 0;
    while (fgets(line, sizeof line, f)) {
        line[strcspn(line, " \t\r\n")] = 0;
        if (!line[0] || line[0] == '#')
            continue;
        if (cnt == cap) {
            cap = cap ? 2 * cap : 1024;
            b = realloc(b, (size_t)cap * MAXN);
        }
        int k = parse_perm(line, b[cnt]);
        if (k < 0 || (*n && k != *n)) {
            fprintf(stderr, "%s: bad bottom '%s'\n", path, line);
            exit(1);
        }
        *n = k;
        cnt++;
    }
    fclose(f);
    *out = b;
    return cnt;
}

int run_search(const char *dir, const char *bfile, const char *patterns,
               const char *prevfile, double hours, int chunk)
{
    if (chunk <= 0)
        chunk = CHUNK;
    int n;
    uint8_t(*B)[MAXN];
    int nb = read_bottoms(bfile, &B, &n);
    if (!nb) {
        fprintf(stderr, "%s: no bottoms\n", bfile);
        return 1;
    }

    uint8_t pat[MAXPAT][MAXN];
    int plen[MAXPAT], npat = 0;
    char pbuf[1024];
    snprintf(pbuf, sizeof pbuf, "%s", patterns);
    for (char *tok = strtok(pbuf, ", "); tok; tok = strtok(NULL, ", ")) {
        if (npat == MAXPAT || (plen[npat] = parse_perm(tok, pat[npat])) < 0) {
            fprintf(stderr, "bad pattern list '%s'\n", patterns);
            return 1;
        }
        npat++;
    }

    char path[4096], line[4096];
    mkdir_p(dir);
    snprintf(path, sizeof path, "%s/chunks", dir);
    mkdir_p(path);
    snprintf(line, sizeof line, "scnp-search v1 n=%d bottoms=%d chunk=%d file=%s patterns=%s\n",
             n, nb, chunk, bfile, patterns);
    snprintf(path, sizeof path, "%s/config", dir);
    if (file_exists(path)) {
        char old[4096] = {0};
        FILE *f = fopen(path, "r");
        if (!f || !fgets(old, sizeof old, f)) {
            fprintf(stderr, "%s unreadable\n", path);
            return 1;
        }
        fclose(f);
        if (strcmp(old, line)) {
            fprintf(stderr, "%s was made with different parameters:\n  %s", path, old);
            return 1;
        }
    } else {
        write_file(path, line, strlen(line));
    }

    double t0 = wall();
    tables_t T;
    tables_build(&T, n);
    uint8_t *cprev;
    uint8_t *cont = containment(n, pat, plen, npat, &cprev);
    uint64_t ncont = 0;
    for (uint32_t v = 0; v < T.nperm; v++)
        ncont += cont[v];
    fprintf(stderr, "# S_%d: tables and pattern containment in %.1f s; %llu of %u permutations contain a pattern\n",
            n, wall() - t0, (unsigned long long)ncont, T.nperm);

    uint8_t *prevbad = NULL;
    uint32_t nprev = factorial(n - 1);
    if (prevfile) {
        size_t nbytes = (nprev + 7) / 8;
        prevbad = malloc(nbytes);
        gzFile g = gzopen(prevfile, "rb");
        if (!g || gzread(g, prevbad, (unsigned)nbytes) != (int)nbytes) {
            fprintf(stderr, "%s: cannot read %zu bytes\n", prevfile, nbytes);
            return 1;
        }
        gzclose(g);
        uint64_t mism = 0, nbad = 0;
        for (uint32_t v = 0; v < nprev; v++) {
            int bad = prevbad[v >> 3] >> (v & 7) & 1;
            nbad += (uint64_t)bad;
            mism += (uint64_t)(cprev && bad != cprev[v]);
        }
        fprintf(stderr, "# S_%d: %llu bad w in %s; bad != contains a pattern for %llu w\n",
                n - 1, (unsigned long long)nbad, prevfile, (unsigned long long)mism);
    }

    size_t nbytes = (T.nperm + 7) / 8;
    uint8_t *badw = calloc(nbytes, 1);
    load_badw(dir, badw, nbytes);

    int nth = omp_get_max_threads();
    sweep_work_t *W = malloc(sizeof(sweep_work_t) * (size_t)nth);
    for (int t = 0; t < nth; t++)
        sweep_work_init(&W[t], &T);
    signal(SIGINT, on_stop);
    signal(SIGTERM, on_stop);

    int nchunks = (nb + chunk - 1) / chunk, done_now = 0;
    for (int c = 0; c < nchunks && !s_stop; c++) {
        snprintf(path, sizeof path, "%s/chunks/%05d.stat", dir, c);
        if (file_exists(path))
            continue;
        if (hours > 0 && wall() - t0 > hours * 3600.0) {
            fprintf(stderr, "# time limit reached\n");
            break;
        }
        int lo = c * chunk, hi = lo + chunk < nb ? lo + chunk : nb;
        bres_t *res = calloc((size_t)chunk, sizeof(bres_t));
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = lo; i < hi; i++) {
            if (s_stop)
                continue;
            bres_t *r = &res[i - lo];
            r->T = &T;
            r->cont = cont;
            r->badw = badw;
            double ts = wall();
            sweep_from(&T, &W[omp_get_thread_num()], perm_rank(B[i], n), T.N,
                       on_search_level, r);
            r->secs = wall() - ts;
        }
        if (s_stop) {
            for (int i = 0; i < hi - lo; i++)
                free(res[i].buf);
            free(res);
            break;
        }
        bres_t all;
        memset(&all, 0, sizeof all);
        uint64_t P = 0, Fl = 0, A = 0;
        double S = 0;
        for (int i = lo; i < hi; i++) {
            bres_t *r = &res[i - lo];
            char a[MAXN + 1];
            format_perm(B[i], n, a);
            appendf(&all, "# %s %llu %llu %llu %.2f\n", a, (unsigned long long)r->pairs,
                    (unsigned long long)r->fails, (unsigned long long)r->avoid, r->secs);
            for (size_t k = 0; k < r->len;) {   /* copy line by line */
                size_t e = k;
                while (r->buf[e] != '\n')
                    e++;
                appendf(&all, "%.*s\n", (int)(e - k), r->buf + k);
                k = e + 1;
            }
            P += r->pairs;
            Fl += r->fails;
            A += r->avoid;
            S += r->secs;
            free(r->buf);
        }
        free(res);
        snprintf(path, sizeof path, "%s/chunks/%05d.txt", dir, c);
        write_file(path, all.buf ? all.buf : "", all.len);
        free(all.buf);
        char tmp[4200];
        snprintf(tmp, sizeof tmp, "%s/badw.bin.tmp", dir);
        snprintf(path, sizeof path, "%s/badw.bin", dir);
        FILE *f = fopen(tmp, "wb");
        if (!f || fwrite(badw, 1, nbytes, f) != nbytes || fclose(f)) {
            perror(tmp);
            exit(1);
        }
        commit_file(tmp, path);
        snprintf(line, sizeof line, "%d %llu %llu %llu %.2f\n", hi - lo, (unsigned long long)P,
                 (unsigned long long)Fl, (unsigned long long)A, S);
        snprintf(path, sizeof path, "%s/chunks/%05d.stat", dir, c);
        write_file(path, line, strlen(line));
        done_now++;
        fprintf(stderr, "# chunk %d/%d done (%d bottoms, %llu fails, %llu avoiding), %.0f s elapsed\n",
                c + 1, nchunks, hi - lo, (unsigned long long)Fl, (unsigned long long)A,
                wall() - t0);
    }

    /* totals over every finished chunk, and the avoiding bad w */
    uint64_t P = 0, Fl = 0, A = 0;
    int cdone = 0, bdone = 0;
    double S = 0;
    for (int c = 0; c < nchunks; c++) {
        snprintf(path, sizeof path, "%s/chunks/%05d.stat", dir, c);
        FILE *f = fopen(path, "r");
        if (!f)
            continue;
        int k;
        unsigned long long p, fl, a;
        double s;
        if (fscanf(f, "%d %llu %llu %llu %lf", &k, &p, &fl, &a, &s) == 5) {
            cdone++;
            bdone += k;
            P += p;
            Fl += fl;
            A += a;
            S += s;
        }
        fclose(f);
    }
    uint64_t nbad = 0, navoid = 0, nminimal = 0;
    snprintf(path, sizeof path, "%s/avoiding_bad_w.txt", dir);
    FILE *out = fopen(path, "w");
    for (uint32_t w = 0; w < T.nperm; w++) {
        if (!(badw[w >> 3] >> (w & 7) & 1))
            continue;
        nbad++;
        if (cont[w])
            continue;
        navoid++;
        int good_deletions = -1;
        if (prevbad) {
            uint8_t q[MAXN];
            good_deletions = 0;
            for (int i = 0; i < n; i++) {
                delete_letter(PERM(&T, w), n, i, q);
                uint32_t r = perm_rank(q, n - 1);
                good_deletions += !(prevbad[r >> 3] >> (r & 7) & 1);
            }
            nminimal += good_deletions == n;
        }
        char b[MAXN + 1];
        format_perm(PERM(&T, w), n, b);
        fprintf(out, "%s %s\n", b, good_deletions < 0 ? "?" : good_deletions == n ? "minimal" : "not-minimal");
    }
    fclose(out);
    fprintf(stderr,
            "# %s: chunks %d/%d (%d of %d bottoms; %d this run), pairs %llu, fails %llu, "
            "fails with avoiding top %llu, thread-seconds %.0f\n"
            "# distinct bad w %llu, avoiding all patterns %llu, of which minimal %llu%s\n",
            dir, cdone, nchunks, bdone, nb, done_now, (unsigned long long)P,
            (unsigned long long)Fl, (unsigned long long)A, S, (unsigned long long)nbad,
            (unsigned long long)navoid, (unsigned long long)nminimal,
            prevbad ? "" : " (no --prev: minimality not tested)");

    for (int t = 0; t < nth; t++)
        sweep_work_free(&W[t]);
    free(W);
    free(badw);
    free(cont);
    free(cprev);
    free(prevbad);
    free(B);
    tables_free(&T);
    return s_stop ? 3 : 0;
}
