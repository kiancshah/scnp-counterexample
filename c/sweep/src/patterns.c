/*
 * patterns.c - compare the set of "bad" w (some u <= w with D_u^w not SCNP,
 * recorded by `scnp sweep` in DIR/badw.bin) with pattern containment.
 *
 * Default patterns: 4231, 541632, 5276143, 5472163 (the conjectured list).
 * The first pattern is also reported on its own.
 */
#include "scnp.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXPAT 16

/* 1 iff w (length n) contains p (length k): choose positions left to right,
 * checking the relative order against every position chosen so far.    */
static int contains_rec(const uint8_t *w, int n, const uint8_t *p, int k,
                        int d, int start, int *pos)
{
    if (d == k)
        return 1;
    for (int i = start; i <= n - (k - d); i++) {
        int ok = 1;
        for (int e = 0; e < d && ok; e++)
            ok = (w[i] > w[pos[e]]) == (p[d] > p[e]);
        if (!ok)
            continue;
        pos[d] = i;
        if (contains_rec(w, n, p, k, d + 1, i + 1, pos))
            return 1;
    }
    return 0;
}

int perm_contains(const uint8_t *w, int n, const uint8_t *p, int k)
{
    int pos[MAXN];
    return k <= n && contains_rec(w, n, p, k, 0, 0, pos);
}

int run_patterns(const char *dir, const char *patterns)
{
    char path[4096];
    snprintf(path, sizeof path, "%s/config", dir);
    FILE *f = fopen(path, "r");
    int n = 0, sym, nchunks, maxlen;
    if (!f || fscanf(f, "scnp-sweep v1 n=%d sym=%d chunks=%d maxlen=%d", &n,
                     &sym, &nchunks, &maxlen) != 4) {
        fprintf(stderr, "%s: missing or unreadable config\n", dir);
        return 1;
    }
    fclose(f);

    uint8_t pat[MAXPAT][MAXN];
    int plen[MAXPAT], npat = 0;
    char buf[1024];
    snprintf(buf, sizeof buf, "%s", patterns);
    for (char *tok = strtok(buf, ", "); tok; tok = strtok(NULL, ", ")) {
        if (npat == MAXPAT || (plen[npat] = parse_perm(tok, pat[npat])) < 0) {
            fprintf(stderr, "bad pattern list: %s\n", patterns);
            return 1;
        }
        npat++;
    }
    if (!npat)
        return 1;

    /* the run must be complete for "not bad" to mean anything */
    int done = 0;
    for (int c = 0; c < nchunks; c++) {
        snprintf(path, sizeof path, "%s/chunks/%05d.stat", dir, c);
        FILE *g = fopen(path, "r");
        if (g) {
            done++;
            fclose(g);
        }
    }

    uint32_t np = 1;
    for (int i = 2; i <= n; i++)
        np *= (uint32_t)i;
    size_t nbytes = (np + 7) / 8;
    uint8_t *bits = calloc(nbytes, 1);
    if (!load_badw(dir, bits, nbytes)) {
        fprintf(stderr, "%s/badw.bin not found (sweep from an older version?)\n", dir);
        return 1;
    }
    uint8_t *cls = malloc(np); /* bit0 bad, bit1 contains pattern 0, bit2 any */
#pragma omp parallel for schedule(dynamic, 4096)
    for (uint32_t v = 0; v < np; v++) {
        uint8_t w[MAXN];
        perm_unrank(v, n, w);
        int first = perm_contains(w, n, pat[0], plen[0]), any = first;
        for (int k = 1; k < npat && !any; k++)
            any = perm_contains(w, n, pat[k], plen[k]);
        cls[v] = (uint8_t)((bits[v >> 3] >> (v & 7) & 1) | first << 1 | any << 2);
    }

    uint64_t nbad = 0, nfirst = 0, nany = 0, unexplained = 0, notbad = 0;
    char p1[4096], p2[4096];
    snprintf(p1, sizeof p1, "%s/bad_avoiding_patterns.txt", dir);
    snprintf(p2, sizeof p2, "%s/pattern_but_not_bad.txt", dir);
    FILE *f1 = fopen(p1, "w"), *f2 = fopen(p2, "w");
    if (!f1 || !f2) {
        perror(dir);
        return 1;
    }
    for (uint32_t v = 0; v < np; v++) {
        int bad = cls[v] & 1, first = cls[v] >> 1 & 1, any = cls[v] >> 2 & 1;
        nbad += (uint64_t)bad;
        nfirst += (uint64_t)first;
        nany += (uint64_t)any;
        if (bad != any) {
            uint8_t w[MAXN];
            char s[MAXN + 1];
            perm_unrank(v, n, w);
            format_perm(w, n, s);
            fprintf(bad ? f1 : f2, "%s\n", s);
            if (bad)
                unexplained++;
            else
                notbad++;
        }
    }
    fclose(f1);
    fclose(f2);
    free(cls);
    free(bits);

    snprintf(path, sizeof path, "%s/patterns.txt", dir);
    FILE *out = fopen(path, "w");
    FILE *dst[2] = {stdout, out};
    for (int k = 0; k < 2; k++) {
        FILE *o = dst[k];
        if (!o)
            continue;
        fprintf(o, "# S_%d: bad w = some u <= w with D_u^w not SCNP\n", n);
        if (done != nchunks)
            fprintf(o, "# ** sweep INCOMPLETE (%d/%d chunks): bad-w counts are lower bounds **\n",
                    done, nchunks);
        fprintf(o, "# patterns: %s\n", patterns);
        fprintf(o, "permutations                        %u\n", np);
        fprintf(o, "bad w                               %llu\n", (unsigned long long)nbad);
        fprintf(o, "w containing %-22s %llu\n", strtok(strcpy(buf, patterns), ", "),
                (unsigned long long)nfirst);
        fprintf(o, "w containing one of the patterns    %llu\n", (unsigned long long)nany);
        fprintf(o, "w avoiding all patterns             %llu\n",
                (unsigned long long)(np - nany));
        fprintf(o, "bad w avoiding all patterns         %llu   (%s)\n",
                (unsigned long long)unexplained, "bad_avoiding_patterns.txt");
        fprintf(o, "w with a pattern but not bad        %llu   (%s)\n",
                (unsigned long long)notbad, "pattern_but_not_bad.txt");
        fprintf(o, "# result: %s\n",
                unexplained || notbad
                    ? "bad w != w containing a pattern"
                    : "bad w == w containing one of the patterns");
    }
    if (out)
        fclose(out);
    return (unexplained || notbad || done != nchunks) ? 2 : 0;
}
