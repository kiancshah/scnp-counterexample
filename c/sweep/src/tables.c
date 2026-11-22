/* tables.c - permutations, Bruhat covers and symmetry maps for S_n. */
#include "scnp.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Lexicographic rank via the Lehmer code, r = sum_i c_i (n-1-i)!. */
uint32_t perm_rank(const uint8_t *p, int n)
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

void perm_unrank(uint32_t r, int n, uint8_t *p)
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

int perm_length(const uint8_t *p, int n)
{
    int l = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            l += p[i] > p[j];
    return l;
}

/* Tableau criterion: u <= w iff #{i<=k : u(i) >= j} <= #{i<=k : w(i) >= j}. */
int bruhat_leq(const uint8_t *u, const uint8_t *w, int n)
{
    int cu[MAXN] = {0}, cw[MAXN] = {0};
    for (int k = 0; k < n; k++) {
        for (int j = 0; j <= u[k]; j++)
            cu[j]++;
        for (int j = 0; j <= w[k]; j++)
            cw[j]++;
        for (int j = 0; j < n; j++)
            if (cu[j] > cw[j])
                return 0;
    }
    return 1;
}

/* Accepts "214365", "21436587A9" ('A' or '0' = 10, 'B' = 11) or "2,1,4,3,6,5". */
int parse_perm(const char *s, uint8_t *p)
{
    int n = 0;
    if (strchr(s, ',')) {
        const char *q = s;
        while (*q) {
            char *end;
            long v = strtol(q, &end, 10);
            if (end == q || v < 1 || v > MAXN || n >= MAXN)
                return -1;
            p[n++] = (uint8_t)(v - 1);
            q = end;
            if (*q == ',')
                q++;
        }
    } else {
        for (; *s; s++) {
            int v;
            if (*s >= '1' && *s <= '9')
                v = *s - '0';
            else if (*s == 'A' || *s == 'a' || *s == '0')
                v = 10;
            else if (*s == 'B' || *s == 'b')
                v = 11;
            else
                return -1;
            if (n >= MAXN)
                return -1;
            p[n++] = (uint8_t)(v - 1);
        }
    }
    unsigned seen = 0;
    for (int i = 0; i < n; i++) {
        if (p[i] >= n || (seen >> p[i] & 1u))
            return -1;
        seen |= 1u << p[i];
    }
    return n > 0 ? n : -1;
}

void format_perm(const uint8_t *p, int n, char *out)
{
    for (int i = 0; i < n; i++)
        out[i] = p[i] < 9 ? (char)('1' + p[i]) : (char)('A' + p[i] - 9);
    out[n] = 0;
}

static int count_covers(const uint8_t *p, int n)
{
    int c = 0;
    for (int a = 0; a < n - 1; a++) {
        unsigned lim = (unsigned)n; /* min value > p[a] strictly between */
        for (int b = a + 1; b < n; b++)
            if (p[b] > p[a] && p[b] < lim) {
                c++;
                lim = p[b];
            }
    }
    return c;
}

void tables_build(tables_t *T, int n)
{
    memset(T, 0, sizeof *T);
    T->n = n;
    T->N = n * (n - 1) / 2;
    uint32_t np = 1;
    for (int i = 2; i <= n; i++)
        np *= (uint32_t)i;
    T->nperm = np;

    int t = 0;
    int8_t tid[MAXN][MAXN];
    for (int a = 0; a < n; a++)
        for (int b = a + 1; b < n; b++) {
            T->ta[t] = (uint8_t)a;
            T->tb[t] = (uint8_t)b;
            tid[a][b] = (int8_t)t++;
        }
    T->ntr = t;
    int I = 0;
    for (int i = 0; i + 1 < n; i++)
        for (int j = i; j + 1 < n; j++) {
            T->ii[I] = (uint8_t)i;
            T->ij[I] = (uint8_t)j;
            I++;
        }
    T->nint = I;
    for (t = 0; t < T->ntr; t++) {
        uint8_t d[VL] = {0};
        for (I = 0; I < T->nint; I++)
            d[I] = T->ii[I] <= T->ta[t] && T->tb[t] - 1 <= T->ij[I];
        memcpy(&T->delta[t], d, VL);
    }

    T->perm = malloc((size_t)np * n);
    T->len = malloc(np);
    T->off = malloc(sizeof(uint32_t) * ((size_t)np + 1));
    T->lw0 = malloc(sizeof(uint32_t) * np);
    T->rw0 = malloc(sizeof(uint32_t) * np);
    T->cw0 = malloc(sizeof(uint32_t) * np);
    if (!T->perm || !T->len || !T->off || !T->lw0 || !T->rw0 || !T->cw0) {
        fprintf(stderr, "out of memory building tables\n");
        exit(1);
    }

#pragma omp parallel for schedule(static)
    for (uint32_t v = 0; v < np; v++) {
        uint8_t *p = T->perm + (size_t)v * n, q[MAXN];
        perm_unrank(v, n, p);
        T->len[v] = (uint8_t)perm_length(p, n);
        T->off[v + 1] = (uint32_t)count_covers(p, n);
        for (int i = 0; i < n; i++)
            q[i] = (uint8_t)(n - 1 - p[i]);
        T->lw0[v] = perm_rank(q, n);
        for (int i = 0; i < n; i++)
            q[i] = p[n - 1 - i];
        T->rw0[v] = perm_rank(q, n);
        for (int i = 0; i < n; i++)
            q[i] = (uint8_t)(n - 1 - p[n - 1 - i]);
        T->cw0[v] = perm_rank(q, n);
    }
    T->off[0] = 0;
    for (uint32_t v = 0; v < np; v++)
        T->off[v + 1] += T->off[v];
    uint32_t E = T->off[np];
    T->cov = malloc(sizeof(uint32_t) * (size_t)E);
    T->ctr = malloc(E ? E : 1);
    if (!T->cov || !T->ctr) {
        fprintf(stderr, "out of memory building cover table\n");
        exit(1);
    }

#pragma omp parallel for schedule(static)
    for (uint32_t v = 0; v < np; v++) {
        const uint8_t *p = T->perm + (size_t)v * n;
        uint8_t q[MAXN];
        uint32_t e = T->off[v];
        for (int a = 0; a < n - 1; a++) {
            unsigned lim = (unsigned)n;
            for (int b = a + 1; b < n; b++)
                if (p[b] > p[a] && p[b] < lim) {
                    lim = p[b];
                    memcpy(q, p, (size_t)n);
                    q[a] = p[b];
                    q[b] = p[a];
                    T->cov[e] = perm_rank(q, n);
                    T->ctr[e] = (uint8_t)tid[a][b];
                    e++;
                }
        }
    }

    uint32_t lev[MAXL] = {0};
    for (uint32_t v = 0; v < np; v++)
        lev[T->len[v]]++;
    for (int l = 0; l <= T->N; l++)
        if (lev[l] > T->maxlevel)
            T->maxlevel = lev[l];
    T->id = 0;
    T->w0 = np - 1;
}

void tables_free(tables_t *T)
{
    free(T->perm);
    free(T->len);
    free(T->off);
    free(T->cov);
    free(T->ctr);
    free(T->lw0);
    free(T->rw0);
    free(T->cw0);
    memset(T, 0, sizeof *T);
}
