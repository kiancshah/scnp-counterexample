/*
 * greedy.c - are greedy chains dominant?  (An-Tung-Zhang, Def. 3.5)
 *
 * A step x' -> x = x' t_ab of a chain from u is greedy unless some x'' >= u
 * with x'' -> x is a step with positions (a, b'), b' > b, or (a', b), a' < a.
 * This depends only on u and the step, so greedy chains are the paths of a
 * subgraph and the recursion of sweep.c carries over: with H and F as there,
 *   gsome(x) = OR  over greedy steps v -> x of gsome(v) && H(v)+delta == H(x)
 *   gall(x)  = AND over greedy steps v -> x reached by a greedy chain of
 *              gall(v) && H(v)+delta == H(x)
 * (a prefix of a chain attaining H(x) attains H(v), as in sweep.c), so
 * gsome = some greedy chain is dominant, gall = every greedy chain is.
 * ATZ prove gall = 1 for u = id and that greedy chains always exist.
 */
#include "scnp.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>



typedef struct {
    uint64_t pairs, good, gsome, gall, unreach;
    uint64_t av, av_good, av_gsome, av_gall;
    uint64_t id_pairs, id_gall;
    int nex;
    uint32_t exu[8], exw[8]; /* good, w avoids patterns, no greedy chain dominant */
} gstat_t;

static int vec_eq(vec_t a, vec_t b)
{
    return !memcmp(&a, &b, VL);
}

static void greedy_from(const tables_t *T, uint32_t u, const uint64_t *badmask,
                        const uint8_t *avoid, int32_t *slot, uint32_t *node[2],
                        vec_t *H[2], uint8_t *F[2], uint8_t *gs[2], uint8_t *ga[2],
                        uint8_t *rc[2], uint64_t *inmask, gstat_t *st, uint8_t *pairgs)
{
    if (pairgs)
        pairgs[(size_t)u * T->nperm + u] = 3; /* empty chain: greedy, dominant */
    int cur = 0;
    uint32_t ncur = 1;
    node[0][0] = u;
    memset(&H[0][0], 0, VL);
    F[0][0] = gs[0][0] = ga[0][0] = rc[0][0] = 1;
    st->pairs++;
    st->good++;
    st->gsome++;
    st->gall++;
    if (avoid[u]) {
        st->av++;
        st->av_good++;
        st->av_gsome++;
        st->av_gall++;
    }
    if (u == T->id) {
        st->id_pairs++;
        st->id_gall++;
    }
    for (int d = 0; ncur; d++) {
        int nx = cur ^ 1;
        uint32_t nn = 0;
        /* phase 1: H, F and the set of incoming step types */
        for (uint32_t i = 0; i < ncur; i++) {
            uint32_t v = node[cur][i];
            for (uint32_t e = T->off[v]; e < T->off[v + 1]; e++) {
                uint32_t x = T->cov[e];
                int t = T->ctr[e];
                vec_t cand = H[cur][i] + T->delta[t];
                int32_t s = slot[x];
                if (s < 0) {
                    s = (int32_t)nn++;
                    slot[x] = s;
                    node[nx][s] = x;
                    H[nx][s] = cand;
                    F[nx][s] = F[cur][i];
                    inmask[s] = 0;
                } else {
                    vec_t h = H[nx][s];
                    vec_t m = h;
                    uint8_t *pm = (uint8_t *)&m;
                    const uint8_t *pc = (const uint8_t *)&cand;
                    for (int k = 0; k < VL; k++)
                        if (pc[k] < pm[k])
                            pm[k] = pc[k];
                    int eh = vec_eq(m, h), ec = vec_eq(m, cand);
                    if (eh && ec)
                        F[nx][s] |= F[cur][i];
                    else if (ec) {
                        H[nx][s] = cand;
                        F[nx][s] = F[cur][i];
                    } else if (!eh) {
                        H[nx][s] = m;
                        F[nx][s] = 0;
                    }
                }
                inmask[s] |= 1ull << t;
            }
        }
        for (uint32_t s = 0; s < nn; s++) {
            gs[nx][s] = 0;
            ga[nx][s] = 1;
            rc[nx][s] = 0;
        }
        /* phase 2: greedy steps */
        for (uint32_t i = 0; i < ncur; i++) {
            if (!rc[cur][i])
                continue;
            uint32_t v = node[cur][i];
            for (uint32_t e = T->off[v]; e < T->off[v + 1]; e++) {
                int t = T->ctr[e];
                int32_t s = slot[T->cov[e]];
                if (inmask[s] & badmask[t])
                    continue; /* not greedy */
                int eq = vec_eq(H[cur][i] + T->delta[t], H[nx][s]);
                rc[nx][s] = 1;
                gs[nx][s] |= (uint8_t)(gs[cur][i] && eq);
                ga[nx][s] &= (uint8_t)(ga[cur][i] && eq);
            }
        }
        for (uint32_t s = 0; s < nn; s++) {
            uint32_t x = node[nx][s];
            slot[x] = -1;
            int g = F[nx][s], s1 = gs[nx][s], sa = ga[nx][s] && rc[nx][s];
            if (pairgs)
                pairgs[(size_t)u * T->nperm + x] = (uint8_t)(1 | s1 << 1);
            st->pairs++;
            st->good += (uint64_t)g;
            st->gsome += (uint64_t)s1;
            st->gall += (uint64_t)sa;
            st->unreach += (uint64_t)!rc[nx][s];
            if (avoid[x]) {
                st->av++;
                st->av_good += (uint64_t)g;
                st->av_gsome += (uint64_t)s1;
                st->av_gall += (uint64_t)sa;
                if (g && !s1 && st->nex < 8) {
                    st->exu[st->nex] = u;
                    st->exw[st->nex] = x;
                    st->nex++;
                }
            }
            if (u == T->id) {
                st->id_pairs++;
                st->id_gall += (uint64_t)sa;
            }
        }
        cur = nx;
        ncur = nn;
        (void)d;
    }
}

int run_greedy(int n, const char *patterns)
{
    tables_t T;
    tables_build(&T, n);
    uint8_t pat[16][MAXN];
    int plen[16], npat = 0;
    char buf[1024];
    snprintf(buf, sizeof buf, "%s", patterns);
    for (char *tok = strtok(buf, ", "); tok && npat < 16; tok = strtok(NULL, ", "))
        if ((plen[npat] = parse_perm(tok, pat[npat])) > 0)
            npat++;
    uint8_t *avoid = malloc(T.nperm);
#pragma omp parallel for schedule(dynamic, 4096)
    for (uint32_t v = 0; v < T.nperm; v++) {
        int any = 0;
        for (int k = 0; k < npat && !any; k++)
            any = perm_contains(PERM(&T, v), n, pat[k], plen[k]);
        avoid[v] = (uint8_t)!any;
    }
    uint64_t badmask[MAXT] = {0};
    for (int t = 0; t < T.ntr; t++)
        for (int s = 0; s < T.ntr; s++)
            if ((T.ta[s] == T.ta[t] && T.tb[s] > T.tb[t]) ||
                (T.tb[s] == T.tb[t] && T.ta[s] < T.ta[t]))
                badmask[t] |= 1ull << s;

    gstat_t tot;
    memset(&tot, 0, sizeof tot);
    /* n <= 7: keep "some greedy chain dominant" per pair, to combine (u,w)
     * with (w0 w, w0 u): a top-down greedy chain of the latter is a chain of
     * the former that is greedy from the bottom (longest swap going up).   */
    uint8_t *pairgs = n <= 7 ? calloc((size_t)T.nperm * T.nperm, 1) : NULL;
#pragma omp parallel
    {
        int32_t *slot = malloc(sizeof(int32_t) * T.nperm);
        memset(slot, 0xff, sizeof(int32_t) * T.nperm);
        uint32_t *node[2];
        vec_t *H[2];
        uint8_t *F[2], *gs[2], *ga[2], *rc[2];
        for (int k = 0; k < 2; k++) {
            node[k] = malloc(sizeof(uint32_t) * T.maxlevel);
            H[k] = aligned_alloc(64, sizeof(vec_t) * T.maxlevel);
            F[k] = malloc(T.maxlevel);
            gs[k] = malloc(T.maxlevel);
            ga[k] = malloc(T.maxlevel);
            rc[k] = malloc(T.maxlevel);
        }
        uint64_t *inmask = malloc(sizeof(uint64_t) * T.maxlevel);
        gstat_t st;
        memset(&st, 0, sizeof st);
#pragma omp for schedule(dynamic, 16)
        for (uint32_t u = 0; u < T.nperm; u++)
            greedy_from(&T, u, badmask, avoid, slot, node, H, F, gs, ga, rc, inmask, &st, pairgs);
#pragma omp critical
        {
            tot.pairs += st.pairs;
            tot.good += st.good;
            tot.gsome += st.gsome;
            tot.gall += st.gall;
            tot.unreach += st.unreach;
            tot.av += st.av;
            tot.av_good += st.av_good;
            tot.av_gsome += st.av_gsome;
            tot.av_gall += st.av_gall;
            tot.id_pairs += st.id_pairs;
            tot.id_gall += st.id_gall;
            for (int i = 0; i < st.nex && tot.nex < 8; i++) {
                tot.exu[tot.nex] = st.exu[i];
                tot.exw[tot.nex] = st.exw[i];
                tot.nex++;
            }
        }
        free(slot);
        free(inmask);
        for (int k = 0; k < 2; k++) {
            free(node[k]);
            free(H[k]);
            free(F[k]);
            free(gs[k]);
            free(ga[k]);
            free(rc[k]);
        }
    }
    printf("# S_%d: greedy chains (An-Tung-Zhang Def. 3.5), patterns %s\n", n, patterns);
    printf("pairs u <= w                                   %llu\n", (unsigned long long)tot.pairs);
    printf("  with SCNP                                    %llu\n", (unsigned long long)tot.good);
    printf("  some greedy chain dominant                   %llu\n", (unsigned long long)tot.gsome);
    printf("  every greedy chain dominant                  %llu\n", (unsigned long long)tot.gall);
    printf("  no greedy chain at all (should be 0)         %llu\n", (unsigned long long)tot.unreach);
    printf("pairs with w avoiding the patterns             %llu\n", (unsigned long long)tot.av);
    printf("  with SCNP                                    %llu\n", (unsigned long long)tot.av_good);
    printf("  some greedy chain dominant                   %llu\n", (unsigned long long)tot.av_gsome);
    printf("  every greedy chain dominant                  %llu\n", (unsigned long long)tot.av_gall);
    printf("pairs with u = id                              %llu\n", (unsigned long long)tot.id_pairs);
    printf("  every greedy chain dominant (ATZ: all)       %llu\n", (unsigned long long)tot.id_gall);
    for (int i = 0; i < tot.nex; i++) {
        char a[MAXN + 1], b[MAXN + 1];
        format_perm(PERM(&T, tot.exu[i]), n, a);
        format_perm(PERM(&T, tot.exw[i]), n, b);
        printf("example: u=%s w=%s  SCNP, w avoids the patterns, no greedy chain dominant\n", a, b);
    }
    if (pairgs) {
        uint64_t av = 0, top = 0, bot = 0, either = 0, neither = 0;
        int shown = 0;
        for (uint32_t u = 0; u < T.nperm; u++)
            for (uint32_t w = 0; w < T.nperm; w++) {
                uint8_t p = pairgs[(size_t)u * T.nperm + w];
                if (!p || !avoid[w])
                    continue;
                uint8_t q = pairgs[(size_t)T.lw0[w] * T.nperm + T.lw0[u]];
                int t = p >> 1 & 1, b = q >> 1 & 1;
                av++;
                top += (uint64_t)t;
                bot += (uint64_t)b;
                either += (uint64_t)(t | b);
                if (!(t | b)) {
                    neither++;
                    if (shown < 4) {
                        char a[MAXN + 1], c[MAXN + 1];
                        format_perm(PERM(&T, u), n, a);
                        format_perm(PERM(&T, w), n, c);
                        printf("example: u=%s w=%s  w avoids the patterns, no greedy chain from either end dominant\n", a, c);
                        shown++;
                    }
                }
            }
        printf("pairs with w avoiding the patterns             %llu\n", (unsigned long long)av);
        printf("  greedy from the top dominant                 %llu\n", (unsigned long long)top);
        printf("  greedy from the bottom dominant              %llu\n", (unsigned long long)bot);
        printf("  greedy from one of the two ends dominant     %llu\n", (unsigned long long)either);
        printf("  neither                                      %llu\n", (unsigned long long)neither);
        free(pairgs);
    }
    free(avoid);
    tables_free(&T);
    return 0;
}
