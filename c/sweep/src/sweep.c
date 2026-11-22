/*
 * sweep.c - the fast SCNP criterion.
 *
 * For a chain C let P_C = sum of the simplices Delta_[a, b-1] over its steps;
 * supp(m_C) = P_C ∩ Z^{n-1}.  P_C = { t : |t| = l, t(I) >= g_C(I) for every
 * interval I } (Postnikov, Permutohedra..., Prop. 6.3; non-interval subsets
 * are implied because g_C is additive over connected components), so
 *     supp(m_C') ⊆ supp(m_C)  <=>  g_C <= g_C'  componentwise on intervals.
 * Since supp(D_{u,w}) is the union of the supp(m_C), SCNP holds iff the
 * vectors g_C have a least element, i.e. iff some chain attains
 * H[w] = componentwise min of g_C.  Such a chain ends with an edge v -> w
 * where its prefix attains H[v] and H[v] + delta = H[w], which gives the
 * recursion implemented in relax().
 */
#include "scnp.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __AVX512BW__
#include <immintrin.h>
#endif

static void *xalloc(size_t align, size_t sz)
{
    sz = (sz + align - 1) / align * align;
    void *p = aligned_alloc(align, sz ? sz : align);
    if (!p) {
        fprintf(stderr, "out of memory (%zu bytes)\n", sz);
        exit(1);
    }
    return p;
}

void sweep_work_init(sweep_work_t *W, const tables_t *T)
{
    W->slot = xalloc(64, sizeof(int32_t) * T->nperm);
    memset(W->slot, 0xff, sizeof(int32_t) * T->nperm);
    for (int k = 0; k < 2; k++) {
        W->node[k] = xalloc(64, sizeof(uint32_t) * T->maxlevel);
        W->H[k] = xalloc(64, sizeof(vec_t) * T->maxlevel);
        W->F[k] = xalloc(64, T->maxlevel);
    }
}

void sweep_work_free(sweep_work_t *W)
{
    free(W->slot);
    for (int k = 0; k < 2; k++) {
        free(W->node[k]);
        free(W->H[k]);
        free(W->F[k]);
    }
}

/* Fold a new predecessor (cand, fc) into the state (h, f) of a node, keeping
 *   h = componentwise min of the candidates seen so far,
 *   f = OR of fc over the candidates that are equal to h.                  */
static inline void relax(vec_t *h, uint8_t *f, vec_t cand, uint8_t fc)
{
#ifdef __AVX512BW__
    __m512i a = (__m512i)*h, c = (__m512i)cand;
    __mmask64 le = _mm512_cmple_epu8_mask(c, a);
    __mmask64 ge = _mm512_cmpge_epu8_mask(c, a);
    if (ge == ~0ULL) {        /* cand >= h */
        if (le == ~0ULL)      /* cand == h */
            *f |= fc;
        return;
    }
    if (le == ~0ULL) {        /* cand < h  */
        *h = cand;
        *f = fc;
        return;
    }
    *h = (vec_t)_mm512_min_epu8(a, c); /* incomparable */
    *f = 0;
#else
    svec_t lt = cand < *h, gt = cand > *h;
    uint64_t l[VL / 8], g[VL / 8], lo = 0, go = 0;
    memcpy(l, &lt, VL);
    memcpy(g, &gt, VL);
    for (int i = 0; i < VL / 8; i++) {
        lo |= l[i];
        go |= g[i];
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
    *h = (*h & ~(vec_t)lt) | (cand & (vec_t)lt);
    *f = 0;
#endif
}

uint64_t sweep_from(const tables_t *T, sweep_work_t *W, uint32_t u, int maxd,
                    level_fn fn, void *ctx)
{
    int cur = 0;
    uint32_t ncur = 1;
    uint64_t pushes = 0;
    int32_t *slot = W->slot;

    W->node[0][0] = u;
    memset(&W->H[0][0], 0, VL);
    W->F[0][0] = 1;
    fn(ctx, u, 0, W->node[0], W->F[0], W->H[0], 1);

    for (int d = 0; d < maxd && ncur; d++) {
        int nx = cur ^ 1;
        uint32_t nn = 0;
        const uint32_t *cn = W->node[cur];
        const vec_t *ch = W->H[cur];
        const uint8_t *cf = W->F[cur];
        uint32_t *xn = W->node[nx];
        vec_t *xh = W->H[nx];
        uint8_t *xf = W->F[nx];

        for (uint32_t i = 0; i < ncur; i++) {
            uint32_t v = cn[i];
            vec_t hv = ch[i];
            uint8_t fv = cf[i];
            uint32_t e0 = T->off[v], e1 = T->off[v + 1];
            for (uint32_t e = e0; e < e1; e++) {
                uint32_t x = T->cov[e];
                vec_t cand = hv + T->delta[T->ctr[e]];
                int32_t s = slot[x];
                if (s < 0) {
                    s = (int32_t)nn++;
                    slot[x] = s;
                    xn[s] = x;
                    xh[s] = cand;
                    xf[s] = fv;
                } else {
                    relax(&xh[s], &xf[s], cand, fv);
                }
            }
            pushes += e1 - e0;
        }
        for (uint32_t s = 0; s < nn; s++)
            slot[xn[s]] = -1;
        cur = nx;
        ncur = nn;
        fn(ctx, u, d + 1, xn, xf, xh, nn);
    }
    return pushes;
}
