/* main.c - command line interface. */
#include "scnp.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void usage(void)
{
    fprintf(stderr,
        "usage:\n"
        "  scnp pair U W [-v|-vv] [--limit K]   report one pair (chains, supports, SCNP)\n"
        "  scnp pair - [--limit K]             same, one line per \"U W ...\" line on stdin\n"
        "  scnp check -n N [--maxlen L] [--snp] brute-force cross-check of every pair in S_N\n"
        "  scnp sweep -n N -o DIR [--chunks K] [--maxlen L] [--nosym]\n"
        "             [--list all|reps|none] [--gzip]\n"
        "                                       checkpointed SCNP sweep of S_N (resumable)\n"
        "  scnp bottom U [--maxlen L] [-q]      all pairs (U, w), w >= U, from one sweep\n"
        "  scnp summary DIR                     aggregate the statistics of a sweep\n"
        "  scnp greedy -n N [--patterns P,..]   are greedy chains (ATZ Def. 3.5) dominant?\n"
        "  scnp patterns DIR [--patterns P,Q..] compare bad w (some u with D_u^w not SCNP)\n"
        "                                       with containing one of the patterns\n"
        "  scnp search -o DIR --bottoms FILE [--patterns P,..] [--prev BADW.bin.gz] [--hours H]\n"
        "              [--chunks K]               (K bottoms per checkpoint, default 16)\n"
        "                                       fast criterion from every bottom in FILE\n"
        "                                       (n <= 11), checkpointed; lists the pairs\n"
        "                                       whose top avoids all the patterns\n"
        "  scnp bigbottom U [--prev BADW10.bin.gz] [--maxlen L] [--tmp DIR [--mem GB]]\n"
        "                                       like bottom, for n <= 16 without tables\n"
        "                                       (letters 1..9, A..G); with --prev, lists the\n"
        "                                       failing pairs whose top avoids the patterns;\n"
        "                                       with --tmp, levels in files under DIR, rounds\n"
        "                                       of at most GB (default 10), all threads,\n"
        "                                       resumable\n"
        "\n"
        "Permutations are one-line, 1-based: 214365, 21436587A9 ('A' or '0' = 10,\n"
        "'B' = 11)\n"
        "or comma separated.  Threads: OMP_NUM_THREADS.\n");
    exit(1);
}

static long num(const char *s)
{
    char *end;
    long v = strtol(s, &end, 10);
    if (!*s || *end)
        usage();
    return v;
}

int main(int argc, char **argv)
{
    if (argc < 2)
        usage();
    const char *cmd = argv[1];
    int n = 0, maxlen = -1, snp = 0, nosym = 0, verbose = 0, chunks = 0;
    int list = 2, gz = 0, quiet = 0;
    size_t limit = 50000000;
    const char *dir = NULL, *pos[2] = {NULL, NULL}, *bottoms = NULL, *prev = NULL, *tmp = NULL;
    double hours = 0, mem = 10;
    const char *patterns = "4231,541632,5276143,5472163";
    int npos = 0;
    for (int i = 2; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-n") && i + 1 < argc)
            n = (int)num(argv[++i]);
        else if (!strcmp(a, "-o") && i + 1 < argc)
            dir = argv[++i];
        else if (!strcmp(a, "--maxlen") && i + 1 < argc)
            maxlen = (int)num(argv[++i]);
        else if (!strcmp(a, "--chunks") && i + 1 < argc)
            chunks = (int)num(argv[++i]);
        else if (!strcmp(a, "--limit") && i + 1 < argc)
            limit = (size_t)num(argv[++i]);
        else if (!strcmp(a, "--patterns") && i + 1 < argc)
            patterns = argv[++i];
        else if (!strcmp(a, "--bottoms") && i + 1 < argc)
            bottoms = argv[++i];
        else if (!strcmp(a, "--prev") && i + 1 < argc)
            prev = argv[++i];
        else if (!strcmp(a, "--hours") && i + 1 < argc)
            hours = atof(argv[++i]);
        else if (!strcmp(a, "--tmp") && i + 1 < argc)
            tmp = argv[++i];
        else if (!strcmp(a, "--mem") && i + 1 < argc)
            mem = atof(argv[++i]);
        else if (!strcmp(a, "--snp"))
            snp = 1;
        else if (!strcmp(a, "--nosym"))
            nosym = 1;
        else if (!strcmp(a, "--gzip"))
            gz = 1;
        else if (!strcmp(a, "-q"))
            quiet = 1;
        else if (!strcmp(a, "--list") && i + 1 < argc) {
            const char *m = argv[++i];
            list = !strcmp(m, "none") ? 0 : !strcmp(m, "reps") ? 1 : !strcmp(m, "all") ? 2 : -1;
            if (list < 0)
                usage();
        }
        else if (!strcmp(a, "-v"))
            verbose = 1;
        else if (!strcmp(a, "-vv"))
            verbose = 2;
        else if (npos < 2 && (a[0] != '-' || !a[1]))
            pos[npos++] = a;
        else
            usage();
    }

    if (!strcmp(cmd, "pair")) {
        if (npos == 1 && !strcmp(pos[0], "-"))
            return run_pairs_stdin(limit);
        if (npos != 2)
            usage();
        uint8_t u[MAXN], w[MAXN];
        int n1 = parse_perm(pos[0], u), n2 = parse_perm(pos[1], w);
        if (n1 < 0 || n2 < 0 || n1 != n2) {
            fprintf(stderr, "U and W must be permutations of the same size <= %d\n", MAXN);
            return 1;
        }
        tables_t T;
        tables_build(&T, n1);
        int rc = run_pair(&T, u, w, verbose ? verbose : 1, limit);
        tables_free(&T);
        return rc;
    }
    if (!strcmp(cmd, "check")) {
        if (n < 1 || n > 7) {
            fprintf(stderr, "check: need 1 <= n <= 7\n");
            return 1;
        }
        return run_check(n, maxlen, snp);
    }
    if (!strcmp(cmd, "sweep")) {
        if (n < 1 || n > MAXN || !dir) {
            fprintf(stderr, "sweep: need -n (1..%d) and -o DIR\n", MAXN);
            return 1;
        }
        if (chunks <= 0)
            chunks = n <= 6 ? 16 : n == 7 ? 64 : n == 8 ? 256 : n == 9 ? 1024 : n == 10 ? 4096 : 16384;
        return run_sweep(dir, n, !nosym, chunks, maxlen, list, gz);
    }
    if (!strcmp(cmd, "bottom")) {
        uint8_t u[MAXN];
        int n1 = npos == 1 ? parse_perm(pos[0], u) : -1;
        if (n1 < 0)
            usage();
        return run_bottom(u, n1, maxlen, quiet);
    }
    if (!strcmp(cmd, "greedy")) {
        if (n < 1 || n > MAXN)
            usage();
        return run_greedy(n, patterns);
    }
    if (!strcmp(cmd, "patterns")) {
        if (npos != 1)
            usage();
        return run_patterns(pos[0], patterns);
    }
    if (!strcmp(cmd, "bigbottom")) {
        if (npos != 1)
            usage();
        if (tmp)
            return run_bigsplit(pos[0], prev, maxlen, tmp, mem);
        return run_bigbottom(pos[0], prev, maxlen);
    }
    if (!strcmp(cmd, "search")) {
        if (!dir || !bottoms)
            usage();
        return run_search(dir, bottoms, patterns, prev, hours, chunks);
    }
    if (!strcmp(cmd, "summary")) {
        if (npos != 1)
            usage();
        return run_summary(pos[0]);
    }
    usage();
    return 1;
}
