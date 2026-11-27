# A counterexample to the SCNP conjecture

Code and data for the note "A counterexample to the conjecture on single-chain Newton polytopes of
Postnikov-Stanley polynomials" (K. C. Shah, 2026).

For `u = 214365` and `w = 541632` the Postnikov-Stanley polynomial `D_u^w` has 81 monomials and no saturated
chain produces more than 73 of them, although `w` avoids `4231`.

## Check it

    python3 python/verify.py          # plain Python 3, a few seconds (Appendix A of the note)
    python3 python/verify_sage.py     # the same with SageMath's Bruhat order (Appendix B)
    python3 python/verify_values.py   # labels on values instead of positions (Remark 2.2)

`verify_sage.py` needs the SageMath library, for example `pip install passagemath-combinat`
(it was run with version 10.8.14). `figures.py` needs `matplotlib`. Everything else is plain Python 3.
The C programs are built with `gcc -O3 -fopenmp` (see the Makefiles).

## Which program backs which statement

| in the note | program | output |
|---|---|---|
| Theorem 2.1, Proposition 2.2 | `python/verify.py`, `python/verify_sage.py` | printed, as in Appendices A and B |
| Remark 2.2, labels on values | `python/verify_values.py` | printed |
| Table 2 for `n <= 10`, Observation 3.1, Table 3 | `c/sweep` | `data/sweep/` |
| Table 2 for `n = 11`, Conjecture 3.3 for `n <= 11` | `c/tops` | `data/tops/` |
| Table 4, Observation 4.5 | `python/moves.py` | printed |
| Observation 4.3, competing chains of the cubes | `cd python && python3 spiral.py --all` | printed |
| bad cubes | `cd python && python3 boolean_intervals.py 4 5 6 7` | `data/boolean/` |
| Figure 6 | `cd python && python3 figures.py ../figures` | `figures/` |

Figures 1 to 5 of the note were drawn separately and are not generated here.

Single pairs take a moment:

    cd c/sweep && make && ./scnp pair 214365 541632
    cd c/tops && make -C fast && fast/scnp_fast top 541632 -v

The complete runs are long. The logs in `data/` record about 10 CPU hours for the sweep of `S_10`
(`data/sweep/summary_S10.txt`) and about 6 CPU hours for the tops of `S_11` avoiding `4231`
(`data/tops/n11_progress.json`).

Licence: MIT (see `LICENSE`).
