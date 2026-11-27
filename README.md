# A counterexample to the SCNP conjecture

Code and data for the note "A counterexample to the conjecture on single-chain Newton polytopes of
Postnikov-Stanley polynomials" (K. C. Shah, 2026).

For `u = 214365` and `w = 541632` the Postnikov-Stanley polynomial `D_u^w` has 81 monomials and no saturated
chain produces more than 73 of them, although `w` avoids `4231`.

## Check it

    python3 python/verify.py          # plain Python, a few seconds (Appendix A of the note)
    python3 python/verify_sage.py     # the same with SageMath's Bruhat order (Appendix B)

## Contents

| path | content |
|---|---|
| `python/verify.py`, `python/verify_sage.py` | the two programs printed in the note |
| `python/core.py` | definitions: covers, chains, supports, the window criterion |
| `python/spiral.py` | atoms, atom graphs, rotations in the cubes |
| `python/boolean_intervals.py` | bad cubes and loops, `n <= 8` |
| `python/figures.py` | spiral and rosette figures |
| `c/sweep/` | C program that decides every pair of `S_n`, `n <= 10`, going up from every bottom |
| `c/tops/` | second C program that decides every pair below a given top, `n <= 11` |
| `data/` | outputs of these programs |

    cd c/sweep && make && ./scnp pair 214365 541632
    cd c/tops && make -C fast && fast/scnp_fast top 541632 -v

Licence: MIT.
