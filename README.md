# Inversion tokens and single-chain Newton polytopes of Postnikov-Stanley polynomials

Code, data and figures for the paper of the same name.

For permutations `u <= w` in Bruhat order, the Postnikov-Stanley polynomial `D_u^w` is a sum over the saturated
chains from `u` to `w`. It has the **single-chain Newton polytope property (SCNP)** if one chain already produces every
monomial. An, Tung and Zhang conjectured that some `D_u^w` fails SCNP exactly when `w` contains the pattern `4231`.

What is here:

* the **counterexample** `u = 214365`, `w = 541632`: `w` avoids `4231`, yet `D_u^w` has 81 monomials and no chain
  produces more than 73;
* the **explicit bottom** `u` that proves the "if" direction: every `w` containing `4231` has a `u` with `D_u^w` not SCNP;
* the complete classification for `n <= 11`: `w` admits a failure of SCNP exactly when it contains one of

      4231   541632   5276143   5472163
      54821A9376   5482A19376   64821A9375   6482A19375        (A = 10)

* the pictures of the competing chains (spirals and rosettes).

## Conventions

Permutations are written in one-line notation, positions and values start at 1, and `A`, `B`, ... stand for 10, 11, ...
A step of a chain swaps two **positions** `a < b`; it is allowed when the entry at `a` is smaller than the entry at `b`
and no entry strictly between the two positions has a value strictly between the two values. Its label is
`x_a + ... + x_(b-1)`. A *window* is an interval of positions `[p, q]`, and `g(p, q)` is the number of steps of a chain
with both positions inside the window. SCNP holds exactly when some chain has the smallest `g` in every window at once.

## Quick start (plain Python 3, nothing to install)

    cd python
    python3 counterexample.py      # everything about 214365 <= 541632, about 10 seconds
    python3 verify_small.py        # all small checks, about 2 minutes
    python3 spiral.py              # the cyclic orders and the missing rotations
    python3 boolean_intervals.py 4 5 6 7   # bad Boolean intervals and loops, a few seconds

`figures.py` needs `matplotlib`:

    python3 figures.py ../figures

## Layout

| path | content |
|---|---|
| `python/core.py` | the definitions, written to be read: steps, Bruhat order, chains, supports, the window criterion |
| `python/tokens.py` | inversions as tokens: follows every inversion along a chain and checks the escape formula |
| `python/construction.py` | the explicit bottom `u` for a top `w` containing `4231` |
| `python/counterexample.py` | the two chains `X`, `Y`, the six window bounds, and the full enumeration of 2688 chains |
| `python/spiral.py` | atoms, atom graph, competing chains as rotations of one cyclic order |
| `python/boolean_intervals.py` | every bad pair contains a bad Boolean pair, and every bad Boolean pair has a loop in its atom graph (`n <= 8`) |
| `python/figures.py` | the figures of the paper |
| `python/verify_small.py` | reruns every statement that fits in a few minutes |
| `c/sweep/` | fast C program: decides every pair of `S_n`, `n <= 10`, bottom by bottom (`make`, then `./scnp`) |
| `c/tops/` | second, independent C program: decides every pair below a given top, `n <= 11` |
| `data/sweep/` | output of the complete sweeps, `n <= 10` |
| `data/tops/` | output of the top-by-top runs, `n <= 11`, and of the partial run at `n = 12` |
| `data/boolean/` | output of `python/boolean_intervals.py`, `n <= 8` |
| `figures/` | the figures as PDF and PNG |
| `RESULTS.md` | all numbers, with the command that produces each of them |

## The two C programs

    cd c/sweep && make && ./scnp pair 214365 541632
    mkdir -p runs && ./scnp sweep -n 8 -o runs/s8 --list none && ./scnp patterns runs/s8

    cd c/tops && make -C fast && fast/scnp_fast top 541632 -v

Both use the same criterion (the least vector `g`), organised differently: `c/sweep` goes up from every bottom,
`c/tops` goes down from every top. `c/sweep` also contains a brute-force mode (`scnp check`, `scnp pair`) that expands
every chain and compares.

## Licence

MIT.
