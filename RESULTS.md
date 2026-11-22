# Results

"Bad pair" means that `D_u^w` fails SCNP. "Bad `w`" means that some `u <= w` gives a bad pair.

## Counts

| n | pairs u <= w | bad pairs | bad w | w containing 4231 | bad w avoiding 4231 |
|---|---:|---:|---:|---:|---:|
| 4 | 213 | 1 | 1 | 1 | 0 |
| 5 | 3 781 | 77 | 17 | 17 | 0 |
| 6 | 98 407 | 5 132 | 208 | 207 | 1 |
| 7 | 3 550 919 | 364 214 | 2 303 | 2 278 | 25 |
| 8 | 170 288 585 | 29 131 995 | 24 905 | 24 527 | 378 |
| 9 | 10 501 351 657 | 2 671 695 123 | 272 641 | 268 104 | 4 537 |
| 10 | 811 107 567 659 | 281 859 859 053 | 3 084 713 | 3 036 850 | 47 863 |
| 11 | | | 36 559 185 | 36 092 688 | 466 497 |

* `n <= 10`: every pair decided. `cd c/sweep && ./scnp sweep -n N -o runs/sN --list none && ./scnp patterns runs/sN`.
  Outputs in `data/sweep/summary_S*.txt`, `data/sweep/patterns_S*.txt`, `data/sweep/badw_S10.bin.gz`.
  Time on 4 cores: 2 s for `n = 8`, 1.5 min for `n = 9`, 2.5 h for `n = 10`.
* `n = 11`: every pair whose top avoids 4231 decided (1 129 772 239 859 pairs, 3 461 289 736 bad), with `c/tops`.
  Tops containing 4231 are bad by the construction in `python/construction.py`. Outputs in `data/tops/`.
* `n <= 6`: also by plain Python, `python3 python/verify_small.py`.

## The eight minimal bad permutations

For `n <= 11`, `w` is bad if and only if it contains one of these. Each avoids the other seven.

| w | a shortest bad u | length | chains | monomials of D | best single chain |
|---|---|---:|---:|---:|---:|
| 4231 | 1324 | 4 | 24 | 6 | 5 |
| 541632 | 214365 | 7 | 2 688 | 81 | 73 |
| 5276143 | 1432765 | 7 | 5 040 | 135 | 126 |
| 5472163 | 3216547 | 7 | 5 040 | 135 | 126 |
| 54821A9376 | 3217564A98 | 10 | 3 628 800 | 2 800 | 2 674 |
| 5482A19376 | 3217654A98 | 10 | 3 628 800 | 5 916 | 5 664 |
| 64821A9375 | 3217654A98 | 10 | 3 628 800 | 2 800 | 2 674 |
| 6482A19375 | 3217654A98 | 11 | 68 428 800 | 14 136 | 13 744 |

Each line: `cd c/sweep && ./scnp pair U W` (brute force over all chains; the lines of size 10 take a few seconds each, the last one longer).

## Further facts

* No bad pair has both `u` avoiding 1324 and `w` avoiding 4231 (`n <= 11`; 9 927 025 944 such pairs in `S_11`).
* Every bad pair has length at least 4 (`n <= 10`); every bad pair whose top avoids 4231 has length at least 7 (`n <= 11`).
* `u = id` never gives a bad pair (`n <= 10`), as proved by An, Tung and Zhang.
* The explicit bottom of `python/construction.py` gives a bad pair for every `w` containing 4231, `n <= 8`.
* If `Inv(u)` is contained in `Inv(w)`, the pair is good and the least vector is `inv_W(w) - inv_W(u)` (`n <= 7`).

## Beyond n = 11 (partial)

* `n = 12`: every `w` containing one of the eight is bad (all 457 889 412 of them); every `w` avoiding the eight with
  at most 30 inversions is good (18 208 989 permutations). `data/tops/n12_if_direction.json`, `data/tops/n11_progress.json`.
* `n = 13`: for every bottom `u` that is a direct sum of three decreasing blocks of sizes at least 2, all tops were
  decided (`./scnp bigbottom U --prev data/sweep/badw_S10.bin.gz --tmp DIR`): 19 657 688 320 pairs, and no bad top
  avoids the eight.

## Independent checks

* The window criterion against expansion of every chain: all pairs of `S_n`, `n <= 6` (`./scnp check -n 6`), and all
  pairs of `S_4`, `S_5` in Python.
* `c/sweep` and `c/tops` agree on every `w` in `S_10` that avoids 4231 (591 950 permutations).
* A further independent implementation reproduced the table above for `n <= 8`, the numbers for `214365 <= 541632`,
  and the badness and minimality of the four permutations of size 10.
