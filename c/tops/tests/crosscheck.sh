#!/usr/bin/env bash
# Cross-checks of the fast checker (fast/scnp_fast) against
#   (1) the brute-force Python reference (polynomial expansion), all pairs n <= 5,
#       plus random pairs of S_6 and S_7,
#   (2) the readable Python criterion, all pairs n <= 6,
#   (3) the failing-pair lists data/scnp_fails_S{4..7}.txt.gz of the main repo,
#   (4) the per-length summaries data/summary_S{6..8}.txt of the main repo,
#   (5) its own second code path (-T: steps generated on the fly, no table).
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
FAST="$HERE/fast/scnp_fast"
REF="$HERE/reference"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
pass() { printf '  ok   %s\n' "$1"; }
fail() { printf '  FAIL %s\n' "$1"; exit 1; }

make -s -C "$HERE/fast" >/dev/null

allperms() { python3 -c "
import itertools,sys
n=int(sys.argv[1])
for p in itertools.permutations(range(1,n+1)): print(''.join('123456789AB'[x-1] for x in p))" "$1"; }

echo "== (1) brute-force reference vs fast checker, every pair of S_4, S_5"
for n in 4 5; do
  allperms $n > "$TMP/p$n"
  "$FAST" list "$TMP/p$n" "$TMP/f$n" -p
  sort "$TMP/f$n.pairs" > "$TMP/fast$n"
  python3 - "$n" > "$TMP/ref$n" <<EOF
import sys, itertools
sys.path.insert(0, "$REF")
from scnp_reference import up_set, is_scnp, length, format_perm
n = int(sys.argv[1])
for u in itertools.permutations(range(1, n + 1)):
    for w in up_set(u):
        if not is_scnp(u, w):
            print(format_perm(u), format_perm(w), length(w) - length(u))
EOF
  sort -o "$TMP/ref$n" "$TMP/ref$n"
  cmp -s "$TMP/fast$n" "$TMP/ref$n" || fail "S_$n bad-pair lists differ"
  pass "S_$n: identical bad-pair lists ($(wc -l < "$TMP/ref$n") pairs)"
done

echo "== (1b) brute-force reference vs fast checker, random pairs of S_6 and S_7"
python3 - > "$TMP/rand" <<EOF
import random, sys
sys.path.insert(0, "$REF")
from scnp_reference import up_set, format_perm, length
random.seed(20261103)
for n, count in ((6, 60), (7, 25)):
    done = 0
    while done < count:
        u = tuple(random.sample(range(1, n + 1), n))
        ups = sorted(up_set(u))
        w = random.choice(ups)
        if length(w) - length(u) > 9:      # keep brute force affordable
            continue
        print(format_perm(u), format_perm(w))
        done += 1
EOF
nr=0
while read -r u w; do
  f=$("$FAST" pair "$u" "$w" | sed 's/.*scnp=//')
  r=$(python3 -c "
import sys; sys.path.insert(0, '$REF')
from scnp_reference import parse_perm, is_scnp
print('yes' if is_scnp(parse_perm('$u'), parse_perm('$w')) else 'NO')")
  [ "$f" = "$r" ] || fail "pair $u $w: fast=$f reference=$r"
  nr=$((nr+1))
done < "$TMP/rand"
pass "$nr random pairs agree (S_6, S_7, length <= 9)"

echo "== (2) readable Python criterion vs fast checker, every pair of S_6"
allperms 6 > "$TMP/p6"
"$FAST" list "$TMP/p6" "$TMP/f6" -p
sort "$TMP/f6.pairs" > "$TMP/fast6"
python3 - > "$TMP/crit6" <<EOF
import sys, itertools
sys.path.insert(0, "$REF")
from scnp_criterion import sweep_down
from scnp_reference import format_perm, length
for w in itertools.permutations(range(1, 7)):
    _, ok = sweep_down(w)
    for v, good in ok.items():
        if not good:
            print(format_perm(v), format_perm(w), length(w) - length(v))
EOF
sort -o "$TMP/crit6" "$TMP/crit6"
cmp -s "$TMP/fast6" "$TMP/crit6" || fail "S_6 criterion lists differ"
pass "S_6: identical bad-pair lists ($(wc -l < "$TMP/crit6") pairs)"

echo "== (3) main-repo failing-pair lists"
for n in 4 5 6 7; do
  [ -f "$TMP/p$n" ] || allperms $n > "$TMP/p$n"
  [ -f "$TMP/f$n.pairs" ] || "$FAST" list "$TMP/p$n" "$TMP/f$n" -p
  sort "$TMP/f$n.pairs" > "$TMP/fast$n"
  zcat "$REPO/data/sweep/scnp_fails_S$n.txt.gz" | sort > "$TMP/repo$n"
  cmp -s "$TMP/fast$n" "$TMP/repo$n" || fail "S_$n differs from data/scnp_fails_S$n"
  pass "S_$n: identical to data/scnp_fails_S$n.txt.gz ($(wc -l < "$TMP/repo$n") pairs)"
done

echo "== (4) main-repo per-length summaries; (5) table vs on-the-fly steps"
for n in 6 7 8; do
  "$FAST" all $n "$TMP/a$n"
  "$FAST" -T all $n "$TMP/b$n"
  awk '$1 ~ /^[0-9]+$/ {print $1, $2, $3}' "$TMP/a$n" > "$TMP/mine$n"
  awk '$1 ~ /^[0-9]+$/ {print $1, $2, $3}' "$REPO/data/sweep/summary_S$n.txt" > "$TMP/theirs$n"
  cmp -s "$TMP/mine$n" "$TMP/theirs$n" || fail "S_$n per-length counts differ"
  cmp -s "$TMP/a$n.bits" "$TMP/b$n.bits" || fail "S_$n -T bitmap differs"
  pass "S_$n: per-length pairs and bad pairs equal data/summary_S$n.txt; -T path identical"
done
echo "all cross-checks passed"
