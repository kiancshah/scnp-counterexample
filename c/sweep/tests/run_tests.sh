#!/usr/bin/env bash
# Regression tests: reproduces the n=6 example and cross-checks every mode.
set -euo pipefail
cd "$(dirname "$0")/.."
BIN=./scnp
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
pass() { printf '  ok   %s\n' "$1"; }
fail() { printf '  FAIL %s\n' "$1"; exit 1; }

echo "== n=6 example u=214365, w=541632"
out=$($BIN pair 214365 541632)
echo "$out" | sed 's/^/     /'
echo "$out" | grep -q '^chains  *2688 saturated' || fail "2688 chains"
echo "$out" | grep -q '^supp D_{u,w}  *81 monomials' || fail "support 81"
echo "$out" | grep -q '^best chain  *73 monomials' || fail "best single chain 73"
echo "$out" | grep -q '^SCNP  *FAILS.*brute force: fails' || fail "SCNP fails (both methods)"
pass "2688 chains, support 81, best single chain 73, SCNP fails"
[ "$(./scnp-generic pair 214365 541632)" = "$out" ] && pass "generic (non-AVX-512) build agrees" || fail "generic build"
if command -v python3 >/dev/null; then
  r=$(python3 ref/scnp_ref.py pair 214365 541632)
  [ "$r" = "chains=2688 supp=81 best_chain=73 scnp=NO" ] && pass "Python reference: $r" || fail "python ref: $r"
fi

echo "== smallest failure u=1324, w=4231"
out=$($BIN pair - <<< "1324 4231")
[ "$out" = "1324 4231 4 chains=24 supp=6 best_chain=5 scnp=NO" ] && pass "$out" || fail "$out"

echo "== brute-force cross-check of the fast criterion on all pairs"
for n in 4 5 6; do
  $BIN check -n $n --snp > "$TMP/check$n"
  grep -q 'agrees with brute force on every pair' "$TMP/check$n" || fail "check -n $n"
  pass "S_$n: $(awk '$1=="all"{print $2" pairs, "$3" SCNP failures, "$5" mismatches, "$6" SNP-uncertified"}' "$TMP/check$n")"
done

echo "== sweep: pair totals (OEIS A007767) and failure counts"
declare -A PAIRS=([4]=213 [5]=3781 [6]=98407 [7]=3550919)
declare -A FAILS=([4]=1 [5]=77 [6]=5132)
for n in 4 5 6 7; do
  $BIN sweep -n $n -o "$TMP/s$n" >/dev/null 2>&1
  $BIN sweep -n $n -o "$TMP/ns$n" --nosym >/dev/null 2>&1
  read -r p f < <(awk '$1=="all"{print $2, $3}' "$TMP/s$n/summary.txt")
  [ "$p" = "${PAIRS[$n]}" ] || fail "S_$n pairs $p"
  [ -z "${FAILS[$n]:-}" ] || [ "$f" = "${FAILS[$n]}" ] || fail "S_$n fails $f"
  cat "$TMP/s$n"/chunks/*.txt | sort > "$TMP/a$n"
  cat "$TMP/ns$n"/chunks/*.txt | sort > "$TMP/b$n"
  cmp -s "$TMP/a$n" "$TMP/b$n" || fail "S_$n symmetric vs unreduced lists differ"
  [ "$(wc -l < "$TMP/a$n")" = "$f" ] || fail "S_$n list length"
  pass "S_$n: $p pairs, $f failures; symmetry-reduced list == unreduced list"
done
if command -v python3 >/dev/null; then
  python3 ref/scnp_ref.py list 5 | sort > "$TMP/py5"
  cmp -s "$TMP/py5" "$TMP/a5" && pass "S_5 failure list == Python reference" || fail "python list 5"
fi
grep -q '^214365 541632 7$' "$TMP/a6" && pass "example pair is in the S_6 list" || fail "example missing"
[ "$(cat "$TMP/check6" | awk '$1=="all"{print $3}')" = "$(wc -l < "$TMP/a6")" ] && pass "check -n 6 and sweep -n 6 agree on the count" || fail "check vs sweep"

echo "== bad w (some u with D_u^w not SCNP) vs the patterns 4231, 541632, 5276143, 5472163"
declare -A BADW=([4]=1 [5]=17 [6]=208 [7]=2303 [8]=24905)
declare -A W4231=([4]=1 [5]=17 [6]=207 [7]=2278 [8]=24527)
declare -A AVOID=([4]=23 [5]=103 [6]=512 [7]=2737 [8]=15415)
$BIN sweep -n 8 -o "$TMP/s8" --list none >/dev/null 2>&1
for n in 4 5 6 7 8; do
  $BIN patterns "$TMP/s$n" > "$TMP/pat$n" || fail "S_$n: bad w != pattern containment"
  b=$(awk '/^bad w   /{print $3}' "$TMP/pat$n")
  c=$(awk '/^w containing 4231/{print $4}' "$TMP/pat$n")
  a=$(awk '/^w avoiding all/{print $5}' "$TMP/pat$n")
  [ "$b $c $a" = "${BADW[$n]} ${W4231[$n]} ${AVOID[$n]}" ] || fail "S_$n: $b $c $a"
  pass "S_$n: $b bad w, $c contain 4231, $a avoid all four; bad w == w containing a pattern"
done
cmp -s "$TMP/s7/badw.bin" "$TMP/ns7/badw.bin" && pass "S_7 bad-w set identical with and without symmetry reduction" || fail "badw sym"
[ "$(zcat ../../data/sweep/scnp_fails_S7.txt.gz | awk '{print $2}' | sort -u | wc -l)" = 2303 ] && pass "S_7 bad w == distinct tops in ../../data/sweep/scnp_fails_S7.txt.gz" || fail "data S7"

echo "== greedy chains (An-Tung-Zhang Def. 3.5)"
for nn in "4 210 147" "5 3567 1733"; do
  set -- $nn
  out=$($BIN greedy -n $1)
  s=$(echo "$out" | awk '/^  some greedy/{print $NF; exit}'); a=$(echo "$out" | awk '/^  every greedy/{print $NF; exit}')
  z=$(echo "$out" | awk '/no greedy chain at all/{print $NF}'); i=$(echo "$out" | awk '/^pairs with u = id/{p=$NF} /ATZ: all/{print ($NF==p)}')
  [ "$s $a $z $i" = "$2 $3 0 1" ] && pass "S_$1: some greedy dominant $s, every greedy dominant $a (Python reference: $2, $3); u = id: all greedy dominant" || fail "greedy $1: $s $a $z $i"
done
[ "$($BIN greedy -n 5 | awk '/^  neither/{print $NF}')" = 3 ] && pass "S_5: greedy from neither end dominant for exactly 3 avoider pairs" || fail "greedy both ends"
if command -v python3 >/dev/null; then
  [ "$(cd ref && python3 greedy_ref.py 4)" = "n=4: pairs 213, scnp 212, some greedy dominant 210, every greedy dominant 147, no greedy chain 0" ] && pass "S_4 greedy totals == Python reference" || fail "greedy python"
fi

echo "== output modes"
$BIN sweep -n 7 -o "$TMP/r7" --list reps --gzip >/dev/null 2>&1
zcat "$TMP/r7"/chunks/*.txt.gz | python3 tests/expand_reps.py | sort > "$TMP/e7"
cmp -s "$TMP/e7" "$TMP/a7" && pass "S_7 orbit representatives expand to the full list" || fail "reps expansion"
$BIN sweep -n 7 -o "$TMP/z7" --list none >/dev/null 2>&1
cmp -s <(tail -n +2 "$TMP/z7/summary.txt" | grep -v -e listed -e "cpu seconds") <(tail -n +2 "$TMP/s7/summary.txt" | grep -v -e listed -e "cpu seconds") && pass "--list none gives identical statistics" || fail "list none"

echo "== checkpoint / resume"
$BIN sweep -n 8 -o "$TMP/k8" --list none --chunks 512 >/dev/null 2>&1 &
pid=$!; sleep 0.7; kill -INT $pid 2>/dev/null || true; wait $pid || true
partial=$(ls "$TMP/k8/chunks/"*.stat 2>/dev/null | wc -l)
$BIN sweep -n 8 -o "$TMP/k8" --list none --chunks 512 >/dev/null 2>&1
$BIN sweep -n 8 -o "$TMP/f8" --list none --chunks 512 >/dev/null 2>&1
cmp -s <(grep -v cpu "$TMP/k8/summary.txt") <(grep -v cpu "$TMP/f8/summary.txt") && pass "interrupted after $partial/512 chunks, resumed run == fresh run" || fail "resume"
if $BIN sweep -n 8 -o "$TMP/k8" --list none --chunks 256 >/dev/null 2>&1; then fail "config mismatch not detected"; else pass "resume with different parameters is refused"; fi

echo "all tests passed"
