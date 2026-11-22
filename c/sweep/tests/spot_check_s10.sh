#!/usr/bin/env bash
# Brute-force spot check in S_10: (1) random pairs from random upward cover
# walks of length 4..9, (2) a sample of pairs the fast criterion rejects.
# Every line must be confirmed by explicit expansion (no BRUTE_DISAGREES).
set -euo pipefail
cd "$(dirname "$0")/.."
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
python3 - > "$TMP/rand" <<'PY'
import random
random.seed(12345)
n = 10
def covers(v):
    out = []
    for a in range(n):
        lim = n + 1
        for b in range(a + 1, n):
            if v[a] < v[b] < lim:
                lim = v[b]; x = list(v); x[a], x[b] = x[b], x[a]; out.append(tuple(x))
    return out
fmt = lambda p: ''.join('A' if x == 10 else str(x) for x in p)
for _ in range(400):
    u = list(range(1, n + 1)); random.shuffle(u); w = u = tuple(u)
    for _ in range(random.randint(4, 9)):
        c = covers(w)
        if not c: break
        w = random.choice(c)
    print(fmt(u), fmt(w))
PY
for u in 21436587A9 1324576A98 3127A45698 1A23456789; do
  ./scnp bottom $u --maxlen 8 2>/dev/null
done | awk 'NR % 97 == 1' | head -150 > "$TMP/fails"
./scnp pair - < "$TMP/rand" > "$TMP/rand.out"
./scnp pair - < "$TMP/fails" > "$TMP/fails.out"
bad=$(cat "$TMP/rand.out" "$TMP/fails.out" | grep -cE 'DISAGREE|INCONSISTENT|supp=\?' || true)
echo "random pairs: $(wc -l < "$TMP/rand.out"), of which SCNP fails: $(grep -c 'scnp=NO' "$TMP/rand.out" || true)"
echo "rejected pairs re-checked: $(wc -l < "$TMP/fails.out"), confirmed failing by brute force: $(grep -c 'scnp=NO' "$TMP/fails.out" || true)"
echo "disagreements: $bad"
[ "$bad" = 0 ]
