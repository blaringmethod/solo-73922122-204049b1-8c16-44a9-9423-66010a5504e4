#!/usr/bin/env bash
# End-to-end demo against a running service (default http://localhost:8000).
# Shows: upload -> progress polling -> download verification, plus cancellation
# and a rejected (malformed) upload.
set -euo pipefail

BASE_URL="${BASE_URL:-http://localhost:8000}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

need() { command -v "$1" >/dev/null || { echo "missing tool: $1" >&2; exit 1; }; }
need curl
need python3

echo "== health =="
curl -fsS "$BASE_URL/health"
echo

mk() { # file id seq side
  # 'I' is Phred+33 score 40, the top of the accepted 0-40 quality window.
  printf '@%s/%s desc\n%s\n+\n%s\n' "$2" "$4" "$3" "$(printf 'I%.0s' $(seq ${#3}))" >> "$1"
}

echo "== generating shuffled R1/R2 (40-bp ACGT reads) =="
R1="$WORK/R1.fastq"; R2="$WORK/R2.fastq"
SEQ1="$(python3 -c 'print("ACGT" * 10, end="")')"   # 40 bp
SEQ2="$(python3 -c 'print("TGCA" * 10, end="")')"   # 40 bp
for id in alpha beta gamma delta epsilon; do
  mk "$R1" "$id" "$SEQ1" 1
done
for id in epsilon gamma alpha delta beta; do   # deliberately shuffled
  mk "$R2" "$id" "$SEQ2" 2
done

echo "== upload =="
TASK=$(curl -fsS -X POST "$BASE_URL/tasks" \
  -F "r1=@$R1;filename=R1.fastq" \
  -F "r2=@$R2;filename=R2.fastq")
echo "$TASK"
ID=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])' <<<"$TASK")

echo "== poll until terminal =="
while :; do
  BODY=$(curl -fsS "$BASE_URL/tasks/$ID")
  STATUS=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["status"])' <<<"$BODY")
  echo "$BODY"
  case "$STATUS" in
    succeeded|failed|cancelled) break ;;
  esac
  sleep 0.5
done

echo "== download and verify order =="
curl -fsS "$BASE_URL/tasks/$ID/download/r1" -o "$WORK/out_R1.fastq"
curl -fsS "$BASE_URL/tasks/$ID/download/r2" -o "$WORK/out_R2.fastq"
python3 - "$WORK/out_R1.fastq" "$WORK/out_R2.fastq" <<'PY'
import sys
ids = []
for path in sys.argv[1:]:
    got = []
    with open(path, "rb") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            got.append(header[1:].split()[0].decode().rsplit("/", 1)[0])
            fh.readline(); fh.readline(); fh.readline()
    ids.append(got)
print("R1 ids:", ids[0])
print("R2 ids:", ids[1])
assert ids[0] == ids[1] == ["alpha", "beta", "gamma", "delta", "epsilon"]
print("OK: paired in R1 order")
PY

echo "== recalibrate the successful task against a reference =="
python3 - "$WORK/manifest.json" <<'PY'
import json, sys
# 60-bp ACGT reference. R1 (40 bp, forward) covers 1..40; R2 (40 bp,
# reverse) covers 21..60. The manifest ids match the uploaded read ids.
ref = ("ACGTACGTAC" * 6)
manifest = {
    "reference": ref,
    "manifest": [
        {"read_id": rid, "r1_start": 1, "r1_strand": "+",
         "r2_start": 21, "r2_strand": "-"}
        for rid in ("alpha", "beta", "gamma", "delta", "epsilon")
    ],
}
with open(sys.argv[1], "w") as fh:
    json.dump(manifest, fh)
PY
curl -fsS -D "$WORK/recal_headers.txt" -X POST "$BASE_URL/tasks/$ID/recalibrate" \
  -H 'Content-Type: application/json' \
  --data @"$WORK/manifest.json" -o "$WORK/recal.zip"
grep -iE 'HTTP/|x-effective|x-mismatch|x-recalibrated|x-score-applied|content-type' \
  "$WORK/recal_headers.txt"
python3 - "$WORK/recal.zip" <<'PY'
import sys, zipfile
zf = zipfile.ZipFile(sys.argv[1])
assert zf.testzip() is None and sorted(zf.namelist()) == ["R1.fastq", "R2.fastq"]
r1 = zf.read("R1.fastq").decode().splitlines()
print("recalibrated R1 quality line:", r1[3])
print("OK: ZIP contains rewritten R1.fastq and R2.fastq")
PY

echo "== malformed input is rejected as a whole task =="
printf '@x\nACGT\n+\n~~\n' > "$WORK/bad_R2.fastq"     # seq 4 != qual 2
BAD=$(curl -fsS -X POST "$BASE_URL/tasks" \
  -F "r1=@$R1;filename=R1.fastq" \
  -F "r2=@$WORK/bad_R2.fastq;filename=R2.fastq")
BAD_ID=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])' <<<"$BAD")
sleep 1
curl -fsS "$BASE_URL/tasks/$BAD_ID"
echo
echo "(download attempt on failed task -> HTTP 409)"
curl -s -o /dev/null -w '%{http_code}\n' "$BASE_URL/tasks/$BAD_ID/download/r1"

echo "== cancellation returns immediately, task never publishes =="
BIG_R1="$WORK/big_R1.fastq"; BIG_R2="$WORK/big_R2.fastq"
python3 - "$BIG_R1" "$BIG_R2" <<'PY'
import sys
n = 200_000
with open(sys.argv[1], "wb") as f1, open(sys.argv[2], "wb") as f2:
    for i in range(n):
        f1.write(b"@r%09d/1\n" % i + b"A" * 40 + b"\n+\n" + b"~" * 40 + b"\n")
    for i in range(n - 1, -1, -1):
        f2.write(b"@r%09d/2\n" % i + b"T" * 40 + b"\n+\n" + b"!" * 40 + b"\n")
PY
CANCEL=$(curl -fsS -X POST "$BASE_URL/tasks" \
  -F "r1=@$BIG_R1;filename=big_R1.fastq" \
  -F "r2=@$BIG_R2;filename=big_R2.fastq")
CANCEL_ID=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])' <<<"$CANCEL")
curl -fsS -X POST "$BASE_URL/tasks/$CANCEL_ID/cancel"
echo
sleep 2
curl -fsS "$BASE_URL/tasks/$CANCEL_ID"
echo

echo "== cancelling the already succeeded task -> HTTP 409 =="
curl -s -o /dev/null -w '%{http_code}\n' -X POST "$BASE_URL/tasks/$ID/cancel"

echo "ALL DEMO STEPS DONE"
