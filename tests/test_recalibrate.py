"""Functional self-tests for control-based quality-score recalibration.

Covers the full HTTP flow (pair task -> recalibrate -> ZIP) plus the
all-or-nothing validation contract: missing/duplicate/unknown manifest ids,
out-of-bounds coordinates, corrupted FASTQ, missing results, alphabet and
quality limits. An independent reimplementation of the counting/scoring rules
is used to verify the outputs.

Run with the project virtualenv:
    .venv/bin/python -I tests/test_recalibrate.py
"""

from __future__ import annotations

import asyncio
import io
import math
import os
import shutil
import sys
import tempfile
import time
import uuid
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import httpx  # noqa: E402

from app import config, db, worker  # noqa: E402
from app.main import app  # noqa: E402

RESULTS: list[tuple[str, bool, str]] = []


def record(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" - {detail}" if detail else ""))


def make_record(rid: str, seq: bytes, quality: bytes,
                suffix: str = "", desc: bytes = b"", crlf: bool = False) -> bytes:
    header = b"@" + rid.encode() + suffix.encode()
    if desc:
        header += b" " + desc
    nl = b"\r\n" if crlf else b"\n"
    return header + nl + seq + nl + b"+" + nl + quality + nl


def make_fastq(rows: list[tuple[str, bytes, bytes]], side: int,
               crlf: bool = False) -> io.BytesIO:
    buf = io.BytesIO()
    for rid, seq, qual in rows:
        buf.write(make_record(rid, seq, qual, suffix=f"/{side}", crlf=crlf))
    buf.seek(0)
    return buf


def phred(score: int) -> bytes:
    return bytes((config.PHRED_OFFSET + score,))


def phreds(scores: list[int]) -> bytes:
    return bytes(config.PHRED_OFFSET + s for s in scores)


async def wait_terminal(client: httpx.AsyncClient, task_id: str,
                        timeout: float = 30.0) -> dict:
    deadline = time.time() + timeout
    while time.time() < deadline:
        body = (await client.get(f"/tasks/{task_id}")).json()
        if body["status"] in ("succeeded", "failed", "cancelled"):
            return body
        await asyncio.sleep(0.05)
    raise AssertionError(f"task {task_id} did not finish in {timeout}s")


async def upload_pair(client: httpx.AsyncClient,
                      r1: io.BytesIO, r2: io.BytesIO) -> str:
    r1.seek(0)
    r2.seek(0)
    r = await client.post(
        "/tasks",
        files={"r1": ("r1.fastq", r1, "application/octet-stream"),
               "r2": ("r2.fastq", r2, "application/octet-stream")},
    )
    assert r.status_code == 202, r.text
    body = await wait_terminal(client, r.json()["id"])
    assert body["status"] == "succeeded", body
    return body["id"]


async def recalibrate(client: httpx.AsyncClient, task_id: str,
                      reference: str, manifest: list[dict],
                      expected_status: int = 200):
    return await client.post(
        f"/tasks/{task_id}/recalibrate",
        json={"reference": reference, "manifest": manifest},
    )


def independent_stats(reference: bytes, pairs: list[dict]) -> tuple[int, int]:
    """Independent reimplementation of the mismatch-counting contract."""
    comp = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")
    total_sites = 0
    total_mm = 0
    for p in pairs:
        coverage: dict[int, tuple[int, int]] = {}
        conflicts: set[int] = set()
        for side in ("r1", "r2"):
            seq = p[side + "_seq"]
            qual = p[side + "_qual"]
            start = p[side + "_start"]
            if p[side + "_strand"] == "-":
                seq = seq.translate(comp)[::-1]
                qual = qual[::-1]
            for i, base in enumerate(seq):
                coord = start + i
                score = qual[i] - config.PHRED_OFFSET
                if coord in conflicts:
                    continue
                if coord in coverage:
                    pb, pq = coverage[coord]
                    if pb == base:
                        if score > pq:
                            coverage[coord] = (base, score)
                    else:
                        conflicts.add(coord)
                        del coverage[coord]
                else:
                    coverage[coord] = (base, score)
        for coord, (base, _q) in coverage.items():
            total_sites += 1
            if base != reference[coord - 1]:
                total_mm += 1
    return total_sites, total_mm


def independent_score(sites: int, mm: int):
    if sites < config.MIN_EFFECTIVE_SITES:
        return None
    raw = -10.0 * math.log10((mm + 1) / (sites + 2))
    return max(0, min(40, math.floor(raw + 0.5)))


def split_records(blob: bytes) -> list[tuple[bytes, bytes, bytes, bytes]]:
    lines = blob.splitlines(keepends=True)
    assert len(lines) % 4 == 0
    return [tuple(lines[i:i + 4]) for i in range(0, len(lines), 4)]


async def main_async() -> int:
    tmp_root = Path(tempfile.mkdtemp(prefix="recal-test-",
                                     dir=str(ROOT / "tests" / "output")))
    data_dir = tmp_root / "data"
    os.environ["FASTQ_DATA_DIR"] = str(data_dir)
    config.DATA_DIR = data_dir
    config.DB_PATH = data_dir / "tasks.db"
    config.TASKS_DIR = data_dir / "tasks"

    db._initialized = False
    db.init_db()
    worker.recover_on_startup()

    transport = httpx.ASGITransport(app=app)
    async with httpx.AsyncClient(transport=transport, base_url="http://test",
                                 timeout=60) as client:
        # ---- 1. happy path: forward reads, independent stats verification
        reference = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"  # 60 bp
        assert len(reference) == 60
        # R1 read 'a': start 1, len 40; R2 read 'a': start 21, len 40.
        # Overlap coordinates 21..40 (20 bp), identical bases, R2 higher qual.
        r1_seq = reference[0:40].encode()
        r2_seq = reference[20:60].encode()
        r1_qual = phreds([20] * 40)
        r2_qual = phreds([30] * 40)
        r1 = make_fastq([("a", r1_seq, r1_qual)], 1)
        r2 = make_fastq([("a", r2_seq, r2_qual)], 2)
        tid = await upload_pair(client, r1, r2)
        manifest = [dict(read_id="a", r1_start=1, r1_strand="+",
                         r2_start=21, r2_strand="+")]
        resp = await recalibrate(client, tid, reference, manifest)
        record("recalibrate returns 200 + zip",
               resp.status_code == 200
               and resp.headers["content-type"] == "application/zip",
               f"{resp.status_code} {resp.text[:200]}")

        zf = zipfile.ZipFile(io.BytesIO(resp.content))
        record("zip contains exactly R1.fastq and R2.fastq",
               sorted(zf.namelist()) == ["R1.fastq", "R2.fastq"],
               str(zf.namelist()))
        bad = zf.testzip()
        record("zip entries are not corrupt", bad is None, str(bad))

        pairs_meta = [dict(r1_seq=r1_seq, r1_qual=r1_qual, r1_start=1, r1_strand="+",
                           r2_seq=r2_seq, r2_qual=r2_qual, r2_start=21, r2_strand="+")]
        exp_sites, exp_mm = independent_stats(reference.encode(), pairs_meta)
        record("overlap merged: 60 effective sites", exp_sites == 60, str(exp_sites))
        record("no mismatches on perfect match", exp_mm == 0, str(exp_mm))
        exp_score = independent_score(exp_sites, exp_mm)
        record("score header matches independent computation",
               resp.headers["X-Recalibrated-Score"] == str(exp_score)
               and resp.headers["X-Effective-Sites"] == "60"
               and resp.headers["X-Mismatches"] == "0",
               dict(resp.headers))

        out_r1 = zf.read("R1.fastq")
        out_r2 = zf.read("R2.fastq")
        recs1 = split_records(out_r1)
        recs2 = split_records(out_r2)
        record("quality lines rewritten uniformly",
               recs1[0][3].strip() == phred(exp_score) * 40
               and recs2[0][3].strip() == phred(exp_score) * 40)
        # All other lines preserved byte-for-byte (compare original task output).
        orig_r1 = (await client.get(f"/tasks/{tid}/download/r1")).content
        orig_r2 = (await client.get(f"/tasks/{tid}/download/r2")).content
        o1 = split_records(orig_r1)[0]
        o2 = split_records(orig_r2)[0]
        record("non-quality lines preserved byte-for-byte",
               recs1[0][0:3] == o1[0:3] and recs2[0][0:3] == o2[0:3])
        record("original task files not modified",
               orig_r1 == (await client.get(f"/tasks/{tid}/download/r1")).content)

        # ---- 2. reverse strand: reverse-complement + reversed quality
        # R1 reverse at start 1 over the first 40 ref bases: its observed
        # bases (after RC) must equal reference[0:40], so its raw sequence is
        # RC(reference[0:40]).
        comp = bytes.maketrans(b"ACGT", b"TGCA")
        raw = reference[0:40].encode().translate(comp)[::-1]
        q_scores = [5 + (i % 30) for i in range(40)]
        rev_qual = phreds(q_scores)
        r1 = make_fastq([("rev", raw, rev_qual)], 1)
        r2 = make_fastq([("rev", r2_seq, r2_qual)], 2)
        tid2 = await upload_pair(client, r1, r2)
        manifest = [dict(read_id="rev", r1_start=1, r1_strand="-",
                         r2_start=21, r2_strand="+")]
        resp = await recalibrate(client, tid2, reference, manifest)
        record("reverse-strand request succeeds", resp.status_code == 200,
               resp.text[:200])
        pairs_meta = [dict(r1_seq=raw, r1_qual=rev_qual, r1_start=1, r1_strand="-",
                           r2_seq=r2_seq, r2_qual=r2_qual, r2_start=21, r2_strand="+")]
        exp_sites, exp_mm = independent_stats(reference.encode(), pairs_meta)
        record("reverse-strand counting matches independent stats",
               resp.headers["X-Effective-Sites"] == str(exp_sites)
               and resp.headers["X-Mismatches"] == str(exp_mm),
               f"{exp_sites}/{exp_mm}")

        # ---- 3. overlap conflict: same coordinate, different base -> excluded
        # R1 covers 1..40 matching ref; R2 covers 21..60 but at coord 30 has a
        # different base. Coord 30 must be excluded; all others counted.
        bad_r2 = bytearray(r2_seq)
        bad_r2[9] = ord("T") if bad_r2[9] != ord("T") else ord("A")  # coord 30
        bad_r2 = bytes(bad_r2)
        r1 = make_fastq([("c", r1_seq, r1_qual)], 1)
        r2 = make_fastq([("c", bad_r2, r2_qual)], 2)
        tid3 = await upload_pair(client, r1, r2)
        manifest = [dict(read_id="c", r1_start=1, r1_strand="+",
                         r2_start=21, r2_strand="+")]
        resp = await recalibrate(client, tid3, reference, manifest)
        record("conflicting overlap request succeeds", resp.status_code == 200,
               resp.text[:200])
        pairs_meta = [dict(r1_seq=r1_seq, r1_qual=r1_qual, r1_start=1, r1_strand="+",
                           r2_seq=bad_r2, r2_qual=r2_qual, r2_start=21, r2_strand="+")]
        exp_sites, exp_mm = independent_stats(reference.encode(), pairs_meta)
        record("conflicting coordinate excluded (59 sites)",
               resp.headers["X-Effective-Sites"] == str(exp_sites) == "59",
               resp.headers.get("X-Effective-Sites", ""))

        # ---- 4. mismatch drives the score down
        # Single pair covering 40 non-overlapping coords; 10 deliberate
        # mismatches against the reference.
        mm_seq = bytearray(reference[0:40].encode())
        flip = {ord("A"): ord("C"), ord("C"): ord("G"),
                ord("G"): ord("T"), ord("T"): ord("A")}
        for i in range(0, 40, 4):
            mm_seq[i] = flip[mm_seq[i]]
        mm_seq = bytes(mm_seq)
        r1 = make_fastq([("m", mm_seq, phreds([40] * 40))], 1)
        # R2 inside the reference: coords 21..60 with exact match.
        r2 = make_fastq([("m", reference[20:60].encode(), phreds([40] * 40))], 2)
        tid4 = await upload_pair(client, r1, r2)
        manifest = [dict(read_id="m", r1_start=1, r1_strand="+",
                         r2_start=21, r2_strand="+")]
        resp = await recalibrate(client, tid4, reference, manifest)
        pairs_meta = [dict(r1_seq=mm_seq, r1_qual=phreds([40] * 40),
                           r1_start=1, r1_strand="+",
                           r2_seq=reference[20:60].encode(),
                           r2_qual=phreds([40] * 40),
                           r2_start=21, r2_strand="+")]
        exp_sites, exp_mm = independent_stats(reference.encode(), pairs_meta)
        exp_score = independent_score(exp_sites, exp_mm)
        record("mismatch score matches independent calc",
               resp.status_code == 200
               and resp.headers["X-Recalibrated-Score"] == str(exp_score)
               and resp.headers["X-Mismatches"] == str(exp_mm),
               f"{resp.status_code} sites={exp_sites} mm={exp_mm} score={exp_score}")
        zf = zipfile.ZipFile(io.BytesIO(resp.content))
        recs = split_records(zf.read("R1.fastq"))
        record("R1 quality rewritten to lower score",
               recs[0][3].strip() == phred(exp_score) * 40)

        # ---- 5. fewer than 20 effective sites -> original scores retained
        # 40-bp mates covering the exact same coordinates with different bases
        # at every position: every coordinate is a conflict, so 0 sites count.
        short_ref = "ACGTACGTACGTACGTACGT"  # 20 bp
        long_seq = (short_ref * 2)[0:40].encode()
        r1 = make_fastq([("s", long_seq, phreds([7] * 40))], 1)
        flip = {ord("A"): ord("C"), ord("C"): ord("G"),
                ord("G"): ord("T"), ord("T"): ord("A")}
        conflict_r2 = bytearray(long_seq)
        for i in range(40):
            conflict_r2[i] = flip[conflict_r2[i]]
        r2 = make_fastq([("s", bytes(conflict_r2), phreds([9] * 40))], 2)
        ref40 = (short_ref * 2)
        tid5 = await upload_pair(client, r1, r2)
        manifest = [dict(read_id="s", r1_start=1, r1_strand="+",
                         r2_start=1, r2_strand="+")]
        resp = await recalibrate(client, tid5, ref40, manifest)
        record("all-conflict overlap -> 0 effective sites",
               resp.status_code == 200
               and resp.headers["X-Effective-Sites"] == "0"
               and resp.headers["X-Score-Applied"] == "0",
               resp.text[:200])
        zf = zipfile.ZipFile(io.BytesIO(resp.content))
        recs1 = split_records(zf.read("R1.fastq"))
        recs2 = split_records(zf.read("R2.fastq"))
        record("original quality retained when sites < 20",
               recs1[0][3].strip() == phreds([7] * 40)
               and recs2[0][3].strip() == phreds([9] * 40))
        record("retained output equals original bytes",
               zf.read("R1.fastq")
               == (await client.get(f"/tasks/{tid5}/download/r1")).content)

        # ---- 6. R1 order preserved and R2 aligned by read id, multi-pair
        ids = ["p1", "p2", "p3"]
        r1_rows = [(rid, reference[0:40].encode(), phreds([20] * 40))
                   for rid in ids]
        r2_rows = [(rid, reference[20:60].encode(), phreds([25] * 40))
                   for rid in reversed(ids)]
        r1 = make_fastq(r1_rows, 1)
        r2 = make_fastq(r2_rows, 2)
        tid6 = await upload_pair(client, r1, r2)
        manifest = [dict(read_id=rid, r1_start=1, r1_strand="+",
                         r2_start=21, r2_strand="+") for rid in ids]
        resp = await recalibrate(client, tid6, reference, manifest)
        record("multi-pair recalibration succeeds", resp.status_code == 200,
               resp.text[:200])
        zf = zipfile.ZipFile(io.BytesIO(resp.content))
        out1 = split_records(zf.read("R1.fastq"))
        out2 = split_records(zf.read("R2.fastq"))
        got_ids1 = [rec[0].split()[0] for rec in out1]
        got_ids2 = [rec[0].split()[0] for rec in out2]
        record("R1 order preserved in ZIP",
               got_ids1 == [b"@p1/1", b"@p2/1", b"@p3/1"], str(got_ids1))
        record("R2 aligned to R1 order in ZIP",
               got_ids2 == [b"@p1/2", b"@p2/2", b"@p3/2"], str(got_ids2))
        record("multi-pair site count header",
               resp.headers["X-Effective-Sites"] == "180",
               resp.headers.get("X-Effective-Sites", ""))

        # ---- 7. error: unknown task id
        r = await client.post(f"/tasks/{uuid.uuid4().hex}/recalibrate",
                              json={"reference": reference, "manifest": manifest})
        record("unknown task -> 404", r.status_code == 404, str(r.status_code))

        # ---- 8. error: task not successful
        bad_upload = io.BytesIO(b"@x\nACGT\n+\n~~\n")
        r = await client.post(
            "/tasks",
            files={"r1": ("r1.fastq", bad_upload, "application/octet-stream"),
                   "r2": ("r2.fastq", make_fastq(
                       [("x", b"A" * 40, phreds([10] * 40))], 2),
                       "application/octet-stream")},
        )
        failed_id = (await wait_terminal(client, r.json()["id"]))["id"]
        r = await recalibrate(client, failed_id, reference, manifest)
        record("non-successful task -> 409", r.status_code == 409,
               str(r.status_code))

        # ---- 9. error: manifest missing an id
        r = await recalibrate(client, tid6, reference, manifest[:2])
        record("missing manifest id -> 422", r.status_code == 422
               and "missing" in r.json()["detail"], r.text[:200])

        # ---- 10. error: duplicate manifest id
        dup_manifest = manifest + [manifest[0]]
        r = await recalibrate(client, tid6, reference, dup_manifest)
        record("duplicate manifest id -> 422", r.status_code == 422
               and "duplicate" in r.json()["detail"], r.text[:200])

        # ---- 11. error: unknown manifest id
        weird = [dict(read_id="ghost", r1_start=1, r1_strand="+",
                      r2_start=21, r2_strand="+")]
        r = await recalibrate(client, tid6, reference, weird)
        record("unknown manifest id -> 422", r.status_code == 422
               and "unknown" in r.json()["detail"], r.text[:200])

        # ---- 12. error: out-of-bounds coordinate
        oob = [dict(read_id=rid, r1_start=41, r1_strand="+",
                    r2_start=21, r2_strand="+") for rid in ids]
        r = await recalibrate(client, tid6, reference, oob)
        record("out-of-bounds start -> 422", r.status_code == 422
               and "out of bounds" in r.json()["detail"], r.text[:200])

        # ---- 13. error: reference too long / non-ACGT
        r = await recalibrate(client, tid6, "A" * 10001, manifest)
        record("oversized reference rejected by schema (422)",
               r.status_code == 422, str(r.status_code))
        r = await recalibrate(client, tid6, reference.replace("A", "N", 1),
                              manifest)
        record("non-ACGT reference -> 422", r.status_code == 422
               and "ACGT" in r.json()["detail"], r.text[:200])

        # ---- 14. error: strand value invalid
        bad_strand = [dict(read_id=rid, r1_start=1, r1_strand="*",
                           r2_start=21, r2_strand="+") for rid in ids]
        r = await recalibrate(client, tid6, reference, bad_strand)
        record("invalid strand -> 422", r.status_code == 422, str(r.status_code))

        # ---- 15. error: extra JSON field forbidden
        r = await client.post(
            f"/tasks/{tid6}/recalibrate",
            json={"reference": reference, "manifest": manifest, "extra": 1},
        )
        record("unknown JSON field -> 422", r.status_code == 422,
               str(r.status_code))

        # ---- 16. error: result file missing -> 409, whole request fails
        result_dir = config.TASKS_DIR / tid6 / "output"
        removed = result_dir / "R1.fastq"
        removed.unlink()
        r = await recalibrate(client, tid6, reference, manifest)
        record("missing result file -> 409, no ZIP",
               r.status_code == 409
               and r.headers.get("content-type") != "application/zip",
               str(r.status_code))

        # ---- 17. error: corrupted FASTQ (trailing junk / truncation) -> 422
        r1 = make_fastq([("j", reference[0:40].encode(), phreds([20] * 40))], 1)
        r2 = make_fastq([("j", reference[20:60].encode(), phreds([20] * 40))], 2)
        tid7 = await upload_pair(client, r1, r2)
        junk_target = config.TASKS_DIR / tid7 / "output" / "R1.fastq"
        with junk_target.open("ab") as fh:
            fh.write(b"@junk\nACGT\n+\n!!!!\n")
        r = await recalibrate(
            client, tid7, reference,
            [dict(read_id="j", r1_start=1, r1_strand="+",
                  r2_start=21, r2_strand="+")])
        record("extra trailing record rejected whole (422, no ZIP)",
               r.status_code == 422
               and r.headers.get("content-type") != "application/zip",
               r.text[:200])

        # A genuinely truncated trailing block must be reported as corruption.
        r1 = make_fastq([("k", reference[0:40].encode(), phreds([20] * 40))], 1)
        r2 = make_fastq([("k", reference[20:60].encode(), phreds([20] * 40))], 2)
        tid7b = await upload_pair(client, r1, r2)
        trunc_target = config.TASKS_DIR / tid7b / "output" / "R1.fastq"
        with trunc_target.open("ab") as fh:
            fh.write(b"@broken\nACGT\n+\n")
        r = await recalibrate(
            client, tid7b, reference,
            [dict(read_id="k", r1_start=1, r1_strand="+",
                  r2_start=21, r2_strand="+")])
        record("truncated trailing block -> 422 corruption, no ZIP",
               r.status_code == 422 and "corrupt" in r.json()["detail"],
               r.text[:200])

        # ---- 18. error: quality outside Phred 0..40
        bad_qual = b"J" + phreds([10] * 39)  # 'J' = score 41, invalid
        r1 = make_fastq([("q", reference[0:40].encode(), bad_qual)], 1)
        r2 = make_fastq([("q", reference[20:60].encode(), phreds([10] * 40))], 2)
        tid8 = await upload_pair(client, r1, r2)
        r = await recalibrate(
            client, tid8, reference,
            [dict(read_id="q", r1_start=1, r1_strand="+",
                  r2_start=21, r2_strand="+")])
        record("Phred score above 40 -> 422", r.status_code == 422
               and "Phred" in r.json()["detail"], r.text[:200])

        # ---- 19. CRLF line endings: other bytes preserved, quality rewritten
        r1 = make_fastq([("w", reference[0:40].encode(), phreds([15] * 40))],
                        1, crlf=True)
        r2 = make_fastq([("w", reference[20:60].encode(), phreds([15] * 40))],
                        2, crlf=True)
        tid9 = await upload_pair(client, r1, r2)
        r = await recalibrate(
            client, tid9, reference,
            [dict(read_id="w", r1_start=1, r1_strand="+",
                  r2_start=21, r2_strand="+")])
        record("CRLF input recalibrates", r.status_code == 200, r.text[:200])
        zf = zipfile.ZipFile(io.BytesIO(r.content))
        rec = split_records(zf.read("R1.fastq"))[0]
        record("CRLF line endings preserved",
               rec[0].endswith(b"\r\n") and rec[3].endswith(b"\r\n"))
        record("CRLF non-quality lines unchanged",
               rec[0:3] == split_records(
                   (await client.get(f"/tasks/{tid9}/download/r1")).content)[0][0:3])

        # ---- 20. score clamps at 40 and never goes below 0
        from app.recalibrate import recalibrated_score
        record("score clamps to 40", recalibrated_score(100000, 0) == 40)
        record("score floors at 0", recalibrated_score(100000, 100000) == 0)

    shutil.rmtree(tmp_root, ignore_errors=True)
    failed = [name for name, ok, _ in RESULTS if not ok]
    print("\n%d/%d checks passed" % (len(RESULTS) - len(failed), len(RESULTS)))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(asyncio.run(main_async()))
