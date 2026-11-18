"""Control-based quality-score recalibration.

Given the paired FASTQ of an already successful pairing task, an ACGT
reference and a manifest of 1-based start positions / strand orientations,
this module:

1. validates the inputs as an all-or-nothing contract (any violation raises
   :class:`RecalError`; no partial output is ever produced);
2. counts reference mismatches using the *original* quality scores, merging
   the two mates where they cover the same reference coordinate;
3. computes one recalibrated Phred score from
   ``-10*log10((mismatches+1)/(sites+2))`` once at least
   ``MIN_EFFECTIVE_SITES`` effective sites exist, otherwise keeps the
   original scores;
4. rewrites only the quality lines, preserving every other original byte and
   the R1 ordering, and packages both FASTQ files into an in-memory ZIP.

The whole input is bounded (max 2000 pairs * 200 bases), so the in-memory
working set is a few megabytes at most and stays well inside the 256 MiB
process cap.
"""

from __future__ import annotations

import io
import math
import zipfile
from dataclasses import dataclass
from pathlib import Path

from . import config

# Sentinel stored in the per-pair coverage map when both mates cover a
# coordinate with different bases: such a coordinate is never counted.
_CONFLICT = object()

_COMPLEMENT = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")
_ACGT = frozenset(b"ACGT")


class RecalError(ValueError):
    """A validation or corruption error; the whole request must fail."""


@dataclass(frozen=True)
class ManifestEntry:
    """Manifest row for one read pair."""

    read_id: str
    r1_start: int
    r1_strand: str  # '+' forward, '-' reverse
    r2_start: int
    r2_strand: str


@dataclass
class FastqRecord:
    """One parsed FASTQ record with its exact four original lines."""

    read_id: str
    header: bytes
    sequence: bytes          # payload without line ending
    plus: bytes              # full line including line ending
    quality: bytes           # payload without line ending
    seq_line_end: bytes
    qual_line_end: bytes

    def original_bytes(self) -> bytes:
        return (
            self.header
            + self.sequence
            + self.seq_line_end
            + self.plus
            + self.quality
            + self.qual_line_end
        )

    def rewritten_bytes(self, qual_payload: bytes) -> bytes:
        return (
            self.header
            + self.sequence
            + self.seq_line_end
            + self.plus
            + qual_payload
            + self.qual_line_end
        )


# ---------------------------------------------------------------------------
# Parsing and structural validation
# ---------------------------------------------------------------------------

def _normalize_read_id(header_payload: bytes, side_label: str, record_no: int) -> str:
    if not header_payload.startswith(b"@"):
        raise RecalError(
            f"{side_label}: record {record_no} header does not start with '@' "
            "(corrupt FASTQ)"
        )
    first_field = (
        header_payload[1:].split(None, 1)[0] if len(header_payload) > 1 else b""
    )
    if not first_field:
        raise RecalError(
            f"{side_label}: record {record_no} has an empty read id (corrupt FASTQ)"
        )
    if (
        len(first_field) >= 2
        and first_field[-2:-1] == b"/"
        and first_field[-1:] in (b"1", b"2")
    ):
        first_field = first_field[:-2]
    if not first_field:
        raise RecalError(
            f"{side_label}: record {record_no} has an empty read id after "
            "stripping '/1' or '/2' (corrupt FASTQ)"
        )
    return first_field.decode("latin-1")


def parse_fastq_file(path: Path, side: int) -> list[FastqRecord]:
    """Parse a whole FASTQ file, enforcing the recalibration contract.

    Returns records in file order. Raises :class:`RecalError` on any format
    violation, length/alphabet/quality breach, or duplicate read id.
    """
    side_label = f"R{side}"
    if not path.is_file():
        raise RecalError(f"{side_label}: result file missing on server")

    records: list[FastqRecord] = []
    seen: set[str] = set()
    record_no = 0
    with path.open("rb") as fh:
        while True:
            lines = [fh.readline(), fh.readline(), fh.readline(), fh.readline()]
            if not lines[0]:
                if any(lines[1:]):
                    raise RecalError(
                        f"{side_label}: truncated final record, expected four "
                        "lines (trailing corruption)"
                    )
                break
            record_no += 1
            header, seq_line, plus_line, qual_line = lines
            if any(not line for line in lines):
                raise RecalError(
                    f"{side_label}: record {record_no} is missing lines, FASTQ "
                    "records must span exactly four lines (corrupt FASTQ)"
                )
            if not header.startswith(b"@"):
                raise RecalError(
                    f"{side_label}: record {record_no} header does not start "
                    "with '@' (corrupt FASTQ)"
                )
            if not plus_line.startswith(b"+"):
                raise RecalError(
                    f"{side_label}: record {record_no} third line does not "
                    "start with '+' (corrupt FASTQ)"
                )

            if seq_line.endswith(b"\r\n"):
                seq_end = b"\r\n"
            elif seq_line.endswith(b"\n"):
                seq_end = b"\n"
            else:
                raise RecalError(
                    f"{side_label}: record {record_no} sequence line has no "
                    "line ending (trailing corruption)"
                )
            if qual_line.endswith(b"\r\n"):
                qual_end = b"\r\n"
            elif qual_line.endswith(b"\n"):
                qual_end = b"\n"
            else:
                raise RecalError(
                    f"{side_label}: record {record_no} quality line has no "
                    "line ending (trailing corruption)"
                )

            sequence = seq_line[: -len(seq_end)]
            quality = qual_line[: -len(qual_end)]
            if not sequence:
                raise RecalError(
                    f"{side_label}: record {record_no} has an empty sequence"
                )
            if not quality:
                raise RecalError(
                    f"{side_label}: record {record_no} has an empty quality string"
                )
            if len(sequence) != len(quality):
                raise RecalError(
                    f"{side_label}: record {record_no} sequence length "
                    f"{len(sequence)} != quality length {len(quality)}"
                )
            read_len = len(sequence)
            if not config.MIN_READ_LEN <= read_len <= config.MAX_READ_LEN:
                raise RecalError(
                    f"{side_label}: record {record_no} length {read_len} is "
                    f"outside the allowed {config.MIN_READ_LEN}-"
                    f"{config.MAX_READ_LEN} base window"
                )
            non_acgt = next((b for b in sequence if b not in _ACGT), None)
            if non_acgt is not None:
                raise RecalError(
                    f"{side_label}: record {record_no} contains non-ACGT base "
                    f"{chr(non_acgt)!r}; only ACGT is accepted"
                )
            bad_q = next(
                (
                    q
                    for q in quality
                    if not config.PHRED_OFFSET
                    <= q
                    <= config.PHRED_OFFSET + config.MAX_PHRED
                ),
                None,
            )
            if bad_q is not None:
                raise RecalError(
                    f"{side_label}: record {record_no} quality byte "
                    f"{chr(bad_q)!r} is outside the Phred+33 range "
                    f"{config.MIN_PHRED}-{config.MAX_PHRED}"
                )

            read_id = _normalize_read_id(
                header.rstrip(b"\r\n"), side_label, record_no
            )
            if read_id in seen:
                raise RecalError(
                    f"{side_label}: duplicate read id '{read_id}' at record "
                    f"{record_no}"
                )
            seen.add(read_id)
            records.append(
                FastqRecord(
                    read_id=read_id,
                    header=header,
                    sequence=sequence,
                    plus=plus_line,
                    quality=quality,
                    seq_line_end=seq_end,
                    qual_line_end=qual_end,
                )
            )

        # Anything after the last complete four-line block is corruption.
        tail = fh.read()
        if tail:
            raise RecalError(
                f"{side_label}: trailing bytes after the last record "
                "(trailing corruption)"
            )

    if record_no == 0:
        raise RecalError(f"{side_label}: file contains no records")
    if record_no > config.MAX_READ_PAIRS:
        raise RecalError(
            f"{side_label}: {record_no} records exceed the limit of "
            f"{config.MAX_READ_PAIRS} read pairs"
        )
    return records


def _index_by_id(records: list[FastqRecord]) -> dict[str, FastqRecord]:
    return {rec.read_id: rec for rec in records}


# ---------------------------------------------------------------------------
# Manifest and reference validation
# ---------------------------------------------------------------------------

def validate_reference(reference: str) -> bytes:
    if not isinstance(reference, str) or not reference:
        raise RecalError("reference must be a non-empty ACGT string")
    ref = reference.encode("ascii")
    if len(ref) > config.MAX_REF_BASES:
        raise RecalError(
            f"reference length {len(ref)} exceeds the "
            f"{config.MAX_REF_BASES}-base limit"
        )
    bad = next((b for b in ref if b not in _ACGT), None)
    if bad is not None:
        raise RecalError(
            f"reference contains non-ACGT base {chr(bad)!r}; only ACGT is "
            "accepted"
        )
    return ref


def validate_manifest(
    manifest: list[ManifestEntry],
    r1_ids: set[str],
    r2_ids: set[str],
    reference_len: int,
    r1_by_id: dict[str, FastqRecord],
    r2_by_id: dict[str, FastqRecord],
) -> None:
    if not manifest:
        raise RecalError("manifest must contain an entry for every read pair")

    seen: set[str] = set()
    for entry in manifest:
        if entry.read_id in seen:
            raise RecalError(
                f"manifest: duplicate read id '{entry.read_id}'"
            )
        seen.add(entry.read_id)
        if entry.read_id not in r1_ids:
            raise RecalError(
                f"manifest: unknown read id '{entry.read_id}' (not present "
                "in R1)"
            )
        if entry.read_id not in r2_ids:
            raise RecalError(
                f"manifest: unknown read id '{entry.read_id}' (not present "
                "in R2)"
            )
        if entry.r1_strand not in ("+", "-") or entry.r2_strand not in ("+", "-"):
            raise RecalError(
                f"manifest: strand for '{entry.read_id}' must be '+' or '-'"
            )
        for side, start in (
            (1, entry.r1_start),
            (2, entry.r2_start),
        ):
            if start < 1:
                raise RecalError(
                    f"manifest: R{side} start for '{entry.read_id}' must be a "
                    "1-based position (>= 1)"
                )
            rec = r1_by_id[entry.read_id] if side == 1 else r2_by_id[entry.read_id]
            # 1-based inclusive window. For reverse-strand mates the start is
            # the leftmost (lowest) reference coordinate the read covers.
            end = start + len(rec.sequence) - 1
            if end > reference_len:
                raise RecalError(
                    f"manifest: R{side} of '{entry.read_id}' covers "
                    f"coordinates {start}-{end} but the reference is only "
                    f"{reference_len} bases long (out of bounds)"
                )

    missing = r1_ids - seen
    if missing:
        raise RecalError(
            f"manifest: missing entry for read id "
            f"'{next(iter(sorted(missing)))}'"
        )
    extra_r2 = r2_ids - r1_ids
    if extra_r2:
        # Pairing should have prevented this, but fail the whole request
        # rather than silently dropping a mate.
        raise RecalError(
            f"R2 contains read id without an R1 mate: "
            f"'{next(iter(sorted(extra_r2)))}'"
        )
    if len(manifest) != len(r1_ids):
        raise RecalError(
            f"manifest has {len(manifest)} entries but the FASTQ holds "
            f"{len(r1_ids)} read pairs"
        )


# ---------------------------------------------------------------------------
# Mismatch statistics
# ---------------------------------------------------------------------------

def _oriented_bases(
    record: FastqRecord, start: int, strand: str
) -> list[tuple[int, int, int]]:
    """Return ``(reference coordinate, base byte, original Phred score)``.

    Coordinates are 1-based. Reverse-strand reads are reverse-complemented and
    their quality string is reversed before comparison.
    """
    if strand == "+":
        seq = record.sequence
        qual = record.quality
    else:
        seq = record.sequence.translate(_COMPLEMENT)[::-1]
        qual = record.quality[::-1]
    return [
        (
            start + i,
            seq[i],
            qual[i] - config.PHRED_OFFSET,
        )
        for i in range(len(seq))
    ]


def count_mismatches(
    manifest: list[ManifestEntry],
    reference: bytes,
    r1_by_id: dict[str, FastqRecord],
    r2_by_id: dict[str, FastqRecord],
) -> tuple[int, int]:
    """Count mismatches against the reference using original quality scores.

    Returns ``(effective_sites, mismatches)``. Where both mates cover the same
    coordinate, equal bases count once at the higher quality, and different
    bases are excluded entirely.
    """
    effective_sites = 0
    mismatches = 0
    for entry in manifest:
        coverage: dict[int, object] = {}
        for coord, base, score in _oriented_bases(
            r1_by_id[entry.read_id], entry.r1_start, entry.r1_strand
        ):
            coverage[coord] = (base, score)
        for coord, base, score in _oriented_bases(
            r2_by_id[entry.read_id], entry.r2_start, entry.r2_strand
        ):
            existing = coverage.get(coord)
            if existing is None:
                coverage[coord] = (base, score)
            elif existing is _CONFLICT:
                continue
            else:
                prev_base, prev_score = existing  # type: ignore[misc]
                if prev_base == base:
                    if score > prev_score:
                        coverage[coord] = (base, score)
                else:
                    coverage[coord] = _CONFLICT
        for coord, value in coverage.items():
            if value is _CONFLICT:
                continue
            base, _score = value
            effective_sites += 1
            # Reference is stored 0-based; coordinates are 1-based.
            if base != reference[coord - 1]:
                mismatches += 1
    return effective_sites, mismatches


# ---------------------------------------------------------------------------
# Score computation and quality-line rewriting
# ---------------------------------------------------------------------------

def recalibrated_score(effective_sites: int, mismatches: int) -> int | None:
    """Return the recalibrated Phred score, or ``None`` to keep originals."""
    if effective_sites < config.MIN_EFFECTIVE_SITES:
        return None
    raw = -10.0 * math.log10((mismatches + 1) / (effective_sites + 2))
    # Half-up rounding (Python's built-in round() uses banker's rounding).
    score = math.floor(raw + 0.5)
    return max(config.MIN_PHRED, min(config.MAX_PHRED, score))


def build_result_zip(
    r1_records: list[FastqRecord],
    r2_by_id: dict[str, FastqRecord],
    score: int | None,
) -> bytes:
    """Rewrite quality lines and return the deliverable ZIP archive bytes.

    R1 keeps its original order; R2 is aligned to the same read ids. When
    ``score`` is ``None`` the original quality bytes are emitted untouched.
    """
    r1_out = io.BytesIO()
    r2_out = io.BytesIO()
    for r1 in r1_records:
        r2 = r2_by_id[r1.read_id]
        if score is None:
            r1_out.write(r1.original_bytes())
            r2_out.write(r2.original_bytes())
        else:
            q_char = (config.PHRED_OFFSET + score).to_bytes(1, "little")
            r1_out.write(r1.rewritten_bytes(q_char * len(r1.quality)))
            r2_out.write(r2.rewritten_bytes(q_char * len(r2.quality)))

    buffer = io.BytesIO()
    # The payloads are small (bounded by the 2000-pair / 200-base contract);
    # ZIP_DEFLATED keeps the archive compact while staying inside the memory
    # cap.
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("R1.fastq", r1_out.getvalue())
        zf.writestr("R2.fastq", r2_out.getvalue())
    return buffer.getvalue()


def run_recalibration(
    r1_path: Path,
    r2_path: Path,
    reference: str,
    manifest: list[ManifestEntry],
) -> tuple[bytes, dict]:
    """End-to-end recalibration used by the API endpoint.

    Validates everything first, then computes statistics and builds the ZIP.
    Any failure raises :class:`RecalError` before a single output byte is
    produced. Returns ``(zip_bytes, stats)``.
    """
    ref = validate_reference(reference)
    r1_records = parse_fastq_file(r1_path, 1)
    r2_records = parse_fastq_file(r2_path, 2)

    r1_by_id = _index_by_id(r1_records)
    r2_by_id = _index_by_id(r2_records)
    if len(r1_records) != len(r2_records):
        raise RecalError(
            f"R1 has {len(r1_records)} records but R2 has {len(r2_records)}; "
            "every read id must have both mates"
        )
    orphans_r2 = sorted(set(r2_by_id) - set(r1_by_id))
    if orphans_r2:
        raise RecalError(
            f"R2 contains read id without an R1 mate: '{orphans_r2[0]}'"
        )
    orphans_r1 = sorted(set(r1_by_id) - set(r2_by_id))
    if orphans_r1:
        raise RecalError(
            f"R1 contains read id without an R2 mate: '{orphans_r1[0]}'"
        )

    validate_manifest(
        manifest,
        set(r1_by_id),
        set(r2_by_id),
        len(ref),
        r1_by_id,
        r2_by_id,
    )

    effective_sites, mismatches = count_mismatches(manifest, ref, r1_by_id, r2_by_id)
    score = recalibrated_score(effective_sites, mismatches)
    zip_bytes = build_result_zip(r1_records, r2_by_id, score)

    stats = {
        "read_pairs": len(r1_records),
        "effective_sites": effective_sites,
        "mismatches": mismatches,
        "recalibrated_score": score,
        "score_applied": score is not None,
    }
    return zip_bytes, stats
