"""Runtime configuration for the paired FASTQ service.

All limits can be overridden through environment variables, but the defaults
match the delivery requirement: each uploaded file <= 256 MiB and the whole
process is expected to live inside a 256 MiB memory cap (see
docker-compose.yml).
"""

from __future__ import annotations

import os
from pathlib import Path

MAX_FILE_BYTES = int(os.environ.get("FASTQ_MAX_FILE_BYTES", 256 * 1024 * 1024))
"""Maximum size of a single uploaded FASTQ file (256 MiB)."""

DATA_DIR = Path(os.environ.get("FASTQ_DATA_DIR", "/data"))
"""Root directory for the SQLite database and all task files."""

DB_PATH = DATA_DIR / "tasks.db"
TASKS_DIR = DATA_DIR / "tasks"

# Uploads are spooled to a temporary file once they exceed this size, so a big
# multipart upload never accumulates in RAM.
UPLOAD_SPOOL_THRESHOLD = int(
    os.environ.get("FASTQ_UPLOAD_SPOOL_THRESHOLD", 1 * 1024 * 1024)
)

# Records are inserted in batches; only one batch is held in memory at a time.
INSERT_BATCH_SIZE = int(os.environ.get("FASTQ_INSERT_BATCH_SIZE", 512))

# SQLite page cache kept deliberately small. The heavy data stays on disk.
SQLITE_CACHE_PAGES = int(os.environ.get("FASTQ_SQLITE_CACHE_PAGES", 2048))
SQLITE_BUSY_TIMEOUT_MS = 30_000

# How often the worker refreshes the progress percentage in the database.
PROGRESS_UPDATE_EVERY_BYTES = 8 * 1024 * 1024

STAGE_UPLOAD = "upload"
STAGE_PARSE = "parse"
STAGE_VALIDATE = "validate"
STAGE_PUBLISH = "publish"
STAGE_DONE = "done"

# ---------------------------------------------------------------------------
# Quality-score recalibration limits (all overridable via environment).
# ---------------------------------------------------------------------------

MAX_READ_PAIRS = int(os.environ.get("FASTQ_MAX_READ_PAIRS", 2000))
"""Maximum number of read pairs accepted by the recalibration endpoint."""

MIN_READ_LEN = int(os.environ.get("FASTQ_MIN_READ_LEN", 40))
MAX_READ_LEN = int(os.environ.get("FASTQ_MAX_READ_LEN", 200))
"""Per-end sequence length window (inclusive, ACGT bases only)."""

MAX_REF_BASES = int(os.environ.get("FASTQ_MAX_REF_BASES", 10_000))
"""Maximum reference length in bases."""

PHRED_OFFSET = 33
MIN_PHRED = 0
MAX_PHRED = 40
"""Accepted/recalibrated quality range for Phred+33 ('!'..'I')."""

MIN_EFFECTIVE_SITES = int(os.environ.get("FASTQ_MIN_EFFECTIVE_SITES", 20))
"""Effective sites required before the recalibrated score replaces originals."""
