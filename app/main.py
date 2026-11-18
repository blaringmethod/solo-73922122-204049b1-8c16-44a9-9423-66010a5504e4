"""FastAPI application: bounded uploads, status, cancellation, downloads."""

from __future__ import annotations

import sqlite3
import tempfile
import uuid
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Literal

from fastapi import FastAPI, File, HTTPException, UploadFile, status
from fastapi.responses import FileResponse, JSONResponse, Response
from pydantic import BaseModel, ConfigDict, Field

from . import config, db, worker
from .recalibrate import ManifestEntry, RecalError, run_recalibration
from .worker import UploadTooLarge, start_worker, task_dir


@asynccontextmanager
async def lifespan(_: FastAPI):
    db.init_db()
    worker.recover_on_startup()
    yield


app = FastAPI(
    title="Paired FASTQ Service",
    version="1.1.0",
    description=(
        "Upload shuffled paired-end FASTQ files (R1/R2, uncompressed, "
        "max 256 MiB each). The service re-pairs them by read id, delivers "
        "both mates in R1 order and preserves the original bytes. A "
        "successful task can then be recalibrated: quality lines are "
        "rewritten from the control/reference mismatch rate and returned "
        "as a ZIP archive."
    ),
    lifespan=lifespan,
)


class ManifestItem(BaseModel):
    model_config = ConfigDict(extra="forbid")

    read_id: str = Field(..., min_length=1, max_length=1024,
                         description="Read id exactly as in the FASTQ headers")
    r1_start: int = Field(..., ge=1, description="1-based R1 start coordinate")
    r1_strand: Literal["+", "-"] = Field(..., description="R1 strand orientation")
    r2_start: int = Field(..., ge=1, description="1-based R2 start coordinate")
    r2_strand: Literal["+", "-"] = Field(..., description="R2 strand orientation")


class RecalibrateRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    reference: str = Field(
        ..., min_length=1, max_length=config.MAX_REF_BASES,
        description="Reference sequence, ACGT only, at most 10000 bases",
    )
    manifest: list[ManifestItem] = Field(
        ..., min_length=1, max_length=config.MAX_READ_PAIRS,
        description="One entry per read pair of the referenced task",
    )


class TaskOut(BaseModel):
    id: str
    status: str
    stage: str
    progress: float
    message: str
    r1_name: str
    r2_name: str
    r1_bytes: int
    r2_bytes: int
    created_at: str
    updated_at: str
    download_r1: str | None = None
    download_r2: str | None = None


def _serialize(row: sqlite3.Row) -> TaskOut:
    download_r1 = download_r2 = None
    if row["status"] == "succeeded" and row["result_r1"]:
        download_r1 = f"/tasks/{row['id']}/download/r1"
        download_r2 = f"/tasks/{row['id']}/download/r2"
    return TaskOut(
        id=row["id"],
        status=row["status"],
        stage=row["stage"],
        progress=round(row["progress"], 4),
        message=row["message"],
        r1_name=row["r1_name"],
        r2_name=row["r2_name"],
        r1_bytes=row["r1_bytes"],
        r2_bytes=row["r2_bytes"],
        created_at=row["created_at"],
        updated_at=row["updated_at"],
        download_r1=download_r1,
        download_r2=download_r2,
    )


async def _spool_bounded(upload: UploadFile, side: str, dest: Path) -> int:
    """Stream an upload to ``dest`` without holding it in memory.

    SpooledTemporaryFile rolls over to a real file after
    ``UPLOAD_SPOOL_THRESHOLD`` and the copy proceeds in fixed chunks; the
    per-file 256 MiB limit is enforced while streaming.
    """
    spool = tempfile.SpooledTemporaryFile(max_size=config.UPLOAD_SPOOL_THRESHOLD)
    total = 0
    try:
        while True:
            chunk = await upload.read(1024 * 1024)
            if not chunk:
                break
            total += len(chunk)
            if total > config.MAX_FILE_BYTES:
                raise UploadTooLarge(side, config.MAX_FILE_BYTES)
            spool.write(chunk)
        spool.rollover()
        spool.seek(0)
        with dest.open("wb") as out:
            while True:
                chunk = spool.read(1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)
    finally:
        spool.close()
    return total


@app.post("/tasks", response_model=TaskOut, status_code=status.HTTP_202_ACCEPTED)
async def create_task(
    r1: UploadFile = File(..., alias="r1", description="R1 FASTQ, uncompressed, <= 256 MiB"),
    r2: UploadFile = File(..., alias="r2", description="R2 FASTQ, uncompressed, <= 256 MiB"),
) -> TaskOut:
    task_id = uuid.uuid4().hex
    base = task_dir(task_id)
    base.mkdir(parents=True, exist_ok=False)
    r1_dest = base / "input_R1.fastq"
    r2_dest = base / "input_R2.fastq"
    try:
        size1 = await _spool_bounded(r1, "R1", r1_dest)
        size2 = await _spool_bounded(r2, "R2", r2_dest)
    except UploadTooLarge as exc:
        worker.cleanup_task_dir(task_id)
        raise HTTPException(
            status_code=getattr(
                status, "HTTP_413_CONTENT_TOO_LARGE",
                status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
            ),
            detail=str(exc),
        )
    except Exception:
        worker.cleanup_task_dir(task_id)
        raise
    finally:
        await r1.close()
        await r2.close()

    conn = db.connect()
    try:
        db.create_task(
            conn,
            task_id,
            r1_name=Path(r1.filename or "r1.fastq").name,
            r2_name=Path(r2.filename or "r2.fastq").name,
            r1_bytes=size1,
            r2_bytes=size2,
        )
        row = db.get_task(conn, task_id)
    finally:
        conn.close()

    worker.register_cancel_event(task_id)
    start_worker(task_id)
    return _serialize(row)


@app.get("/tasks/{task_id}", response_model=TaskOut)
def get_task(task_id: str) -> TaskOut:
    conn = db.connect()
    try:
        row = db.get_task(conn, task_id)
        if row is None:
            raise HTTPException(status_code=404, detail="task not found")
        return _serialize(row)
    finally:
        conn.close()


@app.post("/tasks/{task_id}/cancel", response_model=TaskOut)
def cancel_task(task_id: str) -> TaskOut:
    conn = db.connect()
    try:
        previous = db.try_set_cancelling(conn, task_id)
        if previous is None:
            raise HTTPException(status_code=404, detail="task not found")
        row = db.get_task(conn, task_id)
        if row["status"] in ("succeeded", "failed", "cancelled"):
            # The task reached a terminal state first; tell the caller plainly
            # that it cannot be cancelled.
            raise HTTPException(
                status_code=409,
                detail=f"task is already {row['status']} and cannot be cancelled",
            )
        # Flip the cooperative flag and return without waiting for the worker.
        worker.request_cancel(task_id)
        return _serialize(row)
    finally:
        conn.close()


def _result_path(task_id: str, side: int) -> Path:
    conn = db.connect()
    try:
        row = db.get_task(conn, task_id)
        if row is None:
            raise HTTPException(status_code=404, detail="task not found")
        if row["status"] != "succeeded":
            raise HTTPException(
                status_code=409,
                detail=f"task is {row['status']}; results are not published",
            )
        rel = row["result_r1"] if side == 1 else row["result_r2"]
        if not rel:
            raise HTTPException(status_code=404, detail="result file missing")
    finally:
        conn.close()

    base = task_dir(task_id).resolve()
    path = (base / rel).resolve()
    # Guard against path traversal even though the rel value is server-generated.
    if base not in path.parents or not path.is_file():
        raise HTTPException(status_code=404, detail="result file missing")
    return path


@app.get("/tasks/{task_id}/download/r1")
def download_r1(task_id: str) -> FileResponse:
    path = _result_path(task_id, 1)
    return FileResponse(
        path, media_type="text/plain", filename="R1.fastq",
        headers={"Content-Length": str(path.stat().st_size)},
    )


@app.get("/tasks/{task_id}/download/r2")
def download_r2(task_id: str) -> FileResponse:
    path = _result_path(task_id, 2)
    return FileResponse(
        path, media_type="text/plain", filename="R2.fastq",
        headers={"Content-Length": str(path.stat().st_size)},
    )


def _recalibration_inputs(task_id: str) -> tuple[Path, Path]:
    """Resolve the two published FASTQ files of a successful task."""
    conn = db.connect()
    try:
        row = db.get_task(conn, task_id)
        if row is None:
            raise HTTPException(status_code=404, detail="task not found")
        if row["status"] != "succeeded":
            raise HTTPException(
                status_code=409,
                detail=(
                    f"task is {row['status']}; recalibration requires a "
                    "successful task id"
                ),
            )
        rel1, rel2 = row["result_r1"], row["result_r2"]
        if not rel1 or not rel2:
            raise HTTPException(
                status_code=409,
                detail="result file missing; task has no published FASTQ",
            )
    finally:
        conn.close()

    base = task_dir(task_id).resolve()
    path1 = (base / rel1).resolve()
    path2 = (base / rel2).resolve()
    if base not in path1.parents or base not in path2.parents:
        raise HTTPException(status_code=404, detail="result file missing")
    if not path1.is_file() or not path2.is_file():
        raise HTTPException(
            status_code=409,
            detail="result file missing on disk; recalibration cannot run",
        )
    return path1, path2


@app.post("/tasks/{task_id}/recalibrate")
def recalibrate(task_id: str, payload: RecalibrateRequest) -> Response:
    """Recalibrate quality lines of a successful task and return a ZIP.

    The whole request is all-or-nothing: any missing/duplicate/unknown read
    id, out-of-bounds coordinate, missing result file or corrupted FASTQ
    fails with an error response and no ZIP (or partial output) is produced;
    the original task files are never modified.
    """
    path1, path2 = _recalibration_inputs(task_id)
    manifest = [
        ManifestEntry(
            read_id=item.read_id,
            r1_start=item.r1_start,
            r1_strand=item.r1_strand,
            r2_start=item.r2_start,
            r2_strand=item.r2_strand,
        )
        for item in payload.manifest
    ]
    try:
        zip_bytes, stats = run_recalibration(
            path1, path2, payload.reference, manifest
        )
    except RecalError as exc:
        # Domain validation failures are request errors, not server crashes.
        raise HTTPException(status_code=422, detail=str(exc))

    headers = {
        "Content-Disposition": (
            f'attachment; filename="recalibrated-{task_id}.zip"'
        ),
        "X-Effective-Sites": str(stats["effective_sites"]),
        "X-Mismatches": str(stats["mismatches"]),
        "X-Score-Applied": "1" if stats["score_applied"] else "0",
        "X-Recalibrated-Score": (
            str(stats["recalibrated_score"])
            if stats["recalibrated_score"] is not None
            else ""
        ),
    }
    return Response(
        content=zip_bytes,
        media_type="application/zip",
        headers=headers,
    )


@app.get("/health")
def health() -> JSONResponse:
    return JSONResponse({"status": "ok"})
