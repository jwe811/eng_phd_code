from __future__ import annotations

import json
import sqlite3
from datetime import UTC, datetime
from pathlib import Path
from typing import Protocol

from .config import settings
from .models import JobRecord, JobStatus


def utc_now() -> str:
    return datetime.now(UTC).isoformat()


class JobRepository(Protocol):
    def create(self, record: JobRecord) -> JobRecord: ...
    def get(self, job_id: str) -> JobRecord | None: ...
    def list(self) -> list[JobRecord]: ...
    def update(self, job_id: str, **fields: object) -> JobRecord: ...


class SQLiteJobRepository:
    """Local job persistence; designed to be swappable for Firestore later."""

    def __init__(self, path: Path = settings.sqlite_path) -> None:
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._init_db()

    def _connect(self) -> sqlite3.Connection:
        return sqlite3.connect(self.path)

    def _init_db(self) -> None:
        with self._connect() as db:
            db.execute(
                """
                CREATE TABLE IF NOT EXISTS jobs (
                    id TEXT PRIMARY KEY,
                    tool TEXT NOT NULL,
                    status TEXT NOT NULL,
                    command_json TEXT NOT NULL,
                    created_at TEXT NOT NULL,
                    started_at TEXT,
                    finished_at TEXT,
                    exit_code INTEGER,
                    log_path TEXT,
                    output_paths_json TEXT NOT NULL,
                    error TEXT
                )
                """
            )

    @staticmethod
    def _row_to_record(row: sqlite3.Row | tuple) -> JobRecord:
        return JobRecord(
            id=row[0],
            tool=row[1],
            status=row[2],
            command=json.loads(row[3]),
            created_at=row[4],
            started_at=row[5],
            finished_at=row[6],
            exit_code=row[7],
            log_path=row[8],
            output_paths=json.loads(row[9]),
            error=row[10],
        )

    def create(self, record: JobRecord) -> JobRecord:
        with self._connect() as db:
            db.execute(
                """
                INSERT INTO jobs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    record.id,
                    record.tool,
                    record.status,
                    json.dumps(record.command),
                    record.created_at,
                    record.started_at,
                    record.finished_at,
                    record.exit_code,
                    record.log_path,
                    json.dumps(record.output_paths),
                    record.error,
                ),
            )
        return record

    def get(self, job_id: str) -> JobRecord | None:
        with self._connect() as db:
            row = db.execute("SELECT * FROM jobs WHERE id = ?", (job_id,)).fetchone()
        return self._row_to_record(row) if row else None

    def list(self) -> list[JobRecord]:
        with self._connect() as db:
            rows = db.execute("SELECT * FROM jobs ORDER BY created_at DESC").fetchall()
        return [self._row_to_record(row) for row in rows]

    def update(self, job_id: str, **fields: object) -> JobRecord:
        if not fields:
            found = self.get(job_id)
            if found is None:
                raise KeyError(job_id)
            return found
        allowed = {
            "status",
            "started_at",
            "finished_at",
            "exit_code",
            "log_path",
            "output_paths",
            "error",
        }
        unknown = set(fields) - allowed
        if unknown:
            raise ValueError(f"Unknown job fields: {sorted(unknown)}")
        columns: list[str] = []
        values: list[object] = []
        for key, value in fields.items():
            if key == "output_paths":
                key = "output_paths_json"
                value = json.dumps(value)
            columns.append(f"{key} = ?")
            values.append(value)
        values.append(job_id)
        with self._connect() as db:
            db.execute(f"UPDATE jobs SET {', '.join(columns)} WHERE id = ?", values)
        found = self.get(job_id)
        if found is None:
            raise KeyError(job_id)
        return found

