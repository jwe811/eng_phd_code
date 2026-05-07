from __future__ import annotations

import subprocess
import threading
import uuid
from pathlib import Path
from typing import Protocol

from .config import settings
from .models import JobCreate, JobRecord
from .repositories import JobRepository, utc_now
from .validation import build_command, repo_relative


class JobRunner(Protocol):
    def start(self, request: JobCreate) -> JobRecord: ...
    def cancel(self, job_id: str) -> JobRecord: ...


class LocalSubprocessRunner:
    """Local process runner; future Cloud Run/Tasks execution can replace this."""

    def __init__(self, jobs: JobRepository) -> None:
        self.jobs = jobs
        self.processes: dict[str, subprocess.Popen[str]] = {}
        self.lock = threading.Lock()
        settings.jobs_root.mkdir(parents=True, exist_ok=True)

    def start(self, request: JobCreate) -> JobRecord:
        command = build_command(request)
        job_id = uuid.uuid4().hex[:12]
        job_dir = settings.jobs_root / job_id
        job_dir.mkdir(parents=True, exist_ok=True)
        log_path = job_dir / "run.log"
        record = JobRecord(
            id=job_id,
            tool=request.tool,
            status="queued",
            command=command,
            created_at=utc_now(),
            log_path=repo_relative(log_path),
        )
        self.jobs.create(record)
        thread = threading.Thread(target=self._run, args=(job_id, command, log_path), daemon=True)
        thread.start()
        return self.jobs.get(job_id) or record

    def cancel(self, job_id: str) -> JobRecord:
        with self.lock:
            proc = self.processes.get(job_id)
        current = self.jobs.get(job_id)
        if current is None:
            raise KeyError(job_id)
        if proc and proc.poll() is None:
            proc.terminate()
            return self.jobs.update(job_id, status="cancelled", finished_at=utc_now())
        return current

    def _run(self, job_id: str, command: list[str], log_path: Path) -> None:
        before = self._snapshot_outputs()
        self.jobs.update(job_id, status="running", started_at=utc_now())
        with log_path.open("w", encoding="utf-8") as log:
            log.write("$ " + " ".join(command) + "\n\n")
            log.flush()
            proc = subprocess.Popen(
                command,
                cwd=settings.root,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
            )
            with self.lock:
                self.processes[job_id] = proc
            exit_code = proc.wait()
        with self.lock:
            self.processes.pop(job_id, None)
        after = self._snapshot_outputs()
        changed = {path for path, fingerprint in after.items() if before.get(path) != fingerprint}
        outputs = sorted(repo_relative(path) for path in changed)
        status = "succeeded" if exit_code == 0 else "failed"
        current = self.jobs.get(job_id)
        if current and current.status == "cancelled":
            return
        self.jobs.update(
            job_id,
            status=status,
            finished_at=utc_now(),
            exit_code=exit_code,
            output_paths=outputs,
            error=None if exit_code == 0 else f"Command exited with {exit_code}",
        )

    def _snapshot_outputs(self) -> dict[Path, tuple[int, int]]:
        if not settings.data_root.exists():
            return {}
        return {
            path.resolve(): (path.stat().st_size, path.stat().st_mtime_ns)
            for path in settings.data_root.rglob("*")
            if path.is_file()
        }
