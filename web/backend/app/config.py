from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = ROOT / "data"
BIN_ROOT = ROOT / "bin"
RUNTIME_ROOT = ROOT / "web" / "backend" / "runtime"


@dataclass(frozen=True)
class Settings:
    root: Path = ROOT
    data_root: Path = DATA_ROOT
    bin_root: Path = BIN_ROOT
    runtime_root: Path = RUNTIME_ROOT
    sqlite_path: Path = RUNTIME_ROOT / "jobs.sqlite3"
    jobs_root: Path = RUNTIME_ROOT / "jobs"


settings = Settings()

