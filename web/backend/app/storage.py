from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path
from typing import Protocol

from .config import settings
from .models import DataBrowserItem, DataBrowserListing, ResultFile
from .validation import repo_relative


class ResultStorage(Protocol):
    def list_results(self) -> list[ResultFile]: ...
    def browse(self, relative_path: str) -> DataBrowserListing: ...
    def resolve_data_path(self, relative_path: str) -> Path: ...
    def resolve_text_path(self, relative_path: str) -> Path: ...


class LocalResultStorage:
    """Local data storage; the Cloud Storage equivalent can keep this API."""

    def __init__(self, data_root: Path = settings.data_root) -> None:
        self.data_root = data_root

    def resolve_data_path(self, relative_path: str) -> Path:
        path = (settings.root / relative_path).resolve()
        data_root = self.data_root.resolve()
        if not path.is_file() or data_root not in path.parents:
            raise ValueError("path must be an existing file under data/")
        return path

    def resolve_text_path(self, relative_path: str) -> Path:
        path = (settings.root / relative_path).resolve()
        data_root = self.data_root.resolve()
        runtime_root = settings.runtime_root.resolve()
        allowed = data_root in path.parents or runtime_root in path.parents
        if not path.is_file() or not allowed or path.suffix not in {".txt", ".meta", ".log"}:
            raise ValueError("path must be a text, metadata, or log file under data/ or web/backend/runtime/")
        return path

    def resolve_data_dir(self, relative_path: str) -> Path:
        path = (settings.root / relative_path).resolve()
        data_root = self.data_root.resolve()
        if path == data_root:
            return path
        if not path.is_dir() or data_root not in path.parents:
            raise ValueError("path must be an existing directory under data/")
        return path

    def browse(self, relative_path: str = "data") -> DataBrowserListing:
        directory = self.resolve_data_dir(relative_path)
        data_root = self.data_root.resolve()
        parent = directory.parent if directory != data_root else None
        items: list[DataBrowserItem] = []
        monte_carlo_root = (data_root / "MonteCarlo").resolve()
        monte_carlo_polygon_dirs = {"SAPs", "HamSAPs", "2SAPs", "Ham2SAPs"}
        for child in sorted(directory.iterdir(), key=lambda item: (not item.is_dir(), item.name.lower())):
            if child.name.startswith("."):
                continue
            if child.is_dir():
                items.append(DataBrowserItem(name=child.name, path=repo_relative(child), is_dir=True, kind="directory"))
            elif child.is_file() and child.suffix != ".meta":
                items.append(
                    DataBrowserItem(
                        name=child.name,
                        path=repo_relative(child),
                        is_dir=False,
                        kind=self._kind(child),
                        size_bytes=child.stat().st_size,
                    )
                )
        return DataBrowserListing(
            path=repo_relative(directory),
            parent_path=repo_relative(parent) if parent else None,
            items=items,
        )

    def list_results(self) -> list[ResultFile]:
        if not self.data_root.exists():
            return []
        results: list[ResultFile] = []
        for path in sorted(self.data_root.rglob("*")):
            if not path.is_file() or path.suffix == ".meta":
                continue
            stat = path.stat()
            rel = repo_relative(path)
            meta = path.with_name(path.name + ".meta")
            results.append(
                ResultFile(
                    path=rel,
                    group=self._group(path),
                    size_bytes=stat.st_size,
                    modified_at=datetime.fromtimestamp(stat.st_mtime, UTC).isoformat(),
                    meta_path=repo_relative(meta) if meta.exists() else None,
                    kind=self._kind(path),
                )
            )
        return results

    def _group(self, path: Path) -> str:
        rel = path.relative_to(self.data_root)
        return rel.parts[0] if rel.parts else "data"

    def _kind(self, path: Path) -> str:
        if path.suffix == ".txt":
            try:
                with path.open("r", encoding="utf-8") as fp:
                    first = fp.readline().strip()
                if first == "UofS":
                    return "uofs"
            except UnicodeDecodeError:
                pass
            return "text"
        if path.suffix == ".bin":
            return "binary"
        return path.suffix.lstrip(".") or "file"
