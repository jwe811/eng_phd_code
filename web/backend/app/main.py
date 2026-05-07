from __future__ import annotations

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware

from .config import settings
from .models import DataBrowserListing, JobCreate, JobDetail, ToolStatus
from .analysis_api import run_analysis
from .repositories import SQLiteJobRepository
from .runner import LocalSubprocessRunner
from .storage import LocalResultStorage
from .uofs_api import object_points, summarize_uofs
from .validation import ValidationError


app = FastAPI(title="SAP/2SAP Workbench API")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

jobs = SQLiteJobRepository()
runner = LocalSubprocessRunner(jobs)
storage = LocalResultStorage()


@app.get("/api/health")
def health() -> dict:
    return {"ok": True}


@app.get("/api/tools", response_model=list[ToolStatus])
def tools() -> list[ToolStatus]:
    specs = [
        ("tm", "Transfer Matrix", settings.bin_root / "tm_master"),
        ("mc", "Monte Carlo", settings.bin_root / "mc_master"),
        ("creator", "CreatorAll", settings.bin_root / "creator_all"),
    ]
    return [
        ToolStatus(name=name, label=label, path=str(path), available=path.exists())
        for name, label, path in specs
    ]


@app.post("/api/jobs")
def create_job(request: JobCreate):
    try:
        return runner.start(request)
    except ValidationError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


@app.get("/api/jobs")
def list_jobs():
    return jobs.list()


@app.get("/api/jobs/{job_id}", response_model=JobDetail)
def get_job(job_id: str):
    record = jobs.get(job_id)
    if record is None:
        raise HTTPException(status_code=404, detail="job not found")
    log_tail = ""
    if record.log_path:
        path = (settings.root / record.log_path).resolve()
        if path.exists() and settings.runtime_root.resolve() in path.parents:
            log_tail = path.read_text(encoding="utf-8", errors="replace")[-12000:]
    return JobDetail(**record.model_dump(), log_tail=log_tail)


@app.post("/api/jobs/{job_id}/cancel")
def cancel_job(job_id: str):
    try:
        return runner.cancel(job_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="job not found") from exc


@app.get("/api/results")
def results():
    return storage.list_results()


@app.get("/api/data-browser", response_model=DataBrowserListing)
def data_browser(path: str = "data"):
    try:
        return storage.browse(path)
    except (ValueError, OSError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.get("/api/uofs/summary")
def uofs_summary(path: str = Query(...), polys_per_object: int | None = None):
    try:
        resolved = storage.resolve_data_path(path)
        return summarize_uofs(resolved, polys_per_object)
    except (ValueError, OSError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.get("/api/uofs/object")
def uofs_object(path: str = Query(...), index: int = 0, polys_per_object: int | None = None):
    try:
        resolved = storage.resolve_data_path(path)
        return object_points(resolved, index, polys_per_object)
    except IndexError as exc:
        raise HTTPException(status_code=404, detail="object index not found") from exc
    except (ValueError, OSError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.get("/api/analysis")
def analysis(path: str = Query(...), action: str = Query(...)):
    try:
        resolved = storage.resolve_data_path(path)
        return run_analysis(resolved, action)
    except (ValueError, OSError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.get("/api/file-text")
def file_text(path: str = Query(...)):
    try:
        resolved = storage.resolve_text_path(path)
        return {"path": path, "text": resolved.read_text(encoding="utf-8", errors="replace")[-20000:]}
    except (ValueError, OSError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
