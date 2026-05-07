from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field


ToolName = Literal["tm", "mc", "creator"]
JobStatus = Literal["queued", "running", "succeeded", "failed", "cancelled"]


class JobOptions(BaseModel):
    export_eigenvectors: bool = False
    damping: bool = False
    fugacity: float | None = None
    section_capacity: int | None = Field(default=None, ge=1)
    state_capacity: int | None = Field(default=None, ge=1)
    convergence: float | None = Field(default=None, gt=0)


class JobCreate(BaseModel):
    tool: ToolName
    L: int = Field(ge=0)
    M: int = Field(ge=0)
    mode: int = Field(ge=0, le=3)
    span: int | None = None
    samples: int | None = None
    run: int | None = None
    seed: int | None = None
    options: JobOptions = Field(default_factory=JobOptions)


class JobRecord(BaseModel):
    id: str
    tool: ToolName
    status: JobStatus
    command: list[str]
    created_at: str
    started_at: str | None = None
    finished_at: str | None = None
    exit_code: int | None = None
    log_path: str | None = None
    output_paths: list[str] = Field(default_factory=list)
    error: str | None = None


class JobDetail(JobRecord):
    log_tail: str = ""


class ToolStatus(BaseModel):
    name: ToolName
    label: str
    path: str
    available: bool


class ResultFile(BaseModel):
    path: str
    group: str
    size_bytes: int
    modified_at: str
    meta_path: str | None = None
    kind: str


class DataBrowserItem(BaseModel):
    name: str
    path: str
    is_dir: bool
    kind: str
    size_bytes: int | None = None


class DataBrowserListing(BaseModel):
    path: str
    parent_path: str | None = None
    items: list[DataBrowserItem]
