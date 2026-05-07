from __future__ import annotations

from app import main
from app.validation import ValidationError, build_command
from app.models import JobCreate


def test_health() -> None:
    assert main.health() == {"ok": True}


def test_tools() -> None:
    names = {tool.name for tool in main.tools()}
    assert names == {"tm", "mc", "creator"}


def test_lattice_guardrail_rejects_product_at_ten() -> None:
    request = JobCreate(tool="creator", L=4, M=1, mode=0, span=2)
    try:
        build_command(request)
    except ValidationError as exc:
        assert "(L + 1) * (M + 1)" in str(exc)
    else:
        raise AssertionError("expected guardrail failure")


def test_lattice_guardrail_allows_small_command() -> None:
    command = build_command(JobCreate(tool="creator", L=1, M=1, mode=0, span=2))
    assert command[1:] == ["-L", "1", "-M", "1", "-m", "0", "-s", "2"]


def test_uofs_summary_fixture() -> None:
    body = main.uofs_summary("data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt")
    assert body["objects"] >= 1
    assert body["polys_per_object"] == 1


def test_results_scan() -> None:
    assert isinstance(main.results(), list)


def test_data_browser_lists_data_root() -> None:
    listing = main.data_browser("data")
    assert listing.path == "data"
    names = {item.name for item in listing.items}
    assert {"CreatorAll", "MonteCarlo"}.issubset(names)
    assert "TransferMatrix" not in names


def test_data_browser_hides_monte_carlo_evectors() -> None:
    listing = main.data_browser("data/MonteCarlo")
    names = {item.name for item in listing.items}
    assert "MC_Evectors" not in names
    assert names.issubset({"SAPs", "HamSAPs", "2SAPs", "Ham2SAPs"})


def test_analysis_summary_and_validate_fixture() -> None:
    path = "data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt"
    summary = main.analysis(path, "summary")
    assert summary["metrics"]["objects"] >= 1
    validate = main.analysis(path, "validate")
    assert validate["metrics"]["validated_objects"] >= 1


def test_file_text_reads_metadata_fixture() -> None:
    body = main.file_text("data/CreatorAll/All_SAPs/AllSAPsL1M1span2.summary.meta")
    assert "text" in body


def test_job_rejection_from_api() -> None:
    try:
        main.create_job(JobCreate(tool="mc", L=4, M=1, mode=0, span=2, samples=1, run=1, seed=1))
    except Exception as exc:
        assert getattr(exc, "status_code", None) == 422
    else:
        raise AssertionError("expected API rejection")
