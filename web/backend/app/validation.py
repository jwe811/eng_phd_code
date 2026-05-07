from __future__ import annotations

from pathlib import Path

from .config import settings
from .models import JobCreate


MODE_LABELS = {
    0: "Standard SAP",
    1: "Hamiltonian SAP",
    2: "2SAP",
    3: "2SAP-Hamiltonian",
}


class ValidationError(ValueError):
    pass


def validate_lattice(lattice_l: int, lattice_m: int) -> None:
    if lattice_l < 0 or lattice_m < 0:
        raise ValidationError("L and M must be non-negative.")
    if (lattice_l + 1) * (lattice_m + 1) >= 10:
        raise ValidationError("(L + 1) * (M + 1) must be less than 10 for v1 local runs.")


def validate_job(request: JobCreate) -> None:
    validate_lattice(request.L, request.M)
    if request.mode not in MODE_LABELS:
        raise ValidationError("mode must be one of 0, 1, 2, or 3.")

    if request.tool in {"mc", "creator"}:
        if request.span is None:
            raise ValidationError("span is required for Monte Carlo and CreatorAll runs.")
        if not 2 <= request.span <= 8:
            raise ValidationError("span must be between 2 and 8.")

    if request.tool == "mc":
        if request.samples is None:
            raise ValidationError("samples is required for Monte Carlo runs.")
        if not 1 <= request.samples <= 100:
            raise ValidationError("samples must be between 1 and 100.")
        if request.run is None:
            raise ValidationError("run number is required for Monte Carlo runs.")
        if request.seed is None:
            raise ValidationError("seed is required for Monte Carlo runs.")


def build_command(request: JobCreate) -> list[str]:
    validate_job(request)
    if request.tool == "tm":
        command = [str(settings.bin_root / "tm_master"), "-L", str(request.L), "-M", str(request.M), "-m", str(request.mode)]
        if request.options.fugacity is not None:
            command += ["-x", str(request.options.fugacity)]
        if request.options.export_eigenvectors:
            command.append("-E")
        if request.options.damping:
            command.append("-d")
        if request.options.section_capacity is not None:
            command += ["-S", str(request.options.section_capacity)]
        if request.options.state_capacity is not None:
            command += ["-K", str(request.options.state_capacity)]
        if request.options.convergence is not None:
            command += ["-c", str(request.options.convergence)]
        return command

    if request.tool == "mc":
        return [
            str(settings.bin_root / "mc_master"),
            "-L",
            str(request.L),
            "-M",
            str(request.M),
            "-m",
            str(request.mode),
            "-s",
            str(request.span),
            "-n",
            str(request.samples),
            "-r",
            str(request.run),
            "-S",
            str(request.seed),
        ]

    return [
        str(settings.bin_root / "creator_all"),
        "-L",
        str(request.L),
        "-M",
        str(request.M),
        "-m",
        str(request.mode),
        "-s",
        str(request.span),
    ]


def repo_relative(path: Path) -> str:
    return str(path.resolve().relative_to(settings.root.resolve()))

