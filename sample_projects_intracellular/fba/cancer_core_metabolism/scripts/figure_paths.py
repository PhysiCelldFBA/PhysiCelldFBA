"""Path setup for cancer-tissue figure scripts (pctk + local plot modules)."""
from __future__ import annotations

import os
import sys
from pathlib import Path

_SCRIPTS = Path(__file__).resolve().parent
_PROJECT = _SCRIPTS.parent  # cancer_core_metabolism
_MN5 = _PROJECT.parents[2]  # mn5sync_pdfba
_REPO = _MN5.parent  # pDFBA (local) — may not exist on cluster layout


def _as_cwd_relative(path: Path) -> Path:
    """Return ``path`` expressed relative to the current working directory."""
    return Path(os.path.relpath(path.resolve(), Path.cwd().resolve()))


def find_scripts_dir(start: Path | None = None) -> Path:
    """Locate this scripts/ folder from cwd (or ``start``) via relative candidates."""
    start = Path(start or Path.cwd())
    rel_candidates = [
        Path("."),
        Path("scripts"),
        Path("mn5sync_pdfba/sample_projects_intracellular/fba/cancer_core_metabolism/scripts"),
        Path("../mn5sync_pdfba/sample_projects_intracellular/fba/cancer_core_metabolism/scripts"),
        Path("../../mn5sync_pdfba/sample_projects_intracellular/fba/cancer_core_metabolism/scripts"),
        Path("sample_projects_intracellular/fba/cancer_core_metabolism/scripts"),
        Path("../sample_projects_intracellular/fba/cancer_core_metabolism/scripts"),
    ]
    for rel in rel_candidates:
        cand = (start / rel).resolve()
        if (cand / "figure_paths.py").is_file():
            return _as_cwd_relative(cand)

    # Walk parents looking for the sample-project scripts marker
    for parent in [start.resolve(), *start.resolve().parents]:
        cand = (
            parent
            / "mn5sync_pdfba"
            / "sample_projects_intracellular"
            / "fba"
            / "cancer_core_metabolism"
            / "scripts"
        )
        if (cand / "figure_paths.py").is_file():
            return _as_cwd_relative(cand)
        if parent.name == "cancer_core_metabolism" and (parent / "scripts" / "figure_paths.py").is_file():
            return _as_cwd_relative(parent / "scripts")

    # Last resort: directory containing this file (may be absolute)
    return _as_cwd_relative(_SCRIPTS)


def ensure_paths(scripts_dir: Path | None = None) -> Path:
    """Insert scripts/ and python_libs on sys.path. Return scripts dir."""
    scripts_path = Path(scripts_dir) if scripts_dir is not None else find_scripts_dir()
    # sys.path needs an existing filesystem path; resolve for imports
    scripts_abs = str(scripts_path.resolve())
    if scripts_abs not in sys.path:
        sys.path.insert(0, scripts_abs)

    # python_libs: relative candidates first, then layout next to mn5sync_pdfba
    lib_candidates = [
        Path("python_libs"),
        Path("../python_libs"),
        Path("../../python_libs"),
        Path("mn5sync_pdfba/python_libs"),
        Path("../mn5sync_pdfba/python_libs"),
        Path("../../mn5sync_pdfba/python_libs"),
        _MN5 / "python_libs",
        _REPO / "python_libs",
    ]
    for lib in lib_candidates:
        lib_abs = lib.resolve() if not lib.is_absolute() else lib
        # relative candidates are from cwd
        if not lib.is_absolute():
            lib_abs = (Path.cwd() / lib).resolve()
        if lib_abs.is_dir():
            p = str(lib_abs)
            if p not in sys.path:
                sys.path.insert(0, p)
            break
    return scripts_path


def default_output_dir() -> Path:
    """Best-effort PhysiCell output folder, returned relative to cwd when possible."""
    candidates = [
        Path("output_lit"),
        Path("output"),
        Path("test_parameters_results/cancer_tissue/output_lit"),
        Path("../output_lit"),
        Path("../../test_parameters_results/cancer_tissue/output_lit"),
        Path("../../../test_parameters_results/cancer_tissue/output_lit"),
        Path("../../../../test_parameters_results/cancer_tissue/output_lit"),
        _PROJECT / "output",
        _REPO / "test_parameters_results" / "cancer_tissue" / "output_lit",
    ]
    for c in candidates:
        probe = c if c.is_absolute() else (Path.cwd() / c)
        if probe.is_dir():
            return _as_cwd_relative(probe)
    return candidates[0]
