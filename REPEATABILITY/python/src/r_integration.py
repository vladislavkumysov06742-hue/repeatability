"""Helpers to list and optionally run R scripts bundled in `python/r_scripts/`.

This module does not execute R by default. It provides a safe wrapper that:
- lists available R scripts copied into `python/r_scripts/`.
- can run a script via Rscript if the user explicitly calls `run_r_script` and Rscript is available.

Use this for integration tests or to document mappings between R and Python implementations.
"""
from pathlib import Path
import shutil
import subprocess
from typing import List


R_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "r_scripts"


def list_r_scripts() -> List[Path]:
    """Return list of R script paths in the bundled r_scripts directory."""
    if not R_SCRIPTS_DIR.exists():
        return []
    return sorted([p for p in R_SCRIPTS_DIR.iterdir() if p.suffix in {".R", ".r", ".txt", ".Rmd"}])


def run_r_script(script_path: Path, args: List[str] = None, timeout: int = 300) -> subprocess.CompletedProcess:
    """Run an R script using Rscript. Raises FileNotFoundError if Rscript not found.

    Returns subprocess.CompletedProcess. Caller should handle exceptions.
    """
    if args is None:
        args = []
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise FileNotFoundError("Rscript not found in PATH. Install R and make Rscript available.")
    cmd = [rscript, str(script_path)] + args
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)


def script_mapping() -> dict:
    """Return a simple mapping of R scripts to a short description (manual).

    This is a lightweight 'link' between R assets and Python modules. You can
    expand it to include expected inputs/outputs.
    """
    return {
        "01.RefAltRepeatability.Visualization.R": "Visualization: draws arcs connecting variant to repeats",
        "02.RefAltRepeatability.Functions.R": "Core repeat finding functions (motif extraction, approximate matching)",
        "03.RefAltRepeatability.Function.R": "Pipeline wrapper: runs detection and writes outputs per-variant"
    }
