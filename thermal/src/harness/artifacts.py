# ABOUTME: Failure-only diagnostic emitter (ADR-0008). Dumps the error-vs-h
# table, the fitted rate, and a tricontourf plot of the error field at the
# finest mesh to tests/_artifacts/<test_id>/. Never writes on a passing test;
# the test wiring calls emit_failure_artifacts() from inside a try/except.

from __future__ import annotations

import logging
from pathlib import Path
from typing import Callable

import matplotlib

matplotlib.use("Agg")  # non-interactive backend for CI / pytest.
import matplotlib.pyplot as plt
import numpy as np

logger = logging.getLogger(__name__)


# ============================================================================
# FUNCTIONS
# ============================================================================


def emit_failure_artifacts(
    artifact_dir: Path,
    study,
    exact_solution: Callable[[np.ndarray, np.ndarray], np.ndarray],
) -> None:
    """Write the error-vs-h table, fitted rates, and finest-mesh error plot.

    Idempotent: overwrites existing files in ``artifact_dir``.
    """
    artifact_dir.mkdir(parents=True, exist_ok=True)

    # --- error-vs-h table + fitted rates ------------------------------------
    table_path = artifact_dir / "errors.csv"
    with table_path.open("w", encoding="utf-8") as fh:
        fh.write("h,n_dofs,l2_error,h1_error\n")
        for lvl in study.levels:
            fh.write(
                f"{lvl.mesh_size:.6e},{lvl.n_dofs},"
                f"{lvl.l2_error:.6e},{lvl.h1_error:.6e}\n"
            )

    rates_path = artifact_dir / "rates.txt"
    rates_path.write_text(
        f"l2_rate {study.l2_rate:.6f}\nh1_rate {study.h1_rate:.6f}\n",
        encoding="utf-8",
    )

    # --- finest-mesh error field plot --------------------------------------
    finest = study.levels[-1]
    mesh = finest.basis.mesh
    nodal_exact = exact_solution(mesh.p[0], mesh.p[1])
    err = finest.solution - nodal_exact

    fig, ax = plt.subplots(figsize=(6, 5))
    tcf = ax.tricontourf(mesh.p[0], mesh.p[1], mesh.t.T, err, levels=21, cmap="RdBu_r")
    fig.colorbar(tcf, ax=ax, label="T_h - T_exact")
    ax.set_aspect("equal")
    ax.set_title(
        f"Error field, h={finest.mesh_size:.3g}, "
        f"L2={finest.l2_error:.2e}, H1={finest.h1_error:.2e}"
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    plot_path = artifact_dir / "error_field.png"
    fig.savefig(plot_path, dpi=120, bbox_inches="tight")
    plt.close(fig)

    logger.warning(
        "failure artifacts emitted at %s (csv=%s, plot=%s)",
        artifact_dir,
        table_path.name,
        plot_path.name,
    )
