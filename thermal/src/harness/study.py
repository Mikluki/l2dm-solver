# ABOUTME: Refinement-study driver. Runs the solver at each Problem-declared
# mesh size, computes L^2 and H^1-seminorm errors against the exact solution,
# and fits convergence rates by least-squares on log(error) vs log(h). On
# assertion failure inside a test, the result object carries everything the
# artifact emitter needs to dump a diagnostic bundle.

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from src.harness.norms import h1_seminorm_error, l2_error
from src.solver.solve_scalar import solve_scalar

logger = logging.getLogger(__name__)


# ============================================================================
# DATACLASSES
# ============================================================================


@dataclass(frozen=True)
class LevelResult:
    """One row of the refinement study."""

    mesh_size: float
    n_dofs: int
    pin_dof: int
    l2_error: float
    h1_error: float
    solution: np.ndarray
    basis: Any


@dataclass
class StudyResult:
    """Aggregate of all levels plus fitted rates."""

    levels: list[LevelResult] = field(default_factory=list)
    l2_rate: float = float("nan")
    h1_rate: float = float("nan")


# ============================================================================
# FUNCTIONS
# ============================================================================


def _fit_rate(h: np.ndarray, err: np.ndarray) -> float:
    """Least-squares slope of log(err) vs log(h). Convergence rate = slope."""
    slope = float(np.polyfit(np.log(h), np.log(err), 1)[0])
    return slope


def run_refinement_study(
    problem,
    *,
    mesh_cache_dir: Path | None = None,
) -> StudyResult:
    """Run the solver at each Problem-declared mesh size and fit rates.

    Always returns a :class:`StudyResult`. Artifact emission is *not* the
    study's responsibility - the test wraps its assertions in a try/except
    and invokes the artifact emitter on failure (see ``tests/test_problem_01``).
    """
    sizes = list(problem.mesh_sizes())
    if len(sizes) < 3:
        raise ValueError(
            f"refinement study needs >=3 mesh sizes (verification.md § "
            f"Harness requirements); got {sizes}"
        )

    levels: list[LevelResult] = []
    for h in sizes:
        sr = solve_scalar(problem, h, mesh_cache_dir=mesh_cache_dir)
        err_l2 = l2_error(sr.solution, sr.basis, problem.exact_solution)
        # exact_gradient is not a Protocol method; problems with smooth
        # analytic solutions expose it on the concrete class.
        if not hasattr(problem, "exact_gradient"):
            raise NotImplementedError(
                "H^1 seminorm requires the concrete problem class to expose "
                "an `exact_gradient` callable; problems lacking one need a "
                "future projection-based fallback."
            )
        err_h1 = h1_seminorm_error(sr.solution, sr.basis, problem.exact_gradient)
        logger.info(
            "h=%.4f n_dofs=%d L2=%.3e H1=%.3e pin_dof=%d",
            h,
            sr.basis.N,
            err_l2,
            err_h1,
            sr.pin_dof,
        )
        levels.append(
            LevelResult(
                mesh_size=h,
                n_dofs=sr.basis.N,
                pin_dof=sr.pin_dof,
                l2_error=err_l2,
                h1_error=err_h1,
                solution=sr.solution,
                basis=sr.basis,
            )
        )

    h_arr = np.array([lvl.mesh_size for lvl in levels])
    l2_arr = np.array([lvl.l2_error for lvl in levels])
    h1_arr = np.array([lvl.h1_error for lvl in levels])

    study = StudyResult(
        levels=levels,
        l2_rate=_fit_rate(h_arr, l2_arr),
        h1_rate=_fit_rate(h_arr, h1_arr),
    )
    logger.info("fitted rates: L2=%.3f H1=%.3f", study.l2_rate, study.h1_rate)
    return study
