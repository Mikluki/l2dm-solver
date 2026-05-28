# ABOUTME: End-to-end pytest for Problem 5 (L-shape reentrant corner). The
# harness assertion is INVERTED: a fitted L^2 rate in [1.2, 1.5] is the success
# window, both bounds load-bearing. An auxiliary test verifies the same check
# helper rejects a synthetic smooth-rate StudyResult (submission 0007 §
# Acceptance 4 - same helper as the main test). Acceptance 6 corner-node check
# and Acceptance 7 cache-hit check live alongside.

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from src.harness.artifacts import emit_failure_artifacts
from src.harness.study import StudyResult, run_refinement_study
from src.problems.problem_05_lshape import Problem05LShape
from src.solver.solve_scalar import solve_scalar


# ============================================================================
# CONSTANTS
# ============================================================================

# Inverted-assertion window per submission 0007 § Decisions 2 / Acceptance 2.
# rate > UPPER: rate fitter dishonest (smearing the corner singularity).
# rate < LOWER: rate fitter broken, or finest mesh not in asymptotic regime.
RATE_LOWER = 1.2
RATE_UPPER = 1.5

# Reentrant corner. The OCC point at (1/2, 1/2) is a guaranteed mesh node at
# every refinement; the assert in test_problem_05_corner_is_mesh_node is sharp.
_CORNER = (0.5, 0.5)
_CORNER_TOL = 1e-12


# ============================================================================
# SHARED CHECK HELPER
# ============================================================================


def check_inverted_rate_window(l2_rate: float) -> None:
    """Assert ``l2_rate`` lies in the inverted-success window [LOWER, UPPER].

    Both bounds are load-bearing (submission 0007 § Acceptance 2):
        rate > RATE_UPPER -> rate fitter is dishonest (e.g. the L^2 norm
            regressed, or the coarsest meshes smear the singularity).
        rate < RATE_LOWER -> rate fitter broken or finest mesh too coarse
            to be in the asymptotic regime.

    This helper is the *single* check shared by the main convergence test and
    the forced-failure test that exercises the inversion mechanism itself
    (Acceptance 4). Duplicating the inequality inline would let copy-drift
    defeat the inversion proof.
    """
    if l2_rate < RATE_LOWER:
        raise AssertionError(
            f"L^2 rate {l2_rate:.4f} below inverted window lower bound "
            f"{RATE_LOWER}; either the rate fitter is broken, or the finest "
            f"mesh is not yet in the asymptotic regime for the 4/3 mode."
        )
    if l2_rate > RATE_UPPER:
        raise AssertionError(
            f"L^2 rate {l2_rate:.4f} above inverted window upper bound "
            f"{RATE_UPPER}; the rate fitter is reporting smooth convergence "
            f"on a corner-singular solution (norm regressed or mesh too "
            f"coarse to resolve the corner)."
        )


# ============================================================================
# TESTS
# ============================================================================


def test_problem_05_corner_is_mesh_node(mesh_cache_dir: Path) -> None:
    """Submission 0007 § Acceptance 6: the reentrant corner is a mesh node.

    The OCC builder adds (1/2, 1/2) as an explicit point, so the nearest mesh
    node is exact at every refinement level. The singular solution is
    meaningless if the corner drifts even by a single grid spacing.
    """
    problem = Problem05LShape()
    for h in problem.mesh_sizes():
        sr = solve_scalar(problem, h, mesh_cache_dir=mesh_cache_dir)
        coords = sr.basis.mesh.p
        d_inf = np.maximum(
            np.abs(coords[0] - _CORNER[0]), np.abs(coords[1] - _CORNER[1])
        )
        d_min = float(d_inf.min())
        assert d_min < _CORNER_TOL, (
            f"reentrant-corner node missing at h={h}; nearest node at "
            f"L-inf distance {d_min:.3e} > tol {_CORNER_TOL:.0e}"
        )


def test_problem_05_inverted_rate_window(
    artifact_dir: Path, mesh_cache_dir: Path
) -> None:
    """Submission 0007 § Acceptance 2: 1.2 <= L^2 rate <= 1.5 over >=5 levels.

    Both bounds load-bearing; the artifact emitter records the full error
    table and the finest-mesh error field on either failure. H^1 rate is
    logged for inspection but not asserted (verification.md is silent on the
    expected H^1 rate at the L-shape).
    """
    problem = Problem05LShape()
    study = run_refinement_study(problem, mesh_cache_dir=mesh_cache_dir)
    assert len(study.levels) >= 5, (
        f"Problem 5 requires >=5 refinement levels; got {len(study.levels)}"
    )

    header = f"{'h':>10}  {'n_dofs':>8}  {'L2':>12}  {'H1':>12}"
    print(f"\n{header}")
    for lvl in study.levels:
        print(
            f"{lvl.mesh_size:10.4f}  {lvl.n_dofs:8d}  "
            f"{lvl.l2_error:12.4e}  {lvl.h1_error:12.4e}"
        )
    print(
        f"fitted L2 rate = {study.l2_rate:.4f}  "
        f"H1 rate = {study.h1_rate:.4f}  (H1 logged, not asserted)"
    )

    try:
        # pin_dof must be None at every level (every edge is Dirichlet).
        for lvl in study.levels:
            assert lvl.pin_dof is None, (
                f"pin_dof should be None when Dirichlet present; got "
                f"{lvl.pin_dof} at h={lvl.mesh_size}"
            )
        check_inverted_rate_window(study.l2_rate)
    except AssertionError:
        emit_failure_artifacts(artifact_dir, study, problem.exact_solution)
        raise


def test_problem_05_rate_window_helper_rejects_smooth_rate() -> None:
    """Submission 0007 § Acceptance 4: forced failure on the *mechanism*.

    Build a synthetic ``StudyResult`` with ``l2_rate ~= 2.0`` (the smooth-
    solution rate Problem 1 hits) and assert the SAME helper used by the
    main test rejects it. This catches bugs in the inversion logic itself,
    which the natural Problem 5 run cannot - a passing Problem 5 only
    confirms real corner singularities are detected, not that the inversion
    inequality is the right shape.

    A second case checks the lower-bound branch with rate ~= 0 (the value
    you'd see if the callable-Dirichlet wiring silently zeroed every outer
    edge - Acceptance 5's forced-failure scenario).
    """
    # Synthetic smooth-rate StudyResult (rate >= 2 -> upper-bound branch).
    smooth_study = StudyResult(levels=[], l2_rate=2.0, h1_rate=1.0)
    with pytest.raises(AssertionError, match="above inverted window upper bound"):
        check_inverted_rate_window(smooth_study.l2_rate)

    # Synthetic non-converging StudyResult (rate ~= 0 -> lower-bound branch).
    flat_study = StudyResult(levels=[], l2_rate=0.05, h1_rate=0.0)
    with pytest.raises(AssertionError, match="below inverted window lower bound"):
        check_inverted_rate_window(flat_study.l2_rate)


def test_problem_05_mesh_cache_hits(mesh_cache_dir: Path) -> None:
    """Submission 0007 § Acceptance 7: re-running the study adds zero .msh files.

    The geometry cache key (ADR-0007) is (geometry_name, mesh_size); a warm
    second pass must hit every level. No silent re-mesh.
    """
    problem = Problem05LShape()
    sizes = problem.mesh_sizes()
    # Cold pass to ensure every level is materialised.
    for h in sizes:
        solve_scalar(problem, h, mesh_cache_dir=mesh_cache_dir)
    before = {
        p.name: p.stat().st_mtime_ns
        for p in mesh_cache_dir.glob("l_shape_*.msh")
    }
    # Warm pass - must reuse every cached .msh.
    for h in sizes:
        solve_scalar(problem, h, mesh_cache_dir=mesh_cache_dir)
    after = {
        p.name: p.stat().st_mtime_ns
        for p in mesh_cache_dir.glob("l_shape_*.msh")
    }
    assert after == before, (
        f"Warm re-run produced cache misses: before={sorted(before)}, "
        f"after={sorted(after)} - cache key may include non-geometry params."
    )
