# ABOUTME: End-to-end pytest wiring Problem 1 through the verification harness.
# Asserts the L^2/H^1 convergence rates, the pinned-DOF reproducibility, and
# the small mean-shift between computed and exact T at the finest mesh.
# Diagnostic artifacts are emitted only on assertion failure (ADR-0008).

from __future__ import annotations

from pathlib import Path

import pytest
from skfem import Functional

from src.harness.artifacts import emit_failure_artifacts
from src.harness.study import run_refinement_study
from src.problems.problem_01_manufactured import Problem01Manufactured
from src.solver.solve_scalar import solve_scalar


# ============================================================================
# FORMS
# ============================================================================


@Functional
def _domain_integral(w):
    # Mean shift check (see docs/open-questions.md): integrate the discrete
    # solution over the domain and compare against the analytic mean.
    return w["u"]

# ============================================================================
# CONSTANTS
# ============================================================================

# Rate thresholds per verification.md § Problem 1 acceptance and the
# submission brief. The "within 0.2 of the theoretical value" guard catches
# threshold-skimming (a rate that scrapes 1.8 should be investigated).
L2_RATE_MIN = 1.8
L2_RATE_THEORY = 2.0
H1_RATE_MIN = 0.9
H1_RATE_THEORY = 1.0
RATE_WINDOW = 0.2


# ============================================================================
# TESTS
# ============================================================================


def test_problem_01_converges(artifact_dir: Path, mesh_cache_dir: Path) -> None:
    """Problem 1 must converge at the predicted P1 rates and respect the pin."""
    problem = Problem01Manufactured()

    study = run_refinement_study(problem, mesh_cache_dir=mesh_cache_dir)

    try:
        # --- Convergence rates --------------------------------------------
        assert (
            study.l2_rate >= L2_RATE_MIN
        ), f"L^2 rate {study.l2_rate:.3f} below floor {L2_RATE_MIN}"
        assert abs(study.l2_rate - L2_RATE_THEORY) <= RATE_WINDOW, (
            f"L^2 rate {study.l2_rate:.3f} more than {RATE_WINDOW} from "
            f"theoretical {L2_RATE_THEORY}"
        )
        assert (
            study.h1_rate >= H1_RATE_MIN
        ), f"H^1 rate {study.h1_rate:.3f} below floor {H1_RATE_MIN}"
        assert abs(study.h1_rate - H1_RATE_THEORY) <= RATE_WINDOW, (
            f"H^1 rate {study.h1_rate:.3f} more than {RATE_WINDOW} from "
            f"theoretical {H1_RATE_THEORY}"
        )

        # --- Pinned DOF index is reproducible across refinements -----------
        pin_indices = [r.pin_dof for r in study.levels]
        assert all(
            idx == pin_indices[0] for idx in pin_indices
        ), f"Pinned DOF drifted across refinements: {pin_indices}"
        # And the pin is at (0, 0) per Problem 1's declared pin_point.
        finest = study.levels[-1]
        pin_x = finest.basis.mesh.p[0, finest.pin_dof]
        pin_y = finest.basis.mesh.p[1, finest.pin_dof]
        assert pin_x == pytest.approx(0.0, abs=1e-12)
        assert pin_y == pytest.approx(0.0, abs=1e-12)

        # --- Mean shift: assembly is unbiased -------------------------------
        # Acceptance #4 of the brief, reinterpreted per docs/open-questions.md.
        # mean(T_exact) = 0 on [0,1]^2 for cos(pi x) cos(pi y); the discrete
        # mean integral decays at the FE error rate. Calibrating the tolerance
        # against the finest L^2 error catches sign-flipped contributions that
        # the rate measurement alone cannot.
        mean_h = float(
            _domain_integral.assemble(
                finest.basis, u=finest.basis.interpolate(finest.solution)
            )
        )
        mean_tol = 5.0 * finest.l2_error
        assert abs(mean_h) < mean_tol, (
            f"mean(T_h)={mean_h:.3e} exceeds FE-scale bound {mean_tol:.3e}; "
            f"assembly may be biased or pin location inconsistent"
        )
    except AssertionError:
        # Failure-only artifact emission (ADR-0008, brief acceptance #5).
        emit_failure_artifacts(artifact_dir, study, problem.exact_solution)
        raise


def test_problem_01_solver_smoke(mesh_cache_dir: Path) -> None:
    """Smoke check the solver returns a SolverResult with consistent shapes."""
    problem = Problem01Manufactured()
    result = solve_scalar(
        problem,
        mesh_size=problem.mesh_sizes()[0],
        mesh_cache_dir=mesh_cache_dir,
    )
    assert result.solution.shape == (result.basis.mesh.p.shape[1],)
    assert result.basis.mesh.p.shape[0] == 2


def test_mesh_cache_is_reused(mesh_cache_dir: Path) -> None:
    """Brief acceptance #6: re-running must hit the cache for every mesh size.

    The first pass warms the cache; the second pass must produce no new files
    and no mtime changes. A drift here means ``cache_key`` is non-deterministic
    (object ids, wall-clock, dict iteration order leaked into the hash).
    """
    problem = Problem01Manufactured()
    sizes = problem.mesh_sizes()
    for h in sizes:
        solve_scalar(problem, h, mesh_cache_dir=mesh_cache_dir)
    before = {p.name: p.stat().st_mtime_ns for p in mesh_cache_dir.glob("*.msh")}
    for h in sizes:
        solve_scalar(problem, h, mesh_cache_dir=mesh_cache_dir)
    after = {p.name: p.stat().st_mtime_ns for p in mesh_cache_dir.glob("*.msh")}
    assert after == before, (
        f"mesh cache miss on re-run: before={sorted(before)}, after={sorted(after)}"
    )
