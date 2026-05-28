# ABOUTME: End-to-end pytest wiring Problem 2 (piecewise-constant kappa, 1D
# slab) through the verification harness. Asserts the L^2 convergence rate
# and the kappa_2-independence acceptance signal: at the finest mesh,
# T_h(0.75, h/2) for kappa_2 in {10, 100, 1000} must have spread < 1% and
# each value within 5% of q_0/(8 kappa_1). Also asserts that varying kappa_2
# does not add new .msh files to the cache (geometry-only key).

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from src.harness.artifacts import emit_failure_artifacts
from src.harness.study import run_refinement_study
from src.problems.problem_02_slab import (
    RIGHT_REGION_EXACT_T,
    SLAB_HEIGHT,
    Problem02Slab,
)
from src.solver.solve_scalar import solve_scalar


# ============================================================================
# CONSTANTS
# ============================================================================

# Rate thresholds per verification.md § Problem 2 acceptance ("observed rate
# >= 1.8") and the brief's threshold-skimming guard ("within 0.2 of 2.0").
L2_RATE_MIN = 1.8
L2_RATE_THEORY = 2.0
RATE_WINDOW = 0.2

# kappa_2-sweep acceptance per submission 0003 § Decisions 3.
KAPPA_2_SWEEP = (10.0, 100.0, 1000.0)
SWEEP_SPREAD_TOL = 0.01  # (max - min)/|mean| < 1%
SWEEP_BOUND_TOL = 0.05  # |T_probe - exact| / exact < 5%

# Probe coordinate: x=0.75 sits inside the right subdomain where the exact
# solution is the constant q_0/(8 kappa_1) = 0.125; y is mid-slab.
_PROBE_POINT = np.array([[0.75], [SLAB_HEIGHT / 2.0]])


# ============================================================================
# TESTS
# ============================================================================


def test_problem_02_converges(artifact_dir: Path, mesh_cache_dir: Path) -> None:
    """Convergence rate at the verification.md § Problem 2 default kappa_2."""
    problem = Problem02Slab()  # kappa_2 = 100 (verification.md default)

    study = run_refinement_study(problem, mesh_cache_dir=mesh_cache_dir)

    try:
        # --- Convergence rate (L^2 only per verification.md) -------------
        assert (
            study.l2_rate >= L2_RATE_MIN
        ), f"L^2 rate {study.l2_rate:.3f} below floor {L2_RATE_MIN}"
        assert abs(study.l2_rate - L2_RATE_THEORY) <= RATE_WINDOW, (
            f"L^2 rate {study.l2_rate:.3f} more than {RATE_WINDOW} from "
            f"theoretical {L2_RATE_THEORY}"
        )
        # H^1 rate logged but not asserted (verification.md specifies only L^2).
        # No pin assertion: Problem 2 has Dirichlet, so pin_dof is None for
        # every level - the harness records this and the test trusts it.
        for lvl in study.levels:
            assert lvl.pin_dof is None, (
                f"pin_dof should be None when Dirichlet is present; "
                f"got {lvl.pin_dof} at h={lvl.mesh_size}"
            )
    except AssertionError:
        emit_failure_artifacts(artifact_dir, study, problem.exact_solution)
        raise


def test_problem_02_kappa2_independence(
    artifact_dir: Path, mesh_cache_dir: Path
) -> None:
    """Brief acceptance #3: T_h(0.75, h/2) does not depend on kappa_2.

    Two prongs:
      (a) spread (max - min)/|mean| < 1% across the sweep - the literal
          verification.md signal.
      (b) each value within 5% of the exact constant q_0/(8 kappa_1) - catches
          a sign-flipped or otherwise coherently wrong constant that prong (a)
          alone would miss.
    Both must hold simultaneously.
    """
    finest_h = min(Problem02Slab().mesh_sizes())
    probe_values: dict[float, float] = {}
    for k2 in KAPPA_2_SWEEP:
        problem = Problem02Slab(kappa_2=k2)
        result = solve_scalar(problem, finest_h, mesh_cache_dir=mesh_cache_dir)
        # Faithful FE interpolation via basis.probes - never a nearest-node
        # readoff. probes returns a sparse matrix mapping solution -> point
        # values; the brief flags nearest-node as the foot-gun to avoid.
        probe = result.basis.probes(_PROBE_POINT)
        probe_values[k2] = float((probe @ result.solution)[0])

    values = np.array(list(probe_values.values()))
    spread = (values.max() - values.min()) / abs(values.mean())
    rel_err = np.abs(values - RIGHT_REGION_EXACT_T) / RIGHT_REGION_EXACT_T

    try:
        assert spread < SWEEP_SPREAD_TOL, (
            f"kappa_2 sweep spread {spread:.4f} >= {SWEEP_SPREAD_TOL}; "
            f"probes={probe_values}"
        )
        assert rel_err.max() < SWEEP_BOUND_TOL, (
            f"kappa_2 sweep probe(s) {probe_values} drift from exact "
            f"{RIGHT_REGION_EXACT_T} by > {SWEEP_BOUND_TOL} (max rel err "
            f"{rel_err.max():.4f})"
        )
    except AssertionError:
        # No full study to dump on this test - write a minimal note so the
        # bundle directory still materialises on failure.
        artifact_dir.mkdir(parents=True, exist_ok=True)
        (artifact_dir / "kappa2_sweep.txt").write_text(
            "kappa_2\tT_probe\n"
            + "\n".join(f"{k}\t{v:.6e}" for k, v in probe_values.items())
            + f"\nspread={spread:.6e}\nexact={RIGHT_REGION_EXACT_T}\n",
            encoding="utf-8",
        )
        raise


def test_problem_02_kappa2_sweep_hits_cache(mesh_cache_dir: Path) -> None:
    """Brief acceptance #6: varying kappa_2 must not add new .msh files.

    kappa is not a geometry parameter, so the ADR-0007 cache key must depend
    only on (geometry_name, mesh_size, height). Sweeping kappa_2 across every
    declared mesh size therefore produces zero cache misses after the first
    pass. Walking the full mesh_sizes() list (not just the finest) guards
    against a partial leak that would only fire at one resolution.
    """
    sizes = list(Problem02Slab().mesh_sizes())
    # Warm the cache: one full mesh-size pass at the first kappa_2.
    for h in sizes:
        solve_scalar(
            Problem02Slab(kappa_2=KAPPA_2_SWEEP[0]),
            h,
            mesh_cache_dir=mesh_cache_dir,
        )
    before = {p.name: p.stat().st_mtime_ns for p in mesh_cache_dir.glob("*.msh")}
    # Sweep remaining kappa_2 across every mesh size; none should touch cache.
    for k2 in KAPPA_2_SWEEP[1:]:
        for h in sizes:
            solve_scalar(
                Problem02Slab(kappa_2=k2), h, mesh_cache_dir=mesh_cache_dir
            )
    after = {p.name: p.stat().st_mtime_ns for p in mesh_cache_dir.glob("*.msh")}
    assert after == before, (
        f"kappa_2 sweep produced cache misses: before={sorted(before)}, "
        f"after={sorted(after)} (kappa leaking into cache key?)"
    )
