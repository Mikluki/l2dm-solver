# ABOUTME: End-to-end pytest wiring Problem 3 (radially symmetric disk in
# larger disk) through the verification harness. Asserts the L^2 convergence
# rate and the kappa_1-independence acceptance signal: at the finest mesh,
# T_h(0.6, 0) for kappa_1 in {0.1, 1, 10} must have spread < 1% and each
# value within 5% of T_out(0.6) = q_0 R_0^2/(2 kappa_2) ln(R_out/0.6)
# ≈ 2.299e-3. Also asserts that varying kappa_1 does not add new .msh files
# to the cache (geometry-only key per brief § Decisions 5).

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from src.harness.artifacts import emit_failure_artifacts
from src.harness.study import run_refinement_study
from src.problems.problem_03_disk import (
    OUTER_PROBE_EXACT_T,
    OUTER_PROBE_R,
    Problem03Disk,
)
from src.solver.solve_scalar import solve_scalar


# ============================================================================
# CONSTANTS
# ============================================================================

# Rate thresholds per verification.md § Problem 3 acceptance and the
# brief's threshold-skimming guard (within 0.2 of the theoretical rate).
L2_RATE_MIN = 1.8
L2_RATE_THEORY = 2.0
RATE_WINDOW = 0.2

# kappa_1-sweep acceptance per brief § Decisions 9.
KAPPA_1_SWEEP = (0.1, 1.0, 10.0)
SWEEP_SPREAD_TOL = 0.01  # (max - min)/|mean| < 1%
SWEEP_BOUND_TOL = 0.05   # |T_probe - exact| / exact < 5%

# Probe coordinate: x=0.6 in the outer annulus (R_0 < 0.6 < R_out).
# T_out(0.6) is kappa_1-independent — the load-bearing structural signal.
_PROBE_POINT = np.array([[OUTER_PROBE_R], [0.0]])


# ============================================================================
# TESTS
# ============================================================================


def test_problem_03_converges(artifact_dir: Path, mesh_cache_dir: Path) -> None:
    """Convergence rate at the verification.md § Problem 3 default parameters."""
    problem = Problem03Disk()  # kappa_1 = 1.0 (verification.md default)

    study = run_refinement_study(problem, mesh_cache_dir=mesh_cache_dir)

    # Print convergence table to stdout (visible with -s).
    header = f"{'h':>10}  {'n_dofs':>8}  {'L2':>12}  {'H1':>12}"
    print(f"\n{header}")
    for lvl in study.levels:
        print(
            f"{lvl.mesh_size:10.4f}  {lvl.n_dofs:8d}  "
            f"{lvl.l2_error:12.4e}  {lvl.h1_error:12.4e}"
        )
    print(f"fitted L2 rate = {study.l2_rate:.4f}  H1 rate = {study.h1_rate:.4f}")

    try:
        assert (
            study.l2_rate >= L2_RATE_MIN
        ), f"L^2 rate {study.l2_rate:.3f} below floor {L2_RATE_MIN}"
        assert abs(study.l2_rate - L2_RATE_THEORY) <= RATE_WINDOW, (
            f"L^2 rate {study.l2_rate:.3f} more than {RATE_WINDOW} from "
            f"theoretical {L2_RATE_THEORY}"
        )
        # H^1 rate logged but not asserted (verification.md specifies only L^2).
        # No pin assertion: Problem 3 has Dirichlet, so pin_dof is None for
        # every level - the harness records this and the test trusts it.
        for lvl in study.levels:
            assert lvl.pin_dof is None, (
                f"pin_dof should be None when Dirichlet is present; "
                f"got {lvl.pin_dof} at h={lvl.mesh_size}"
            )
    except AssertionError:
        emit_failure_artifacts(artifact_dir, study, problem.exact_solution)
        raise


def test_problem_03_kappa1_independence(
    artifact_dir: Path, mesh_cache_dir: Path
) -> None:
    """Brief acceptance: T_h(0.6, 0) does not depend on kappa_1.

    Two prongs:
      (a) spread (max - min)/|mean| < 1% across the sweep.
      (b) each value within 5% of T_out(0.6) — catches a coherently wrong
          constant that prong (a) alone would miss.
    Both must hold simultaneously.
    """
    finest_h = min(Problem03Disk().mesh_sizes())
    probe_values: dict[float, float] = {}
    for k1 in KAPPA_1_SWEEP:
        problem = Problem03Disk(kappa_1=k1)
        result = solve_scalar(problem, finest_h, mesh_cache_dir=mesh_cache_dir)
        # Faithful FE interpolation via basis.probes — never a nearest-node readoff.
        probe = result.basis.probes(_PROBE_POINT)
        probe_values[k1] = float((probe @ result.solution)[0])

    values = np.array(list(probe_values.values()))
    spread = (values.max() - values.min()) / abs(values.mean())
    rel_err = np.abs(values - OUTER_PROBE_EXACT_T) / OUTER_PROBE_EXACT_T

    print(f"\nkappa_1-independence probe at r={OUTER_PROBE_R}:")
    print(f"  T_out(0.6) exact = {OUTER_PROBE_EXACT_T:.6e}")
    for k1, v in probe_values.items():
        print(f"  kappa_1={k1:5.1f}: T_h = {v:.6e}  rel_err = {abs(v - OUTER_PROBE_EXACT_T)/OUTER_PROBE_EXACT_T:.4f}")
    print(f"  spread = {spread:.4e}  (tol {SWEEP_SPREAD_TOL})")

    try:
        assert spread < SWEEP_SPREAD_TOL, (
            f"kappa_1 sweep spread {spread:.4f} >= {SWEEP_SPREAD_TOL}; "
            f"probes={probe_values}"
        )
        assert rel_err.max() < SWEEP_BOUND_TOL, (
            f"kappa_1 sweep probe(s) {probe_values} drift from exact "
            f"{OUTER_PROBE_EXACT_T:.4e} by > {SWEEP_BOUND_TOL} (max rel err "
            f"{rel_err.max():.4f})"
        )
    except AssertionError:
        artifact_dir.mkdir(parents=True, exist_ok=True)
        (artifact_dir / "kappa1_sweep.txt").write_text(
            "kappa_1\tT_probe\n"
            + "\n".join(f"{k}\t{v:.6e}" for k, v in probe_values.items())
            + f"\nspread={spread:.6e}\nexact={OUTER_PROBE_EXACT_T:.6e}\n",
            encoding="utf-8",
        )
        raise


def test_problem_03_kappa1_sweep_hits_cache(mesh_cache_dir: Path) -> None:
    """Brief acceptance: varying kappa_1 must not add new .msh files.

    kappa_1 is not a geometry parameter; the ADR-0007 cache key depends only
    on (geometry_name, R_inner, R_outer, mesh_size). Sweeping kappa_1 at fixed
    h therefore produces zero cache misses after the first instantiation.
    """
    finest_h = min(Problem03Disk().mesh_sizes())
    # Warm the cache with one kappa_1 value.
    solve_scalar(Problem03Disk(kappa_1=KAPPA_1_SWEEP[0]), finest_h, mesh_cache_dir=mesh_cache_dir)
    before = {p.name: p.stat().st_mtime_ns for p in mesh_cache_dir.glob("disk_in_disk_*.msh")}
    # Sweep the remaining values; none should touch the .msh cache.
    for k1 in KAPPA_1_SWEEP[1:]:
        solve_scalar(Problem03Disk(kappa_1=k1), finest_h, mesh_cache_dir=mesh_cache_dir)
    after = {p.name: p.stat().st_mtime_ns for p in mesh_cache_dir.glob("disk_in_disk_*.msh")}
    assert after == before, (
        f"kappa_1 sweep produced cache misses: before={sorted(before)}, "
        f"after={sorted(after)} (kappa_1 leaking into cache key?)"
    )
