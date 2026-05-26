# ABOUTME: Regression test for the failure-only artifact emitter (ADR-0008,
# brief 0001 acceptance #5). Builds a synthetic StudyResult, calls
# emit_failure_artifacts, and asserts the three-file diagnostic bundle
# materialises. Keeps the failure path from rotting silently.

from __future__ import annotations

from pathlib import Path

import numpy as np
from skfem import Basis, ElementTriP1, MeshTri

from src.harness.artifacts import emit_failure_artifacts
from src.harness.study import LevelResult, StudyResult


def test_emit_failure_artifacts_writes_bundle(tmp_path: Path) -> None:
    """All three artifacts (csv, rates, plot) must land on the failure path."""
    mesh = MeshTri()
    basis = Basis(mesh, ElementTriP1())
    level = LevelResult(
        mesh_size=1.0,
        n_dofs=basis.N,
        pin_dof=0,
        l2_error=1.0e-2,
        h1_error=1.0e-1,
        solution=np.zeros(basis.N),
        basis=basis,
    )
    study = StudyResult(levels=[level], l2_rate=0.5, h1_rate=0.25)

    emit_failure_artifacts(tmp_path, study, lambda x, _y: np.zeros_like(x))

    assert (tmp_path / "errors.csv").exists()
    assert (tmp_path / "rates.txt").exists()
    assert (tmp_path / "error_field.png").exists()
