# ABOUTME: Direct harness-norm regression tests that do not involve the solver,
# gmsh, or convergence fitting. These keep norm bugs from hiding behind rates.

from __future__ import annotations

import numpy as np
from skfem import Basis, ElementTriP1, MeshTri

from src.harness.norms import l2_error


def test_l2_error_uses_exact_solution_at_quadrature_points() -> None:
    """Nodal interpolation of a quadratic must not have zero L2 error."""
    mesh = MeshTri.init_tensor(
        np.linspace(0.0, 1.0, 3),
        np.linspace(0.0, 1.0, 3),
    )
    basis = Basis(mesh, ElementTriP1())

    def exact(x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return x**2

    solution = exact(basis.mesh.p[0], basis.mesh.p[1])

    err = l2_error(solution, basis, exact)

    assert 1.0e-2 < err < 1.0e-1
