# ABOUTME: SolverResult dataclass returned by solve_scalar(). Carries the
# nodal solution vector, the Basis it was computed in, the mesh, and the
# pinned-DOF index (for nullspace-handled problems). No methods - this is the
# interchange struct that survives unchanged into Part 2's integral form.

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

# ============================================================================
# DATACLASS
# ============================================================================


@dataclass(frozen=True)
class SolverResult:
    """Discrete solution + the scikit-fem context that produced it.

    Attributes:
        solution: nodal DOF values, shape ``(n_dofs,)``.
        basis: the scikit-fem Basis the solution lives in.
        mesh: the MeshTri the basis was built on (convenience handle).
        pin_dof: the global DOF index pinned to the exact value, or ``None``
            if the problem had Dirichlet BCs and required no pin.
    """

    solution: np.ndarray
    basis: Any
    mesh: Any
    pin_dof: int | None
