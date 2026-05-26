# ABOUTME: solve_scalar(problem, mesh_size) - the Part 1 scalar PDE solver.
# Materialises the cached mesh, builds a P1 basis, assembles the stiffness
# and load forms, applies the node-pin from Problem.pin_point() per ADR-0005,
# and runs a direct sparse solve (ADR-0006). For Problem 1 specifically:
# single subdomain, uniform kappa=1, smooth source - no P0 indirection
# needed yet (deferred to Problem 2's submission).

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import meshio
import numpy as np
from skfem import (
    Basis,
    BilinearForm,
    ElementTriP1,
    LinearForm,
    condense,
    solve,
)
from skfem.helpers import dot, grad
from skfem.io.meshio import from_meshio

from src.geometry.unit_square import (
    UnitSquareSpec,
    materialise as materialise_unit_square,
)
from src.solver.result import SolverResult

logger = logging.getLogger(__name__)


# ============================================================================
# CONFIG
# ============================================================================

# Default mesh cache root. Tests override via the public solve_scalar() arg.
_DEFAULT_CACHE = Path(__file__).resolve().parents[2] / "tests" / "_mesh_cache"


# ============================================================================
# FORMS
# ============================================================================


@BilinearForm
def _stiffness(u, v, _w):
    # Single subdomain, kappa = 1 - the P0-field-per-tag indirection is the
    # Problem 2 pattern (architecture.md § Coefficient handling) and is
    # explicitly deferred for this submission.
    return dot(grad(u), grad(v))


def _make_load(source):
    """Build a LinearForm that calls ``source(x, y)`` at quadrature points.

    Sampling at quadrature points - rather than projecting nodal values via
    ``basis.interpolate`` - avoids an O(h^2) projection error that otherwise
    contaminates the L^2 rate measurement. The "right rate, wrong constant"
    failure mode in architecture.md § Coefficient handling applies to source
    fields just as much as to kappa.
    """

    @LinearForm
    def _load(v, w):
        return source(w.x[0], w.x[1]) * v

    return _load


# ============================================================================
# HELPERS
# ============================================================================


def _materialise_geometry(spec: Any, cache_dir: Path) -> Path:
    """Dispatch the geometry-spec object to the right cache-backed builder.

    Kept tiny on purpose: Problem 1 only needs UnitSquareSpec; further shapes
    extend this dispatch in their own submissions.
    """
    if isinstance(spec, UnitSquareSpec):
        return materialise_unit_square(spec, cache_dir)
    raise TypeError(f"unsupported geometry spec: {type(spec).__name__}")


def _load_mesh(msh_path: Path):
    """Read a `.msh` and convert to a scikit-fem MeshTri, preserving tags."""
    return from_meshio(meshio.read(msh_path))


def _nearest_node_dof(basis: Basis, point: tuple[float, float]) -> int:
    """Return the global DOF index of the P1 node nearest ``point``.

    For P1 elements there is one DOF per node so the node index is the DOF
    index. Stable choice in ties: ``argmin`` returns the lowest index.
    """
    px, py = point
    coords = basis.mesh.p
    d2 = (coords[0] - px) ** 2 + (coords[1] - py) ** 2
    return int(np.argmin(d2))


# ============================================================================
# SOLVER
# ============================================================================


def solve_scalar(
    problem,
    mesh_size: float,
    *,
    mesh_cache_dir: Path | None = None,
) -> SolverResult:
    """Solve the Part 1 scalar PDE for ``problem`` at ``mesh_size``.

    The pipeline mirrors architecture.md § Pipeline: geometry -> cached mesh
    -> P1 basis -> assembly -> BC/pin -> direct sparse solve.
    """
    cache_dir = mesh_cache_dir if mesh_cache_dir is not None else _DEFAULT_CACHE

    # --- Geometry / mesh ----------------------------------------------------
    geometry_builder = problem.geometry()
    spec = geometry_builder(mesh_size)
    msh_path = _materialise_geometry(spec, cache_dir)
    mesh = _load_mesh(msh_path)
    logger.debug("loaded mesh: %d nodes, %d elements", mesh.p.shape[1], mesh.t.shape[1])

    # --- Basis + assembly ---------------------------------------------------
    basis = Basis(mesh, ElementTriP1())
    A = _stiffness.assemble(basis)
    # Source evaluated directly at quadrature points - see _make_load.
    b = _make_load(problem.source).assemble(basis)

    # --- BCs ---------------------------------------------------------------
    # For Problem 1 there are no Dirichlet tags; zero-flux Neumann is the
    # natural BC and contributes nothing to b. When future problems introduce
    # Dirichlet, extend boundary_conditions() handling here.
    bcs = problem.boundary_conditions()
    if bcs:
        raise NotImplementedError(
            "Dirichlet BC handling is deferred to Problem 2's submission; "
            f"problem returned tags {list(bcs)}"
        )

    # --- Nullspace pin ------------------------------------------------------
    # architecture.md § Nullspace handling: caller declares pin_point(); the
    # solver pins the nearest node to exact_solution(*pin_point). No default,
    # no fallback - returning None means "Dirichlet present, no pin needed".
    pin_point = problem.pin_point()
    if pin_point is None:
        raise ValueError(
            "Problem 1 has no Dirichlet BCs and requires a pin_point(); got None"
        )
    pin_dof = _nearest_node_dof(basis, pin_point)
    pin_value = float(
        problem.exact_solution(np.array([pin_point[0]]), np.array([pin_point[1]]))[0]
    )
    x = np.zeros(basis.N)
    x[pin_dof] = pin_value
    logger.debug("pinning DOF %d to T=%.6f at %s", pin_dof, pin_value, pin_point)

    # --- Solve --------------------------------------------------------------
    solution = solve(*condense(A, b, x=x, D=np.array([pin_dof])))

    return SolverResult(
        solution=solution,
        basis=basis,
        mesh=mesh,
        pin_dof=pin_dof,
    )
