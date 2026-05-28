# ABOUTME: solve_scalar(problem, mesh_size) - the Part 1 scalar PDE solver.
# Materialises the cached mesh, builds a P1 basis and a derived P0 basis for
# the per-element kappa field (architecture.md § Coefficient handling),
# assembles the stiffness and load forms, applies Dirichlet DOFs from
# problem.boundary_conditions() and/or the node-pin from problem.pin_point()
# (ADR-0005, nullspace handling rule 1: at least one Dirichlet => no pin), and
# runs a direct sparse solve (ADR-0006). Subdomain assignment is by gmsh
# physical-surface name propagated via meshio (ADR-0003) - never by
# element-centroid coordinate.

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import meshio
import numpy as np
from skfem import (
    Basis,
    BilinearForm,
    ElementTriP0,
    ElementTriP1,
    LinearForm,
    condense,
    solve,
)
from skfem.helpers import dot, grad
from skfem.io.meshio import from_meshio

from src.geometry.disk_in_disk import (
    DiskInDiskSpec,
    materialise as materialise_disk_in_disk,
)
from src.geometry.rectangle_split import (
    RectangleSplitSpec,
    materialise as materialise_rectangle_split,
)
from src.geometry.unit_square import (
    UnitSquareSpec,
    materialise as materialise_unit_square,
)
from src.problems.protocol import DirichletBC
from src.solver.result import SolverResult

logger = logging.getLogger(__name__)


# ============================================================================
# CONFIG
# ============================================================================

# Default mesh cache root. Tests override via the public solve_scalar() arg.
_DEFAULT_CACHE = Path(__file__).resolve().parents[2] / "tests" / "_mesh_cache"

# Subdomain names emitted by gmsh+meshio that are not real material regions
# (e.g. the auxiliary group meshio synthesises for line-bounded entities).
# We filter these out before building the per-element kappa table so the
# Problem.kappa() lookup is never asked for a name it doesn't know about.
_RESERVED_SUBDOMAIN_NAMES = {"gmsh:bounding_entities"}


# ============================================================================
# FORMS
# ============================================================================


@BilinearForm
def _stiffness_kappa(u, v, w):
    """Per-element kappa via the w.kappa indirection (P0 field).

    The form sees a different kappa value in each element; passing the raw
    P0 numpy array (rather than basis_p0.interpolate(...)) would silently
    broadcast and give "right rate, wrong constant" (architecture.md §
    Coefficient handling). The Problem-2 acceptance check exists to catch
    exactly that footgun.
    """
    return w.kappa * dot(grad(u), grad(v))


def _make_load(source):
    """Build a LinearForm that calls ``source(x, y)`` at quadrature points.

    Sampling at quadrature points - rather than projecting nodal values via
    ``basis.interpolate`` - avoids an O(h^2) projection error that otherwise
    contaminates the L^2 rate measurement. The "right rate, wrong constant"
    failure mode in architecture.md § Coefficient handling applies to source
    fields just as much as to kappa. For mesh-aligned discontinuous sources
    (Problem 2: q_0 inside x<1/2, 0 outside) the alignment guarantees each
    element's quadrature points all live on the same side of the interface,
    so per-element evaluation is exact.
    """

    @LinearForm
    def _load(v, w):
        return source(w.x[0], w.x[1]) * v

    return _load


# ============================================================================
# HELPERS
# ============================================================================


def _materialise_geometry(spec: Any, cache_dir: Path) -> Path:
    """Dispatch the geometry-spec object to the right cache-backed builder."""
    if isinstance(spec, UnitSquareSpec):
        return materialise_unit_square(spec, cache_dir)
    if isinstance(spec, RectangleSplitSpec):
        return materialise_rectangle_split(spec, cache_dir)
    if isinstance(spec, DiskInDiskSpec):
        return materialise_disk_in_disk(spec, cache_dir)
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


def _build_kappa_p0_values(mesh, problem) -> np.ndarray:
    """Per-element kappa array, indexed by mesh element id.

    Reads ``mesh.subdomains`` (a dict from physical-group name to element-
    index ndarray) and asks ``problem.kappa(name)`` for each material region.
    Reserved gmsh-internal names are skipped. Every element must end up in
    exactly one Problem-known subdomain; an uncovered element is a tagging
    bug (e.g. interface not synchronised with the surfaces) and we want it
    loud, not silent.
    """
    if not mesh.subdomains:
        # No tags exported - fall back to a single Problem.kappa(name) call
        # with the empty string so single-subdomain problems (Problem 1)
        # keep working without re-tagging their existing mesh.
        return np.full(mesh.t.shape[1], problem.kappa(""))

    n_elem = mesh.t.shape[1]
    kappa = np.full(n_elem, np.nan)
    material_names = [
        name for name in mesh.subdomains if name not in _RESERVED_SUBDOMAIN_NAMES
    ]
    for name in material_names:
        elem_ids = np.asarray(mesh.subdomains[name])
        kappa[elem_ids] = problem.kappa(name)

    if np.isnan(kappa).any():
        n_missing = int(np.isnan(kappa).sum())
        raise ValueError(
            f"{n_missing} of {n_elem} elements not covered by any Problem-known "
            f"subdomain; checked names={material_names}. Subdomain tagging is "
            f"by gmsh physical-surface name (ADR-0003); coordinate-based "
            f"assignment is forbidden."
        )
    return kappa


def _resolve_dirichlet_dofs(
    basis: Basis,
    bcs: dict[str, DirichletBC],
) -> tuple[np.ndarray, np.ndarray]:
    """Return (dirichlet_dof_indices, values_for_those_dofs).

    Each boundary name in ``bcs`` is looked up in ``mesh.boundaries`` and
    converted to a DOF-index array via ``basis.get_dofs(facets=...)``. Only
    ``DirichletBC(value: float)`` is supported in Part 1; callable Dirichlet
    and inhomogeneous Neumann are deferred (Protocol contract).
    """
    if not bcs:
        return np.array([], dtype=np.int64), np.array([], dtype=float)

    mesh = basis.mesh
    if not mesh.boundaries:
        raise ValueError(
            f"Problem declared Dirichlet BCs {list(bcs)} but mesh.boundaries "
            f"is empty - the geometry builder must tag boundary edges as "
            f"named gmsh physical curves."
        )

    dof_chunks: list[np.ndarray] = []
    value_chunks: list[np.ndarray] = []
    for name, spec in bcs.items():
        if not isinstance(spec, DirichletBC):
            raise NotImplementedError(
                f"BC spec {type(spec).__name__} for boundary {name!r} not "
                f"supported in Part 1; only DirichletBC is wired."
            )
        if name not in mesh.boundaries:
            raise KeyError(
                f"Dirichlet boundary {name!r} not in mesh.boundaries "
                f"(available: {list(mesh.boundaries)})"
            )
        dofs = basis.get_dofs(facets=mesh.boundaries[name]).flatten()
        dof_chunks.append(np.asarray(dofs, dtype=np.int64))
        value_chunks.append(np.full(dofs.shape, float(spec.value)))

    all_dofs = np.concatenate(dof_chunks)
    all_values = np.concatenate(value_chunks)
    # Deduplicate corner DOFs that belong to two named edges; silently keep
    # the first value seen. Moot for the only Part-1 case (homogeneous
    # Dirichlet, where every chunk's value is identical anyway). If an
    # inhomogeneous Dirichlet problem is added later, this needs to grow a
    # consistency check that raises when two edges disagree at a shared corner.
    unique_dofs, first_idx = np.unique(all_dofs, return_index=True)
    return unique_dofs, all_values[first_idx]


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

    Pipeline (architecture.md § Pipeline):
        geometry -> cached mesh -> P1 basis -> P0 kappa table -> assemble
        stiffness (w.kappa) and load -> Dirichlet DOFs from boundary tags +
        optional pin -> direct sparse solve.
    """
    cache_dir = mesh_cache_dir if mesh_cache_dir is not None else _DEFAULT_CACHE

    # --- Geometry / mesh ----------------------------------------------------
    geometry_builder = problem.geometry()
    spec = geometry_builder(mesh_size)
    msh_path = _materialise_geometry(spec, cache_dir)
    mesh = _load_mesh(msh_path)
    logger.debug(
        "loaded mesh: %d nodes, %d elements", mesh.p.shape[1], mesh.t.shape[1]
    )

    # --- Basis + per-element kappa ------------------------------------------
    basis = Basis(mesh, ElementTriP1())
    basis_p0 = basis.with_element(ElementTriP0())
    kappa_values = _build_kappa_p0_values(mesh, problem)
    kappa_field = basis_p0.interpolate(kappa_values)

    # --- Assembly -----------------------------------------------------------
    A = _stiffness_kappa.assemble(basis, kappa=kappa_field)
    b = _make_load(problem.source).assemble(basis)

    # --- BCs (Dirichlet) ----------------------------------------------------
    bcs = problem.boundary_conditions()
    dirichlet_dofs, dirichlet_values = _resolve_dirichlet_dofs(basis, bcs)
    has_dirichlet = dirichlet_dofs.size > 0

    # --- Nullspace pin (only if no Dirichlet) -------------------------------
    # architecture.md § Nullspace handling rule 1: pin iff no Dirichlet
    # exists. With Dirichlet present, the operator is non-singular on its own.
    x = np.zeros(basis.N)
    constrained_dofs: list[np.ndarray] = []
    pin_dof: int | None = None
    if dirichlet_dofs.size:
        x[dirichlet_dofs] = dirichlet_values
        constrained_dofs.append(dirichlet_dofs)

    if not has_dirichlet:
        pin_point = problem.pin_point()
        if pin_point is None:
            raise ValueError(
                "Pure-Neumann problem requires a pin_point() to remove the "
                "constant nullspace; problem returned None."
            )
        pin_dof = _nearest_node_dof(basis, pin_point)
        pin_value = float(
            problem.exact_solution(
                np.array([pin_point[0]]), np.array([pin_point[1]])
            )[0]
        )
        x[pin_dof] = pin_value
        constrained_dofs.append(np.array([pin_dof]))
        logger.debug(
            "pinning DOF %d to T=%.6f at %s", pin_dof, pin_value, pin_point
        )
    else:
        # Defensive sanity check: if pin_point() returns a coordinate while
        # Dirichlet is present, the Problem is internally inconsistent.
        if problem.pin_point() is not None:
            raise ValueError(
                "Problem declared Dirichlet BCs AND a pin_point(); rule 1 "
                "says pin iff no Dirichlet present. Pick one."
            )

    D = np.concatenate(constrained_dofs) if constrained_dofs else np.array([], dtype=np.int64)

    # --- Solve --------------------------------------------------------------
    solution = solve(*condense(A, b, x=x, D=D))

    return SolverResult(
        solution=solution,
        basis=basis,
        mesh=mesh,
        pin_dof=pin_dof,
    )
