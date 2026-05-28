# ABOUTME: Structural typing.Protocol that every verification problem implements.
# Defines the cross-document contract from verification.md § Problem definition
# interface: geometry, kappa, source, boundary_conditions, exact_solution,
# mesh_sizes, expected_rate, and the nullspace pin_point accessor. Subdomains
# are keyed by string name (matching scikit-fem's mesh.subdomains dict); BCs
# are keyed by boundary string name and carry a DirichletBC dataclass payload.

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Protocol

import numpy as np

# ============================================================================
# BC SPEC
# ============================================================================


@dataclass(frozen=True)
class DirichletBC:
    """Homogeneous-or-constant Dirichlet boundary condition.

    Only the scalar ``value`` form is supported in Part 1; callable-valued
    Dirichlet (needed by Problem 5) and inhomogeneous Neumann are deferred
    until a problem actually exercises them. Boundary tags absent from
    ``Problem.boundary_conditions()`` are natural zero-flux Neumann.
    """

    value: float


# ============================================================================
# PROTOCOL
# ============================================================================


class Problem(Protocol):
    """Structural contract for a verification problem.

    Implementations are pure data + pure functions. They must not reference
    scikit-fem; the solver and harness consume them through this interface.
    """

    def geometry(self) -> Callable[[float], object]:
        """Return a callable ``build(mesh_size) -> gmsh-aware geometry handle``.

        The returned callable is what the geometry-cache layer keys on; the
        opaque handle is consumed by the solver's mesh-loading step.
        """
        ...

    def kappa(self, subdomain_name: str) -> float:
        """Conductivity for the named subdomain (W/(m K)).

        The name is the gmsh physical-group name propagated through meshio to
        ``mesh.subdomains``; subdomain assignment is never by element-centroid
        coordinate (ADR-0003, verification.md § Problem 2 failure diagnostic).
        """
        ...

    def source(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Vectorized source field Q(x, y)."""
        ...

    def boundary_conditions(self) -> dict[str, DirichletBC]:
        """Map boundary name -> BC specification.

        Unlisted boundary names are natural zero-flux Neumann (no
        contribution to the load vector). Only ``DirichletBC`` is currently
        supported; future inhomogeneous Neumann would extend this union.
        """
        ...

    def exact_solution(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Vectorized analytic solution T(x, y)."""
        ...

    def mesh_sizes(self) -> list[float]:
        """Target characteristic mesh sizes for the refinement study (>=3)."""
        ...

    def expected_rate(self) -> float:
        """Expected L^2 convergence rate for P1 elements."""
        ...

    def pin_point(self) -> tuple[float, float] | None:
        """Nullspace pin location.

        Returns ``(x, y)`` to pin the nearest-node DOF to ``exact_solution(x, y)``;
        returns ``None`` when at least one Dirichlet DOF is present and no pin is
        required. The solver has no default - the choice is explicit per ADR-0005
        and architecture.md § Nullspace handling.
        """
        ...
