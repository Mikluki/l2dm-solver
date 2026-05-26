# ABOUTME: Structural typing.Protocol that every verification problem implements.
# Defines the cross-document contract from verification.md § Problem definition
# interface: geometry, kappa, source, boundary_conditions, exact_solution,
# mesh_sizes, expected_rate, and the nullspace pin_point accessor.

from __future__ import annotations

from typing import Callable, Protocol

import numpy as np

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

    def kappa(self, subdomain_tag: int) -> float:
        """Conductivity for the given subdomain tag (W/(m K))."""
        ...

    def source(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Vectorized source field Q(x, y)."""
        ...

    def boundary_conditions(self) -> dict[int, object]:
        """Map boundary tag -> BC specification (Dirichlet value or Neumann flux)."""
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
