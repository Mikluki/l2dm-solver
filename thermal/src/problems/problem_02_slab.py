# ABOUTME: Verification Problem 2 (verification.md § Problem 2). Piecewise-
# constant kappa, 1D slab on a rectangle split at x=1/2. Left subdomain has
# kappa_1=1 and uniform source q_0=1; right subdomain has variable kappa_2 and
# zero source. Dirichlet T=0 on the left edge, zero-flux Neumann everywhere
# else; no nullspace pin (Dirichlet present, architecture.md § Nullspace
# handling rule 1). The exact solution is parabolic in the left region and
# constant T = q_0/(8 kappa_1) = 0.125 in the right region - independent of
# kappa_2 (algebraic confirmation in submission 0002 § Problem 2).

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from src.geometry.rectangle_split import (
    LEFT_SUBDOMAIN_NAME,
    RIGHT_SUBDOMAIN_NAME,
    SLAB_HEIGHT,
    build_rectangle_split,
)
from src.problems.protocol import DirichletBC

# ============================================================================
# CONSTANTS
# ============================================================================

# Defaults tied to verification.md § Problem 2. kappa_2 is a constructor arg
# (varied in the kappa_2-sweep acceptance check); the other values are fixed.
_KAPPA_1 = 1.0
_KAPPA_2_DEFAULT = 100.0
_Q_0 = 1.0
_X_INTERFACE = 0.5
_LEFT_EDGE_NAME = "left"

# Exact right-region constant T(x>1/2) = q_0 / (8 kappa_1). Pre-computed for
# the kappa_2-sweep acceptance test as the bound-vs-exact reference.
RIGHT_REGION_EXACT_T = _Q_0 / (8.0 * _KAPPA_1)


# ============================================================================
# PROBLEM
# ============================================================================


@dataclass(frozen=True)
class Problem02Slab:
    """Piecewise-constant kappa on a 1D slab; varies kappa_2 by construction.

    The construct holds kappa_2 explicitly so the kappa_2-sweep test can
    instantiate variants without touching geometry; mesh-cache hits across
    the sweep are the consequence (kappa is not a geometry parameter, so
    ADR-0007 keys do not depend on it).
    """

    kappa_2: float = _KAPPA_2_DEFAULT
    name: str = field(default="problem_02_slab")

    # --- Protocol surface ---------------------------------------------------

    def geometry(self) -> Callable[[float], object]:
        return build_rectangle_split

    def kappa(self, subdomain_name: str) -> float:
        # Subdomain assignment is by gmsh physical-surface name (the planner
        # decision in submission 0003 § Decisions resolved) - never by
        # element-centroid coordinate (ADR-0003).
        if subdomain_name == LEFT_SUBDOMAIN_NAME:
            return _KAPPA_1
        if subdomain_name == RIGHT_SUBDOMAIN_NAME:
            return self.kappa_2
        raise KeyError(
            f"unknown subdomain {subdomain_name!r}; "
            f"expected {LEFT_SUBDOMAIN_NAME!r} or {RIGHT_SUBDOMAIN_NAME!r}"
        )

    def source(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        # q_0 in the left subdomain, 0 in the right. The mesh aligns with
        # x=1/2 (ADR-0003) so quadrature points never straddle the
        # discontinuity - each element is wholly in one subdomain.
        return np.where(x < _X_INTERFACE, _Q_0, 0.0)

    def boundary_conditions(self) -> dict[str, DirichletBC]:
        # Dirichlet T=0 on the left edge; the other three sides are natural
        # zero-flux Neumann (unlisted = Neumann zero, per Protocol contract).
        return {_LEFT_EDGE_NAME: DirichletBC(value=0.0)}

    def exact_solution(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        # Piecewise: parabolic for x<1/2, constant for x>=1/2. The constant
        # is q_0/(8 kappa_1) and explicitly does not depend on kappa_2 - the
        # load-bearing structural property the kappa_2-sweep asserts (proof
        # in submission 0002 § Problem 2).
        x = np.asarray(x, dtype=float)
        return np.where(
            x < _X_INTERFACE,
            _Q_0 * x * (1.0 - x) / (2.0 * _KAPPA_1),
            RIGHT_REGION_EXACT_T,
        )

    def exact_gradient(
        self, x: np.ndarray, y: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Analytic gradient of ``exact_solution``.

        Used by the harness's H^1-seminorm rate diagnostic. The right region
        is constant, so dT/dx vanishes there; the left region has
        dT/dx = q_0 (1 - 2x) / (2 kappa_1). dT/dy is identically zero (the
        framing is 1D in x). Per-quadrature-point evaluation is exact off
        the interface because the mesh is aligned with x=1/2; values on the
        measure-zero interface itself are arbitrary.
        """
        x = np.asarray(x, dtype=float)
        dT_dx = np.where(
            x < _X_INTERFACE,
            _Q_0 * (1.0 - 2.0 * x) / (2.0 * _KAPPA_1),
            0.0,
        )
        dT_dy = np.zeros_like(dT_dx)
        return dT_dx, dT_dy

    def mesh_sizes(self) -> list[float]:
        # Five log-spaced sizes. Lower bound chosen large enough that the
        # slab height h=0.1 (SLAB_HEIGHT) is resolved (>=2 elements through
        # the thickness); upper bound chosen small enough to give an honest
        # asymptotic fit. The bilinear-form rate is independent of the
        # framing direction, so this list need not be as long as Problem 1.
        return [0.05, 0.035, 0.025, 0.017, 0.012]

    def expected_rate(self) -> float:
        return 2.0

    def pin_point(self) -> tuple[float, float] | None:
        # Dirichlet on the left edge kills the nullspace; no pin per
        # architecture.md § Nullspace handling rule 1.
        return None


# Re-exported for tests that want to probe the right-region constant without
# importing the module-private symbol.
__all__ = ["Problem02Slab", "RIGHT_REGION_EXACT_T", "SLAB_HEIGHT"]
