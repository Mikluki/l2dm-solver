# ABOUTME: Verification Problem 1 (verification.md § Problem 1). Smooth
# manufactured solution T = cos(pi x) cos(pi y) on the unit square with pure
# zero-flux Neumann BCs and kappa = 1. Pin point is the (0, 0) corner where
# T = 1 - a guaranteed mesh node, well-conditioned, reproducible.

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from src.geometry.unit_square import build_unit_square
from src.problems.protocol import DirichletBC

# ============================================================================
# CONSTANTS
# ============================================================================

_PI = np.pi


# ============================================================================
# PROBLEM
# ============================================================================


@dataclass(frozen=True)
class Problem01Manufactured:
    """Smooth manufactured solution on the unit square, pure Neumann.

    See verification.md § Problem 1. Single subdomain (tag 1), zero-flux on
    every edge, kappa = 1. Source chosen so T(x, y) = cos(pi x) cos(pi y)
    solves the PDE exactly; mean(T) = 0 on the unit square, so the
    compatibility condition int Q dA = 0 is satisfied.
    """

    name: str = "problem_01_manufactured"

    def geometry(self) -> Callable[[float], object]:
        return build_unit_square

    def kappa(self, _subdomain_name: str) -> float:
        # Single subdomain, uniform conductivity. Argument unused; the name
        # parameter is kept for Protocol conformance (verification.md §
        # Problem definition interface).
        return 1.0

    def source(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return 2.0 * _PI * _PI * np.cos(_PI * x) * np.cos(_PI * y)

    def source_integral(self) -> float:
        """Analytic value of $\\int_\\Omega Q\\,dA$ for the Layer 2 source check.

        Closed form via separability + cosine symmetry on $[0,1]$:
        $\\int_0^1 \\cos(\\pi x)\\,dx = 0$, so the whole product integral vanishes.
        See artifacts/inspect/CONVENTIONS.md § Technique 2 — this must be a
        one-liner, not a re-derivation, or a matching FE bug hides forever.
        """
        return 0.0

    def boundary_conditions(self) -> dict[str, DirichletBC]:
        # Pure Neumann, zero flux everywhere. Empty dict signals "no Dirichlet
        # tags"; the zero-flux Neumann condition is the natural BC for the
        # bilinear form and requires no explicit contribution.
        return {}

    def exact_solution(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return np.cos(_PI * x) * np.cos(_PI * y)

    def exact_gradient(
        self, x: np.ndarray, y: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Analytic gradient of ``exact_solution``.

        Not part of the Protocol (verification.md § Problem definition
        interface fixes that contract); supplied here so the harness can
        compute an honest H^1-seminorm rate without projection pollution.
        Problems without a closed-form gradient can omit this and the harness
        falls back to a higher-order projection.
        """
        dT_dx = -_PI * np.sin(_PI * x) * np.cos(_PI * y)
        dT_dy = -_PI * np.cos(_PI * x) * np.sin(_PI * y)
        return dT_dx, dT_dy

    def mesh_sizes(self) -> list[float]:
        # Log-spaced sequence with seven levels. gmsh's unstructured triangle
        # quality varies enough between sizes that a four-level pure-halving
        # sequence produces noisy pairwise rates; the wider sequence smooths
        # the least-squares fit toward the asymptotic theoretical value.
        return [0.2, 0.14, 0.1, 0.07, 0.05, 0.035, 0.025]

    def expected_rate(self) -> float:
        return 2.0

    def pin_point(self) -> tuple[float, float] | None:
        # Geometry corner -> guaranteed mesh node at every refinement; exact
        # value T(0, 0) = 1, away from the saddle at (1/2, 1/2).
        return (0.0, 0.0)
