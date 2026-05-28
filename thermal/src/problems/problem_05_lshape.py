# ABOUTME: Verification Problem 5 (verification.md § Problem 5). L-shape
# [0,1]^2 \ [1/2,1]^2 with a reentrant corner at (1/2, 1/2), uniform kappa=1,
# zero source, singular exact solution T(r, theta) = r^{2/3} sin(2 theta/3) in
# polar coords centred on the corner. The four outer edges (south, east_lower,
# north_left, west) carry the singular solution as a callable Dirichlet BC -
# this is the protocol extension named in submission 0007 § Decisions 1. The
# two cut edges meeting at the corner are scalar T=0 (the singular solution
# vanishes identically there), keeping the callable path exercised only by
# nonzero edges. Inverted L^2 rate window [1.2, 1.5] confirms the rate fitter
# honestly reports the corner-degraded 4/3 convergence.

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from src.geometry.l_shape import (
    CUT_EAST_NAME,
    CUT_NORTH_NAME,
    EAST_LOWER_NAME,
    INTERIOR_NAME,
    NORTH_LEFT_NAME,
    SOUTH_NAME,
    WEST_NAME,
    build_l_shape,
)
from src.problems.protocol import DirichletBC

# ============================================================================
# CONSTANTS
# ============================================================================

# Reentrant corner at (1/2, 1/2); polar coords are centred here per
# submission 0007 § Decisions 3.
_CORNER_X: float = 0.5
_CORNER_Y: float = 0.5

# Exponent of the singular solution. The L-shape interior subtends 3pi/2,
# so the leading singular mode is r^{pi/(3pi/2)} = r^{2/3}.
_ALPHA: float = 2.0 / 3.0


# ============================================================================
# HELPERS
# ============================================================================


def _polar(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (r, theta_cw) centred on the reentrant corner.

    theta_cw is the clockwise angle from the (1/2,1/2) -> (1, 1/2) edge,
    valued in [0, 2*pi). Submission 0007 § Decisions 3 specifies a piecewise
    closed form; the modular form ``(-atan2) mod 2*pi`` is mathematically
    equivalent everywhere except at the negative-x-axis edge case
    (xhat<0, yhat=0) where the piecewise form returns -pi instead of +pi.
    That edge case happens to coincide with the midpoint DOF (0, 1/2) of the
    west boundary, so the discrepancy would mis-set one Dirichlet value and
    pollute the solution. The modular form returns the geometric limit at
    that point and matches the piecewise form everywhere else.
    """
    xhat = x - _CORNER_X
    yhat = y - _CORNER_Y
    r = np.sqrt(xhat**2 + yhat**2)
    theta = np.mod(-np.arctan2(yhat, xhat), 2.0 * np.pi)
    return r, theta


# ============================================================================
# PROBLEM
# ============================================================================


@dataclass(frozen=True)
class Problem05LShape:
    """L-shape reentrant-corner verification (verification.md § Problem 5).

    Single subdomain, kappa=1, Q=0, Dirichlet on every edge. The L^2 rate
    is asserted to lie in [1.2, 1.5] - both bounds load-bearing (submission
    0007 § Decisions 2). H^1 rate is logged but not asserted (~2/3 expected).
    """

    name: str = "problem_05_lshape"

    # --- Protocol surface ---------------------------------------------------

    def geometry(self) -> Callable[[float], object]:
        return build_l_shape

    def kappa(self, subdomain_name: str) -> float:
        # Single subdomain; kappa = 1 uniformly.
        if subdomain_name == INTERIOR_NAME:
            return 1.0
        raise KeyError(
            f"unknown subdomain {subdomain_name!r}; expected {INTERIOR_NAME!r}"
        )

    def source(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return np.zeros_like(x, dtype=float)

    def boundary_conditions(self) -> dict[str, DirichletBC]:
        # Two cut edges: T identically zero -> scalar form (keeps the callable
        # code path exercised only by the four nonzero edges, per submission
        # 0007 § Decisions 4).
        # Four outer edges: callable form, evaluates the singular exact
        # solution at each boundary DOF.
        return {
            CUT_EAST_NAME: DirichletBC(value=0.0),
            CUT_NORTH_NAME: DirichletBC(value=0.0),
            SOUTH_NAME: DirichletBC(value=self.exact_solution),
            EAST_LOWER_NAME: DirichletBC(value=self.exact_solution),
            NORTH_LEFT_NAME: DirichletBC(value=self.exact_solution),
            WEST_NAME: DirichletBC(value=self.exact_solution),
        }

    def exact_solution(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        r, theta = _polar(x, y)
        return r**_ALPHA * np.sin(_ALPHA * theta)

    def exact_gradient(
        self, x: np.ndarray, y: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Analytic gradient of ``exact_solution``.

        Diverges like r^{-1/3} at the corner; np.where guards r=0 to keep
        evaluation finite (submission 0007 § Decisions-left-to-worker:
        no quadrature scheme places a node *at* the corner for any mesh
        in the Decision-6 range, so the guarded value is never load-bearing).

        Derivation (sketch). With xh = x-1/2, yh = y-1/2, r = sqrt(xh^2+yh^2),
        theta = (-atan2(yh, xh)) mod 2*pi:
            dr/dx = xh/r,   dr/dy = yh/r
            dtheta/dx = yh/r^2,  dtheta/dy = -xh/r^2
        Then alpha = 2/3 gives
            dT/dx = alpha * r^{-4/3} * (xh * sin(alpha*theta) + yh * cos(alpha*theta))
            dT/dy = alpha * r^{-4/3} * (yh * sin(alpha*theta) - xh * cos(alpha*theta))
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        xhat = x - _CORNER_X
        yhat = y - _CORNER_Y
        r = np.sqrt(xhat**2 + yhat**2)
        # r=0 only at the reentrant corner; the guard keeps r_safe**(-4/3)
        # finite at that single point (gradient is set to zero there by the
        # np.where below - the actual exact gradient diverges, but no
        # quadrature point lands on the corner).
        r_safe = np.where(r > 0, r, 1.0)
        _, theta = _polar(x, y)
        s = np.sin(_ALPHA * theta)
        c = np.cos(_ALPHA * theta)
        prefactor = _ALPHA * r_safe ** (-4.0 / 3.0)
        gx = np.where(r > 0, prefactor * (xhat * s + yhat * c), 0.0)
        gy = np.where(r > 0, prefactor * (yhat * s - xhat * c), 0.0)
        return gx, gy

    def mesh_sizes(self) -> list[float]:
        # Five levels spanning ~9x; submission 0007 § Decisions 6 gives the
        # range. Slow 4/3 decay needs more headroom than rate-2 problems for
        # the least-squares fit to land in [1.2, 1.5].
        return [0.12, 0.07, 0.040, 0.022, 0.013]

    def expected_rate(self) -> float:
        return 4.0 / 3.0

    def pin_point(self) -> tuple[float, float] | None:
        # Every edge is Dirichlet -> the operator is non-singular on its own
        # (architecture.md § Nullspace handling rule 1).
        return None


__all__ = ["Problem05LShape"]
