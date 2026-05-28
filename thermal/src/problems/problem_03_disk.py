# ABOUTME: Verification Problem 3 (verification.md § Problem 3). Radially
# symmetric disk in disk: inner disk radius R_0=0.3 with kappa_1, uniform
# source q_0=1; outer annulus radius R_out=1.0 with kappa_2=10, zero source.
# Dirichlet T=0 on the outer circle; no nullspace pin (Dirichlet present,
# architecture.md § Nullspace handling rule 1). kappa_1 is a constructor
# parameter so the kappa_1-independence acceptance sweep can vary it without
# changing geometry (mesh-cache hits across the sweep are the consequence;
# kappa does not enter the cache key per brief § Decisions 5). The exact outer
# solution T_out(r) depends only on kappa_2, confirming the structural property
# that the outer region is kappa_1-independent (algebraic proof in
# docs/derivations/algebraic-verification.md § Problem 3).

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from src.geometry.disk_in_disk import (
    INNER_DISK_NAME,
    OUTER_ANNULUS_NAME,
    OUTER_BOUNDARY_NAME,
    build_disk_in_disk,
)
from src.problems.protocol import DirichletBC

# ============================================================================
# CONSTANTS
# ============================================================================

# Geometry and physics per verification.md § Problem 3. kappa_1 is a
# constructor arg (varied in the kappa_1-sweep acceptance check). kappa_2 and
# geometry parameters are fixed per brief § Decisions 4.
_R_INNER: float = 0.3
_R_OUTER: float = 1.0
_KAPPA_1_DEFAULT: float = 1.0
_KAPPA_2: float = 10.0
_Q_0: float = 1.0

# T_out at r=0.6 (outer annulus), independent of kappa_1. Pre-computed for the
# kappa_1-sweep bound check in the acceptance test. Derivation in
# algebraic-verification.md § Problem 3: T_out(r) = q_0 R_0^2/(2 kappa_2) ln(R_out/r).
OUTER_PROBE_R: float = 0.6
OUTER_PROBE_EXACT_T: float = (
    _Q_0 * _R_INNER**2 / (2.0 * _KAPPA_2) * np.log(_R_OUTER / OUTER_PROBE_R)
)


# ============================================================================
# PROBLEM
# ============================================================================


@dataclass(frozen=True)
class Problem03Disk:
    """Radially symmetric disk-in-disk; kappa_1 is varied in the sweep test.

    kappa_2 and the outer radius are fixed per verification.md. R_inner is
    constructor-parameterised so Problem 4's per-config FE-error proxy
    (submission 0008 § Decisions 2) can evaluate the same closed form at
    smaller inner radii without re-tagging or duplicating the problem class.
    The dataclass holds kappa_1 explicitly so sweep instantiation never
    touches geometry; cache hits are the consequence (kappa is not a geometry
    parameter).
    """

    kappa_1: float = _KAPPA_1_DEFAULT
    R_inner: float = _R_INNER
    name: str = field(default="problem_03_disk")

    # --- Protocol surface ---------------------------------------------------

    def geometry(self) -> Callable[[float], object]:
        return build_disk_in_disk(self.R_inner, _R_OUTER)

    def kappa(self, subdomain_name: str) -> float:
        # Subdomain assignment is by gmsh physical-surface name (brief § Decisions 6).
        # Never by element-centroid coordinate (ADR-0003).
        if subdomain_name == INNER_DISK_NAME:
            return self.kappa_1
        if subdomain_name == OUTER_ANNULUS_NAME:
            return _KAPPA_2
        raise KeyError(
            f"unknown subdomain {subdomain_name!r}; "
            f"expected {INNER_DISK_NAME!r} or {OUTER_ANNULUS_NAME!r}"
        )

    def source(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        # q_0 inside the inner disk, 0 in the outer annulus. Mesh alignment
        # with r=R_inner (via OCC fragment) ensures each element's quadrature
        # points are wholly on one side, so the discontinuity is never straddled.
        r = np.sqrt(x**2 + y**2)
        return np.where(r < self.R_inner, _Q_0, 0.0)

    def boundary_conditions(self) -> dict[str, DirichletBC]:
        # Dirichlet T=0 on the outer circle; outer circle kills the nullspace.
        # Inner circle is a material interface — no BC (brief § Decisions 7–8).
        return {OUTER_BOUNDARY_NAME: DirichletBC(value=0.0)}

    def exact_solution(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        # Piecewise radial solution from verification.md § Problem 3.
        # Inner: T_in(r) = q_0(R_0^2 - r^2)/(4 kappa_1)
        #                 + q_0 R_0^2/(2 kappa_2) ln(R_out/R_0)
        # Outer: T_out(r) = q_0 R_0^2/(2 kappa_2) ln(R_out/r)
        # T_out depends on kappa_2 only — the load-bearing structural property
        # the kappa_1-sweep asserts (proof in algebraic-verification.md § P3).
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        r = np.sqrt(x**2 + y**2)
        inside = r < self.R_inner

        t_in = (
            _Q_0 * (self.R_inner**2 - r**2) / (4.0 * self.kappa_1)
            + _Q_0 * self.R_inner**2 / (2.0 * _KAPPA_2) * np.log(_R_OUTER / self.R_inner)
        )
        # Guard r=0 (in inside region) from log(R_out / 0). np.where evaluates
        # both branches; substitute r=1 for inside points so log stays finite.
        r_safe = np.where(inside, 1.0, r)
        t_out = _Q_0 * self.R_inner**2 / (2.0 * _KAPPA_2) * np.log(_R_OUTER / r_safe)

        return np.where(inside, t_in, t_out)

    def exact_gradient(
        self, x: np.ndarray, y: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Analytic gradient of exact_solution.

        Inner: dT_in/dx = -q_0 x / (2 kappa_1), dT_in/dy = -q_0 y / (2 kappa_1).
        Outer: dT_out/dx = -q_0 R_0^2 x / (2 kappa_2 r^2).
        At r=0 (inside region): gradient = 0 by L'Hopital. np.where guard used
        to avoid the r^2 denominator blowing up on the inside branch evaluation.
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        r_sq = x**2 + y**2
        inside = r_sq < self.R_inner**2

        # Inside gradient: well-defined everywhere (= 0 at origin by symmetry)
        gx_in = -_Q_0 * x / (2.0 * self.kappa_1)
        gy_in = -_Q_0 * y / (2.0 * self.kappa_1)

        # Outside gradient: guard r^2=0 for inside points where r<R_inner
        r_sq_safe = np.where(inside, 1.0, r_sq)
        factor = -_Q_0 * self.R_inner**2 / (2.0 * _KAPPA_2 * r_sq_safe)
        gx_out = factor * x
        gy_out = factor * y

        return np.where(inside, gx_in, gx_out), np.where(inside, gy_in, gy_out)

    def mesh_sizes(self) -> list[float]:
        # Five sizes spanning a factor of ~7. Lower bound: at h=0.015 the inner
        # circle (circumference ~1.88) has ~125 nodes — well-resolved. Upper bound
        # h=0.10 gives ~19 nodes around inner circle — adequate for honest rate fit.
        # Finest bound is comparable to Problem 2's finest h≈0.012.
        return [0.10, 0.065, 0.040, 0.025, 0.015]

    def expected_rate(self) -> float:
        return 2.0

    def pin_point(self) -> tuple[float, float] | None:
        # Dirichlet on the outer circle kills the nullspace; no pin per
        # architecture.md § Nullspace handling rule 1.
        return None


__all__ = ["Problem03Disk", "OUTER_PROBE_EXACT_T", "OUTER_PROBE_R"]
