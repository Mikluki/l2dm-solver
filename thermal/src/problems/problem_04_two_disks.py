# ABOUTME: Verification Problem 4 (verification.md § Problem 4). Two
# disconnected inner disks of radius R_inner centred at (+/- d_sep/2, 0)
# inside a larger annulus of radius R_outer. kappa_1=1 in both inner disks,
# kappa_2=10 in the outer annulus; uniform source q_0=1 in both inner disks,
# 0 outside. Dirichlet T=0 on r=R_outer. Exact solution is approximated by
# superposition of two shifted Problem-3 closed forms (submission 0008 §
# Decisions 1); the two stacked approximation errors — O(R_0/d_sep)
# finite-separation and O((d_sep/(2 R_out))^2) joint-boundary residual — are
# intentional and bounded by the 10% acceptance. The Protocol surface is
# implemented unchanged: expected_rate() returns NaN and mesh_sizes() a
# single-element list because Problem 4's convergence variable is the
# geometric ratio R_0/d_sep, not h (the sweep is driven by the test, not by
# run_refinement_study; submission 0008 § Decisions 2).

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from src.geometry.two_disks_in_disk import (
    INNER_DISK_A_NAME,
    INNER_DISK_B_NAME,
    OUTER_ANNULUS_NAME,
    OUTER_BOUNDARY_NAME,
    build_two_disks_in_disk,
)
from src.problems.protocol import DirichletBC

# ============================================================================
# CONSTANTS
# ============================================================================

# Materials per verification.md § Problem 4 (also matches Problem 3 so the
# superposed reference is literally the Problem 3 closed form at each centre).
_KAPPA_1: float = 1.0
_KAPPA_2: float = 10.0
_Q_0: float = 1.0


# ============================================================================
# HELPERS
# ============================================================================


def _t_single(
    r: np.ndarray, R_inner: float, R_outer: float
) -> np.ndarray:
    """Problem-3 closed form T_single(r; R_inner, kappa_1=1, kappa_2=10, R_outer).

    Returns the radial Problem 3 solution at distance r from a disk centre.
    Used to build the superposed reference; carries the same kappa_1/kappa_2
    constants Problem 4 fixes. Vectorized in r.

    The two stacked approximation errors (finite-separation + common-boundary
    residual) live in the *superposition*, not here — this routine is exact
    for the corresponding single-disk problem.
    """
    r = np.asarray(r, dtype=float)
    inside = r < R_inner

    t_in = (
        _Q_0 * (R_inner**2 - r**2) / (4.0 * _KAPPA_1)
        + _Q_0 * R_inner**2 / (2.0 * _KAPPA_2) * np.log(R_outer / R_inner)
    )
    # Guard r=0 (inside region) from log(R_out/0). np.where evaluates both
    # branches; substitute a safe r for inside points so log stays finite.
    r_safe = np.where(inside, 1.0, r)
    t_out = _Q_0 * R_inner**2 / (2.0 * _KAPPA_2) * np.log(R_outer / r_safe)

    return np.where(inside, t_in, t_out)


# ============================================================================
# PROBLEM
# ============================================================================


@dataclass(frozen=True)
class Problem04TwoDisks:
    """Two well-separated disks in a common annulus.

    Geometry is parametrised by R_inner, d_sep, R_outer; the inner-disk
    centres sit at (+/- d_sep/2, 0). The test exercises a sweep of geometric
    configurations at a fixed mesh size (submission 0008 § Decisions 2),
    not a refinement study; expected_rate() therefore returns NaN and
    mesh_sizes() returns a single-element list.
    """

    R_inner: float
    d_sep: float
    R_outer: float
    mesh_size: float
    kappa_1: float = _KAPPA_1
    kappa_2: float = _KAPPA_2
    # Overrides for the disk-B kappa and source value. When None the disk
    # inherits the symmetric defaults (kappa_1 and q_0) so T_h is mirror-
    # symmetric. The forced-failure check in test_problem_04 sets q_inner_B=0
    # to break mirror symmetry to O(1) at annulus probes (kappa_inner_B alone
    # cannot do this because Problem 3's outer solution is kappa_1-independent
    # — submission 0008 § Decisions 5 anticipated kappa as the mechanism, but
    # the kappa_1-independence theorem from Problem 3 makes it bite only at
    # second order in the two-disk coupling; q is the load-bearing knob).
    kappa_inner_B: float | None = None
    q_inner_B: float | None = None
    name: str = field(default="problem_04_two_disks")

    # --- Geometry handles ---------------------------------------------------

    @property
    def center_a(self) -> tuple[float, float]:
        return (-0.5 * self.d_sep, 0.0)

    @property
    def center_b(self) -> tuple[float, float]:
        return (+0.5 * self.d_sep, 0.0)

    # --- Protocol surface ---------------------------------------------------

    def geometry(self) -> Callable[[float], object]:
        return build_two_disks_in_disk(self.R_inner, self.d_sep, self.R_outer)

    def kappa(self, subdomain_name: str) -> float:
        # Subdomain assignment is by gmsh physical-surface name (ADR-0003).
        # Distinct A/B tags both map to kappa_1 by default; the forced-failure
        # check in the test sets kappa_inner_B asymmetrically to break mirror
        # symmetry of T_h (submission 0008 § Decisions 6).
        if subdomain_name == INNER_DISK_A_NAME:
            return self.kappa_1
        if subdomain_name == INNER_DISK_B_NAME:
            return (
                self.kappa_1 if self.kappa_inner_B is None else self.kappa_inner_B
            )
        if subdomain_name == OUTER_ANNULUS_NAME:
            return self.kappa_2
        raise KeyError(
            f"unknown subdomain {subdomain_name!r}; expected one of "
            f"{INNER_DISK_A_NAME!r}, {INNER_DISK_B_NAME!r}, "
            f"{OUTER_ANNULUS_NAME!r}"
        )

    def source(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        # q_0 inside either inner disk, 0 in the outer annulus. Mesh alignment
        # with both r=R_inner circles (via OCC fragment) ensures each element's
        # quadrature points are wholly on one side, so the discontinuity is
        # never straddled. q_inner_B overrides the source value inside disk B
        # only; in the symmetric case it equals q_0 (None default) and the
        # field is even in x.
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        xa, ya = self.center_a
        xb, yb = self.center_b
        in_a = (x - xa) ** 2 + (y - ya) ** 2 < self.R_inner**2
        in_b = (x - xb) ** 2 + (y - yb) ** 2 < self.R_inner**2
        q_b = _Q_0 if self.q_inner_B is None else self.q_inner_B
        out = np.zeros_like(x)
        out = np.where(in_a, _Q_0, out)
        out = np.where(in_b, q_b, out)
        return out

    def boundary_conditions(self) -> dict[str, DirichletBC]:
        # Dirichlet T=0 on the outer circle kills the nullspace. Both inner
        # circles are interior material interfaces — no BC (submission 0008
        # § Decisions 8).
        return {OUTER_BOUNDARY_NAME: DirichletBC(value=0.0)}

    def exact_solution(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Superposition reference per submission 0008 § Decisions 1.

        T_ref(r) = T_single(|r - r_A|) + T_single(|r - r_B|), each evaluated
        with the *joint* outer radius R_outer. Approximation, not exact PDE
        solution; the 10% acceptance accommodates the stacked errors.
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        xa, ya = self.center_a
        xb, yb = self.center_b
        ra = np.sqrt((x - xa) ** 2 + (y - ya) ** 2)
        rb = np.sqrt((x - xb) ** 2 + (y - yb) ** 2)
        return _t_single(ra, self.R_inner, self.R_outer) + _t_single(
            rb, self.R_inner, self.R_outer
        )

    def mesh_sizes(self) -> list[float]:
        # Single fixed h per submission 0008 § Decisions 2 — the convergence
        # variable here is geometric (R_inner/d_sep), not h.
        return [self.mesh_size]

    def expected_rate(self) -> float:
        # No h-refinement rate is defined for Problem 4 (submission 0008 §
        # Decisions 2); NaN signals "not applicable".
        return float("nan")

    def pin_point(self) -> tuple[float, float] | None:
        # Dirichlet on the outer circle kills the nullspace; no pin per
        # architecture.md § Nullspace handling rule 1.
        return None


__all__ = ["Problem04TwoDisks"]
