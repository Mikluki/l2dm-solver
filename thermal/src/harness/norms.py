# ABOUTME: L^2 and H^1-seminorm error functionals built on scikit-fem's
# Functional, so the integration quadrature matches assembly. Each takes a
# SolverResult and the problem's vectorized exact_solution callable.

from __future__ import annotations

from typing import Callable

import numpy as np
from skfem import Functional
from skfem.helpers import dot, grad

# ============================================================================
# FORMS
# ============================================================================


@Functional
def _l2_squared(w):
    # (T_h - T_exact)^2 integrated over the domain.
    return (w["uh"] - w["uex"]) ** 2


@Functional
def _h1_semi_squared(w):
    # |grad(T_h) - grad(T_exact)|^2 integrated. H^1 *seminorm* per
    # verification.md § Scope - it is what the doc asserts and is the right
    # rate measurement for P1 elements.
    diff = grad(w["uh"]) - w["gex"]
    return dot(diff, diff)


# ============================================================================
# FUNCTIONS
# ============================================================================


def l2_error(
    solution: np.ndarray,
    basis,
    exact: Callable[[np.ndarray, np.ndarray], np.ndarray],
) -> float:
    """L^2 error between the discrete solution and ``exact``.

    The analytic exact solution is evaluated directly at the basis quadrature
    points. This is an FE quadrature norm, matching the assembly quadrature,
    without first projecting ``exact`` into the P1 space.
    """
    qpts = basis.global_coordinates()
    qpts_arr = np.asarray(qpts)
    uh = basis.interpolate(solution)
    uex = exact(qpts_arr[0], qpts_arr[1])
    val = _l2_squared.assemble(basis, uh=uh, uex=uex)
    return float(np.sqrt(val))


def h1_seminorm_error(
    solution: np.ndarray,
    basis,
    exact_gradient: Callable[[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]],
) -> float:
    """H^1 seminorm error using an analytic gradient of the exact solution.

    Using the analytic gradient (rather than a nodal projection followed by
    discrete grad) keeps the rate measurement honest at P1 - the projection's
    gradient is one order lower and would dominate the residual.
    """
    qpts = basis.global_coordinates()  # (2, n_elem, n_qp) DiscreteField-like
    qpts_arr = np.asarray(qpts)
    gx_arr, gy_arr = exact_gradient(qpts_arr[0], qpts_arr[1])
    gex = np.stack([gx_arr, gy_arr], axis=0)
    uh = basis.interpolate(solution)
    val = _h1_semi_squared.assemble(basis, uh=uh, gex=gex)
    return float(np.sqrt(val))
