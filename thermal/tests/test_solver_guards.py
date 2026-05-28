# ABOUTME: Direct unit tests for solver-internal guards that aren't exercised
# end-to-end by any verification problem. Covers solve_scalar's pin/Dirichlet
# exclusivity rule (architecture.md § Nullspace handling rule 1): declaring
# both Dirichlet BCs and a pin_point() over-constrains the nullspace and can
# silently mask wrong-edge BC bugs by clamping the discrete solution to the
# correct value at the pinned node.

from __future__ import annotations

import pytest

from src.problems.problem_02_slab import Problem02Slab
from src.solver.solve_scalar import solve_scalar


class _DirichletPlusPin:
    """Wrap a Dirichlet-equipped Problem and force pin_point() to a coordinate.

    __getattr__ delegates every other Protocol method to the wrapped problem,
    so the composite still satisfies the Problem protocol. Only pin_point is
    overridden — that's the exact misconfiguration the solver guard rejects.
    """

    def __init__(self, inner, pin_point):
        self._inner = inner
        self._pin_point = pin_point

    def __getattr__(self, name):
        return getattr(self._inner, name)

    def pin_point(self):
        return self._pin_point


def test_solver_refuses_dirichlet_and_pin(mesh_cache_dir):
    # Problem 2 has Dirichlet on the left edge; layering a pin on top is the
    # internal-inconsistency the guard exists to catch.
    problem = _DirichletPlusPin(Problem02Slab(), pin_point=(0.5, 0.05))

    # Coarsest mesh from Problem 2's refinement list — cheapest mesh that the
    # rectangle_split geometry builder can produce; the assembly happens but
    # the guard fires before any solve.
    mesh_size = Problem02Slab().mesh_sizes()[0]

    with pytest.raises(ValueError, match=r"Dirichlet BCs AND a pin_point"):
        solve_scalar(problem, mesh_size=mesh_size, mesh_cache_dir=mesh_cache_dir)
