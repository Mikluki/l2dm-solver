# Submission 0001 — Problem 1 end-to-end vertical slice

**Status:** accepted
**Predecessors:** none (first submission)
**Successors:** 0003 (Problem 2 — piecewise-constant κ).

## Goal

Stand up the smallest vertical slice — geometry → mesh → solver → harness → pytest — that runs `verification.md` Problem 1 (smooth manufactured solution on a unit square, pure Neumann) and converges at the expected rates.

## Acceptance — outcome

1. **Test passes.** `uv run pytest tests/test_problem_01.py -v` ✓
2. **Convergence honest.** ≥ 3 mesh sizes; L² rate ≥ 1.8 and within 0.2 of 2.0; H¹ rate ≥ 0.9 and within 0.2 of 1.0. ✓
3. **Node-pinning convention applied as specified.** `Problem.pin_point()` returns the pin location; solver pins nearest mesh DOF to the exact solution; pin-DOF index reproducible across refinements. ✓
4. **Mean integral decays at FE-error scale.** Test asserts $|\int T_h\,dA| < 5\cdot\|T_h - T\|_{L^2}$ — orthogonal signal to the rate, catches sign-flipped assembly that still converges. ✓ Note: superseded by 0004 once the L² norm was corrected to use the analytic exact at quadrature points; the assertion shape is unchanged.
5. **Failure mode is visible.** Forced-failure (multiply assembled κ by 2): test fails on rate; `tests/_artifacts/test_problem_01/` populated with error table + rate + finest-mesh error plot; passing run leaves no artifact dir. ✓
6. **Mesh cache exercised.** Rerun produces no new `.msh` files. ✓
7. **Predicted first-run failure.** Initial run failed with `ModuleNotFoundError` on the not-yet-written solver — the expected diagnostic. ✓

## What shipped

```
src/problems/{__init__.py, protocol.py, problem_01_manufactured.py}
src/geometry/{__init__.py, cache.py, unit_square.py}
src/solver/{__init__.py, result.py, solve_scalar.py}
src/harness/{__init__.py, norms.py, study.py, artifacts.py}
tests/{conftest.py, test_problem_01.py}
.gitignore — adds tests/_artifacts/, tests/_mesh_cache/
```

`Problem` is a `typing.Protocol` (structural typing, no inheritance) with `pin_point()` declared by each Problem; solver has no default. Subtract-the-mean is not used anywhere — node-pinning is the project-wide convention, recorded in `architecture.md` § Nullspace handling.
