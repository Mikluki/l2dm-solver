# Submission 0007 — Problem 5 (L-shape reentrant corner)

**Status:** accepted
**Predecessors:** 0006 (Problem 3 — OCC pattern, Dirichlet-DOF lookup); 0004 (L² norm at quadrature points — the rate-fitter this brief inverts an assertion against).
**Successors:** any future verification problem with inhomogeneous Dirichlet; Part 2 problems carrying non-constant BCs.

## Goal

Implement `verification.md` § Problem 5 and the tests that exercise it. The submission is a **harness-correctness** check, not a solver one: the inverted-assertion rate window $1.2 \le r_{L^2} \le 1.5$ confirms the rate-fitter honestly reports the corner-degraded $4/3$ convergence rather than masking it. The L-shape geometry is a single subdomain with $\kappa=1$, $Q=0$; the solver-side novelty is **callable Dirichlet** — the protocol extension the `DirichletBC` docstring already named Problem 5 as the trigger for.

## Acceptance — outcome

1. **Tests pass — Problem 5 and Problems 1–3.** Existing tests not regressed. The `DirichletBC.value` extension does not break the scalar-value path. ✓
2. **Inverted-rate window holds:** $1.2 \le r_{L^2} \le 1.5$ over 5 mesh refinements. Both bounds load-bearing — upper-bound failure means the rate-fitter is dishonestly reporting smooth convergence; lower-bound failure means the fitter is broken or the mesh is outside the asymptotic regime. Artifact bundle records which bound blew up. ✓
3. **H¹ rate logged, not asserted.** ✓ (rate ~2/3, as expected from $\nabla T \sim r^{-1/3}$).
4. **Forced-failure on the inversion mechanism.** A synthetic `StudyResult` with $r_{L^2} \approx 2.0$ run through the rate-window check raises. The main test and the forced-failure test call the *same* extracted helper, so copy-drift cannot defeat the inversion logic. ✓
5. **Callable-Dirichlet wiring is loud-fail.** Replacing the four nonzero-edge callables with `DirichletBC(value=0.0)` made the unique solution $T \equiv 0$; the L² discrepancy stopped decaying; rate dropped to ≈0; bundle emitted. ✓
6. **Reentrant corner is a guaranteed mesh node** at every refinement: $\|p_{\text{nearest}} - (1/2, 1/2)\|_\infty < 10^{-12}$. ✓
7. **Mesh cache hits across re-runs.** ✓

## What shipped

```
src/geometry/l_shape.py
src/problems/problem_05_lshape.py
src/problems/protocol.py — DirichletBC.value extended to float | Callable[[x, y], np.ndarray]
src/solver/solve_scalar.py — _resolve_dirichlet_dofs type-dispatches on the value field
tests/test_problem_05.py
```

Key conventions baked in: polar angle uses the modular form `θ = (-atan2(ŷ, x̂)) mod 2π` (not piecewise on $\hat y$ — the piecewise form returns $-\pi$ at the west-edge midpoint and mis-sets one Dirichlet value); the two cut edges meeting at the reentrant corner use scalar `DirichletBC(value=0.0)` while the four perimeter edges use the callable form, so the callable code path is exercised only where the boundary value is non-zero; mesh-size list spans ~8× (`[0.12, 0.07, 0.04, 0.022, 0.013]`) to give the slow $4/3$ rate enough signal. **Brief-specific anti-creep:** no graded / corner-refined meshes on Problem 5 — they would recover rate 2 and defeat the inverted assertion that is the whole point. The callable-Dirichlet extension does not pull inhomogeneous-Neumann along; that stays deferred.
