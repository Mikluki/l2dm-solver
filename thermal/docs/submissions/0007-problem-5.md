# Submission 0007 — Problem 5 (L-shape reentrant corner)

**Status:** accepted
**Predecessors:** 0006 (Problem 3 — OCC pattern, Dirichlet-DOF lookup); 0004 (L² norm at quadrature points — load-bearing for the rate-fitter this brief inverts an assertion against).
**Successors:** any future verification problem with inhomogeneous Dirichlet; Part 2 problems carrying non-constant BCs.

## Goal

Implement `verification.md` § Problem 5 and the tests that exercise it. The submission is a **harness-correctness** check, not a solver one: the inverted-assertion rate window $1.2 \le r_{L^2} \le 1.5$ confirms the rate-fitter honestly reports the corner-degraded $4/3$ convergence rather than masking it. The L-shape geometry is a single subdomain with $\kappa=1$, $Q=0$; the only solver-side novelty is **callable Dirichlet** — the protocol extension the `DirichletBC` docstring already names Problem 5 as the trigger for.

## Relevant core-doc sections

- `verification.md` § Problem 5 — geometry, exact solution, inverted-assertion semantics.
- `src/problems/protocol.py` § `DirichletBC` — the deferred-extension note that this submission discharges.
- `architecture.md` § Out of scope — *no AMR* and *no alternative element orders* are binding; this submission must not reach for either as an escape from a failing rate.
- 0004 § Why — the corrected norm is the mechanism by which an honest rate is even measurable; reverting it would mask Problem 5's signal.
- 0006 § Decisions 1, 6, 7 — OCC pattern, name-keyed subdomain assignment, named physical curves for boundary lookup. Mirror exactly.

## Decisions resolved before implementation

1. **Callable-Dirichlet protocol extension.** Extend `DirichletBC.value` to `float | Callable[[np.ndarray, np.ndarray], np.ndarray]`. The solver branches on type: scalar broadcasts as before; callable is evaluated at the boundary DOFs' physical coordinates and the resulting array is written into the constraint vector `x`. One dataclass, one field, type-dispatched in `_resolve_dirichlet_dofs`. Problems 1–3 keep working unchanged (their `value` is `float`).

2. **Inverted-assertion bounds: $1.2 \le r_{L^2} \le 1.5$, both sides asserted.** Upper-bound failure is the load-bearing inverted case (the rate-fitter is dishonestly reporting smooth convergence). Lower-bound failure means either the fitter is broken or the mesh is too coarse to be in the asymptotic regime — either of which invalidates Problem 5 as a harness check. Artifact bundle records *which* bound blew up so the diagnosis is one-step.

3. **Singular-solution polar convention.** Reentrant corner at $(1/2, 1/2)$. Define $\hat{x}=x-1/2$, $\hat{y}=y-1/2$, $r=\sqrt{\hat{x}^2+\hat{y}^2}$. Angle $\theta$ measured **clockwise** through the L-shape interior, $\theta=0$ along the edge $(1/2,1/2)\to(1,1/2)$, $\theta=3\pi/2$ along $(1/2,1/2)\to(1/2,1)$. Closed form:
   $$\theta = \begin{cases} -\operatorname{atan2}(\hat{y}, \hat{x}), & \hat{y} \le 0 \\ 2\pi - \operatorname{atan2}(\hat{y}, \hat{x}), & \hat{y} > 0 \end{cases}$$
   $T(x,y) = r^{2/3}\sin(2\theta/3)$ then vanishes on both edges meeting at the reentrant corner and grows like $r^{2/3}$ elsewhere.

4. **Boundary tagging.** All six edges as named physical curves. Two of them — the cut edges meeting at the reentrant corner — get `DirichletBC(value=0.0)` (scalar, because the exact solution is identically zero there). The remaining four — the outer perimeter of the L — get `DirichletBC(value=<callable>)`, where the callable returns $T(x,y)$ evaluated pointwise. Using the scalar form on the two cut edges (rather than a callable that happens to return zero) keeps the *callable* code path exercised only by the four nonzero edges, sharpening diagnosis when the new pathway misbehaves.

5. **Geometry builder: OCC polygon.** New `src/geometry/l_shape.py`. Six points → six lines → curve loop → plane surface. The reentrant corner is a polygon vertex, guaranteed a mesh node at every refinement (acceptance check #6). Single subdomain named `"interior"`, single named curve per edge (`"south"`, `"east_lower"`, `"cut_east"`, `"cut_north"`, `"west"`, `"north_left"` or similar — worker names them; only the two corner-vanishing edges and the four nonzero edges need distinguishing).

6. **Mesh sizes: 5 levels spanning ~8×.** Slow $4/3$ decay needs more levels and a wider range than rate-2 problems to give the fitter signal-over-noise. Representative list: `[0.12, 0.07, 0.040, 0.022, 0.013]`. Worker may tune within these bounds (finest comparable to Problem 2's $h \approx 0.012$). If a tuned list shifts the fitted rate outside $[1.2, 1.5]$ for honest reasons, surface — do not retune to chase the window (`CLAUDE.md` § "When a decision needs making mid-task").

7. **No nullspace pin.** Every edge carries Dirichlet → `pin_point()` returns `None`. The pin/Dirichlet exclusivity guard in the solver remains the loud check.

8. **Single subdomain, $\kappa=1$, $Q=0$.** No coefficient-indirection novelty. `kappa("interior") = 1.0`, `source(x, y) = 0`.

## Decisions left to the worker

- **Mesh-size list within the Decision 6 bounds.** If the rate-fit at the finest level looks anomalous, surface; do not silently retune.
- **`exact_gradient` for the H¹ norm.** $\nabla T$ diverges at the corner like $r^{-1/3}$ but is finite at every quadrature point (no quadrature scheme places a node *at* the corner for any mesh in Decision 6's range). Implement it with an `np.where(r > 0, ..., 0.0)` guard; the H¹ rate will be ~$2/3$, logged-not-asserted (Acceptance #3). Alternative: omit `exact_gradient` and have the study skip H¹. Implementation is the lower-risk path — the gradient is informative even if not asserted, and omitting it would mean special-casing the harness's hasattr-driven H¹ path.
- **Shared-corner consistency.** The existing `_resolve_dirichlet_dofs` "silently take first value at a shared corner" comment becomes load-bearing here. For Problem 5 every edge's callable returns the same $T(x,y)$ at any geometric corner, so the silent-take-first remains fine — but adding a `max-disagreement < 1e-10` assertion at this point would catch a future inhomogeneous-BC bug cheaply. Worker may add it; not load-bearing for this brief's acceptance.

## Acceptance

1. **Tests pass — Problem 5 and Problems 1–3.** Existing tests not regressed. The `DirichletBC.value` extension must not break the scalar-value path.

2. **Inverted-rate window holds:** $1.2 \le r_{L^2} \le 1.5$ over ≥ 5 mesh refinements. **Both bounds load-bearing.** Failure-mode mapping:
   - $r_{L^2} > 1.5$: rate-fitter is dishonest (likely the norm regressed, or the corner is being smeared by an over-coarse mesh that hides the singularity).
   - $r_{L^2} < 1.2$: rate-fitter is broken, or the finest mesh is not in the asymptotic regime.
   Artifact bundle records which bound blew up, the full error-vs-$h$ table, and the finest-mesh error field.

3. **H¹ rate logged, not asserted.** Expected ~$2/3$. `verification.md` is silent on H¹ here; asserting would bind the harness on an unspecified rate.

4. **Forced-failure on the inverted-assertion *mechanism*.** Build a synthetic `StudyResult` with `l2_rate ≈ 2.0` (e.g., copy Problem 1's). Assert that running it through the rate-window check raises. **The main test and this forced-failure test must call the *same* check function** (extracted helper, not inline-duplicated logic) — otherwise the forced-failure verifies a parallel copy of the check, and copy-drift defeats the whole point. Without this acceptance item, Problem 5 only catches *real* corner singularities, not bugs in the inversion logic itself. Single dedicated test in `tests/test_problem_05.py`.

5. **Callable-Dirichlet wiring is loud-fail.** Forced-failure (run once before commit): replace the four nonzero-edge callables with `DirichletBC(value=0.0)`. The solution now satisfies a *different* PDE — zero Dirichlet on every edge with $Q=0$, whose unique solution is $T \equiv 0$. The L² *discrepancy* against the singular reference no longer decays at any rate; the rate-window assertion fails on the lower bound (rate ≈ 0). Artifact bundle emitted; clean rerun produces no artifact directory.

6. **Reentrant corner is a guaranteed mesh node.** For each refinement level, locate the node closest to $(1/2, 1/2)$ and assert $\|p - (1/2, 1/2)\|_\infty < 10^{-12}$. The L-shape singular solution is meaningless if the corner drifts. Cheap, sharp.

7. **Mesh cache hits across re-runs.** Per the geometry-only-key convention (ADR-0007 / 0006 § Decisions 5). No new `.msh` files after the first cold pass.

8. **Predicted first-run failure recorded.** Most likely: `ModuleNotFoundError` on `problem_05_lshape` before the file exists, or `TypeError` from `float(spec.value)` in `_resolve_dirichlet_dofs` the first time a callable BC reaches it (the existing code does not yet dispatch on type). A test that passes on first run is suspicious.

## Pre-implementation checkpoint

- **Codebase prerequisite:** 0006 accepted; Problems 1–3 pass their full suites.
- **Protocol extension scope.** Decision 1 is the only change to `protocol.py`. Confirm the worker treats it as a small extension (one dataclass field's type union), not a refactor.
- **No AMR, no P2.** `architecture.md` § Out of scope binds — Problem 5 will *not* recover rate 2 on a uniform P1 mesh, and that is the point.

## Out of scope

- **Graded / refined-toward-corner meshes.** Would recover rate 2 and defeat the inverted assertion.
- **AMR / P2 / iterative solvers / GPU / parallelism.** `architecture.md` § Out of scope.
- **Other singular-solution benchmarks** (Motz, slit domain). One L-shape suffices for the inversion proof.
- **Problem 4.** Tracked separately; multi-disk composition is independent of corner-singularity rate honesty.
- **Inhomogeneous Neumann.** Still deferred — no problem in the queue exercises it. The callable-Dirichlet extension does not bring it along.
- **Shared-corner Dirichlet consistency check as a load-bearing assertion** — see § Decisions left to the worker. Not required for this submission.

## Done

Acceptance passes; status moves to `accepted` per `_conventions.md` § Post-accept compaction (head matter retained; "What shipped" added; acceptance lines marked ✓ or superseded; forced-failure logs and convergence tables moved into the commit that flipped status).
