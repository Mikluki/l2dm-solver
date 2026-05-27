# Submission 0001 — Problem 1 end-to-end vertical slice

**Status:** accepted
**Predecessors:** none (first submission)
**Successors:** 0002 (Problem 2 — piecewise-constant κ) once this is accepted.

## Goal

Stand up the smallest vertical slice — geometry → mesh → solver → harness → pytest — that runs `verification.md` Problem 1 (smooth manufactured solution on a unit square, pure Neumann) and converges at the expected rates. The first time this test passes for the right reasons is the moment Part 1 starts paying off.

## Relevant core-doc sections

- `physics.md` § Equations → "What Part 1 solves" — the scalar PDE and why Part 1 solves it.
- `verification.md` § Harness requirements, § Problem definition interface, § Problem 1.
- `architecture.md` § Pipeline, § Key abstractions, § Module structure, § Coefficient handling, § Nullspace handling.
- ADR-0001 (stack), ADR-0002 (P1 elements), ADR-0004 (P0 κ), ADR-0005 (node-pinning), ADR-0006 (direct solve), ADR-0007 (mesh cache).

## Decisions resolved before implementation

These were settled in the planning round; they are not open for the implementer to revisit silently.

1. **Nullspace handling: node-pinning uniformly.** Per ADR-0005, `architecture.md` § Nullspace handling, and `verification.md` § Problem 1 (now consistent). Subtract-the-mean is not used.
   - The `Problem` Protocol exposes `pin_point() -> tuple[float, float] | None`. The solver requires it; there is no solver-side default.
   - For Problem 1 specifically: `pin_point()` returns `(0.0, 0.0)`. The exact solution there is $T = \cos(0)\cos(0) = 1$ — a guaranteed mesh node (geometry corner), reproducible across refinements, well-conditioned (nonzero, away from the saddle at the centroid).
2. **`Problem` interface as `typing.Protocol`.** Structural typing, no inheritance burden, plays cleanly with the eventual Part 2 integral-form variant. One Protocol in `src/problems/protocol.py`.
3. **Failure-only diagnostic artifacts wired now.** Error-vs-$h$ table + fitted rate + finest-mesh error-field plot, emitted to `tests/_artifacts/{test_id}/` only when the test fails. Suppressed on pass.
4. **Mesh cache stub wired now.** SHA-256 of geometry parameters → `.msh` in `tests/_mesh_cache/` per ADR-0007. The unit-square mesh is cheap to build, but the cache layer must exist so Problems 2–4 inherit it without retrofitting.
5. **$H^1$ semi-norm acceptance threshold ≥ 0.9** per `verification.md`. Pinned in the test code, not negotiable in this submission.

## Decisions deferred (not this submission)

- The Dirichlet BC code path. Problem 1 is pure-Neumann; Problem 2 introduces Dirichlet. Implement just enough BC handling for Problem 1; do not pre-build a Dirichlet abstraction.
- Per-element subdomain coefficient indirection (`w.kappa` over P0 field). Problem 1 is single-subdomain, $\kappa = 1$. A uniform scalar suffices here. Problem 2 introduces the P0-κ-per-tag pattern.
- Geometry builders for shapes other than the unit square.
- Artifact retention/cleanup policy beyond `.gitignore`.

## Deliverable

### Files created

```
src/problems/__init__.py
src/problems/protocol.py                    # typing.Protocol for Problem
src/problems/problem_01_manufactured.py     # Problem 1 instance
src/geometry/__init__.py
src/geometry/cache.py                       # SHA-256-keyed .msh cache (ADR-0007)
src/geometry/unit_square.py                 # gmsh model builder, cache-backed
src/solver/__init__.py
src/solver/result.py                        # SolverResult dataclass (solution, basis, mesh)
src/solver/solve_scalar.py                  # solve_scalar(problem, mesh_size) -> SolverResult
src/harness/__init__.py
src/harness/norms.py                        # l2_error, h1_error via skfem forms
src/harness/study.py                        # run_refinement_study(problem) -> StudyResult
src/harness/artifacts.py                    # failure-only emitter (table, rate, plot)
tests/conftest.py                           # pytest fixture: per-test artifact dir
tests/test_problem_01.py                    # wires Problem 1 through the harness
.gitignore                                  # adds tests/_artifacts/, tests/_mesh_cache/
```

Each source file begins with the `# ABOUTME:` header per CLAUDE.md.

### Files modified

- `pyproject.toml` — no expected changes; flag if any dep is unexpectedly missing.

### Artifact

A passing `uv run pytest tests/test_problem_01.py -v` reporting observed $L^2$ rate close to 2 and $H^1$ semi-norm rate close to 1.

## Acceptance

All must hold simultaneously.

1. **Test passes.** `uv run pytest tests/test_problem_01.py -v` succeeds.
2. **Convergence is honest, not threshold-skimming.** The refinement study uses ≥ 3 mesh sizes; fitted rates satisfy
   - $L^2$ rate ≥ 1.8 *and* within 0.2 of 2.0,
   - $H^1$ semi-norm rate ≥ 0.9 *and* within 0.2 of 1.0.
   A rate that scrapes 1.8 should be treated as suspect, not as success — refer back to `verification.md` § Problem 1 failure diagnostic.
3. **Node-pinning convention applied as specified.** The `Problem` exposes `pin_point()`; the solver pins the nearest-node DOF to the exact solution at that point; pinned-DOF index is reproducible across refinements (verified by an assertion in the test).
4. **Mean integral decays at FE-error scale.** At the finest mesh, $\left|\int_\Omega T_h\,dA - \int_\Omega T_{\text{exact}}\,dA\right|$ is bounded by a constant times the L² error. Concretely, the test asserts $|\int T_h\,dA| < 5\cdot\|T_h - T\|_{L^2}$ (note $\int T_{\text{exact}}\,dA = 0$ for Problem 1). This is an *orthogonal* signal to the L² rate — it catches sign-flipped assembly contributions that still converge — and it confirms both that the pin landed correctly and that the assembly is unbiased.

   Rationale: an earlier draft of this brief asked for "agreement on mean to within solver tolerance ($\sim 10^{-10}$ relative)". That phrasing was wrong: $\int T_h\,dA$ is a discretization quantity that decays at the L² rate (∼ $h^2$), not at the linear-solver tolerance. Asserting against 1e-10 would either always fail (any reasonable mesh) or be replaced by a solver-side tautology (writing the pin value back and asserting we wrote it). The implemented check is the right scale.
5. **Failure mode is visible.** Forced-failure check (run once, removed before commit): temporarily multiply the assembled $\kappa$ by 2 inside the solver. Confirm:
   - the test fails on the rate assertion,
   - `tests/_artifacts/test_problem_01/` is populated with the error-vs-$h$ table, fitted-rate value, and finest-mesh error-field plot,
   - no artifact directory is created on the passing run.
6. **Mesh cache is exercised.** Re-running the test produces no new entries in `tests/_mesh_cache/` (cache hit on every mesh size), confirmed by file mtimes.
7. **Predicted first-run failure (test-first cadence per CLAUDE.md).** Before any solver code is written, `pytest tests/test_problem_01.py` should fail with `ModuleNotFoundError` (or equivalent) referencing the not-yet-written `src/solver/solve_scalar.py`. If the first run fails for any *other* reason, the test wiring is wrong and must be fixed before implementing the solver. A first-run pass is a red flag and must be investigated.

## Out of scope for this submission

- Problems 2–5 and their geometry builders.
- The subdomain-tagged P0 κ indirection (deferred to Problem 2's submission).
- Dirichlet BC handling beyond what node-pinning needs.
- Iterative solvers (ADR-0006: direct via `scipy.sparse.linalg.spsolve`).
- Adaptive mesh refinement, parallelism, GPU, alternative element orders (`architecture.md` § Out of scope).
- Benchmarking or timing reports.
- Validation against $\langle\kappa\rangle$ (Part 2).
- Any structural change to `verification.md` or `physics.md` (human-owned per CLAUDE.md).

## Pre-implementation checkpoint

Before the next worker begins coding:

- Confirm the file list above is the minimum needed. If a file in the list turns out unnecessary during implementation, drop it and note the change in this brief's acceptance discussion — do not silently add files not listed.
- Confirm that `architecture.md` § Nullspace handling is the authoritative reference for the pin convention; the implementer reads it before writing `solve_scalar`.
- Confirm that the `Problem` Protocol matches `verification.md` § Problem definition interface field-for-field (including the new `pin_point()` accessor introduced by the nullspace section).

## Done definition

This submission is done when:

1. All acceptance criteria above pass.
2. `docs/architecture.md` is updated only if the implementation revealed an ambiguity in it (with a separate ADR if structural).
3. `docs/submissions/0001-problem-1.md` (this file) status changes from `proposed` to `accepted`.
4. No new `.md` files at root; no per-script docs; no findings docs.
