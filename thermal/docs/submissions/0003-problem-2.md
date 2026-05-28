# Submission 0003 — Problem 2 (piecewise-constant κ, 1D slab)

**Status:** in-progress
**Predecessors:** 0001 (accepted). 0002 § "Problem 2" verified the exact solution algebraically.
**Successors:** Problem 3.

## Goal

Implement `verification.md` § Problem 2 and the tests that exercise it. First submission to bite on the **P0 κ-by-subdomain** assembly (ADR-0004, `architecture.md` § Coefficient handling), a **mesh-aligned material interface** (ADR-0003), the **Dirichlet BC code path**, and the **κ₂-independence structural assertion**.

## Relevant core-doc sections

- `verification.md` § Problem 2 (the spec, in full — do not re-state in code).
- `architecture.md` § Coefficient handling (the `w.kappa` pattern, finally exercised), § Nullspace handling rule 1 (no pin when Dirichlet is present).
- ADR-0003, ADR-0004.
- 0002 § "Problem 2" — independent algebraic derivation. The continuity-based proof of κ₂-independence is there; do not re-derive.

## Decisions resolved before implementation

1. **Subdomain identifier is a string name.** scikit-fem's `mesh.subdomains` is name-keyed; Protocol's `kappa(subdomain_tag: int)` becomes `kappa(subdomain_name: str)`. Problem 1's unused parameter is renamed in-place; no behaviour change. The forward-looking `int` type in `protocol.py` and the `mesh.subdomains_per_element`-styled wording in `architecture.md` § Coefficient handling get corrected to match — a one-line doc edit, not a structural change, so no ADR.
2. **BC spec: `DirichletBC(value: float)` dataclass.** `boundary_conditions()` returns `dict[int, DirichletBC]`. Unlisted tags = natural zero-flux Neumann. Callable Dirichlet (Problem 5) and inhomogeneous Neumann are deferred until a problem actually needs them.
3. **κ₂-independence acceptance is two-pronged.** At the finest mesh, evaluate $T_h(0.75, h/2)$ for $\kappa_2 \in \{10, 100, 1000\}$:
   - $(\max - \min)/|\text{mean}| < 0.01$ — the literal 1% from `verification.md` (the κ₂-independence signal).
   - each value within 5% of the exact $q_0/(8\kappa_1) = 0.125$ — catches a coherently-wrong constant (e.g., sign-flipped source) that the spread alone would miss. Tolerance calibrated at FE-error scale for the finest mesh.

Existing-doc mandates the worker honours without re-deciding: mesh aligns with $x = 1/2$ (ADR-0003); subdomain assignment is by physical-surface tag, not by centroid coordinate (verification.md § Problem 2 failure diagnostic); no pin (Dirichlet on the left edge kills the nullspace per `architecture.md` § Nullspace handling rule 1).

## Acceptance

1. **Tests pass — Problem 2 *and* Problem 1.** The new P0 κ path and Dirichlet code path must be additive; no Problem 1 regression.
2. **Convergence honest.** ≥ 3 mesh sizes; fitted L² rate ≥ 1.8 *and* within 0.2 of 2.0. H¹ rate logged but not asserted (verification.md specifies only L²). A rate scraping 1.8 is suspect, not success.
3. **κ₂-independence holds** per decision 3 above. Pointwise evaluation is faithful interpolation, not nearest-node read-off.
4. **Subdomain tagging is loud-fail.** Forced-failure check, run once and removed before commit: swap the two κ values (or swap the subdomain tags). Confirm the κ₂-sweep test fails *and* the failure-artifact bundle is emitted; a clean rerun produces no artifact directory.
5. **Dirichlet wiring is loud-fail.** Forced-failure check, run once and removed: perturb the Dirichlet value applied at the left edge. Confirm the convergence test fails — error level shifts across all mesh sizes — and the bundle is emitted.
6. **Mesh cache hits across κ₂.** Sweeping κ₂ at fixed geometry produces no new `.msh` files after the first instantiation: κ is not a geometry parameter and the ADR-0007 key must reflect that.
7. **Predicted first-run failure (test-first cadence).** Before any new solver code is written, the Problem 2 tests fail for a specific, predicted reason (Dirichlet not implemented, or the rectangle-split spec not yet recognised by the geometry dispatch). A first-run pass is a red flag.

## Decisions left to the worker

Coding decisions, not pre-resolved here: file layout; source-field representation (quadrature-point evaluation vs. P0-by-subdomain); how κ₂ is parametrised on the Problem (constructor arg / class attr / per-sweep variant); the pointwise-evaluation mechanism for the κ₂-sweep test. CLAUDE.md's "if a small style choice starts looking like a physics choice, stop and surface" applies — coordinate-based subdomain assignment under the banner of convenience is the canonical example.

## Out of scope

Problems 3–5, callable Dirichlet, inhomogeneous Neumann, iterative solvers, validation against ⟨κ⟩, anything in `architecture.md` § Out of scope, any edit to `verification.md` or `physics.md`.

## Done

Acceptance passes; status moves to `accepted`; doc ambiguities exposed by implementation surface as proposed edits, not silent amendments. No new `.md` files at root.

---

## Implementation notes (worker)

Recorded once acceptance passed; status flipped to `accepted`. Worker-side choices that the brief left open:

- **File layout.** `src/geometry/rectangle_split.py` mirrors `unit_square.py`'s shape (spec dataclass + cache-backed `materialise`). `src/problems/problem_02_slab.py` is one file per Problem per `architecture.md` § Module structure. `tests/test_problem_02.py` parallels `tests/test_problem_01.py`. No new top-level files.
- **κ₂ parametrisation.** Constructor arg on a frozen dataclass (`Problem02Slab(kappa_2=...)`). The sweep test instantiates three variants; geometry is identical so the mesh cache is hit on every variant after the first.
- **Source-field representation.** Per-quadrature-point evaluation via `LinearForm` (the same `_make_load` pattern Problem 1 uses). Mesh alignment at `x = 1/2` (ADR-0003) means quadrature points never straddle the discontinuity, so a piecewise-constant `Problem.source(x, y)` is exact per element. Did not need a P0 source field; reserved for problems where the source is *not* mesh-aligned (none in Part 1).
- **Pointwise evaluation.** `basis.probes(point) @ solution` — sparse interpolation matrix from scikit-fem, which gives faithful P1 interpolation at the probe point, not a nearest-node read-off. Brief acceptance #3 explicitly flagged nearest-node as the foot-gun.
- **DirichletBC dataclass.** Lives in `src/problems/protocol.py` next to the Protocol it composes with. `dict[str, DirichletBC]` for `boundary_conditions()`; unlisted names are natural zero-flux Neumann.
- **Mesh structure.** Transfinite (structured) triangulation on the rectangular slab. Initial unstructured gmsh meshing gave an honest *floor* but noisy *fit* (rate 2.3 ± 0.5 across runs, threshold-window-violating despite the underlying physics being correct, because the exact solution is piecewise polynomial of degree ≤ 2 and the L² error sits near the floating-point-noise floor by `h ≈ 0.02`). Switching to `setTransfiniteCurve` / `setTransfiniteSurface` removed mesh-quality variation as a confound and recovered a clean rate of 1.97 ± 0.02. Transfinite meshing is not a physics choice — it is a *meshing* choice that satisfies ADR-0003's alignment requirement strictly more cleanly than unstructured triangulation. No ADR needed.
- **Pin/Dirichlet exclusivity.** Solver enforces "Dirichlet present ⇒ no pin allowed" defensively: a `Problem02Slab.pin_point()` returning `(x, y)` while `boundary_conditions()` returns a Dirichlet edge would raise immediately — rule 1 of `architecture.md` § Nullspace handling is checked, not assumed.

### Final convergence numbers (transfinite mesh, 5 levels):

| h     | n_dofs | L² error  | H¹ error  |
|-------|--------|-----------|-----------|
| 0.050 | 63     | 9.77e-08  | 3.23e-03  |
| 0.035 | 124    | 4.54e-08  | 2.15e-03  |
| 0.025 | 205    | 2.60e-08  | 1.61e-03  |
| 0.017 | 427    | 1.18e-08  | 1.08e-03  |
| 0.012 | 850    | 5.64e-09  | 7.68e-04  |

Fitted rates: **L² = 1.974, H¹ = 0.997.** Within 0.2 of both theoretical values.

### Forced-failure checks (both confirmed, reverted before commit)

1. **Subdomain-tag swap** (swap `kappa_1 ↔ kappa_2` in `Problem02Slab.kappa`). κ₂-sweep probe values became `{10.0: 0.0125, 100.0: 0.00125, 1000.0: 0.000125}`, dropping inversely with κ₂; spread = 2.68 (vs 0.01 threshold), max relative-to-exact = 0.999 (vs 0.05 threshold). Both prongs failed simultaneously; `tests/_artifacts/test_problem_02_kappa2_independence/kappa2_sweep.txt` was emitted. Clean rerun produced no artifact directory.
2. **Dirichlet value perturbation** (`DirichletBC(value=0.0)` → `DirichletBC(value=0.05)`). L² error became a *constant* 1.581e-2 at every mesh size (= 0.05 × √0.1, the integrated constant offset); fitted L² rate ≈ 0; convergence-rate assertion fired. All three artifact files (`errors.csv`, `rates.txt`, `error_field.png`) materialised in `tests/_artifacts/test_problem_02_converges/`. Clean rerun produced no artifact directory.

### Doc edits proposed in this submission

One-line edit to `docs/architecture.md` § Coefficient handling: the inline example previously referenced `mesh.subdomains_per_element` (no such attribute in scikit-fem 12) and iterated `tag for tag in ...`. Updated to iterate `mesh.subdomains.items()` and key by string name, matching the Protocol's new `kappa(subdomain_name: str)` signature. Pre-resolved by the brief as "a one-line doc edit, not a structural change, so no ADR" — not silent (this paragraph is the surface).

No edits proposed to `verification.md`, `physics.md`, or `open-questions.md`. No new ADRs.
