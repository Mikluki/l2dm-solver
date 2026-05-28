# Submission 0006 — Problem 3 (radially symmetric disk in disk)

**Status:** accepted
**Predecessors:** 0004 done
**Successors:** Problem 4 (multi-disk extends the disk-in-disk OCC pattern); any Part 2 curved-substrate geometry.

## Goal

Implement `verification.md` § Problem 3 and the tests that exercise it. First submission to exercise a **curved material interface** ($r = R_0$ as a mesh feature edge), a **curved Dirichlet boundary** ($r = R_{\text{out}}$ with $T = 0$), and the **gmsh OCC kernel + `fragment`** pattern that all subsequent multi-region geometry builders will follow. The structural assertion mirrors Problem 2: probe the outer region at a fixed point and confirm the value is independent of $\kappa_1$.

## Relevant core-doc sections

- `verification.md` § Problem 3 — geometry, materials, BCs, source, exact solution, acceptance.
- `docs/derivations/algebraic-verification.md` § Problem 3 — independent re-derivation, clean.
- `architecture.md` § Coefficient handling, § Nullspace handling, § Key decisions (proposed addition: geometry kernel).
- 0003 § What shipped — patterns to mirror: `basis.probes(point) @ solution` for pointwise probing; pin/Dirichlet exclusivity enforced in the solver; geometry-only mesh-cache key.

## Decisions resolved before implementation

1. **Geometry kernel: gmsh OCC.** `addDisk(outer)` + `addDisk(inner)` + `occ.fragment([(2, outer)], [(2, inner)])`. `synchronize()` is called before any `addPhysicalGroup`; post-boolean surface tags are recovered from the `outDimTagsMap` returned by `fragment`. Hard-coded integer tags carried through the boolean are forbidden — they renumber. The geo kernel is **not** acceptable here; the `rectangle_split` → OCC migration consolidates the codebase on a single kernel as a precursor to this submission (handed off without a brief).
2. **Proposed addition to `architecture.md` § Key decisions:**
   > **Geometry kernel.** gmsh OCC kernel for all geometry builders. The geo kernel is no longer used; `rectangle_split` migrated to OCC as a precursor to 0006. Curved-boundary geometries rely on OCC `fragment` to produce conforming shared edges (ADR-0003 mesh-alignment). Tags after `fragment` are queried via the returned map, never hard-coded.

   Planner approves this addition before this brief moves to `accepted`.
3. **No structured/transfinite meshing on the disk.** Unstructured Delaunay (gmsh default) is acceptable. If the L² rate is noisy under the 0004-corrected norm, surface — do not engineer geometry to chase a rate. Transfinite-on-a-disk requires quadrilateral blocking that has no precedent here.
4. **Geometry parameters per verification.md:** $R_0 = 0.3$, $R_{\text{out}} = 1.0$, $\kappa_1 = 1.0$, $\kappa_2 = 10.0$, $q_0 = 1.0$. No tuning.
5. **Mesh-cache key:** `{R_inner, R_outer, mesh_size}`. $\kappa_1$ and $\kappa_2$ must not enter the key (mirror 0003 acceptance #6).
6. **Subdomain assignment** by gmsh physical-surface **name** (`"inner_disk"`, `"outer_annulus"`), never by element-centroid coordinate (ADR-0003 / 0003 § Decisions).
7. **Boundary tagging.** Outer circle as named physical curve `"outer_boundary"` → Dirichlet $T = 0$. The inner circle is an interior interface — tagged for inspection if convenient, no BC attached.
8. **No nullspace pin.** Dirichlet on $r = R_{\text{out}}$ removes the nullspace; `pin_point()` returns `None`. The pin/Dirichlet exclusivity check in the solver remains the loud guard.
9. **Structural assertion: κ₁-independence at $r = 0.6$.** Sweep $\kappa_1 \in \{0.1, 1, 10\}$ at fixed $\kappa_2 = 10$. Acceptance:
   - **Spread** $(\max - \min)/|\text{mean}| < 1\%$ at the finest mesh.
   - **Bound** — each probed value within 5% of $T_{\text{out}}(0.6) = (q_0 R_0^2 / 2\kappa_2)\ln(R_{\text{out}}/0.6) \approx 2.299\cdot10^{-3}$.

   Two prongs because spread alone misses uniformly-sign-flipped or scaled outputs; bound alone misses a κ₁-leak that happens to cancel between sweep values. If the finest-mesh probe is dominated by discretization noise rather than the structural signal, surface for retuning — do not silently widen.

## Decisions left to the worker

- **Mesh sizes list.** ≥ 3 sizes; lower bound set by perimeter resolution of $R_0 = 0.3$ (enough nodes around the inner circle for an honest rate fit); upper bound chosen so the list spans the asymptotic regime. If the list turns unusually fine (lower bound substantially below Problem 2's $h \approx 0.012$), surface — that signals a geometry issue, not a refinement decision.
- **`exact_gradient` at $r = 0$.** $\nabla T(0) = 0$ analytically; use `np.where` with a radius guard or an explicit branch. Either is fine.
- **OCC tag recovery.** Whether to recover post-fragment surface tags via the returned `outDimTagsMap` or via `getEntities(dim=2)` filtered by a single point on each disk. Map is cleaner; either is acceptable provided it does not degenerate into point-in-polygon classification (the ADR-0003 anti-pattern).
- **Radial-symmetry diagnostic** (recommended, not load-bearing). Recording max-over-θ of $T_h(r) - \overline{T_h}(r)$ at the finest mesh as a failure-only artifact catches mesh asymmetry as a signal separable from the convergence rate.

## Acceptance

1. **Tests pass — Problem 3 and Problems 1–2.** Existing tests are not regressed.
2. **Convergence honest.** ≥ 3 mesh refinements; L² rate $\geq 1.8$ and within $0.2$ of $2.0$. H¹ rate logged, not asserted (`verification.md` § Problem 3 specifies L² only).
3. **κ₁-independence holds.** Spread and bound prongs both pass at the finest mesh per § Decisions 9.
4. **Subdomain tagging is loud-fail.** Forced-failure (swap inner/outer κ tags): rate breaks; artifact bundle emitted; clean rerun produces no artifact directory.
5. **Outer-Dirichlet wiring is loud-fail.** Forced-failure (perturb outer-boundary Dirichlet by an O(1) amount): L² rate collapses; bundle emitted.
6. **Inner-vs-outer Dirichlet is loud-fail.** Forced-failure (attach Dirichlet to the inner circle by mistake): convergence fails in an identifiable way; bundle emitted. New for Problem 3 — both circles are named physical curves, easy to confuse.
7. **Mesh cache hits across κ₁ sweep.** $\kappa_1$ is not geometry; the sweep produces zero new `.msh` files after the first instantiation.
8. **Predicted first-run failure** recorded in the convergence-table acceptance log. Most likely `ModuleNotFoundError` on `problem_03_disk` before the file exists, or an unsynchronised-`fragment` error from gmsh. A test that passes on first run is suspicious.

## Pre-implementation checkpoint

- **Codebase prerequisite:** `rectangle_split.py` migrated to the OCC kernel; Problem 2 still passes its full acceptance suite at the same rate/magnitudes as 0003's recorded run. Without this, the § Decisions 2 architecture-doc addition is dishonest.
- **0004 status flip:** `proposed → accepted` on 0004's frontmatter (code change has already landed; norm-correctness check verified at `src/harness/norms.py`).
- **Architecture-decision approval:** § Decisions 2 above approved by the planner.

## Out of scope

- Problems 4–5.
- Kapitza coupling, in-plane graphene conduction, validation against ⟨κ⟩.
- P2 / isoparametric curved-boundary elements (locked out by `architecture.md` § Out of scope until Part 2 demands it).
- Iterative solvers; AMR; parallelism; GPU.
- The radial-symmetry diagnostic as a load-bearing assertion — it is a diagnostic artifact only.

## Done

Acceptance passes; status moves to `accepted` per `_conventions.md` § Post-accept compaction (head matter retained; "What shipped" added; acceptance lines marked ✓ or superseded; forced-failure logs and convergence tables moved into the commit that flipped status).
