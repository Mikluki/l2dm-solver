# Submission 0008 — Problem 4 (two well-separated disks)

**Status:** proposed
**Predecessors:** 0006 (Problem 3 — disk-in-disk OCC pattern; the Problem 3 closed-form solution is what gets superposed); 0007 (Problem 5 — callable Dirichlet extension; not used here but the protocol surface is the same).
**Successors:** Part 2 multi-patch geometries (graphene + metal contacts, periodic metasurface) — both reuse the OCC three-disk `fragment` pattern this brief introduces and the convergence-in-geometric-parameter harness pattern it pioneers.

## Goal

Implement `verification.md` § Problem 4 and the tests that exercise it. The new code path under verification is **multi-region subdomain composition**: two disconnected $\kappa_1$-regions inside a $\kappa_2$ annulus, mesh-aligned with both inner circles, each region carrying its own source. The harness gains a **convergence-in-geometric-parameter** pattern that Part 2 will reuse for effective-conductivity-vs-separation tests — Problem 4's value to the harness is establishing that pattern cleanly, on top of the solver-verification value of confirming the composition path.

## Relevant core-doc sections

- `verification.md` § Problem 4 — geometry, materials, BCs, source, the two-mechanism approximation analysis, acceptance.
- `verification.md` § Problem 3 — closed-form solution being superposed; the structural identity $T_{\text{out}}(r) \propto \ln(R_{\text{out}}/r)$ that the superposition inherits.
- 0006 § Decisions 1, 6, 7 — OCC `fragment` pattern, name-keyed subdomain assignment, named physical curves. Mirror exactly.
- `architecture.md` § Coefficient handling — the P0-kappa indirection. Two inner subdomains with the same $\kappa_1$ exercise the per-element table at distinct subdomain *names*.
- `architecture.md` § Geometry kernel (Key decisions) — OCC `fragment` with conforming shared edges; tags recovered from the returned map, never hard-coded.

## Decisions resolved before implementation

1. **Reference approximation: shifted single-disk superposition** (verbatim per `verification.md` § Problem 4). The reference is
   $$T_{\text{ref}}(\mathbf{r}) = T_{\text{single}}(|\mathbf{r}-\mathbf{r}_A|; R_0, \kappa_1, \kappa_2, R_{\text{out}}) + T_{\text{single}}(|\mathbf{r}-\mathbf{r}_B|; R_0, \kappa_1, \kappa_2, R_{\text{out}})$$
   where $T_{\text{single}}$ is the Problem 3 closed form and **$R_{\text{out}}$ is the joint outer radius**, not a per-disk one. The two stacked approximation errors — $O(R_0/d_{\text{sep}})$ finite-separation and $O((d_{\text{sep}}/2R_{\text{out}})^2)$ joint-boundary residual — are intentional and the loose 10% acceptance accommodates them. No monopole+dipole correction, no numerical fine-mesh reference. Adding either is solver-cleverness creep (`CLAUDE.md` § "Don't add features").

2. **Harness routing: parameter sweep, not h-refinement.** Problem 4's convergence variable is the geometric ratio $R_0/d_{\text{sep}}$, not $h$. The test calls `solve_scalar` directly at one fixed `mesh_size` per configuration and iterates over a sweep of geometric configurations. `expected_rate()` returns `float("nan")`; `mesh_sizes()` returns a single-element list with the fixed FE resolution. **`run_refinement_study` is not used** — at a fixed Problem 4 config, the discrepancy $\|T_h - T_{\text{ref}}\|_{L^2}$ *saturates* (rate-2 in $h$ until FE error drops below the approximation floor, then rate $\approx 0$), so a single-slope least-squares fit would average the two regimes and report a meaningless number. Two-regime detection or Richardson extrapolation are real harness work with no second user in Part 1; if Part 2 needs a parameter-sweep driver, that's the second user that earns the abstraction.

   **Per-config FE-error proxy logged alongside the discrepancy.** The discrepancy alone conflates FE error and approximation error. At each config, also solve a *single-disk* Problem 3 at the same $R_0$ and same $h$, against its closed-form analytic solution, and log its relative $L^2$ error. That number is the FE-error scale at that $(R_0, h)$; the Problem 4 discrepancy minus that scale isolates the approximation contribution. Diagnostic only — no assertion on it. Requires `Problem03Disk` to accept `R_inner` as a constructor argument (currently hardcoded at $0.3$); see § Decisions left to the worker.

3. **Structural assertion: mirror symmetry of $T_h$ about the $y$-axis.** With identical disks at $(\pm a, 0)$ and identical $\kappa/Q$ in both, the exact PDE solution is even in $x$. The discrete solution should be mirror-symmetric to within a tolerance set by mesh asymmetry. Implementation: probe $T_h$ at 4 symmetric pairs (e.g. $(\pm 0.5, 0)$, $(\pm 0.5, 0.5)$, $(\pm 1.0, 0.3)$, $(\pm 1.5, 0)$ — worker picks the exact list, all interior to the outer disk and outside both inner disks) and assert $|T_h(+x_i,y_i) - T_h(-x_i,y_i)| / \max_i |T_h(\pm x_i, y_i)| < 1\%$ for each pair.

   **Why this is the load-bearing structural check.** It tests the *new code path* — multi-region composition — directly. A "skipped one subdomain" bug, a "tag collapse picks the wrong inner disk's source", or a wrong-side assembly indexing all break mirror symmetry to $O(1)$, while honest mesh asymmetry is $O(h^2)$ relative to $T$. The factor-100 separation is what makes 1% a sharp threshold.

4. **Sweep family: hold $d_{\text{sep}}$ and $R_{\text{out}}$ fixed, shrink $R_0$.**
   - Config A: $R_0 = 0.4,\ d_{\text{sep}} = 2,\ R_{\text{out}} = 8.\ R_0/d_{\text{sep}} = 0.20$, $d_{\text{sep}}/R_{\text{out}} = 0.25$.
   - Config B: $R_0 = 0.2,\ d_{\text{sep}} = 2,\ R_{\text{out}} = 8.\ R_0/d_{\text{sep}} = 0.10$, $d_{\text{sep}}/R_{\text{out}} = 0.25$. ← verification.md representative.
   - Config C: $R_0 = 0.1,\ d_{\text{sep}} = 2,\ R_{\text{out}} = 8.\ R_0/d_{\text{sep}} = 0.05$, $d_{\text{sep}}/R_{\text{out}} = 0.25$.

   $d_{\text{sep}}/R_{\text{out}}$ fixed at $0.25$ across the sweep so the joint-boundary residual stays approximately constant; the $R_0/d_{\text{sep}}$ trend is what the test measures. Holding $d_{\text{sep}},R_{\text{out}}$ literally fixed (rather than scaling them together with $R_0$) keeps the compute small without sacrificing the controlled trend.

5. **Geometry builder: OCC three-disk fragment.** New `src/geometry/two_disks_in_disk.py`. Three `addDisk` calls (one outer + two inner) → `occ.fragment([(2, outer)], [(2, innerA), (2, innerB)])` → conforming shared edges at both inner circles. Surface tags recovered from `outDimTagsMap`. Naming: two distinct inner-disk subdomains `"inner_disk_A"`, `"inner_disk_B"` (not collapsed to one tag), outer annulus `"outer_annulus"`, outer curve `"outer_boundary"`, inner curves `"inner_boundary_A"`, `"inner_boundary_B"` (tagged for inspection; no BCs attached).

6. **Distinct inner-disk tags, both κ₁.** `kappa("inner_disk_A") = kappa("inner_disk_B") = κ₁`. Distinct tags rather than collapsing both into one is the more future-proof choice — a future variant with asymmetric materials slots in without re-tagging the mesh. It also makes the **mirror-symmetry forced-failure check sharp**: change `kappa("inner_disk_B")` alone and the symmetry breaks loudly.

7. **Mesh-cache key:** `{R_inner, d_sep, R_outer, mesh_size}`. Per config, the geometry changes (`R_inner` varies), so the three configs produce three cached meshes. Within a config, re-runs hit cache. $\kappa$ is not geometry and stays out of the key (mirror 0006 § Decisions 5).

8. **Boundary tagging.** Only `"outer_boundary"` carries a BC: `DirichletBC(value=0.0)`. The two inner circles are interior material interfaces; no BC. No nullspace pin (`pin_point()` returns `None`). Mirror Problem 3's pattern.

9. **Materials per verification.md:** $\kappa_1 = 1$, $\kappa_2 = 10$, $q_0 = 1$ — same as Problem 3 so the superposed $T_{\text{single}}$ is literally Problem 3's solution at each disk center.

10. **Fixed FE mesh size across the sweep.** $h$ is fixed at a value fine enough that the FE error is much smaller than the 10% approximation-error threshold. Recommend $h = 0.05$ (informed by Problem 3's table: at $h \approx 0.025$ the relative L² error is well under 1%; at $h=0.05$ it is still in the few-percent range, comfortably below 10%). Worker may tune within $\{0.04, 0.05, 0.06\}$. The single-disk FE-error reference per config (§ Decisions 2) replaces what was previously a pre-commit one-shot sanity check — it is now logged as part of the test's diagnostic output.

## Decisions left to the worker

- **Probe-pair coordinates for the symmetry assertion.** Four pairs, all interior to the outer disk and outside both inner disks. The exact picks are worker's call; one near each disk's exterior + one in the deep annulus gives good coverage.
- **`exact_gradient` for the H¹ norm.** Optional — the harness's `run_refinement_study` is not used here, so the H¹ machinery is unexercised. Worker may omit `exact_gradient` entirely (no H¹ rate to log) or implement it as the sum of two shifted Problem-3 gradients. Omission is the lower-effort path and acceptable.

- **Parameterising `Problem03Disk.R_inner`.** Decision 2's per-config FE-error proxy needs Problem 3 evaluated at $R_0 \in \{0.4, 0.2, 0.1\}$, but `Problem03Disk._R_INNER` is currently hardcoded at $0.3$. Extend the dataclass with `R_inner: float = 0.3`, mirroring the existing `kappa_1` constructor pattern; the underlying `build_disk_in_disk(R_inner, R_outer)` already accepts the parameter, so the change is one field and threading it through `geometry()` and the cached `OUTER_PROBE_*` constants. Worker must verify that Problem 3's existing tests still pass at the default — the constants and the κ₁-independence acceptance must keep using $R_0=0.3$, $R_{\text{out}}=1.0$ unchanged.
- **OCC tag-recovery mechanism.** Map-based (preferred per 0006) or single-point classification — either acceptable provided it does not degenerate into point-in-polygon centroid classification (ADR-0003 anti-pattern).
- **L² discrepancy metric.** Either `l2_error` from the harness against `Problem04TwoDisks.exact_solution` (returning the superposition), or a one-off relative-L² computation. The harness path is preferred — it reuses the corrected norm from 0004 — but the test must remember that what's being measured is approximation+FE error, not FE error alone.
- **How to compute the relative L².** Normalize by $\|T_{\text{ref}}\|_{L^2}$ over the FE domain. The integral over the full annulus + inner disks gives a stable denominator that doesn't vanish in the sweep.

## Acceptance

1. **Tests pass — Problem 4 and Problems 1–3, 5.** Existing tests not regressed.

2. **Monotone-decreasing L² discrepancy across the sweep.** Relative L² of $T_h$ against $T_{\text{ref}}$ at fixed $h$: $\text{disc}_A > \text{disc}_B > \text{disc}_C$, strictly. Monotone alone is a necessary condition, not sufficient (a uniformly-too-loose constant could still pass).

3. **Finest config meets the absolute threshold.** $\text{disc}_C \leq 10\%$ per `verification.md` § Problem 4 acceptance. Both #2 and #3 load-bearing — #2 catches a stale or non-controlled approximation, #3 catches an approximation that's controlled but biased by a constant.

4. **Mirror symmetry holds at fixed config B.** $|T_h(+x_i, y_i) - T_h(-x_i, y_i)| / \max_i |T_h(\pm x_i, y_i)| < 1\%$ at each of 4 symmetric probe pairs. Honest mesh asymmetry is $O(h^2)$; a multi-region composition bug breaks this to $O(1)$. **This is the load-bearing structural assertion — see § Decisions 3.**

5. **Forced-failure on multi-region composition is loud.** Run once before commit: set `kappa("inner_disk_B") = 100` (asymmetric κ), confirm the symmetry assertion in #4 fails by $\gg$ 1%. Artifact bundle emitted; clean rerun produces no artifact dir. This validates that the symmetry assertion bites on the *new code path*.

6. **Mesh cache hits within a config.** Re-running the same config produces zero new `.msh` files. The three configs produce three distinct cached meshes (different `R_inner`).

7. **Per-config diagnostic table printed on every run** (not just on failure). For each config, log: $R_0/d_{\text{sep}}$, $d_{\text{sep}}/R_{\text{out}}$, $h$, $n_{\text{dofs}}$, Problem 4 discrepancy (relative $L^2$ vs $T_{\text{ref}}$), and the **single-disk Problem 3 FE-error proxy at the same $(R_0, h)$** (per § Decisions 2). The FE-error column lets a reader distinguish "approximation error dominates" (Problem 4 discrepancy $\gg$ Problem 3 proxy — the expected case) from "FE error dominates" ($h$ too coarse). On assertion failure, the table also lands in the artifact bundle alongside $\|T_h\|_{L^2}$ and $\|T_{\text{ref}}\|_{L^2}$ for the magnitude check.

8. **Predicted first-run failure recorded.** Most likely `ModuleNotFoundError` on `problem_04_two_disks` or `two_disks_in_disk` before the files exist; or an `out_map` shape surprise from `occ.fragment` with three input disks (the inner-disk count was 1 in 0006, is 2 here, and the per-input fragment count needs verifying — Problem 6 of the OCC `fragment` semantics, not a bug). A test that passes on first run is suspicious.

## Pre-implementation checkpoint

- **Codebase prerequisite:** 0006 accepted, 0007 accepted. Problems 1–3, 5 pass their full suites.
- **No protocol extension.** Problem 4 implements the existing `Problem` protocol with `expected_rate() = float("nan")` and `mesh_sizes() = [<single h>]`. No new methods. If the worker finds themselves wanting to add one, surface — that signals an architectural decision, not a worker call.
- **Sweep is *not* a refinement study.** `run_refinement_study` must not be invoked. If the worker is tempted to reuse it by abuse (e.g. passing the configuration index as a fake mesh size), stop and surface.

## Out of scope

- **Monopole + dipole correction to the reference.** Decision 1 commits to the superposition. Adding higher-order corrections is a derivation with no current user.
- **Numerical fine-mesh reference.** Decision 1 rejects it; it's circular (the solver under test would be the reference).
- **A `parameter_sweep` harness driver.** No second user in Part 1. If Part 2's first effective-κ test wants one, that's when it earns the abstraction.
- **Three or more inner disks.** Problem 4 specifies two. A future "$N$-patch composition" verification problem is conceivable but out of scope here.
- **AMR, P2, iterative solvers, GPU, parallelism, time dependence** — `architecture.md` § Out of scope binds.
- **Kapitza coupling, in-plane graphene conduction, validation against $\langle\kappa\rangle$** — Part 2 territory.

## Done

Acceptance passes; status moves to `accepted` per `_conventions.md` § Post-accept compaction (head matter retained; "What shipped" added; acceptance lines marked ✓ or superseded; the per-config discrepancy table and the forced-failure log moved into the commit that flipped status).
