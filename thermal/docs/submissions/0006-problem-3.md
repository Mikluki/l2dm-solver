# Submission 0006 — Problem 3 (radially symmetric disk in disk)

**Status:** done
**Predecessors:** 0004 done.
**Successors:** 0008 (Problem 4 — multi-disk extends the disk-in-disk OCC pattern); 0007 (Problem 5 — OCC polygon + named-curve Dirichlet lookup); any Part 2 curved-substrate geometry.

## Goal

Implement `verification.md` § Problem 3 and the tests that exercise it. First submission to exercise a **curved material interface** ($r = R_0$ as a mesh feature edge), a **curved Dirichlet boundary** ($r = R_{\text{out}}$ with $T = 0$), and the **gmsh OCC kernel + `fragment`** pattern that all subsequent multi-region geometry builders follow. Structural assertion mirrors Problem 2: probe the outer region at a fixed point and confirm the value is independent of $\kappa_1$.

## Acceptance — outcome

1. **Tests pass — Problem 3 and Problems 1–2.** Existing tests not regressed. ✓
2. **Convergence honest.** ≥ 3 mesh refinements; L² rate within [1.8, 2.2]. H¹ rate logged, not asserted. ✓
3. **κ₁-independence holds at $r = 0.6$.** Spread $(\max - \min)/|\text{mean}| < 1\%$ across $\kappa_1 \in \{0.1, 1, 10\}$, and each probed value within 5% of $T_{\text{out}}(0.6) = (q_0 R_0^2 / 2\kappa_2)\ln(R_{\text{out}}/0.6) \approx 2.299\cdot10^{-3}$. ✓
4. **Subdomain tagging is loud-fail** (swap inner/outer κ tags → rate breaks; bundle emitted). ✓
5. **Outer-Dirichlet wiring is loud-fail** (perturb outer-boundary Dirichlet by O(1) → rate collapses). ✓
6. **Inner-vs-outer Dirichlet confusion is loud-fail** (attach Dirichlet to inner circle → identifiable failure). ✓
7. **Mesh cache hits across the κ₁ sweep.** κ is not geometry; cache key is `{R_inner, R_outer, mesh_size}` only. ✓
8. **Predicted first-run failure recorded** (`ModuleNotFoundError` on the not-yet-written problem; an unsynchronised-`fragment` error from gmsh was the secondary risk). ✓

## What shipped

```
src/geometry/disk_in_disk.py
src/problems/problem_03_disk.py
tests/test_problem_03.py
src/geometry/rectangle_split.py — migrated to OCC kernel as a precursor (handed off without a brief)
docs/architecture.md § Key decisions — "Geometry kernel" entry added (OCC throughout; tags from outDimTagsMap, never hard-coded)
```

Worker conventions baked in: post-`fragment` surface tags recovered from `outDimTagsMap`, never via centroid classification (ADR-0003 anti-pattern); subdomain assignment by gmsh physical-surface **name** (`"inner_disk"`, `"outer_annulus"`); outer circle as a named physical curve `"outer_boundary"`; inner circle is an interior interface, no BC attached; no nullspace pin (Dirichlet removes the kernel). Unstructured Delaunay mesh — no transfinite engineering on the disk. The radial-symmetry diagnostic (max-over-θ of $T_h(r) - \overline{T_h}(r)$) stays a diagnostic artifact, not a load-bearing assertion — its job is to separate mesh-asymmetry signal from convergence-rate signal, not to gate acceptance.
