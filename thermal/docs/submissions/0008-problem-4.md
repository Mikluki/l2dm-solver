# Submission 0008 — Problem 4 (two well-separated disks)

**Status:** done
**Predecessors:** 0006 (Problem 3 — disk-in-disk OCC pattern; the Problem 3 closed-form solution is what gets superposed); 0007 (Problem 5 — callable Dirichlet extension; protocol surface mirrored).
**Successors:** Part 2 multi-patch geometries (graphene + metal contacts, periodic metasurface) — both reuse the OCC three-disk `fragment` pattern and the convergence-in-geometric-parameter harness pattern this brief pioneers.

## Goal

Implement `verification.md` § Problem 4 and the tests that exercise it. Verifies **multi-region subdomain composition** — two disconnected $\kappa_1$-regions inside a $\kappa_2$ annulus, mesh-aligned with both inner circles, each region carrying its own source. The harness gains a **convergence-in-geometric-parameter** sweep pattern that Part 2's effective-conductivity-vs-separation tests will reuse.

## Acceptance — outcome

1. **Tests pass — Problem 4 and Problems 1–3, 5.** ✓
2. **Monotone-decreasing L² discrepancy across the sweep** (Config A → B → C, with $R_0/d_{\text{sep}} = 0.20, 0.10, 0.05$ at fixed $d_{\text{sep}}/R_{\text{out}} = 0.25$). ✓
3. **Finest config meets the 10% absolute threshold.** Config C discrepancy ≤ 10% per `verification.md` § Problem 4. ✓
4. **Mirror symmetry holds at fixed config B.** $|T_h(+x_i,y_i) - T_h(-x_i,y_i)| / \max_i |T_h(\pm x_i, y_i)| < 1\%$ at 4 symmetric probe pairs. Load-bearing structural assertion: honest mesh asymmetry is $O(h^2)$ while multi-region composition bugs break it to $O(1)$. ✓
5. **Forced-failure on composition is loud.** Setting $\kappa(\text{inner\_disk\_B}) = 100$ broke mirror symmetry by $\gg$ 1%; artifact bundle emitted; clean rerun left no artifact dir. ✓
6. **Mesh cache hits within a config.** Three configs → three distinct cached meshes (different `R_inner`); re-runs hit cache. ✓
7. **Per-config diagnostic table printed on every run.** Logs $R_0/d_{\text{sep}}$, $d_{\text{sep}}/R_{\text{out}}$, $h$, $n_{\text{dofs}}$, Problem 4 discrepancy, and a single-disk Problem 3 FE-error proxy at the same $(R_0, h)$ so a reader can separate FE error from approximation error. ✓

## What shipped

```
src/geometry/two_disks_in_disk.py
src/problems/problem_04_two_disks.py
src/problems/problem_03_disk.py — R_inner parameterised (default 0.3 preserves Problem 3 tests)
tests/test_problem_04.py
```

Key decisions baked in: reference approximation is the shifted single-disk superposition from Problem 3 (no monopole+dipole correction, no numerical fine-mesh reference — both are solver-cleverness creep against a deliberately loose 10% acceptance); the harness's `run_refinement_study` is **not** used because the discrepancy saturates against the approximation floor — three configs are tested directly at one fixed $h = 0.05$ each; two distinct inner-disk subdomain tags (`inner_disk_A`, `inner_disk_B`, both κ₁) so the mirror-symmetry forced-failure check is sharp; OCC three-disk `fragment` with surface tags recovered from `outDimTagsMap`. **Brief-specific anti-creep:** no general `parameter_sweep` harness driver until a second user shows up — Part 2's effective-κ-vs-separation test is the candidate; if and when it lands, that's the earned abstraction.
