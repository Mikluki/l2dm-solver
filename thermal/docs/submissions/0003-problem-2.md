# Submission 0003 — Problem 2 (piecewise-constant κ, 1D slab)

**Status:** done
**Predecessors:** 0001 (accepted). 0002 § Problem 2 verified the exact solution algebraically (now in `docs/derivations/algebraic-verification.md`).
**Successors:** Problem 3; 0004 (harness L² norm correctness — surfaced by this submission).

## Goal

Implement `verification.md` § Problem 2 and the tests that exercise it. First submission to exercise the **P0 κ-by-subdomain** assembly (see `architecture.md` § Coefficient handling), a **mesh-aligned material interface**, the **Dirichlet BC code path**, and the **κ₂-independence structural assertion**.

## Acceptance — outcome

1. **Tests pass — Problem 2 and Problem 1.** ✓
2. **Convergence honest.** ≥ 3 mesh sizes; L² rate ≥ 1.8 and within 0.2 of 2.0. ✓ on transfinite mesh — fitted L² = 1.974, H¹ = 0.997 (convergence table in commit history). Caveat: the L² magnitude is dominated by nodal super-convergence under the then-current `l2_error`, which projected the exact solution onto P1 before subtracting. 0004 fixes the norm.
3. **κ₂-independence holds.** $T_h(0.75, h/2)$ for $\kappa_2 \in \{10, 100, 1000\}$: spread $< 1\%$ ✓, each within 5% of exact $0.125$ ✓.
4. **Subdomain tagging is loud-fail.** Forced-failure (swap κ tags): both prongs failed; artifact bundle emitted; clean rerun produced no artifact dir. ✓
5. **Dirichlet wiring is loud-fail.** Forced-failure (perturb left-edge value): L² became a constant 1.581e-2 across all mesh sizes; rate dropped to ≈0; bundle emitted. ✓
6. **Mesh cache hits across κ₂.** κ is not a geometry parameter; cache key reflects that. ✓
7. **Predicted first-run failure.** Initial run failed because Dirichlet was not yet implemented — expected. ✓

## What shipped

```
src/geometry/rectangle_split.py
src/problems/problem_02_slab.py
tests/test_problem_02.py
src/problems/protocol.py — DirichletBC dataclass added; kappa(subdomain_name: str)
```

Worker decisions worth knowing:
- **Transfinite (structured) mesh** on the slab. Unstructured gmsh meshing put L² error at the floating-point floor by `h ≈ 0.02`, producing a noisy fitted rate. Transfinite removed mesh-quality variation as a confound. Mesh-side choice, not physics. The deeper reason — that the L² norm was measuring a near-zero quantity — is what 0004 addresses.
- **Pointwise probe via `basis.probes(point) @ solution`** — faithful P1 interpolation, not nearest-node.
- **Pin/Dirichlet exclusivity** enforced defensively in the solver — `Problem.pin_point()` returning a value while a Dirichlet edge exists raises immediately.

One-line edit to `architecture.md` § Coefficient handling: code example updated to iterate `mesh.subdomains.items()` and key by string name (was referencing a non-existent `mesh.subdomains_per_element` attribute).

