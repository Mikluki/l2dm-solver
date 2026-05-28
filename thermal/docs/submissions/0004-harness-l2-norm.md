# Submission 0004 — Harness L² norm correctness

**Status:** done
**Predecessors:** 0003.
**Successors:** every subsequent verification problem; the corrected norm is what they measure against.

## Goal

Fix `src/harness/norms.py:l2_error` to integrate $(T_h(x,y) - T_{\text{exact}}(x,y))^2$ at quadrature points using the analytic `exact_solution` callable, **not** its P1 nodal projection. The harness is the project's deliverable; this was a harness correctness defect 0003 surfaced and did not resolve.

## Why it mattered

`l2_error` had been measuring $\|T_h - I_h T_{\text{exact}}\|$, not $\|T_h - T_{\text{exact}}\|$. For Problem 2 on a structured mesh the 2D Galerkin system collapses to a 1D 3-point stencil by y-symmetry; 1D P1 against a piecewise-quadratic source is **nodal-exact** (classical super-convergence), so $T_h \equiv I_h T_{\text{exact}}$ and the measured "error" was round-off (5.6e-9 at h=0.012, three orders below the honest `~0.020·h² ≈ 3e-6`). 0003's rate of 1.974 was satisfied accidentally — no future problem's L² rate could be trusted on any mesh where nodal super-convergence might apply.

## Acceptance — outcome

1. **Tests pass under the corrected norm.** Problem 1 and Problem 2 both pass; rates honest. ✓
2. **Magnitudes physically plausible.** Problem 2's L² error sits in the $10^{-6}$–$10^{-5}$ range at fine $h$, matching the back-of-envelope $0.020 \cdot h^2$. Problem 1's magnitudes shift only within constants (no super-convergence on a smooth exact). ✓
3. **Forced-failure on the harness itself.** Reverting `l2_error` to the nodal-projection form snapped Problem 2's L² magnitudes back to round-off scale (~1e-8 at h=0.05); reverting again restored the honest measurement. ✓
4. **No regressions in `h1_seminorm_error`** (its analytic-gradient path was already correct and untouched). ✓

## What shipped

```
src/harness/norms.py — l2_error mirrors h1_seminorm_error's idiom: basis.global_coordinates()
                       → call analytic `exact` at quadrature points → assemble (uh - exact_at_qp)²
```

No new Protocol method, no new abstraction; the assembly form change is local to `norms.py`. Acceptance thresholds in `verification.md` (rate ≥ 1.8, within 0.2 of 2.0) bind unchanged — the fix changes the *measurement*, not the *target*. 0003's `accepted` status was left in place; the prior acceptance reflected the defective norm and is now historical.
