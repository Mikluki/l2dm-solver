# Submission 0004 — Harness L² norm correctness

**Status:** proposed
**Predecessors:** 0003 (accepted, but its acceptance #2 "convergence honest" did not hold — see § Why).
**Successors:** every future verification problem; the corrected norm is what they measure against.

## Goal

Fix `src/harness/norms.py:l2_error` to integrate `(T_h(x,y) − T_exact(x,y))²` at quadrature points using the analytic `exact_solution` callable, **not** its P1 nodal projection. Re-verify Problem 1 and Problem 2 against the corrected norm.

The harness is the project's deliverable (`CLAUDE.md` § Project disposition). This is a harness correctness defect that 0003 surfaced and did not resolve.

## Why

`l2_error` currently does:

```python
nodal_exact = exact(basis.mesh.p[0], basis.mesh.p[1])
uex = basis.interpolate(nodal_exact)         # P1 interpolant of T_exact
val = _l2_squared.assemble(basis, uh=uh, uex=uex)
```

It measures `‖T_h − I_h T_exact‖`, not `‖T_h − T_exact‖`. The two differ by an O(h²) term whose constant depends on the problem.

For Problem 1 (smooth cosine) they agree to within constants and the rate is honest. For Problem 2 on a structured mesh, the 2D Galerkin system collapses to the 1D 3-point stencil by y-symmetry, and 1D P1 against a piecewise-quadratic source gives **nodal-exact** `T_h` (classical super-convergence). Therefore `T_h ≡ I_h T_exact` and the measured "error" is round-off. The reported L² of 5.6e-9 at h=0.012 is three orders of magnitude below the honest `~0.020 · h² ≈ 3e-6` from `‖T_exact − I_h T_exact‖` on a quadratic. 0003's rate of 1.974 was satisfied accidentally.

Without the fix, no future verification problem's L² rate can be trusted on any mesh where nodal super-convergence may apply.

## Relevant core-doc sections

- `verification.md` § Harness requirements — the contract this submission upholds.
- 0003 § Implementation notes — the convergence table whose magnitudes are the symptom.
- `src/harness/norms.py:h1_seminorm_error` — already uses the analytic-callable-at-quadrature-points pattern; the L² fix mirrors it.

## Decisions resolved before implementation

1. **Mechanism: mirror `h1_seminorm_error`.** That function already uses `basis.global_coordinates()` to get quadrature-point coordinates, calls the analytic callable there, and assembles the form. The L² fix uses the same idiom; no new Protocol method, no new abstraction.
2. **Thresholds unchanged.** `verification.md` § Problem 1 / § Problem 2 acceptance (rate ≥ 1.8, within 0.2 of 2.0) binds as written. The fix changes the *measurement*, not the *target*. If an honest measurement no longer satisfies the threshold, surface it; do not retune.
3. **Re-measurement is part of acceptance.** Record the new L² error tables and fitted rates for Problem 1 and Problem 2 in this brief's implementation notes (as 0003 did). If a rate falls outside the window for honest reasons (e.g. the structured-transfinite mesh + super-convergence still gives a faster-than-2 rate against the *correct* norm), that goes to `docs/open-questions.md` — do not paper over.

## Acceptance

1. **Tests pass under the corrected norm.** Problem 1 and Problem 2 tests pass. A rate that *was* passing under the old norm and *now* fails is the right outcome — surface, do not retune.
2. **Magnitudes are physically plausible.** Problem 2's L² error at h=0.012 sits in the 10⁻⁶–10⁻⁵ range (matches the back-of-envelope `0.020 · h²` from § Why). Problem 1's magnitudes shift within constants — no super-convergence on a smooth exact, so the new and old numbers should agree at the same order.
3. **Forced-failure on the harness itself.** Run once before commit: revert `l2_error` to the old nodal-projection form temporarily. Confirm Problem 2's L² magnitudes snap back to round-off scale (~1e-8 at h=0.05). Revert. The point: show the corrected norm bites on a signal the old one didn't.
4. **No regressions in untouched code.** `h1_seminorm_error` is untouched (its analytic-gradient path was already correct).

## Decisions left to the worker

The assembly shape of the corrected form: one new `Functional` that takes both `w.uh` and an `exact_at_qp` array, or an extension of the existing `_l2_squared`. Either is acceptable provided the quadrature scheme matches `h1_seminorm_error`.

## Out of scope

- Test threshold *values*.
- The H¹ path (already correct).
- The small 0003 follow-ups: `LevelResult.pin_dof` annotation, unused `y` params, unused `pytest` import, retroactive `int → str` amendment of 0003 § Decisions resolved 2. Bundle separately or roll into this submission's cleanup as the worker prefers — they are surgical and do not affect this submission's load-bearing acceptance.
- Anything in `architecture.md` § Out of scope.

## Done

Acceptance passes; status moves to `accepted`. If Problem 1 or Problem 2 fails the rate window under the corrected norm, status stays `in-progress` and the failure is logged in `docs/open-questions.md` with the corrected error table attached. 0003's status (currently `accepted`) reflects the prior acceptance under the defective norm; whether to revert it to `in-progress` is the planner's call and not gated by this submission.
