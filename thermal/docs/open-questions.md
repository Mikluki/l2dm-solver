# Open questions

Live questions about the project, with the resolution path. Closed questions are deleted, not archived.

## Acceptance #4 of submission 0001 — "mean shift" check is ambiguous

The brief reads: *"Computed solution and exact solution agree on mean to within
solver tolerance (~10⁻¹⁰ relative) at the finest mesh — confirms the pin is in
the right place and the assembly is unbiased."*

Two problems surfaced during the 0001 review:

1. **The original test did not check the mean.** It asserted
   `solution[pin_dof] ≈ exact(pin_point)` to 1e-10. The solver writes
   `x[pin_dof] = pin_value` before `condense(...)`, scikit-fem reinserts that
   exact value, and the diff is 0 by construction — the assertion could not
   fail.

2. **The brief's 1e-10 tolerance is unattainable for the mean.** `∫T_h dA` is
   a discretization quantity (O(h²·|u|_H²), measured ~1.2e-3 at h=0.05). No
   choice of pin or assembly puts it at solver tolerance in finite-precision
   FE.

**Resolution applied:** (a) — the test now asserts `|∫T_h dA - ∫T_exact dA|`
falls at FE-error scale (calibrated against the finest-level L² error), which
is an orthogonal signal to the L² rate (catches sign-flipped contributions
that still converge). The corner-location assertion (`pin_x == 0.0`) survives
unchanged — it remains the real check that `_nearest_node_dof` worked.

**To close this question:** the brief's acceptance #4 text needs rewording to
match the implemented check. That edit is human-owned (the brief is
`accepted`), so this remains open until the brief is updated.
