# Inspect — conventions and agent guidance

Read it before extending `scripts/inspect.py` for new problems. The rules below are not preferences; they're rules-from-incident.

For the visual legend itself (what each marker means), see the inline comments in `scripts/inspect.py` and the panel titles on the dashboard. This file is about *the principles* behind those choices.

---

## Rule 1 — Don't conflate solver tolerance with discretization scale

A check at "~1e-10 relative" is **solver-tolerance** scale. Quantities like `mean(T_h)`, `‖T_h − T‖_L²`, `|T_h − T|_H¹` are **discretization** quantities; they decay at `O(h^p)` with `p` the convergence rate, nowhere near 1e-10 for any reasonable mesh.

Symptom of getting this wrong: an agent writes `assert solution[pin_dof] ≈ pin_value` to satisfy a "small mean" acceptance. That assertion is a tautology — the solver wrote `pin_value` into that index before `condense()` and scikit-fem reinserts it on return.

**Pattern to use** when checking "assembly is unbiased":

```python
mean_h = float(_domain_integral.assemble(basis, u=basis.interpolate(solution)))
assert abs(mean_h - mean_exact) < 5 * finest_l2_error
```

The `5 *` is a looseness margin; tighten as you observe more problem instances. The bound calibrates against the L² error so it scales correctly with `h`.

---

## Rule 2 — Part 1 verification problems are non-dimensional

`physics.md` § Symbol glossary, `Q` row, says verbatim: *"W/m³ (formally; in practice dimensionless in the warm-up tests)"*. 

This means:
- `x`, `y`, `h` → `[-]`
- `T`, `T_h`, `T − T_h` → `[-]`
- `κ` → `[-]`
- `Q` → `[-]`

Do **not** label axes with `[m]` or `[K]` — that's a category error for the verification context. The unit-square domain has no physical length scale.

**On the convergence plot** the y-axis label stays generic (`error norm [-]`) because L² and H¹-seminorm have *different* physical dimensions in general (K·m vs K). Even though both happen to be dimensionless in Problem 1, labeling the y-axis with a specific norm's units would mislead the moment Problem 3+ ships with real SI κ. Keep the legend distinguishing the norms with math notation.

---

## Rule 3 — Every panel needs four pieces of text, not one

A panel with only an axis title is "random letters" to a first-time reader. The minimum:

1. **Title** — what the quantity is, with math notation (e.g. `source $Q(x, y)$`).
2. **Colorbar label** — the symbol + units (e.g. `$Q$  $[-]$`).
3. **x-axis label** — `$x$  $[-]$`.
4. **y-axis label** — `$y$  $[-]$`.

Use the `_label_axes(ax)` helper in `scripts/inspect.py` for (3) and (4) so the convention stays consistent. New panels should pass `title=` and `cbar_label=` to `_panel_field` rather than relying on defaults.

---

## Rule 4 — Uniform fields need a locked colorbar

`tripcolor`/`tricontourf` on data where `min == max` defaults to a ~±10% auto-range. The constant field then renders as a fake gradient and the colorbar reads as if there's variation.

**Pattern** (already implemented in `_panel_kappa`):

```python
is_uniform = np.isclose(kmin, kmax)
if is_uniform:
    pad = max(abs(kmin) * 1e-3, 1e-12)
    vmin, vmax = kmin - pad, kmin + pad
else:
    vmin, vmax = kmin, kmax
# ...pass vmin/vmax to tripcolor.
# When uniform, also set:
cbar.set_ticks([kmin])
cbar.set_ticklabels([f"{kmin:.4g}"])
```

This makes a uniform field look uniform. When Problem 2's κ contrast (1 vs 100) lands, the `else` branch takes over and shows the real range.

---

## Rule 5 — Theory-confirming patterns are not bugs

On a smooth manufactured solution, pointwise FE error tracks `|∇²T| · h²`. For Problem 1's `T = cos(πx)cos(πy)`, `|∇²T|² = 2π⁴ cos²(πx)cos²(πy)` — **maximal at corners, zero on the cross `x = 1/2 ∪ y = 1/2`**. The error plot reflects exactly this pattern. It is not a bug. Do not "fix" it.

Similarly, the pin at `(0,0)` enforces `error = 0` there and breaks the otherwise-symmetric error field across the diagonal. That asymmetry is **the pin doing its job**, not a defect.

If a future agent sees these patterns and proposes to "investigate" them, that's wasted effort. The corrective action is to make the patterns *read as expected* by annotating them (see Rule 6).

---

## Rule 6 — Annotate features that look surprising but aren't

Concrete, in the error panel:

```python
# Pin: error is 0 here by construction.
err_panel.plot(mesh.p[0, pin_dof], mesh.p[1, pin_dof],
               marker="*", markersize=14,
               color="white", markeredgecolor="black",
               linestyle="none", label=f"pin (DOF {pin_dof}, err=0)")

# |error| max location.
imax = int(np.argmax(np.abs(err_vals)))
err_panel.plot(mesh.p[0, imax], mesh.p[1, imax],
               marker="x", markersize=10, color="black", markeredgewidth=2,
               linestyle="none", label=f"max |err| = {abs(err_vals[imax]):.2e}")
```

For Problem 2+, add similar markers for whatever feature would otherwise read as suspicious — interface line, Dirichlet edge, anywhere a forced constraint creates visible asymmetry in a field.

---

## Rule 7 — Pyright "Import could not be resolved" is environmental noise

`scripts/inspect.py`, `tests/*.py`, and `src/*.py` will trigger Pyright `reportMissingImports` for `matplotlib`, `pytest`, `skfem`, etc. This is because Pyright doesn't see the `uv`-managed `.venv`. **Do not** add imports, type stubs, or `# type: ignore` comments to silence these. They are not real bugs; the code runs fine.

---

## Rule 8 — Tests and the inspector are different surfaces

Tests in `tests/` enforce correctness assertions; they fail loudly when something is wrong. The inspector in `scripts/inspect.py` is for human comprehension; it never asserts and never raises on "interesting" patterns.

Do not:
- Add asserts to the inspector (move them to a test).
- Add visualization to the test harness on a passing run (the
  `tests/_artifacts/` bundle is failure-only per `architecture.md` § Key decisions).
- Conflate the `artifacts/inspect/` tree (this tool) with the
  `tests/_artifacts/` tree (failure diagnostics).

---

# Layer 2 — solver invariant checks

Layer 1 looks at the problem (inputs / outputs). Layer 2 looks at the solver — whether assembly, condense, and solve preserved the invariants the math guarantees. It catches bugs that *pass the rate test*: "right rate, wrong constant" (verification.md § Problem 2 failure diagnostic).

Principles:
- **Numbers, not pictures.** Output is `internals.md`: tick-or-fail invariants with one number backing each verdict.
- **No exact solution required.** Layer 2 must work with only the Problem definition. That's the regime where the harness goes blind.
- **Silent on healthy.** Clean report = all `✓`. Anything else halts investigation.
- **Forward-loaded.** Problem 1 reports clean across the board; the value lands at Problem 2+.

Not in scope: condition number, spectral plots, per-DOF residual, quadrature-point dumps. None catch a bug the three techniques below don't already cover for direct-solve P1 FEM.

## Technique 1 — Subdomain integrity

**Catches.** (a) `gmsh:bounding_entities` metadata leak treated as a real subdomain — wrong-length κ array. (b) Wrong-side-of-interface κ in Problem 2: rate passes, constant is off by ~κ₂/κ₁. (c) Collapsed inclusion tags in Problem 4: two disks reported as one.

**Invariant.** `sum(area per real subdomain) == domain_area` to roundoff, summing only over names that don't start with `gmsh:`. Every iterated tag must have a `problem.kappa(tag)` value.

**Report.** One markdown table per problem:

| tag/name | n_elem | area | fraction | κ | metadata? |

plus a `Coverage: <sum> / <expected>  (Δ = …)` line. Domain area is Problem-side: unit square → 1.0, rectangle 0.1×1 → 0.1, disk → π R².

## Technique 2 — Source verification (`∫Q dA` two ways)

**Catches.** `LinearForm` integrating something other than the declared source. The canonical right-rate-wrong-constant load bug.

**Invariant.** FE-quadrature `∫_Ω Q dA` matches `problem.source_integral()` (new optional method, mirrors how `exact_gradient` is optional) to within quadrature scale (~`h⁴` for smooth Q under default P1 quadrature).

**Gotcha.** `source_integral()` must be a *one-liner*: a symmetry argument or a closed form, never a re-derivation. If the Problem author can't write one obviously-correct line, drop the check rather than fake it — a wrong analytic value silently agrees with a wrong FE value and the bug hides forever.

**Compute.** FE side: `Functional` over `problem.source(w.x[0], w.x[1])` on the P1 basis (same quadrature `_make_load` uses, so the check is an honest shadow).

## Technique 3 — Linear residual and matrix symmetry

**Catches.** Asymmetric assembly (BilinearForm typo); wrong-condense reinsert; silent `spsolve` failure.

**Invariants.**
- `‖K − Kᵀ‖_F / ‖K‖_F  <  ~1e-14`  (Laplacian and piecewise-κ Laplacian are self-adjoint when κ is per-element scalar).
- `‖A · u_h − b‖_∞      <  ~1e-10`  on the **enforced rows only** (see gotcha below).

**Gotcha — exclude substituted rows.** The pin row and any Dirichlet rows are *replaced* by `u_i = value` during `condense`; the original PDE equation `K[i, :] @ u = b[i]` at those rows is never solved. The residual at substituted rows lives at problem scale (~1e-9 for Problem 1), not solver scale, so including them in `‖A·u − b‖_∞` makes the threshold meaningless. Mask them out before taking the max.

**Compute.**
```python
sym = float(np.sqrt((K - K.T).multiply(K - K.T).sum())
            / np.sqrt(K.multiply(K).sum()))

enforced = np.ones(basis.N, dtype=bool)
if pin_dof is not None:
    enforced[pin_dof] = False
enforced[dirichlet_dofs] = False  # if any
res = float(np.max(np.abs((K @ u_h - b)[enforced])))
```
Sparse-aware, no dense allocation. On the enforced rows, the residual drops to machine epsilon × condition-number scale (Problem 1 at h=0.025: ~1.4e-15).

## Output

```
artifacts/inspect/<problem>/internals.md
```

Single file. Three sections — subdomain integrity, source verification, linear system invariants — each ticks or fails with the numeric backing. No extra PNGs at this layer.

## Where it lives (TBD)

Pick at implementation: new `scripts/diagnose.py` (separate CLI) or `--layer2` on `scripts/inspect.py` (one CLI, more flag-driven). Both write to the same per-problem directory.

---

## Adding a new problem to the inspector

1. Implement the `Problem` per `verification.md` § Problem definition    interface and the protocol in `src/problems/protocol.py`.
2. Register it in `scripts/inspect.py`:
   ```python
   PROBLEMS: dict[str, type] = {
       "problem_01": Problem01Manufactured,
       "problem_02": Problem02PiecewiseKappa,   # new
   }
   ```
3. The geometry-spec dispatch in `src/solver/solve_scalar.py:_materialise_geometry`    needs an `isinstance` branch for the new geometry's spec class.
4. The `_panel_kappa` panel reads `mesh.subdomains` and skips    `gmsh:*` metadata leaks — Problem 2's real two-subdomain map should    surface here automatically. Verify visually: the dashboard's κ panel    should show a clean two-color split, not a fake gradient.
5. If the Problem exposes an `exact_gradient` callable, the H¹-seminorm    path in `src/harness/study.py` picks it up. Otherwise it raises; that's    intentional (no silent projection fallback).
6. Run `uv run python -m scripts.inspect problem_NN` and confirm all six    panels render. Then write a pytest like `tests/test_problem_NN.py`.

---

## When this file is wrong

These rules are derived from review feedback on Problem 1. If a future problem's verification setup makes one of them counterproductive, **do not just break the rule** — raise it in `docs/open-questions.md` first with the specific scenario, propose a resolution, and update this file when the question closes. Silent rule-breaking is what created the loop this file exists to prevent.
