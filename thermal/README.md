# thermal

Scalar steady-state heat-equation solver + verification harness. See `docs/` for the source of truth — `physics.md`, `verification.md`, `architecture.md`, `CLAUDE.md`.

## What runs

Verification problems live under `tests/test_problem_NN.py`, one file per problem; submission briefs are in `docs/submissions/`.

| # | Problem | Brief |
| - | ------- | ----- |
| 1 | Smooth manufactured solution, pure Neumann | `0001-problem-1.md` |
| 2 | Piecewise-constant κ, 1D slab | `0003-problem-2.md` |
| 3 | Radially symmetric disk in disk | `0006-problem-3.md` |
| 4 | Two well-separated disks in a common annulus | `0008-problem-4.md` |
| 5 | Reentrant corner (L-shape), inverted-rate assertion | `0007-problem-5.md` |

```bash
# Full suite (~23 s with warm cache; runs in parallel by default via pytest-xdist).
uv run pytest tests/

# Just one problem (Problem 4 dominates wall time; skip it for fast TDD).
uv run pytest tests/test_problem_01.py tests/test_problem_02.py \
    tests/test_problem_03.py tests/test_problem_05.py

# Stream the per-level errors and fitted rates; serial run so logs don't interleave.
uv run pytest tests/test_problem_03.py -s -n0 --log-cli-level=INFO

# Failure-artifact emitter regression.
uv run pytest tests/test_artifacts.py
```

Failure-only diagnostics land in `tests/_artifacts/<test_id>/`. Mesh artifacts cache in `tests/_mesh_cache/`. Both are gitignored.

## Comprehension — looking at your work

The inspector is the human-facing tool. Tests catch what they assert on; the
inspector renders the things that bugs hide in.

```bash
# 6-panel dashboard (mesh + BC tags + pin, source Q, kappa, exact T,
# computed T_h, signed error) + convergence plot.
uv run python -m scripts.inspect problem_01

# Single mesh size, skip the refinement study.
uv run python -m scripts.inspect problem_01 --mesh-size 0.05 --no-convergence
```

Outputs land in `artifacts/inspect/<problem_name>/dashboard.png` and `convergence.png` (gitignored, regeneratable).

## Solver invariant checks — `diagnose`

Layer 2: structural sanity on the *solver*, independent of the exact solution. Catches "right rate, wrong constant" bugs that pass the rate test (`verification.md` § Problem 2 failure diagnostic).

```bash
# Three checks: subdomain integrity, source verification, K symmetry + residual.
uv run python -m scripts.diagnose problem_01

# At a specific mesh size (default: finest from problem.mesh_sizes()).
uv run python -m scripts.diagnose problem_01 --mesh-size 0.1
```

Writes `artifacts/inspect/<problem_name>/internals.md`. Healthy = all three sections ✓. Any ✗ is a halt-and-investigate signal. Read `docs/inspector.md` § Layer 2 before extending.
