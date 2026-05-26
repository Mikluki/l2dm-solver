# thermal

Scalar steady-state heat-equation solver + verification harness. See `docs/` for the source of truth — `physics.md`, `verification.md`, `architecture.md`, `CLAUDE.md`.

## Problem 1 — what runs

Smooth manufactured solution on the unit square, pure Neumann, P1 elements. Submission brief: `docs/submissions/0001-problem-1.md`.

```bash
# Full suite (4 tests, ~0.7 s with warm cache).
uv run pytest tests/ -v

# Just the convergence study (L^2 rate ≈ 2, H^1 rate ≈ 1).
uv run pytest tests/test_problem_01.py::test_problem_01_converges -v

# Solver smoke: one solve at the coarsest mesh size, shape checks only.
uv run pytest tests/test_problem_01.py::test_problem_01_solver_smoke -v

# Mesh cache determinism: re-running must not regenerate any .msh files.
uv run pytest tests/test_problem_01.py::test_mesh_cache_is_reused -v

# Failure-artifact emitter regression: csv + rates.txt + error_field.png land
# in the per-test directory whenever the harness raises.
uv run pytest tests/test_artifacts.py -v

# Stream the per-level errors and fitted rates.
uv run pytest tests/test_problem_01.py::test_problem_01_converges -s --log-cli-level=INFO
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
