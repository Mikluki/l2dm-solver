# Architecture

This document describes how the code is organized: the stack, the module layout, the key abstractions, the decisions already made, and what is deliberately out of scope. It is the cross-reference for `verification.md` and the contract the agent works against.

The physics lives in `physics.md`. The test problems live in `verification.md`. Coding conventions live in `CLAUDE.md`. This document covers everything else.

## Stack

- **scikit-fem** — finite element assembly, basis functions, DOF management. Pure Python on NumPy/SciPy, no compiled dependencies, no install friction. See ADR-0001.
- **gmsh** — mesh generation via the Python API. Subdomains and boundaries as physical groups, propagated through to scikit-fem.
- **meshio** — `.msh` ↔ scikit-fem mesh bridge.
- **pytest** — verification harness runner.
- **numpy, scipy** — pulled in by the above; direct use only for sparse linear algebra and array operations.
- **matplotlib** — diagnostic artifact plotting.

Versions are pinned in `pyproject.toml`. No additional runtime dependencies without an ADR.

## Pipeline

The data flow for any verification problem run:

1. **Geometry build** — `geometry/` builders construct a gmsh model with physical groups for subdomains and boundaries. Cached on disk by content hash; see ADR-0007.
2. **Mesh** — gmsh generates the `.msh` file at the requested characteristic size.
3. **Load** — `meshio.read` → `skfem.MeshTri.load`, with subdomain and boundary tags preserved as `mesh.subdomains` and `mesh.boundaries`.
4. **Basis** — `Basis(mesh, ElementTriP1())` for the trial/test space; a derived P0 basis via `basis.with_element(ElementTriP0())` for piecewise-constant fields ($\kappa$, source).
5. **Coefficients** — $\kappa$ built as a P0 array indexed by subdomain tag; source $Q$ either as a P0 array (piecewise-constant case) or evaluated at quadrature points inside the form (smooth case).
6. **Assemble** — bilinear form for the stiffness, linear form for the load. Coefficients passed via the `w` dictionary; the form uses `w.kappa` and `w.source` indirection.
7. **BCs and nullspace** — Dirichlet DOFs identified from boundary tags. For pure-Neumann problems, one node is pinned (ADR-0005).
8. **Solve** — `solve(*condense(A, b, D=dirichlet_dofs))`. Direct sparse solver via SciPy (ADR-0006).
9. **Post-process** — error norms against the exact solution, diagnostic artifacts to `tests/_artifacts/{test_id}/`.

The pipeline is the same for every problem. New problems add a Problem definition (see below); they do not modify the pipeline.

## Key abstractions

Four, deliberately. Each owns one concern and exposes a narrow interface.

### `Problem`

Defined by the contract in `verification.md` § Problem definition interface. A `Problem` is a data class (or dataclass) carrying:

- `geometry()` — gmsh-model builder, parametric in target mesh size.
- `kappa(tag)` — conductivity per subdomain tag.
- `source(x, y)` — vectorized $Q$ evaluator.
- `boundary_conditions()` — tag → BC spec dict.
- `exact_solution(x, y)` — vectorized $T$ evaluator.
- `mesh_sizes()` — list of characteristic sizes for the refinement study.
- `expected_rate()` — expected $L^2$ convergence rate.

A `Problem` is **pure data and pure functions**. It does not assemble, solve, or know about scikit-fem. This boundary is what lets Part 2 problems plug into the same harness — a Part 2 problem just returns a different operator builder instead of a different bilinear form, and the harness is none the wiser.

One Problem per file in `problems/`.

### `Solver`

Stateless functions consuming a `Problem` and a target mesh size, returning a discrete solution as a numpy array of nodal DOF values, plus the `Basis` used.

```python
def solve_scalar(problem: Problem, mesh_size: float) -> SolverResult:
    ...
```

`SolverResult` is a small dataclass: `solution`, `basis`, `mesh`. No methods.

The solver owns: geometry-to-mesh build, basis construction, coefficient assembly into P0 fields, bilinear/linear form assembly, BC application, nullspace pinning, linear solve. It does **not** own: error computation, plotting, refinement-study orchestration.

### `Harness`

The pytest-facing layer. Owns the refinement-study driver, error norms (`l2_error`, `h1_error`), rate fitting, and artifact emission.

```python
def run_refinement_study(problem: Problem) -> StudyResult:
    ...
```

`StudyResult` carries per-level errors, fitted rates, and paths to emitted artifacts. The harness does not interpret `StudyResult` — assertions live in test functions.

The harness is the part that **does not change between Part 1 and Part 2**. It must not know about scikit-fem specifics beyond the `SolverResult` shape; the integral-form solver in Part 2 returns the same shape and plugs in unchanged.

### `Geometry builders`

Helper functions in `geometry/` for the common shapes used across verification problems:

- `unit_square(mesh_size)` — Problem 1.
- `rectangle_split(mesh_size)` — Problem 2.
- `disk_in_disk(R_inner, R_outer, mesh_size)` — Problem 3.
- `two_disks_in_disk(...)` — Problem 4.
- `l_shape(mesh_size)` — Problem 5.

Each returns a gmsh model with documented physical group tags. Builders are cached on disk; see ADR-0007.

## Module structure

```
src/
  solver/          — assembly, BC, solve. Pure functions of (Problem, mesh_size).
  geometry/        — gmsh model builders. One function per shape.
  problems/        — Problem instances. One file per verification problem.
  harness/         — refinement-study driver, error norms, artifact emitters.
tests/
  test_*.py        — pytest test functions wiring Problems through the harness.
  _artifacts/      — diagnostic outputs. Gitignored.
  _mesh_cache/     — cached `.msh` files keyed by geometry hash. Gitignored.
docs/
  physics.md
  verification.md
  architecture.md
  open-questions.md
  decisions/       — ADRs.
  submissions/     — per-submission briefs.
CLAUDE.md
pyproject.toml
```

Five source modules. If any grows past one file's worth of code, that's its own decision recorded as an ADR.

## Decisions

Each is summarized here in one line; ADR files in `docs/decisions/` carry the longer rationale.

- **ADR-0001:** Stack is scikit-fem + gmsh + meshio + pytest. Pure-Python, no compiled deps, suits agentic dev introspection.
- **ADR-0002:** P1 Lagrange elements throughout Part 1. P2 deferred to Part 2 if singular-quadrature accuracy demands it.
- **ADR-0003:** All material interfaces are mesh-aligned (gmsh feature edges enforce alignment).
- **ADR-0004:** $\kappa$ represented as a P0 field, indexed per element from the subdomain tag.
- **ADR-0005:** Neumann nullspace handled by **node-pinning** — one DOF clamped to the exact value (or zero) before solve. Loud failure mode (singular matrix if forgotten) preferred over silent mean-offset.
- **ADR-0006:** Direct sparse solver (`scipy.sparse.linalg.spsolve`) for Part 1. Iterative solvers deferred to Part 2 where dense kernel matrices may require them.
- **ADR-0007:** Mesh artifacts cached on disk in `tests/_mesh_cache/` keyed by SHA-256 of geometry parameters. Re-meshing on refinement studies is the dominant cost otherwise.
- **ADR-0008:** Diagnostic artifacts emitted to `tests/_artifacts/{test_id}/` on test failure. Gitignored.

## Coefficient handling: the one example worth inlining

The form-with-coefficient pattern is unusual enough that one example is worth carrying here. Coefficients enter the form via the `w` dictionary, populated by `interpolate()` on a P0 basis:

```python
from skfem import BilinearForm, LinearForm, Basis, ElementTriP0, ElementTriP1
from skfem.helpers import grad, dot

basis = Basis(mesh, ElementTriP1())
basis_p0 = basis.with_element(ElementTriP0())

# kappa as a P0 field, one value per element
kappa_values = np.array([problem.kappa(tag) for tag in mesh.subdomains_per_element])

@BilinearForm
def stiffness(u, v, w):
    return w.kappa * dot(grad(u), grad(v))

A = stiffness.assemble(basis, kappa=basis_p0.interpolate(kappa_values))
```

The `w.kappa` indirection inside the form is what makes coefficients piecewise-constant per element. Forgetting `basis_p0.interpolate(...)` and passing the raw array silently broadcasts; the symptom is the "right rate, wrong constant" bug Problem 2 catches.

This is the only inlined example in this doc. All other API specifics live in code.

## Part 2 (out of scope, sketched)

When Part 2 begins, the integral form requires:

- Dense matrix assembly with a translation-invariant kernel $g_T(\mathbf{r}-\mathbf{r}')$.
- Singular quadrature for the $|q|^{-1}$ short-range behavior of the kernel.
- FFT-accelerated matvecs for the periodic metasurface case; iterative solver (GMRES) for the matvec-only regime.
- A `Solver`-shaped function with the same `SolverResult` output, so the harness is unchanged.

This is a placeholder, not a design. It exists so the agent knows the harness contract has to survive the transition.

## Out of scope (with motives)

The following are deliberately excluded from Part 1. Each entry says **why excluded** and **when to revisit**.

### No parallelism (single-threaded)

**Why:** Verification problems are small — finest mesh ~10⁴–10⁵ DOFs, seconds per run. Parallelism (scikit-fem's `nthreads`, MPI, etc.) adds complexity for no measured speedup at this scale.

**When to revisit:** When a single problem run exceeds ~30 s wall-clock and profiling identifies assembly or solve as the bottleneck. Then measure first, parallelize second.

### No adaptive mesh refinement

**Why:** The harness measures convergence by halving uniform mesh size and watching error decay. That only works when "mesh size" is a single number. AMR would make the rate measurement ambiguous and would actively defeat Problem 5, whose purpose is to confirm we measure the degraded corner-singularity rate honestly.

**When to revisit:** Possibly in Part 2 if singular-quadrature error concentrates around patch edges. Uniform refinement first; AMR only with explicit justification.

### No time-dependent problems

**Why:** The physics is steady-state — the substrate paper derives $-\nabla\cdot(\kappa\nabla T) = Q$, not the heat equation with $\partial_t T$. Adding time would require a time-stepping loop, stability analysis, and a parallel verification track for transient solutions — none of which serves the goal of computing $\langle\kappa\rangle$.

**When to revisit:** Only if the experimental question shifts to transient detector response. Different project.

### No GPU

**Why:** Same as parallelism, more painful. scikit-fem is CPU-only; GPU acceleration would mean either porting to a different library (JAX/CuPy-based FEM) or hand-rolling kernels. Problem sizes don't justify it. The agentic-dev tax is real — GPU failure modes (synchronization, memory placement) are much harder to debug in an agent loop.

**When to revisit:** When problem size exceeds ~10⁶ DOFs *and* profiling shows dense linear algebra as the bottleneck. Likely never in this project.

### No alternative element orders mid-Part-1

**Why:** P1 elements deliver $L^2$ rate 2, which is what all required verification problems expect. Switching to P2 mid-Part-1 would:

1. Change expected convergence rates, invalidating acceptance thresholds.
2. Complicate the P0-coefficient pattern, since P2 has edge DOFs.
3. Tempt the agent to bump element order to "fix" a failing convergence test that's actually a setup bug.

The third point is the load-bearing one. Locking element order during Part 1 closes off the worst escape route from a real verification failure.

**When to revisit:** When Part 2's integral-form kernel accuracy demonstrably plateaus with P1 — i.e., a measured need, not an aesthetic preference.

### No integral form yet (Part 2)

**Why:** Part 1's deliverable is the harness plus the scalar PDE solver. The integral form has fundamentally different machinery — dense matrices, singular kernels, no infinity BCs — that does not carry over from the scalar setting. Building both at once would either bake scalar-PDE assumptions into the harness, or force premature generalization before either operator is well-understood.

**When to revisit:** When all required (R) verification problems pass and harness artifacts are demonstrably useful. The "Future verification problems" stub in `verification.md` becomes real problems, and Part 2 begins.

---

Common theme: most exclusions exist to protect the harness from premature complexity. The harness is the deliverable. Anything that makes it harder to write, test, or trust does not earn its way in until measured pressure demands it.
