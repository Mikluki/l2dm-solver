# Verification

This document catalogs the test problems used to verify the scalar heat-equation solver, and specifies the test harness that runs them.

The primary deliverable of Part 1 is the **harness**, not the solver. The solver is one component the harness exercises; the harness is what gets reused — unchanged — when Part 2 introduces the integral-equation form. Every choice in this document is made with that in mind.

## Scope and conventions

**Verification vs. validation.** Verification asks whether the code solves the equations it claims to solve. Validation asks whether those equations describe reality. This document is verification only. Validation against the document's $\langle\kappa\rangle$ curves belongs elsewhere and must not be used to verify components, on pain of circular reasoning.

**Equation under test.** All problems verify the scalar PDE
$$-\nabla\cdot(\kappa(\mathbf{r})\,\nabla T) = Q(\mathbf{r})$$
on a 2D domain $\Omega$. Boundary conditions are stated per problem.

**Stack.** Solver built on scikit-fem (assembly, basis, DOFs) with gmsh for mesh generation. Verification harness in pytest. Subdomain-tagged element-aligned material interfaces throughout — all problems assume the mesh respects every $\kappa$-discontinuity.

**Convergence rates.** For P1 Lagrange elements on smooth solutions, the expected rate in the $L^2$ norm is 2 and in the $H^1$ semi-norm is 1. Degraded rates from corner singularities are stated per problem.

**What failure means here.** scikit-fem is a tested library; we are not verifying its assembly. We are verifying *our* problem setup: subdomain tagging, coefficient indirection, source-field definition, BC application, post-processing. Failure diagnostics are framed accordingly.

**Acceptance.** Each problem states a qualitative acceptance criterion (e.g., "observed $L^2$ convergence rate $\geq 1.8$"). Numerical thresholds live in the test code, not here, so they can be tuned without doc churn.

**Required vs. stretch.** Problems are tagged **R** (required for the scalar solver to be considered verified), **+** (recommended), or **S** (stretch).

---

## Harness requirements

The harness is the part of Part 1 that lives longest. It must:

1. **Run as `pytest`.** Each verification problem is one or more test functions. CI-able.
2. **Support mesh refinement studies.** Each problem can be run at multiple mesh sizes; the harness computes errors at each level, fits a convergence rate, and asserts on the rate. Mesh sizes are problem-defined, not hard-coded in the harness.
3. **Compute $L^2$ and $H^1$ error norms against analytic solutions.** Implemented once, reused across problems. Uses scikit-fem's form-based integration so the norms are consistent with the assembly.
4. **Separate problem definition from test execution.** A problem is a data structure (geometry, $\kappa$, $Q$, BCs, exact solution) and a small set of functions. The harness consumes problems uniformly. New problems are added by writing a new definition, not by editing the harness.
5. **Emit diagnostic artifacts on failure.** When a convergence test fails, the harness produces: the error-vs-$h$ table, a fitted-rate value, and a plot of the error field on the finest mesh. These go to a per-test output directory, not the terminal.
6. **Be re-usable for the integral form.** A Part 2 integral-equation problem should plug into the same harness by providing the same problem-definition interface. The harness must not bake in PDE-specific assumptions beyond "there is a 2D field $T$ with a known reference."

The harness is *not* a benchmarking tool. Timing and memory are out of scope.

---

## Problem definition interface

Each verification problem is a Python object exposing:

- `geometry()` — returns a gmsh model (or a callable that builds one given a target mesh size). Subdomains tagged by physical group; boundaries tagged by physical curve.
- `kappa(subdomain_tag)` — returns the conductivity for a given tag.
- `source(x, y)` — returns $Q(x, y)$. Vectorized.
- `boundary_conditions()` — returns a dict mapping boundary tag to BC specification (Dirichlet value or Neumann flux).
- `exact_solution(x, y)` — returns $T(x, y)$ as a callable. Vectorized.
- `mesh_sizes()` — returns the list of target element sizes for the refinement study.
- `expected_rate()` — returns the expected $L^2$ convergence rate.

This interface is the cross-document contract between `verification.md` and `architecture.md`. Changes to it are decisions, not refactors.

---

## Problem 1 — Smooth manufactured solution **R**

**Purpose:** Verify that the harness end-to-end works on the simplest possible setting. The first time this passes is the moment Part 1 starts paying off.

**Domain:** Unit square $\Omega = [0,1]^2$. One subdomain.

**Material:** $\kappa = 1$ everywhere.

**Boundary conditions:** Zero-flux on all four edges.

**Source:** $Q(x,y) = 2\pi^2\cos(\pi x)\cos(\pi y)$.

**Exact solution:** $T(x,y) = \cos(\pi x)\cos(\pi y)$. Satisfies the zero-flux BC exactly and has zero mean over $\Omega$, so $\int Q\,dA = 0$ as required by pure Neumann.

**Nullspace handling:** The pure-Neumann problem defines $T$ only up to a constant. Resolved by node-pinning per `architecture.md` § Nullspace handling — pin one DOF to the exact value at a point declared by the `Problem`. For Problem 1: `pin_point() = (0, 0)`, where $T = 1$. A geometry corner is a guaranteed mesh node at every refinement, keeping the pinned DOF index reproducible across the convergence study. Note in the test code which one is used.

**Failure diagnostic:** If this fails, the harness itself is wrong — wiring between problem definition, scikit-fem, and the error norm. Until this passes, no other problem is meaningful. Do not move on.

**Expected convergence rate:** $L^2$ rate 2, $H^1$ rate 1.

**Acceptance:** Observed $L^2$ rate $\geq 1.8$ over at least three mesh refinements; $H^1$ rate $\geq 0.9$.

**scikit-fem notes:** Single-subdomain mesh, no tagging needed. `Basis(mesh, ElementTriP1())`. Standard `solve(*condense(A, b))` after applying the nullspace fix.

---

## Problem 2 — Piecewise-constant $\kappa$, 1D slab **R**

**Purpose:** Verify that subdomain tagging, coefficient indirection, and the assembly's use of per-element $\kappa$ all line up correctly. This is the test most likely to pass at the right convergence rate but the wrong constant if any link in the chain is wrong.

**Domain:** Rectangle $\Omega = [0,1]\times[0,h]$ with $h = 0.1$. Geometry is effectively 1D; the $y$-direction is present only because the framework is 2D.

**Material:** Two subdomains split at $x = 1/2$:
$$\kappa(x,y) = \begin{cases} \kappa_1 = 1, & x < 1/2,\\ \kappa_2 = 100, & x > 1/2.\end{cases}$$
Large contrast is deliberate — it makes wrong-side-of-interface bugs loud.

**Boundary conditions:** Dirichlet $T = 0$ on the left edge ($x = 0$); zero-flux on the other three edges. The Dirichlet condition removes the nullspace cleanly.

**Source:** Uniform in the left subdomain, zero in the right:
$$Q(x,y) = \begin{cases} q_0 = 1, & x < 1/2,\\ 0, & x > 1/2.\end{cases}$$

**Exact solution:** From the 1D ODE $-(\kappa T')' = Q$ with $T(0)=0$ and $T'(1)=0$:

For $x < 1/2$: $T(x) = q_0 x (1 - x) / (2\kappa_1)$.

For $x > 1/2$: $T(x) = q_0 / (8\kappa_1) = $ constant.

Continuity of $T$ and of $\kappa T'$ at $x = 1/2$ is built in.

Key sanity check during development: the answer for $x > 1/2$ **does not depend on $\kappa_2$**. Heat that crosses the interface has nowhere to go (zero flux on the right), so the right region equilibrates to a constant. If the test result depends on $\kappa_2$, the bug is found.

**Failure diagnostic:** If Problem 1 passes and this fails, the bug is in *our* problem setup:
- Wrong subdomain tag → wrong $\kappa$ in assembly (most common).
- gmsh physical groups not propagating through meshio to scikit-fem (check `mesh.subdomains`).
- Coefficient passed as a uniform scalar instead of a P0 element-wise field.
- $Q$ field defined on wrong subdomain.

If the $L^2$ error is right at coarse meshes but degrades at fine ones, the mesh is not aligning with $x = 1/2$ — gmsh needs the line $x = 1/2$ as a feature edge.

**Expected convergence rate:** $L^2$ rate 2.

**Acceptance:** Observed rate $\geq 1.8$. Additionally, **assert that the computed $T$ at $x = 0.75$ does not depend on $\kappa_2$ to within 1% when $\kappa_2$ is varied between 10 and 1000.** This catches bugs that get the rate right but the constant wrong.

**scikit-fem notes:** Two subdomains via `mesh.with_subdomains({...})` or gmsh physical groups. $\kappa$ as a P0 (`ElementTriP0`) field, evaluated per element from subdomain tags. Standard `condense(A, b, D=dirichlet_dofs)`.

---

## Problem 3 — Radially symmetric disk in larger disk **R**

**Purpose:** Verify curved-boundary handling. Catches issues that the axis-aligned interface in Problem 2 cannot.

**Domain:** Disk of radius $R_{\text{out}} = 1$, with an inner disk of radius $R_0 = 0.3$ centered at the origin. Two subdomains.

Note: we use a disk outer boundary rather than a square because that gives a closed-form radial solution. A square outer boundary breaks radial symmetry and the exact solution becomes a reference-solution comparison instead of a closed-form one — that's a strictly worse test.

**Material:** $\kappa_1 = 1$ inside, $\kappa_2 = 10$ outside.

**Boundary conditions:** Dirichlet $T = 0$ on $r = R_{\text{out}}$.

**Source:** $Q = q_0 = 1$ inside the disk, $0$ outside.

**Exact solution:** From the 1D radial ODE $-r^{-1}(r\kappa T')' = Q$ with $T(R_{\text{out}}) = 0$:

For $r < R_0$:
$$T(r) = \frac{q_0(R_0^2 - r^2)}{4\kappa_1} + \frac{q_0 R_0^2}{2\kappa_2}\ln\frac{R_{\text{out}}}{R_0}.$$

For $R_0 < r < R_{\text{out}}$:
$$T(r) = \frac{q_0 R_0^2}{2\kappa_2}\ln\frac{R_{\text{out}}}{r}.$$

Continuity of $T$ and $\kappa T'$ at $r = R_0$ is built in; total flux through $r = R_{\text{out}}$ equals $\pi R_0^2 q_0$ as required.

**Mesh requirement:** The mesh must respect both circles ($r = R_0$ and $r = R_{\text{out}}$). In gmsh, both as physical curves with an enforced characteristic length.

**Failure diagnostic:** If Problems 1–2 pass and this fails:
- Curved boundary not being captured (check that the mesh has edges along $r = R_0$).
- Subdomain assignment by point-in-circle test failing on near-boundary elements — should be by gmsh physical surface tag, not by element-centroid coordinates.
- Error concentrated in a thin annular shell around $r = R_0$ → mesh-to-circle approximation is the bottleneck, not a bug.

**Expected convergence rate:** $L^2$ rate 2. The constant is somewhat worse than Problem 2 because the curved interface is piecewise-linearly approximated; the rate is still 2.

**Acceptance:** Observed rate $\geq 1.8$. Additionally, plot the error field at the finest mesh and confirm visually (recorded as an artifact, not asserted) that the error is not dominated by a sharp annular shell.

**scikit-fem notes:** Two subdomains from gmsh physical surfaces. The outer disk boundary as a physical curve for Dirichlet DOF lookup. P0 $\kappa$ field as in Problem 2.

---

## Problem 4 — Two well-separated disks **+**

**Purpose:** Verify that multiple disconnected subdomains compose correctly. This is the closest scalar-setting analog to the eventual multi-patch (graphene + metal contacts, or metasurface) case.

**Domain:** Larger disk of radius $R_{\text{out}}$ containing two inner disks of radius $R_0$, centered at $(\pm a, 0)$ with center-to-center separation $d_{\text{sep}} = 2a$. The acceptance family must keep both $R_0/d_{\text{sep}}$ and $d_{\text{sep}}/R_{\text{out}}$ controlled; a representative geometry is $d_{\text{sep}} = 2$, $R_0 = 0.2$, and $R_{\text{out}} \geq 8$.

**Material:** $\kappa_1 = 1$ inside both inner disks, $\kappa_2 = 10$ outside.

**Boundary conditions:** Dirichlet $T = 0$ on $r = R_{\text{out}}$.

**Source:** $Q = q_0 = 1$ inside both inner disks, $0$ outside.

**Exact solution:** No closed form. Acceptance is **superposition-based**, but the approximation has two distinct error mechanisms that must not be conflated.

The reference approximation is the sum of two single-disk solutions from Problem 3, each shifted to its disk center:
$$T_{\text{two}}(\mathbf{r}) \approx T_{\text{single}}(|\mathbf{r} - \mathbf{r}_A|) + T_{\text{single}}(|\mathbf{r} - \mathbf{r}_B|).$$

Finite-separation/interface corrections scale as $O(R_0/d_{\text{sep}})$. The reason is local: near disk A, the field from disk B has gradient scale $O(q_0R_0^2/(\kappa_2 d_{\text{sep}}))$. Applied across disk A's $\kappa_1/\kappa_2$ interface, that smooth field creates a relative flux-continuity defect of order $R_0/d_{\text{sep}}$ compared with the self flux scale $O(q_0R_0)$.

The shifted single-disk superposition also fails the **common outer Dirichlet boundary**. Each shifted single-disk solution is zero on a circle centered at its own disk, not on the joint circle $r = R_{\text{out}}$. For symmetric centers $(\pm a, 0)$, the leading boundary residual is quadrupolar and scales as $O((a/R_{\text{out}})^2) = O((d_{\text{sep}}/(2R_{\text{out}}))^2)$. Therefore the approximation becomes controlled only when both $R_0/d_{\text{sep}} \to 0$ and $d_{\text{sep}}/R_{\text{out}} \to 0$ (or when the reference approximation is changed to use a common-domain Green's function).

**Failure diagnostic:** If Problems 1–3 pass and this fails:
- gmsh failing to honor both embedded circles as separate physical surfaces.
- Subdomain tagging collapsing both inner disks into one tag (legitimate, but then $\kappa$ assignment needs to handle it).
- Assembly looping over elements but skipping one subdomain due to indexing.

**Expected behavior:** This is a *convergence-in-geometric-parameter* test, not a convergence-in-$h$ test. Fix the mesh resolution and vary $R_0/d_{\text{sep}}$ while keeping $d_{\text{sep}}/R_{\text{out}}$ small and approximately fixed; equivalently, vary $d_{\text{sep}}$ and $R_{\text{out}}$ together so the outer-boundary residual does not become the dominant changing error. Do **not** infer convergence by moving the disks farther apart inside a fixed outer disk: that decreases $R_0/d_{\text{sep}}$ but increases the boundary/image error through $d_{\text{sep}}/R_{\text{out}}$.

**Acceptance:** With $R_0/d_{\text{sep}} \leq 0.1$ and $d_{\text{sep}}/R_{\text{out}} \leq 0.25$, relative discrepancy below 10% in $L^2$. Loose threshold because the approximation has both finite-separation and boundary/image error; the controlled trend matters more than the number.

**scikit-fem notes:** Two physical surfaces in gmsh, both tagged with the same subdomain marker (since both are "inner") — or tagged distinctly if we later want different $\kappa$ in each. Verify by inspection that the mesh has both circles resolved.

---

## Problem 5 — Reentrant corner **S**

**Purpose:** Verify that the convergence-rate measurement in the harness is honest in the presence of solution singularities. If we ever claim a high convergence rate on a problem with a corner, this test confirms whether the measurement is real.

**Domain:** L-shape — unit square with the upper-right quadrant removed: $\Omega = [0,1]^2 \setminus [1/2, 1]^2$. Reentrant corner at $(1/2, 1/2)$.

**Material:** $\kappa = 1$ uniform.

**Boundary conditions:** Dirichlet, with values chosen so that the exact solution in polar coordinates centered at the reentrant corner is
$$T(r,\theta) = r^{2/3}\sin(2\theta/3),$$
$\theta$ measured from one of the boundary edges. Harmonic, singular at the corner.

**Source:** $Q = 0$.

**Exact solution:** As above; standard L-shape benchmark.

**Failure diagnostic:** Not really a "bug" failure mode. If this passes but reports an $L^2$ rate of 2, the rate-fitting code is wrong (probably masking the corner by not refining around it).

**Expected convergence rate:** $L^2$ rate $4/3$ on uniform meshes (not 2).

**Acceptance:** Observed rate in $[1.2, 1.5]$. **A rate $\geq 1.8$ here fails the test** — that's the inverted-assertion case.

**scikit-fem notes:** Single subdomain. L-shape can be built in gmsh as a polygon. Dirichlet BCs applied as a callable that evaluates the exact solution on boundary DOFs.

---

## Future verification problems

To be added when the corresponding features come into scope:

- **Kapitza coupling:** scalar PDE with a linear sink term $-G_K T$ representing leakage. Tests the bilinear-form modification.
- **Integral form, single isolated patch:** verifies the kernel evaluation and singular quadrature. Reference solution via the document's effective-conductivity formulas for a single disk on a substrate.
- **Periodic metasurface:** verifies the regularized Green's function and the zero-mode subtraction.

These define the Part 2 deliverable. Their interfaces should match the problem-definition contract above; if they don't, the contract needs revision.
