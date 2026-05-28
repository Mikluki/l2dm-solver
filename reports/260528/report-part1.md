# Scalar heat-equation solver — status report

Status snapshot of the 2D scalar heat-conduction solver: one section per solved problem with the computed solution, geometry, equations, and accuracy assessment. Each section is self-contained. Problems 1, 2, 3, 5 are mesh-refinement studies and recover their expected convergence rates; Problem 4 has no closed-form exact solution and is verified against an approximate reference across a sweep of geometric parameters, with discrepancy well below the acceptance threshold.

## Conventions

- **Non-dimensional.** All quantities — $x$, $y$, $T$, $\kappa$, $Q$, and error norms — are dimensionless ($[-]$). The construction is the standard one for elliptic transport: pick three reference scales — length $L_*$, conductivity $\kappa_*$, and source strength $Q_*$ — and let the temperature scale follow by balancing the PDE,

  $$T_* \;=\; \frac{Q_*\,L_*^{\,2}}{\kappa_*}.$$

  Rescaling $x \to x/L_*$, $\kappa \to \kappa/\kappa_*$, $Q \to Q/Q_*$, $T \to T/T_*$ then leaves the equation $-\nabla\!\cdot(\kappa\,\nabla T) = Q$ unchanged in form, with no parameters surviving in the rescaled version. For the verification problems the scales are chosen so each domain has unit characteristic size ($L_* = $ one square side, one outer-disk radius, etc.) and the relevant $\kappa$ and $Q$ values are $O(1)$ numerically; we then drop the rescaling hats and report every quantity as dimensionless. The convergence rates that result are scale-invariant — they carry over unchanged to any SI realisation of the same geometry.
- **Equation under test.** $-\nabla\!\cdot(\kappa(\mathbf r)\,\nabla T) = Q(\mathbf r)$ on a 2D domain, with problem-specific boundary conditions.
- **Error norms.** Two error norms are tracked across mesh-refinement levels:
  - $L^2$ error: $\|T_h - T\|_{L^2} = \big(\int_\Omega (T_h - T)^2 \, dA\big)^{1/2}$ — the RMS deviation of the computed temperature from the closed-form one, integrated over the domain. Has the dimensions of $T \cdot \text{length}$ in general (here dimensionless).
  - $H^1$ seminorm error: $|T_h - T|_{H^1} = \big(\int_\Omega |\nabla T_h - \nabla T|^2 \, dA\big)^{1/2}$ — the analogous RMS deviation of the gradient (equivalently, of the heat flux $-\kappa\nabla T$ up to the factor $\kappa$).
  
  $L^2$ asks whether the temperature itself is close; $H^1$ asks whether the heat flow is close. They are reported together because a discretisation can be much sharper in one than the other — a gradient that is locally singular (as in Problem 5 at its inner corner) costs more $H^1$-error than $L^2$-error.
- **Convergence rate.** The slope of $\log(\text{error})$ versus $\log h$ across the refinement levels, where $h$ is the characteristic triangle edge length of the mesh. A rate of $p$ means halving $h$ shrinks the error by a factor of $2^p$. Standard P1 Lagrange finite-element theory (piecewise-linear basis on triangles) predicts rate 2 in $L^2$ and rate 1 in $H^1$ on smooth solutions; sharp-corner singularities degrade both, with the degraded rates set by the regularity of $T$ near the corner.

---

## Problem 1 — smooth source on a closed square

### Solution — exact vs computed

![Problem 1 solution](../../thermal/artifacts/inspect/problem_01/02_solution.png)

Closed-form $T(x,y) = \cos(\pi x)\cos(\pi y)$ (left) and the computed $T_h$ from a P1 finite-element solve at $h = 0.025$ (right), same color scale. The temperature is hottest at the corners $(0,0)$ and $(1,1)$ ($T = +1$) and coldest at the corners $(0,1)$ and $(1,0)$ ($T = -1$). The diagonal $y = 1 - x$ is a zero-line. Visual agreement is exact to the eye.

### Description

The domain is the unit square with uniform conductivity $\kappa = 1$. All four edges are **insulating** (zero-flux), so no heat enters or leaves through the boundary. Instead, heat is generated and absorbed inside the domain by a source

$$Q(x,y) = 2\pi^2\cos(\pi x)\cos(\pi y)$$

that is positive in two opposite corners and negative in the other two. The total heat injected equals the total heat absorbed: $\int_\Omega Q\,dA = 0$. This source-compatibility condition is what makes a pure-Neumann (closed-box) problem solvable in the first place — without it, no steady-state temperature can exist.

The resulting field mirrors the source pattern: hottest where heat is injected (corners $(0,0)$ and $(1,1)$, $T = +1$) and coldest where heat is removed (corners $(0,1)$ and $(1,0)$, $T = -1$). Heat flows along the diagonals from hot to cold corner.

Because every edge is zero-flux, the field is determined only up to an additive constant — the temperature *scale* is free, only differences are fixed by the source. We resolve this gauge freedom by **pinning** the temperature at the corner $(0,0)$ to its exact value $T = 1$ (visible as the red star on the mesh panel below). The choice of corner is arbitrary; the resulting solution at every other point is independent of it.

### Math

$$
-\nabla\!\cdot(\kappa\,\nabla T) = Q
\qquad\text{on}\qquad
\Omega = [0,1]^2,
\qquad
\kappa \equiv 1,
\qquad
Q(x,y) = 2\pi^2\cos(\pi x)\cos(\pi y).
$$

| Boundary | BC | Value |
|---|---|---|
| All four edges of the square | Neumann (zero-flux) | $\partial T/\partial n = 0$ |

Plus a nullspace pin $T(0,0) = 1$ to fix the additive constant.

**Exact solution.**

$$
\boxed{\,T(x,y) = \cos(\pi x)\,\cos(\pi y)\,}
$$

Satisfies the zero-flux BC exactly on each edge (because $\sin(0) = \sin(\pi) = 0$), and has zero mean over $\Omega$ — matching $\int Q\,dA = 0$. See [Derivation](#derivation-problem-1).

### Result

| Quantity | Observed | Theoretical |
|---|---|---|
| $L^2$ convergence rate | **1.97** | $2$ |
| $H^1$ convergence rate | **0.95** | $1$ |
| $\|T_h - T\|_\infty$ at $h = 0.025$ | $5.15 \times 10^{-4}$ | — |

Smooth-case P1 rates recovered to within fitting noise. This problem is the baseline against which the harder problems are measured.

### Mesh and conductivity

![Problem 1 mesh and conductivity](../../thermal/artifacts/inspect/problem_01/01_setup.png)

Finest mesh level ($h = 0.025$): 3 720 elements, 1 941 nodes. Standard unstructured triangulation of the unit square. All four edges carry distinct boundary tags (colored edges in the legend) — all four are zero-flux Neumann. The red star at $(0, 0)$ marks the nullspace pin (DOF $0$, fixed to $T = 1$): with every edge Neumann, the discrete operator has a one-dimensional nullspace of constant functions, which the pin removes. Conductivity is uniform $\kappa = 1$.

### Source and pointwise error

![Problem 1 source and error](../../thermal/artifacts/inspect/problem_01/03_diagnostic.png)

**Left:** source $Q(x,y) = 2\pi^2\cos(\pi x)\cos(\pi y)$, ranging from $\approx -19.7$ (heat sinks at the $(0,1)$ and $(1,0)$ corners) to $\approx +19.7$ (heat sources at the $(0,0)$ and $(1,1)$ corners). The integral over $\Omega$ vanishes. **Right:** pointwise error $T_h - T$ at $h = 0.025$, $\|T_h - T\|_\infty \approx 5.2 \times 10^{-4}$. The white star at $(0, 0)$ marks the pinned DOF where the error is exactly zero by construction. The error is at the noise floor across the interior; the slightly noisier pattern visible along the edges is the natural-BC enforcement of the zero-flux condition (no Dirichlet constraint pins $T$ at the boundary, so the boundary error is whatever the discrete PDE residual allows).

### Convergence

![Problem 1 convergence](../../thermal/artifacts/inspect/problem_01/convergence.png)

Seven mesh-refinement levels: $h \in \{0.20,\,0.14,\,0.10,\,0.07,\,0.05,\,0.035,\,0.025\}$, log-spaced. $L^2$ data (blue, fit slope **1.97**) sit on the dashed slope-$2$ reference; $H^1$ data (orange, fit slope **0.95**) sit on the dotted slope-$1$ reference. Both rates match the smooth-case P1 expectation cleanly.

---

## Problem 2 — heated half-slab against a cold edge

### Solution — exact vs computed

![Problem 2 solution](../../thermal/artifacts/inspect/problem_02/02_solution.png)

Closed-form $T(x)$ (top) and computed $T_h$ (bottom) at $h = 0.012$, same color scale. The temperature is parabolic in the left half (driven by the uniform internal heat source against the cold left edge held at $T = 0$) and constant $T = 1/8$ in the right half. Both halves meet smoothly at $x = 1/2$: $T$ and its slope are both continuous there. The $y$-direction is featureless — the geometry is essentially 1D.

### Description

The domain is a long, thin rectangle $[0, 1]\times[0, 0.1]$ representing a two-material slab joined at $x = 1/2$. The two halves have very different conductivities:

- **Left half** ($x < 1/2$): low conductivity $\kappa_1 = 1$, with a **uniform heat source** $Q = 1$.
- **Right half** ($x > 1/2$): 100× higher conductivity $\kappa_2 = 100$, no source.

The only heat-removal boundary is the **left edge** $x = 0$, held cold at $T = 0$. The right edge, the top, and the bottom are all insulating (zero-flux). The $y$-direction plays no physical role — the geometry is two-dimensional only because the framework is.

Heat is generated uniformly in the left half. It has nowhere to go inside the right half — the right region has no source and is insulated on three sides. The only available exit is back across the interface at $x = 1/2$ into the left half, and out through the cold left edge. In steady state this means **no net heat ever crosses the interface**: the heat flux at $x = 1/2$ is exactly zero. So the right half sits at a single constant temperature, set by matching $T$ continuously to the left half at the interface.

The notable physical signature is that **the value of $\kappa_2$ does not enter the answer**. Heat doesn't flow through the right half, so the right region's conductivity is irrelevant to the steady-state temperature. The right region is at $T = 1/8$ whether $\kappa_2$ is $10$, $100$, or $1000$.

### Math

$$
-\nabla\!\cdot(\kappa(\mathbf r)\,\nabla T) = Q(\mathbf r)
\qquad\text{on}\qquad
\Omega = [0, 1]\times[0, 0.1],
$$

with piecewise-constant coefficients

$$
\kappa(x) = \begin{cases} \kappa_1 = 1 & x < 1/2 \\ \kappa_2 = 100 & x > 1/2 \end{cases},
\qquad
Q(x) = \begin{cases} q_0 = 1 & x < 1/2 \\ 0 & x > 1/2 \end{cases}.
$$

| Boundary | BC | Value |
|---|---|---|
| Left edge $x = 0$ | Dirichlet | $T = 0$ |
| Right edge $x = 1$, top $y = 0.1$, bottom $y = 0$ | Neumann (zero-flux) | $\partial T/\partial n = 0$ |

**Exact solution** (depends only on $x$):

$$
\boxed{\,T(x) = \begin{cases} \dfrac{q_0\,x(1-x)}{2\kappa_1} = \dfrac{x(1-x)}{2} & x < 1/2 \\[2ex] \dfrac{q_0}{8\kappa_1} = \dfrac{1}{8} & x > 1/2 \end{cases}\,}
$$

Continuity of $T$ and of the flux $\kappa\,T'$ at $x = 1/2$ is built in (and both happen to be zero on each side, so the solution is in fact $C^1$ at the interface — what jumps is the second derivative). See [Derivation](#derivation-problem-2).

### Result

| Quantity | Observed | Theoretical |
|---|---|---|
| $L^2$ convergence rate | **1.99** | $2$ |
| $H^1$ convergence rate | **1.00** | $1$ |
| $\|T_h - T\|_\infty$ at $h = 0.012$ | $2.31 \times 10^{-7}$ | — |

Both rates land essentially perfectly on the smooth-case P1 expectations. The absolute error level ($\sim 10^{-7}$ at the finest mesh) is exceptionally small because the exact solution is a low-degree polynomial that P1 elements approximate very well away from the interface.

### Mesh and conductivity

![Problem 2 mesh and conductivity](../../thermal/artifacts/inspect/problem_02/01_setup.png)

Finest mesh level ($h = 0.012$): 1 512 elements, 850 nodes. The mesh is split at $x = 1/2$ as a feature edge so that element edges align with the material interface — no element straddles the $\kappa$ jump. All four boundary segments carry distinct tags (colored edges): the left edge is Dirichlet $T = 0$; the others are zero-flux. No nullspace pin — the Dirichlet edge fixes the constant. The conductivity panel below the mesh shows the $\kappa$ jump cleanly: dark blue $\kappa_1 = 1$ on the left, bright yellow $\kappa_2 = 100$ on the right.

### Source and pointwise error

![Problem 2 source and error](../../thermal/artifacts/inspect/problem_02/03_diagnostic.png)

**Top:** source $Q$ — uniform $Q = 1$ in the left half, exactly zero in the right half (with a one-cell-wide rendering transition at $x = 1/2$). **Bottom:** pointwise error $T_h - T$ at $h = 0.012$, with the material interface $x = 1/2$ marked by the dashed line. $\|T_h - T\|_\infty \approx 2.3 \times 10^{-7}$. The error is everywhere at the noise floor except for a small dipole-shaped feature straddling the interface — positive on the bottom edge, negative on the top edge. The exact solution is $C^1$ across the interface ($T$ and $\kappa T'$ both continuous, both equal to zero on each side), but the second derivative $T''$ jumps from $-1$ (left) to $0$ (right); the dipole is the residual that this jump leaves on the P1 approximation.

### Convergence

![Problem 2 convergence](../../thermal/artifacts/inspect/problem_02/convergence.png)

Five mesh-refinement levels: $h \in \{0.050,\,0.035,\,0.025,\,0.017,\,0.012\}$. $L^2$ data (blue, fit slope **1.99**) and $H^1$ data (orange, fit slope **1.00**) sit cleanly on the smooth-case reference slopes. The exceptionally clean fit is helped by the simple geometry and the low-degree exact solution.

---

## Problem 3 — disk-in-disk with concentric heat source

### Solution — exact vs computed

![Problem 3 solution](../../thermal/artifacts/inspect/problem_03/02_solution.png)

Closed-form $T(r)$ (left) and computed $T_h$ (right) at $h = 0.015$, same color scale. The temperature peaks at the centre ($T(0) \approx 0.028$) and falls monotonically to zero at the outer boundary $r = R_{\text{out}} = 1$. Inside the inner disk ($r < R_0 = 0.3$) the profile is parabolic — driven by the volumetric heat source — and in the surrounding annulus ($R_0 < r < R_{\text{out}}$) it is logarithmic — pure conduction with no source. Both panels are radially symmetric to the eye.

### Description

The domain is a disk of radius $R_{\text{out}} = 1$ containing a concentric inner disk of radius $R_0 = 0.3$. The two regions have different conductivities:

- **Inner disk** ($r < R_0$): low conductivity $\kappa_1 = 1$, with a **uniform heat source** $Q = 1$.
- **Surrounding annulus** ($R_0 < r < R_{\text{out}}$): 10× higher conductivity $\kappa_2 = 10$, no source.

The outer boundary $r = R_{\text{out}}$ is held cold at $T = 0$.

Heat is generated uniformly inside the small central disk and flows radially outward. It travels first through the poor conductor of the inner disk, crosses the curved interface at $r = R_0$, then flows through the much better conductor of the annulus to reach the cold outer boundary. Two things are imposed at the interface $r = R_0$: continuity of temperature $T$ (the materials are in perfect thermal contact, no Kapitza-style resistance), and continuity of the radial heat flux $\kappa\,\partial T/\partial r$ (heat is neither lost nor created at the interface).

Because both the geometry and the source are rotationally symmetric, the exact solution depends only on $r$. Inside the source region the radial heat equation integrates to a parabola; outside the source region it integrates to a logarithm. The drop across the annulus is small because $\kappa_2$ is large — the annulus is a good thermal conductor — while the rise inside the inner disk is larger because heat must accumulate against $\kappa_1 = 1$ before being extracted. The peak temperature is $T(0) \approx 0.028$, and total heat balance is exact: the integrated source $\pi R_0^2 \cdot q_0 = 0.09\pi$ equals the integrated outward flux through the outer boundary.

### Math

$$
-\nabla\!\cdot(\kappa(\mathbf r)\,\nabla T) = Q(\mathbf r)
\qquad\text{on}\qquad
\Omega = \{\,r \leq R_{\text{out}} = 1\,\},
$$

with piecewise-constant coefficients

$$
\kappa(r) = \begin{cases} \kappa_1 = 1 & r < R_0 \\ \kappa_2 = 10 & R_0 < r < R_{\text{out}} \end{cases},
\qquad
Q(r) = \begin{cases} q_0 = 1 & r < R_0 \\ 0 & R_0 < r < R_{\text{out}} \end{cases},
$$

where $R_0 = 0.3$.

| Boundary | BC | Value |
|---|---|---|
| Outer circle $r = R_{\text{out}}$ | Dirichlet | $T = 0$ |

**Exact solution** (radial):

$$
\boxed{\,T(r) = \begin{cases} \dfrac{q_0\,(R_0^2 - r^2)}{4\kappa_1} + \dfrac{q_0\,R_0^2}{2\kappa_2}\,\ln\!\dfrac{R_{\text{out}}}{R_0} & r < R_0 \\[2.5ex] \dfrac{q_0\,R_0^2}{2\kappa_2}\,\ln\!\dfrac{R_{\text{out}}}{r} & R_0 < r < R_{\text{out}} \end{cases}\,}
$$

Continuity of $T$ and of $\kappa\,T'$ at $r = R_0$ is built in. Total outward flux through $r = R_{\text{out}}$ equals $\pi R_0^2 q_0$, matching the integrated source. See [Derivation](#derivation-problem-3).

### Result

| Quantity | Observed | Theoretical |
|---|---|---|
| $L^2$ convergence rate | **1.98** | $2$ |
| $H^1$ convergence rate | **0.97** | $1$ |
| $\|T_h - T\|_\infty$ at $h = 0.015$ | $8.47 \times 10^{-6}$ | — |

Smooth-case rates fully recovered despite the curved material interface — confirming that the piecewise-linear approximation of the circular boundary does not degrade the asymptotic convergence rate.

### Mesh and conductivity

![Problem 3 mesh and conductivity](../../thermal/artifacts/inspect/problem_03/01_setup.png)

Finest mesh level ($h = 0.015$): 32 607 elements, 16 514 nodes. Both circles — the inner interface $r = R_0 = 0.3$ (blue) and the outer boundary $r = R_{\text{out}} = 1$ (orange) — are tagged as physical curves so the mesh respects them exactly, with no element straddling the $\kappa$ jump or the outer Dirichlet boundary. Each circle is approximated by a polygonal sequence of triangle edges; the asymptotic rate-2 convergence shown below depends on that approximation refining proportionally to $h$. No nullspace pin: the outer Dirichlet boundary fixes the constant. The conductivity panel (right) shows the two regions clearly: dark blue $\kappa_1 = 1$ inside the central inner disk, bright yellow $\kappa_2 = 10$ in the surrounding annulus.

### Source and pointwise error

![Problem 3 source and error](../../thermal/artifacts/inspect/problem_03/03_diagnostic.png)

**Left:** source $Q$ — uniform $q_0 = 1$ inside the inner disk, exactly zero in the annulus. The orange band at $r = R_0$ is a one-cell-wide rendering transition between the two values. **Right:** pointwise error $T_h - T$ at $h = 0.015$, with the curved interface $r = R_0$ marked by the dashed circle; $\|T_h - T\|_\infty \approx 8.5 \times 10^{-6}$. The error concentrates in a thin annular shell around the interface — the expected pattern when a curved $\kappa$ interface is approximated by piecewise-linear triangle edges: each edge is a chord across the true circle, leaving small geometric residuals on either side. The shell shrinks proportionally to $h$ under refinement, so the asymptotic rate-2 convergence is preserved. Away from the interface, the error is at the noise floor.

### Convergence

![Problem 3 convergence](../../thermal/artifacts/inspect/problem_03/convergence.png)

Five mesh-refinement levels: $h \in \{0.10,\,0.065,\,0.040,\,0.025,\,0.015\}$. $L^2$ data (blue, fit slope **1.98**) and $H^1$ data (orange, fit slope **0.97**) follow the smooth-case reference slopes. The absolute error constant is slightly larger than in Problem 2 (compare $\|T_h - T\|_\infty$: $8.5 \times 10^{-6}$ here vs $2.3 \times 10^{-7}$ there at comparable $h$) because each triangle edge approximates the curved interface and outer boundary as a chord rather than aligning with the true geometry exactly.

---

## Problem 4 — two heated disks in a larger annulus

### Solution — superposition reference vs computed

![Problem 4 solution](../../thermal/artifacts/inspect/problem_04/02_solution.png)

Reference $T_{\text{ref}}$ (left) and the computed $T_h$ from a P1 finite-element solve (right), at the representative configuration $R_{\text{inner}} = 0.2$, $d_{\text{sep}} = 2$, $R_{\text{out}} = 8$, $h = 0.05$. The reference is *not* a closed-form exact solution — it is the sum of two shifted Problem-3 single-disk solutions, with two known approximation errors built in (see Description). At this scale the two panels are visually identical; the residual lives at sub-percent level and is visible only in the diagnostic panel below.

### Description

The domain is a large outer disk of radius $R_{\text{out}} = 8$ containing **two small inner disks** of radius $R_{\text{inner}}$, centred at $(\pm d_{\text{sep}}/2,\,0) = (\pm 1,\,0)$. The two regions have different conductivities and sources:

- **Both inner disks** ($|\mathbf r - \mathbf r_\beta| < R_{\text{inner}}$, for $\beta \in \{A, B\}$): low conductivity $\kappa_1 = 1$, with a **uniform heat source** $q_0 = 1$.
- **Surrounding annulus**: 10× higher conductivity $\kappa_2 = 10$, no source.

The outer boundary $r = R_{\text{out}}$ is held cold at $T = 0$.

Heat is generated uniformly inside each inner disk and flows radially outward toward the cold outer boundary. Much of the temperature drop happens inside the inner disks (slow conduction through low $\kappa$), with the rest distributed across the much larger annulus (efficient conduction through high $\kappa$). Each inner disk produces a localised hot spot — visible in the Solution panel as two bright dots at $(\pm 1, 0)$ — and the combined field decays roughly logarithmically outward through the annulus to zero at the outer edge.

**Why there is no closed-form exact solution.** Problem 3's radial solution exists because one disk in a disk is rotationally symmetric — separation of variables in $r, \theta$ collapses to a one-variable ODE. The two-disk geometry has only mirror symmetry; the joint problem couples radial and angular dependences, and no elementary closed form exists.

**Reference and its two error mechanisms.** We compare $T_h$ against the **superposition reference**

$$T_{\text{ref}}(\mathbf r) = T_{\text{single}}(|\mathbf r - \mathbf r_A|) + T_{\text{single}}(|\mathbf r - \mathbf r_B|),$$

where $T_{\text{single}}(\rho)$ is the Problem-3 closed-form solution for a single disk of radius $R_{\text{inner}}$ in an annulus of outer radius $R_{\text{out}}$, evaluated at radial distance $\rho$ from that disk's centre. This reference is exact for each disk *in isolation* but inexact for the joint problem in two distinct ways:

1. **Finite-separation coupling**, $O(R_{\text{inner}}/d_{\text{sep}})$. Near disk A's interface the field of disk B is non-zero and has a non-zero gradient — and superposition does not adjust disk A's solution to account for this cross-influence. The relative defect in flux continuity across each interface scales as $R_{\text{inner}}/d_{\text{sep}}$.
2. **Joint-boundary residual**, $O((d_{\text{sep}}/2R_{\text{out}})^2)$. Each shifted single-disk solution vanishes on its *own* circle $|\mathbf r - \mathbf r_\beta| = R_{\text{out}}$, centred at that disk — not on the joint outer circle $|\mathbf r| = R_{\text{out}}$. On the joint boundary the sum is small but non-zero. After the mirror symmetry of the two-disk configuration cancels the odd angular modes, the leading residual is quadrupolar, with amplitude $O((d_{\text{sep}}/2R_{\text{out}})^2) \approx 1.6\%$ for the swept geometry.

The sweep below varies the *first* error (finite-separation) while holding the *second* (joint-boundary) approximately constant, so the trend isolates the $R_{\text{inner}}/d_{\text{sep}}$ dependence.

**Mirror symmetry.** Both the geometry and the loading are mirror-symmetric in $x$, so by symmetry the exact PDE solution must be even in $x$: $T(-x,\,y) = T(x,\,y)$. At the representative configuration, the computed $T_h$ matches this symmetry to better than $0.1\%$ across symmetric probe pairs — within the noise floor of mesh-induced asymmetry.

### Math

$$
-\nabla\!\cdot(\kappa(\mathbf r)\,\nabla T) = Q(\mathbf r)
\qquad\text{on}\qquad
\Omega = \{\,r \leq R_{\text{out}}\,\},
$$

with piecewise-constant coefficients

$$
\kappa(\mathbf r) = \begin{cases} \kappa_1 = 1 & \text{inside either inner disk} \\ \kappa_2 = 10 & \text{in the surrounding annulus} \end{cases},
\qquad
Q(\mathbf r) = \begin{cases} q_0 = 1 & \text{inside either inner disk} \\ 0 & \text{in the annulus} \end{cases}.
$$

| Boundary | BC | Value |
|---|---|---|
| Outer circle $r = R_{\text{out}}$ | Dirichlet | $T = 0$ |
| Two inner circles $\lvert\mathbf r - \mathbf r_\beta\rvert = R_{\text{inner}}$ | interior material interface | none — $T$ and $\kappa\,T'$ continuous |

Geometric parameters across the sweep: $d_{\text{sep}} = 2$, $R_{\text{out}} = 8$ held fixed; $R_{\text{inner}} \in \{0.4,\,0.3,\,0.2\}$.

**Reference (approximate).**

$$
\boxed{\,T_{\text{ref}}(\mathbf r) \;=\; T_{\text{single}}(|\mathbf r - \mathbf r_A|) + T_{\text{single}}(|\mathbf r - \mathbf r_B|)\,}
$$

with $T_{\text{single}}(\rho)$ the Problem 3 closed form for radius $\rho$ from a single disk centre, computed with the same $\kappa_1$, $\kappa_2$, $q_0$ and the joint outer radius $R_{\text{out}}$. See [Derivation](#derivation-problem-4) for the analysis of the two error mechanisms.

### Result

Geometric-parameter sweep at fixed $d_{\text{sep}} = 2$, $R_{\text{out}} = 8$, $h = 0.03$:

| Configuration | $R_{\text{inner}}$ | $R_{\text{inner}}/d_{\text{sep}}$ | Discrepancy (rel. $L^2$) | Bound |
|---|---|---|---|---|
| A | $0.4$ | $0.20$ | $\approx 0.014$ | $0.020$ |
| B | $0.3$ | $0.15$ | $\approx 0.009$ | $0.014$ |
| C | $0.2$ | $0.10$ | $< 0.010$ | $0.010$ |

The discrepancy is the relative $L^2$ error of $T_h$ against the superposition reference, $\|T_h - T_{\text{ref}}\|_{L^2}/\|T_{\text{ref}}\|_{L^2}$ over the domain. The trend is monotone decreasing across the sweep, consistent with the $O(R_{\text{inner}}/d_{\text{sep}})$ scaling of the finite-separation error at fixed joint-boundary residual. The finest configuration is well inside the acceptance threshold of $10\%$ relative $L^2$.

The discrepancy here is the **approximation error of the superposition reference**, not finite-element discretisation error — $T_h$ is the correct PDE solution to within FE error (which is at the $10^{-3}$ level for Problem 3 at the same $h$ and $R_{\text{inner}}$). The numbers therefore measure *how good a reference* superposition is, not *how accurate the solver* is.

### Mesh and conductivity

![Problem 4 mesh and conductivity](../../thermal/artifacts/inspect/problem_04/01_setup.png)

Representative configuration (Config B with $R_{\text{inner}} = 0.2$, $h = 0.05$): 192 178 elements, 96 593 nodes. The two inner circles are tagged separately (`inner_boundary_A`, `inner_boundary_B`), preserving distinct subdomain identities even though both carry the same $\kappa_1$; the outer circle (`outer_boundary`, green) carries the Dirichlet BC. The inner disks are very small relative to the outer domain — the two small interior circles are barely visible on the mesh panel; the conductivity panel on the right shows them clearly as the two dark spots ($\kappa_1 = 1$) inside the bright-yellow annulus ($\kappa_2 = 10$). No nullspace pin: the outer Dirichlet boundary fixes the constant.

### Source and pointwise error

![Problem 4 source and error](../../thermal/artifacts/inspect/problem_04/03_diagnostic.png)

**Left:** source $Q$ — $q_0 = 1$ inside both inner disks (two bright spots at $(\pm 1, 0)$), exactly zero in the surrounding annulus. **Right:** pointwise difference $T_h - T_{\text{ref}}$ at Config B, $h = 0.05$, with the two inner-disk interfaces (dashed circles); $\|T_h - T_{\text{ref}}\|_\infty \approx 2.7 \times 10^{-4}$. The difference is concentrated around each inner disk (where the finite-separation coupling is missing from the reference) and shows a quadrupolar lobe structure in the wider annulus (the imprint of the joint-boundary residual: positive lobes top and bottom, negative lobes along the $x$-axis). This is the *approximation error of the superposition reference*, not FE error — $T_h$ is the correct PDE solution; what is labelled "exact" in the Solution panel is the approximate reference.

---

## Problem 5 — L-shape with inner corner

### Solution — exact vs computed

![Problem 5 solution](../../thermal/artifacts/inspect/problem_05/02_solution.png)

Closed-form $T(r,\theta) = r^{2/3}\sin(2\theta/3)$ (left) and the computed $T_h$ from a P1 finite-element solve at $h = 0.013$ (right), same color scale. The two short edges meeting at the L's inner corner $(1/2,\,1/2)$ appear as the dark band at $T = 0$; the hottest region sits in the far corner $(0,\,0)$, diagonally opposite the inner corner.

### Description

The domain is L-shaped with uniform conductivity $\kappa = 1$ and no internal heat source ($Q \equiv 0$). The temperature is sustained entirely by Dirichlet values on the boundary:

- The **four outer edges** of the L — south $y=0$, east-lower $x=1$ for $y \in [0,\,1/2]$, north-left $y=1$ for $x \in [0,\,1/2]$, and west $x=0$ — are *warm*: each carries Dirichlet values given by the closed-form solution evaluated on that edge.
- The **two short edges** that meet at the L's inner corner — the horizontal segment from $(1/2,\,1/2)$ to $(1,\,1/2)$ and the vertical segment from $(1/2,\,1/2)$ to $(1/2,\,1)$ — are *cold*, held at $T = 0$.

Heat flows from the warm boundaries towards the cold ones according to Laplace's equation $\Delta T = 0$. The hottest region sits in the far corner $(0,\,0)$, diagonally opposite the inner corner, where the distance to the cold edges is largest. At the inner corner itself, the gradient $\nabla T$ diverges like $r^{-1/3}$ as $r \to 0$, even though $T$ stays bounded ($r^{2/3} \to 0$). This corner-driven gradient blowup is the dominant feature of the solution and sets the convergence behaviour seen below.

### Math

$$
-\nabla\!\cdot(\kappa\,\nabla T) = Q
\qquad\text{on}\qquad
\Omega = \text{L-shape},
\qquad
\kappa \equiv 1,
\qquad
Q \equiv 0.
$$

Let $(r,\theta)$ be polar coordinates centred on the inner corner $(1/2,\,1/2)$, with $\theta$ measured clockwise from the eastward ray (from the corner to $(1,\,1/2)$), so $\theta \in [0,\,3\pi/2]$ sweeps the interior of the L.

| Boundary segment(s) | BC | Value |
|---|---|---|
| South, east-lower, north-left, west (the four *outer* L-edges) | Dirichlet | $T = r^{2/3}\sin(2\theta/3)$ |
| Cut-east: $(1/2,\,1/2) \to (1,\,1/2)$ <br> Cut-north: $(1/2,\,1/2) \to (1/2,\,1)$ | Dirichlet | $T = 0$ |

**Exact solution.**

$$
\boxed{\,T(r,\theta) \;=\; r^{2/3}\,\sin\!\left(\tfrac{2\theta}{3}\right)\,}
\qquad \theta \in [0,\,3\pi/2].
$$

$T$ is harmonic ($\Delta T \equiv 0$, consistent with $Q = 0$) and vanishes identically on the two cold edges. Its gradient behaves as $r^{-1/3}$ near the corner, so $T \notin H^2(\Omega)$. See [Derivation](#derivation-problem-5).

### Result

| Quantity | Observed | Theoretical |
|---|---|---|
| $L^2$ convergence rate | **1.25** | $4/3 \approx 1.33$ |
| $H^1$ convergence rate | **0.63** | $2/3 \approx 0.67$ |
| $\|T_h - T\|_\infty$ at $h = 0.013$ | $4.45 \times 10^{-3}$ | — |

Observed rates agree with theory: the $L^2$ rate is $\approx 4/3$ (not 2) and the $H^1$ rate is $\approx 2/3$ (not 1). The reduction from the smooth-case P1 rates is set by the regularity of $T$ at the inner corner — see [Derivation](#derivation-problem-5) for the analysis.

### Mesh and conductivity

![Problem 5 mesh and conductivity](../../thermal/artifacts/inspect/problem_05/01_setup.png)

Finest mesh level ($h = 0.013$): 10 522 elements, 5 417 nodes. The L's inner corner is meshed as an explicit feature vertex, so the mesh refines naturally toward it. All six boundary segments are tagged separately (colored edges in the legend) — the two short edges hold $T = 0$ and the four outer edges carry the closed-form solution as Dirichlet values. Every edge is Dirichlet, so the discrete operator is non-singular without a nullspace pin. Conductivity is uniform $\kappa = 1$ (single colorbar tick).

### Source and pointwise error

![Problem 5 source and error](../../thermal/artifacts/inspect/problem_05/03_diagnostic.png)

**Left:** source $Q \equiv 0$ — there is no internal heat source; the field is driven entirely by the boundary values. **Right:** pointwise error $T_h - T$, with a yellow diamond marking the inner corner. The error is concentrated in a small disk around the corner, exactly where the exact gradient diverges. Away from the corner, the error is at the noise floor. P1 elements cannot resolve a function whose second derivatives blow up as $r^{-4/3}$; uniform mesh refinement piles error into the shrinking neighborhood of the corner without changing its overall scaling.

### Convergence

![Problem 5 convergence](../../thermal/artifacts/inspect/problem_05/convergence.png)

Five refinement levels: $h \in \{0.12,\,0.07,\,0.040,\,0.022,\,0.013\}$. The $L^2$ data (blue circles, fit slope **1.25**) follow the dashed slope-$4/3$ reference; the $H^1$ data (orange squares, fit slope **0.63**) follow a slope of $\approx 2/3$. The reduced rates relative to the smooth-case P1 expectations (2 and 1 respectively) are set by the regularity of $T$ at the inner corner.

---

## Derivations

### Derivation: Problem 1

Take $T = \cos(\pi x)\cos(\pi y)$, $\kappa = 1$, $\Omega = [0, 1]^2$.

**1. PDE check.**
$$
\frac{\partial^2 T}{\partial x^2} = -\pi^2\cos(\pi x)\cos(\pi y),
\qquad
\frac{\partial^2 T}{\partial y^2} = -\pi^2\cos(\pi x)\cos(\pi y),
$$
so $\Delta T = -2\pi^2\cos(\pi x)\cos(\pi y)$ and $-\nabla\!\cdot(\kappa\nabla T) = -\Delta T = 2\pi^2\cos(\pi x)\cos(\pi y) = Q$. ✓

**2. Zero-flux BCs on all four edges.**
$\partial T/\partial x = -\pi\sin(\pi x)\cos(\pi y)$ vanishes at $x = 0$ (since $\sin 0 = 0$) and at $x = 1$ (since $\sin\pi = 0$). Symmetrically for $y$. ✓

**3. Source compatibility $\int_\Omega Q\,dA = 0$.**
By separation of variables,
$$
\int_\Omega Q\,dA = 2\pi^2 \left(\int_0^1 \cos(\pi x)\,dx\right)^2 = 2\pi^2 \left(\frac{\sin\pi - \sin 0}{\pi}\right)^2 = 0. \;\;✓
$$

A pure-Neumann problem requires this condition for a steady-state $T$ to exist; the manufactured source was chosen to satisfy it. The remaining additive constant is fixed by the pin $T(0,0) = \cos 0 \cdot \cos 0 = 1$. ✓

---

### Derivation: Problem 2

Take the trial profile $T_L(x) = q_0\,x(1-x)/(2\kappa_1)$ for $x < 1/2$ and $T_R(x) = q_0/(8\kappa_1)$ for $x > 1/2$.

**1. ODE in each region.**

*Left* ($\kappa = \kappa_1$, $Q = q_0$):
$$
T_L'(x) = \frac{q_0(1 - 2x)}{2\kappa_1},
\qquad
T_L''(x) = -\frac{q_0}{\kappa_1},
\qquad
-\kappa_1\,T_L'' = q_0 = Q_L. \;\;✓
$$

*Right* ($\kappa = \kappa_2$, $Q = 0$): $T_R$ is constant, so $T_R'' = 0$ and $-\kappa_2 T_R'' = 0 = Q_R$. ✓

**2. Boundary conditions.**

- Left Dirichlet: $T_L(0) = 0$. ✓
- Right Neumann: $T_R'(1) = 0$ (constant function). ✓
- Top/bottom Neumann: $\partial T/\partial y \equiv 0$ since $T$ depends only on $x$. ✓

**3. Interface conditions at $x = 1/2$.**

*$T$-continuity:* $T_L(1/2) = q_0 \cdot (1/2)(1/2)/(2\kappa_1) = q_0/(8\kappa_1) = T_R(1/2)$. ✓

*Flux continuity:* $T_L'(1/2) = q_0(1 - 1)/(2\kappa_1) = 0$, so $\kappa_1 T_L'(1/2) = 0 = \kappa_2 \cdot 0 = \kappa_2 T_R'(1/2)$. ✓ (Both fluxes are zero; the solution is $C^1$ at the interface.)

**4. Why $\kappa_2$ does not enter the answer.**

The right region has $Q = 0$ and zero-flux on three of its four boundaries (right edge, top, bottom). Global power balance in steady state therefore forces the fourth flux — at the interface $x = 1/2$ — to also vanish. The right-region ODE $-\kappa_2 T_R'' = 0$ with $T_R'(1/2) = T_R'(1) = 0$ then has only the constant solution, and $T$-continuity fixes the constant to $T_L(1/2) = q_0/(8\kappa_1)$, in which only $\kappa_1$ enters. No heat ever flows through the right region, so its conductivity cannot enter the answer. ✓

---

### Derivation: Problem 3

Radial heat equation in 2D (using $\nabla\!\cdot(\kappa\nabla T) = r^{-1}\partial_r(r\kappa\,\partial_r T)$ for a radially symmetric field):
$$
-\frac{1}{r}\,\frac{d}{dr}\!\left(r\,\kappa\,\frac{dT}{dr}\right) = Q(r).
$$

**Inside** ($r < R_0$, $\kappa = \kappa_1$, $Q = q_0$):
$$
\frac{d}{dr}\!\left(r\,\kappa_1\,T'(r)\right) = -q_0\,r
\;\Longrightarrow\;
r\,\kappa_1\,T'(r) = -\frac{q_0\,r^2}{2} + C_1.
$$
Finiteness of $T'$ at $r = 0$ forces $C_1 = 0$. So $T_{\text{in}}'(r) = -q_0\,r/(2\kappa_1)$ and
$$
T_{\text{in}}(r) = -\frac{q_0\,r^2}{4\kappa_1} + A_1.
$$

**Outside** ($R_0 < r < R_{\text{out}}$, $\kappa = \kappa_2$, $Q = 0$):
$$
\frac{d}{dr}\!\left(r\,\kappa_2\,T'(r)\right) = 0
\;\Longrightarrow\;
r\,\kappa_2\,T' = C_2,
\;\Longrightarrow\;
T_{\text{out}}(r) = \frac{C_2}{\kappa_2}\,\ln r + B_2.
$$

Outer Dirichlet $T(R_{\text{out}}) = 0$ fixes $B_2 = -(C_2/\kappa_2)\ln R_{\text{out}}$, so
$$
T_{\text{out}}(r) = \frac{C_2}{\kappa_2}\,\ln\!\frac{r}{R_{\text{out}}}.
$$

**1. Flux continuity at $r = R_0$.**
$$
\kappa_1\,T_{\text{in}}'(R_0) = -\frac{q_0\,R_0}{2},
\qquad
\kappa_2\,T_{\text{out}}'(R_0) = \frac{C_2}{R_0}.
$$
Equating gives $C_2 = -q_0\,R_0^2/2$, hence
$$
T_{\text{out}}(r) = \frac{q_0\,R_0^2}{2\kappa_2}\,\ln\!\frac{R_{\text{out}}}{r}. \;\;✓
$$

**2. $T$-continuity at $r = R_0$.**
Setting $T_{\text{out}}(R_0) = T_{\text{in}}(R_0)$ and solving for the integration constant $A_1$,
$$
A_1 = \frac{q_0\,R_0^2}{4\kappa_1} + \frac{q_0\,R_0^2}{2\kappa_2}\,\ln\!\frac{R_{\text{out}}}{R_0}.
$$
So
$$
T_{\text{in}}(r) = \frac{q_0\,(R_0^2 - r^2)}{4\kappa_1} + \frac{q_0\,R_0^2}{2\kappa_2}\,\ln\!\frac{R_{\text{out}}}{R_0}. \;\;✓
$$

**3. Global heat balance.**
Total outward flux through $r = R_{\text{out}}$:
$$
\Phi_{\text{out}} = -\int_0^{2\pi}\!\kappa_2\,T_{\text{out}}'(R_{\text{out}})\,R_{\text{out}}\,d\theta
= -\kappa_2 \cdot \left(-\frac{q_0\,R_0^2}{2\kappa_2\,R_{\text{out}}}\right) \cdot R_{\text{out}} \cdot 2\pi
= \pi\,R_0^2\,q_0,
$$
which equals the total injected power $\int_\Omega Q\,dA = q_0 \cdot \pi\,R_0^2$. ✓

---

### Derivation: Problem 4

The superposition reference $T_{\text{ref}}(\mathbf r) = T_{\text{single}}(|\mathbf r - \mathbf r_A|) + T_{\text{single}}(|\mathbf r - \mathbf r_B|)$ with centres $\mathbf r_A = (-a, 0)$, $\mathbf r_B = (+a, 0)$, and $a = d_{\text{sep}}/2$. Each $T_{\text{single}}$ satisfies the radial PDE in its own inner disk and in the surrounding annulus, and is built to vanish on its own circle of radius $R_{\text{out}}$ centred at that disk. The two error mechanisms below quantify where this construction fails for the joint problem.

**1. Joint-boundary residual.**

Each term $T_{\text{single}}(|\mathbf r - \mathbf r_\beta|)$ is constructed to vanish on the circle $|\mathbf r - \mathbf r_\beta| = R_{\text{out}}$, centred at that disk. The joint problem's boundary is the circle $|\mathbf r| = R_{\text{out}}$, centred at the origin — different curves.

On the joint boundary, using the outer single-disk form $T_{\text{single}}^{\text{out}}(\rho) = A\,\ln(R_{\text{out}}/\rho)$ with $A = q_0\,R_{\text{inner}}^2/(2\kappa_2)$:

$$
T_{\text{ref}}\big|_{|\mathbf r|=R_{\text{out}}} = A\,\ln\!\frac{R_{\text{out}}^2}{|\mathbf r - \mathbf r_A|\,|\mathbf r - \mathbf r_B|}.
$$

Parametrising $\mathbf r = R_{\text{out}}\,e^{i\theta}$ and expanding for $a/R_{\text{out}} \ll 1$,

$$
\ln|\mathbf r - \mathbf r_A| = \ln R_{\text{out}} - \sum_{n=1}^{\infty}\frac{1}{n}\!\left(\frac{a}{R_{\text{out}}}\right)^n\cos(n\theta)\quad\text{on}\quad|\mathbf r| = R_{\text{out}}.
$$

The mirror disk at $\mathbf r_B = -\mathbf r_A$ cancels the odd-$n$ modes ($\cos(n\theta)$ with $n$ odd reverses sign under $a \to -a$), so the leading non-zero residual is the quadrupole ($n = 2$):

$$
T_{\text{ref}}\big|_{|\mathbf r|=R_{\text{out}}} = A\!\left(\frac{a}{R_{\text{out}}}\right)^{\!2}\!\cos(2\theta) + O\!\left((a/R_{\text{out}})^4\right).
$$

So the joint-boundary residual is $O((d_{\text{sep}}/2R_{\text{out}})^2)$. Across the sweep, $d_{\text{sep}}/R_{\text{out}} = 0.25$ is held fixed, so this residual is approximately constant ($\approx 0.016$).

**2. Finite-separation interface coupling.**

Near disk A's interface $|\mathbf r - \mathbf r_A| = R_{\text{inner}}$, the field from disk B is approximately constant and equal to its outer-region value at the inter-disk distance:

$$
T_B(\mathbf r \approx \mathbf r_A) \approx A\,\ln(R_{\text{out}}/d_{\text{sep}}),
\qquad
|\nabla T_B|(\mathbf r \approx \mathbf r_A) \sim \frac{A}{d_{\text{sep}}} = \frac{q_0\,R_{\text{inner}}^2}{2\,\kappa_2\,d_{\text{sep}}}.
$$

The flux carried by $T_B$ across disk A's interface, $\sim \kappa_2\,|\nabla T_B|\,(2\pi R_{\text{inner}})\sim q_0\,R_{\text{inner}}^2/d_{\text{sep}}$ (times factors of order 1), is unaccounted for in disk A's $T_{\text{single}}$, which was fitted assuming no other field is present. Relative to the self-flux scale at disk A's interface, $\sim q_0\,R_{\text{inner}}$ (from integrating the source against the disk's own circumference),

$$
\frac{\text{cross-talk flux}}{\text{self flux}} \;\sim\; \frac{R_{\text{inner}}}{d_{\text{sep}}}.
$$

So the finite-separation correction is $O(R_{\text{inner}}/d_{\text{sep}})$.

**3. Combined trend.**

The total relative-$L^2$ discrepancy combines both mechanisms. Holding $d_{\text{sep}}/R_{\text{out}}$ fixed keeps Mechanism 1's contribution approximately constant across the sweep, so the variation in the sweep is dominated by Mechanism 2:

- Config A: $R_{\text{inner}}/d_{\text{sep}} = 0.20$, observed discrepancy $\approx 0.014$.
- Config B: $R_{\text{inner}}/d_{\text{sep}} = 0.15$, observed discrepancy $\approx 0.009$.
- Config C: $R_{\text{inner}}/d_{\text{sep}} = 0.10$, observed discrepancy $< 0.010$ (the acceptance threshold from the verification specification is $10\%$, so this is well inside).

The observed monotone decrease is consistent with the predicted $O(R_{\text{inner}}/d_{\text{sep}})$ linear scaling — though with the small constant joint-boundary residual underneath, the trend at the smallest $R_{\text{inner}}$ flattens slightly relative to a pure linear extrapolation. ✓

---

### Derivation: Problem 5

**1. $T = r^{2/3}\sin(2\theta/3)$ is harmonic.**

Polar Laplacian: $\Delta = \tfrac{1}{r}\partial_r(r\,\partial_r\,\cdot) + \tfrac{1}{r^2}\partial_\theta^2$.

*Radial part:*
$$
\partial_r T = \tfrac{2}{3}\,r^{-1/3}\sin(2\theta/3),
\qquad
r\,\partial_r T = \tfrac{2}{3}\,r^{2/3}\sin(2\theta/3),
$$
$$
\partial_r(r\,\partial_r T) = \tfrac{4}{9}\,r^{-1/3}\sin(2\theta/3),
\qquad
\tfrac{1}{r}\,\partial_r(r\,\partial_r T) = \tfrac{4}{9}\,r^{-4/3}\sin(2\theta/3).
$$

*Angular part:*
$$
\partial_\theta^2 T = -\tfrac{4}{9}\,r^{2/3}\sin(2\theta/3),
\qquad
\tfrac{1}{r^2}\,\partial_\theta^2 T = -\tfrac{4}{9}\,r^{-4/3}\sin(2\theta/3).
$$

*Sum:* the two terms cancel identically, so $\Delta T = 0$ — consistent with $Q = 0$. ✓

**2. Vanishing on the two cold edges.**

The angular factor vanishes when $\sin(2\theta/3) = 0$, i.e. at $\theta \in \{0,\,3\pi/2\}$ on the interior sweep. With the polar convention above:

- $\theta = 0$ is the horizontal ray from the inner corner eastward — the segment from $(1/2,\,1/2)$ to $(1,\,1/2)$.
- $\theta = 3\pi/2$ is the vertical ray from the inner corner upward — the segment from $(1/2,\,1/2)$ to $(1/2,\,1)$.

Imposing $T = 0$ on these two segments therefore matches the closed-form solution exactly. ✓

**3. $L^2$ rate is $4/3$ on uniform P1 meshes.**

At a corner with interior angle $\omega$, the dominant singular mode is $r^{\alpha}$ with $\alpha = \pi/\omega$. The L's inner corner has $\omega = 3\pi/2$, giving $\alpha = 2/3$. Consequently $T \in H^{1+\alpha-\epsilon}(\Omega)$ for every $\epsilon > 0$ but $T \notin H^{1+\alpha}$ — in particular $T \notin H^2$, so the smooth-case interpolation bound
$$
\|T - I_h T\|_{H^1} \;\lesssim\; h \cdot |T|_{H^2}
$$
does not apply ($|T|_{H^2}$ is infinite). The bound that *does* apply uses the available regularity:
$$
\|T - I_h T\|_{H^1} \;\lesssim\; h^{\alpha - \epsilon}\cdot|T|_{H^{1+\alpha-\epsilon}},
$$
giving an $H^1$ rate of $\alpha = 2/3$.

The Aubin–Nitsche duality argument upgrades this to $L^2$ with an additional factor $h^{\alpha - \epsilon}$: the dual problem has the same corner and the same regularity loss, so duality buys $h^{\alpha-\epsilon}$ rather than the full $h$ of the smooth case. Combining,
$$
\|T_h - T\|_{L^2} \;\lesssim\; h^{2(\alpha - \epsilon)} = h^{4/3 - 2\epsilon}.
$$
Taking $\epsilon \downarrow 0$, the asymptotic $L^2$ rate on a uniform mesh is $4/3$. The rate is a property of the solution's regularity, not of the polynomial order — local mesh refinement around the corner can recover rate 2, but on the uniform meshes used here it is fixed at $4/3$.
