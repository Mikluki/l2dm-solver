# Submission 0002 — Algebraic verification of exact solutions

**Status:** in-progress (Problem 4 discrepancy pending adjudication; see § Verification results → Problem 4)
**Predecessors:** none. Runs independently of 0001 (no shared files, no shared decisions).
**Successors:** every subsequent verification-problem submission (Problem 2 onward) depends on the conclusions here being signed off.

## Goal

Independently re-derive each exact solution in `verification.md` Problems 1–5 and confirm it satisfies the stated PDE, boundary conditions, source compatibility, and (where applicable) interface and flux continuity — before any of those solutions get baked into a `Problem` definition.

## Why now

Every verification problem is only as trustworthy as its analytic answer. A doc-level algebra error becomes a silent test failure that looks like a solver bug: the implementer has no way to distinguish "my code is wrong" from "the doc was wrong", and the harness's whole point — measuring whether the solver does what it claims — collapses. This pass catches doc-side errors before they reach code.

The work is **pure derivation**. It touches no Python files and no `src/` or `tests/` directory. It can run fully in parallel with the 0001 implementation worker; the two never collide.

## Relevant core-doc sections

- `docs/verification.md` Problems 1–5 — exact solutions, BCs, sources, acceptance criteria to be re-derived.
- `docs/physics.md` § Equations → "What Part 1 solves" — the operator under test, $-\nabla\cdot(\kappa\nabla T) = Q$, and the symbol convention.
- `CLAUDE.md` § Doc-editing rules — `verification.md` is **human-owned**. If a discrepancy is found, raise it in `open-questions.md` and **stop**. The worker does not edit `verification.md` under any circumstance.

## Decisions resolved before this submission

These are settled; not open for the worker to revisit silently.

1. **`verification.md` is human-owned and read-only for this worker.** Per CLAUDE.md doc-editing rules. Discrepancies are logged, not patched.
2. **Discrepancy log goes to `docs/open-questions.md`.** Append, do not overwrite. The file may be empty or may already contain entries from elsewhere; either is fine. Each discrepancy is its own section with a clear title, the source-doc location (Problem N), the claim being checked, the worker's derivation, and the proposed resolution (for human adjudication).
3. **Derivations are recorded in this brief itself.** Append a `## Verification results` section at the bottom. One subsection per Problem. The brief becomes the durable artifact; no new `.md` file is created elsewhere (per CLAUDE.md doc discipline).
4. **Symbolic computation tools (SymPy, pen and paper, etc.) are allowed as scratch.** Their session output does not land in the repo. What lands is the human-readable LaTeX/math in this brief.

## Per-problem checks required

For each problem, **show the work**. "Verified" or "looks correct" without explicit derivation does not count.

### Problem 1 — smooth manufactured solution

1. Compute $-\nabla\cdot(\kappa\nabla T)$ with $\kappa=1$ and $T = \cos(\pi x)\cos(\pi y)$. Confirm it equals the stated source $Q = 2\pi^2\cos(\pi x)\cos(\pi y)$.
2. Evaluate $\partial T/\partial n$ on each of the four edges of the unit square. Confirm zero on all four (this is what makes the zero-flux BC compatible).
3. Compute $\int_\Omega Q\,dA$. Confirm zero — the pure-Neumann source-compatibility condition.
4. Compute $\int_\Omega T\,dA$. Confirm zero — required for "subtract the mean" nullspace option (a) to be self-consistent with the stated exact solution. (Note: 0001 uses node-pinning, not option (a). This check is still necessary because it verifies the doc is internally consistent regardless of implementation choice.)

### Problem 2 — piecewise-constant κ, 1D slab

1. Verify the piecewise $T$ satisfies $-(\kappa T')' = Q$ separately on each subdomain ($x<1/2$ with $\kappa_1, q_0$; $x>1/2$ with $\kappa_2, 0$).
2. Check $T(0) = 0$ (left Dirichlet) and $T'(1) = 0$ (right Neumann).
3. Check continuity of $T$ and continuity of the flux $\kappa T'$ at $x = 1/2$.
4. **Confirm explicitly that the right-region formula $T(x>1/2) = q_0/(8\kappa_1)$ has no $\kappa_2$ dependence.** This is the load-bearing structural property the acceptance check asserts; if the doc's formula secretly depends on $\kappa_2$, the κ₂-sweep acceptance test will fail at the doc level, not the code level. Show the derivation from continuity, not just substitute and stare.

### Problem 3 — radially symmetric disk in disk

1. Verify the stated $T(r)$ satisfies the radial form $-r^{-1}\,d/dr(r\kappa\,dT/dr) = Q$ in each region.
2. Check $T(R_\text{out}) = 0$.
3. Check continuity of $T$ and of the radial flux $\kappa\,dT/dr$ at $r = R_0$.
4. Compute the total outward flux at $r = R_\text{out}$: $\Phi = -\int_0^{2\pi} \kappa_2 (dT/dr)|_{R_\text{out}} R_\text{out}\,d\theta$. Confirm $\Phi = \pi R_0^2 q_0$ (total source). This is the global flux-balance check; it is independent of the spatial profile and catches sign or factor errors that point-checks miss.

### Problem 4 — two well-separated disks

No closed-form exact solution; this problem is a **convergence-in-geometric-parameter** test. The derivation work is therefore about justifying the form of the approximation, not checking a formula.

1. State why the naive superposition $T_\text{two}(\mathbf{r}) = T_\text{single}(|\mathbf{r}-\mathbf{r}_A|) + T_\text{single}(|\mathbf{r}-\mathbf{r}_B|)$ **fails to exactly satisfy** $T = 0$ on $r = R_\text{out}$, and identify the asymptotic limit in which it becomes exact ($R_0/d \to 0$ and $R_0/R_\text{out} \to 0$).
2. Justify the doc's claim that the leading correction to superposition is $O(R_0/d)$ (rather than e.g. $O((R_0/d)^2)$). An image-charge / multipole argument is fine; the goal is to verify the order, not compute the prefactor. If the worker's analysis suggests a different leading order, that is a discrepancy and goes to `open-questions.md`.

### Problem 5 — reentrant corner

1. Confirm $T = r^{2/3}\sin(2\theta/3)$ is harmonic ($\Delta T = 0$) wherever it is defined. Use the polar Laplacian; the calculation is short.
2. Identify which two edges of the L-shape this $T$ vanishes on naturally (so the Dirichlet boundary values $T|_{\partial\Omega}$ are nontrivial only on the remaining edges). This matters because the test wiring needs to evaluate the exact solution on the right subset of boundary DOFs.
3. **Derive (briefly) why uniform-mesh P1 elements give $L^2$ rate $4/3$, not 2.** Relate the corner-singularity exponent $\alpha = 2/3$ to the rate. The standard argument: $T \in H^{1+\alpha-\epsilon}$, so best P1 approximation in $L^2$ on a uniform mesh of size $h$ is $O(h^{2\alpha}) = O(h^{4/3})$. This is the test for whether the worker understands what Problem 5 is *measuring*, not just what it computes.

## Deliverable

- This brief, updated in-place with a `## Verification results` section appended at the bottom containing one subsection per problem and a short LaTeX/math derivation per check.
- If any discrepancy is found: a new entry in `docs/open-questions.md` for each one, with enough detail (Problem N, claim, derivation, proposed resolution) that a human can adjudicate without re-deriving from scratch.
- **No new files**, no edits to `verification.md`, no edits to `physics.md`, no edits to `architecture.md`, no `.py` files created or modified.

## Acceptance

All must hold simultaneously.

1. **All five problems have a recorded derivation in the `## Verification results` section of this brief.** A bare "verified" without showing the math is not acceptance.
2. **Every check listed in "Per-problem checks required" above is explicitly addressed**, either with the derivation confirming the doc, or with an entry in `open-questions.md` flagging a discrepancy. Skipping a check is not allowed.
3. **No edits to `verification.md`.** If a worker finds verification.md needs to be corrected, the discrepancy goes to `open-questions.md` and this submission's status stays `in-progress` until a human resolves it.
4. **Status transitions:**
   - If all 5 problems verify cleanly: status → `done`.
   - If any discrepancy is logged: status stays `in-progress`. The submission is not done until the human resolves the discrepancy and (if needed) `verification.md` is updated. The worker does not close the submission unilaterally on discrepancy.

## Out of scope for this submission

- Any Python code. No `Problem` objects, no harness code, no tests, no `src/` or `tests/` edits.
- Numerical verification (plugging a value into a calculator) as a *substitute* for symbolic derivation. It is fine as a sanity cross-check, but the brief records the symbolic work.
- Recommending changes to the verification-problem set itself (e.g., "Problem 4 should be replaced with..."). That is a scope decision, not a verification task. Out of scope here; raise in `open-questions.md` as a separate item if the worker has such a concern.
- Editing `verification.md` or `physics.md` (per CLAUDE.md doc-editing rules — both are human-owned).
- Editing `architecture.md` or adding ADRs.
- The "harness machinery unit tests" recommendation (l2_error returns 0 on the exact solution, rate-fitter returns 2 on synthetic $h^2$ errors). That is a separate planned submission (one option for 0003) and should not be smuggled in here.

## Done definition

This submission is done when:

1. All acceptance criteria above hold.
2. This file's status is updated from `proposed` to either `done` or `in-progress` per the acceptance rules.
3. The `## Verification results` section is present, complete, and self-contained (a human can read it without re-deriving).
4. If any discrepancies were found, the corresponding entries exist in `docs/open-questions.md` and are clearly cross-referenced from this brief.

> **Discrepancy-log location override (this submission only).** CLAUDE.md was updated mid-submission to mark `docs/open-questions.md` as planner-only and off-limits to workers. Per direct user instruction, the Problem 4 discrepancy below is recorded **in this brief only**, not appended to `open-questions.md`. The done-definition clause above is satisfied by the in-brief record; cross-reference to `open-questions.md` is intentionally absent.

---

## Verification results

Problems 1, 2, 3, 5 verify cleanly against the PDE, BCs, source compatibility, and (where applicable) interface and flux continuity stated in `verification.md`.

Problem 4 has a discrepancy in the stated asymptotic justification. The doc's $O(R_0/d)$ scaling is a real finite-separation/interface mechanism, but the documented bounded-domain superposition also has a common outer-boundary residual controlled by the disk-center offset relative to $R_{\text{out}}$. Details in the Problem 4 subsection below; this is what holds the submission at `in-progress`.

### Problem 1 — smooth manufactured solution

$T = \cos(\pi x)\cos(\pi y)$, $\kappa = 1$, $\Omega = [0,1]^2$.

**1. PDE check.** Componentwise,
$$
\frac{\partial^2 T}{\partial x^2} = -\pi^2\cos(\pi x)\cos(\pi y), \qquad \frac{\partial^2 T}{\partial y^2} = -\pi^2\cos(\pi x)\cos(\pi y),
$$
so $\Delta T = -2\pi^2\cos(\pi x)\cos(\pi y)$ and $-\nabla\cdot(\kappa\nabla T) = -\Delta T = 2\pi^2\cos(\pi x)\cos(\pi y) = Q$. ✓

**2. Zero-flux BC on each edge.**
- $x=0$: $\partial T/\partial x = -\pi\sin(0)\cos(\pi y) = 0$. ✓
- $x=1$: $\partial T/\partial x = -\pi\sin(\pi)\cos(\pi y) = 0$. ✓
- $y=0$: $\partial T/\partial y = -\pi\cos(\pi x)\sin(0) = 0$. ✓
- $y=1$: $\partial T/\partial y = -\pi\cos(\pi x)\sin(\pi) = 0$. ✓

**3. Source-compatibility $\int_\Omega Q\,dA = 0$.**
$$
\int_\Omega Q\,dA = 2\pi^2\Bigl[\textstyle\int_0^1 \cos(\pi x)\,dx\Bigr]^2 = 2\pi^2\Bigl[\frac{\sin\pi - \sin 0}{\pi}\Bigr]^2 = 0. \;\;✓
$$

**4. Zero-mean solution $\int_\Omega T\,dA = 0$.**
$$
\int_\Omega T\,dA = \Bigl[\textstyle\int_0^1 \cos(\pi x)\,dx\Bigr]^2 = 0. \;\;✓
$$
(Consistent with the doc regardless of which nullspace fix the implementation uses.)

### Problem 2 — piecewise-constant κ, 1D slab

Let $T_L(x) = q_0 x(1-x)/(2\kappa_1)$ for $x < 1/2$, $T_R(x) = q_0/(8\kappa_1)$ for $x > 1/2$.

**1. ODE in each region.**
- Left: $T_L' = q_0(1-2x)/(2\kappa_1)$, $T_L'' = -q_0/\kappa_1$, so $-\kappa_1 T_L'' = q_0 = Q_L$. ✓
- Right: $T_R$ is constant, so $T_R'' = 0$ and $-\kappa_2 T_R'' = 0 = Q_R$. ✓

**2. Boundary conditions.**
- Left Dirichlet: $T_L(0) = 0$. ✓
- Right Neumann: $T_R'(1) = 0$ (constant). ✓

**3. Interface conditions at $x = 1/2$.**
- $T$-continuity: $T_L(1/2) = q_0 \cdot (1/2) \cdot (1/2) / (2\kappa_1) = q_0/(8\kappa_1) = T_R(1/2)$. ✓
- Flux continuity: $T_L'(1/2) = q_0(1-1)/(2\kappa_1) = 0$, so $\kappa_1 T_L'(1/2) = 0 = \kappa_2 \cdot 0 = \kappa_2 T_R'(1/2)$. ✓

**4. κ₂-independence of the right region (load-bearing).** Derived from continuity, not by substitution:

The right region has $Q = 0$, Neumann $T_R'(1) = 0$, and an unknown flux at $x = 1/2$. The 1D slab framing gives three of four right-region boundaries as zero-flux (the two long edges by the problem framing, plus $x = 1$). With no sources in the right region, global power balance forces the fourth flux (at $x = 1/2$) to also vanish:
$$
\kappa_2 T_R'(1/2) = 0 \;\Longrightarrow\; T_R'(1/2) = 0.
$$
The right-region ODE $-\kappa_2 T_R'' = 0$ with $T_R'(1/2) = T_R'(1) = 0$ has only the constant solution. $T$-continuity then fixes the constant to $T_L(1/2) = q_0/(8\kappa_1)$, in which **only $\kappa_1$ enters**. ✓ This is precisely the structural property the κ₂-sweep acceptance check probes.

### Problem 3 — radially symmetric disk in disk

Radial ODE: $-\frac{1}{r}\,\frac{d}{dr}\!\left(r\kappa\,\frac{dT}{dr}\right) = Q$. Derive from scratch to check the doc's formulas.

**Inside** ($r < R_0$, $\kappa = \kappa_1$, $Q = q_0$): integrating $\frac{d}{dr}(r\kappa_1 T') = -q_0 r$ gives $r\kappa_1 T'(r) = -q_0 r^2/2 + C_1$. Finiteness at $r = 0$ forces $C_1 = 0$, so $T'(r) = -q_0 r/(2\kappa_1)$ and $T(r) = -q_0 r^2/(4\kappa_1) + A_1$.

**Outside** ($R_0 < r < R_{\text{out}}$, $\kappa = \kappa_2$, $Q = 0$): $r\kappa_2 T' = C_2$ (constant), so $T(r) = (C_2/\kappa_2)\ln r + B_2$. Outer Dirichlet $T(R_{\text{out}}) = 0$ gives $B_2 = -(C_2/\kappa_2)\ln R_{\text{out}}$, hence
$$
T_{\text{out}}(r) = \frac{C_2}{\kappa_2}\ln\!\frac{r}{R_{\text{out}}}.
$$

**1. Flux continuity at $r = R_0$.**
$$
\kappa_1 T'(R_0^-) = -\frac{q_0 R_0}{2}, \qquad \kappa_2 T'(R_0^+) = \frac{C_2}{R_0}.
$$
Equating: $C_2 = -q_0 R_0^2/2$. Substituting:
$$
T_{\text{out}}(r) = \frac{q_0 R_0^2}{2\kappa_2}\,\ln\!\frac{R_{\text{out}}}{r}. \;\;✓\;\text{matches doc}
$$

**2. $T$-continuity at $r = R_0$.**
$$
T_{\text{out}}(R_0) = \frac{q_0 R_0^2}{2\kappa_2}\ln\!\frac{R_{\text{out}}}{R_0} \;=\; T_{\text{in}}(R_0) = -\frac{q_0 R_0^2}{4\kappa_1} + A_1
$$
gives $A_1 = \frac{q_0 R_0^2}{4\kappa_1} + \frac{q_0 R_0^2}{2\kappa_2}\ln\!\frac{R_{\text{out}}}{R_0}$, so
$$
T_{\text{in}}(r) = \frac{q_0(R_0^2 - r^2)}{4\kappa_1} + \frac{q_0 R_0^2}{2\kappa_2}\ln\!\frac{R_{\text{out}}}{R_0}. \;\;✓\;\text{matches doc}
$$

**3. $T(R_{\text{out}}) = 0$.** $T_{\text{out}}(R_{\text{out}}) = (q_0 R_0^2/2\kappa_2)\ln 1 = 0$. ✓

**4. Global flux balance at $r = R_{\text{out}}$.**
$$
\frac{dT_{\text{out}}}{dr}\bigg|_{R_{\text{out}}} = \frac{q_0 R_0^2}{2\kappa_2}\cdot\left(-\frac{1}{R_{\text{out}}}\right).
$$
Total outward flux:
$$
\Phi = -\int_0^{2\pi}\kappa_2\,\frac{dT}{dr}\bigg|_{R_{\text{out}}}\!R_{\text{out}}\,d\theta = -\kappa_2 \cdot \left(-\frac{q_0 R_0^2}{2\kappa_2 R_{\text{out}}}\right)\cdot R_{\text{out}} \cdot 2\pi = \pi R_0^2 q_0. \;\;✓
$$
Equals $\int_\Omega Q\,dA = q_0 \cdot \pi R_0^2$ as required.

### Problem 4 — two well-separated disks ⚠ discrepancy

No closed-form exact solution; the verification is qualitative + asymptotic. The recheck found two distinct approximation defects that need to be separated.

**1. Why naive superposition fails the outer Dirichlet BC.**

Each shifted single-disk solution $T_{\text{single}}(|\mathbf{r}-\mathbf{r}_\beta|)$ is constructed to satisfy $T = 0$ on the circle $|\mathbf{r}-\mathbf{r}_\beta| = R_{\text{out}}$, centered at that disk. The joint problem's boundary is instead the circle $|\mathbf{r}|=R_{\text{out}}$ centered at the origin. For centers $\mathbf{r}_A=(a,0)$ and $\mathbf{r}_B=(-a,0)$, the exterior contribution from each disk is
$$
T_c(x)=A\ln\frac{R_{\text{out}}}{|x-c|},
\qquad
A=\frac{q_0R_0^2}{2\kappa_2}.
$$
On the actual boundary $x=R_{\text{out}}e^{i\theta}$,
$$
T_A+T_B
=A\ln\frac{R_{\text{out}}^2}
{|R_{\text{out}}e^{i\theta}-a|\,|R_{\text{out}}e^{i\theta}+a|},
$$
which is not identically zero. For the documented geometry $a=1$, $R_{\text{out}}=2$, it is positive at $\theta=0$ because the product of distances is $1\cdot3=3<R_{\text{out}}^2$, and negative at $\theta=\pi/2$ because the product is $\sqrt5\sqrt5=5>R_{\text{out}}^2$.

Expanding for $a/R_{\text{out}}\ll1$,
$$
\ln|\mathbf{r}-\mathbf{r}_A|
=\ln R_{\text{out}}
-\sum_{n=1}^{\infty}\frac{(a/R_{\text{out}})^n}{n}\cos(n\theta)
\quad\text{on }|\mathbf{r}|=R_{\text{out}}.
$$
The symmetric disk at $-\mathbf{r}_A$ cancels odd modes, so the leading boundary residual is the quadrupole $n=2$, with amplitude $O((a/R_{\text{out}})^2)$. Thus the shifted-single-disk superposition does not become exact from $R_0/d\to0$ and $R_0/R_{\text{out}}\to0$ alone unless the common-boundary effect is also controlled, e.g. by $d/R_{\text{out}}\to0$ or by using a common-domain Green's function.

**Conclusion:** the doc is incomplete on the asymptotic limit for the total bounded-domain discrepancy.

**2. Leading finite-separation/interface correction.**

The doc's $O(R_0/d)$ statement is nevertheless a real finite-separation mechanism. Near disk $A$, the field from disk $B$ is smooth:
$$
T_B(x)=T_B(c_A)+\nabla T_B(c_A)\cdot(x-c_A)+\cdots,
\qquad
|\nabla T_B(c_A)|\sim \frac{A}{d}
=\frac{q_0R_0^2}{2\kappa_2d}.
$$
The constant term does not affect flux continuity at $\partial A$, but the linear term does. The same scalar field is being added on both sides of disk $A$, while physical flux multiplies its normal derivative by $\kappa_1$ inside and $\kappa_2$ outside. Therefore the other disk's contribution creates an interface-flux jump
$$
[\kappa\partial_n T_B]_{\partial A}
=(\kappa_1-\kappa_2)\partial_nT_B
=O\!\left(\frac{q_0R_0^2}{d}\right).
$$
The self flux scale at disk $A$ is
$$
|\kappa_1T_A'(R_0)|=\frac{q_0R_0}{2},
$$
so the relative inter-disk interface defect is
$$
O\!\left(\frac{q_0R_0^2/d}{q_0R_0}\right)
=O(R_0/d).
$$
This confirms the doc's finite-separation scaling as an interface-flux estimate; it is not the full story for the bounded-domain $L^2$ discrepancy because the outer-boundary residual above can dominate or vary in the opposite direction when $d$ changes at fixed $R_{\text{out}}$.

**Discrepancy summary.** The statements "finite-separation corrections are $O(R_0/d)$" and "the discrepancy should decrease as $d/R_0$ grows" conflate two effects. The $O(R_0/d)$ interface defect is plausible and confirmed. The shifted-superposition approximation in the documented outer disk also has a boundary/image correction controlled by $a/R_{\text{out}}=d/(2R_{\text{out}})$, with leading symmetric residual $O((a/R_{\text{out}})^2)$. If the acceptance metric is total relative $L^2$ error against the numerical bounded-domain solution, the trend cannot be justified by $R_0/d$ alone.

**Proposed resolution.** Clarify `verification.md` Problem 4 before implementation: either define the approximation using the common outer boundary so the test isolates the $O(R_0/d)$ interface interaction, or state that the measured discrepancy includes both $O(R_0/d)$ interface error and boundary/image error controlled by $d/R_{\text{out}}$.

This submission stays at status `in-progress` until the `verification.md` discrepancy is resolved, per the brief's done-definition rules.

### Problem 5 — reentrant corner

**1. $T = r^{2/3}\sin(2\theta/3)$ is harmonic.**

Polar Laplacian $\Delta T = \frac{1}{r}\partial_r(r\,\partial_r T) + \frac{1}{r^2}\partial_\theta^2 T$.

Radial part:
$$
\partial_r T = \tfrac{2}{3}r^{-1/3}\sin(2\theta/3), \quad r\,\partial_r T = \tfrac{2}{3}r^{2/3}\sin(2\theta/3),
$$
$$
\partial_r(r\,\partial_r T) = \tfrac{4}{9}r^{-1/3}\sin(2\theta/3), \quad \tfrac{1}{r}\partial_r(r\,\partial_r T) = \tfrac{4}{9}r^{-4/3}\sin(2\theta/3).
$$

Angular part:
$$
\partial_\theta T = \tfrac{2}{3}r^{2/3}\cos(2\theta/3), \quad \partial_\theta^2 T = -\tfrac{4}{9}r^{2/3}\sin(2\theta/3), \quad \tfrac{1}{r^2}\partial_\theta^2 T = -\tfrac{4}{9}r^{-4/3}\sin(2\theta/3).
$$

Sum:
$$
\Delta T = \tfrac{4}{9}r^{-4/3}\sin(2\theta/3) \;-\; \tfrac{4}{9}r^{-4/3}\sin(2\theta/3) = 0. \;\;✓
$$
Consistent with $Q = 0$.

**2. Edges where $T$ vanishes naturally.** Zeros of $\sin(2\theta/3)$: $2\theta/3 = n\pi$, i.e., $\theta = 3n\pi/2$. With $\theta \in [0, 3\pi/2]$ sweeping the interior angle at the reentrant corner $(1/2, 1/2)$, $T = 0$ at $\theta = 0$ and $\theta = 3\pi/2$. These two rays from the reentrant corner are the two L-shape boundary edges meeting at it:
- the horizontal segment from $(1/2, 1/2)$ to $(1, 1/2)$ (part of $y = 1/2$),
- the vertical segment from $(1/2, 1/2)$ to $(1/2, 1)$ (part of $x = 1/2$).

The remaining four outer L-shape edges ($y = 0$ across the bottom; $x = 1$ from $y=0$ to $1/2$; $x = 0$ across the left; $y = 1$ from $x = 0$ to $1/2$) carry nontrivial Dirichlet values from evaluating $r^{2/3}\sin(2\theta/3)$ on them. Test wiring needs to identify the "zero-edges" by physical-group tag, not by "all of $\partial\Omega$ minus the nonzero set" — the exact-solution evaluator returns zero on the corner-adjacent segments by construction.

**3. Why $L^2$ rate is $4/3$, not $2$, on uniform P1 meshes.**

The reentrant corner has interior angle $\omega = 3\pi/2$, giving the dominant singularity exponent $\alpha = \pi/\omega = 2/3$. The solution behaves as $r^\alpha$ times an angular function near the corner; its second derivatives blow up as $r^{\alpha-2}$, which lie in $L^2$ only after losing fractional regularity. Quantitatively, $T \in H^{1+\alpha-\epsilon}(\Omega)$ for any $\epsilon > 0$ but not in $H^{1+\alpha}$.

Standard interpolation estimate on a quasi-uniform mesh of size $h$:
$$
\|T - I_h T\|_{H^1} \;\le\; C\, h^{\alpha - \epsilon}\,|T|_{H^{1+\alpha - \epsilon}}.
$$
The Aubin–Nitsche duality argument upgrades the $L^2$-error by another factor of $h^{\alpha-\epsilon}$ — the dual problem has the same singular geometry and hence the same regularity loss, so duality buys $h^{\alpha-\epsilon}$ rather than the full $h$ of the smooth case. Combined:
$$
\|T - T_h\|_{L^2} \;\le\; C\, h^{2(\alpha - \epsilon)} = C\, h^{4/3 - 2\epsilon}.
$$
Sending $\epsilon \downarrow 0$, the asymptotic uniform-mesh $L^2$ rate is $4/3$. The optimal rate $2$ from the smooth case is lost because the standard interpolation argument requires $|T|_{H^2}$, which is unavailable here.

This rate is "honest" in the sense that uniform refinement cannot recover it — it is a property of the solution's regularity, not a defect of P1 elements. The inverted assertion in `verification.md` Problem 5 (rate $\ge 1.8$ here *fails* the test) catches setups that secretly mask the corner singularity.
