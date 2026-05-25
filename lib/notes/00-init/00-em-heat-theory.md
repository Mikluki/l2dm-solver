# Layered Photodetector Solver — Meeting Notes

## 1. Task Formulation

We are building a finite-element solver for layered photodetector structures. A
typical structure is:

```
    Metal (2D) / Graphene (2D)          ← active layer
    ─────────── ε₁ ───────────          ← dielectric stack
    ─────────── ε₂ ───────────             (contrasting permittivities,
    ─────────── ε₃ ───────────              arbitrary number of layers)
    ─────────── substrate ────
```

Two coupled physical problems must be solved on the same geometry:

1. **Electrodynamics** — given an incident wave `E₀(x, y)`, find the surface
   current `j(x, y)` and the electric field on the 2D layer.
2. **Heat transport** — given the absorbed power from (1) as a source, find the
   temperature distribution `T_surf(x, y)` on the 2D layer, coupled to the 3D
   bulk temperature of the substrate.

In both problems the 2D active layer lies on a **contrasting substrate**
(layered dielectric / contrasting thermal conductivity), which is the main
complication versus the free-space case.

The implementation strategy is the method of finite elements applied to an
**integral** formulation of each problem, where the substrate has been
eliminated analytically via a Green's function / transfer-matrix construction.

---

## 2. Electrodynamic Equation

### 2.1 Integral equation on the 2D layer

The unknown is the surface current `j(r)` (equivalently, the tangential field
`E(r)` at `z = 0`). The equation closes on the 2D plane:

$$
\vec{E}(\vec{r}) \;=\; \vec{E}_0(\vec{r}) \;+\; \int
\bigl(k_0^2 \,-\, \operatorname{grad\,div}\bigr)\,
G(\vec{r} - \vec{r}')\,\vec{j}(\vec{r}')\, d^2\vec{r}'
$$

closed by the local constitutive relation

$$
\vec{E}(\vec{r}) \;=\; \frac{\vec{j}(\vec{r})}{\sigma(\vec{r})}.
$$

**Term meanings:**

| Symbol | Meaning |
|---|---|
| `E(r)` | Total in-plane electric field at point `r = (x, y)` on the 2D layer. The primary unknown when σ is finite. |
| `E₀(r)` | Incident (driving) electromagnetic wave evaluated on the plane `z = 0`. Known input. |
| `j(r')` | Surface current density (A/m) on the 2D layer at source point `r'`. Alternative primary unknown; preferred when the layer contains perfect conductors (see below). |
| `G(r − r')` | Scalar Green's function of the substrate stack — the field response at `r` to a point current at `r'` with the full layered dielectric present. Generalises the free-space Green's function `G₀ = e^{ik₀ R}/(4π R)`, where `R = \|r − r'\|`. |
| `k₀²` term | Radiation / "mass" term from the vector wave equation — current oscillating at ω₀ radiates through the photon wavenumber `k₀ = ω/c`. |
| `grad div` term | Longitudinal correction: accounts for the scalar potential produced by charge accumulation `ρ = i(∇·j)/ω`. Without it the equation would only describe transverse radiation. |
| `σ(r)` | Local 2D conductivity of the active material (graphene, 2DEG, metal). Spatially varying: different in graphene regions, metal contacts, gaps. |

**Why formulate in terms of `j` rather than `E`:** metal contacts have `σ → ∞`,
so `E = j/σ → 0` locally. Rearranging the equation as `j/σ = E₀ + ∫(…)j`, the
term `j/σ` simply vanishes inside metals. The same equation therefore covers
graphene regions and perfect-metal contacts uniformly.

### 2.2 Fourier / q-space representation

Because the substrate is translation-invariant in `(x, y)`, the Green's
function is a convolution kernel, and the equation becomes a pointwise product
in the in-plane wavevector `q`:

$$
\vec{E}(\vec{q}) \;=\; \hat{G}(\vec{q}) \, \vec{j}(\vec{q}).
$$

`Ĝ(q)` is what the transfer-matrix method delivers.

### 2.3 Helmholtz decomposition and polarization split

The transfer matrix has **different** expressions for p- and s-polarized waves,
so the vector current must be split before applying it. Any vector field
decomposes as

$$
\vec{j}(\vec{q}) \;=\;
\underbrace{\frac{(\vec{q}\cdot\vec{j})\,\vec{q}}{q^2}}_{\text{longitudinal}}
\;+\;
\underbrace{\left(\vec{j} - \frac{(\vec{q}\cdot\vec{j})\,\vec{q}}{q^2}\right)}_{\text{transverse}},
$$

or in tensor form

$$
j_i \;=\; \frac{q_i q_j}{q^2}\, j_j \;+\; \left(\delta_{ij} - \frac{q_i q_j}{q^2}\right) j_j.
$$

**Term meanings and physical mapping:**

| Piece | Geometry | Polarization | Transfer matrix to use |
|---|---|---|---|
| Longitudinal: `j ∥ q` | Current oscillates along its propagation direction → charge accumulation → scalar potential → field has a component in the plane of incidence. | **p-polarization** | p-wave transfer matrix. |
| Transverse: `j ⊥ q` | Current oscillates perpendicular to `q` → no charge pile-up → purely solenoidal field. | **s-polarization** | s-wave transfer matrix. |

For simple geometries (e.g. a slit + plane wave) the polarization is known *a
priori* and the split is trivial. For arbitrary in-plane geometries (triangle,
asymmetric grating) the current has mixed components and the decomposition
must be applied numerically in q-space before multiplying by `Ĝ(q)`.

### 2.4 Solver skeleton (pseudocode)

```
INPUT:  substrate stack { ε_l, d_l }_{l=1..N}
        incident wave   E₀(x, y)
        geometry of active layer + conductivity map σ(x, y)

1. Build Ĝ_p(q), Ĝ_s(q) by transfer matrix over the stack.
2. Triangulate the active 2D region; choose basis functions on triangles.
3. Assemble the matrix form of the integral equation using the
   Helmholtz split: for each q, project j onto q̂ and q̂_⊥, apply the
   corresponding Ĝ, recombine.
4. Solve the resulting dense linear system for j(x, y).
5. Inverse-transform j(q) if fields off the plane are needed.
6. Output: E(x, y), j(x, y), absorbed power density P(x, y) = 1/2 Re(j·E*).
```

---

## 3. Heat-Transport Equation

### 3.1 Coupled surface + bulk equations

Define

- `T_surf(R)` — temperature of the 2D layer, `R = (x, y)`.
- `T_bulk(R, z)` — temperature in the 3D substrate.
- `T_surf(R) = T_bulk(R, z → 0⁺)` only up to the Kapitza-like jump below.

**2D surface heat equation:**

$$
\nabla_{2D}\!\cdot\!\bigl(\kappa_{\text{surf}}(\vec{R})\, \nabla_{2D} T_{\text{surf}}\bigr)
= - Q_{\text{surf}}(\vec{R})
+ \frac{C}{\tau_e}\bigl(T_{\text{surf}}(\vec{R}) - T_{\text{bulk}}(\vec{R}, z\to 0^+)\bigr).
$$

**3D bulk heat equation:**

$$
\nabla_{3D}\!\cdot\!\bigl(\kappa_{\text{bulk}}(\vec{r})\, \nabla_{3D} T_{\text{bulk}}(\vec{r})\bigr)
= 0.
$$

**Coupling boundary condition at `z = 0`:**

$$
-\,\kappa_{\text{bulk}}\, \nabla T_{\text{bulk}} \cdot \hat{n}
\;=\; \frac{C}{\tau_e}\bigl(T_{\text{surf}} - T_{\text{bulk}}|_{z\to 0^+}\bigr).
$$

**Term meanings:**

| Symbol | Meaning |
|---|---|
| `κ_surf(R)` | In-plane (sheet) thermal conductivity. Spatially varying — high under metal contacts, different under graphene, zero in gaps. |
| `κ_bulk(r)` | 3D thermal conductivity of the substrate. Varies with layer: e.g. SiO₂ (low κ) then Si (high κ). |
| `Q_surf(R)` | Absorbed power per area on the 2D layer. **Input from the EM solver** (P(x, y) computed in §2). |
| `C / τ_e` | Surface-to-bulk thermal coupling. `C` = sheet heat capacity of the 2D layer, `τ_e` = energy relaxation time to the substrate. Mechanistically close to inverse Kapitza resistance: finite `(T_surf − T_bulk)` is required to push heat across the interface. |
| `T_bulk(R, z → 0⁺)` | Substrate temperature just below the 2D layer. In general `≠ T_surf`. |
| Right-hand-side sign | Positive `Q_surf` heats the layer; positive `(T_surf − T_bulk)` cools it by leaking heat into the bulk. |

In non-stationary problems add `c_p ∂T/∂t` on the left-hand side (bulk heat
capacity × warming rate). In this project we start stationary.

### 3.2 Reducing the bulk: pure surface integral equation

The bulk equation with the coupling BC is linear — solve it "in the head" via
the heat Green's function and substitute back. The bulk Green's function of
the half-space Laplacian is Coulomb-like:

$$
G_{\text{heat}}(\vec{r}, \vec{r}') \;\propto\; \frac{1}{4\pi\,\|\vec{r}-\vec{r}'\|},
$$

with the half-space handled by the method of images. The elimination yields a
**2D integral equation** for `T_surf(R)` alone — structurally analogous to the
EM equation in §2, but scalar (simpler).

**Order-of-magnitude estimate** for a hot spot of radius `R₀` on a
semi-infinite substrate:

$$
T_{\text{surf}} \;\sim\; \frac{Q_{\text{surf}}}{\kappa_{\text{substrate}} / R_0},
$$

which is the thermal analogue of Ohm's law with `ΔT / R₀` replacing `∇T`. The
Green's-function method supplies the O(1) prefactor exactly.

### 3.3 Boundary conditions on the finite computational domain

- **Substrate bottom / sides (infinity):** fixed at the cryostat / ambient
  temperature. In integral form this is absorbed automatically by choosing the
  Green's function that decays at infinity.
- **Edges of finite metal / graphene patches:** zero normal heat flux,
  `q · n̂ = 0`. Heat leaves through the substrate, not through the zero-area
  edge of a 2D film.
- **Periodic metasurfaces:** periodic BCs on the unit cell.
- Special case: when a metal electrode is itself a large heat sink
  comparable to the substrate (e.g. long wires to the chip boundary), the
  zero-flux BC must be replaced with a sink BC on the metal — this is
  physically important for real detectors.

### 3.4 Note on 1D reduction

Unlike the EM problem, the 1D case of the heat equation is **pathological** —
a line source on a 1D substrate has no transverse dimension to sink heat into
and diverges. So the heat solver must be exercised in ≥ 2D. The 2D surface
problem with 3D substrate stays finite because the 3D substrate is a working
heat sink.

### 3.5 Kapitza-resistance generalization (transfer-matrix analogue)

Just as the EM transfer matrix was extended to treat a 2D conducting sheet by
imposing a jump in `H` at the interface, the thermal version can carry a jump
in `T` governed by the Kapitza resistance `R_K`:

$$
T_{\text{above}} - T_{\text{below}} \;=\; R_K \cdot q_{\text{flux}},
$$

so layered-substrate Green's functions with finite inter-layer thermal
resistance are a straightforward extension.

---

## 4. Conversation Story (Q & A)

Filtered to substantive questions and answers; interpersonal exchanges
omitted.

### 4.1 Framing

- **Goal.** Build a finite-element solver for 2D active layers (metal /
  graphene) on contrasting layered substrates. First priority:
  electrodynamics. Parallel priority: thermal conductivity, because
  detectors operate on a thermal (bolometric) principle.
- **Why integral form.** A fully 3D `(x, y, z)` formulation is expensive in
  `z`. Reformulating as a 2D surface equation with the substrate eliminated
  via Green's function / transfer matrix keeps the unknown on the plane.

### 4.2 Electrodynamics — questions raised and resolved

**Q: The operator `(k₀² − grad div)` acts on a vector `j`; are the tensor
ranks consistent?**
A: Yes. `k₀² j` is vectorial. `div j` is scalar, but `grad(div j)` is
vectorial again. Split and written as two separate terms the structure is
obviously vector-valued.

**Q: What happens at perfect metals where `σ → ∞`?**
A: Formulate the integral equation with `j/σ` on the left. In metal regions
`j/σ → 0` and the term drops out automatically. The same equation covers
graphene and metal contacts without special cases — this is the main
argument for the current-based formulation.

**Q: The transfer-matrix method looks different for p and s waves. How do we
know which to apply at a given point?**
A: Project the current in q-space via Helmholtz decomposition:
longitudinal part (`j ∥ q`) → p-polarization transfer matrix; transverse
part (`j ⊥ q`) → s-polarization transfer matrix. For simple geometries
polarization is known in advance and no split is needed; for arbitrary
in-plane geometries both components are generically present and must be
handled per-`q`.

**Q: Does the decomposition need to be performed in real space?**
A: No — stay in Fourier space throughout. The transfer matrix itself is
applied in q-space, and Helmholtz decomposition is algebraically trivial
there. Real-space inversion is only needed at output.

### 4.3 Validation / division of labour (EM side)

- **1D first.** Existing spectral-method result for a single strip on GaAs
  (D. Mylnikov's code, already benchmarked against experiment by Muravyov)
  will serve as the reference for a 1D finite-element port. Implementing 1D
  first is to nail down the basis-function / triangulation machinery before
  handling vector-valued 2D problems.
- **Responsibilities.** Kirill: pseudocode + derivation of the EM integral
  equation & project lead. Vlad: 1D EM FEM implementation, validated against the strip
  reference. CST Microwave Studio reference solutions (square scatterer,
  tooth scatterer, strip-on-stack) to be produced by Ilya / Sasha. Mikhail:
  aggregation.

### 4.4 Heat equation — questions raised and resolved

**Q: Why does the bulk equation have zero on the right-hand side — isn't
heat being injected?**
A: Heat injection happens on the 2D surface, represented either as a
boundary condition or as a `δ(z)` source. In the stationary case there is
no volumetric source inside the substrate (barring dielectric losses, see
extension). Non-stationarity adds `c_p ∂T/∂t`.

**Q: Can dielectric layers in the substrate themselves heat up?**
A: Yes — add `Q_bulk(R, z)` on the right-hand side for any volumetric
source. The integral treatment generalises with essentially the same
machinery.

**Q: Does the transfer-matrix approach extend to the heat problem if there
is Kapitza (thermal contact) resistance between layers?**
A: Yes. The EM transfer matrix handled a 2D sheet by jumping `H` at the
interface; the thermal version handles inter-layer thermal resistance by
jumping `T` at the interface, with `ΔT = R_K · q_flux`. Structure of the
layered Green's function is unchanged.

**Q: Where does heat physically leave a finite structure?**
A: Not through the zero-area edges of a 2D film — through the substrate
bottom (heat sink at infinity) and through the leads (metal wires to the
chip periphery). These are the two real heat sinks.

**Q: Can we reduce to 1D for testing, as we do in EM?**
A: No — a 1D heat problem has nowhere to sink heat and temperature
diverges. The heat solver must be exercised in at least 2D surface + 3D
bulk. This is fine: the 3D substrate is a genuine heat sink.

**Q: Bounding-box / BC sensitivity (CST experience shows results depend on
domain padding)?**
A: Integral formulation removes this: the Green's function encodes the
correct decay at infinity. Finite-domain BCs are needed only when a real
second heat sink exists (e.g. highly conductive leads), in which case a
sink BC on those leads is imposed explicitly.

### 4.5 Heat-solver responsibilities

- Svintsov: derivation of the surface integral heat equation eliminating
  the bulk.
- Mikhail: first implement a **simpler scalar differential** problem —
  `∇·(κ(R) ∇T) = Q` on a 2D square with zero-flux BCs and spatially
  varying `κ` imitating metal/graphene patches — to develop the
  triangulation and basis-function infrastructure. Then migrate to the
  integral form once derived.
- Kirill: aggregate programmatic side (FEM assembly, meshing via
  off-the-shelf triangulation libraries).

### 4.6 Open qualitative question (for order-of-magnitude work on paper)

In a real detector there are two competing heat sinks: (a) the substrate
(bottom-infinity), (b) the contact wires leading to the chip edge. Which
dominates?

- Wire geometry: width ~ tens of μm, thickness ~ hundreds of nm, length to
  periphery ~ 1 mm, κ of Au/Cu known.
- Substrate: heat spreading from a small spot into a half-space, `T ~ Q
  R₀ / κ_sub`.

Compute both thermal resistances, compare. This answer guides optimization
(thinner/longer leads vs. low-κ substrate vs. thinner substrate) **before**
the solver is finished (solver ETA ≥ 1 month).

Current community heuristic (to be confirmed numerically): substrate
dominates — most reported bolometer work tunes the substrate aggressively
while leaving leads alone.

### 4.7 Benchmarks to establish

| # | Structure | Reference |
|---|---|---|
| 1 | Strip on layered substrate | Mylnikov's spectral code + experiment |
| 2 | Square 2DEG patch | CST Microwave Studio |
| 3 | Tooth / slot scatterer | CST Microwave Studio |
| 4 | Asymmetric dual grating gate (future) | CST Microwave Studio |

Target: results from the new FEM solver within known tolerance of these
references, then exercise on geometries the spectral method cannot handle
(asymmetric, arbitrary in-plane shapes).

---

## 5. Summary of Decisions

1. **Two solvers**, sharing geometry and FEM infrastructure:
   EM (vector, q-space transfer-matrix Green's function) →
   heat (scalar, half-space heat Green's function + Kapitza-like coupling).
2. **Formulate EM in surface currents**, not fields, so perfect metals need
   no special treatment.
3. **Eliminate the bulk analytically** in both problems before discretizing,
   keeping the unknown on the 2D plane.
4. **Validate in 1D** on the EM side against existing spectral results;
   **skip 1D** on the heat side (divergent) and start directly at 2D+3D.
5. **Pipeline.** EM solver outputs `Q_surf(R)` → heat solver inputs it.
   Each solver also has independent value (heat solver can ingest existing
   CST absorbed-power maps).
