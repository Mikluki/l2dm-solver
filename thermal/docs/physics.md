# Physics

This document states the problem the code computes, the equations involved, and the regime of physical validity. It is the shared mental model between you, future you, and the agent.

For the derivations, see the source document: `The_heat_conduction_equation_graphene_substrate_EN.md`. This document compresses, it does not re-derive.

## Problem

A graphene sheet on a Si/SiO₂ substrate, surrounded by metallic contacts, is heated by illumination (typically IR). We want its steady-state temperature rise, and specifically the **effective substrate thermal conductivity** $\langle\kappa\rangle$ that characterizes how efficiently heat leaks from the graphene into the substrate.

The deliverable is to reproduce the source document's $\langle\kappa\rangle$ curves (Figures A and B) as a function of sample size $L$ (square) or $R$ (disk), for SiO₂ thicknesses of 5 nm, 50 nm, 300 nm, and 600 nm.

## Physical setup

The stack, from top to bottom:

- **Graphene** — 2D sheet, in-plane conductivity $\sigma_{2d}$ (units W/K, since it's 2D).
- **Metallic contacts** — surrounding the graphene, thickness $t_{\text{met}}$, conductivity $\kappa_{\text{met}}$.
- **Graphene–substrate interface** — Kapitza thermal boundary conductance $G_K$, in W/(m²·K).
- **SiO₂ layer** — thickness $d$, conductivity $\kappa_1 \approx 1$ W/(m·K).
- **Si bulk** — conductivity $\kappa_2 \approx 150$ W/(m·K).

Three channels carry heat away from the graphene: vertical flow through the substrate (Kapitza + SiO₂ + Si), lateral flow along the metallic contacts, and in-plane conduction within the graphene itself.

For typical IR detector geometries (sample radius ~30 μm, room temperature), the source document estimates the resistances:

| Channel | Resistance |
|---|---|
| Kapitza interface | ~7 K/W |
| Si substrate (alone) | ~60 K/W |
| SiO₂ substrate (alone) | ~6000–9000 K/W |
| Metallic contacts | ~9000 K/W |

Vertical substrate flow dominates. The Kapitza resistance is small enough that **graphene and substrate top surface can be taken as isothermal** for typical geometries; this assumption breaks at cryogenic temperatures (where $G_K$ drops) and for sub-5 μm spots.

## Equations

### What we ultimately solve (Part 2)

The physically correct equation for the graphene temperature, in the isothermal-graphene-substrate limit and neglecting in-plane conduction:

$$T_{2d}(\mathbf{r}) = \int d\mathbf{r}'\,g_T(\mathbf{r}-\mathbf{r}')\,p_{em}(\mathbf{r}'),$$

where $p_{em}$ is the absorbed power density (W/m²) and $g_T$ is the substrate Green's function — a closed-form Fourier integral encoding the SiO₂/Si stack:

$$g_T(\mathbf{r}-\mathbf{r}') = \int\frac{d\mathbf{q}}{(2\pi)^2}\,\frac{e^{i\mathbf{q}(\mathbf{r}-\mathbf{r}')}}{|q|\kappa_1}\,\frac{\kappa_2\tanh(|q|d) + \kappa_1}{\kappa_1\tanh(|q|d) + \kappa_2}.$$

This is tier 1 of the source document's four-tier hierarchy (the simplest non-trivial case). The relevant features:

- The kernel is **translation-invariant** but not separable.
- It has a $1/|q|$ short-range singularity, meaning $g_T(\mathbf{r}) \sim 1/|\mathbf{r}|$ as $|\mathbf{r}| \to 0$.
- It interpolates correctly between $\kappa_{\text{SiO}_2}$ (for $|q|d \gg 1$, i.e., small samples) and $\kappa_{\text{Si}}$ (for $|q|d \ll 1$, i.e., large samples).

### What Part 1 solves (the warm-up)

Part 1 does **not** solve the equation above. It solves the scalar PDE

$$-\nabla\cdot(\kappa(\mathbf{r})\,\nabla T) = Q(\mathbf{r})$$

on a 2D domain with mixed Dirichlet/Neumann boundary conditions. This is a different operator (differential, not integral) and a different physical problem (heat conduction in a 2D slab with piecewise material properties).

**This is a warm-up, not physics.** The scalar PDE is solved because:

1. It exercises the mesh, basis, assembly, BC handling, subdomain tagging, and verification harness — all of which Part 2 reuses unchanged.
2. It has closed-form solutions for simple geometries, enabling the verification problems in `verification.md`.
3. It is the same shape of equation (elliptic, 2D, piecewise coefficients) without the kernel-evaluation and singular-quadrature complexity.

The scalar PDE does not correspond to any limit of the substrate problem. It is not a "tier 0" of the hierarchy — it is infrastructure development. Calling it physics would be a category error.

## Symbol glossary

| Symbol | Meaning | Units |
|---|---|---|
| $T(\mathbf{r}, z)$ | substrate temperature | K |
| $T_{2d}(\mathbf{r})$ | graphene temperature | K |
| $T_0$ | ambient temperature | K |
| $\kappa(\mathbf{r}, z)$ | bulk thermal conductivity | W/(m·K) |
| $\kappa_1$ | SiO₂ thermal conductivity ($\approx$ 1) | W/(m·K) |
| $\kappa_2$ | Si thermal conductivity ($\approx$ 150) | W/(m·K) |
| $\langle\kappa\rangle$ | effective substrate conductivity | W/(m·K) |
| $\sigma_{2d}(\mathbf{r})$ | graphene 2D thermal conductivity | W/K |
| $\kappa_{\text{met}}$ | metallic contact bulk conductivity | W/(m·K) |
| $t_{\text{met}}$ | metallic contact thickness | m |
| $d$ | SiO₂ layer thickness | m |
| $R$ | disk sample radius | m |
| $L$ | square sample side length | m |
| $G_K(\mathbf{r})$ | Kapitza thermal boundary conductance | W/(m²·K) |
| $g_T(\mathbf{r}-\mathbf{r}')$ | substrate Green's function | K·m²/W |
| $\mathbf{q}$ | in-plane wave vector | 1/m |
| $p_{em}(\mathbf{r})$ | absorbed electromagnetic power density | W/m² |
| $p(\mathbf{r})$ | heat flux from graphene into substrate | W/m² |
| $Q(\mathbf{r})$ | source term in the scalar warm-up PDE | W/m³ (formally; in practice dimensionless in the warm-up tests) |
| $Z_T$ | thermal resistance (with various subscripts) | K/W |

Note the dimensional distinction between $\sigma_{2d}$ (2D, W/K) and $\kappa$ (3D bulk, W/(m·K)). Conflating them is the most likely unit bug in this domain.

## Regime and assumptions

Active assumptions across both Part 1 and Part 2:

- **Steady state.** $\partial_t T = 0$. No time-dependent or transient physics.
- **Linear materials.** $\kappa$ does not depend on $T$. Valid for the small temperature rises typical in detector operation.
- **Classical heat conduction.** Fourier's law applies; no ballistic effects, no phonon-hydrodynamic regime, no quantum corrections.
- **Graphene-substrate isothermal.** $T_{2d}(\mathbf{r}) \approx T(\mathbf{r}, z=0)$ at the substrate top surface, justified by the Kapitza resistance being much smaller than the substrate resistance. Breaks at low temperatures and small spot sizes.
- **In-plane graphene conduction neglected.** Lateral heat flow within the graphene sheet is small compared to vertical flow into the substrate, for typical geometries.
- **Metallic contacts neglected as a heat-flow channel.** Justified for thin contacts ($t_{\text{met}} \ll R$). Worsens for thick oxide where substrate resistance approaches contact resistance.
- **Substrate as a layered 1D stack in the vertical direction.** No lateral structure to $\kappa(z)$. Reduces the 3D heat equation to a Fourier-domain ODE in $z$ and gives the closed-form $g_T$.
- **Infinite substrate extent in-plane.** The Green's function above assumes the substrate extends to infinity in $x$ and $y$. For the isolated-antenna case (Part 2) this is the right form; for the metasurface case, a regularized periodic version replaces it.

Assumptions that are *not* baked in (i.e., the code should handle these):

- Sample shape is arbitrary (disk, square, polygon, multi-patch).
- Material distribution can be arbitrary in-plane (subdomain-tagged regions with distinct $\kappa$).
- Source $p_{em}$ can be arbitrary in-plane.

## What "done" looks like

The validation target — distinct from the verification problems — is reproduction of the source document's $\langle\kappa\rangle$ curves:

- Square samples (side $L$) and circular samples (radius $R$), sizes from 1 nm to 1 mm, log-spaced.
- SiO₂ thicknesses of 5 nm, 50 nm, 300 nm, 600 nm.
- $\kappa_1 = 1$ W/(m·K), $\kappa_2 = 150$ W/(m·K).
- Curves should sigmoidally interpolate between $\kappa_1$ (small samples) and $\kappa_2$ (large samples), with crossover at sample size of order $(\kappa_2/\kappa_1)\,d \approx 150\,d$.

This is a Part 2 deliverable. Part 1 produces no physics output; it produces the harness and the scalar solver that Part 2 will reuse.
