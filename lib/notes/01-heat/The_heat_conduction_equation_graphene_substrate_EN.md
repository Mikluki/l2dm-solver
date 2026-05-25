# Derivation of the heat conduction equation for the graphene–contacts–substrate system

## Substrate temperature for a prescribed surface heat source

The temperature of a Si/SiO$_2$ substrate is obtained by solving the heat conduction equation

$$-\nabla\bigl[\kappa(z)\,\nabla T\bigr] = 0,$$

with a prescribed heat flux on the upper surface

$$-\kappa(0)\,\frac{dT}{dz} = p(x,y).$$

The dissipated power $p(x,y)$ (in W/m$^2$) is treated as given for the time being.

We take a Fourier transform in the in-plane coordinates $x$ and $y$,

$$T(\mathbf{q},z) = \int dx\,dy\;e^{-i\mathbf{q}\mathbf{r}_\parallel}\,T(\mathbf{r}_\parallel,z),\qquad p(\mathbf{q}) = \int dx\,dy\;e^{-i\mathbf{q}\mathbf{r}_\parallel}\,p(\mathbf{r}_\parallel),$$

so that the heat equation reduces to an ordinary differential equation in $z$:

$$q^2\kappa(z)\,T(\mathbf{q},z) - \frac{d}{dz}\!\left[\kappa(z)\,\frac{dT(\mathbf{q},z)}{dz}\right] = 0,\qquad -\kappa(0)\,\frac{dT(\mathbf{q},z)}{dz}\bigg|_{z=0} = p(\mathbf{q}).$$

At each internal layer interface (if any) we impose continuity of temperature and of heat flux. For instance, at the Si/SiO$_2$ interface located at $z=-d$ below the graphene,

$$T(\mathbf{q},z=-d_+) = T(\mathbf{q},z=-d_-),\qquad \kappa_1 T'(\mathbf{q},z=-d_+) = \kappa_2 T'(\mathbf{q},z=-d_-),$$

where $\kappa_1$ is the thermal conductivity of the upper layer (SiO$_2$) and $\kappa_2$ that of the lower layer (Si). Once $T(\mathbf{q},z=0)$ is found, the graphene temperature in real space follows by inverse Fourier transform; the mean graphene temperature is obtained by an additional area average over the sheet.

The Fourier-space solution at the top surface reads

$$T(\mathbf{q},z=0) = \frac{p(\mathbf{q})}{|q|\kappa_1}\,\frac{\kappa_2\tanh(|q|d) + \kappa_1}{\kappa_1\tanh(|q|d) + \kappa_2}.$$

As a sanity check: at $d=0$ the upper-layer conductivity $\kappa_1$ (SiO$_2$) drops out, and at $d\to\infty$ the lower-layer conductivity $\kappa_2$ (Si) drops out. Inverting the Fourier transform yields a Green's-function representation,

$$T(\mathbf{r},z=0) = \int d\mathbf{r}'\,p(\mathbf{r}')\,g_T(\mathbf{r}-\mathbf{r}'),$$

$$g_T(\mathbf{r}-\mathbf{r}') = \int\frac{d\mathbf{q}}{(2\pi)^2}\,\frac{e^{i\mathbf{q}(\mathbf{r}-\mathbf{r}')}}{|q|\kappa_1}\,\frac{\kappa_2\tanh(|q|d) + \kappa_1}{\kappa_1\tanh(|q|d) + \kappa_2}.$$

## Coupling the substrate-surface heat flux to the graphene temperature

So far we have computed the substrate-surface temperature assuming a known heat flux $p(\mathbf{r}')$ from the graphene. We now allow the graphene temperature $T_{2d}(\mathbf{r})$ to differ from the substrate top-surface temperature $T(\mathbf{r},z=0)$ on account of the Kapitza resistance. The *surface* heat conduction equation governing $T_{2d}(\mathbf{r})$ is

$$-\nabla\bigl[\sigma_{2d}(\mathbf{r})\,\nabla T_{2d}(\mathbf{r})\bigr] = p_{em}(\mathbf{r}) - p(\mathbf{r}),$$

where $p(\mathbf{r})$ accounts for heat leaking from the graphene into the substrate, $p_{em}(\mathbf{r})$ is the electromagnetic absorption density, and $\sigma_{2d}(\mathbf{r})$ is the position-dependent 2D thermal conductivity of the graphene (the metallic-contact contribution can also be folded into this distribution). We deliberately use a different symbol for the 2D conductivity to highlight the change in dimensions: W/K for the 2D quantity versus W/(m·K) for the bulk one. Heat outflow into the substrate is taken to be linear, with proportionality coefficient set by the Kapitza conductance $G_K$:

$$p(\mathbf{r}) = G_K(\mathbf{r})\bigl[T_{2d}(\mathbf{r}) - T(\mathbf{r},z=0)\bigr].$$

The result is a coupled system for the graphene temperature and the substrate top-surface temperature:

$$\begin{cases}
-\nabla\bigl[\sigma_{2d}(\mathbf{r})\,\nabla T_{2d}(\mathbf{r})\bigr] = p_{em}(\mathbf{r}) - G_K(\mathbf{r})\bigl[T_{2d}(\mathbf{r}) - T(\mathbf{r},z=0)\bigr],\\[4pt]
T(\mathbf{r},z=0) = \displaystyle\int d\mathbf{r}'\,G_K(\mathbf{r}')\,g_T(\mathbf{r}-\mathbf{r}')\bigl[T_{2d}(\mathbf{r}') - T(\mathbf{r}',z=0)\bigr].
\end{cases}$$

Two simplifying limits are worth noting:

(a) If the 2D-system properties are independent of $\mathbf{r}$, the planar Fourier transform handles the whole problem in closed form. We won't pursue this further here.

(b) If the Kapitza conductance is the largest scale in the problem, the graphene and the substrate top surface can be taken at the same temperature, $T_{2d}(\mathbf{r}) \approx T(\mathbf{r},z=0)$. The system then collapses to a single integral equation:

$$T_{2d}(\mathbf{r}) = \int d\mathbf{r}'\,g_T(\mathbf{r}-\mathbf{r}')\Bigl[p_{em}(\mathbf{r}') + \nabla\bigl[\sigma_{2d}(\mathbf{r}')\,\nabla T_{2d}(\mathbf{r}')\bigr]\Bigr].$$

## Numerical estimates and limiting cases

We now make some order-of-magnitude estimates to decide which heat-leakage channels matter. Take a graphene sample with characteristic size $R$ dissipating a total power $P$. Balancing supply and outflow through the substrate gives

$$P = \pi R^2\,\kappa_{sub}\,\frac{dT}{dr} \approx \pi R^2\,\kappa_{sub}\,\frac{T_{2d}-T_0}{R},$$

so the substrate thermal resistance is of order

$$Z_{T,sub} \sim \frac{1}{\pi R\,\kappa_{sub}}.$$

The exact Fourier calculation outlined above yields

$$Z_{T,sub} = \frac{8}{3\pi^2\,R\,\kappa_{sub}}.$$

Taking $\kappa_{sub} = 150$ W/(m·K) for silicon and $R = 30$ μm — a typical spot radius in IR experiments — we get

$$Z_{T,\,\mathrm{Si\,sub}} = 60 \text{ K/W}.$$

If we instead use the conductivity of silicon oxide, the answer is roughly 100–150 times larger:

$$Z_{T,\,\mathrm{SiO}_2\,\mathrm{sub}} = 6000\ldots 9000 \text{ K/W}.$$

The effective conductivity of the actual composite substrate is recovered by averaging the above expression over wave vectors:

$$\frac{1}{\langle\kappa\rangle} = \frac{3\pi}{8}\,\frac{1}{\kappa_1}\int_{-\infty}^{+\infty}\frac{d\tau}{2}\left(\frac{2J_1(\tau)}{\tau}\right)^{\!2}\,\frac{\kappa_1 + \kappa_2\tanh(\tau d/R)}{\kappa_2 + \kappa_1\tanh(\tau d/R)}.$$

Mak, Lui and Heinz [K. F. Mak, C. H. Lui, and T. F. Heinz, "Measurement of the thermal conductance of the graphene/SiO2 interface," *Appl. Phys. Lett.*, vol. 97, no. 22, Nov. 2010, doi: 10.1063/1.3511537] report a graphene–SiO$_2$ Kapitza conductance $G_K = 5\times 10^7$ W/(m$^2$·K). For a circular sample of radius $R = 30$ μm this yields a Kapitza resistance

$$Z_K = \frac{1}{\pi R^2\,G_K} = 7 \text{ K/W},$$

which is much smaller than the substrate resistance. This vindicates our assumption that the graphene and the substrate top surface are essentially isothermal. The assumption can break down at low temperatures, where $G_K$ drops sharply, and for smaller spot radii (already around 5 μm) $Z_K$ and $Z_{T,sub}$ may become comparable.

Lastly, we estimate the contribution of the metallic contacts $Z_{met}$. Surround the graphene by an annular metallic contact of inner radius $R = 30$ μm and outer radius $R_{max}\gg R$. Heat leakage along the contact surface is equivalent to a Corbino-disk resistance problem, giving

$$Z_{met} = \frac{\ln(R_{max}/R)}{2\pi\,\kappa_{met}\,t_{met}},$$

where $t_{met}$ is the vertical thickness of the contact. With $R_{max} = 1$ mm, $\kappa_{met} = 300$ W/(m·K), $t_{met} = 200$ nm we obtain

$$Z_{met} \approx 9\times 10^3 \text{ K/W}.$$

The contact resistance is therefore very large, and lateral heat conduction along the contacts can safely be neglected next to the substrate channel — ultimately because of the small contact thickness, $t_{met}\ll R$. The approximation is less safe for thick SiO$_2$ layers, where substrate and contact conductances become comparable. In any case, in-plane conduction through the graphene itself can be neglected with confidence, as it is far smaller than through the metal contacts.

This gives the following hierarchy of approximations for the graphene temperature under irradiation:

1. The simplest limit ignores the Kapitza resistance and lateral conduction through both the metals and the graphene. The surface temperature is then obtained directly from

   $$T_{2d}(\mathbf{r}) = \int d\mathbf{r}'\,g_T(\mathbf{r}-\mathbf{r}')\,p_{em}(\mathbf{r}').$$

   Strictly speaking this is no longer an equation but an explicit formula expressing the surface temperature in terms of the known electromagnetic power.

2. For strongly insulating substrates (e.g., a thick SiO$_2$ layer of order 1 μm) or thick metallic contacts (also of order 1 μm), lateral conduction through the contacts must be kept and one has to solve the integral equation

   $$T_{2d}(\mathbf{r}) = \int d\mathbf{r}'\,g_T(\mathbf{r}-\mathbf{r}')\Bigl[p_{em}(\mathbf{r}') + \nabla\bigl[\sigma_{2d}(\mathbf{r}')\,\nabla T_{2d}(\mathbf{r}')\bigr]\Bigr].$$

3. At cryogenic temperatures or for very small samples, the finite Kapitza resistance must also be retained, and one solves the full system

   $$\begin{cases}
   -\nabla\bigl[\sigma_{2d}(\mathbf{r})\,\nabla T_{2d}(\mathbf{r})\bigr] = p_{em}(\mathbf{r}) - G_K(\mathbf{r})\bigl[T_{2d}(\mathbf{r}) - T(\mathbf{r},z=0)\bigr],\\[4pt]
   T(\mathbf{r},z=0) = \displaystyle\int d\mathbf{r}'\,G_K(\mathbf{r}')\,g_T(\mathbf{r}-\mathbf{r}')\bigl[T_{2d}(\mathbf{r}') - T(\mathbf{r}',z=0)\bigr].
   \end{cases}$$

   In this regime the substrate plays a less prominent role than the surface conduction.

4. At ultra-low temperatures (strong decoupling of electrons from the substrate), the substrate can simply be assumed to stay at ambient temperature, and

   $$-\nabla\bigl[\sigma_{2d}(\mathbf{r})\,\nabla T_{2d}(\mathbf{r})\bigr] = p_{em}(\mathbf{r}) - G_K(\mathbf{r})\,T_{2d}(\mathbf{r}).$$

## Implementation: isolated antennas versus metasurfaces

The Green's function of the heat equation takes a slightly different form for an isolated antenna-coupled detector and for a "periodic metasurface" detector.

For an isolated antenna-coupled detector, the integration runs over the area of all heat-conducting elements on the surface — i.e., the metals and the graphene:

$$T_{2d}(\mathbf{r}) = \int_{\mathrm{Gr,Me}} d\mathbf{r}'\,g_T(\mathbf{r}-\mathbf{r}')\Bigl[p_{em}(\mathbf{r}') + \nabla\bigl[\sigma_{2d}(\mathbf{r}')\,\nabla T_{2d}(\mathbf{r}')\bigr]\Bigr],$$

$$\sigma_{2d}(\mathbf{r}') = \begin{cases} \sigma_{Gr}, & \mathbf{r}'\in\mathrm{Gr},\\ \kappa_{Me}\,t_{Me}, & \mathbf{r}'\in\mathrm{Me}.\end{cases}$$

The integral equation thus lives on a bounded region, with no need to impose boundary conditions at infinity. The Green's function takes its standard form,

$$g_T(\mathbf{r}-\mathbf{r}') = \int\frac{d\mathbf{q}}{(2\pi)^2}\,\frac{e^{i\mathbf{q}(\mathbf{r}-\mathbf{r}')}}{|q|\kappa_1}\,\frac{\kappa_2\tanh(|q|d) + \kappa_1}{\kappa_1\tanh(|q|d) + \kappa_2}.$$

For a metasurface, the equation is solved inside a single unit cell,

$$T_{2d}(\mathbf{r}) = \int_{\mathrm{Unit\,cell}} d\mathbf{r}'\,\tilde g_T(\mathbf{r}-\mathbf{r}')\Bigl[p_{em}(\mathbf{r}') + \nabla\bigl[\sigma_{2d}(\mathbf{r}')\,\nabla T_{2d}(\mathbf{r}')\bigr]\Bigr],$$

but the Green's function is now a Fourier *series* — over wave vectors that are multiples of the reciprocal-lattice vector — rather than a Fourier integral:

$$\tilde g_T(\mathbf{r}-\mathbf{r}') = \frac{1}{L_x L_y}\sum_{n_x,n_y=-\infty}^{+\infty}\frac{e^{i\mathbf{q}(\mathbf{r}-\mathbf{r}')}}{|q|\kappa_1}\,\frac{\kappa_2\tanh(|q|d) + \kappa_1}{\kappa_1\tanh(|q|d) + \kappa_2},\qquad \mathbf{q} = \frac{2\pi n_x}{L_x}\mathbf{e}_x + \frac{2\pi n_y}{L_y}\mathbf{e}_y.$$

A problem now appears: when both indices vanish, $n_x = n_y = 0$, the corresponding term in the sum is strictly infinite. This is not accidental — heating an infinite planar structure at a constant power density must produce an infinite temperature rise. But what we actually care about is not the overall temperature rise, only its gradients within a unit cell. The infinite term can therefore be dropped, and we define a Green's function with the zero mode removed:

$$\tilde g_{T,\mathrm{reg}}(\mathbf{r}-\mathbf{r}') = \frac{1}{L_x L_y}\sum_{\substack{n_x,n_y=-\infty\\ n_x\,\cup\,n_y\neq 0}}^{+\infty}\frac{e^{i\mathbf{q}(\mathbf{r}-\mathbf{r}')}}{|q|\kappa_1}\,\frac{\kappa_2\tanh(|q|d) + \kappa_1}{\kappa_1\tanh(|q|d) + \kappa_2},\qquad \mathbf{q} = \frac{2\pi n_x}{L_x}\mathbf{e}_x + \frac{2\pi n_y}{L_y}\mathbf{e}_y.$$

To make the subtraction procedure more rigorous, consider a metasurface illuminated by a laser beam whose intensity decays slowly from center to edges. For the metasurface as a whole (not for a single cell) we have a "large-scale" equation

$$T_{2d}(\mathbf{r}) = \int_{\mathrm{metasurface}} d\mathbf{r}'\,g_T(\mathbf{r}-\mathbf{r}')\,q(\mathbf{r}').$$

Break the integral into a sum over unit cells:

$$T_{2d}(\mathbf{r}) = \left(\int_{\mathrm{cell\,1}} + \int_{\mathrm{cell\,2}} + \ldots + \int_{\mathrm{cell\,N}}\right)\!\Bigl\{d\mathbf{r}'\,g_T(\mathbf{r}-\mathbf{r}')\,q(\mathbf{r}')\Bigr\}.$$

We now reduce the integration to a single unit cell. A radius-vector shift transforms the Green's function as

$$g_T(\mathbf{r}-\mathbf{r}'+\mathbf{L}\times n) = \int\frac{d\mathbf{q}}{(2\pi)^2}\,g_T(\mathbf{q})\,e^{i\mathbf{q}(\mathbf{r}-\mathbf{r}'+\mathbf{L}\times n)},$$

i.e., a real-space shift multiplies its Fourier components by a phase factor $e^{i\mathbf{q}\mathbf{L}\times n}$. The same shift makes the dissipation function nearly periodic:

$$q(\mathbf{r}+\mathbf{L}\times n) = q(\mathbf{r})\,\exp(-\gamma n).$$

For an infinite laser spot on an infinite metasurface we have exact periodicity, $\gamma = 0$. For a finite spot, $\gamma \ll 1$; one can estimate $\gamma\sim L/R_{beam}$, where $L$ is the unit-cell size and $R_{beam}$ is the spot size. Substituting the periodicity rules gives

$$T_{2d}(\mathbf{r}) = \int_{\mathrm{cell\,1}} d\mathbf{r}'\,G_T(\mathbf{r}-\mathbf{r}')\,q(\mathbf{r}'),$$

$$G_T(\mathbf{r}-\mathbf{r}') = \sum_{n_x,n_y}\int\frac{d\mathbf{q}}{(2\pi)^2}\,e^{i\mathbf{q}(\mathbf{r}-\mathbf{r}')}\,g_T(\mathbf{q})\,\exp\bigl(iq_x n_x L_x - \gamma|n_x|\bigr)\,\exp\bigl(iq_y n_y L_y - \gamma|n_y|\bigr).$$

The discrete sum evaluates by the identity

$$2\sum_{n=0}^{\infty}\cos(qnL)\,e^{-\gamma n} = 1 + \frac{\sinh\gamma}{\cosh\gamma - \cos qL},$$

so the metasurface Green's function becomes

$$G_T(\mathbf{r}-\mathbf{r}') = \int\frac{d\mathbf{q}}{(2\pi)^2}\,g_T(\mathbf{q})\,e^{i\mathbf{q}(\mathbf{r}-\mathbf{r}')}\,\frac{\sinh\gamma}{\cosh\gamma - \cos q_x L_x}\cdot\frac{\sinh\gamma}{\cosh\gamma - \cos q_y L_y}.$$

In the limit $\gamma\ll 1$, each factor develops a narrow (delta-like) peak at $q_x L_x \approx 2\pi n_x$, $q_y L_y \approx 2\pi n_y$:

$$\frac{\sinh\gamma}{\cosh\gamma - \cos q_x L_x} \approx \frac{2\gamma}{\gamma^2 + (q_x L_x - 2\pi n_x)^2} \approx 2\pi\,\delta(q_x L_x - 2\pi n_x).$$

Carrying out the $\mathbf{q}$-integration via the delta functions yields

$$G_T(\mathbf{r}-\mathbf{r}') = \frac{1}{L_x L_y}\sum_{n_x,n_y}g_T(\mathbf{q}_{n_x n_y})\,e^{i\mathbf{q}_{n_x n_y}(\mathbf{r}-\mathbf{r}')},\qquad \mathbf{q}_{n_x n_y} = \frac{2\pi n_x}{L_x}\mathbf{e}_x + \frac{2\pi n_y}{L_y}\mathbf{e}_y.$$

The delta-function approximation fails at $n_x = n_y = 0$, however, since $g_T(\mathbf{q})\approx (\kappa|\mathbf{q}|)^{-1}$ as $|\mathbf{q}|\to 0$. The regularized formula, by contrast, lets us evaluate the small-$|\mathbf{q}|$ contribution explicitly:

$$\begin{aligned}
\delta G_T(\mathbf{r}-\mathbf{r}') &= \int\frac{dq_x\,dq_y}{(2\pi)^2}\,\frac{1}{\kappa\sqrt{q_x^2 + q_y^2}}\,e^{i\mathbf{q}(\mathbf{r}-\mathbf{r}')}\,\frac{2\gamma}{\gamma^2 + (q_x L_x)^2}\,\frac{2\gamma}{\gamma^2 + (q_y L_y)^2} \\
&\approx \int\frac{dq_y}{(2\pi)^2}\,\frac{1}{\kappa}\,\frac{2\gamma}{\gamma^2 + (q_y L_y)^2}\,\frac{4\,\operatorname{arccosh}\!\bigl(\gamma/(q_y L_y)\bigr)}{\sqrt{\gamma^2 + (q_y L_y)^2}} \\
&= \frac{1}{L_y\gamma\kappa}\int_{-\infty}^{+\infty}\frac{dt}{(2\pi)^2}\,\frac{8\,\operatorname{arccosh}(1/t)}{1+t^2} \approx \frac{0.8}{L_y\kappa}.
\end{aligned}$$

The upshot: the divergent piece of the sum — the $n_x = n_y = 0$ term — has no spatial dependence, and after regularization is simply a large constant added to the Green's function. This term corresponds to the overall temperature rise of the metasurface, but does not affect gradients or temperature differences within a unit cell.

## Average thermal conductivity of the SiO$_2$/Si composite substrate

We collect here the results for the mean temperature of square (side $L$) and circular (radius $R$) samples driven by a uniform power density.

For the disk:

$$p(\mathbf{r}) = \begin{cases} p_0, & r < R,\\ 0, & r > R.\end{cases}$$

For the square:

$$p(\mathbf{r}) = \begin{cases} p_0, & |x| < L/2,\ |y| < L/2,\\ 0, & \text{otherwise}.\end{cases}$$

$$\langle T_\square\rangle = \frac{Lp_0}{\kappa_1}\int_{-\infty}^{+\infty}\frac{d\xi\,d\eta}{(2\pi)^2}\,\frac{1}{\sqrt{\xi^2+\eta^2}}\left(\frac{\sin\xi/2}{\xi/2}\right)^{\!2}\left(\frac{\sin\eta/2}{\eta/2}\right)^{\!2}\,\frac{\kappa_1 + \kappa_2\tanh\!\bigl(d\sqrt{\xi^2+\eta^2}/L\bigr)}{\kappa_2 + \kappa_1\tanh\!\bigl(d\sqrt{\xi^2+\eta^2}/L\bigr)},$$

$$\langle T_\circ\rangle = \frac{Rp_0}{\kappa_1}\int_{-\infty}^{+\infty}\frac{d\tau}{2}\left(\frac{2J_1(\tau)}{\tau}\right)^{\!2}\,\frac{\kappa_1 + \kappa_2\tanh(\tau d/R)}{\kappa_2 + \kappa_1\tanh(\tau d/R)}.$$

Now introduce an effective thermal conductivity $\langle\kappa\rangle$ for the two-layer substrate, defined so that the single-layer problem with this effective conductivity reproduces the true two-layer solution. For the square we obtain

$$\frac{1}{\langle\kappa\rangle_\square}\,J_\square = \frac{1}{\kappa_1}\int_{-\infty}^{+\infty}\frac{d\xi\,d\eta}{(2\pi)^2}\,\frac{1}{\sqrt{\xi^2+\eta^2}}\left(\frac{\sin\xi/2}{\xi/2}\right)^{\!2}\left(\frac{\sin\eta/2}{\eta/2}\right)^{\!2}\,\frac{\kappa_1 + \kappa_2\tanh\!\bigl(d\sqrt{\xi^2+\eta^2}/L\bigr)}{\kappa_2 + \kappa_1\tanh\!\bigl(d\sqrt{\xi^2+\eta^2}/L\bigr)},$$

$$J_\square = \int_{-\infty}^{+\infty}\frac{d\xi\,d\eta}{(2\pi)^2}\,\frac{1}{\sqrt{\xi^2+\eta^2}}\left(\frac{\sin\xi/2}{\xi/2}\right)^{\!2}\left(\frac{\sin\eta/2}{\eta/2}\right)^{\!2} \approx 0.47,$$

and for the disk

$$\frac{J_\circ}{\langle\kappa_\circ\rangle} = \frac{1}{\kappa_1}\int_{-\infty}^{+\infty}\frac{d\tau}{2}\left(\frac{2J_1(\tau)}{\tau}\right)^{\!2}\,\frac{\kappa_1 + \kappa_2\tanh(\tau d/R)}{\kappa_2 + \kappa_1\tanh(\tau d/R)},$$

$$J_\circ = \int_{-\infty}^{+\infty}\frac{d\tau}{2}\left(\frac{2J_1(\tau)}{\tau}\right)^{\!2} = \frac{8}{3\pi}.$$

We plot the resulting effective conductivities for realistic SiO$_2$ thicknesses — 5 nm to 600 nm — and across a range of device sizes. Silicon is taken at 150 W/(m·K) and silicon oxide at 1 W/(m·K).

![Square sample](chart-square.png)

**(A)**

![Circular sample](chart-circle.png)

**(B)**

*Computed effective thermal conductivity of the Si/SiO2 composite substrate at various SiO2 thicknesses (color-coded) and for various graphene-device sizes. A — square sample, B — circular sample.*

The computed curves show three notable features:

- they interpolate correctly between the SiO$_2$ conductivity for small devices and the Si conductivity for large devices;
- when the device size becomes comparable to the SiO$_2$ thickness, the effective conductivity is still well below silicon's, because of the large contrast between the silicon and oxide conductivities;
- the effective conductivity reaches roughly half of silicon's only when the device size is of order $(\kappa_2/\kappa_1)\,d \approx 150\,d$. For instance, with 300 nm of oxide and a 100 μm square device the effective conductivity is 75 W/(m·K). Even for sizable devices, in other words, the oxide produces a noticeable reduction in the effective conductivity.

# Goal:

**Part 1.** Implement a **simpler scalar differential** problem — `∇·(κ(R) ∇T) = Q` on a 2D square with zero-flux BCs and spatially varying `κ` imitating metal/graphene patches — to develop the triangulation and basis-function infrastructure. 

**Part 2.** Migrate to the integral form once derived.

