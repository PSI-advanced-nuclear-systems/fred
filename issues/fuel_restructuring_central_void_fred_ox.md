# Feature: Implementation of fuel restructuring model and central void formation within FRED-OX

## Summary

At the elevated temperatures typical of sodium fast reactor (SFR) oxide fuel pins, the steep
radial temperature gradient drives thermally-activated pore migration and grain growth, partitioning
the fuel pellet into distinct radial zones with different microstructure and density.  The innermost
zone develops a **central void** as lenticular pores migrate up the temperature gradient and
coalesce near the pellet axis.  FRED-OX currently applies the Waltar-Reynolds FGR zone model but
does not simulate the physical redistribution of fuel and void formation that establishes those
zones.  Implementing a fuel restructuring model would correct the conductivity field and power
density used in the heat conduction solve, provide self-consistent zone boundaries for the FGR
model, and update the pellet inner boundary condition once a central void forms.

---

## Motivation: Phenomena Not Captured by the Current Model

FRED-OX applies empirical temperature threshold criteria to identify the CRSZ and equiaxed zone
for the Waltar-Reynolds FGR model, but:

- The fuel geometry remains a **solid cylinder** ($r_{fi} = 0$) throughout irradiation.
- The thermal conductivity porosity correction uses the **as-fabricated** porosity uniformly
  across all radial nodes; no redistribution of porosity between zones is modelled.
- The volumetric power density is **radially uniform** within the pellet, ignoring the density
  variation between restructured zones.
- The FGR zone boundaries are driven by the current temperature profile but are not tied to
  a density or microstructure field that evolves with burnup.

---

## Physical Model — Walter & Reynolds (1981), §9-2B

### Nomenclature

Figure 9-4 of Walter & Reynolds identifies four concentric regions after restructuring (from
centre outward), with radii $R_0, R_1, R_2, R_3 = R_f$ and temperatures $T_0, T_1, T_2, T_3 = T_f$:

| Label | Region | Inner radius | Outer radius | Density | Inner T | Outer T |
|---|---|---|---|---|---|---|
| — | Central void | $0$ | $R_0$ | $0$ | — | $T_0 \approx T_{solidus}$ |
| Region 1 | Columnar grain | $R_0$ | $R_1$ | $\rho_1$ | $T_0$ | $T_1 = T_{cg}$ |
| Region 2 | Equiaxed grain | $R_1$ | $R_2$ | $\rho_2$ | $T_1$ | $T_2 = T_{eg}$ |
| Region 3 | As-fabricated | $R_2$ | $R_f$ | $\rho_3$ | $T_2$ | $T_f$ |

Temperature thresholds and zone densities (Walter & Reynolds Fig. 9-5, MOX, $\rho_3 = 85\%$ TD):

| Parameter | Value | Notes |
|---|---|---|
| $T_{cg}$ (columnar grain threshold) | $1800\,°\text{C}$ (2073 K) | outer boundary of Region 1 |
| $T_{eg}$ (equiaxed grain threshold) | $1600\,°\text{C}$ (1873 K) | outer boundary of Region 2 |
| $\rho_1$ | $0.98\,\rho_{TD}$ | pores driven to void; columnar zone near full density |
| $\rho_2$ | $0.95\,\rho_{TD}$ | grain growth closes some pores |
| $\rho_3$ | as-fabricated | typically $0.85$–$0.93\,\rho_{TD}$ |

### §1  Mass balance (Eq. 9-11)

Regions 1 and 2 contain material that was originally at density $\rho_3$ inside radius $R_2$
before restructuring.  Conserving mass per unit length inside $R_2$:

$$\boxed{\rho_1\!\left(R_1^2 - R_0^2\right) + \rho_2\!\left(R_2^2 - R_1^2\right) = \rho_3 R_2^2}$$

This gives the central void radius $R_0$ once $R_1$ (columnar outer boundary) and $R_2$
(equiaxed outer boundary) are known from the temperature profile:

$$R_0 = \sqrt{R_1^2 - \frac{\rho_3 R_2^2 - \rho_2(R_2^2 - R_1^2)}{\rho_1}}
      = \sqrt{\frac{\rho_1 R_1^2 + \rho_2(R_2^2 - R_1^2) - \rho_3 R_2^2}{\rho_1}}\cdot\frac{1}{1}$$

In the limiting case where the equiaxed zone is negligible ($R_1 \approx R_2$):

$$R_0 \approx R_1\sqrt{1 - \frac{\rho_3}{\rho_1}}$$

### §2  Volumetric heat generation rates (Eqs. 9-9a, 9-9b, 9-9c)

The linear power $X$ [W/m] is conserved across restructuring.  Since the fission rate per unit
volume scales with local atom (fissile) number density, hence with local mass density, the
heat generation rate in each zone is:

$$Q_1 = \frac{X}{\pi R_f^2}\,\frac{\rho_1}{\rho_3}, \qquad
  Q_2 = \frac{X}{\pi R_f^2}\,\frac{\rho_2}{\rho_3}, \qquad
  Q_3 = \frac{X}{\pi R_f^2}$$

These are self-consistent with mass conservation: substituting into
$X = Q_1\pi(R_1^2-R_0^2) + Q_2\pi(R_2^2-R_1^2) + Q_3\pi(R_f^2-R_2^2)$
and applying Eq. (9-11) recovers $X = X$.

Because $\rho_1 > \rho_2 > \rho_3$, the innermost (columnar) zone generates more heat per unit
volume than the unrestructured outer annulus.  This is partly offset by the higher conductivity
in the denser inner zones.

### §3  Temperature distribution via the conductivity integral (Eqs. 9-12, 9-13, 9-14)

Walter & Reynolds solve the heat conduction equation analytically for each zone without any
spatial discretisation.  The approach integrates the steady-state cylindrical energy balance
across each zone to express the temperature drop as a **conductivity integral** $\int k\,dT$.

For a zone bounded by $[r_a, r_b]$ with volumetric power $Q$ and a total outward linear power
$\chi(r)$ at radius $r$, the heat equation gives:

$$-k\frac{dT}{dr} = \frac{\chi(r)}{2\pi r}$$

with $\chi(r) = $ (total power generated from void surface to $r$, per unit length).

Integrating across each zone from outer boundary (lower T, known) to inner boundary (higher T):

**Region 3 — unrestructured ($R_2 \leq r \leq R_f$)** (Eq. 9-12):

$$\int_{T_f}^{T_2} k_3\,dT = \frac{X}{4\pi}\left[1 - \left(\frac{R_2}{R_f}\right)^2\right]$$

This is the same as the solid-cylinder conductivity integral formula because $Q_3 = X/(\pi R_f^2)$
is independent of the inner zones' restructuring.

**Region 2 — equiaxed ($R_1 \leq r \leq R_2$)** (Eq. 9-13):

Using mass conservation (Eq. 9-11) to eliminate $R_0$, the integral across the equiaxed zone is:

$$\int_{T_2}^{T_1} k_2\,dT = \frac{X}{4\pi\,\rho_3 R_f^2}
  \left[\rho_2\!\left(R_2^2 - R_1^2\right)
       + 2\!\left(\rho_3 - \rho_2\right)R_2^2\ln\frac{R_2}{R_1}\right]$$

or equivalently:

$$= \frac{X}{4\pi}\left[\frac{\rho_2}{\rho_3}\frac{R_2^2-R_1^2}{R_f^2}
  + 2\!\left(1 - \frac{\rho_2}{\rho_3}\right)\frac{R_2^2}{R_f^2}\ln\frac{R_2}{R_1}\right]$$

Since $\rho_2 > \rho_3$, the factor $(1 - \rho_2/\rho_3) < 0$, and the logarithm term reduces
the temperature drop relative to an equivalent solid annulus — the equiaxed zone "benefits"
from the extra heat delivered from the denser inner zone.

**Region 1 — columnar ($R_0 \leq r \leq R_1$)** (Eq. 9-14):

$$\int_{T_1}^{T_0} k_1\,dT = \frac{X}{4\pi}\frac{\rho_1}{\rho_3}\frac{1}{R_f^2}
  \left[\left(R_1^2 - R_0^2\right) - 2R_0^2\ln\frac{R_1}{R_0}\right]$$

The void surface temperature $T_0$ can be found from this equation if $T_1 = T_{cg}$ and $R_0$
are known, or it can be prescribed as $T_0 = T_{solidus}$ to close the system.

### §4  Applying the W&R model (solving from outside in)

Given the linear power $X$ and fuel surface temperature $T_f$ (from gap conductance):

1. **Find $R_2$**: solve Eq. (9-12) for $R_2$ such that $\int_{T_f}^{T_{eg}} k_3\,dT = \frac{X}{4\pi}[1-(R_2/R_f)^2]$ — a
   one-dimensional root-find in $R_2$.

2. **Find $R_1$**: with $R_2$ known, solve Eq. (9-13) for $R_1$ such that the integral evaluates
   to $T_1 = T_{cg}$.

3. **Find $R_0$**: apply Eq. (9-11) (mass balance) to compute $R_0$ from $R_1$ and $R_2$.

4. **Find $T_0$**: evaluate Eq. (9-14) to get the void surface temperature.

The entire temperature profile is then specified by four radii and four temperatures — no radial
mesh required.

**The example in W&R Fig. 9-5** ($X = 50$ kW/m, $T_c = 1000°C$, $\rho_3 = 85\%$ TD) shows
a ~400°C reduction in peak fuel temperature after restructuring compared to the uniform-density
solid cylinder, due to the higher conductivity in the denser inner zones.

---

## Key Architectural Mismatch: Analytical Model vs FRED-OX Discretised Solver

The W&R model is **fundamentally continuous**: zone boundaries $R_0, R_1, R_2$ are
algebraic unknowns determined by temperature isotherms and mass conservation.  The temperature
profile is derived analytically once the boundaries are known — no spatial mesh is used.

FRED-OX solves for temperatures on a **fixed radial mesh** of $n_f$ nodes using a finite-volume
residual system.  Zone boundaries in general do not coincide with node positions.  This creates
a two-sided mismatch:

| Aspect | W&R analytical | FRED-OX discretised |
|---|---|---|
| Temperature unknowns | Zone boundary temperatures $T_1, T_2, T_0$ (3 scalars per layer) | Node temperatures $T[0\ldots n_f-1]$ ($n_f$ unknowns per layer) |
| Zone boundaries | Primary unknowns: $R_1, R_2, R_0$ solved from isotherms | Derived: found by isotherm interpolation in $T[i]$ after each step |
| Zone assignment | Exact by construction | Approximate: each node assigned to one zone; boundary straddles a cell |
| Conductivity | Piecewise constant per zone, exact | Step function in $i$; O($\Delta r$) boundary error |
| Power density | Piecewise constant per zone | Step function in $i$ |
| Inner BC | Moves to $R_0$ continuously | Innermost active node is always a fixed grid point |

### Two implementation approaches

The options described below reduce to two distinct strategies.  The choice should be made
once and committed to in the new app (see §"Recommended app: fred-ox-fr"):

| | Approach 1 — Fixed mesh, void node treatment | Approach 2 — Moving mesh (redistribution) |
|---|---|---|
| Nodes | All $n_f$ active; void nodes assigned $Q=0$, $\rho\to0$ | All $n_f$ active; mesh always spans $[R_0, R_f]$ |
| Inner BC (thermal) | Zero-flux at $r=0$ (wrong when void forms) | Zero-flux at `rad0[0]`$=R_0$ (correct automatically) |
| Inner BC (mechanical) | Symmetry at $r=0$ (wrong when void forms) | $\sigma_{rr}(R_0) = -p_{gas}$ (single BC change needed) |
| T–stress consistency | Broken: void nodes have high $T$ but no material | Consistent: all nodes are physical fuel |
| Remapping cost | None | Interpolation of state vars when $R_0$ advances |
| Recommended for | Phase 1 (zone coefficients, no void) | Phase 2 (full void geometry) |

### Proposed detail per approach

#### Option A — Post-step isotherm detection (recommended first step)

After each accepted SUNDIALS step, scan the node temperature array to find zone boundaries
by linear interpolation:

```
r_cg = interpolate(r, T, T_cg)   // find r where T[i] crosses T_cg
r_eg = interpolate(r, T, T_eg)   // find r where T[i] crosses T_eg
```

Each node $i$ is assigned to a zone by comparing $r[i]$ against $r_{cg}$ and $r_{eg}$.
The central void radius $R_0$ is then computed from mass conservation (Eq. 9-11).
Node densities `rho_zone[i]` and `k_fuel_per_node[i]` are updated accordingly, using the
Path B infrastructure already in `FuelRodState.hpp`.  These updated coefficients are frozen
for the next SUNDIALS Newton solve, just as `k_fuel_per_node` is in FRED-M-Na.

**Accuracy**: zone boundary is located to within one cell width $\Delta r = R_f/n_f$.  With
$n_f = 10$ nodes and $R_f \approx 3$ mm, $\Delta r \approx 0.3$ mm — comparable to the typical
central void radius at moderate linear power.

#### Option B — Sub-cell volume fractions at straddling nodes (higher accuracy, O(Δr²))

For nodes whose cell straddles a zone boundary, compute the fraction of the cell area inside
each zone using the interpolated boundary position:

```
f_cg[i] = fraction of cell i inside columnar zone  [0,1]
f_eg[i] = fraction of cell i inside equiaxed zone  [0,1]
f_0[i]  = fraction of cell i inside as-fabricated zone [0,1]
```

Assign a volume-fraction-weighted density and conductivity:

$$\rho_i = f_{cg}^{(i)}\rho_1 + f_{eg}^{(i)}\rho_2 + f_0^{(i)}\rho_3$$

This eliminates the staircase artifact at zone boundaries without requiring additional nodes
or a moving mesh.  Recommended for cells near $r_{cg}$ and $r_{eg}$; nodes far from boundaries
use the binary assignment of Option A.

#### Option C — W&R analytical solution as a zone-boundary predictor

Implement the W&R conductivity integral equations (9-12), (9-13) as a fast standalone
procedure that takes the current $X$ and $T_f$ and solves for $R_1, R_2, R_0$ by numerical
root-finding (bisection in each zone working outside-in).  The conductivity integral
$\int_{T_a}^{T_b} k(T)\,dT$ is evaluated by adaptive quadrature over the existing
`thermalConductivity(T)` function.

These analytical zone radii are then used in place of the interpolation in Option A to assign
nodes, giving accurate zone boundaries regardless of $\Delta r$.

**When to use**: when $n_f$ is small (e.g. $n_f = 5$–$7$), the O($\Delta r$) error of Option A
can misplace the columnar boundary by a significant fraction of the void radius.  Option C
gives sub-cell accuracy at negligible computational cost ($\sim 3$ root-finds per axial layer
per accepted step).

**Limitation**: the W&R model assumes steady-state conduction and a power profile proportional
to density.  During transients with rapidly changing $X$, the instantaneous $T_f$ and $X$
may not represent the actual thermal state inside the pellet; Option A tracking the actual
SUNDIALS temperature field is more appropriate in that case.

#### Option D — Mesh contraction (recommended for Phase 2)

When $R_0 > 0$, update the reference geometry so that the $n_f$ fuel nodes always span the
active fuel domain $[R_0,\, R_f]$:

$$r_0^{(i)} = R_0 + i\,\frac{R_f - R_0}{n_f - 1}, \qquad i = 0,\ldots,n_f-1$$

Node spacing contracts to $\Delta r(t) = (R_f - R_0(t))/(n_f-1)$ as the void grows.  All
$n_f$ nodes remain in the solve as active fuel material — there are no void nodes.

**Why this resolves the thermal-mechanical consistency problem** (see §"Consistency problem"
below): the thermal inner BC (zero-flux at `rad0[0]`) automatically relocates to the void
surface $r = R_0$ because that is now the innermost node.  The stress-strain inner BC changes
from the solid-pellet symmetry condition ($u_r = 0$ at $r = 0$) to the annular traction
condition ($\sigma_{rr}(R_0) = -p_{gas}$), but this is a single BC change in the stress-strain
residual, not an index-mapping restructuring.  Since both the thermal and mechanical domains
now start at $R_0$, the two solves remain coupled consistently.

**State variable remapping**: when `rad0[i]` changes, each node represents a different material
point.  All nodal state must be interpolated to the new positions:
- $T[i]$, $\varepsilon_{th}[i]$, $\varepsilon_{h}[i]$ — smooth fields, linear interpolation adequate.
- Creep and swelling strain history $\varepsilon_c[i]$, $\varepsilon_s[i]$ — Lagrangian history
  variables attached to material points; spatial interpolation introduces physical error.

**Energy balance during remapping**: naively interpolating $T[i]$ to new node positions without
accounting for the change in ring volume would violate thermal energy conservation.  The thermal
energy per node is $E_i = \rho_i c_p T_i A_i \Delta z$, and both $\rho_i$ and $A_i$ change
when the mesh moves.  However, the mass conservation equation (Eq. 9-11) guarantees:

$$\sum_i^{\text{new}} \rho_i^{new} A_i^{new} \approx \sum_i^{\text{old}} \rho_i^{old} A_i^{old}$$

That is, the density increase in the restructured zones exactly cancels the reduction in ring
area due to mesh contraction.  If the temperature field is approximately constant during the
small remap event ($\Delta T \ll T$), then $\sum \rho A T$ is conserved to first order — the
density update itself carries the energy balance.  A formal energy-conservative interpolation
(weight the $T$ field by $\rho c_p A$ on both meshes before and after) is only needed if $R_0$
jumps by more than $\sim \Delta r$ in a single step, which is controlled by the time-step
constraint discussed in §"Time-step constraint" below.

The error is negligible in practice because $R_0$ grows slowly: a practical guard is to remap
only when $R_0$ has advanced by more than $\Delta r / 4$ since the last remap, so the
interpolation perturbation per remap is small.

**Geometry updates required**: `FuelRodGeometry::rad0[i]`, `drf0`, `area0[i]`, and `rad_half[i]`
must all be recomputed after each remap.  These are currently set once at initialisation in
`FuelRodGeometry`; making them mutable (or storing them in the layer state) is the primary
structural change.

#### What Options A/B/C do and do not reproduce from W&R Fig. 9-5

Restructuring affects the FRED-OX solve through two paths:

**Path 1 — Direct thermal (dominant, captured by A/B/C for the non-void zones):**
zone assignment → `k_fuel_per_node[i]`, `Q_node[i]` per zone → `HeatConduction::computeResiduals` → $T[i]$

The ~400 °C drop in peak fuel temperature shown in Fig. 9-5 is almost entirely from this path:
the columnar zone at $\rho_1 = 98\%$ TD has substantially higher conductivity than as-fabricated
85% TD material.  Updating `k_fuel_per_node` for the non-void zones (columnar and equiaxed)
captures this dominant effect.

**Path 2 — Indirect mechanical (secondary):**
$R_0$ → annular pellet mechanics → pellet outer radius → gap width → $h_{gap}$ → $T_f$ → $T[i]$

Without the annular stress model, the gap mechanics are wrong, which errors $T_f$, but this is
a smaller correction than Path 1.

#### The thermal-mechanical consistency problem with void nodes

If void nodes ($r < R_0$) are left in the solve and assigned $Q = 0$ (no heat source) to suppress
fission heating:

- **Thermal**: the steady-state solution in the void region becomes $\nabla^2 T = 0$ with
  zero-flux at $r = 0$, giving $T_{void} \approx T(R_0) \approx T_{solidus} \approx 1800$–$2000\,°C$.
  This is approximately correct for an isothermal gas void.

- **Stress-strain**: the stress-strain equations evaluate thermal expansion as
  $\varepsilon_{th}[i] = \alpha(T[i])\,(T[i] - T_{ref})$ at void nodes.  With $T_{void} \approx 2000\,°C$,
  this gives a large thermal strain as if solid hot fuel is present — the model does not know
  the material has been removed.  This incorrect expansion propagates into the pellet outer
  radius, changing the gap and $T_f$.

Setting $\rho_{void} \to 0$ removes the thermal mass ($\rho c_p \dot{T} \to 0$) but does **not**
fix the stress-strain, because the elastic modulus $E(T)$ and expansion coefficient $\alpha(T)$
in FRED are evaluated from temperature only, not density.  Stiffness matrix entries for void
nodes remain non-zero and produce unphysical thermal stresses.

**Conclusion**: the thermal and mechanical treatments of the void must be changed together.
Implementing Option D (thermal inner BC) without the annular stress model — or vice versa —
creates a different inconsistency.  The two changes are a mandatory pair.

#### Safe Phase 1 scope: zone coefficients only, no void

To avoid the consistency problem, Phase 1 should update zone-resolved coefficients for the
**non-void zones only** and make no attempt to model the central void:

- Detect $r_{cg}$ and $r_{eg}$ by isotherm interpolation of the current $T[i]$ (Option A).
- Update `k_fuel_per_node[i]` and `Q_node[i]` for nodes with $r[i] \geq r_{cg}$ (equiaxed
  and as-fabricated zones).  Nodes inside $r_{cg}$ (columnar zone and future void) receive
  the columnar-zone conductivity.  **Do not attempt to zero out any nodes.**
- Compute $r_{cv}$ from mass conservation and store it as a diagnostic, but do not use it to
  modify the thermal or mechanical solve.
- Use $r_{cg}$ and $r_{eg}$ as CRSZ and equiaxed boundaries for the FGR model.

This is self-consistent: all nodes remain in both the thermal and mechanical solves with valid
material properties.  The dominant conductivity improvement is captured.  The central void
geometry is deferred entirely to Phase 2.

#### Recommended phasing

| Phase | Content | What it fixes |
|---|---|---|
| 1a | Option A isotherm scan; zone-resolved $k_i$, $Q_i$ for non-void zones only; FGR zone sync; $R_0$ stored as diagnostic | Dominant ~400 °C peak $T$ reduction; self-consistent FGR zones; no consistency risk |
| 1b | Option B (sub-cell fractions) | Removes staircase artefact at zone boundaries |
| 1c | Option C (W&R predictor) | Accurate $r_{cg}$, $r_{eg}$ with coarse $n_f$ |
| 2  | Option D (mesh contraction to $[R_0, R_f]$) + annular inner BC for stress-strain ($\sigma_r(R_0) = -p_{gas}$), implemented together | Consistent thermal and mechanical treatment of void; nodes compress naturally as void grows; no void-node inconsistency |
| 3  | Pore migration ODE | Kinetic $R_0(t)$ during transients |

---

## Recommended App: fred-ox-fr

The restructuring model involves changes to FRED-OX that are too invasive to carry inside the
existing `fred_ox` app without risk of breaking it:

- `FuelRodGeometry` changes from immutable to mutable (moving mesh).
- `HeatConduction::computeResiduals` receives per-node power density and conductivity (Path B).
- The stress-strain inner BC changes from symmetry to traction at a runtime-determined radius.
- New layer state fields (`r_cv`, `r_cg`, `r_eg`, `rho_zone`) are added.

The recommended approach is a **new app `fred-ox-fr`** (`fr` = fuel restructuring), living at
`src/apps/fred_ox_fr/`.  It inherits the full FRED-OX physics stack (same platform residuals,
same gap model, same FGR correlations, same SUNDIALS setup) and adds a `FuelRestructuring`
module that runs in `afterAcceptedStep`.  FRED-OX itself remains unchanged and continues to
be the validated, stable reference for non-restructuring cases.

```
src/apps/
  fred_ox/          ← unchanged; no restructuring
  fred_ox_fr/       ← new; inherits fred_ox residuals, adds restructuring layer
    FredOxFrSolver.cpp        # afterAcceptedStep calls FuelRestructuring::update(...)
    FredOxFrLayerState.hpp    # extends AxialLayerState with r_cv, r_cg, r_eg, rho_zone
    FuelRestructuring.hpp/.cpp # zone detection + either void-node or moving-mesh approach
```

### Approach selection for fred-ox-fr

Only one approach should be implemented in fred-ox-fr.  Based on the analysis above:

- **Phase 1**: implement Approach 1 (fixed mesh, zone coefficients only, no void).  This is
  the safe first step: zone-resolved $k$ and $Q$, FGR zone sync, $R_0$ stored as diagnostic.
  No void geometry is modified.

- **Phase 2**: replace with Approach 2 (moving mesh).  The mesh contraction is added together
  with the annular inner BC for stress-strain.  Phase 1 code is replaced, not extended, because
  the two approaches are architecturally incompatible (Approach 1 keeps void nodes active;
  Approach 2 removes them from the domain).

---

## Post-step Coupling and Time-step Constraint

### Restructuring runs outside the SUNDIALS solve

The restructuring model (zone boundary detection, $R_0$ computation, coefficient update, and
mesh remap in Phase 2) runs entirely in `afterAcceptedStep`, after SUNDIALS has converged.
The updated coefficients — `k_fuel_per_node`, `rho_zone`, `rad0` — are then **frozen** for
the duration of the next SUNDIALS Newton solve, identical in structure to how FRED-M-Na handles
`k_fuel_per_node` and GRSIS state.

This introduces a one-step lag: if SUNDIALS takes a long adaptive step, the zone boundaries
and void radius used in that step may be stale relative to the temperature profile at the end
of the step.

### Error from the one-step lag

The error in the peak fuel temperature due to a stale $R_0$ is approximately:

$$\Delta T_{peak} \approx \left.\frac{\partial T_{peak}}{\partial R_0}\right\vert_{X,T_f} \cdot \dot{R}_0 \cdot \Delta t$$

where $\dot{R}_0 = dR_0/dt$ is the void growth rate.  To keep $\Delta T_{peak}$ below a
fraction $\varepsilon$ of the peak temperature:

$$\Delta t \leq \varepsilon \cdot \frac{R_0}{\dot{R}_0} \equiv \varepsilon \cdot \tau_{void}$$

The same argument applies to the zone radii $R_1$ and $R_2$: if these move significantly
during a step, the zone-resolved conductivity and power density used in that step are wrong.

### Time-step ceiling via IDASetMaxStep

After each accepted step, estimate the void growth rate:

```cpp
const double dR0_dt = (R0_new - R0_old) / dt;
if (dR0_dt > 1e-15) {
    const double tau_void = R0_new / dR0_dt;
    IDASetMaxStep(m_ida_mem, restructure_step_fraction * tau_void);
}
```

with `restructure_step_fraction` $\sim 0.05$–$0.1$.  This is self-regulating: when
restructuring is complete and $\dot{R}_0 \approx 0$, $\tau_{void} \to \infty$ and the
constraint is non-binding.  During the initial restructuring transient (when the void first
forms and grows rapidly), $\tau_{void}$ is short and SUNDIALS is forced to take small steps.

### Short-term hard ceiling

During the initial restructuring period (first occurrence of $R_0 > 0$, or after a power
increase that drives a step change in zone boundaries), a hard ceiling of **1–7 days** via
`IDASetMaxStep` is a safe, simple bound before the adaptive estimate above is available.
This mirrors the approach recommended in the SCIANTIX integration issue for FGR lag.

---

## Impact on FRED-OX

### 1.  Modified inner boundary condition

Before restructuring: zero-flux symmetry at $r = 0$, enforced by the stencil in
`HeatConduction.cpp`.

After void formation: when $R_0 > 0$, the inner boundary moves to $r = R_0$.  The void surface
temperature is approximately $T_{solidus}$ (since the columnar zone forms at the melting point).
Until Option D is implemented, a conservative approximation is to retain zero-flux at the
innermost node and accept a small error in the inner temperature.

### 2.  Mesh contraction and annular inner BC (Phase 2, must ship together)

When $R_0$ first exceeds $\Delta r / 4$, remap all fuel nodes onto $[R_0, R_f]$, update
`FuelRodGeometry::rad0`, `drf0`, `area0`, `rad_half`, and interpolate the current node
state to the new positions.  Simultaneously change the stress-strain inner boundary condition
from symmetry ($u_r = 0$ at $r = 0$) to the annular traction condition
($\sigma_{rr}(R_0) = -p_{gas}$).  These two changes must be made together: the thermal solve
inner BC relocates automatically because `rad0[0] = R_0`, while the mechanical BC requires
an explicit change in the stress-strain residual.

### 3.  Per-node conductivity and density (Path B)

The restructuring model requires per-node conductivity $k_{fuel\_per\_node}[i]$, which FRED-OX
currently does not populate (it uses Path A with scalar `k_irr_factor`).  Enabling Path B for
FRED-OX is architecturally straightforward: the infrastructure is already in `FuelRodState.hpp`
and `HeatConduction.cpp`.  The restructuring module fills `k_fuel_per_node[i]` and `rho_zone[i]`
in `afterAcceptedStep`; from that point the Newton solve uses frozen per-node conductivity for
the next step.

### 3.  Self-consistent FGR zone boundaries

The Waltar-Reynolds FGR model uses the CRSZ outer boundary ($T > 2000$ K) and equiaxed
boundary.  Once restructuring is tracked, these boundaries should be taken directly from
$R_1$ and $R_2$ (the zone radii computed by the restructuring model) rather than re-derived
from in-line temperature threshold tests.

### 4.  Stress-strain (deferred)

The central void converts the pellet from a solid to an annular cylinder.  The current
`FuelStressStrain` model assumes $r_{fi} = 0$.  Correcting this requires an annular fuel
stress formulation with a traction-free inner surface at $r = R_0$, deferred to Phase 2.

---

## Implementation Plan

### Phase 1 — fred-ox-fr, Approach 1 (fixed mesh, zone coefficients)

| Component | Change |
|---|---|
| New app `src/apps/fred_ox_fr/` | Inherits fred_ox residuals; adds `FuelRestructuring` post-step call |
| New `FuelRestructuring.hpp/.cpp` | Isotherm scan (Option A); sub-cell fractions (Option B); W&R predictor (Option C); mass conservation for $R_0$ (diagnostic only) |
| `FredOxFrLayerState` | Extends `AxialLayerState` with `r_cv`, `r_cg`, `r_eg`, `rho_zone[nf]` |
| `FredOxFrSolver::afterAcceptedStep` | Call `FuelRestructuring::update(...)`; fill `k_fuel_per_node` (Path B); update per-node `qqv`; update `IDASetMaxStep` from $\tau_{void}$ |
| `FredOxFrResiduals` (power density) | Use `rho_zone[i] / rho_avg * qqv_layer` per node |
| `FredOxFrResiduals` (FGR) | Take CRSZ / equiaxed boundaries from `r_cg`, `r_eg` in layer state |
| Verification | Reproduce W&R Fig. 9-5 temperature profile; check convergence with $n_f$; compare $R_0$ diagnostic against analytical mass-balance result |

### Phase 2 — fred-ox-fr, Approach 2 (moving mesh, full void geometry)

| Component | Change |
|---|---|
| `FuelRodGeometry` | Make `rad0`, `drf0`, `area0`, `rad_half` mutable; add `updateMesh(R0)` method |
| `FuelRestructuring` | After computing $R_0$: call `geom.updateMesh(R0)` when $R_0 > R_0^{prev} + \Delta r/4$; conservative $T[i]$ remap (weight by $\rho c_p A$) |
| Stress-strain residual | Change inner BC from symmetry at $r=0$ to $\sigma_{rr}(r_{fi}) = -p_{gas}$ when `rad0[0]` $> 0$ |
| `fred_ox` | Unchanged — remains the stable non-restructuring reference |

---

## References

- Walter, A.E. & Reynolds, A.B. (1981). *Fast Breeder Reactors*. Pergamon Press.
  §9-2B: Effects of Restructuring, pp. 318–321.  Source of Eqs. (9-9)–(9-14), Fig. 9-4
  (zone geometry), Fig. 9-5 ($T_{cg} = 1800°C$, $T_{eg} = 1600°C$, $\rho_1 = 98\%$,
  $\rho_2 = 95\%$ TD).

- Karahan, A. (2009). *Modelling of Thermo-Mechanical and Irradiation Behaviour of Metallic
  and Oxide Fuels for Sodium Fast Reactors*. PhD Thesis, MIT.

- Timpano, D. (2024). Master's Thesis, EPFL.
  [Local PDF](../../timpano_fred/01.new/04.MANUAL/Timpano_MasterThesis_2024_v6.pdf)
