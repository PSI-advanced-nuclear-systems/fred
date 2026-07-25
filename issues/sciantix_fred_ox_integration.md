# Feature: Incorporate SCIANTIX fission gas release model into FRED-OX

## Summary

The theory manual section on FRED-OX details the Waltar-Reynolds fission gas release (FGR) model, which partitions the fuel pellet into empirical zone-based regions (central restructured, unrestructured, as-fabricated) and assigns a release fraction to each. While this is computationally efficient, it is physically insufficient for scenarios where the microstructure evolution, grain boundary saturation, and intra-granular diffusion kinetics determine the release. Incorporating the open-source mechanistic code **SCIANTIX** as the FGR sub-model in FRED-OX would provide physically coherent modelling grounded in atomic-scale phenomena.

---

## Motivation: Why the Current Model is Insufficient

### Current model — Waltar-Reynolds

The current FGR model (see `fred_ox.tex`, §Fission Gas Generation and Release) is:

- **Zone-based and empirical.** Three radial zones are identified by temperature threshold ($T > 1773\,\text{K}$ for the CRSZ) and burnup threshold ($B > 3.5\,\text{MWd/kgU}$ for partial release from the unrestructured zone). The zone boundaries and release fractions are algebraic correlations, not ODEs tracking physical state.
- **No microstructure state.** The model carries no information about grain size, grain boundary gas inventory, intragranular bubble density, or grain boundary saturation fraction. All of this physics is folded into the empirical constants calibrated to a narrow set of EBR-II/Halden data.
- **Release is a fraction of generation, not a flux.** The global release ODE is $d\mu_{rel}/dt = \overline{\text{FGR}}_j \cdot d\mu_{gen}/dt$, which means release is instantaneously proportional to generation. This misses the time-delay between gas production and venting: gas must first diffuse to grain boundaries, saturate them, and then nucleate inter-granular porosity before release occurs.
- **No irradiation-enhanced diffusion.** The diffusion coefficient of Xe in UO2 has three contributions — thermal, irradiation-enhanced, and athermal (ballistic recoil) — that depend on fission rate and temperature separately. The Waltar-Reynolds model does not distinguish them.
- **Not predictive outside its calibration range.** The model is not physically predictive for high-burnup MOX ($B > 60\,\text{MWd/kgU}$), transient conditions, or pellet geometries significantly different from the calibration dataset. Extrapolation is unsupported.

### Phenomena the current model cannot represent

| Phenomenon | Physical mechanism | Status in current model |
|---|---|---|
| Intra-granular gas diffusion | Fick diffusion of Xe/Kr atoms in UO₂ lattice; $D_{eff}(T, \dot{F})$ | Absent |
| Grain boundary saturation | Gas sweeps to grain boundaries; saturation triggers inter-granular bubble formation | Absent |
| Inter-granular bubble growth/venting | Bubble coalescence and connection to open porosity; venting to plenum | Absent |
| Irradiation-enhanced diffusion | Additional athermal/ballistic diffusion at high fission rate | Absent |
| Burst release during power ramp | Rapid de-saturation and venting when grain boundary coverage threshold exceeded | Absent |
| Grain size evolution | Grain growth changes the diffusion path length and hence release lag | Absent |

---

## SCIANTIX: What it provides

[SCIANTIX](https://github.com/sciantix/sciantix) (Pizzocri et al., 2020, *Journal of Nuclear Materials*) is an open-source, standalone C++ fission gas behaviour module for nuclear fuel. It solves a system of coupled ODEs tracking:

1. **Intra-granular gas concentration** $c_g$ [at/m³]: gas produced by fission minus gas that reaches grain boundaries.
2. **Intra-granular bubble density** $N_b$ [1/m³] and mean radius $r_b$: nucleation and growth under fission rate and temperature.
3. **Grain boundary gas inventory** $c_{gb}$ [at/m²]: flux arriving from grain interior.
4. **Grain boundary bubble coverage fraction** $F_{cov}$ [-]: bubble growth until saturation threshold $F_{sat} \approx 0.5$.
5. **Released gas fraction** $f_{rel}$ [-]: once $F_{cov} \geq F_{sat}$, excess gas vents to open porosity.

The key ODEs are:

**Intra-granular diffusion** (modified Speight model):
$$\frac{dc_g}{dt} = \dot{F}\,y_{Xe} - \frac{D_{eff}}{a^2}\,S\,c_g$$

where $y_{Xe}$ is the fission yield of Xe+Kr, $a$ the grain radius, $S$ a shape factor, and $D_{eff}$ the effective diffusivity including thermal, irradiation-enhanced, and athermal contributions.

**Effective diffusivity** (Turnbull et al.):
$$D_{eff} = D_1\,e^{-Q_1/RT} + D_2\,\dot{F}^{0.5}\,e^{-Q_2/RT} + D_3\,\dot{F}$$

**Grain boundary gas balance:**
$$\frac{dc_{gb}}{dt} = \frac{D_{eff}}{a^2}\,S\,c_g\,\frac{4\pi a^2}{3\,N_{gb}} - \dot{c}_{rel}$$

where $N_{gb}$ is the number of grain boundary sites per unit volume and $\dot{c}_{rel}$ is the venting flux once grain boundary saturation is reached.

SCIANTIX is designed as a **standalone module** with a clean C++ interface: given the current temperature $T$, fission rate $\dot{F}$, time step $\Delta t$, and grain size $a$, it advances all internal state variables and returns the updated gas release fraction $f_{rel}$ and intra-granular swelling contribution. It manages its own internal sub-stepping for numerical stability.

---

## Proposed Integration Architecture

### What couples FGR to the rest of FRED-OX

FGR affects the FRED-OX solve through two channels:

1. **Gas pressure** `gpres`: $p_{gas} = (\mu_0 + \mu_{rel})\,R\,/\,V_t$ (Eq. `gpres`). `gpres` is an algebraic DAE variable that closes on the momentum balance.
2. **Gap conductance** `hgap`: the released Xe/Kr fraction reduces $k_{mix}$ via the geometric-mean mixing rule (Eq. `kmix`), directly affecting the radial temperature distribution.

Both of these are coupled into the SUNDIALS residual vector. Any change to `fgrel` must flow through them consistently.

### Should SCIANTIX ODEs enter the SUNDIALS residual vector?

This is the key architectural question. Two approaches are possible:

**Option A: Operator-split (recommended)**

Keep SCIANTIX outside the SUNDIALS y-vector. Advance SCIANTIX once per accepted SUNDIALS step, using the accepted temperature and fission-rate profile, and inject the resulting `fgrel` into the y-vector before the next step.

```
SUNDIALS BDF step (thermo-mechanical DAE + fggen ODE + fgrel algebraic)
    └─ converges at t_{n+1} with frozen fgrel from step n
After acceptance:
    SCIANTIX.advance(T[i], Fdot[i], dt)  // per layer, per node
    mu_rel_new = SCIANTIX.gasRelease()   // mol
    inject mu_rel_new → y-vector         // IDAReInit or direct write
    update gpres, gap conductance
```

This is identical in structure to how GRSIS bubble physics are handled in FRED-M-Na (post-step explicit advance, one-step lag into the next implicit solve). The one-step lag is acceptable because the timescale for FGR to change gas pressure significantly is much longer than a single SUNDIALS step.

**Advantages of Option A:**
- SCIANTIX manages its own internal sub-stepping; no interaction with SUNDIALS step-size control.
- No increase in DAE system dimension (SCIANTIX state — typically 10–15 scalars per node — does not enter the SUNDIALS y-vector).
- SUNDIALS convergence rate is unaffected; the Jacobian structure is unchanged.
- Clean separation: SCIANTIX is a black-box sub-model called once per accepted step.

**Option B: Fully implicit (SCIANTIX ODEs inside SUNDIALS)**

Add SCIANTIX internal variables ($c_g$, $c_{gb}$, $F_{cov}$ per node) to the SUNDIALS y-vector as additional ODE rows. The SCIANTIX residuals would then be evaluated inside `computeLayerResiduals` at every Newton iterate.

**Disadvantages of Option B:**
- The intra-granular diffusion timescale ($\sim a^2/D_{eff}$) can be much shorter than the thermo-mechanical timescale, requiring SUNDIALS to take very small steps — precisely what SCIANTIX's own sub-stepping is designed to absorb.
- Increases the DAE system by ~15 unknowns per radial node per layer (potentially hundreds of new variables per pin), degrading Newton convergence.
- The Jacobian coupling between $c_g$ and $T$ is non-trivial and would need finite-difference columns for all new variables, increasing the cost of the dense solve in each layer.

**Recommendation: Option A**, consistent with how GRSIS is handled in FRED-M-Na. The one-step lag on `fgrel` is well-justified given the slow evolution of cumulative gas release relative to the thermo-mechanical time step.

### Structural changes required

| Component | Change |
|---|---|
| `FredOxResiduals.cpp` | Replace algebraic FGR correlation call with identity (`fgrel` advances from post-step SCIANTIX output, not from current-step FGR fraction) |
| `FredOxSolver.cpp::afterAcceptedStep` | Add SCIANTIX advance call per layer; inject updated `fgrel` |
| `FredOxLayerState` | Add SCIANTIX internal state per layer (`SciantixNodeState` struct wrapping SCIANTIX's C++ objects) |
| `FredOxResiduals::syncAuxLayerState` | Copy SCIANTIX state alongside existing irradiation accumulators |
| `Makefile` | Add SCIANTIX as a dependency (header-only or static library) |
| Input deck / Python bindings | Expose grain size `a`, SCIANTIX model flags to user |

### Per-node vs per-layer resolution

SCIANTIX operates on a single point (single $T$, $\dot{F}$). Since FRED-OX resolves $n_f$ radial nodes per axial layer, SCIANTIX should be called once per radial node per axial layer — i.e., $n_z \times n_f$ times per accepted step. The released gas from all nodes is then summed to update the global $\mu_{rel}$. This is consistent with the radial temperature gradient being the primary driver of the spatially variable FGR.

---

## Time-step constraint imposed by the operator-split lag

The operator-split approach introduces a one-step lag: `fgrel` is frozen at its value from the end of step $n$ for the entire duration of SUNDIALS step $n+1$. The gas pressure used inside the Newton solve for step $n+1$ is therefore:

$$p_{gas}(t) \approx \frac{(\mu_0 + \mu_{rel}^{(n)})\,R}{V_t}, \qquad t \in [t_n,\, t_{n+1}]$$

If SUNDIALS takes a long adaptive step, the true `fgrel` has evolved during that interval and the error in the predicted plenum pressure is:

$$\Delta p_{error} \approx \frac{\dot{\mu}_{rel} \cdot \Delta t \cdot R}{V_t}$$

For this to remain below a fraction $\varepsilon$ of $p_{gas}$, the step must satisfy:

$$\Delta t \leq \varepsilon \cdot \frac{p_{gas}\,V_t}{\dot{\mu}_{rel}\,R} = \varepsilon \cdot \frac{\mu_0 + \mu_{rel}}{\dot{\mu}_{rel}} \equiv \varepsilon \cdot \tau_{FGR}$$

**The operator-split is only valid when $\Delta t \ll \tau_{FGR}$.** This is not guaranteed by SUNDIALS' own adaptive step-size control, which is driven by the thermo-mechanical DAE truncation error — not by the accuracy of the gas pressure.

### Why this matters in practice

During steady-state base irradiation at constant power, SUNDIALS may select $\Delta t$ of days or weeks (the thermal/mechanical state is near-stationary). $\tau_{FGR}$ can be comparably short, particularly late in life when accumulated gas inventory grows rapidly and grain boundary saturation drives burst release. In that regime, a long SUNDIALS step would use a stale `fgrel` to compute `gpres`, causing the gap conductance model to see the wrong He/Xe/Kr mixture and the cladding stress state to reflect the wrong internal pressure.

### Short-term fix: hard maximum step size

The simplest immediate safeguard is to impose a hard ceiling on the SUNDIALS step via `IDASetMaxStep` set once at solver initialisation, on the order of **1–7 days** (86 400–604 800 s). This bounds the lag error to at most one week's worth of FGR accumulation regardless of how slowly the thermo-mechanical state is evolving. It is a conservative, easy-to-audit measure that can be applied before SCIANTIX is integrated and requires no runtime estimation of $\dot{\mu}_{rel}$. The appropriate value should be guided by the expected FGR rate at the operating condition: for fast-reactor MOX at nominal power, a 1-day ceiling is a safe default.

### Long-term fix: FGR-aware maximum step via `IDASetMaxStep`

After each accepted step, estimate $\dot{\mu}_{rel}$ from the difference in SCIANTIX output between the current and previous step, then enforce:

$$h_{max}^{FGR} = \varepsilon \cdot \frac{\mu_0 + \mu_{rel}}{\dot{\mu}_{rel}}$$

with $\varepsilon \sim 0.05$–$0.1$. Pass this ceiling to SUNDIALS immediately after each accepted step:

```cpp
// afterAcceptedStep, after SCIANTIX advance:
const double dfgrel_dt = (fgrel_new - fgrel_old) / dt;
if (dfgrel_dt > 1e-20) {
    const double tau_fgr = (m_state.mu0 + fgrel_new) / dfgrel_dt;
    IDASetMaxStep(m_ida_mem, fgr_step_fraction * tau_fgr);
}
```

This is self-regulating: when FGR is negligible (early life), $\dot{\mu}_{rel} \approx 0$, $\tau_{FGR} \to \infty$, and the constraint is non-binding — SUNDIALS runs freely. During burst release or a power ramp, $\dot{\mu}_{rel}$ rises sharply, $\tau_{FGR}$ shortens, and SUNDIALS is automatically forced to take smaller steps — precisely when pressure accuracy matters most. The user's output cadence `dtout` is unaffected; SUNDIALS sub-steps internally.

### Why not sub-cycle SCIANTIX inside the SUNDIALS step?

A more accurate but invasive alternative is to call SCIANTIX multiple times within a single SUNDIALS step, updating `fgrel` and hence `gpres` at each sub-cycle. This requires either `IDAReInit` mid-step (expensive, discards BDF history) or treating `fgrel` as a non-DAE accumulator decoupled from the y-vector mid-step. Neither is recommended: the `IDASetMaxStep` approach achieves the same accuracy bound with far less complexity.

---

## References

- Pizzocri, D. et al. (2020). *SCIANTIX: A new open source multi-scale code for fission gas behaviour modelling in LWR fuel*. Journal of Nuclear Materials, 532, 152047.
- Waltar, A.E. & Reynolds, A.B. (1981). *Fast Breeder Reactors*. Pergamon Press. (current FRED-OX FGR model basis)
- Turnbull, J.A. et al. (1982). *The diffusion coefficients of gaseous and volatile species during the irradiation of uranium dioxide*. Journal of Nuclear Materials, 107, 168–184.
- Speight, M.V. (1969). *A calculation on the migration of fission gas in nuclear fuel*. Nuclear Science and Engineering, 37, 180–185.
- SCIANTIX repository: https://github.com/sciantix/sciantix
