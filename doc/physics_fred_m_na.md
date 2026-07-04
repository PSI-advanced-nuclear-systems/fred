# FRED-M-Na: Physics Documentation

## Overview

FRED-M-Na is a FRED platform application for simulating U-Pu-Zr metallic fuel rods in
sodium-cooled fast reactors. It uses the same core thermal conduction and stress-strain
equations as FRED-ROD, and adds irradiation physics specific to metallic fuel:

1. **Sodium infiltration into fuel pores** — affects thermal conductivity
2. **Zirconium redistribution** — drives Zr from hot centre to cooler periphery
3. **Fission gas release (FGR)** — empirical model for upuzr (Karahan 2009)
4. **Fuel swelling** — solid fission product (FP) contribution
5. **Cladding wastage** — lanthanide diffusion into HT-9 cladding

All correlations are ported from the legacy Fortran FRED_M_OCT24.SRC
(Timpano/EPFL 2018–2024, building on Mikityuk/EPFL).

---

## 1. Fuel Material: U-Pu-Zr

### 1.1 Phase Diagram  (Fphase.for — Karahan 2007)

The ternary U-Pu-Zr phase diagram is approximated by linear interpolation between
the binary U-Zr boundaries (0% Pu) and U-19%Pu-Zr boundaries (19% Pu).

Six phase regions are identified:
- **alpha** (BCC + orthorhombic mixture, low T, low Zr)
- **alpha+delta** (two-phase, low T, intermediate Zr)
- **delta** (high Zr, low T)
- **beta** (BCC, intermediate T, low Zr)
- **beta+gamma** (two-phase)
- **gamma** (high T, fully BCC)

The phase boundaries are parameterised by six temperature lines (line1–line6)
each interpolated as a function of T and Pu content.

### 1.2 Thermal Conductivity  (Flamb.for — Aydin/Karahan 2009)

#### Fresh fuel

k = (A + B·T + C·T²)  [W/(m·K)]

where the coefficients depend on composition (pu, zr weight fractions):

| Coefficient | Formula |
|-------------|---------|
| A | 17.5 · [(1 – 2.23·zr)/(1 + 1.61·zr) – 2.62·pu] |
| B | 1.54×10⁻² · (1 + 0.061·zr)/(1 + 1.61·zr) |
| C | 9.38×10⁻⁶ · (1 – 2.7·pu) |

#### Irradiation correction (detailed model with Na infiltration)

The irradiation-corrected thermal conductivity accounts for sodium filling the fuel
pores and reduces the gas-filled porosity effect:

```
pna  = psod × (ε_tot – ε_gas)    [sodium-filled porosity]
pgas = ε_gas + (ε_tot – ε_gas) × (1 – psod)    [remaining gas porosity]
pc   = 1 – 3·pna/(1–pgas) × (1–k_Na/k_U) / (1.163 + 1.837·k_Na/k_U)
k    = k_fresh × pc × (1–pgas)^1.5
```

where:
- psod: sodium infiltration fraction (from Sinf model, see §2)
- k_Na: sodium thermal conductivity [W/(m·K)] (MFUEL correlation)
- ε_tot, ε_gas: total and gas-filled porosities (from fuel swelling model)

#### Empirical correction (simplified)

An empirical correction based on burnup (Karahan 2009):

| Burnup (FIMA) | k_correction |
|---------------|-------------|
| ≤ 2% | (1 – 0.135·B) / (1 + 1.7·0.135·B) |
| 2–5% | 0.5 + 0.0667·(B–2) |
| > 5% | 0.7 (constant) |

where B = bup_FIMA × 100 [at%].

### 1.3 Thermal Expansion  (Ftexp.for — Karahan 2009)

Linear thermal expansion strain referenced to 298.15 K:

```
if T ≤ 877 K:   ε = 1.76×10⁻⁵ × (T – 298.15)
if T ≥ 936 K:   ε = 2.01×10⁻⁵ × (T – 298.15)
877 < T < 936:  linear interpolation
```

The slope change at 877–936 K corresponds to the alpha→beta phase transition.

### 1.4 Elastic Modulus  (Felmod.for — Karahan 2009, Hofman 1985)

```
E = 1.6×10¹¹ × (1 – 1.2·por) × [(1+0.17·zr)/(1+1.34·zr) – pu]
              × [1 – 1.06·(T–588)/1405]   [Pa]
```

where por = 1 – ρ/ρ₀ is the current porosity.

### 1.5 Poisson Ratio  (Fpoir.for)

```
ν = 0.24 × (1 – 0.8·por) × [(1+3.4·zr)/(1+1.9·zr)]
           × [1 + 1.2·(T–588)/1405]
```

### 1.6 Heat Capacity  (Fcp.for — Karahan 2009, Hales 2016 correction)

```
if Tc ≤ 600°C:  cp = (26.58 + 0.027·Tc) / ma    [J/(kg·K)]
if Tc ≥ 650°C:  cp = (15.84 + 0.026·Tc) / ma
600 < Tc < 650: linear interpolation
```

where ma = average atomic mass of the U-Pu-Zr mixture [kg/mol]:
ma = 238.02891·xU + 244.06·xPu + 91.22·xZr  (weight fraction weighted)

### 1.7 Creep Rate  (Fcreep.for — Karahan 2009, Hofman 1985, Gruber 1986)

Two temperature regimes:

**Low T (< 923 K):**
```
ε̇_creep = [5000·σ + 6·σ^4.5] × exp(–26170/T)
         + 7.7×10⁻²³ × ḟ × 10⁻⁶ × σ     [1/s]
```

**High T (≥ 923 K):**
```
ε̇_creep = 0.08 × σ^3 × exp(–14350/T)    [1/s]
```

where σ = equivalent (von Mises) stress [MPa], ḟ = fission rate density [fiss/(m³·s)].

---

## 2. Sodium Infiltration Model  (Sinf.for)

Sodium from the coolant infiltrates the fuel porosity after fuel-clad gap closure.
The sodium infiltration fraction psod (0–1) is applied only in:
- alpha or alpha_delta phase nodes
- Nodes at r > 0.6 × r_max (outer 40% of fuel)

The infiltration fraction depends on gap contact state:

| Gap state | psod |
|-----------|------|
| open | 0.6 |
| soft (fuel touches cladding via relocation) | 0.6 |
| clos (hard contact) | max(0.3, 0.6 – 5·(Bu – Bu_hard)) |

where Bu and Bu_hard are the current and hard-contact burnups in FIMA.

---

## 3. Zirconium Redistribution  (Zrdist.for — SAS4A model)

Zirconium migrates radially driven by concentration gradients (diffusion) and
temperature gradients (thermotransport). The governing equation is a Fick-type
diffusion with a Soret transport term:

```
∂c_Zr/∂t = (1/r) ∂/∂r [r × (–D·∂c_Zr/∂r – D·Q/(R·T²) × c_Zr × ∂T/∂r)]
```

where:
- c_Zr: Zr atomic density [atom/m³]
- D: effective inter-diffusion coefficient [m²/s]
- Q: heat of transport [kJ/mol]
- R: universal gas constant [kJ/(mol·K)]

The diffusion coefficients D and heat of transport Q depend on phase and Pu content:

| Phase | D₀ (pu < 7 at%) | Q_h (pu < 7 at%) |
|-------|-----------------|------------------|
| alpha | 2×10⁻⁷ m²/s | 0 kJ/mol |
| delta | 2×10⁻⁷ m²/s | 0 kJ/mol |
| beta  | 5.7×10⁻⁶ m²/s | 0 kJ/mol |
| gamma | D₀(Zr) | –80 kJ/mol |

For pu ≥ 7 at%, D₀ values are 20× larger and Q_h values change sign.

Solubility limits:
- alpha phase: x_Zr < 0.05 (5 at%)
- beta/gamma: x_Zr < 0.05 
- delta: x_Zr > 0.40

At solubility limits, thermotransport currents are suppressed.

---

## 4. Fission Gas Release  (FGR.for — Karahan 2009, empirical)

### 4.1 Gas generation

25 atoms of fission gas (Kr, Xe) produced per 100 fissions:

```
fggen_rate = 0.25 × (qqv / 200 MeV) / N_A    [mol/(m³·s)]
```

### 4.2 Release fraction (empirical Karahan 2009)

Based on EBR-II irradiation data for U-Pu-Zr fuel:

```
if Bu_avg < 0.8 at%:  fgr = 0
else:  fgr = 0.8 × (1 – exp(–Bu_local[at%] / 1.8))
```

where Bu_avg and Bu_local are burnup in atom percent (at%).

The threshold of 0.8 at% corresponds to the onset of interconnected porosity in
metallic fuel, which allows gas transport to the free volume.

### 4.3 Gas pressure

```
gpres = (mu₀ + n_fgrel) × R / V_T × 10⁻⁶   [MPa]
```

where V_T = Σ_j [π·dz·(r_ci² – r_fo²)/T_gap,j] + V_plenum/T_plenum is the
thermal-volume factor summed over all axial layers.

---

## 5. Fuel Swelling  (simplified Ogata 1999)

For metallic fuel, swelling is dominated by solid fission products:

```
d(ε_vol)/dt = 0.015 × d(Bu[at%])/dt   (1.5% per at%)
```

The GRSIS model (Karahan 2009) that tracks gas bubble populations (closed and open
bubbles separately) is not implemented in this version. The 1.5%/at% solid FP model
is the simpler baseline from Ogata (1999).

The volumetric swelling is decomposed into linear strain per direction,
partitioned according to the gap-contact state (`flag`, Section 8) — this is
the Baseir.for anisotropy partition, implemented in
`apps/fred_m_na/FredMNaGapBehavior.cpp` and consumed by
`apps/fred_m_na/FredMNaStressStrain.cpp`:

```
open / clos:  Δε_z = Δε_h = Δε_r = max(0, ΔSw_tot / 3)      (isotropic)
soft:         Δε_z = 0,  Δε_h = Δε_r = 0.4995 × ΔSw_tot      (radial/hoop only)
```

Under `soft` contact the fuel is constrained axially by localized contact
before the bulk gap closes, so new swelling is redirected entirely into
hoop/radial growth rather than axial elongation (0.4995 + 0.4995 = 0.999, not
1.0 — an intentional ~0.1% volume loss/roundoff already present in the legacy
Fortran, not a transcription error). The volumetric ODE state `efs` (used for
porosity/density bookkeeping, e.g. `nd.poros_tot`) still tracks `Sw_tot/3`
unconditionally regardless of contact state; only the *directional*
accumulators (`efsz`/`efsh`/`efsr`, `FredMNaLayerState`-only, not part of the
shared platform `AxialLayerState`) follow the table above.

---

## 6. Cladding: HT-9 Ferritic-Martensitic Steel

### 6.1 Thermal Conductivity  (Clamb.for)

```
k = 29.65 – 0.06668·T + 2.184×10⁻⁴·T² – 2.527×10⁻⁷·T³ + 9.621×10⁻¹¹·T⁴   [W/(m·K)]
```

### 6.2 Thermal Expansion  (Ctexp.for — Karahan 2007)

```
ε = 10⁻² × (–0.2191 + 5.678×10⁻⁴·T + 8.111×10⁻⁷·T² – 2.576×10⁻¹⁰·T³)
```

### 6.3 Elastic Modulus  (Celmod.for)

```
E = 2.137×10¹¹ – 1.0274×10⁸·T   [Pa]   (capped at Tc ≤ 800°C)
```

### 6.4 Poisson Ratio  (Cpoir.for)

```
ν = 0.5 × (2.137×10⁵ – 102.74·Tc) / (8.964×10⁴ – 53.78·Tc) – 1.0
```

### 6.5 Creep  (Ccreep.for)

HT-9 creep has both irradiation (in-pile) and thermal components.

**Irradiation creep:**
```
ε̇_irr = (1.83×10⁻⁴ + 2.59×10¹⁴·exp(–7.3×10⁴/(R·T))) × Φ × σ^1.3   [%/h]
```

**Thermal creep (primary, secondary, tertiary contributions):**
See `Ccreep.for` for the full seven-term expression.

Neutron flux proxy: Φ ≈ qqv × 8.4×10⁵ × 10⁻²² [10²² n/cm²/s]

### 6.6 Void Swelling  (Cswel.for — Hofmann 1985)

```
s₀ = 0.01·R × [Φ_eq + (1/0.75)·ln((1+exp(0.75·(14.2–Φ_eq))) / (1+exp(0.75·14.2)))]
D  = 0.01 × 0.15 × (1 – exp(–0.1·Φ_eq))
ε_swel = s₀ + D
```

where R = 0.085·exp(–10⁻⁴·(T[°C]–400)²) and Φ_eq is the neutron fluence [10²² n/cm²].

---

## 7. Cladding Wastage  (Clanth.for — Karahan 2007)

Lanthanide fission products (La, Ce, Nd, Pr) from the fuel diffuse into the
HT-9 cladding and form low-melting eutectic phases, causing chemical wastage.

**Lanthanide diffusion coefficient in HT-9:**
```
D_l = D₀_l × exp(–Q_l / (R·T_ci))   [m²/s]
```
where T_ci is the cladding inner surface temperature.

**Wastage rate (square-root of time growth):**
```
d(x_waste)/dt = 0.5 × (c_l – c_sol,fuel) / (c_sol,clad – c_sol,fuel)
              × sqrt(D_l / t)
```

Applied only when gap is in 'soft' or 'clos' contact state.

Parameters (fitted to SAS4A data):
- Q_l = 333962 J/mol
- D₀_l = 1.6×10⁶ m²/s
- Fission yield of lanthanides: Y_l = 15.29×10⁹ J⁻¹

---

## 8. Gap Contact States for U-Pu-Zr

Unlike oxide fuels, U-Pu-Zr metallic fuel can swell via relocation and achieve
'soft contact' before hard pellet-clad mechanical interaction:

| State | Condition | Physics |
|-------|-----------|---------|
| open | r_ci – r_fo > ruff + rufc | Gas conduction + radiation |
| soft | r_ci – r_loc ≤ ruff + rufc | Sodium infiltration begins; Na conductivity enhancement; wastage begins |
| clos | mechanical constraint | Hard contact; pfc > 0; axial strain coupling |

The 'soft' state is driven by fuel relocation (Fueloc model, fanis factor). The
effective relocated fuel outer radius is:
```
r_loc = r_fo + f_anis × gap₀
```
where f_anis depends on the linear power and Pu content (Fanis.for).

The `flag` transition is a **monotonic ratchet** (`open → soft → clos`, never
reverses for U-Pu-Zr fuel — legacy's reopening check is gated on
`fmat.ne.'upuzr'`) computed once per accepted IDA step in
`FredMNaSolver::afterAcceptedStep` (no SUNDIALS root-finder is used for this
model, unlike FRED-ROD/FRED-OX's two-state `gapOpen` — see
`doc/architecture.md`'s "FRED-M-Na: three-state gap-contact model" section).

`flag` affects five distinct pieces of physics, summarized here (formulas
elsewhere in this document):

| # | Physics | open | soft | clos |
|---|---------|------|------|------|
| 1 | Swelling anisotropy (Section 5) | isotropic | radial/hoop only | isotropic |
| 2 | Fuel axial force balance | fuel-only + gas pressure | combined fuel+clad | combined fuel+clad |
| 3 | Fuel outer radial BC | `σ_fr = -p_gas` | strain continuity (`ε_fz - ε_z = offset`) | strain continuity |
| 4 | Cladding axial/interface BC | clad-only equilibrium | traction continuity | traction continuity |
| 5 | Cladding inner radial BC | `σ_r = -p_gas` | `σ_r = -p_gas` (gas still surrounds clad ID) | geometric contact (`gap = roughness`) |

Items 2–5 are implemented in `apps/fred_m_na/FredMNaStressStrain.cpp`, a
`.cpp`-only fork of the shared platform `StressStrain.cpp` (composition, not
subclassing, since the platform class's `computeResiduals` is non-virtual).
The shared platform mechanics/gap-state files
(`FuelRodState.hpp`, `StressStrain.{hpp,cpp}`, `GapStateManager.{hpp,cpp}`,
`RodResiduals.{hpp,cpp}`) are unmodified; FRED-ROD/FRED-OX are unaffected.

---

## 9. Burnup Conventions

Two burnup metrics are tracked:

| Variable | Units | Formula |
|----------|-------|---------|
| bup | MWd/kgU | d(bup)/dt = qqv / (rof₀ × 8.64×10⁴) |
| bup_FIMA | FIMA (fraction of initial metal atoms) | d(bup_FIMA)/dt = qqv / (3.204354×10⁻¹¹ × n_HM) |

where n_HM = rof₀ × N_A / M_HM is the heavy metal atom density.

The conversion is approximately: 1 at% ≈ 9.4 MWd/kgU for U-19%Pu-10%Zr.

---

## 10. References

1. **Karahan, A.** (2009). *Modelling of thermo-mechanical and irradiation behavior of
   metallic and oxide fuels for sodium fast reactors*. PhD thesis, MIT.

2. **Hofman, G.L.** (1985). *Metallic Fuels Handbook*. ANL-AF-1 report, Argonne National
   Laboratory.

3. **Ogata, T.** (1999). Development and validation of ALFUS: an irradiation behavior
   analysis code for metallic fast reactor fuels. *Nucl. Technol.*, 128, 113.

4. **SAS4A Users' Manual** (2008). Argonne National Laboratory, ANL/NE-08/1.

5. **Timpano, D.** (2018-2024). FRED_M_OCT24.SRC — FRED metallic fuel branch.
   EPFL, Laboratory for Reactor Physics.

6. **Aydin, K.** (PhD Thesis). U-Pu-Zr fuel thermal conductivity correlation.
   Referenced in Karahan (2009).
