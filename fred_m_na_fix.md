# FRED-M-Na fixes: Zr-redistribution phase/conductivity coupling

## Background

After Zr redistribution runs each accepted step, each fuel node holds an updated
`zr_wf` (Zr weight fraction) that differs from the as-fabricated nominal value
stored in the `UPuZr` object (`m_zr`).  Three downstream quantities depended on
composition but were not updated to use the redistributed profile:

1. **Phase** (`nd.phase`) — determines which phase-field diffusion coefficient is
   used in the next redistribution sub-step, and feeds `sodiumInfiltration` which
   computes the sodium infiltration fraction `nd.psod`.
2. **Thermal conductivity** — the radial heat conduction equation needs `k(r)` that
   varies with `r` as Zr redistributes.  The old code used a single scalar
   `k_irr_factor` evaluated at one midpoint node with the nominal `m_zr`, which
   misses both the per-node composition variation and the post-redistribution profile.

The radial conduction equation is:

```
q''' = -1/r · d/dr( r · k(r) · dT/dr )
```

A scalar k is physically incorrect when `k` varies radially.  `upuzrFreshConductivity`
is nonlinear in `zr`, so volume-averaging `zr` and evaluating once is also wrong.
The correct fix is per-interval conductivity.

---

## Changes

### 1. `UPuZr.hpp` / `UPuZr.cpp` — per-node conductivity overloads

Two new methods accept explicit local composition `(pu_wf, zr_wf)` instead of
using the stored nominal `m_pu`/`m_zr`:

```cpp
double thermalConductivityLocal(double T, double pu_wf, double zr_wf) const;
double thermalConductivityIrradiatedLocal(double T_K, double pu_wf, double zr_wf,
                                           double bup_FIMA,
                                           double poros_tot, double poros_gas,
                                           double psod) const;
```

All three conductivity model branches (`DetailedNaSodium`, `EmpiricalBurnup`,
`EsfrSimple`) are forwarded with the local composition.

---

### 2. `platform/FuelRodState.hpp` — `k_fuel_per_node` field

Added to `AxialLayerState`:

```cpp
std::vector<double> k_fuel_per_node;  // size nf when populated; empty → scalar path
```

`HeatConduction` uses this when non-empty; FRED-ROD and FRED-OX leave it empty and
continue using the scalar `k_irr_factor` path unchanged.

---

### 3. `platform/HeatConduction.cpp` — per-interval conductivity dispatch

The fuel-fuel boundary loop now branches on whether `k_fuel_per_node` is populated:

```cpp
for (int i = 0; i < nf-1; ++i) {
    double kf;
    if (!s.k_fuel_per_node.empty()) {
        // Average adjacent node values for the interface.
        kf = 0.5 * (s.k_fuel_per_node[i] + s.k_fuel_per_node[i+1]);
    } else {
        double T_half = 0.5 * (s.T[i] + s.T[i+1]);
        kf = m_fuel.thermalConductivity(T_half) * s.k_irr_factor;
    }
    qf[i] = kf * (s.T[i] - s.T[i+1]) / dr * 2.0 * PI * m_geom.rad_half[i];
}
```

FRED-ROD and FRED-OX are unaffected (`k_fuel_per_node` remains empty for them).

---

### 4. `FredMNaSolver.cpp` — `afterAcceptedStep`

#### Bug fix A: initial phase loop uses local `zr_wf`

```cpp
// before:
auto ph = upuzrPhase(s.T[i], pu, zr);       // zr = nominal m_zr for all nodes
// after:
auto ph = upuzrPhase(s.T[i], pu, nd.zr_wf); // local redistributed value
```

#### Change B: `fillKFuelPerNode` replaces scalar `computeKirrFactor`

```cpp
auto fillKFuelPerNode = [&]() {
    s.k_fuel_per_node.resize(nf);
    for (int i = 0; i < nf; ++i) {
        const auto& nd = s.nodes[i];
        s.k_fuel_per_node[i] = m_fuel.thermalConductivityIrradiatedLocal(
            s.T[i], nd.pu_wf, nd.zr_wf,
            s.bup_FIMA, nd.poros_tot, nd.poros_gas, nd.psod);
    }
};
```

Each node stores `k` in W/(m·K) directly.  `HeatConduction` then uses
`0.5*(k[i]+k[i+1])` for interval `[i, i+1]`, giving a spatially-resolved
conductivity profile that tracks the radial Zr gradient.

#### Bug fix C: phase/psod/conductivity re-evaluation after redistribution

After `upuzrZirconiumRedistribution` writes back and the outermost node is
mirrored, a second pass recomputes phase, psod, and `k_fuel_per_node`:

```cpp
for (int i = 0; i < nf; ++i) {
    auto& nd = s.nodes[i];
    auto ph = upuzrPhase(s.T[i], pu, nd.zr_wf);
    nd.phase = ph.phase; nd.pfrac = ph.pfrac;
    nd.psod  = sodiumInfiltration(..., nd.phase, s.flag);
}
fillKFuelPerNode();
```

GRSIS (which runs after) therefore sees the phase corresponding to the current
step's Zr profile, and `HeatConduction` in the next Newton solve uses the
just-redistributed conductivity profile.

---

### 5. `FredMNaResiduals.cpp` — `syncAuxLayerState`

`k_fuel_per_node` is copied alongside `k_irr_factor` so the residual object's
layer state always holds the latest conductivity profile:

```cpp
dst.k_fuel_per_node = src.k_fuel_per_node;
```

---

## Effect summary

| Quantity | Before | After |
|---|---|---|
| `nd.phase` (initial loop) | nominal `zr` for all nodes | local `nd.zr_wf` |
| `nd.phase` after redistribution | never updated | re-evaluated with new `zr_wf` |
| `nd.psod` after redistribution | stale (previous step phase) | updated with new phase |
| Conductivity representation | scalar `k_irr_factor` at one midpoint | per-node `k_fuel_per_node[i]` for all nf nodes |
| Conductivity composition | global nominal `m_zr` | local `nd.zr_wf`, `nd.pu_wf` per node |
| Conductivity spatial resolution | single value for whole layer | per-interval: `0.5*(k[i]+k[i+1])` |
| Conductivity timing | before redistribution only | recomputed after redistribution |
| FRED-ROD / FRED-OX | unchanged | unchanged (scalar path, `k_fuel_per_node` empty) |

The `aak`, `bbk`, `cck` coefficients in `upuzrFreshConductivity` are now evaluated
at the actual local `zr_wf` per node, so the inner (Zr-depleted) and outer
(Zr-enriched) zones have self-consistent conductivities, and the FD flux between
any two nodes uses the arithmetic mean of those values.
