# Advanced Prototyping: Python → C++ workflow

This example demonstrates a two-step development pattern for adding new physics
to FRED-ROD. The worked example is a **fission gas release model for gap conductance**:
as fuel burns up, fission gas (Xe, Kr) is released into the helium fill gas,
lowering gap thermal conductivity and raising the peak fuel temperature.

The two-step workflow is:

1. **Prototype in Python** — subclass `fred.GapMaterial` in pure Python.
   Iterate on the correlations with no build step.  Check physical trends,
   run parameter studies, compare against reference data.

2. **Compile in C++** — translate the validated model to a C++ header file,
   add a one-block pybind11 binding, rebuild.  The compiled class appears in
   Python as `fred.HeKrXe` with identical behaviour and no Python overhead
   in the inner loop.

---

## Physics: fission gas release and gap conductance

### The problem

Fresh fuel pins are filled with helium to a modest pressure (typically 0.1 MPa).
Helium has excellent thermal conductivity (~0.25 W/mK at 700 K), making the
gap an efficient conductor despite its small width (~100 µm as-fabricated).

During irradiation, a fraction of the fission products are noble gases — primarily
xenon (88.5 %) and krypton (7.7 %), with a small helium component (3.8 %).
These gases are produced in the fuel matrix and, at sufficient temperature and
burnup, diffuse to grain boundaries and are released into the free gas volume.

Xe and Kr have conductivities roughly 10–20× lower than He at typical gap
temperatures.  Even a moderate fission gas release fraction significantly
degrades gap conductance.

### Fission gas release fraction

The Waltar-Reynolds model gives the fraction of generated gas released to the
free volume for the restructured (high-temperature) zone:

```
bup_MWd = bup_atpct * 9.4           (conversion: 1 at% ≈ 9.4 MWd/kgHM for oxide)
a        = 4.7 / bup_MWd * (1 - exp(-bup_MWd / 5.9))
FGR      = max(0,  1 - a)
```

This is an upper-bound estimate (restructured zone only).  At low burnup
(< 1 MWd/kgHM) the formula saturates at FGR ≈ 0; at high burnup (> 50 MWd/kgHM)
FGR approaches 1.

### Molar balance

The total moles of gap gas are the sum of the initial helium fill and the
released fission gas.  The fission gas production rate is approximately 0.25
gas atoms per fission (stoichiometric fixed yield for oxide fuel).  For
typical pin parameters (10 g of fuel, M_HM = 238 g/mol):

```
n_He0   = P_fill * V_free / (R * T_ref)     # ideal gas law at fabrication
          Typical: P=0.1 MPa, V=1.5 cm³, T=293 K  →  n_He0 ≈ 6.15e-5 mol

fggen   = 0.25 * (bup_atpct / 100) * m_fuel / M_HM   # total gas generated
fgrel   = FGR * fggen                                  # released to gap

n_Xe    = 0.8846 * fgrel          # FRED.f90 gaphtc composition
n_Kr    = 0.0769 * fgrel
n_He_FG = 0.0385 * fgrel

n_total = n_He0 + n_He_FG + n_Kr + n_Xe

y_He = (n_He0 + n_He_FG) / n_total
y_Kr = n_Kr / n_total
y_Xe = n_Xe / n_total
```

### Mixture conductivity

The geometric mean rule (used in FRED.f90 `gaphtc`) gives the mixture
conductivity from individual species conductivities and mole fractions:

```
k_He(T) = 2.639e-3 * T^0.7085   [W/(m·K)]
k_Kr(T) = 8.247e-5 * T^0.8363
k_Xe(T) = 4.351e-5 * T^0.8618

k_mix(T) = k_He^y_He * k_Kr^y_Kr * k_Xe^y_Xe
```

This is used directly by the FRED-ROD platform's Lanning-Hann gap conductance
model, which computes the total gap conductance from the medium conductivity
and the actual gap width (gas-jump correction included).

---

## Step 1: Python prototyping

### The GapMaterial interface

`fred.GapMaterial` is an abstract base class with one mandatory method:

```python
class MyGapMaterial(fred.GapMaterial):
    def gapConductivity(self, T: float) -> float:
        """Gap-medium thermal conductivity [W/(m·K)] at temperature T [K]."""
        ...
```

The platform calls `gapConductivity(T)` at each solver step to evaluate the
gap medium conductivity.  Gap conductance (accounting for gas-jump resistances
and the actual gap width) is computed by the platform from this value — the
material class does not need to know the gap width.

The `isGasBond()` method (default `True`) controls how the solver converts
conductivity to a gap conductance:

- **Gas bond** (`isGasBond() = True`): the solver applies the Lanning-Hann
  gas-jump model to `k` and adds linearised radiation — `h = f(k, gap_width, P, T)`.
- **Liquid bond** (`isGasBond() = False`): the solver uses `h = k / gap_eff`
  directly, with no gas-jump correction.  Override `clampConductance()` if the
  material needs to enforce min/max bounds on contact (e.g. sodium forced to
  1e6 W/m²K when the gap closes).

For both bond types, `gapConductivity` is the only method that must be
overridden.

### Running the Python example

```bash
cd examples/fred_rod/advanced_prototyping
python run_python.py
```

The script:
- Prints a table of gap-gas composition and conductance vs burnup.
- Runs three FRED-ROD simulations (bup = 0, 1, 5 at%) and prints the
  steady-state peak fuel temperature for each case.

No build step is required.  Any changes to the Python class take effect
immediately on the next run.

### Iterating on the Python model

During prototyping you can freely:
- Swap the FGR correlation (e.g. a temperature-dependent model, or an
  empirical tabulation from experiment).
- Change the gas composition fractions.
- Add new species (e.g. Ar fill gas — replace y_He by y_He and y_Ar).
- Drive `bup_atpct` from a time-varying burnup array and update the gap
  material between solver steps.

None of these require touching any C++ file.

---

## Step 2: Compiled C++ implementation

Once the Python model is validated, translating it to C++ has three steps:
(1) write the header, (2) add a pybind11 binding, (3) rebuild.

`HeKrXe.hpp` (already in the repository) is the compiled result of the model
prototyped in `run_python.py`.  The instructions below describe the general
process so you can follow them for your own model.

### 2a. Write the C++ header

Create a header file in the gap-material subfolder:

```
src/apps/fred_rod/gapmaterial/MyGapMaterial.hpp
```

Minimal template:

```cpp
#pragma once
#include "platform/GapMaterial.hpp"
#include <cmath>
#include <algorithm>

namespace fred {

class MyGapMaterial : public GapMaterial {
public:
    explicit MyGapMaterial(double bup_atpct = 0.0) {
        // ... compute m_y_He, m_y_Kr, m_y_Xe from bup_atpct ...
    }

    double gapConductivity(double T) const override {
        // ... return k_mix(T) ...
    }

private:
    double m_y_He, m_y_Kr, m_y_Xe;
};

} // namespace fred
```

Rules:
- Inherit from `fred::GapMaterial`.
- Override `gapConductivity(double T) const` — this is the only pure-virtual
  method you must implement.
- Keep the class header-only if there is no non-trivial state (no `.cpp` needed).
- If you need a `.cpp` file, add it to `ROD_LIB_SRCS` in the `Makefile` (see below).

`HeKrXe.hpp` is a header-only class, so no Makefile source-list change is needed
for it.

### 2b. Add a Python binding

Open `src/apps/fred_rod/bindings.cpp` and make two edits:

**Add the include** (near the top, with the other gapmaterial includes):

```cpp
#include "apps/fred_rod/gapmaterial/MyGapMaterial.hpp"
```

**Add the binding block** (inside `bind_fred_rod`, after the He binding):

```cpp
py::class_<MyGapMaterial, GapMaterial>(m, "MyGapMaterial")
    .def(py::init<double>(),
         py::arg("bup_atpct") = 0.0,
         "One-line description of what the material models.\n"
         "bup_atpct: local fuel burnup in atomic percent.");
```

The `py::class_<MyGapMaterial, GapMaterial>` template argument tells pybind11
that `MyGapMaterial` inherits from the already-registered `GapMaterial` base.
This is required so Python code can pass a `MyGapMaterial` instance wherever
a `fred.GapMaterial` is expected (e.g. the `gap` argument to `FredRodSolver`).

### 2c. Makefile — header-only class (no change needed)

For a header-only class the Makefile does not need to change: the header is
compiled as part of `bindings.cpp` when it includes it.

If you need a `.cpp` file (e.g. to separate a long implementation from the
header), add it to the `ROD_LIB_SRCS` list in the Makefile:

```make
ROD_LIB_SRCS = \
    ...
    $(SRC)/apps/fred_rod/gapmaterial/MyGapMaterial.cpp \
    ...
```

The compile rule `$(BUILD)/%.o: $(SRC)/%.cpp` and the `@mkdir -p $(dir $@)`
guard handle any new subdirectory automatically — no further Makefile changes
are needed.

### 2d. Rebuild

From the repository root:

```bash
make fred-rod
```

This recompiles only the changed translation units and relinks `fred_rod.*.so`.
A typical incremental rebuild (bindings.cpp only) takes a few seconds.

After a successful build, the class is immediately available:

```python
import fred_rod as fred
gap = fred.MyGapMaterial(bup_atpct=2.0)
print(gap.gap_conductivity(700.0))   # W/(m·K)
```

### Running the C++ example

```bash
cd examples/fred_rod/advanced_prototyping
python run_cpp.py
```

The script runs the same three simulations as `run_python.py` using the
compiled `fred.HeKrXe` and verifies that conductance values agree with the
Python reference to machine precision.

---

## Files in this example

```
advanced_prototyping/
  README.md          — this file
  run_python.py      — Step 1: HeKrXePy class, prototype simulation
  run_cpp.py         — Step 2: fred.HeKrXe, compiled equivalent

Related source files:
  src/apps/fred_rod/gapmaterial/HeKrXe.hpp   — C++ class
  src/apps/fred_rod/bindings.cpp             — pybind11 binding for HeKrXe
```

---

## Extending this pattern

The same workflow applies to fuel pellet and cladding materials:

| Interface class          | Abstract method(s) to override             |
|--------------------------|--------------------------------------------|
| `fred.GapMaterial`       | `gapConductivity(T)`                       |
| `fred.FuelPelletMaterial`| `thermalConductivity`, `heatCapacity`, `thermalExpansionStrain`, `youngsModulus`, `poissonRatio`, `referenceDensity`, `theoreticalDensity`, `meltingTemperature` |
| `fred.CladdingMaterial`  | `thermalConductivity`, `heatCapacity`, `thermalExpansionStrain`, `youngsModulus`, `poissonRatio`, `meyerHardness`, `referenceDensity` |

See `examples/fred_rod/python_materials/run.py` for a full fuel + cladding
prototyping example.

A similar workflow to develop specific correlations can also be applied to the other applications
such as FRED-OX and is encouraged when modifications involve only changes to variable values without
introducing new physics. 
