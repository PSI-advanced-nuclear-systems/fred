#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>

#include "FredRodSolver.hpp"
#include "platform/FuelPelletMaterial.hpp"
#include "platform/CladdingMaterial.hpp"
#include "platform/GapMaterial.hpp"
#include "apps/fred_rod/fuelpelletmaterial/UO2.hpp"
#include "apps/fred_rod/fuelpelletmaterial/MOX.hpp"
#include "apps/fred_rod/claddingmaterial/AIM1.hpp"
#include "apps/fred_rod/claddingmaterial/T91.hpp"
#include "apps/fred_rod/fuelpelletmaterial/DummyFuelPellet.hpp"
#include "apps/fred_rod/claddingmaterial/DummyCladding.hpp"
#include "apps/fred_rod/gapmaterial/DummyGapMaterial.hpp"
#include "apps/fred_rod/gapmaterial/He.hpp"
#include "apps/fred_rod/gapmaterial/HeKrXe.hpp"
#include "platform/FuelRodGeometry.hpp"

namespace py = pybind11;
using namespace fred;

// -----------------------------------------------------------------------
// Trampoline classes — allow Python subclasses to override virtual methods.
// These are only needed for FRED-ROD where users define custom material
// correlations on the Python side.
// -----------------------------------------------------------------------
namespace {

class PyGapMaterial : public GapMaterial {
public:
    double gapConductivity(double T) const override {
        PYBIND11_OVERRIDE_PURE(double, GapMaterial, gapConductivity, T);
    }
    bool isGasBond() const override {
        PYBIND11_OVERRIDE(bool, GapMaterial, isGasBond);
    }
    double clampConductance(double h_gap, bool gap_closed) const override {
        PYBIND11_OVERRIDE(double, GapMaterial, clampConductance, h_gap, gap_closed);
    }
    // Unconditionally false: a Python subclass may exist (this trampoline
    // was instantiated), so calling into it from a non-GIL-holding OpenMP
    // worker thread would be unsafe regardless of what it overrides.
    bool isThreadSafe() const override { return false; }
};

class PyFuelPelletMaterial : public FuelPelletMaterial {
public:
    double thermalConductivity(double T) const override {
        PYBIND11_OVERRIDE_PURE(double, FuelPelletMaterial, thermalConductivity, T);
    }
    double heatCapacity(double T) const override {
        PYBIND11_OVERRIDE_PURE(double, FuelPelletMaterial, heatCapacity, T);
    }
    double thermalExpansionStrain(double T) const override {
        PYBIND11_OVERRIDE_PURE(double, FuelPelletMaterial, thermalExpansionStrain, T);
    }
    double youngsModulus(double T, double density) const override {
        PYBIND11_OVERRIDE_PURE(double, FuelPelletMaterial, youngsModulus, T, density);
    }
    double poissonRatio() const override {
        PYBIND11_OVERRIDE_PURE(double, FuelPelletMaterial, poissonRatio);
    }
    double referenceDensity() const override {
        PYBIND11_OVERRIDE_PURE(double, FuelPelletMaterial, referenceDensity);
    }
    double theoreticalDensity() const override {
        PYBIND11_OVERRIDE_PURE(double, FuelPelletMaterial, theoreticalDensity);
    }
    double meltingTemperature() const override {
        PYBIND11_OVERRIDE_PURE(double, FuelPelletMaterial, meltingTemperature);
    }
    bool isThreadSafe() const override { return false; }
};

class PyCladdingMaterial : public CladdingMaterial {
public:
    double thermalConductivity(double T) const override {
        PYBIND11_OVERRIDE_PURE(double, CladdingMaterial, thermalConductivity, T);
    }
    double heatCapacity(double T) const override {
        PYBIND11_OVERRIDE_PURE(double, CladdingMaterial, heatCapacity, T);
    }
    double thermalExpansionStrain(double T) const override {
        PYBIND11_OVERRIDE_PURE(double, CladdingMaterial, thermalExpansionStrain, T);
    }
    double youngsModulus(double T) const override {
        PYBIND11_OVERRIDE_PURE(double, CladdingMaterial, youngsModulus, T);
    }
    double poissonRatio() const override {
        PYBIND11_OVERRIDE_PURE(double, CladdingMaterial, poissonRatio);
    }
    double meyerHardness(double T) const override {
        PYBIND11_OVERRIDE_PURE(double, CladdingMaterial, meyerHardness, T);
    }
    double referenceDensity() const override {
        PYBIND11_OVERRIDE_PURE(double, CladdingMaterial, referenceDensity);
    }
    bool isThreadSafe() const override { return false; }
};

} // namespace

// -----------------------------------------------------------------------
// bind_fred_rod — registers platform abstract types (with trampolines),
// shared geometry, FRED-ROD concrete materials, and FredRodSolver.
// Must be called BEFORE bind_fred_ox / bind_fred_m_na because OX and M-Na
// concrete material classes inherit from the abstract types registered here.
// -----------------------------------------------------------------------
void bind_fred_rod(py::module_& m) {

    // ---- Platform abstract material classes (trampolines allow Python subclasses) ----
    py::class_<FuelPelletMaterial, PyFuelPelletMaterial>(m, "FuelPelletMaterial")
        .def(py::init<>())
        .def("thermal_conductivity",    &FuelPelletMaterial::thermalConductivity,
             py::arg("T"), "Thermal conductivity [W/(m·K)]")
        .def("heat_capacity",           &FuelPelletMaterial::heatCapacity,
             py::arg("T"), "Specific heat [J/(kg·K)]")
        .def("thermal_expansion_strain",&FuelPelletMaterial::thermalExpansionStrain,
             py::arg("T"), "Linear thermal expansion strain relative to 293.15 K [-]")
        .def("youngs_modulus",          &FuelPelletMaterial::youngsModulus,
             py::arg("T"), py::arg("density"), "Young's modulus [MPa]")
        .def("poisson_ratio",           &FuelPelletMaterial::poissonRatio,
             "Poisson's ratio [-]")
        .def("reference_density",       &FuelPelletMaterial::referenceDensity,
             "As-fabricated density [kg/m3]")
        .def("theoretical_density",     &FuelPelletMaterial::theoreticalDensity,
             "Fully dense (theoretical) density [kg/m3]")
        .def("melting_temperature",     &FuelPelletMaterial::meltingTemperature,
             "Melting temperature [K]");

    py::class_<CladdingMaterial, PyCladdingMaterial>(m, "CladdingMaterial")
        .def(py::init<>())
        .def("thermal_conductivity",    &CladdingMaterial::thermalConductivity, py::arg("T"))
        .def("heat_capacity",           &CladdingMaterial::heatCapacity,        py::arg("T"))
        .def("thermal_expansion_strain",&CladdingMaterial::thermalExpansionStrain, py::arg("T"))
        .def("youngs_modulus",          &CladdingMaterial::youngsModulus,       py::arg("T"))
        .def("poisson_ratio",           &CladdingMaterial::poissonRatio)
        .def("meyer_hardness",          &CladdingMaterial::meyerHardness,       py::arg("T"))
        .def("reference_density",       &CladdingMaterial::referenceDensity);

    py::class_<GapMaterial, PyGapMaterial>(m, "GapMaterial")
        .def(py::init<>())
        .def("gap_conductivity", &GapMaterial::gapConductivity, py::arg("T"),
             "Gap-medium thermal conductivity k [W/(m·K)] at temperature T [K].\n"
             "For gas bond: the platform converts k to h via the Lanning-Hann model.\n"
             "For liquid bond (is_gas_bond=False): the platform computes h = k / max(gap, roughness).")
        .def("is_gas_bond", &GapMaterial::isGasBond,
             "Return True for gas-filled gaps (default); False for liquid-bond gaps (Na, Pb, ...).\n"
             "Override in Python to select the liquid-bond conductance path.")
        .def("clamp_conductance", &GapMaterial::clampConductance,
             py::arg("h_gap"), py::arg("gap_closed"),
             "Optional post-calculation clamp on h_gap [W/(m²·K)].\n"
             "Called after h = k/gap_eff for liquid-bond materials.\n"
             "Override to enforce application-specific min/max bounds (e.g. sodium floor on closure).\n"
             "Default implementation is identity (no clamping).");

    // ---- Shared geometry ----
    py::class_<FuelRodGeometry>(m, "FuelRodGeometry")
        .def(py::init<>())
        .def_readwrite("nf",   &FuelRodGeometry::nf,   "Number of fuel radial nodes")
        .def_readwrite("nc",   &FuelRodGeometry::nc,   "Number of cladding radial nodes")
        .def_readwrite("nz",   &FuelRodGeometry::nz,   "Number of axial layers")
        .def_readwrite("rfi0", &FuelRodGeometry::rfi0, "Initial inner fuel radius [m] (0 for solid)")
        .def_readwrite("rfo0", &FuelRodGeometry::rfo0, "Initial outer fuel radius [m]")
        .def_readwrite("rci0", &FuelRodGeometry::rci0, "Initial inner cladding radius [m]")
        .def_readwrite("rco0", &FuelRodGeometry::rco0, "Initial outer cladding radius [m]")
        .def_readwrite("dz0",  &FuelRodGeometry::dz0,  "Axial layer heights [m], list of length nz")
        .def_readwrite("vgp",  &FuelRodGeometry::vgp,  "Gas plenum volume [m3]")
        .def_readwrite("ruff", &FuelRodGeometry::ruff, "Fuel surface roughness [m]")
        .def_readwrite("rufc", &FuelRodGeometry::rufc, "Cladding inner surface roughness [m]")
        .def("build", &FuelRodGeometry::build,
             "Compute all derived geometry arrays (call after setting all parameters)");

    // ---- FRED-ROD concrete materials ----
    py::class_<UO2, FuelPelletMaterial>(m, "UO2")
        .def(py::init<double>(), py::arg("reference_density") = 10400.0,
             "UO2 fuel pellet.  reference_density: as-fabricated density [kg/m3]")
        .def("thermal_conductivity",    &UO2::thermalConductivity,    py::arg("T"))
        .def("heat_capacity",           &UO2::heatCapacity,           py::arg("T"))
        .def("thermal_expansion_strain",&UO2::thermalExpansionStrain, py::arg("T"))
        .def("youngs_modulus",          &UO2::youngsModulus,          py::arg("T"), py::arg("density"))
        .def("poisson_ratio",           &UO2::poissonRatio)
        .def("reference_density",       &UO2::referenceDensity)
        .def("theoretical_density",     &UO2::theoreticalDensity)
        .def("melting_temperature",     &UO2::meltingTemperature);

    py::class_<MOX, FuelPelletMaterial>(m, "MOX")
        .def(py::init<double, double>(),
             py::arg("pu_content"), py::arg("reference_density") = -1.0,
             "MOX fuel pellet (Philipponneau 1992 conductivity, Popov ORNL/TM-2000/351 Cp,\n"
             "MATPRO thermal expansion and Young's modulus).\n"
             "pu_content      : plutonium mole fraction [-], e.g. 0.15 (low) or 0.30 (high)\n"
             "reference_density: as-fabricated density [kg/m3]; <=0 → 95% of theoretical density")
        .def("thermal_conductivity",    &MOX::thermalConductivity,    py::arg("T"))
        .def("heat_capacity",           &MOX::heatCapacity,           py::arg("T"))
        .def("thermal_expansion_strain",&MOX::thermalExpansionStrain, py::arg("T"))
        .def("youngs_modulus",          &MOX::youngsModulus,          py::arg("T"), py::arg("density"))
        .def("poisson_ratio",           &MOX::poissonRatio)
        .def("reference_density",       &MOX::referenceDensity)
        .def("theoretical_density",     &MOX::theoreticalDensity)
        .def("melting_temperature",     &MOX::meltingTemperature);

    py::class_<AIM1, CladdingMaterial>(m, "AIM1")
        .def(py::init<double>(), py::arg("reference_density") = 7900.0,
             "AIM1 austenitic stainless steel cladding.")
        .def("thermal_conductivity",    &AIM1::thermalConductivity,    py::arg("T"))
        .def("heat_capacity",           &AIM1::heatCapacity,           py::arg("T"))
        .def("thermal_expansion_strain",&AIM1::thermalExpansionStrain, py::arg("T"))
        .def("youngs_modulus",          &AIM1::youngsModulus,          py::arg("T"))
        .def("poisson_ratio",           &AIM1::poissonRatio)
        .def("meyer_hardness",          &AIM1::meyerHardness,          py::arg("T"))
        .def("reference_density",       &AIM1::referenceDensity);

    py::class_<T91, CladdingMaterial>(m, "T91")
        .def(py::init<double>(), py::arg("reference_density") = 7750.0,
             "T91 (9Cr-1Mo-V) ferritic-martensitic steel cladding (no irradiation).")
        .def("thermal_conductivity",    &T91::thermalConductivity,    py::arg("T"))
        .def("heat_capacity",           &T91::heatCapacity,           py::arg("T"))
        .def("thermal_expansion_strain",&T91::thermalExpansionStrain, py::arg("T"))
        .def("youngs_modulus",          &T91::youngsModulus,          py::arg("T"))
        .def("poisson_ratio",           &T91::poissonRatio)
        .def("meyer_hardness",          &T91::meyerHardness,          py::arg("T"))
        .def("reference_density",       &T91::referenceDensity);

    py::class_<DummyFuelPellet, FuelPelletMaterial>(m, "DummyFuelPellet")
        .def(py::init<>(),
             "Dummy fuel pellet: k=10 W/mK, rho=10000 kg/m3, cp=100 J/kgK, "
             "CTE=10e-6/K, E=100 GPa, nu=0.3")
        .def("thermal_conductivity",    &DummyFuelPellet::thermalConductivity,    py::arg("T"))
        .def("heat_capacity",           &DummyFuelPellet::heatCapacity,           py::arg("T"))
        .def("thermal_expansion_strain",&DummyFuelPellet::thermalExpansionStrain, py::arg("T"))
        .def("youngs_modulus",          &DummyFuelPellet::youngsModulus,          py::arg("T"), py::arg("density"))
        .def("poisson_ratio",           &DummyFuelPellet::poissonRatio)
        .def("reference_density",       &DummyFuelPellet::referenceDensity)
        .def("theoretical_density",     &DummyFuelPellet::theoreticalDensity)
        .def("melting_temperature",     &DummyFuelPellet::meltingTemperature);

    py::class_<DummyCladding, CladdingMaterial>(m, "DummyCladding")
        .def(py::init<>(),
             "Dummy cladding: k=10 W/mK, rho=10000 kg/m3, cp=100 J/kgK, "
             "CTE=10e-6/K, E=100 GPa, nu=0.3, H_M=1 GPa")
        .def("thermal_conductivity",    &DummyCladding::thermalConductivity,    py::arg("T"))
        .def("heat_capacity",           &DummyCladding::heatCapacity,           py::arg("T"))
        .def("thermal_expansion_strain",&DummyCladding::thermalExpansionStrain, py::arg("T"))
        .def("youngs_modulus",          &DummyCladding::youngsModulus,          py::arg("T"))
        .def("poisson_ratio",           &DummyCladding::poissonRatio)
        .def("meyer_hardness",          &DummyCladding::meyerHardness,          py::arg("T"))
        .def("reference_density",       &DummyCladding::referenceDensity);

    py::class_<DummyGapMaterial, GapMaterial>(m, "DummyGapMaterial")
        .def(py::init<>(),
             "Gas bond: He medium conductivity.  The platform computes gap conductance\n"
             "from the actual gap width using the Lanning-Hann gas-jump model plus\n"
             "linearised radiation (matches legacy FRED gaphtc for dummy material).");

    py::class_<He, GapMaterial>(m, "He")
        .def(py::init<>(),
             "Pure helium medium conductivity for as-fabricated gap fills.\n"
             "Gas bond: the platform calls gapConductivity(T) and applies the\n"
             "Lanning-Hann gas-jump model internally.");

    py::class_<HeKrXe, GapMaterial>(m, "HeKrXe")
        .def(py::init<double>(),
             py::arg("bup_atpct") = 1.0,
             "He-Kr-Xe mixed-gas gap material for base-irradiated gap conditions.\n"
             "bup_atpct: local fuel burnup in atomic percent (at%) of HM atoms fissioned.\n"
             "Mole fractions are derived from a molar balance: ideal-gas He fill inventory\n"
             "plus released fission gas (Waltar-Reynolds FGR, 88.46% Xe/7.69% Kr/3.85% He).\n"
             "  0 at%  → pure He fill gas (as-fabricated)\n"
             "  1 at%  (≈  9.4 MWd/kgHM) → moderate base irradiation, y_Xe ≈ 0.45\n"
             "  5 at%  (≈ 47  MWd/kgHM)  → high burnup,               y_Xe ≈ 0.78\n"
             "Mixture conductivity: geometric mean rule (FRED.f90 gaphtc).")
        ;

    // ---- FredRodSolver ----
    py::class_<FredRodSolver>(m, "FredRodSolver")
        .def(py::init<const FuelRodGeometry&, const FuelPelletMaterial&,
                      const CladdingMaterial&, const GapMaterial&>(),
             py::keep_alive<1, 2>(),  // solver keeps geometry alive
             py::keep_alive<1, 3>(),  // solver keeps fuel alive
             py::keep_alive<1, 4>(),  // solver keeps clad alive
             py::keep_alive<1, 5>(),  // solver keeps gap alive
             py::arg("geometry"), py::arg("fuel"), py::arg("clad"), py::arg("gap"),
             "Pure thermo-mechanical fuel rod solver (FRED-ROD).\n"
             "No irradiation model — suitable for out-of-pile thermal tests.\n"
             "gap: GapMaterial instance (e.g. DummyGapMaterial())")
        .def("set_power_density_history", &FredRodSolver::setPowerDensityHistory,
             py::arg("times"), py::arg("qv_W_m3"),
             "Uniform power density history [W/m3] broadcast to all axial layers.")
        .def("set_power_density_history_per_layer",
             &FredRodSolver::setPowerDensityHistoryPerLayer,
             py::arg("times"), py::arg("qqv_per_layer"),
             "Per-layer power density histories [W/m3]. "
             "qqv_per_layer: list of nz value arrays, one per axial layer.")
        .def("set_coolant_temperature", &FredRodSolver::setCoolantTemperature,
             py::arg("times"), py::arg("T_K"),
             "Uniform coolant temperature history [K] broadcast to all axial layers.")
        .def("set_coolant_temperature_per_layer",
             &FredRodSolver::setCoolantTemperaturePerLayer,
             py::arg("times"), py::arg("T_per_layer"),
             "Per-layer coolant temperature histories [K]. "
             "T_per_layer: list of nz value arrays, one per axial layer.")
        .def("set_coolant_pressure",    &FredRodSolver::setCoolantPressure,
             py::arg("pcool_MPa"), "Constant coolant pressure [MPa].")
        .def("set_coolant_pressure_history", &FredRodSolver::setCoolantPressureHistory,
             py::arg("times"), py::arg("pcool_MPa"),
             "Time-varying coolant pressure history [MPa].")
        .def("set_initial_temperature",  &FredRodSolver::setInitialTemperature,
             py::arg("T0_K"), "Uniform initial temperature [K].")
        .def("set_initial_gas_pressure", &FredRodSolver::setInitialGasPressure,
             py::arg("gpres0_MPa"), "Initial fill-gas pressure [MPa] (must be called before run()).")
        .def("set_output_file", &FredRodSolver::setOutputFile,
             py::arg("filename"),
             "Path to HDF5 output file. When set, each output step is written to the file\n"
             "immediately (crash-safe) and in-memory vectors are trimmed to 1 entry.")
        .def("set_plenum_temperature_history", &FredRodSolver::setPlenumTemperatureHistory,
             py::arg("times"), py::arg("T_K"),
             "Gas plenum temperature vs time [K].")
        .def("set_tolerances",          &FredRodSolver::setTolerances,
             py::arg("rtol"), py::arg("atol"), "Relative and absolute tolerances for IDA.")
        .def("set_max_step",            &FredRodSolver::setMaxStep,
             py::arg("hmax"), "Maximum IDA step size [s] (IDASetMaxStep). <0 = no limit.")
        .def("set_init_step",           &FredRodSolver::setInitStep,
             py::arg("hinit"), "Initial IDA step size [s] (IDASetInitStep). <0 = auto.")
        .def("set_max_nonlin_iters",    &FredRodSolver::setMaxNonlinIters,
             py::arg("n"), "Max Newton iterations per IDA step (IDASetMaxNonlinIters, default 4).")
        .def("set_max_order",           &FredRodSolver::setMaxOrd,
             py::arg("order"),
             "Maximum BDF order (IDASetMaxOrd, default IDA max 5). Set to 1 to force\n"
             "backward Euler; can reduce Newton convergence failures at the cost of\n"
             "some local accuracy.")
        .def("set_verbosity",           &FredRodSolver::setVerbosity,
             py::arg("level"),
             "Console diagnostics level (0=silent, default 3=base time-loop status line).")
        .def("set_hot_start",           &FredRodSolver::setHotStart,
             py::arg("enable"),
             "Enable (True) or disable (False, default) a hot steady-state start:\n"
             "before the real time integration, march the rod at its t=0 boundary\n"
             "conditions (irradiation physics off) until temperatures/stresses stop\n"
             "changing, then use that converged state as t=0 instead of a cold start\n"
             "at set_initial_temperature(). Must be called before run(); ignored on\n"
             "checkpoint/snapshot restarts.")
        .def("set_enable_heat_conduction", &FredRodSolver::setEnableHeatConduction,
             py::arg("enable"),
             "Enable (True) or disable (False) the heat conduction block.\n"
             "When disabled, all temperatures are pinned to set_initial_temperature().")
        .def("set_enable_stress_strain",   &FredRodSolver::setEnableStressStrain,
             py::arg("enable"),
             "Enable (True) or disable (False) the stress-strain block.\n"
             "When disabled, all mechanical state variables are pinned to zero.")
        .def("run", &FredRodSolver::run,
             py::arg("tend"), py::arg("dtout"), py::arg("all_steps") = false,
             py::arg("threads") = 1,
             "Run simulation from t=0 to tend [s], output every dtout [s]. "
             "all_steps=True records every internal IDA step. "
             "threads: number of OpenMP threads for the per-axial-layer "
             "residual loop (default 1 = serial); silently clamped to 1 if "
             "a Python-defined material subclass is in use.")
        .def("time_points",           &FredRodSolver::timePoints, "Output time points [s].")
        .def("temperatures",          [](const FredRodSolver& s) {
                 auto& v  = s.temperatures();
                 auto& tp = s.timePoints();
                 int ns     = (int)tp.size();
                 int stride = (ns > 0 && !v.empty()) ? (int)(v.size() / ns) : 0;
                 return py::array_t<double>({ns, stride}, v.data());
             }, "Temperatures at all output times — shape (n_times, nz*(nf+nc)) [K].")
        .def("peak_fuel_temperature",  &FredRodSolver::peakFuelTemperature,
             "Peak fuel temperature at each output time [K].")
        .def("clad_outer_hoop_stress", &FredRodSolver::cladOuterHoopStress,
             "Outer cladding hoop stress [MPa] at each output time (axial average).")
        .def("clad_outer_radius",      &FredRodSolver::cladOuterRadius,
             "Outer cladding radius [m] at each output time (axial average).")
        .def("gap_width",              &FredRodSolver::gapWidth,
             "Radial gap width [m] at each output time (axial average).")
        .def("contact_pressure",       &FredRodSolver::contactPressure,
             "Pellet-cladding contact pressure [MPa] at each output time (axial average).")
        .def("fuel_outer_radius",      &FredRodSolver::fuelOuterRadius,
             "Fuel outer radius [m] at each output time (axial average).")
        .def("clad_inner_radius",      &FredRodSolver::cladInnerRadius,
             "Cladding inner radius [m] at each output time (axial average).")
        .def("neq_total",              &FredRodSolver::neqTotal,
             "Total number of IDA equations (0 for quasi-static path).")
        .def("y_out",  [](const FredRodSolver& s) {
                 auto& v  = s.yOut();
                 auto& tp = s.timePoints();
                 int ns  = (int)tp.size();
                 int neq = (ns > 0 && !v.empty()) ? (int)(v.size() / ns) : 0;
                 return py::array_t<double>({ns, neq}, v.data());
             }, "IDA y-vector at each output step — shape (n_steps, neq). "
                "Empty array if quasi-static path was used (heat off).")
        .def("yp_out", [](const FredRodSolver& s) {
                 auto& v  = s.ypOut();
                 auto& tp = s.timePoints();
                 int ns  = (int)tp.size();
                 int neq = (ns > 0 && !v.empty()) ? (int)(v.size() / ns) : 0;
                 return py::array_t<double>({ns, neq}, v.data());
             }, "IDA yp-vector at each output step — shape (n_steps, neq). "
                "Empty array if quasi-static path was used (heat off).")
        // Checkpoint API (fault recovery — single overwriting file)
        .def("set_checkpoint_prefix", &FredRodSolver::setCheckpointPrefix,
             py::arg("prefix"),
             "Enable fault-recovery checkpointing. Prefix for '<prefix>_checkpoint.chk'.\n"
             "Checkpoint is overwritten at every output step. Empty string disables (default).")
        .def("load_checkpoint",       &FredRodSolver::loadCheckpoint,
             py::arg("filename"),
             "Load a checkpoint file. Call before run(); simulation resumes at the saved time.")
        .def("restart_time",          &FredRodSolver::restartTime,
             "Return the restart time loaded from the last load_checkpoint/load_snapshot (<0 if none).")
        // Snapshot API (permanent restart states — time resets to 0 on load)
        .def("set_snapshot_prefix",
             [](FredRodSolver& s, const std::string& prefix,
                std::vector<double> times) {
                 s.setSnapshotPrefix(prefix, std::move(times));
             },
             py::arg("prefix"), py::arg("snapshot_times") = std::vector<double>{},
             "Enable snapshot saving. Prefix for '<prefix>_frameN.snapshot' files.\n"
             "snapshot_times: additional times (besides end-of-run) to save snapshots.\n"
             "The final simulation state is always saved as a snapshot.")
        .def("load_snapshot",         &FredRodSolver::loadSnapshot,
             py::arg("filename"),
             "Load a snapshot file. Physical state is restored but simulation time resets to 0.\n"
             "Use for starting a new transient from a pre-irradiated fuel state.");
}
