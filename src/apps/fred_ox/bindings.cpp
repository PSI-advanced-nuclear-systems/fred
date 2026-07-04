#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "apps/fred_ox/fuelpelletmaterial/FredOxMOX.hpp"
#include "apps/fred_ox/claddingmaterial/FredOxAIM1.hpp"
#include "apps/fred_ox/claddingmaterial/FredOxT91.hpp"
#include "apps/fred_ox/gapmaterial/FredOxGapMaterial.hpp"
#include "apps/fred_ox/FredOxSolver.hpp"
#include "platform/FuelPelletMaterial.hpp"
#include "platform/CladdingMaterial.hpp"
#include "platform/GapMaterial.hpp"

namespace py = pybind11;
using namespace fred;

void bind_fred_ox(py::module_& m) {

    py::class_<FredOxMOX, FuelPelletMaterial>(m, "FredOxMOX")
        .def(py::init<double, double, double>(),
             py::arg("pu_content"), py::arg("rof0"), py::arg("sto0") = 1.97,
             "MOX fuel for FRED-OX with burnup-dependent thermal conductivity.\n"
             "pu_content: Pu mole fraction [-]\n"
             "rof0: as-fabricated density [kg/m3]\n"
             "sto0: initial stoichiometry (default 1.97)")
        .def("set_burnup",       &FredOxMOX::setBurnup,        py::arg("bup_MWdkgU"))
        .def("set_stoichiometry",&FredOxMOX::setStoichiometry,  py::arg("sto"))
        .def("burnup",           &FredOxMOX::burnup)
        .def("pu_content",       &FredOxMOX::puContent)
        .def("thermal_conductivity",    &FredOxMOX::thermalConductivity,    py::arg("T"))
        .def("heat_capacity",           &FredOxMOX::heatCapacity,           py::arg("T"))
        .def("thermal_expansion_strain",&FredOxMOX::thermalExpansionStrain, py::arg("T"))
        .def("youngs_modulus",          &FredOxMOX::youngsModulus,          py::arg("T"), py::arg("density"))
        .def("poisson_ratio",           &FredOxMOX::poissonRatio)
        .def("reference_density",       &FredOxMOX::referenceDensity)
        .def("theoretical_density",     &FredOxMOX::theoreticalDensity)
        .def("melting_temperature",     &FredOxMOX::meltingTemperature);

    py::class_<FredOxAIM1, CladdingMaterial>(m, "FredOxAIM1")
        .def(py::init<double>(), py::arg("reference_density") = 7900.0,
             "AIM1 cladding for FRED-OX with irradiation creep model (Luzzi 2013).")
        .def("thermal_conductivity",    &FredOxAIM1::thermalConductivity,    py::arg("T"))
        .def("heat_capacity",           &FredOxAIM1::heatCapacity,           py::arg("T"))
        .def("thermal_expansion_strain",&FredOxAIM1::thermalExpansionStrain, py::arg("T"))
        .def("youngs_modulus",          &FredOxAIM1::youngsModulus,          py::arg("T"))
        .def("poisson_ratio",           &FredOxAIM1::poissonRatio)
        .def("meyer_hardness",          &FredOxAIM1::meyerHardness,          py::arg("T"))
        .def("reference_density",       &FredOxAIM1::referenceDensity)
        .def("creep_rate",              &FredOxAIM1::creepRate,
             py::arg("T"), py::arg("sigma_MPa"),
             "Effective creep rate [1/s] at T [K] and effective stress sigma [MPa].");

    py::class_<FredOxT91, CladdingMaterial>(m, "FredOxT91")
        .def(py::init<double>(), py::arg("reference_density") = 7750.0,
             "T91 cladding for FRED-OX (no irradiation creep — zero creep rate).")
        .def("thermal_conductivity",    &FredOxT91::thermalConductivity,    py::arg("T"))
        .def("heat_capacity",           &FredOxT91::heatCapacity,           py::arg("T"))
        .def("thermal_expansion_strain",&FredOxT91::thermalExpansionStrain, py::arg("T"))
        .def("youngs_modulus",          &FredOxT91::youngsModulus,          py::arg("T"))
        .def("poisson_ratio",           &FredOxT91::poissonRatio)
        .def("meyer_hardness",          &FredOxT91::meyerHardness,          py::arg("T"))
        .def("reference_density",       &FredOxT91::referenceDensity);

    py::class_<FredOxGapMaterial, GapMaterial>(m, "FredOxGapMaterial")
        .def(py::init<>(),
             "Gap material for FRED-OX with fission gas mixture (He/Xe/Kr) conductivity.")
        .def("set_gas_inventory",       &FredOxGapMaterial::setGasInventory,      py::arg("mu0_mol"))
        .def("set_fission_gas_release", &FredOxGapMaterial::setFissionGasRelease, py::arg("fgrel_mol"))
        .def("set_gas_pressure",        &FredOxGapMaterial::setGasPressure,       py::arg("gpres_MPa"))
        .def("set_burnup",              &FredOxGapMaterial::setBurnup,            py::arg("bup_MWdkgU"))
        .def("set_linear_power",        &FredOxGapMaterial::setLinearPower,       py::arg("ql_Wm"))
        .def("set_enable_relocation",   &FredOxGapMaterial::setEnableRelocation,  py::arg("enable"))
        .def("gap_conductivity",        &FredOxGapMaterial::gapConductivity,      py::arg("T"));

    py::class_<FredOxSolver>(m, "FredOxSolver")
        .def(py::init<const FuelRodGeometry&, FredOxMOX&,
                      const CladdingMaterial&, FredOxGapMaterial&>(),
             py::keep_alive<1, 2>(),  // solver keeps geometry alive
             py::keep_alive<1, 3>(),  // solver keeps mox alive
             py::keep_alive<1, 4>(),  // solver keeps clad alive
             py::keep_alive<1, 5>(),  // solver keeps gap alive
             py::arg("geometry"), py::arg("mox"), py::arg("clad"), py::arg("gap"),
             "FRED-OX solver: burnup accumulation, fission gas, swelling, burnup-dependent MOX.")
        .def("set_power_density_history", &FredOxSolver::setPowerDensityHistory,
             py::arg("times"), py::arg("qqv_W_m3"),
             "Uniform power density [W/m3] vs time (broadcast to all layers).")
        .def("set_power_density_history_per_layer",
             &FredOxSolver::setPowerDensityHistoryPerLayer,
             py::arg("times"), py::arg("qqv_per_layer"),
             "Per-layer power density. qqv_per_layer[j] = power vs time for layer j.")
        .def("set_coolant_temperature", &FredOxSolver::setCoolantTemperature,
             py::arg("times"), py::arg("T_K"),
             "Outer cladding temperature vs time (uniform for all layers).")
        .def("set_coolant_temperature_per_layer",
             &FredOxSolver::setCoolantTemperaturePerLayer,
             py::arg("times"), py::arg("T_per_layer"),
             "Per-layer outer cladding temperature.")
        .def("set_plenum_temperature_history", &FredOxSolver::setPlenumTemperatureHistory,
             py::arg("times"), py::arg("T_K"),
             "Gas plenum temperature vs time [K].")
        .def("set_coolant_pressure",     &FredOxSolver::setCoolantPressure,
             py::arg("pcool_MPa"), "Constant coolant pressure [MPa].")
        .def("set_coolant_pressure_history", &FredOxSolver::setCoolantPressureHistory,
             py::arg("times"), py::arg("pcool_MPa"),
             "Time-varying coolant pressure history [MPa].")
        .def("set_initial_temperature",  &FredOxSolver::setInitialTemperature,  py::arg("T0_K"))
        .def("set_initial_gas_pressure", &FredOxSolver::setInitialGasPressure,  py::arg("gpres0_MPa"))
        .def("set_swelling_multiplier",  &FredOxSolver::setSwellingMultiplier,
             py::arg("fswelmlt"), "Fuel swelling multiplier (1.0 = MATPRO full; 0.8 = legacy 10d case).")
        .def("set_tolerances",           &FredOxSolver::setTolerances,          py::arg("rtol"), py::arg("atol"))
        .def("set_max_step",             &FredOxSolver::setMaxStep,             py::arg("hmax"),
             "Maximum IDA step size [s] (IDASetMaxStep). <0 = no limit.")
        .def("set_init_step",            &FredOxSolver::setInitStep,            py::arg("hinit"),
             "Initial IDA step size [s] (IDASetInitStep). <0 = auto.")
        .def("set_max_nonlin_iters",     &FredOxSolver::setMaxNonlinIters,      py::arg("n"),
             "Max Newton iterations per IDA step (IDASetMaxNonlinIters, default 4).")
        .def("set_max_order",            &FredOxSolver::setMaxOrd,              py::arg("order"),
             "Maximum BDF order (IDASetMaxOrd, default IDA max 5). Set to 1 to force\n"
             "backward Euler; can reduce Newton convergence failures at the cost of\n"
             "some local accuracy.")
        .def("set_verbosity",            &FredOxSolver::setVerbosity,           py::arg("level"),
             "Console diagnostics level (0=silent, default 3=base time-loop status line).")
        .def("set_hot_start",            &FredOxSolver::setHotStart,            py::arg("enable"),
             "Enable (True) or disable (False, default) a hot steady-state start:\n"
             "before the real time integration, march the rod at its t=0 boundary\n"
             "conditions (irradiation physics off) until temperatures/stresses stop\n"
             "changing, then use that converged state as t=0 instead of a cold start\n"
             "at set_initial_temperature(). Must be called before run(); ignored on\n"
             "checkpoint/snapshot restarts.")
        .def("set_enable_heat_conduction", &FredOxSolver::setEnableHeatConduction, py::arg("enable"))
        .def("set_enable_stress_strain",   &FredOxSolver::setEnableStressStrain,   py::arg("enable"))
        .def("set_output_file", &FredOxSolver::setOutputFile,
             py::arg("filename"),
             "Path to HDF5 output file. When set, each output step is written to the file\n"
             "immediately (crash-safe) and in-memory vectors are trimmed to 1 entry.")
        .def("run", &FredOxSolver::run,
             py::arg("tend"), py::arg("dtout"), py::arg("all_steps") = false,
             py::arg("threads") = 1,
             "Run to tend [s], output every dtout [s]. all_steps=True records every IDA step. "
             "threads: number of OpenMP threads for the per-axial-layer residual loop "
             "(default 1 = serial).")
        .def("time_points",          &FredOxSolver::timePoints)
        .def("temperatures",         [](const FredOxSolver& s) {
                 auto& v  = s.temperatures();
                 auto& tp = s.timePoints();
                 int ns     = (int)tp.size();
                 int stride = (ns > 0 && !v.empty()) ? (int)(v.size() / ns) : 0;
                 return py::array_t<double>({ns, stride}, v.data());
             }, "Temperatures shape [n_times, nz*(nf+nc)] [K].")
        .def("peak_fuel_temperature",  &FredOxSolver::peakFuelTemperature)
        .def("gas_pressure",           &FredOxSolver::gasPressure,  "Internal gas pressure [MPa] vs time.")
        .def("fg_generated",           &FredOxSolver::fgGenerated,  "Fission gas generated [mol] vs time.")
        .def("fg_released",            &FredOxSolver::fgReleased,   "Fission gas released [mol] vs time.")
        .def("gap_width",              &FredOxSolver::gapWidth,     "Axial-avg gap width [m] vs time.")
        .def("burnup",                 &FredOxSolver::burnup,       "Axial-avg burnup [MWd/kgU] vs time.")
        .def("neq_total",              &FredOxSolver::neqTotal,     "Total number of IDA equations.")
        .def("y_out",  [](const FredOxSolver& s) {
                 auto& v  = s.yOut();
                 auto& tp = s.timePoints();
                 int ns  = (int)tp.size();
                 int neq = (ns > 0 && !v.empty()) ? (int)(v.size() / ns) : 0;
                 return py::array_t<double>({ns, neq}, v.data());
             }, "IDA y-vector at each output step — shape (n_steps, neq).")
        .def("yp_out", [](const FredOxSolver& s) {
                 auto& v  = s.ypOut();
                 auto& tp = s.timePoints();
                 int ns  = (int)tp.size();
                 int neq = (ns > 0 && !v.empty()) ? (int)(v.size() / ns) : 0;
                 return py::array_t<double>({ns, neq}, v.data());
             }, "IDA yp-vector at each output step — shape (n_steps, neq).")
        // Checkpoint API (fault recovery — single overwriting file)
        .def("set_checkpoint_prefix", &FredOxSolver::setCheckpointPrefix,
             py::arg("prefix"),
             "Enable fault-recovery checkpointing. Prefix for '<prefix>_checkpoint.chk'.\n"
             "Checkpoint is overwritten at every output step. Empty string disables (default).")
        .def("load_checkpoint",       &FredOxSolver::loadCheckpoint,
             py::arg("filename"),
             "Load a checkpoint file. Call before run(); simulation resumes at the saved time.")
        .def("restart_time",          &FredOxSolver::restartTime,
             "Return the restart time loaded from the last load_checkpoint/load_snapshot (<0 if none).")
        // Snapshot API (permanent restart states — time resets to 0 on load)
        .def("set_snapshot_prefix",
             [](FredOxSolver& s, const std::string& prefix,
                std::vector<double> times) {
                 s.setSnapshotPrefix(prefix, std::move(times));
             },
             py::arg("prefix"), py::arg("snapshot_times") = std::vector<double>{},
             "Enable snapshot saving. Prefix for '<prefix>_frameN.snapshot' files.\n"
             "snapshot_times: additional times (besides end-of-run) to save snapshots.\n"
             "The final simulation state is always saved as a snapshot.")
        .def("load_snapshot",         &FredOxSolver::loadSnapshot,
             py::arg("filename"),
             "Load a snapshot file. Physical state is restored but simulation time resets to 0.\n"
             "Use for starting a new transient from a pre-irradiated fuel state.");
}
