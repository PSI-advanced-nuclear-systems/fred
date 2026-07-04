#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "apps/fred_m_na/fuelpelletmaterial/UPuZr.hpp"
#include "apps/fred_m_na/claddingmaterial/HT9.hpp"
#include "apps/fred_m_na/gapmaterial/Sodium.hpp"
#include "apps/fred_m_na/FredMNaSolver.hpp"
#include "apps/fred_m_na/FredMNaGrsis.hpp"
#include "apps/fred_m_na/FredMNaSubchannelMode.hpp"
#include "platform/FuelPelletMaterial.hpp"
#include "platform/CladdingMaterial.hpp"

namespace py = pybind11;
using namespace fred;
using fred::ConductivityModel;

void bind_fred_m_na(py::module_& m) {

    py::enum_<ConductivityModel>(m, "ConductivityModel",
        "Irradiation correction model for U-Pu-Zr thermal conductivity (Flamb.for).")
        .value("DetailedNaSodium", ConductivityModel::DetailedNaSodium,
               "f=1: Maxwell-Eucken with Na infiltration (Karahan 2009 + MFUEL). "
               "Requires per-node porosity split and sodium infiltration fraction.")
        .value("EmpiricalBurnup", ConductivityModel::EmpiricalBurnup,
               "f=2: Piecewise empirical burnup degradation. Requires only bup_FIMA.")
        .value("EsfrSimple", ConductivityModel::EsfrSimple,
               "f=3: ESFR-SIMPLE sigmoid fit in burnup (at%). Requires only bup_FIMA.")
        .export_values();

    py::enum_<GrsisDataMode>(m, "GrsisDataMode",
        "GRSIS bubble model parameter set (FEAST or original GRSIS).")
        .value("FEAST", GrsisDataMode::FEAST,
               "FEAST-METAL calibrated parameters (default in Fortran).")
        .value("GRSIS", GrsisDataMode::GRSIS,
               "Original GRSIS model parameters (Lee 1999).")
        .export_values();

    py::enum_<SodiumMode>(m, "SodiumMode",
        "Sodium gap conductance mode (igap in Gaphtc.for).")
        .value("TDependent", SodiumMode::TDependent,
               "igap=1: T-dependent polynomial k_Na(T) [W/(m·K)].")
        .value("Constant",   SodiumMode::Constant,
               "igap=2: constant 62.9 W/(m·K) with conductance caps [1e5,1e6] W/(m²·K).")
        .export_values();

    py::class_<UPuZr, FuelPelletMaterial>(m, "UPuZr")
        .def(py::init<double, double, double>(),
             py::arg("pu_weight_frac"), py::arg("zr_weight_frac"),
             py::arg("reference_density") = -1.0,
             "U-Pu-Zr metallic fuel pellet (Karahan 2009 correlations).\n"
             "pu_weight_frac: Pu weight fraction [-] (e.g. 0.19 for 19 wt%)\n"
             "zr_weight_frac: Zr weight fraction [-] (e.g. 0.10 for 10 wt%)\n"
             "reference_density: as-fabricated density [kg/m3]; <=0 → 75% of theoretical")
        .def("thermal_conductivity",    &UPuZr::thermalConductivity,    py::arg("T"))
        .def("heat_capacity",           &UPuZr::heatCapacity,           py::arg("T"))
        .def("thermal_expansion_strain",&UPuZr::thermalExpansionStrain, py::arg("T"))
        .def("youngs_modulus",          &UPuZr::youngsModulus,          py::arg("T"), py::arg("density"))
        .def("poisson_ratio",           &UPuZr::poissonRatio)
        .def("reference_density",       &UPuZr::referenceDensity)
        .def("theoretical_density",     &UPuZr::theoreticalDensity)
        .def("melting_temperature",     &UPuZr::meltingTemperature,
             "Solidus temperature [K] for the as-fabricated Zr content (Mmelt.for 25-point table).")
        .def("solidus_temperature",     &UPuZr::solidusTemperature, py::arg("zr_wf"),
             "Solidus temperature [K] for arbitrary Zr weight fraction (nearest-neighbour lookup).")
        .def("pu_content",              &UPuZr::puContent, "Pu weight fraction [-]")
        .def("zr_content",              &UPuZr::zrContent, "Zr weight fraction [-]");

    py::class_<HT9, CladdingMaterial>(m, "HT9")
        .def(py::init<double>(), py::arg("reference_density") = 7750.0,
             "HT-9 ferritic-martensitic steel cladding (Hofmann 1985, Karahan 2007).")
        .def("thermal_conductivity",    &HT9::thermalConductivity,    py::arg("T"))
        .def("heat_capacity",           &HT9::heatCapacity,           py::arg("T"))
        .def("thermal_expansion_strain",&HT9::thermalExpansionStrain, py::arg("T"))
        .def("youngs_modulus",          &HT9::youngsModulus,          py::arg("T"))
        .def("poisson_ratio",           &HT9::poissonRatio)
        .def("meyer_hardness",          &HT9::meyerHardness,          py::arg("T"))
        .def("yield_stress",            &HT9::yieldStress,            py::arg("T"),
             "Cladding yield stress [Pa] (Csigy.for table interpolation).")
        .def("burst_stress",            &HT9::burstStress,
             py::arg("T"), py::arg("ssy0"),
             "Burst stress [Pa] = ssy0*(1.1 - 0.1*tanh((T_C-200)/200)) (Csigb.for, bug-corrected).")
        .def("reference_density",       &HT9::referenceDensity);

    py::enum_<FredMNaSubchannelMode::HtcCorrelation>(m, "HtcCorrelation",
        "Heat transfer coefficient correlation for sodium subchannel cooling.")
        .value("Mikityuk", FredMNaSubchannelMode::HtcCorrelation::kMikityuk,
               "Mikityuk liquid-metal rod-bundle correlation (default).")
        .value("Subbotin", FredMNaSubchannelMode::HtcCorrelation::kSubbotin,
               "Subbotin simplified liquid-metal correlation.")
        .export_values();

    py::class_<FredMNaSolver>(m, "FredMNaSolver")
        .def(py::init<const FuelRodGeometry&, UPuZr&, HT9&>(),
             py::keep_alive<1, 2>(),  // solver keeps geometry alive
             py::keep_alive<1, 3>(),  // solver keeps fuel alive
             py::keep_alive<1, 4>(),  // solver keeps clad alive
             py::arg("geometry"), py::arg("fuel"), py::arg("clad"),
             "FRED-M-Na solver: U-Pu-Zr fuel / HT-9 cladding / sodium fast reactor.\n"
             "Physics: burnup, FGR (Karahan empirical), fuel swelling (solid FP),\n"
             "sodium infiltration, Zr redistribution, cladding wastage.")
        .def("set_power_density_history", &FredMNaSolver::setPowerDensityHistory,
             py::arg("times"), py::arg("qqv_W_m3"))
        .def("set_power_density_history_per_layer",
             &FredMNaSolver::setPowerDensityHistoryPerLayer,
             py::arg("times"), py::arg("qqv_per_layer"),
             "Set time-varying, axially-varying power density [W/m3].\n"
             "times: vector of times [s]\n"
             "qqv_per_layer: list of nz vectors, each with len(times) values")
        .def("set_coolant_temperature",
             [](FredMNaSolver&, py::args, py::kwargs) {
                 throw py::attribute_error(
                     "FredMNaSolver uses subchannel calculation for coolant temperature and pressure\n"
                     "with a Robin boundary condition for the cladding outer surface heat flux.\n\n"
                     "Example: solver.set_coolant_channel(dhyd, xarea, flowr, T_inlet_times, T_inlet_vals)\n\n"
                     "See FredMNaSubchannelMode.cpp for implementation details.");
             })
        .def("set_coolant_temperature_per_layer",
             [](FredMNaSolver&, py::args, py::kwargs) {
                 throw py::attribute_error(
                     "FredMNaSolver uses subchannel calculation for coolant temperature and pressure\n"
                     "with a Robin boundary condition for the cladding outer surface heat flux.\n\n"
                     "Example: solver.set_coolant_channel(dhyd, xarea, flowr, T_inlet_times, T_inlet_vals)\n\n"
                     "See FredMNaSubchannelMode.cpp for implementation details.");
             })
        .def("set_coolant_channel",
             [](FredMNaSolver& s, double dhyd, double xarea, double flowr,
                std::vector<double> times, std::vector<double> vals,
                FredMNaSubchannelMode::HtcCorrelation corr) {
                 s.setCoolantChannel(dhyd, xarea, flowr, std::move(times), std::move(vals), corr);
             },
             py::arg("dhyd"), py::arg("xarea"), py::arg("flowr"),
             py::arg("T_inlet_times"), py::arg("T_inlet_vals"),
             py::arg("corr") = FredMNaSubchannelMode::HtcCorrelation::kMikityuk,
             "Set subchannel cooling: axial energy balance + Peclet-number HTC.\n"
             "dhyd: hydraulic diameter [m]\n"
             "xarea: flow area [m2]\n"
             "flowr: mass flow rate [kg/s]\n"
             "T_inlet_times: time vector [s]\n"
             "T_inlet_vals: inlet coolant temperature vector [K]\n"
             "corr: HtcCorrelation.Mikityuk (default) or HtcCorrelation.Subbotin")
        .def("set_plenum_temperature_history", &FredMNaSolver::setPlenumTemperatureHistory,
             py::arg("times"), py::arg("T_K"))
        .def("set_coolant_pressure",      &FredMNaSolver::setCoolantPressure,    py::arg("pcool_MPa"))
        .def("set_initial_temperature",   &FredMNaSolver::setInitialTemperature, py::arg("T0_K"))
        .def("set_conductivity_model",    &FredMNaSolver::setConductivityModel,
             py::arg("model"),
             "Select thermal conductivity correction model (default: DetailedNaSodium).")
        .def("set_initial_gas_pressure",  &FredMNaSolver::setInitialGasPressure,
             py::arg("p_MPa"),
             "Initial fill-gas pressure [MPa] at T_ref; sets the gas inventory (mu0).")
        .def("set_bond_pressure",         &FredMNaSolver::setBondPressure,
             py::arg("p_MPa"),
             "Legacy alias for set_initial_gas_pressure.")
        .def("set_step_size",              &FredMNaSolver::setStepSize,           py::arg("dt"),
             "Fixed internal time step [s] for the one-step backward-Euler integrator\n"
             "(distinct from dtout, the output-writing cadence passed to run()). <=0\n"
             "(default): use dtout itself as the step size. FRED-M-Na does not use\n"
             "SUNDIALS IDA (no adaptive step control) — this is the actual, literal\n"
             "step size used for every backward-Euler solve, matching legacy FRED-M's\n"
             "fixed dt=dtout scheme; backward Euler is unconditionally stable, so this\n"
             "is a truncation-accuracy/performance trade-off, not a stability one.")
        .def("set_verbosity",             &FredMNaSolver::setVerbosity,          py::arg("level"),
             "Console diagnostics level: 0=silent, default 3=per-output-step status\n"
             "line, 4=+per-fixed-step outer-Picard-sweep diagnostics, 5=+afterAcceptedStep\n"
             "wall-clock timing.")
        .def("set_hot_start",             &FredMNaSolver::setHotStart,           py::arg("enable"),
             "Enable (True) or disable (False, default) a hot steady-state start:\n"
             "before the real time integration, march the rod at its t=0 boundary\n"
             "conditions (irradiation physics off — GRSIS/Zr redistribution/cladding\n"
             "wastage/burnup/fission gas all skipped) until temperatures/stresses stop\n"
             "changing, then use that converged state as t=0 instead of a cold start\n"
             "at set_initial_temperature(). Must be called before run(); ignored on\n"
             "checkpoint/snapshot restarts.")
        .def("set_enable_heat_conduction",   &FredMNaSolver::setEnableHeatConduction,   py::arg("enable"))
        .def("set_enable_stress_strain",     &FredMNaSolver::setEnableStressStrain,     py::arg("enable"))
        .def("set_enable_zr_redistribution", &FredMNaSolver::setEnableZrRedistribution, py::arg("enable"))
        .def("set_enable_clad_wastage",      &FredMNaSolver::setEnableCladWastage,      py::arg("enable"))
        .def("set_enable_grsis",             &FredMNaSolver::setEnableGrsis,            py::arg("enable"),
             "Enable/disable GRSIS bubble swelling model (default: enabled).")
        .def("set_grsis_data_mode",          &FredMNaSolver::setGrsisDataMode,
             py::arg("mode"),
             "Select GRSIS parameter set: GrsisDataMode.FEAST (default) or .GRSIS.")
        .def("set_sodium_mode",              &FredMNaSolver::setSodiumMode,
             py::arg("mode"),
             "Select sodium gap conductance mode: SodiumMode.TDependent (default) or .Constant.")
        .def("set_output_file", &FredMNaSolver::setOutputFile,
             py::arg("filename"),
             "Path to HDF5 output file. When set, each output step is written to the file\n"
             "immediately (crash-safe) and in-memory vectors are trimmed to 1 entry.")
        .def("run", &FredMNaSolver::run,
             py::arg("tend"), py::arg("dtout"), py::arg("all_steps") = false,
             py::arg("threads") = 1,
             "threads: number of OpenMP threads for the per-axial-layer Newton "
             "solve (default 1 = serial).")
        .def("times",               &FredMNaSolver::times)
        .def("gas_pressure",        &FredMNaSolver::gasPressure)
        .def("fg_generated",        &FredMNaSolver::fgGen)
        .def("fg_released",         &FredMNaSolver::fgRel)
        .def("gap_width",           &FredMNaSolver::gapWidth)
        .def("burnup",              &FredMNaSolver::burnup)
        .def("peak_fuel_temperature",&FredMNaSolver::peakFuelTemperature)
        .def("clad_wastage",        &FredMNaSolver::cladWastage,
             "Max cladding wastage thickness [m] across layers at each output time.")
        .def("grsis_swelling_total", &FredMNaSolver::grisSwellingTotal,
             "Spatially-averaged total GRSIS swelling fraction [-] at each output time.")
        .def("grsis_swelling_open",  &FredMNaSolver::grisSwellingOpen,
             "Spatially-averaged open-bubble swelling fraction [-] at each output time.")
        .def("burst_margin",        &FredMNaSolver::burstMargin,
             "Max cladding burst criterion (sigma/sigma_burst) [-] at each output time.")
        .def("melt_margin",         &FredMNaSolver::meltMargin,
             "Max fuel melt criterion (T/T_solidus) [-] at each output time.")
        .def("coolant_temperature_per_layer",
             [](const FredMNaSolver& s) {
                 auto& v  = s.coolantTemperaturePerLayer();
                 auto& tp = s.times();
                 int ns  = (int)tp.size();
                 int nz  = (ns > 0 && !v.empty()) ? (int)(v.size() / ns) : 0;
                 return py::array_t<double>({ns, nz}, v.data());
             }, "Coolant temperature per layer shape [n_times, nz] [K].")
        .def("gap_width_per_layer",
             [](const FredMNaSolver& s) {
                 auto& v  = s.gapWidthPerLayer();
                 auto& tp = s.times();
                 int ns  = (int)tp.size();
                 int nz  = (ns > 0 && !v.empty()) ? (int)(v.size() / ns) : 0;
                 return py::array_t<double>({ns, nz}, v.data());
             }, "Gap width per layer shape [n_times, nz] [m].")
        .def("hgap_per_layer",
             [](const FredMNaSolver& s) {
                 auto& v  = s.hgapPerLayer();
                 auto& tp = s.times();
                 int ns  = (int)tp.size();
                 int nz  = (ns > 0 && !v.empty()) ? (int)(v.size() / ns) : 0;
                 return py::array_t<double>({ns, nz}, v.data());
             }, "Gap conductance per layer shape [n_times, nz] [W/(m2*K)].")
        .def("temperatures",        [](const FredMNaSolver& s) {
                 auto& v  = s.temperatures();
                 auto& tp = s.times();
                 int ns     = (int)tp.size();
                 int stride = (ns > 0 && !v.empty()) ? (int)(v.size() / ns) : 0;
                 return py::array_t<double>({ns, stride}, v.data());
             }, "Temperatures shape [n_times, nz*(nf+nc)] [K].")
        .def("neq_total", &FredMNaSolver::neqTotal)
        .def("y_out",  [](const FredMNaSolver& s) {
                 auto& v  = s.y_out();
                 auto& tp = s.times();
                 int ns  = (int)tp.size();
                 int neq = (ns > 0 && !v.empty()) ? (int)(v.size() / ns) : 0;
                 return py::array_t<double>({ns, neq}, v.data());
             })
        .def("yp_out", [](const FredMNaSolver& s) {
                 auto& v  = s.yp_out();
                 auto& tp = s.times();
                 int ns  = (int)tp.size();
                 int neq = (ns > 0 && !v.empty()) ? (int)(v.size() / ns) : 0;
                 return py::array_t<double>({ns, neq}, v.data());
             })
        // Checkpoint API (fault recovery — single overwriting file)
        .def("set_checkpoint_prefix", &FredMNaSolver::setCheckpointPrefix,
             py::arg("prefix"),
             "Enable fault-recovery checkpointing. Prefix for '<prefix>_checkpoint.chk'.\n"
             "Checkpoint is overwritten at every output step. Empty string disables (default).")
        .def("load_checkpoint",       &FredMNaSolver::loadCheckpoint,
             py::arg("filename"),
             "Load a checkpoint file. Call before run(); simulation resumes at the saved time.")
        .def("restart_time",          &FredMNaSolver::restartTime,
             "Return the restart time loaded from the last load_checkpoint/load_snapshot (<0 if none).")
        // Snapshot API (permanent restart states — time resets to 0 on load)
        .def("set_snapshot_prefix",
             [](FredMNaSolver& s, const std::string& prefix,
                std::vector<double> times) {
                 s.setSnapshotPrefix(prefix, std::move(times));
             },
             py::arg("prefix"), py::arg("snapshot_times") = std::vector<double>{},
             "Enable snapshot saving. Prefix for '<prefix>_frameN.snapshot' files.\n"
             "snapshot_times: additional times (besides end-of-run) to save snapshots.\n"
             "The final simulation state is always saved as a snapshot.")
        .def("load_snapshot",         &FredMNaSolver::loadSnapshot,
             py::arg("filename"),
             "Load a snapshot file. Physical state is restored but simulation time resets to 0.\n"
             "Accumulated irradiation state (burnup, Zr redistribution, etc.) is preserved;\n"
             "only the time coordinate resets to allow a fresh transient simulation.");
}
