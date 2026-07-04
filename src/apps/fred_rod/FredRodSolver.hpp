#pragma once
#include "FredRodResiduals.hpp"
#include "../../platform/FredIdaSolverBase.hpp"
#include "../../platform/GapPressureModel.hpp"
#include "../../platform/Constants.hpp"
#include "platform/GapMaterial.hpp"
#include <vector>
#include <memory>

namespace fred {

class FuelPelletMaterial;
class CladdingMaterial;

// \coderef{fred_rod_solver}
// -----------------------------------------------------------------------
// FredRodSolver — IDA-based time integration driver for FRED-ROD.
//
// Inherits shared SUNDIALS infrastructure from FredIdaSolverBase and adds:
//   • Gas pressure model management (GasBondPressure or prescribed)
//   • Quasi-static (heat-off) path
//   • Gap closure/reopening root detection with mechanical IC re-solve
//   • Result accessors specific to thermo-mechanical output
// -----------------------------------------------------------------------
class FredRodSolver : public FredIdaSolverBase {
public:
    FredRodSolver(const FuelRodGeometry&    geom,
                  const FuelPelletMaterial& fuel,
                  const CladdingMaterial&   clad,
                  const GapMaterial&        gap);

    // Set the initial fill-gas pressure [MPa] (must be called before run()).
    void setInitialGasPressure(double gpres0_MPa);

    // \coderef{fred_rod_run}
    // threads: number of OpenMP threads for the per-axial-layer residual
    // loop (default 1 = serial, matching pre-threading behavior). Values
    // >1 are silently clamped back to 1 if a Python-subclassed (hence
    // not thread-safe) material was supplied to this solver.
    void run(double tend, double dtout, bool all_steps = false, int threads = 1);

    int neqTotal() const override { return m_res.neqTotal(); }

    // ROD-specific result accessors.
    std::vector<double> peakFuelTemperature() const { return FredSolverBase::peakFuelTemperature(); }
    const std::vector<double>& cladOuterHoopStress() const { return m_sigh_outer_out; }
    const std::vector<double>& cladOuterRadius()     const { return m_rco_out; }
    const std::vector<double>& gapWidth()            const { return m_gap_out; }
    const std::vector<double>& contactPressure()     const { return m_pfc_out; }
    const std::vector<double>& fuelOuterRadius()     const { return m_rfo_out; }
    const std::vector<double>& cladInnerRadius()     const { return m_rci_out; }

protected:
    // FredSolverBase boundary-condition hooks.
    RodResiduals& rodResiduals()                                             override { return m_res; }
    void onPowerDensityPerLayerSet      (const std::vector<TimeTable>& tabs) override;
    void onCoolantTemperaturePerLayerSet(const std::vector<TimeTable>& tabs) override;
    void onCoolantPressureSet   (const TimeTable& tbl) override;
    void onInitialTemperatureSet(double T0)            override;
    void onPlenumTemperatureSet (const TimeTable& tbl) override;

    // FredSolverBase setup hooks.
    int      getSolverNeq()                            const override { return m_res.neqTotal(); }
    IDAResFn getSolverResFn()                          const override { return RodResiduals::idaResidual; }
    void*    getSolverUserData()                             override { return &m_res; }
    void     fillSolverIdArray(double* id, int neq)    const override;
    void     packSolverState  (double* y)              const override { m_res.pack(m_state, y); }

    // FredSolverBase runTimeLoop hooks.
    void unpackCurrentState()                                override;
    void storeOutput       (double t)                        override;
    void afterGapEventsHandled(double tret, double dtout)    override;

private:
    const FuelPelletMaterial& m_fuel;
    const CladdingMaterial&   m_clad;
    const GapMaterial&        m_gap;
    // Pressure model declared before m_res so it is constructed first.
    std::unique_ptr<GapPressureModel> m_pressure_model;
    FredRodResiduals m_res;
    RodState         m_state;

    // ROD-specific output vectors.
    std::vector<double> m_sigh_outer_out;
    std::vector<double> m_rco_out;
    std::vector<double> m_gap_out;
    std::vector<double> m_pfc_out;
    std::vector<double> m_rfo_out;
    std::vector<double> m_rci_out;

    // HDF5 dataset handles for ROD-specific quantities (int64_t = hid_t).
    struct RodH5Ctx {
        int64_t grp_mech = -1;
        int64_t ds_gap   = -1;
        int64_t ds_sigh  = -1;
        int64_t ds_rco   = -1;
        int64_t ds_pfc   = -1;
        int64_t ds_rfo   = -1;
        int64_t ds_rci   = -1;
    };
    RodH5Ctx m_h5rod;

    void openAppH5Datasets   () override;
    void appendAppH5Row      () override;
    void trimAppOutputVectors() override;
    void closeAppH5Datasets  () override;

    void initState();
};

} // namespace fred
