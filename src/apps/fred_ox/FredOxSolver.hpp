#pragma once
#include "FredOxResiduals.hpp"
#include "FredOxState.hpp"
#include "platform/FredIdaSolverBase.hpp"
#include "platform/Constants.hpp"
#include <vector>
#include <memory>
#include <functional>

namespace fred {

class FredOxMOX;
class CladdingMaterial;

// Per-axial-layer power density table: nz TimeTables, one per layer.
struct LayerPowerTable {
    std::vector<TimeTable> layers;
    int nz = 0;
};

// \coderef{fred_ox_solver}
// -----------------------------------------------------------------------
// FredOxSolver — IDA-based time integration driver for FRED-OX.
//
// Inherits shared SUNDIALS infrastructure from FredIdaSolverBase and adds:
//   • Per-layer power density and coolant temperature time tables
//   • Fission gas production / release ODEs
//   • Fuel swelling ODEs (MATPRO, optional multiplier)
//   • MOX burnup-dependent thermal conductivity
//   • Gap conductance with fission gas mixture
// -----------------------------------------------------------------------
class FredOxSolver : public FredIdaSolverBase {
public:
    FredOxSolver(const FuelRodGeometry& geom,
                 FredOxMOX&             mox,
                 const CladdingMaterial& clad,
                 FredOxGapMaterial&     gap);

    // --- OX-only boundary condition setters ---
    void setInitialGasPressure(double gpres0_MPa);
    void setSwellingMultiplier(double fswelmlt);

    // threads: number of OpenMP threads for the per-axial-layer residual
    // loop (default 1 = serial, matching pre-threading behavior).
    void run(double tend, double dtout, bool all_steps = false, int threads = 1);

    int neqTotal() const override { return m_res.neqTotal(); }

    // --- Result accessors ---
    const std::vector<double>& gasPressure() const { return m_gpres_out; }
    const std::vector<double>& fgGenerated() const { return m_fggen_out; }
    const std::vector<double>& fgReleased()  const { return m_fgrel_out; }
    const std::vector<double>& gapWidth()    const { return m_gap_out; }
    const std::vector<double>& burnup()      const { return m_bup_out; }

protected:
    // FredSolverBase boundary-condition hooks.
    RodResiduals& rodResiduals()                                             override { return m_res; }
    void onPowerDensityPerLayerSet      (const std::vector<TimeTable>& tabs) override;
    void onCoolantTemperaturePerLayerSet(const std::vector<TimeTable>& tabs) override;
    void onCoolantPressureSet   (const TimeTable& tbl) override;
    void onInitialTemperatureSet(double T0)            override;
    void onPlenumTemperatureSet (const TimeTable& tbl) override;

    // FredSolverBase setup hooks.
    int      getSolverNeq()                          const override { return m_res.neqTotal(); }
    IDAResFn getSolverResFn()                        const override { return RodResiduals::idaResidual; }
    void*    getSolverUserData()                           override { return &m_res; }
    void     fillSolverIdArray(double* id, int neq)  const override { m_res.fillIdArray(id); }
    void     packSolverState  (double* y)            const override { m_res.pack(m_state, y); }

    // FredSolverBase runTimeLoop hooks.
    void unpackCurrentState()                              override;
    void storeOutput       (double t)                      override;
    void logStepOutput     (double tret, double dt_next)   override;

private:
    FredOxMOX&             m_mox;
    const CladdingMaterial& m_clad;
    FredOxGapMaterial&     m_gap;
    FredOxResiduals        m_res;
    FredOxRodState         m_state;

    std::vector<TimeTable> m_powerTab;
    std::vector<TimeTable> m_coolantTab;
    TimeTable              m_plenumTab;

    // OX-specific output vectors.
    std::vector<double> m_gpres_out, m_fggen_out, m_fgrel_out;
    std::vector<double> m_gap_out, m_bup_out;

    // HDF5 dataset handles for OX-specific quantities (int64_t = hid_t).
    struct OxH5Ctx {
        int64_t grp_burnup = -1;
        int64_t ds_gap     = -1;
        int64_t ds_gpres   = -1;
        int64_t ds_fggen   = -1;
        int64_t ds_fgrel   = -1;
        int64_t ds_bup     = -1;
    };
    OxH5Ctx m_h5ox;

    void openAppH5Datasets   () override;
    void appendAppH5Row      () override;
    void trimAppOutputVectors() override;
    void closeAppH5Datasets  () override;

    void initState();
};

} // namespace fred
