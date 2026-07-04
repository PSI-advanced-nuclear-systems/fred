#include "RodResiduals.hpp"
#include "FuelPelletMaterial.hpp"
#include "CladdingMaterial.hpp"
#include <nvector/nvector_serial.h>
#include <algorithm>
#include <limits>
#include <vector>
#include <omp.h>

namespace fred {

RodResiduals::RodResiduals(const FuelRodGeometry&    geom,
                            const FuelPelletMaterial& fuel,
                            const CladdingMaterial&   clad,
                            const GapMaterial&        gap_mat,
                            double                    coolant_pressure_MPa,
                            int                       neq_irr)
    : m_geom    (geom),
      m_clad    (clad),
      m_heat    (geom, fuel, clad, gap_mat),
      m_mech    (geom, fuel, clad),
      m_neq_th  (m_heat.neq()),
      m_neq_irr (neq_irr),
      m_neq_mech(m_mech.neq()),
      m_neq_j   (m_heat.neq() + neq_irr + m_mech.neq()),
      m_pcool   (coolant_pressure_MPa),
      m_gapMgr  (geom),
      m_materials_thread_safe(fuel.isThreadSafe() && clad.isThreadSafe() && gap_mat.isThreadSafe())
{}

// ---------------------------------------------------------------------------
// Per-layer BC lookup helpers
// ---------------------------------------------------------------------------
double RodResiduals::layerPower(int j, double t) const {
    return (!m_layerPowerFns.empty() && j < (int)m_layerPowerFns.size())
           ? m_layerPowerFns[j](t)
           : (m_powerFn ? m_powerFn(t) : 0.0);
}

double RodResiduals::layerCoolant(int j, double t) const {
    if (m_subchannel != nullptr) return m_subchannel->T_co(j);
    return (!m_layerCoolantFns.empty() && j < (int)m_layerCoolantFns.size())
           ? m_layerCoolantFns[j](t)
           : (m_coolantTFn ? m_coolantTFn(t) : m_T_init);
}

double RodResiduals::layerPcool(double t) const {
    return m_pcoolFn ? m_pcoolFn(t) : m_pcool;
}

double RodResiduals::layerHTC(int j) const {
    if (m_subchannel != nullptr) return m_subchannel->htc(j);
    return std::numeric_limits<double>::infinity();
}

// ---------------------------------------------------------------------------
// Default freeze: zero all thermal residuals (IDA leaves vars unchanged).
// ROD overrides to explicitly pin AEs to as-fabricated values.
// ---------------------------------------------------------------------------
void RodResiduals::freezeThermalResiduals(int /*j*/, const double* /*yj*/, double* rj) const {
    for (int k = 0; k < m_neq_th; ++k) rj[k] = 0.0;
}

// ---------------------------------------------------------------------------
// Shared IDA residual entry point (NVI)
// ---------------------------------------------------------------------------
int RodResiduals::computeResiduals(double t, const double* y, const double* yp, double* r) {
    unpackState(y, yp);

    // Hot-start pseudo-time march (see setBcTimeOverride): all boundary-
    // condition lookups downstream of this point use t_bc instead of the
    // integrator's own advancing t, so power/coolant/plenum stay pinned at
    // t_fixed while the DAE itself is still marched forward in real pseudo-
    // time. Subclasses' computeGlobalResiduals/computeThermalResiduals
    // overrides inherit this automatically since they only ever see t_bc.
    const double t_bc = bcTime(t);

    const int    ng    = globalOffset();
    const double gpres = globalPressure();
    const double pcool = layerPcool(t_bc);

    computeGlobalResiduals(t_bc, yp, r);

    if (m_subchannel != nullptr) {
        std::vector<double> qqv_all(m_geom.nz);
        for (int j = 0; j < m_geom.nz; ++j)
            qqv_all[j] = layerPower(j, t_bc);
        m_subchannel->updateCoolantField(t_bc, qqv_all.data(), m_geom.dz0.data());
    }

    const int nz = m_geom.nz;
    prepareThreadLocalState(m_num_threads);
    #pragma omp parallel for schedule(static) num_threads(m_num_threads) if(m_num_threads > 1)
    for (int j = 0; j < nz; ++j) {
        prepareLayer(j);
        const double* yj  = y  + ng + j * m_neq_j;
        const double* ypj = yp + ng + j * m_neq_j;
        double*       rj  = r  + ng + j * m_neq_j;

        if (m_heat_on)
            computeThermalResiduals(j, t_bc, yj, ypj, rj);
        else
            freezeThermalResiduals(j, yj, rj);

        if (m_irr_on)
            computeIrradiationResiduals(j, ypj, rj);
        else
            freezeIrradiationResiduals(j, ypj, rj);

        if (m_mech_on) {
            std::vector<double>& r_mech = threadLocalMechBuffer();
            m_mech.computeResiduals(layerState(j), gpres, pcool, r_mech);
            std::copy(r_mech.begin(), r_mech.end(), rj + m_neq_th + m_neq_irr);
        } else {
            // Pin each mechanical slot to its initial value (0) so the Jacobian
            // has a well-conditioned diagonal block instead of all-zero columns,
            // which would make IDA's dense linear solver singular.
            for (int k = m_neq_th + m_neq_irr; k < m_neq_j; ++k) rj[k] = yj[k];
        }
    }
    return 0;
}

// ---------------------------------------------------------------------------
// IDA static callbacks
// ---------------------------------------------------------------------------
int RodResiduals::idaResidual(double t, N_Vector y, N_Vector yp, N_Vector r, void* user_data) {
    auto* self = static_cast<RodResiduals*>(user_data);
    return self->computeResiduals(
        t,
        N_VGetArrayPointer(y),
        N_VGetArrayPointer(yp),
        N_VGetArrayPointer(r));
}

// gapRoot — IDARootFn shared by all RodResiduals subclasses.
//
// Uses the y-vector layout shared across ROD and OX:
//   gpres at y[globalOffset()-1]
//   per-layer block at y[globalOffset() + j*neq_j]:
//     [0] = qqv/rci   [1] = gap   [2] = pfc   [3] = hgap   [4..] = T
//
// Open  gap: gout[j] = (gap - thresh) / thresh    (crosses 0 downward at closure)
// Closed gap: gout[j] = (gpres - pfc) / (pfc + gpres + eps)  (crosses 0 upward at reopening)
int RodResiduals::gapRoot(double /*t*/, N_Vector y, N_Vector /*yp*/,
                           double* gout, void* user_data)
{
    auto*        self  = static_cast<RodResiduals*>(user_data);
    const double* yv   = N_VGetArrayPointer(y);
    const int    ng    = self->globalOffset();
    const double gpres = yv[ng - 1];
    const double thresh = self->m_geom.ruff + self->m_geom.rufc;

    for (int j = 0; j < self->m_geom.nz; ++j) {
        const double* yj = yv + ng + j * self->m_neq_j;
        const double  gap = yj[1];
        const double  pfc = yj[2];
        if (self->m_gapMgr.isGapOpen(j)) {
            gout[j] = (gap - thresh) / thresh;
        } else {
            gout[j] = (gpres - pfc) / (pfc + gpres + 1.0e-10);
        }
    }
    return 0;
}

// ---------------------------------------------------------------------------
// setGapClosed — NVI: unpack state then dispatch to typed onGapClosed hook.
// ---------------------------------------------------------------------------
void RodResiduals::setGapClosed(int layer, const double* y_data, const double* yp_data) {
    unpackState(y_data, yp_data);
    onGapClosed(layer);
}

} // namespace fred
