#include "GapPressureModel.hpp"
#include <algorithm>

namespace fred {

GapPressureModel::GapPressureModel(double gpres0_MPa, const FuelRodGeometry& geom,
                                    const GapMaterial& gapmat)
    : m_gpres0(gpres0_MPa), m_Tplenum(T_REF), m_geom(&geom), m_gapmat(&gapmat)
{
    m_h_total = 0.0;
    for (int j = 0; j < geom.nz; ++j)
        m_h_total += geom.dz0[j];

    double vfree0;
    if (gapmat.isGasBond()) {
        // Gas bond: fill gas occupies annular gap + central void + plenum.
        vfree0 = PI * (geom.rci0*geom.rci0 - geom.rfo0*geom.rfo0
                       + geom.rfi0*geom.rfi0) * m_h_total + geom.vgp;
    } else {
        // Liquid bond: gas is confined to the plenum.
        vfree0 = geom.vgp;
    }
    m_mu0 = gpres0_MPa * 1.0e6 * vfree0 / (R_GAS * T_REF);
}

void GapPressureModel::setInitialGasPressure(double gpres0_MPa) {
    m_gpres0 = gpres0_MPa;
    double vfree0;
    if (m_gapmat->isGasBond()) {
        vfree0 = PI * (m_geom->rci0*m_geom->rci0 - m_geom->rfo0*m_geom->rfo0
                       + m_geom->rfi0*m_geom->rfi0) * m_h_total + m_geom->vgp;
    } else {
        vfree0 = m_geom->vgp;
    }
    m_mu0 = gpres0_MPa * 1.0e6 * vfree0 / (R_GAS * T_REF);
}

double GapPressureModel::residual(double t, double gpres,
                                   const RodState& state,
                                   const FuelRodGeometry& geom) const
{
    if (m_plenumTFn) m_Tplenum = m_plenumTFn(t);

    if (m_mu0 < 1e-30)
        return gpres - m_gpres0;

    double vt = 0.0;
    if (m_gapmat->isGasBond()) {
        vt = computeGasBondVT(geom, state.layers, m_Tplenum);
    } else {
        // Liquid bond: v_T = (V_plenum + delta annular gap volume) / T_plenum
        const double gap_area0 = PI * (geom.rci0*geom.rci0 - geom.rfo0*geom.rfo0);
        vt = geom.vgp;
        for (int j = 0; j < geom.nz; ++j) {
            const auto& s = state.layers[j];
            vt += PI * geom.dz0[j] * (s.rci*s.rci - s.rfo*s.rfo);
        }
        vt -= gap_area0 * m_h_total;
        vt  = vt / std::max(m_Tplenum, 1.0);
    }

    if (vt < 1e-30)
        return gpres - m_gpres0;

    return gpres - m_mu0 * R_GAS / vt * 1.0e-6;
}

double GapPressureModel::initialPressure() const {
    return m_gpres0;
}

void GapPressureModel::setInitialTemperature(double T0) {
    m_Tplenum = T0;
}

void GapPressureModel::setPlenumTemperature(std::function<double(double)> fn) {
    m_plenumTFn = std::move(fn);
}

} // namespace fred
