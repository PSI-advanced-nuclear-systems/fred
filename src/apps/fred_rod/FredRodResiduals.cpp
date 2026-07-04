#include "FredRodResiduals.hpp"
#include "../../platform/Constants.hpp"
#include "../../platform/GapModel.hpp"
#include "platform/FuelPelletMaterial.hpp"
#include "platform/CladdingMaterial.hpp"
#include <cassert>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace {

double computeRodGapConductance(
    const fred::AxialLayerState& s,
    const fred::FuelRodGeometry& geom,
    const fred::GapMaterial& gap_mat,
    double gpres)
{
    const int nf = geom.nf;
    const double T_fuel = s.T[nf - 1];
    const double T_clad = s.T[nf];

    if (gap_mat.isGasBond()) {
        // Gas bond: use actual gap geometry with Lanning-Hann gas-jump correction
        // (matches legacy FRED gaphtc) plus linearised radiation.
        const double k_gap = gap_mat.gapConductivity(0.5 * (T_fuel + T_clad));
        const double h_gas = fred::computeGapMediumConductance(
            T_fuel, T_clad, s.gap, gpres, k_gap, geom.ruff, geom.rufc);
        const double h_rad = fred::computeRadiationConductance(
            T_fuel, T_clad, geom.rfo0, geom.rci0);
        return h_gas + h_rad;
    } else {
        // Liquid bond: k / gap_eff, no gas-jump, no radiation.
        const double T_avg   = 0.5 * (T_fuel + T_clad);
        const double k_gap   = gap_mat.gapConductivity(T_avg);
        const double rough   = std::sqrt(geom.ruff * geom.ruff + geom.rufc * geom.rufc);
        const double gap_eff = std::max(s.gap, rough);
        return gap_mat.clampConductance(k_gap / gap_eff, !s.gapOpen);
    }
}

}

namespace fred {

FredRodResiduals::FredRodResiduals(const FuelRodGeometry&    geom,
                                   const FuelPelletMaterial& fuel,
                                   const CladdingMaterial&   clad,
                                   const GapMaterial&        gap_mat,
                                   double coolant_pressure_MPa)
    : RodResiduals(geom, fuel, clad, gap_mat, coolant_pressure_MPa, /*neq_irr=*/0),
      m_gap_mat(gap_mat),
      m_state  (geom.nf, geom.nc, geom.nz)
{}

void FredRodResiduals::initAlgebraicState(RodState& state,
                                           const FuelPelletMaterial& fuel,
                                           const CladdingMaterial& clad) const {
    const int nf = m_geom.nf;
    const int nc = m_geom.nc;

    // Seed the gas pressure algebraic variable.
    // IDACalcIC will refine it to the consistent initial condition.
    state.gpres = m_pressure_model ? m_pressure_model->initialPressure() : 0.1;

    double pcool  = m_pcoolFn ? m_pcoolFn(0.0) : m_pcool;

    for (int j = 0; j < m_geom.nz; ++j) {
        auto& s = state.layers[j];

        // ---- Thermal block ----
        s.qqv = (!m_layerPowerFns.empty() && j < (int)m_layerPowerFns.size())
                ? m_layerPowerFns[j](0.0)
                : (m_powerFn ? m_powerFn(0.0) : 0.0);

        double gap0 = m_geom.rci0 - m_geom.rfo0;
        s.gap = gap0;
        // Initialize surface radii so pack() stores the correct as-fabricated gap
        // even when stress-strain is OFF and updateGeometry() is never called.
        s.rfo = m_geom.rfo0;
        s.rci = m_geom.rci0;
        m_gapMgr.applyToLayer(j, s);
        s.pfc = 0.0;

        double T_cool = (!m_layerCoolantFns.empty() && j < (int)m_layerCoolantFns.size())
                        ? m_layerCoolantFns[j](0.0)
                        : (m_coolantTFn ? m_coolantTFn(0.0) : T_REF);
        s.T[nf + nc - 1] = T_cool;

        s.hgap = computeRodGapConductance(s, m_geom, m_gap_mat, state.gpres);

        // ---- Mechanical block: Lamé initial equilibrium ----
        // Pre-solve the elastic equilibrium so IDACalcIC starts from a
        // consistent state and does not need to take large Newton steps.
        if (m_mech_on && m_gapMgr.isGapOpen(j)) {
            // Fuel: solid cylinder exposed to fill-gas pressure → hydrostatic
            for (int i = 0; i < nf; ++i) {
                double Ef  = fuel.youngsModulus(s.T[i], fuel.referenceDensity());
                double nuF = fuel.poissonRatio();
                s.eft[i]   = fuel.thermalExpansionStrain(s.T[i]);
                s.sigfr[i] = -state.gpres;
                s.sigfh[i] = -state.gpres;
                s.sigfz[i] = -state.gpres;
                double e_mech = -state.gpres * (1.0 - 2.0*nuF) / Ef;
                s.efh[i] = s.eft[i] + e_mech;
                s.efr[i] = s.eft[i] + e_mech;
                s.efz    = s.eft[i] + e_mech;
            }

            // Cladding: Lamé thick-wall cylinder solution
            // sigr(r) = A + C/r^2,  sigh(r) = A − C/r^2
            // A = (gpres*rci^2 − pcool*rco^2) / (rco^2 − rci^2)
            // C = (pcool − gpres)*rci^2*rco^2 / (rco^2 − rci^2)
            // sigz = A  (axial force balance, open gap, uniform)
            double a = m_geom.rci0, b = m_geom.rco0;
            double a2 = a*a, b2 = b*b, denom = b2 - a2;
            double A_c = (state.gpres * a2 - pcool * b2) / denom;
            double C_c = (pcool - state.gpres) * a2 * b2 / denom;

            for (int i = 0; i < nc; ++i) {
                double ri  = m_geom.rad0[nf + i];
                double ri2 = ri * ri;
                double Ec  = clad.youngsModulus(s.T[nf + i]);
                double nuC = clad.poissonRatio();
                s.et[i]  = clad.thermalExpansionStrain(s.T[nf + i]);
                s.sigr[i] = A_c + C_c / ri2;
                s.sigh[i] = A_c - C_c / ri2;
                s.sigz[i] = A_c;
                s.eh[i] = s.et[i] + (s.sigh[i] - nuC*(s.sigr[i] + s.sigz[i])) / Ec;
                s.er[i] = s.et[i] + (s.sigr[i] - nuC*(s.sigh[i] + s.sigz[i])) / Ec;
                s.ez    = s.et[i] + (s.sigz[i] - nuC*(s.sigh[i] + s.sigr[i])) / Ec;
            }

            // Keep geometry-consistent fields in sync with the seeded strain state
            // before packing to the IDA vector.
            double dzf, dzc;
            s.updateGeometry(m_geom.rad0, m_geom.drf0, m_geom.drc0,
                             m_geom.dz0[j], dzf, dzc);
            s.gap = s.rci - s.rfo;
        }
    }
}

// ---------------------------------------------------------------------------
// solveMechanicalIC — Newton on the true mechanical unknowns (not the full
// y-vector). The y-vector has placeholder-zero slots for constraint rows whose
// Jacobian columns are identically zero, making Newton on the full vector
// impossible. Here we work on the 6nf+6nc+2 actual stress/strain unknowns
// directly, calling StressStrain::computeResiduals for the residual.
// ---------------------------------------------------------------------------

// Map mechanical struct fields → flat vector (matches StressStrain eq. order)
static void extractMech(const AxialLayerState& s, int nf, int nc, double* x) {
    int k = 0;
    for (int i = 0; i < nf; ++i) x[k++] = s.eft[i];
    for (int i = 0; i < nf; ++i) x[k++] = s.efh[i];
    for (int i = 0; i < nf; ++i) x[k++] = s.efr[i];
    x[k++] = s.efz;
    for (int i = 0; i < nf; ++i) x[k++] = s.sigfh[i];
    for (int i = 0; i < nf; ++i) x[k++] = s.sigfr[i];
    for (int i = 0; i < nf; ++i) x[k++] = s.sigfz[i];
    for (int i = 0; i < nc; ++i) x[k++] = s.et[i];
    for (int i = 0; i < nc; ++i) x[k++] = s.eh[i];
    for (int i = 0; i < nc; ++i) x[k++] = s.er[i];
    x[k++] = s.ez;
    for (int i = 0; i < nc; ++i) x[k++] = s.sigh[i];
    for (int i = 0; i < nc; ++i) x[k++] = s.sigr[i];
    for (int i = 0; i < nc; ++i) x[k++] = s.sigz[i];
}

static void injectMech(const double* x, int nf, int nc, AxialLayerState& s) {
    int k = 0;
    for (int i = 0; i < nf; ++i) s.eft[i] = x[k++];
    for (int i = 0; i < nf; ++i) s.efh[i] = x[k++];
    for (int i = 0; i < nf; ++i) s.efr[i] = x[k++];
    s.efz = x[k++];
    for (int i = 0; i < nf; ++i) s.sigfh[i] = x[k++];
    for (int i = 0; i < nf; ++i) s.sigfr[i] = x[k++];
    for (int i = 0; i < nf; ++i) s.sigfz[i] = x[k++];
    for (int i = 0; i < nc; ++i) s.et[i] = x[k++];
    for (int i = 0; i < nc; ++i) s.eh[i] = x[k++];
    for (int i = 0; i < nc; ++i) s.er[i] = x[k++];
    s.ez = x[k++];
    for (int i = 0; i < nc; ++i) s.sigh[i] = x[k++];
    for (int i = 0; i < nc; ++i) s.sigr[i] = x[k++];
    for (int i = 0; i < nc; ++i) s.sigz[i] = x[k++];
}

void FredRodResiduals::solveMechanicalIC(RodState& state, double t, double rtol, double atol) const {
    const int nf = m_geom.nf;
    const int nc = m_geom.nc;
    const int nm = m_mech.neq();  // = 6nf+6nc+2

    const double EPS_FD   = 1.0e-6;
    const int    MAX_ITER = 50;
    const double WRMS_TOL = 1.0e-6;

    double pcool = m_pcoolFn ? m_pcoolFn(t) : m_pcool;

    std::vector<double> x(nm), r(nm), rp(nm);
    std::vector<double> J((size_t)nm * nm), rhs(nm);

    for (int lay = 0; lay < m_geom.nz; ++lay) {
        auto& s = state.layers[lay];
        extractMech(s, nf, nc, x.data());

        for (int iter = 0; iter < MAX_ITER; ++iter) {
            injectMech(x.data(), nf, nc, s);
            double dzf, dzc;
            s.updateGeometry(m_geom.rad0, m_geom.drf0, m_geom.drc0,
                             m_geom.dz0[lay], dzf, dzc);

            m_mech.computeResiduals(s, state.gpres, pcool, r);

            double wrms = 0.0;
            for (int i = 0; i < nm; ++i) {
                double w = rtol * std::abs(x[i]) + atol;
                if (w < 1e-30) w = 1e-30;
                wrms += (r[i] / w) * (r[i] / w);
            }
            wrms = std::sqrt(wrms / nm);
            if (wrms < WRMS_TOL) {
                std::cout << "  Mech IC Newton (layer " << lay << "): "
                          << iter << " iter(s), WRMS=" << wrms << "\n";
                break;
            }
            if (iter == MAX_ITER - 1) {
                std::cerr << "  Mech IC Newton (layer " << lay
                          << "): not converged, WRMS=" << wrms << "\n";
                break;
            }

            // Numerical Jacobian: column-wise forward differences
            for (int jc = 0; jc < nm; ++jc) {
                double x0 = x[jc];
                double h = EPS_FD * (std::abs(x0) > 1e-15 ? std::abs(x0) : 1.0);
                x[jc] = x0 + h;
                injectMech(x.data(), nf, nc, s);
                s.updateGeometry(m_geom.rad0, m_geom.drf0, m_geom.drc0,
                                 m_geom.dz0[lay], dzf, dzc);
                m_mech.computeResiduals(s, state.gpres, pcool, rp);
                for (int i = 0; i < nm; ++i)
                    J[i*nm + jc] = (rp[i] - r[i]) / h;
                x[jc] = x0;
            }

            // Gaussian elimination with partial pivoting: J * dx = -r
            for (int i = 0; i < nm; ++i) rhs[i] = -r[i];
            for (int k = 0; k < nm; ++k) {
                int imax = k;
                for (int i = k+1; i < nm; ++i)
                    if (std::abs(J[i*nm+k]) > std::abs(J[imax*nm+k])) imax = i;
                if (imax != k) {
                    for (int jc = 0; jc < nm; ++jc) std::swap(J[k*nm+jc], J[imax*nm+jc]);
                    std::swap(rhs[k], rhs[imax]);
                }
                double piv = J[k*nm+k];
                if (std::abs(piv) < 1e-30) continue;
                for (int i = k+1; i < nm; ++i) {
                    double m = J[i*nm+k] / piv;
                    for (int jc = k; jc < nm; ++jc) J[i*nm+jc] -= m * J[k*nm+jc];
                    rhs[i] -= m * rhs[k];
                }
            }
            for (int k = nm-1; k >= 0; --k) {
                for (int jc = k+1; jc < nm; ++jc) rhs[k] -= J[k*nm+jc] * rhs[jc];
                double d = J[k*nm+k];
                rhs[k] = (std::abs(d) > 1e-30) ? rhs[k] / d : 0.0;
            }

            for (int i = 0; i < nm; ++i) x[i] += rhs[i];
        }

        // Write converged solution back to the layer state
        injectMech(x.data(), nf, nc, s);
        double dzf, dzc;
        s.updateGeometry(m_geom.rad0, m_geom.drf0, m_geom.drc0,
                         m_geom.dz0[lay], dzf, dzc);
    }
}

// ---------------------------------------------------------------------------
// pack / unpack  (RodState ↔ flat SUNDIALS vector)
// ---------------------------------------------------------------------------
void FredRodResiduals::pack(const RodState& state, double* y) const {
    const int nf = m_geom.nf;
    const int nc = m_geom.nc;

    y[0] = state.gpres;  // global gas pressure AE

    for (int j = 0; j < m_geom.nz; ++j) {
        double* yj = y + 1 + j * m_neq_j;
        const auto& s = state.layers[j];

        // Thermal block
        yj[0] = s.qqv;
        yj[1] = s.rci - s.rfo;
        yj[2] = s.pfc;
        yj[3] = s.hgap;
        for (int i = 0; i < nf+nc; ++i) yj[4 + i] = s.T[i];

        // Mechanics block — layout matches legacy FRED write_to_y order
        double* ym = yj + m_neq_th;
        int off = 0;
        for (int i = 0; i < nf; ++i) ym[off+i] = s.eft[i];  off += nf;   // eft[nf]
        for (int i = 0; i < nf; ++i) ym[off+i] = s.efh[i];  off += nf;   // efh[nf]
        for (int i = 0; i < nf; ++i) ym[off+i] = s.efr[i];  off += nf;   // efr[nf]
        ym[off] = s.efz;                                       off += 1;   // efz (scalar)
        for (int i = 0; i < nf; ++i) ym[off+i] = s.sigfh[i]; off += nf;  // sigfh[nf]
        for (int i = 0; i < nf; ++i) ym[off+i] = s.sigfr[i]; off += nf;  // sigfr[nf]
        for (int i = 0; i < nf; ++i) ym[off+i] = s.sigfz[i]; off += nf;  // sigfz[nf]
        for (int i = 0; i < nc; ++i) ym[off+i] = s.et[i];    off += nc;  // et[nc]
        for (int i = 0; i < nc; ++i) ym[off+i] = s.eh[i];    off += nc;  // eh[nc]
        for (int i = 0; i < nc; ++i) ym[off+i] = s.er[i];    off += nc;  // er[nc]
        ym[off] = s.ez;                                        off += 1;  // ez (scalar)
        for (int i = 0; i < nc; ++i) ym[off+i] = s.sigh[i];  off += nc;  // sigh[nc]
        for (int i = 0; i < nc; ++i) ym[off+i] = s.sigr[i];  off += nc;  // sigr[nc]
        for (int i = 0; i < nc; ++i) ym[off+i] = s.sigz[i];  off += nc;  // sigz[nc]
    }
}

void FredRodResiduals::unpack(const double* y, const double* yp, RodState& state) const {
    const int nf = m_geom.nf;
    const int nc = m_geom.nc;

    state.gpres = y[0];  // global gas pressure AE

    for (int j = 0; j < m_geom.nz; ++j) {
        const double* yj  = y  + 1 + j * m_neq_j;
        const double* ypj = yp + 1 + j * m_neq_j;
        auto& s = state.layers[j];

        s.qqv  = yj[0];
        s.gap  = yj[1];
        s.pfc  = yj[2];
        s.hgap = yj[3];
        for (int i = 0; i < nf+nc; ++i) {
            s.T[i]  = yj[4 + i];
            s.dT[i] = ypj[4 + i];
        }
        // Gap state is tracked by m_gapMgr (authoritative); copy into layer.
        m_gapMgr.applyToLayer(j, s);

        // Mechanics block — layout matches legacy FRED read_from_y order
        const double* ym = yj + m_neq_th;
        int off = 0;
        for (int i = 0; i < nf; ++i) s.eft[i]  = ym[off+i]; off += nf;
        for (int i = 0; i < nf; ++i) s.efh[i]  = ym[off+i]; off += nf;
        for (int i = 0; i < nf; ++i) s.efr[i]  = ym[off+i]; off += nf;
        s.efz = ym[off];                                      off += 1;
        s.dEfz = ypj ? ypj[m_neq_th + off - 1] : 0.0;
        for (int i = 0; i < nf; ++i) s.sigfh[i] = ym[off+i]; off += nf;
        for (int i = 0; i < nf; ++i) s.sigfr[i] = ym[off+i]; off += nf;
        for (int i = 0; i < nf; ++i) s.sigfz[i] = ym[off+i]; off += nf;
        for (int i = 0; i < nc; ++i) s.et[i]   = ym[off+i]; off += nc;
        for (int i = 0; i < nc; ++i) s.eh[i]   = ym[off+i]; off += nc;
        for (int i = 0; i < nc; ++i) s.er[i]   = ym[off+i]; off += nc;
        s.ez = ym[off];                                       off += 1;
        s.dEz = ypj ? ypj[m_neq_th + off - 1] : 0.0;
        for (int i = 0; i < nc; ++i) s.sigh[i]  = ym[off+i]; off += nc;
        for (int i = 0; i < nc; ++i) s.sigr[i]  = ym[off+i]; off += nc;
        for (int i = 0; i < nc; ++i) s.sigz[i]  = ym[off+i]; off += nc;

        // Update deformed geometry from unpacked strains
        double dzf_unused, dzc_unused;
        s.updateGeometry(m_geom.rad0, m_geom.drf0, m_geom.drc0,
                         m_geom.dz0[j], dzf_unused, dzc_unused);
    }
}

// ---------------------------------------------------------------------------
// RodResiduals virtual hook implementations
// ---------------------------------------------------------------------------
void FredRodResiduals::unpackState(const double* y, const double* yp) const {
    unpack(y, yp, m_state);
}

void FredRodResiduals::computeGlobalResiduals(double t, const double* /*yp*/, double* r) const {
    const double gpres = m_state.gpres;
    r[0] = m_pressure_model
         ? m_pressure_model->residual(t, gpres, m_state, m_geom)
         : gpres - 0.1;
}

void FredRodResiduals::computeThermalResiduals(int j, double t,
                                                const double* /*yj*/, const double* /*ypj*/,
                                                double* rj) const
{
    const auto& s        = m_state.layers[j];
    const double gpres   = m_state.gpres;
    const double qqv     = layerPower(j, t);
    const double T_cool  = layerCoolant(j, t);

    std::vector<double>& r_th = threadLocalThermalBuffer();
    std::fill(r_th.begin(), r_th.end(), 0.0);
    m_heat.computeResiduals(s, T_cool, gpres, qqv, r_th);
    r_th[3] = s.hgap - computeRodGapConductance(s, m_geom, m_gap_mat, gpres);
    std::copy(r_th.begin(), r_th.end(), rj);
}

void FredRodResiduals::freezeThermalResiduals(int j, const double* yj, double* rj) const {
    // Pin AEs to as-fabricated values; pin T ODEs to initial temperature.
    const double gap0 = m_geom.rci0 - m_geom.rfo0;
    rj[0] = yj[0];
    rj[1] = yj[1] - gap0;
    rj[2] = yj[2];
    rj[3] = yj[3];
    for (int i = 0; i < m_geom.nf + m_geom.nc; ++i)
        rj[4 + i] = yj[4 + i] - m_T_init;
    (void)j;
}

void FredRodResiduals::setGapOpen(int layer, bool open) {
    m_gapMgr.setGapOpen(layer, open);
    m_gapMgr.applyToLayer(layer, m_state.layers[layer]);
}

void FredRodResiduals::onGapClosed(int layer) {
    m_gapMgr.applyGapClosed(layer, m_state.layers[layer].efz, m_state.layers[layer].ez);
    m_gapMgr.applyToLayer(layer, m_state.layers[layer]);
}

} // namespace fred
