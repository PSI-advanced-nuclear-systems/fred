#include "HeatConduction.hpp"
#include "FuelPelletMaterial.hpp"
#include "CladdingMaterial.hpp"
#include "Constants.hpp"
#include <cmath>
#include <cassert>

namespace fred {

HeatConduction::HeatConduction(const FuelRodGeometry& geom,
                               const FuelPelletMaterial& fuel,
                               const CladdingMaterial& clad,
                               const GapMaterial& gap_mat)
    : m_geom(geom), m_fuel(fuel), m_clad(clad), m_gap_mat(gap_mat)
{}

void HeatConduction::computeResiduals(const AxialLayerState& s,
                                      double T_coolant,
                                      double gpres,
                                      double power_density,
                                      std::vector<double>& r,
                                      bool   isCladTempOuterBoundaryPrescribed,
                                      double h_cool) const
{
    (void)gpres;
    const int nf = m_geom.nf;
    const int nc = m_geom.nc;
    const int nr = nf + nc;
    assert((int)r.size() >= neq());

    // --- Equation 0: power density AE ---
    r[0] = s.qqv - power_density;

    // --- Equation 1: gap width AE ---
    //   gap = rci - rfo  (deformed geometry)
    r[1] = s.gap - (s.rci - s.rfo);

    // --- Equation 2: contact pressure AE ---
    //   open or liquid bond: pfc = 0 (no solid contact)
    //   gas bond closed:     pfc = |sigr[0]| (cladding inner surface stress magnitude)
    if (!s.gapOpen) {
        r[2] = s.pfc - std::abs(s.sigr[0]);
    } else {
        r[2] = s.pfc;
    }

    // --- Equation 3: gap conductance AE ---
    // The caller supplies the application-specific gap conductance residual.
    r[3] = 0.0;

    // --- Compute heat fluxes at all internal boundaries (W/m — per unit length) ---
    // We store flux * 2*pi*r at each half-node boundary.
    // Index convention:
    //   qf[i]   for i=0..nf-2  : fuel-fuel boundary between node i and i+1
    //   qf[nf-1]               : fuel-clad gap
    //   qf[nf..nf+nc-2]        : clad-clad boundary
    std::vector<double> qf(nr - 1, 0.0);

    // Fuel-fuel boundaries — conductivity dispatch (see AxialLayerState header).
    //
    // Path A (k_fuel_per_node empty): kf re-evaluated at the current T_half every
    // residual call, so the Newton Jacobian sees the live temperature dependence.
    // Used by FRED-ROD and FRED-OX, where k(T) is the primary variation and the
    // irradiation correction is a simple step-lagged scalar (k_irr_factor).
    //
    // Path B (k_fuel_per_node non-empty): kf is frozen for the entire Newton solve
    // at values pre-computed in afterAcceptedStep from the local post-redistribution
    // composition (zr_wf, pu_wf) and porosity at each node.  The one-step lag on T
    // is the trade-off for capturing the radial k(zr(r)) profile, which requires
    // per-node irradiation state that HeatConduction cannot access directly.
    // Used by FRED-M-Na once Zr redistribution has run at least one accepted step.
    for (int i = 0; i < nf-1; ++i) {
        double kf;
        if (!s.k_fuel_per_node.empty()) {
            // Path B: average adjacent pre-computed node values for the interface.
            kf = 0.5 * (s.k_fuel_per_node[i] + s.k_fuel_per_node[i+1]);
        } else {
            // Path A: evaluate live at current Newton iterate temperature.
            double T_half = 0.5 * (s.T[i] + s.T[i+1]);
            kf = m_fuel.thermalConductivity(T_half) * s.k_irr_factor;
        }
        double dr = m_geom.drf0;  // uniform fuel spacing
        qf[i] = kf * (s.T[i] - s.T[i+1]) / dr * 2.0 * PI * m_geom.rad_half[i];
    }

    // Fuel-clad gap (index nf-1): heat flux through gap [W/m]
    qf[nf-1] = s.hgap * (s.T[nf-1] - s.T[nf])
                * 2.0 * PI * m_geom.rfo0;  // use initial rfo as reference area

    // Clad-clad boundaries
    for (int i = nf; i < nr-1; ++i) {
        double T_half = 0.5 * (s.T[i] + s.T[i+1]);
        double kc = m_clad.thermalConductivity(T_half);
        double dr = m_geom.drc0;
        qf[i] = kc * (s.T[i] - s.T[i+1]) / dr * 2.0 * PI * m_geom.rad_half[i];
    }

    // --- Equations 4..4+nf-1: fuel temperature ODEs ---
    // Energy balance for ring i:
    //   rho * cp * A_i * dT/dt = qqv * A_i  +  q_in  -  q_out
    // where q_in/q_out are the boundary fluxes (W/m, i.e. per unit axial length).
    for (int i = 0; i < nf; ++i) {
        double cp  = m_fuel.heatCapacity(s.T[i]);
        double rho = m_fuel.referenceDensity();
        double Ai  = m_geom.area0[i]; // ring cross-sectional area [m2]

        double q_in  = (i == 0)    ? 0.0 : qf[i-1]; // zero heat at centre (symmetry)
        double q_out = qf[i];                          // nf-1 is the gap flux

        // r[n] = 0  ⟺  rho*cp*A*dT/dt = qqv*A + q_in - q_out
        r[4 + i] = s.qqv * Ai + q_in - q_out - rho * cp * Ai * s.dT[i];
    }

    // --- Equations 4+nf..4+nf+nc-2: inner cladding temperature ODEs ---
    // No volumetric heat generation in the cladding.
    for (int i = 0; i < nc-1; ++i) {
        int ii = nf + i; // global node index
        double cp  = m_clad.heatCapacity(s.T[ii]);
        double rho = m_clad.referenceDensity();
        double Ai  = m_geom.area0[ii];

        double q_in  = qf[ii - 1]; // from inner boundary
        double q_out = qf[ii];     // to outer boundary

        r[4 + nf + i] = q_in - q_out - rho * cp * Ai * s.dT[ii];
    }

    // --- Equation 4+nf+nc-1: outer cladding temperature AE ---
    if (isCladTempOuterBoundaryPrescribed) {
        // Dirichlet: outer surface temperature fixed by coolant BC
        r[4 + nf + nc - 1] = s.T[nr-1] - T_coolant;
    } else {
        // Robin: conductive flux in = convective flux out
        // qf[nf+nc-2] is the flux [W/m] from second-to-last clad ring to outer ring
        // convective loss [W/m] = h_cool * 2*pi*rco * (T_outer - T_coolant)
        const double rco_eff = m_geom.rco0;
        const double q_conv  = h_cool * 2.0 * M_PI * rco_eff * (s.T[nr-1] - T_coolant);
        r[4 + nf + nc - 1]  = qf[nf + nc - 2] - q_conv;
    }
}

} // namespace fred
