#pragma once
#include "FuelRodGeometry.hpp"
#include "FuelRodState.hpp"
#include "GapMaterial.hpp"
#include <vector>

namespace fred {

class FuelPelletMaterial;
class CladdingMaterial;

// Heat-conduction residuals for a single axial layer.
//
// Finite-volume energy balance in cylindrical geometry.
// The fuel and cladding are each discretised into rings; heat flows
// radially through conduction inside each zone and through the gap
// (gas + radiation + contact) at the fuel–clad interface.
//
// Equation layout for axial layer j (0-indexed), neq_thermal equations:
//
//   n=0                : qqv AE       — power density interpolation
//   n=1                : gap AE       — gap = rci - rfo
//   n=2                : pfc AE       — contact pressure (open: pfc=0; closed: pfc=-sigfr[0])
//   n=3                : hgap AE      — caller-supplied gap conductance residual
//   n=4..4+nf-1        : fuel T ODE   — energy balance, each fuel ring
//   n=4+nf..4+nf+nc-2  : clad T ODE  — energy balance, inner clad rings
//   n=4+nf+nc-1        : outer clad T AE  — prescribed by coolant BC
//
// neq_thermal = 4 + nf + nc
class HeatConduction {
public:
    HeatConduction(const FuelRodGeometry& geom,
                   const FuelPelletMaterial& fuel,
                   const CladdingMaterial& clad,
                   const GapMaterial& gap_mat);

    // Number of equations contributed by the thermal block per axial layer.
    // \coderef{heat_conduction_neq}
    int neq() const { return 4 + m_geom.nf + m_geom.nc; }

    // Fill residual vector r[0..neq()-1] for axial layer with state `s`.
    //
    // Additional inputs:
    //   T_coolant     : prescribed outer cladding surface temperature [K] (Dirichlet)
    //                   or coolant bulk temperature (Robin)
    //   gpres         : fill-gas pressure [MPa]
    //   power_density : time-interpolated volumetric power [W/m3]
    //   isCladTempOuterBoundaryPrescribed : true -> Dirichlet (default, all existing apps)
    //                                      false -> Robin BC using h_cool
    //   h_cool        : cladding-to-coolant HTC [W/(m2*K)] (only used when Robin)
    // \coderef{heat_conduction_residual}
    void computeResiduals(const AxialLayerState& s,
                          double T_coolant,
                          double gpres,
                          double power_density,
                          std::vector<double>& r,
                          bool   isCladTempOuterBoundaryPrescribed = true,
                          double h_cool = 0.0) const;

private:
    const FuelRodGeometry&    m_geom;
    const FuelPelletMaterial& m_fuel;
    const CladdingMaterial&   m_clad;
    const GapMaterial&        m_gap_mat;
};

} // namespace fred
