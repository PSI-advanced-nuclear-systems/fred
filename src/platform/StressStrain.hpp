#pragma once
#include "FuelRodGeometry.hpp"
#include "FuelRodState.hpp"
#include <vector>

namespace fred {

class FuelPelletMaterial;
class CladdingMaterial;

// Thermo-elastic + irradiation stress-strain residuals for a single axial layer.
//
// Fuel: includes isotropic swelling via AxialLayerState::efs (0 for FRED-ROD).
// Cladding: includes hoop creep via AxialLayerState::ec (0 for FRED-ROD).
//   Volume-preserving decomposition: ec_hoop = +ec, ec_radial = ec_axial = -ec/2.
//   The ODE d(ec)/dt = creepRate(T, sigh) is integrated by FRED-OX irradiation residuals.
//
// Equation layout (continuing from the thermal block, base index 0):
//
// FUEL (nf nodes):
//   [0..nf-1]           eft AE     : fuel linear thermal expansion AE
//                                    r = eft[i] - (ftexp(T[i]) - ftexp(T_REF))
//   [nf..2nf-1]         efh AE     : Hooke's law, hoop component
//   [2nf..3nf-1]        efr AE     : Hooke's law, radial component
//   [3nf..4nf-1]        efz AE     : Hooke's law, axial component
//   [4nf..5nf-2]        compat     : strain compatibility (nf-1 equations)
//   [5nf-1..6nf-3]      equil      : stress equilibrium   (nf-1 equations)
//   [6nf-2]             BC1_fuel   : axial force balance
//   [6nf-1]             BC2_fuel   : inner radial BC (symmetry or central-hole)
//   [6nf]               BC3_fuel   : outer radial BC (gap open → sigfr=-gpres;
//                                                     gap closed → defz_fuel=defz_clad)
// CLADDING (nc nodes):
//   [6nf+1..6nf+nc]     et AE      : clad thermal expansion AE
//   [6nf+nc+1..6nf+2nc] eh AE      : Hooke's law, hoop
//   [6nf+2nc+1..6nf+3nc]er AE      : Hooke's law, radial
//   [6nf+3nc+1..6nf+4nc]ez AE      : Hooke's law, axial
//   [6nf+4nc+1..6nf+5nc-1] compat  : strain compatibility (nc-1)
//   [6nf+5nc-1..6nf+6nc-3] equil   : stress equilibrium   (nc-1)
//   [6nf+6nc-2]         BC1_clad   : axial force or interface stress
//   [6nf+6nc-1]         BC2_clad   : inner radial stress
//   [6nf+6nc]           BC3_clad   : outer radial stress = -pcool
//
// neq_mech = 7*(nf + nc)  — verify: 3 BCs each side + (nf-1+nf-1) + 4nf + (nc-1+nc-1) + 4nc
//           = 3 + 2(nf-1) + 4nf + 3 + 2(nc-1) + 4nc
//           = 6 + 2nf - 2 + 4nf + 2nc - 2 + 4nc - 2 + 4
//           = 6 + 6nf + 6nc
// Corrected: 4nf + (nf-1) + (nf-1) + 3 + 4nc + (nc-1) + (nc-1) + 3
//          = 4nf + 2nf - 2 + 3 + 4nc + 2nc - 2 + 3 = 6nf + 6nc + 2 -- so neq_mech = 6nf+6nc+2
// (combined neq per layer = (4+nf+nc) + (6nf+6nc+2) = 6 + 7nf + 7nc -- matches FuelRodState)
class StressStrain {
public:
    StressStrain(const FuelRodGeometry& geom,
                 const FuelPelletMaterial& fuel,
                 const CladdingMaterial& clad);

    // Number of equations contributed by the mechanics block per axial layer.
    // \coderef{stress_strain_neq}
    int neq() const;

    // Fill residual vector r[0..neq()-1] for axial layer with state `s`.
    //
    //  gpres  : fill-gas pressure [MPa]
    //  pcool  : coolant pressure [MPa]  (= outer clad BC: sigr[nc-1] = -pcool)
    // \coderef{stress_strain_residual}
    void computeResiduals(const AxialLayerState& s,
                          double gpres, double pcool,
                          std::vector<double>& r) const;

private:
    const FuelRodGeometry&    m_geom;
    const FuelPelletMaterial& m_fuel;
    const CladdingMaterial&   m_clad;
};

} // namespace fred
