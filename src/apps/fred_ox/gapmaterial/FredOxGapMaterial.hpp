#pragma once
#include "platform/GapMaterial.hpp"

namespace fred {

// Gap conductance for FRED-OX.
//
// Implements the full legacy FRED gaphtc model for 'he' fill-gas with
// fission gas release.  Accounts for:
//   - He/Xe/Kr gas mixture conductivity (geometric mean rule)
//   - Lanning-Hann temperature-jump model
//   - Radiation conductance (KIT emissivities)
//
// State setters are called by FredOxSolver before each residual evaluation:
//   setFissionGasRelease() — from current fgrel [mol]
//   setGasInventory()      — fill-gas amount [mol]
//   setGasPressure()       — current gas pressure [MPa]
//   setBurnup()            — layer burnup [MWd/kgU]   (for future use)
//   setLinearPower()       — layer linear power [W/m]  (for fuel relocation)
//   setEnableRelocation()  — toggle FUEL_RELOC model
//
// The fission gas composition from legacy FRED gaphtc:
//   88.46% Xe,  7.69% Kr,  3.85% He  (mol fractions of released FG).
class FredOxGapMaterial : public GapMaterial {
public:
    // State setters — call before computeResiduals for each layer
    void setGasInventory(double mu0_mol)        { m_mu0    = mu0_mol; }
    void setFissionGasRelease(double fgrel_mol) { m_fgrel  = fgrel_mol; }
    void setGasPressure(double gpres_MPa)       { m_gpres  = gpres_MPa; }
    void setBurnup(double bu_MWd_kgU)           { m_burnup = bu_MWd_kgU; }
    void setLinearPower(double ql_Wm)           { m_ql     = ql_Wm; }
    void setEnableRelocation(bool en)           { m_reloc  = en; }

    // Accessors
    double gasInventory()   const { return m_mu0; }
    double fgRelease()      const { return m_fgrel; }
    double gasPressure()    const { return m_gpres; }

    // GapMaterial interface.
    double gapConductivity(double T) const override;
    double relocationFraction() const;

private:
    double m_mu0    = 0.0;   // He fill-gas inventory [mol]
    double m_fgrel  = 0.0;   // released fission gas [mol]
    double m_gpres  = 0.1;   // internal gas pressure [MPa]
    double m_burnup = 0.0;   // local burnup [MWd/kgU]
    double m_ql     = 0.0;   // local linear power [W/m]
    bool   m_reloc  = false; // enable FRAPCON fuel relocation model

    // Gas mixture thermal conductivity [W/(m·K)].
    // Geometric mean rule for He/Ar/Kr/Xe (FRED.f90 gaphtc).
    //   mu0_he  : moles of He fill gas
    //   mu0_ar  : moles of Ar fill gas (usually 0)
    //   fgrel   : moles of released fission gas
    //             composition: 88.46% Xe, 7.69% Kr, 3.85% He
    static double gasMixtureConductivity(double mu0_he, double mu0_ar,
                                         double fgrel, double T_avg);
};

} // namespace fred
