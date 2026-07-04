#pragma once
#include "platform/FuelRodState.hpp"
#include <vector>

namespace fred {

// Extended axial-layer state for FRED-OX.
//
// Adds irradiation state fields on top of the base AxialLayerState:
// \coderef{fred_ox_layer_state}
//   bup  : local burnup [MWd/kgU]
//   defs : fuel volumetric swelling strain rate per node [1/s], size nf
//   dec  : cladding hoop creep strain rate per ring [1/s], size nc
//
// Note: efs (fuel swelling) and ec (cladding creep) are defined in the base
// AxialLayerState so that StressStrain can apply them generically.
//
// The bup_y field is the y-vector representation:
//   bup_y = bup * 8.64e10   [J/kgU]
struct FredOxLayerState : public AxialLayerState {
    double bup  = 0.0; // local burnup [MWd/kgU]
    double dbup = 0.0; // burnup rate [MWd/kgU/s]

    std::vector<double> defs; // fuel swelling rate [1/s], size nf
    std::vector<double> dec;  // cladding hoop creep rate [1/s], size nc

    explicit FredOxLayerState(int nf_, int nc_)
        : AxialLayerState(nf_, nc_),
          defs(nf_, 0.0),
          dec (nc_, 0.0)
    {}
};

// Full rod state for FRED-OX.
// \coderef{fred_ox_rod_state}
struct FredOxRodState {
    int nf, nc, nz;

    // Global gas state
    double gpres  = 0.1; // internal gas pressure [MPa]
    double mu0    = 0.0; // initial fill-gas inventory [mol]
    double fggen  = 0.0; // total fission gas generated [mol]
    double fgrel  = 0.0; // total fission gas released [mol]
    double dfggen = 0.0; // generation rate [mol/s]
    double dfgrel = 0.0; // release rate [mol/s]

    std::vector<FredOxLayerState> layers;

    FredOxRodState(int nf_, int nc_, int nz_)
        : nf(nf_), nc(nc_), nz(nz_),
          layers(nz_, FredOxLayerState(nf_, nc_))
    {}
};

} // namespace fred
