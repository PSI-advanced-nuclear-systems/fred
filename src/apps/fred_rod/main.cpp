#include "FredRodSolver.hpp"
#include "apps/fred_rod/fuelpelletmaterial/UO2.hpp"
#include "apps/fred_rod/claddingmaterial/AIM1.hpp"
#include "apps/fred_rod/gapmaterial/DummyGapMaterial.hpp"
#include <iostream>

// Minimal stand-alone driver for FRED-ROD.
// Use the Python interface to set up the problem. 

int main() {
    using namespace fred;

    // ---- Geometry --------------------------------------------------------
    FuelRodGeometry geom;
    geom.nf   = 5;       // fuel radial nodes
    geom.nc   = 3;       // cladding radial nodes
    geom.nz   = 1;       // axial layers (single-layer test)
    geom.rfi0 = 0.0;     // solid pellet (no central hole)
    geom.rfo0 = 4.1e-3;  // fuel outer radius [m]
    geom.rci0 = 4.2e-3;  // cladding inner radius [m] (gap = 0.1 mm)
    geom.rco0 = 4.75e-3; // cladding outer radius [m]
    geom.dz0  = {0.01};  // axial layer height [m]
    geom.vgp  = 1.0e-5;  // gas plenum volume [m3]
    geom.build();

    // ---- Materials -------------------------------------------------------
    UO2  fuel(10400.0);   // as-fabricated density 10400 kg/m3 (95% TD)
    AIM1 clad(7900.0);

    // ---- Solver setup ----------------------------------------------------
    DummyGapMaterial gap;
    FredRodSolver solver(geom, fuel, clad, gap);

    // Ramp power to 200 W/cm (linear), hold for 1000 s
    double qv_max = 2.0e8; // W/m3
    solver.setPowerDensityHistory({0.0, 100.0, 1000.0},
                                  {0.0, qv_max, qv_max});

    // Coolant temperature (sodium at 700 K)
    solver.setCoolantTemperature({0.0, 1000.0}, {700.0, 700.0});

    solver.setCoolantPressure(0.5);       // 0.5 MPa coolant pressure
    solver.setInitialTemperature(700.0);  // uniform initial temperature
    solver.setTolerances(1.0e-6, 1.0e-8);

    // ---- Run -------------------------------------------------------------
    solver.run(1000.0, 100.0);

    // ---- Report ----------------------------------------------------------
    auto peak = solver.peakFuelTemperature();
    std::cout << "\nTime [s]   Peak fuel T [K]\n";
    for (size_t k = 0; k < solver.timePoints().size(); ++k)
        std::cout << "  " << solver.timePoints()[k] << "   " << peak[k] << "\n";

    return 0;
}
