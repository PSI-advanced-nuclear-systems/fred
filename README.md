# FRED 2.0 Platform

FRED is a code package to (1) facilitate the development of 2D fuel performance codes 
and (2) provide ready-made applications for (a) generic fuel rods (b) MOX fuel and (c) 
metallic fuel in a sodium environment. 

This work grew as a result of a desire improving the user experience in running the code and 
specifically, to improve the developer experience in modifying and creating new applications 
built on shared libraries with modern object oriented and encapsulation methodologies. 

Currently, FRED comes with three ready-to-run applications: 
1. FRED-ROD: A generic fuel pin which solves thermal conduction and stress-strain fields 
2. FRED-OX: MOX fuel performance code based on the legacy FRED developed primarily by Konstantin Mikityuk
3. FRED-M-Na: Metallic fuel in sodium based on the work carried out in the ESFR-SIMPLE and ESFR-SMART frameworks

### Citations

Please cite the following publications if FRED is used:

\begin{verbatim}
    @software{FRED,
    author  = {PSI Advanced Nuclear Systems group},
    title   = {FRED},
    version = {2.0.0},
    year    = {2026},
    url     = {https://github.com/your-org/fred-ox}
    }
\end{verbatim}

For usage of FRED-OX: 

\begin{verbatim}
@article{Mikityuk2011FRED,
  author    = {Konstantin Mikityuk and Andrei Shestopalov},
  title     = {FRED fuel behaviour code: Main models and analysis of Halden IFA-503.2 tests},
  journal   = {Nuclear Engineering and Design},
  year      = {2011},
  volume    = {241},
  number    = {7},
  pages     = {2455--2461},
  doi       = {10.1016/j.nucengdes.2011.04.033}
}

@inproceedings{Kriventsev2015TopFuel,
  author    = {Vladimir Kriventsev and Andrei Rineiski and Werner Pfrang and
               Sara Perez-Martin and Konstantin Mikityuk and
               Grigori Khvostov},
  title     = {Benchmark on Behavior of MOX Fuel Pin Under Irradiation at Nominal Power in Sodium Fast Reactor},
  booktitle = {Proceedings of TopFuel 2015},
  year      = {2015},
  address   = {Zurich, Switzerland}
}
\end{verbatim}

For usage of FRED-M-Na:
\begin{verbatim}
@mastersthesis{Timpano2024FREDM,
  author      = {Timpano, Daniele},
  title       = {Metallic Fuel Behaviour Modeling for the European Sodium Fast Reactor},
  school      = {ETH Zurich},
  year        = {2024},
  month       = aug,
  type        = {Master's Thesis},
  doi         = {10.3929/ethz-b-000700694},
  url         = {https://doi.org/10.3929/ethz-b-000700694}
}
@techreport{Timpano2024ANL,
  author      = {Timpano, Daniele and Karahan, Aydin and Stauff, Nicolas E. and Mikityuk, Konstantin and Jim{\'e}nez, Antonio},
  title       = {Metallic Fuel Performance Analysis for the European Sodium Fast Reactor (ESFR-SIMPLE): Analysis of Metallic Fuel Performance Using SAS4A/SASSYS-1--MFUEL},
  institution = {Argonne National Laboratory},
  year        = {2024},
  number      = {ANL/NSE-24/24},
  doi         = {10.2172/2371706},
  url         = {https://doi.org/10.2172/2371706}
}

@inproceedings{Timpano2025TopFuel,
  author    = {Timpano, Daniele and Bubelis, Evaldas and Farina, Luca and Jim{\'e}nez-Carrascosa, Antonio and Karahan, Aydin and Mikityuk, Konstantin and Perez-Martin, Sara and Scolaro, Alessandro and Stauff, Nicolas E.},
  title     = {Fuel Performance Benchmark for the European Sodium Fast Reactor with Metallic Fuel},
  booktitle = {Proceedings of TopFuel 2025: Nuclear Reactor Fuel Performance Conference},
  year      = {2025},
  address   = {Nashville, Tennessee, USA},
  pages     = {1598--1607},
  doi       = {10.13182/TOPFUEL25-48153}
}

\end{verbatim}


## Prerequisites

### System tools

- `g++` (C++17) — OpenMP (`libgomp`) is bundled with GCC; no separate install needed
- `make`

On Ubuntu/Debian:

```bash
sudo apt install build-essential
```

### Python environment

The build system uses a conda environment called `fred-dev` by default. Create it
and install the required Python packages:

```bash
conda create -n fred-dev python=3.14
conda activate fred-dev
conda install -c conda-forge cmake sundials
pip install pybind11 numpy matplotlib h5py
```

If your conda base is not at `~/miniforge3`, override it at build time:

```bash
make fred-m-na CONDA_BASE=/path/to/conda
```

To use a different environment name:

```bash
make fred-m-na CONDA_ENV_NAME=myenv
```

## Building

### SUNDIALS

SUNDIALS is provided by the conda environment (installed in the step above via
`conda install -c conda-forge sundials`). Point the build at the conda prefix:

```bash
make fred-m-na SUNDIALS_PREFIX=~/miniforge3/envs/fred-dev
```

If your conda base is not at `~/miniforge3`:

```bash
make fred-m-na SUNDIALS_PREFIX=/path/to/conda/envs/fred-dev
```

> **Note:** the patched SUNDIALS from the legacy Fortran FRED installation is not
> compatible. Use only a clean (unpatched) build.

### Build a single application

```bash
make fred-rod      # builds build/fred_rod.<ext>.so
make fred-ox       # builds build/fred_ox.<ext>.so
make fred-m-na     # builds build/fred_m_na.<ext>.so
```

### Build all applications

```bash
make               # equivalent to: make all
make all
```

### Clean

```bash
make clean         # removes build/
```

### Install Python modules system-wide

```bash
make install       # copies all three .so files into the active Python site-packages
```

## Running examples

Each application has examples under `examples/<app>/`. The built modules are loaded
from the `build/` directory automatically by the example scripts.

Run using the `fred-dev` Python interpreter (not the system `python`):

```bash
# FRED-ROD examples
/home/chan_y/miniforge3/envs/fred-dev/bin/python \
    examples/fred_rod/heat_conduction_only/run.py

/home/chan_y/miniforge3/envs/fred-dev/bin/python \
    examples/fred_rod/heat_conduction_stress_strain/run.py

/home/chan_y/miniforge3/envs/fred-dev/bin/python \
    examples/fred_rod/stress_strain_only/run.py

/home/chan_y/miniforge3/envs/fred-dev/bin/python \
    examples/fred_rod/gapclose_reopen/run.py

# FRED-OX examples
/home/chan_y/miniforge3/envs/fred-dev/bin/python \
    examples/fred_ox/base_irradiation/irradiate_10d/run.py

# FRED-M-Na examples
/home/chan_y/miniforge3/envs/fred-dev/bin/python \
    examples/fred_m_na/simple_irradiation/run.py
```

Alternatively, activate the environment first:

```bash
conda activate fred-dev
python examples/fred_m_na/simple_irradiation/run.py
```

## Source layout

```
src/
  platform/          # shared solver core (IDA DAE, gap, heat, stress-strain)
  apps/
    fred_rod/        # FRED-ROD app + material library + Python bindings
    fred_ox/         # FRED-OX app + material library + Python bindings
    fred_m_na/       # FRED-M-Na app + material library + Python bindings
examples/
  fred_rod/
  fred_ox/
  fred_m_na/
python/              # fred_input.py, fred_output.py helpers (copied to build/)
doc/                 # LaTeX documentation
build/               # compiled outputs (gitignored)
```

## SUNDIALS dependency

FRED 2.0 uses the [SUNDIALS IDA](https://github.com/LLNL/sundials) DAE solver.
The conda-forge `sundials` package is used (installed as part of the `fred-dev`
environment). The IDA solver, OpenMP NVector, dense linear solver, and Newton
nonlinear solver are the components required at link time.
