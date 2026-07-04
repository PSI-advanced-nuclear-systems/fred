SHELL = sh

# ===========================================================================
# SUNDIALS — defaults to the active conda environment
#
# Install via: conda install -c conda-forge sundials
#
# Override if your conda env is elsewhere:
#   make fred-m-na SUNDIALS_PREFIX=/path/to/conda/envs/fred-dev
#
# NOTE: the patched SUNDIALS from the legacy Fortran FRED installation is not
# compatible — use only a clean conda or source build.
# ===========================================================================

SUNDIALS_PREFIX     ?= $(CONDA_ENV)
SUNDIALS_INC         = $(SUNDIALS_PREFIX)/include
# Some installations put libs in lib64/ instead of lib/.
SUNDIALS_LIB        := $(or $(wildcard $(SUNDIALS_PREFIX)/lib64),$(SUNDIALS_PREFIX)/lib)

# ===========================================================================
# Python / pybind11 — default conda environment is fred-dev (overridable)
# Example override: make fred-m-na CONDA_ENV_NAME=calc
# ===========================================================================

CONDA_BASE      ?= /home/chan_y/miniforge3
CONDA_ENV_NAME  ?= fred-dev
CONDA_ENV       ?= $(CONDA_BASE)/envs/$(CONDA_ENV_NAME)
PYTHON          = $(CONDA_ENV)/bin/python3
CMAKE           = $(CONDA_ENV)/bin/cmake
PYBIND11_INC   := $(shell $(PYTHON) -c "import pybind11; print(pybind11.get_include())")
PYTHON_INC     := $(shell $(PYTHON)-config --includes)
PYTHON_LDFLAGS := $(shell $(PYTHON)-config --ldflags)
EXT_SUFFIX     := $(shell $(PYTHON) -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")
PYTHON_SITEPKG := $(shell $(PYTHON) -c "import site; print(site.getsitepackages()[0])")

# Source tree root
SRC = src

# Build directory
BUILD = build

# ===========================================================================
# Compiler flags
# ===========================================================================

CXX       = g++
CXXSTD    = -std=c++17
HDF5_INC  = $(CONDA_ENV)/include
HDF5_LIB  = $(CONDA_ENV)/lib

CXXFLAGS  = $(CXXSTD) -O2 -Wall -Wextra -MMD -MP -fopenmp -I$(SUNDIALS_INC) -I$(SRC) -I$(HDF5_INC)
CXXFLAGS_PY = $(CXXFLAGS) $(PYTHON_INC) -I$(PYBIND11_INC) -fPIC

SUNDIALS_STATIC = \
    $(SUNDIALS_LIB)/libsundials_ida.a \
    $(SUNDIALS_LIB)/libsundials_core.a \
    $(SUNDIALS_LIB)/libsundials_nvecserial.a \
    $(SUNDIALS_LIB)/libsundials_nvecopenmp.a \
    $(SUNDIALS_LIB)/libsundials_sunmatrixdense.a \
    $(SUNDIALS_LIB)/libsundials_sunlinsoldense.a \
    $(SUNDIALS_LIB)/libsundials_sunnonlinsolnewton.a

LDFLAGS   = -Wl,--start-group $(SUNDIALS_STATIC) -Wl,--end-group -lm -no-pie \
            -L$(HDF5_LIB) -lhdf5 -Wl,-rpath,$(HDF5_LIB) -fopenmp
LDFLAGS_PY= -Wl,--start-group $(SUNDIALS_STATIC) -Wl,--end-group -lm \
            $(PYTHON_LDFLAGS) -shared \
            -L$(HDF5_LIB) -lhdf5 -Wl,-rpath,$(HDF5_LIB) -fopenmp

# ===========================================================================
# Source file groups
# ===========================================================================

# --- Platform (shared by all apps) ----------------------------------------
PLAT_SRCS = \
    $(SRC)/platform/TimeTable.cpp \
    $(SRC)/platform/FredSolverBase.cpp \
    $(SRC)/platform/FredIdaSolverBase.cpp \
    $(SRC)/platform/GapModel.cpp \
    $(SRC)/platform/HeatConduction.cpp \
    $(SRC)/platform/StressStrain.cpp \
    $(SRC)/platform/GapPressureModel.cpp \
    $(SRC)/platform/GapStateManager.cpp \
    $(SRC)/platform/RodResiduals.cpp \
    $(SRC)/platform/SubchannelMode.cpp

PLAT_OBJS = $(patsubst $(SRC)/%.cpp,$(BUILD)/%.o,$(PLAT_SRCS))

# --- FRED-ROD -------------------------------------------------------------
ROD_LIB_SRCS = \
    $(SRC)/apps/fred_rod/fuelpelletmaterial/UO2.cpp \
    $(SRC)/apps/fred_rod/fuelpelletmaterial/MOX.cpp \
    $(SRC)/apps/fred_rod/claddingmaterial/AIM1.cpp \
    $(SRC)/apps/fred_rod/claddingmaterial/T91.cpp \
    $(SRC)/apps/fred_rod/fuelpelletmaterial/DummyFuelPellet.cpp \
    $(SRC)/apps/fred_rod/claddingmaterial/DummyCladding.cpp \
    $(SRC)/apps/fred_rod/FredRodResiduals.cpp \
    $(SRC)/apps/fred_rod/FredRodSolver.cpp

ROD_LIB_OBJS  = $(patsubst $(SRC)/%.cpp,$(BUILD)/%.o,$(ROD_LIB_SRCS))
ROD_BIND_OBJ  = $(BUILD)/apps/fred_rod/bindings.o
ROD_MOD_OBJ   = $(BUILD)/apps/fred_rod/module.o
ROD_MAIN_OBJ  = $(BUILD)/apps/fred_rod/main.o

# --- FRED-OX --------------------------------------------------------------
OX_LIB_SRCS = \
    $(SRC)/apps/fred_ox/fuelpelletmaterial/FredOxMOX.cpp \
    $(SRC)/apps/fred_ox/claddingmaterial/FredOxAIM1.cpp \
    $(SRC)/apps/fred_ox/claddingmaterial/FredOxT91.cpp \
    $(SRC)/apps/fred_ox/gapmaterial/FredOxGapMaterial.cpp \
    $(SRC)/apps/fred_ox/FredOxResiduals.cpp \
    $(SRC)/apps/fred_ox/FredOxSolver.cpp

OX_LIB_OBJS  = $(patsubst $(SRC)/%.cpp,$(BUILD)/%.o,$(OX_LIB_SRCS))
OX_BIND_OBJ  = $(BUILD)/apps/fred_ox/bindings.o
OX_MOD_OBJ   = $(BUILD)/apps/fred_ox/module.o

# --- FRED-M-Na ------------------------------------------------------------
MNA_LIB_SRCS = \
    $(SRC)/apps/fred_m_na/fuelpelletmaterial/UPuZr.cpp \
    $(SRC)/apps/fred_m_na/claddingmaterial/HT9.cpp \
    $(SRC)/apps/fred_m_na/FredMNaStressStrain.cpp \
    $(SRC)/apps/fred_m_na/FredMNaGapBehavior.cpp \
    $(SRC)/apps/fred_m_na/FredMNaResiduals.cpp \
    $(SRC)/apps/fred_m_na/FredMNaSolver.cpp \
    $(SRC)/apps/fred_m_na/FredMNaSubchannelMode.cpp

MNA_LIB_OBJS = $(patsubst $(SRC)/%.cpp,$(BUILD)/%.o,$(MNA_LIB_SRCS))
MNA_BIND_OBJ = $(BUILD)/apps/fred_m_na/bindings.o
MNA_MOD_OBJ  = $(BUILD)/apps/fred_m_na/module.o

# ===========================================================================
# Output targets
# ===========================================================================

EXE         = $(BUILD)/fred_rod.x

PYMOD_FRED  = $(BUILD)/fred$(EXT_SUFFIX)
PYMOD_ROD   = $(BUILD)/_fred_rod$(EXT_SUFFIX)
PYMOD_OX    = $(BUILD)/_fred_ox$(EXT_SUFFIX)
PYMOD_MNA   = $(BUILD)/_fred_m_na$(EXT_SUFFIX)

FRED_BASE_OBJ = $(BUILD)/bindings.o

# User-facing Python modules (one per app + shared infrastructure)
PYMOD_ROD_PY  = $(BUILD)/fred_rod.py
PYMOD_OX_PY   = $(BUILD)/fred_ox.py
PYMOD_MNA_PY  = $(BUILD)/fred_m_na.py
PYMOD_GEOM_PY = $(BUILD)/fred_geometry.py
PYMOD_MAT_PY  = $(BUILD)/fred_materials.py
PYMOD_BASE_PY = $(BUILD)/fred_solver_base.py
PYMOD_OUT_PY  = $(BUILD)/fred_output.py

PY_MODULES = $(PYMOD_ROD_PY) $(PYMOD_OX_PY) $(PYMOD_MNA_PY) \
             $(PYMOD_GEOM_PY) $(PYMOD_MAT_PY) $(PYMOD_BASE_PY) $(PYMOD_OUT_PY)

# ===========================================================================
# Top-level targets
# ===========================================================================

.PHONY: all fred fred-rod fred-ox fred-m-na exe install doc clean

all: fred fred-rod fred-ox fred-m-na

fred:      $(PYMOD_FRED) $(PY_MODULES)
fred-rod:  $(PYMOD_ROD)  $(PY_MODULES)
fred-ox:   $(PYMOD_OX)   $(PY_MODULES)
fred-m-na: $(PYMOD_MNA)  $(PY_MODULES)

exe: $(EXE)

# ===========================================================================
# Link rules — each app module is self-contained
#
# fred_ox and fred_m_na include ROD objects because bind_fred_rod() registers
# the shared platform abstract types (FuelRodGeometry, FuelPelletMaterial,
# CladdingMaterial, GapMaterial) that the app-specific concrete types inherit.
# ===========================================================================

$(EXE): $(PLAT_OBJS) $(ROD_LIB_OBJS) $(ROD_MAIN_OBJ)
	$(CXX) $(CXXSTD) -o $@ $^ $(LDFLAGS)
	@echo "Built: $@"

$(PYMOD_FRED): $(PLAT_OBJS) $(ROD_LIB_OBJS) $(ROD_BIND_OBJ) \
               $(OX_LIB_OBJS) $(OX_BIND_OBJ) \
               $(MNA_LIB_OBJS) $(MNA_BIND_OBJ) \
               $(FRED_BASE_OBJ)
	$(CXX) $(CXXSTD) -o $@ $^ $(LDFLAGS_PY)
	@echo "Built: $@"

$(PYMOD_ROD): $(PLAT_OBJS) $(ROD_LIB_OBJS) $(ROD_BIND_OBJ) $(ROD_MOD_OBJ)
	$(CXX) $(CXXSTD) -o $@ $^ $(LDFLAGS_PY)
	@echo "Built: $@"

$(PYMOD_OX): $(PLAT_OBJS) $(ROD_LIB_OBJS) $(ROD_BIND_OBJ) \
             $(OX_LIB_OBJS) $(OX_BIND_OBJ) $(OX_MOD_OBJ)
	$(CXX) $(CXXSTD) -o $@ $^ $(LDFLAGS_PY)
	@echo "Built: $@"

$(PYMOD_MNA): $(PLAT_OBJS) $(ROD_LIB_OBJS) $(ROD_BIND_OBJ) \
              $(MNA_LIB_OBJS) $(MNA_BIND_OBJ) $(MNA_MOD_OBJ)
	$(CXX) $(CXXSTD) -o $@ $^ $(LDFLAGS_PY)
	@echo "Built: $@"

# ===========================================================================
# Compile rules
# ===========================================================================

$(BUILD)/%.o: $(SRC)/%.cpp | $(BUILD)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS_PY) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/fred_rod.py:         python/fred_rod.py         | $(BUILD)
	cp $< $@
	@echo "Copied: $@"

$(BUILD)/fred_ox.py:          python/fred_ox.py          | $(BUILD)
	cp $< $@
	@echo "Copied: $@"

$(BUILD)/fred_m_na.py:        python/fred_m_na.py        | $(BUILD)
	cp $< $@
	@echo "Copied: $@"

$(BUILD)/fred_geometry.py:    python/fred_geometry.py    | $(BUILD)
	cp $< $@
	@echo "Copied: $@"

$(BUILD)/fred_materials.py:   python/fred_materials.py   | $(BUILD)
	cp $< $@
	@echo "Copied: $@"

$(BUILD)/fred_solver_base.py: python/fred_solver_base.py | $(BUILD)
	cp $< $@
	@echo "Copied: $@"

$(BUILD)/fred_output.py:      python/fred_output.py      | $(BUILD)
	cp $< $@
	@echo "Copied: $@"

# ===========================================================================
# Install / doc / clean
# ===========================================================================

install: all
	cp $(PYMOD_FRED) $(PYMOD_ROD) $(PYMOD_OX) $(PYMOD_MNA) $(PYTHON_SITEPKG)/
	cp $(PY_MODULES) $(PYTHON_SITEPKG)/
	@echo "Installed _fred_rod, _fred_ox, _fred_m_na + Python wrappers → $(PYTHON_SITEPKG)/"

doc:
	pdflatex -output-directory doc doc/fred_platform.tex
	pdflatex -output-directory doc doc/fred_platform.tex
	@echo "Documentation: doc/fred_platform.pdf"

clean:
	rm -rf $(BUILD)
	@echo "Cleaned build directory"

# ===========================================================================
# Dependency auto-generation
# ===========================================================================

ALL_OBJS = $(PLAT_OBJS) $(ROD_LIB_OBJS) $(ROD_BIND_OBJ) $(ROD_MOD_OBJ) $(ROD_MAIN_OBJ) \
           $(OX_LIB_OBJS)  $(OX_BIND_OBJ)  $(OX_MOD_OBJ) \
           $(MNA_LIB_OBJS) $(MNA_BIND_OBJ) $(MNA_MOD_OBJ)
-include $(ALL_OBJS:.o=.d)
