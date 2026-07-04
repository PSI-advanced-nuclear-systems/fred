# FRED 2.0

## overview
FRED 2.0 is a platform to develop fuel performance codes. 

The legacy FRED code was not object oriented and incorporated empirical correlations specific to MOX. The FRED platform takes on a more modular and object oriented approach, where we decouple the actual equations being solved: thermal conductivit/ stress-strain, with base irradiation type of correlations. 

FRED-ROD: solves thermal conductivity and stress strain equations, including gap closing/reopening logic. New material correlations can be added to FRED-ROD for integration into testing/ verification pipelines. 

FRED-OX: Now taking the logic from FRED-ROD and adding base irradiation effects. (1) specific empirical correlations for base irradiation on some materials (2) gas mixture conductivity from the release of empirical fission gas model. 

## history 

FRED is an old code developed by Konstantin Mikityuk. 

In 2018, a version of this (alias: metal branch) was forked by Daniele Timpano to incorporate metallic fuel properties. Path to source: /mnt/c/Users/chan_y/Documents/timpano_fred/01.new/07.FRED_M_OCT24.SRC 

In 2026, Konstantin and Yi Meng Chan incorporated the SUNDIALS solver and fixed gap-closing/ reopening logic (alias: SUNDIALS branch), path to source: /home/chan_y/Documents/fred

In 2026, Yi Meng Chan is rewriting the code to improve modularity and object oriented approach, with the aim of combining `metal branch` and `Sundials branch`. The structure is now that core thermal conductance and stress-strain equations are included in a `FRED platform` layer. Applications are developed on top of this platform. 

`FRED-ROD`: solves thermal conduction and stress-strain. User encouraged to implement their own material properties for e.g. for specifying base irradiated material properties etc. No irradiation behaviour is modelled in the code. This approach of decoupling the actual physics equations from base irradiation effects is deliberate to allow for experimental workflows. 

`FRED-OX` : uses the same equations as `FRED-ROD` and implements additional physics applied to MOX fuel: (1) base irradiation effects on fuel, cladding properties from empircal correlations (2) fission gas production which affects gas pressure, gap conductance properties etc. 

`FRED-M-Na` : uses the same equations as `FRED-ROD` but implements additional physics related to metallic fuel (U-Pu-Zr) in sodium, incorporating many of the correlations in SAS4A. 