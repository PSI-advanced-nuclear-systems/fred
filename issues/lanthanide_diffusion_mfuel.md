# Lanthanide Diffusion in Metallic SFR Fuel: Theory and Numerical Implementation for FRED

Currently, FRED-M-Na estimates the cladding wastage thickness based on the precipitation kinetics (PK) model, where the coefficients of the PK model were optimized based on simulations of SAS4A/MFUEL on ESFR-M pin types. This work was conducted by D. Timpano and documented partially in his master's thesis, inspired by the PK model implemented by Karahan's Phd thesis (2007), see the theory manual for more elaboration. However, it is not certain if these coefficients can be generalized to any fuel pin geometry, hence there is a motivation to make a copy of the more physically correct MFRED Lanthanide diffusion model which (1) updates the diffusion coefficient by local fuel node temperature and calculates the (2) lanthanide concentrations and most importantly the Lanthanide flux at the fuel-cladding interface, which determines the growth of the cladding wastage layer. 


**Source of model parameters:** SAS4A/SASSYS-1 v5.7 Theory Manual, Appendix 9.3 ([link](https://sas-doc.nse.anl.gov/5.7/Part03/Ch09/clad_param.html)).

---

## 1. Governing Equation

The lanthanide atom number density $u(r,t)$ [#/m³] is governed by the 1D radial diffusion equation with a uniform fission source in cylindrical geometry:

$$\frac{\partial u}{\partial t} = \frac{1}{r}\frac{\partial}{\partial r} \left(r\,D_{\mathrm{La}}(T)\,\frac{\partial u}{\partial r}\right) + S$$

where:

| Symbol | Definition | Value / units |
|---|---|---|
| $u(r,t)$ | Lanthanide atom number density | #/m³ |
| $r$ | Radial coordinate from pin centreline | m |
| $t$ | Time | s |
| $D_{\mathrm{La}}(T)$ | Effective lanthanide diffusion coefficient | m²/s - see §1.2 |
| $S = Y_{\mathrm{La}}\,q'''$ | Volumetric fission source of lanthanides | #/(m³·s) |
| $Y_{\mathrm{La}}$ | Lanthanide fission yield per unit energy | $15.29 \times 10^{9}$ #/J |
| $q'''$ | Local volumetric power density | W/m³ |

The domain spans $r \in [0, R_{\mathrm{fc}}]$ where $R_{\mathrm{fc}}$ is the fuel outer radius (fuel–cladding interface).

### 1.1 Boundary conditions

**Centreline** - symmetry (no flux through the axis):

$$\left.\frac{\partial u}{\partial r}\right|_{r=0} = 0$$

**Fuel–cladding interface** - Dirichlet sink (lanthanides precipitate instantly upon arrival):

$$u(R_{\mathrm{fc}},t) = 0$$

**Initial condition:**

$$u(r,0) = 0$$

The zero initial condition is valid assuming no lanthanide pre-exists in a fresh fuel pin. For a restart from a partially irradiated state, $u(r,0)$ is taken from the stored radial profile at the previous shutdown.

### 1.2 Diffusion coefficient - Arrhenius form with irradiation floor

(MFUEL Appendix 9.3, Table 9.8.5)

$$D_{\mathrm{La}}(T) = D_0 \exp\!\left(-\frac{Q_0}{R_g T^{\ast}}\right), \quad T^{\ast} = \max(T, T_{\mathrm{floor}})$$

| Parameter | Symbol | Value |
|---|---|---|
| Pre-exponential factor | $D_0$ | $5.0$ m²/s |
| Activation energy | $Q_0$ | $250.0 \times 10^3$ J/mol |
| Universal gas constant | $R_g$ | $8.314$ J/(mol·K) |
| Irradiation floor temperature | $T_{\mathrm{floor}}$ | $800$ K |

The floor at 800 K prevents $D_{\mathrm{La}}$ from collapsing to zero at low temperatures. Physically, it accounts for irradiation-enhanced diffusivity: the fast neutron flux maintains a supersaturation of vacancies that keeps lanthanides mobile even when purely thermal diffusion would be negligible.

**Numerical values at key temperatures:**

| $T$ (K) | $T^*$ (K) | $D_{\mathrm{La}}$ (m²/s) |
|---|---|---|
| $\leq 800$ | 800 | $5.0\exp(-250000/(8.314\times800)) \approx 2.3\times10^{-16}$ |
| 850 | 850 | $\approx 5.6\times10^{-15}$ |
| 900 | 900 | $\approx 7.0\times10^{-14}$ |
| 950 | 950 | $\approx 5.1\times10^{-13}$ |
| 1000 | 1000 | $\approx 2.3\times10^{-12}$ |
| 1100 | 1100 | $\approx 1.7\times10^{-10}$ |

The four-order-of-magnitude variation between 800 K and 1100 K means that the radial temperature profile drives the diffusion strongly: a ring of high-temperature fuel near the centre is far more effective at driving La outward than simple concentration-gradient arguments would suggest. This is why a non-uniform $D(r,T)$ cannot be handled analytically and requires the numerical solution (§3).

### 1.3 Cladding wastage thickness

Since no diffusion of lanthanides is modelled within the cladding (it acts as a pure sink), the wastage depth is obtained by integrating the arriving lanthanide flux at the fuel–cladding interface:

$$\boxed{\delta(t) = \frac{1}{u_{\mathrm{sol}}}\int_0^t J_{\mathrm{La}}(R_{\mathrm{fc}},\,\tau)\,d\tau}$$

where $J_{\mathrm{La}} = -D_{\mathrm{La}}(T)\ \partial u/\partial r\big|_{r=R_{\mathrm{fc}}}$ is the outward lanthanide flux [#/(m²·s)] and $u_{\mathrm{sol}}$ is the solubility limit of the cladding steel (MFUEL Appendix 9.3, Table 9.8.5):

| Cladding grade | $u_{\mathrm{sol}}$ (#/m³) |
|---|---|
| HT9 (ferritic-martensitic) | $6.4\times10^{27}$ |
| D9 (austenitic) | $0.5\times10^{27}$ |

Physically, $u_{\mathrm{sol}} \cdot \delta$ is the total number of lanthanide atoms per unit interfacial area that the wastage layer has absorbed. The layer grows until cladding failure, declared when:

$$\delta(t) \geq 0.5\,t_{\mathrm{clad,0}}$$

where $t_{\mathrm{clad,0}}$ is the initial cladding wall thickness (MFUEL §9.2.11). 

The load-bearing cladding thickness passed to the MFUEL mechanical model at each time step is $t_{\mathrm{clad,0}} - \delta(t)$. Currently, there are no plans to have apply feedback on the cladding wastage layer thickness to the stress-strain model in FRED, so the wastage thickness is solely used to assess operating limits of the fuel pin under base irradiation. 