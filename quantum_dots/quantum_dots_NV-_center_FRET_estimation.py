
# coding: utf-8

# # Estimation of FRET efficiency for emission from WS2 and absorption by diamond NV- center.
# 
# Based on [Medintz & Hildebrandt, 2014](https://books.google.com/books/about/FRET_F%C3%B6rster_Resonance_Energy_Transfer.html?id=GXWAAQAAQBAJ), pages 23-31, and [Wu & Brand et al., 1994](https://doi.org/10.1006/abio.1994.1134), page 3.
# 
# FRET efficiency, $E$, is
# $$
# E = \frac{R_0^6}{R_0^6 + r_{DA}^6}
# $$
# where
# $r_{DA}$ is the distance between donor and acceptor
# and 
# $R_0$ is the Förster distance.
# This is also expressed as
# $$
# E = \frac{k_T}{k_T + 1/\tau_D}
# $$
# where
# $k_T$ is the rate of energy transfer
# and $t_D$ is the lifetime of the donor excited state in the absence of acceptor.
# The formula for $k_T$ is
# $$
# k_T = \frac{1}{\tau_D} \left( \frac{R_0^6}{r_{DA}^6} \right)
# $$

# The Förster distance $R_0^6$ is a key quantity to compute,
# and is given as
# $$
# R_0^6 = \frac{9 (\ln 10) \kappa^2 \Phi_D J }{128\pi^5 n^4 N_A}
# $$
# where $\kappa^2$ is the orientation factor,
# $\Phi_D$ is the quantum yield of the donor fluorescence in the absence of acceptor,
# $J$ is the overlap integral,
# $\pi = 3.14159...$,
# $n$ is the index of refraction of the medium,
# and $N_A$ is the Avogadro number $6.0221415 \times 10^{23}$ per mol.

# The overlap integral $J$, is computed differently depending on whether the spectra are given in terms of wavelength, wave number, or frequency. For wavelength,
# $$
# J^\lambda = J = \int f_D(\lambda) \epsilon_A(\lambda) \lambda^4 \mathrm{d} \lambda
# $$
# where $f_D(\lambda)$ is the fluorescence spectrum of the donor,
# $\epsilon_A(\lambda)$ is the molar extinction coefficient (a.k.a. molar absorptivity) of the acceptor (usually in units of $\mathrm{M}^{-1} \mathrm{cm}^{-1}$), and $\lambda$ is the wavelength (usually in nm).

# As a discrete sum, $J$ is computed as
# $$
# J = \frac{\sum_i F_D(\lambda_i) \epsilon_A(\lambda_i) \lambda_i^4}{\sum_i F_D(\lambda_i)}
# $$
# where $F_D(\lambda_i)$ is the intensity of the donor spectrum at $\lambda_i$ (it need not be normalized),
# $\epsilon_A(\lambda_i)$ is the molar attenuation coefficient at $\lambda_i$, and $\lambda_i$ is the wavelength.

# In[83]:


import math
import numpy as np
import matplotlib.pyplot as plt
import IPython.display


# ## WS2 emission spectrum

# In[84]:


WS2_eV, WS2_intensity, _, _, _, _ = np.loadtxt("WS2_emission.txt", unpack=True)
plt.clf()
plt.scatter(WS2_eV, WS2_intensity,color="blue")
plt.xlabel("Photon energy (eV)")
plt.ylabel("Intensity (Counts)");


# In[85]:


hc = 1239.8 # eV nm
WS2_nm = hc/WS2_eV
plt.clf()
plt.scatter(WS2_nm, WS2_intensity,color="tab:blue")
plt.title("Emission spectrum of $\mathrm{WS_2}$")
plt.xlabel("Photon wavelength (nm)")
plt.ylabel("Intensity (Counts)");


# ## NV- center absorption spectrum

# In[86]:


nv_center_emission_nm, nv_center_emission_intensity, _, _, _, _ = np.loadtxt("NV_center_emission.txt", unpack=True)
plt.clf()
plt.scatter(nv_center_emission_nm, nv_center_emission_intensity,color="tab:red")
plt.title("Emission spectrum of NV- center")
plt.xlabel("Photon wavelength (nm)")
plt.ylabel("Emission intensity (a.u.)");


# Compare to Figure 7 of Wee et al., 2007:

# In[87]:


IPython.display.Image(filename="nv-center-excitation-spectrum/fig7.png")


# > Figure 7. Comparison of one-photon and two-photon excited fluorescence spectra of nitrogen-vacancy centers in type Ib diamond treated with 3 MeV proton irradiation.

# > Figure 7 displays the emission spectra of the (N-V)⁻ centers excited by 532 and 1064 nm laser pulses at room temperature for a sample prepared with the 3 MeV proton irradiation. Both spectra were collected at an excitation time of 4 s using a 40× microscope objective with an incident laser power of 0.08 μW and 4.6 mW for the one-photon and two-photon excitations, respectively. In addition to the ZPL of the (N-V)⁻ center at 637 nm, a sharp ZPL derived from the (N-V)⁰ center can be found at 575 nm (2.156 eV) in both spectra. Interestingly, this center cannot be identified in the absorption spectrum even at 80 K (Figure 6b) but reveals itself clearly in the emission profile.

# Since this is the emission spectrum of the diamond NV- center, not the absorption / excitation spectrum, we will need to flip the spectrum around the zero phonon line at 637 nm to approximate the absorption spectrum.

# In[88]:


nv_center_absorption_nm = 637 - (nv_center_emission_nm - 637)


# In[89]:


plt.clf()
plt.scatter(nv_center_emission_nm, nv_center_emission_intensity, label="emission", color="tab:red")
plt.scatter(nv_center_absorption_nm, nv_center_emission_intensity, label="mirror of emission around 637 nm", color="tab:orange")
plt.vlines(x=637, ymin=0, ymax=max(nv_center_emission_intensity), color="red", label="zero phonon line (637 nm)")
plt.vlines(x=532, ymin=0, ymax=max(nv_center_emission_intensity), color="green", label="excitation wavelength (532 nm)")
plt.legend(loc="center right", bbox_to_anchor=(1.65,0.8))
plt.xlabel("Photon wavelength (nm)")
plt.ylabel("Intensity (a.u.)")
plt.savefig("NV_emission_flipped_around_ZPL.eps", bbox_inches="tight")
plt.savefig("NV_emission_flipped_around_ZPL.png", dpi=300, bbox_inches="tight");


# Compare the [image](https://commons.wikimedia.org/wiki/File:NVple.JPG) from the [Wikipedia page for NV centers](https://en.wikipedia.org/wiki/Nitrogen-vacancy_center):

# In[90]:


IPython.display.Image(filename="NVple.JPG")


# The NV- spectrum has a one-photon absorption cross section specified at 532 nm, so we should normalize the spectrum so that it has the correct absorption at that value.

# In[91]:


i_closest_to_532nm = np.abs(532 - nv_center_absorption_nm).argmin()
intensity_closest_to_532nm = nv_center_emission_intensity[i_closest_to_532nm]
nv_center_absorption_intensity = nv_center_emission_intensity/intensity_closest_to_532nm


# In[92]:


plt.clf()
plt.hlines(y=1.0, xmin=min(nv_center_absorption_nm), xmax=max(nv_center_absorption_nm), color="gray", linestyles="dotted")
plt.scatter(nv_center_absorption_nm, nv_center_absorption_intensity, color="tab:orange", label="NV- center absorption normalized to 532nm")
plt.vlines(x=637, ymin=0, ymax=max(nv_center_absorption_intensity), color="red", label="zero phonon line (637 nm)")
plt.vlines(x=532, ymin=0, ymax=max(nv_center_absorption_intensity), color="green", label="excitation wavelength (532 nm)")
plt.legend(loc="center right", bbox_to_anchor=(1.8,0.9))
plt.xlabel("Photon wavelength (nm)")
plt.ylabel("Intensity (a.u.)")
plt.savefig("NV_absorption_normalized.eps", bbox_inches="tight")
plt.savefig("NV_absorption_normalized.png", dpi=300, bbox_inches="tight");


# Now let's see what the $\mathrm{WS_2}$ emission looks like overlaid on the NV- center absorption spectrum.

# In[93]:


plt.clf()
plt.scatter(WS2_nm, WS2_intensity/WS2_intensity.max(), label="$\mathrm{WS}_2$ emission", color="tab:blue")
plt.scatter(nv_center_absorption_nm, nv_center_absorption_intensity/nv_center_absorption_intensity.max(), color="tab:orange", label="NV- center absorptivity")
plt.xlabel("Photon wavelength (nm)")
plt.legend(loc="best", bbox_to_anchor=(1.5,0.9));
plt.ylabel("Normalized intensity");
plt.savefig("NV_WS2_both_spectra.eps", bbox_inches="tight")
plt.savefig("NV_WS2_both_spectra.png", dpi=300, bbox_inches="tight");


# ## Overlap integral

# In[94]:


def overlap_OLI(
        donor_lambda,
        donor_fluorescence,
        acceptor_lambda,
        acceptor_extinction,
        molar_attenuation_coefficient):
    # In this method, we run over donor lambdas.
    assert len(donor_lambda) == len(donor_fluorescence)
    assert len(acceptor_lambda) == len(acceptor_extinction)
    # Peak-normalize donor fluorescence.
    donor_norm = donor_fluorescence/donor_fluorescence.max()
    # Peak-normalize acceptor.
    acceptor_norm = acceptor_extinction/acceptor_extinction.max()
    # Only go over the overlapping lambda range.
    lambda_min = max(donor_lambda.min(),acceptor_lambda.min())
    lambda_max = min(donor_lambda.max(),acceptor_lambda.max())
    wavelengths = []
    J_raw = []
    for i, wavelength_nm in enumerate(donor_lambda):
        if wavelength_nm < lambda_min:
            continue
        elif wavelength_nm > lambda_max:
            continue
        wavelengths.append(wavelength_nm)
        f_D = donor_norm[i]
        # Find index of closest corresponding wavelength in acceptor_lambda.
        # TODO: issue a warning if the acceptor data is too sparse.
        j = np.abs(acceptor_lambda - wavelength_nm).argmin()
        # eps_A has units of 1/(M cm)
        eps_A = acceptor_norm[j]*molar_attenuation_coefficient
        # OLI is 10^14 M^-1 cm^-1 nm^4
        J_val = 1e-14 * f_D * eps_A * (wavelength_nm**4)
        J_raw.append(J_val)
    lambda_J = np.array(wavelengths)
    J = np.array(J_raw) / donor_norm.sum()
    return lambda_J, J


# We need to relate the absorption cross section $\sigma$ of an NV- center to the molar attenuation coefficient $\epsilon_A$.
# 
# $$\epsilon_A = \frac{N_A}{\ln 10} \sigma$$

# In[95]:


nv_center_absorption_cross_section = 3.1e-17 # cm^2/defect
N_A = 6.02214076e23 # defect/mol
molar_attenuation_coefficient_raw = N_A*nv_center_absorption_cross_section/math.log(10) # cm^2/mol
cm2_per_mol_to_inverse_molar_cm = 1/1000.
molar_attenuation_coefficient = molar_attenuation_coefficient_raw*cm2_per_mol_to_inverse_molar_cm # 1/(M cm)
print("{:.3g} M^-1 cm^-1".format(molar_attenuation_coefficient))


# In[96]:


lambda_J, J_OLI = overlap_OLI(
    WS2_nm,
    WS2_intensity,
    nv_center_absorption_nm,
    nv_center_absorption_intensity,
    molar_attenuation_coefficient)


# In[97]:


J_sum_OLI = J_OLI.sum()
print("J = {:.4g} * 10^14 M^-1 cm^-1 nm^4 (OLI)".format(J_sum_OLI))


# In[98]:


plt.clf()
plt.scatter(lambda_J, J_OLI, label="overlap",color="tab:green")
plt.legend()
plt.xlabel("Photon wavelength (nm)");


# In[99]:


plt.clf()
plt.scatter(WS2_nm, WS2_intensity/WS2_intensity.max(), color="tab:blue", label="$\mathrm{WS}_2$ emission")
plt.scatter(nv_center_absorption_nm, nv_center_absorption_intensity/nv_center_absorption_intensity.max(), color="tab:orange", label="NV- center absorptivity")
plt.scatter(lambda_J, J_OLI/J_OLI.max(), label="overlap", color="tab:green")
plt.legend()
plt.xlim(min(lambda_J),max(lambda_J))
plt.xlabel("Photon wavelength (nm)")
plt.savefig("WS2_NV_overlap.eps", bbox_inches="tight")
plt.savefig("WS2_NV_overlap.png", dpi=300, bbox_inches="tight");


# In[100]:


plt.clf()
plt.scatter(WS2_nm, WS2_intensity/WS2_intensity.max(), color="tab:blue", label="$\mathrm{WS}_2$ emission")
plt.scatter(lambda_J, J_OLI/J_OLI.max(), label="overlap", color="tab:green")
plt.legend()
plt.xlabel("Photon wavelength (nm)");


# ## Förster distance
# 
# Now that we have the overlap integral, we can calculate the Förster distance from a few other parameters of the host material.

# In[101]:


def forster_distance_OLI(kappa2, Phi_D, J, n):
    R6_0 = kappa2*Phi_D*J/(n**4)
    R_0 = 4.542*(R6_0)**(1./6.) # nm
    return R_0


# Recall
# $$
# R_0^6 = \frac{9 (\ln 10) \kappa^2 \Phi_D J }{128\pi^5 n^4 N_A}
# $$
# which for OLI units is
# $$
# R_0 = 4.542 \left(\frac{\kappa^2 \Phi_D J }{n^4}\right)^{1/6}
# $$

# [Peimyoo et al. 2018](https://doi.org/10.1021/nn4046002): "The PL quantum yield of our 1L-WS2 on quartz is ∼2.0%, which is much larger than that (0.42%) of suspended 1L-MoS2 exfoliated from the bulk MoS2 crystal."

# In[102]:


kappa2 = 2./3 # orientation factor for dipole interaction
# "The orientation factor takes on a value of κ² = 2/3 for a dynamic and
# isotropic distribution of donor and acceptor orientations"
# "Varies between 0 and 4."

n_diamond = 2.417 # index of refraction of diamond
# "All refractive index values in the literature are in the 1.33 - 1.6 range.
# The values 1.34 and 1.6 are the ones used most frequently."
# TODO: check if this correct at 532 nm.

Phi_D = 0.02 # quantum yield of WS2.
# Using the 2% value.
# "varies between 0 and 1"

R_0_nm = forster_distance_OLI(kappa2, Phi_D, J = J_sum_OLI, n = n_diamond)

print("R₀ = {:.4} nm".format(R_0_nm))


# ## FRET efficiency

# In[103]:


def FRET_efficiency(R_0, r_DA):
    # Need both in same units.
    R6_0 = R_0**6
    E = R6_0 / (R6_0 + r_DA**6)
    return E


# We're interested in a range of donor-acceptor distances.

# In[104]:


r_DA_range_nm = np.linspace(0.1, 10, num=100)


# In[105]:


FRET_efficiencies_list = []
for r_DA_nm in r_DA_range_nm:
    E = FRET_efficiency(R_0_nm, r_DA_nm)
    FRET_efficiencies_list.append(E)
FRET_efficiencies = np.array(FRET_efficiencies_list)


# In[106]:


plt.clf()
plt.hlines(y=0.5, xmin=min(r_DA_range_nm), xmax=max(r_DA_range_nm), color="gray", linestyles="dotted")
plt.vlines(x=R_0_nm, ymin=0, ymax=1, color="black", linestyles="dashed", label="Förster distance ({:.3g} nm)".format(R_0_nm))
plt.scatter(r_DA_range_nm, FRET_efficiencies, color="tab:purple")
plt.xlabel("$r_{DA}$ (nm)")
plt.legend()
plt.ylabel("FRET efficiency")
plt.savefig("FRET_efficiency.eps", bbox_inches="tight")
plt.savefig("FRET_efficiency.png", dpi=300, bbox_inches="tight");


# # Break-even point
# 
# Finally, we wish to compare direct excitation of the two-photon absorption of the NV center at 1064 nm with the indirect excitation of the NV center via FRET from the $\mathrm{WS_2}$ two-photon absorption. The $\mathrm{WS_2}$ two-photon absorption cross section is about four orders of magnitude larger than that of the NV center.
# 
# The FRET efficiency drops off as a function of distance, and once we reach the distance where the FRET efficiency is the ratio of the absorption coefficients, i.e. $\approx 1 \times 10^{-4}$, then we have reached the "break-even" point where FRET becomes equally as inefficient as direct excitation.

# To do this, we will need $r_{DA}$ in terms of $E$ and $R_0$, i.e.
# 
# $$
# r_{DA} = R_0 \left(\frac{1-E}{E}\right)^{1/6}
# $$

# We will also need $E$ from the ratio of the two-photon cross sections,
# 
# $$
# E = \frac{ \sigma^{(2)}_{A}}
#          { \sigma^{(2)}_{B}}
# $$
# where $\sigma^{(2)}_{A}$ and $\sigma^{(2)}_{B}$ are the two-photon absorption cross-sections of the NV-center and the $\mathrm{WS_2}$ film respectively, generally in units of $\mathrm{cm^4 s/photon}$. Note that $\sigma^{(2)}$ depends on the wavelength.

# From [Wee et al. 2007](https://doi.org/10.1021/jp073938o) we have $\sigma^{(2)}_{A} = (0.45 \pm 0.23) \times 10^{-50} \mathrm{cm^4} \cdot \mathrm{s / photon}$ for the NV center, measured at a wavelengths of 1064 nm against a Rhodamine B standard.

# From [Zhang et al. 2015](https://doi.org/10.1021/acsnano.5b03480) we have $\mathrm{WS_2}$, we have the two-photon absorption coefficient $\beta = 1.0 \pm 0.8 \times 10^4 \mathrm{cm/GW}$, measured at 1030 nm via fitting Z-scan curves. To transform this into an absorption cross section, we use
# $$
# \sigma^{(2)} = \frac{\varepsilon}{N} \beta
# $$
# where $\sigma^{(2)}$ is the two-photon absorption cross-section at a given wavelength, $\varepsilon = \frac{h c}{\lambda}$ is the energy of a single photon at that same wavelength, $N$ is the number density of active atoms, and $\beta$ is the two-photon absorption coefficient (also at the same wavelength).

# From [Zhou et al. 2018](https://doi.org/10.1364/OE.26.016093) we have $N = 1.88 \times 10^{22} \, \mathrm{cm^{-3}}$, so we find
# \begin{align*}
# \sigma^{(2)} &= \frac{h c}{\lambda N} \beta \\
#              &= \frac{h c}{1030 nm \cdot 1.88 \times 10^{22} \mathrm{cm}^{-3}} 1.0 \pm 0.8 \times 10^4 \mathrm{cm/GW} \\
#              &= 1.03 \times 10^{-46} \mathrm{cm^4 s/photon}
# \end{align*}

# In[107]:


def get_sigma_TPA(beta_cm_per_GW, N_per_cm3, wavelength_nm):
    h = 6.62607015e-34 # J s
    m_to_cm = 1e2
    c = 299792458*m_to_cm # cm/s
    nm_to_cm = 1e-7
    E = h*c/(wavelength_nm*nm_to_cm) # J
    GW_to_W = 1e-9
    beta = beta_cm_per_GW*GW_to_W # cm/W
    # N is already in cm^{-3}
    sigma_TPA = E * beta / N_per_cm3 # cm^4 s
    return sigma_TPA

sigma_TPA_NV = 0.45e-50 # cm^4 s/photon
sigma_TPA_WS2 = get_sigma_TPA(beta_cm_per_GW = 1.0e4, N_per_cm3 = 1.88e22, wavelength_nm = 1030)
print("σ_TPA_NV- = {:.2e} cm⁴ s/photon".format(sigma_TPA_NV))
print("σ_TPA_WS2 = {:.2e} cm⁴ s/photon".format(sigma_TPA_WS2))


# In[108]:


E_ratio = sigma_TPA_NV / sigma_TPA_WS2

def get_r_DA(E, R_0):
    return R_0*((1-E)/E)**(1./6.)

r_DA_breakeven = get_r_DA(E_ratio, R_0_nm)

print("E = {:.4g}".format(E_ratio))
print("r_DA = {:.1f} nm (break-even point)".format(r_DA_breakeven))


# In[109]:


Phi_D_range = np.linspace(0.01, 1.0, num=100, endpoint=True)
R_0_nm_range = forster_distance_OLI(kappa2, Phi_D_range, J = J_sum_OLI, n = n_diamond)
r_DA_breakeven = get_r_DA(E_ratio, R_0_nm_range)


# In[110]:


plt.clf()
plt.scatter(Phi_D_range, r_DA_breakeven, color="tab:pink")
plt.title("Dependence of break-even distance on quantum yield of $\mathrm{WS_2}$")
plt.xlabel("$\Phi_D$ (dimensionless)")
plt.ylabel("$r_{DA}$ (nm)")
plt.savefig("break-even_distance.eps", bbox_inches="tight")
plt.savefig("break-even_distance.png", dpi=300, bbox_inches="tight");


# # Appendix

# Public git repository for this notebook: <https://github.com/nbeaver/ws2-diamond-fret-calculation>
