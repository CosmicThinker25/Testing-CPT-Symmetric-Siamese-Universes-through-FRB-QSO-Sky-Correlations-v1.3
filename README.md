# 🌀 FRB–Siamese Universes v1.3: Final CPT Axis Analysis

**Author:** [CosmicThinker](https://github.com/CosmicThinker25)  
**Affiliation:** Independent Researcher  
**Date:** October 2025  
**Version:** 1.3 (Final CPT)  
**Status:** Submitted to Zenodo (DOI pending)

---

## 🌌 Abstract

This repository contains the data products, figures, and analysis scripts accompanying  
the research note **“Testing CPT–Symmetric Siamese Universes through FRB–QSO Sky Correlations (v1.3)”**.  

The study performs a cross-correlation analysis between **Fast Radio Bursts (FRBs)** from the CHIME/FRB Catalog 1  
and **Quasars (QSOs)** from the DESI DR1 release, searching for directional anisotropies  
consistent with a *CPT-symmetric mirror universe framework*.  

After applying $|b|>20^\circ$ and $\mathrm{DM}\ge800$ filters to 123 FRBs and $\sim4.7\times10^5$ QSOs  
($z<3$), the sky appears globally isotropic ($R^2\simeq10^{-4}$, $p\simeq0.995$),  
yet exhibits a **localized excess of supra-threshold clusters ($Z=7.0\pm1.1$)**  
within $15^\circ$ of the proposed CPT axis $(\alpha_0,\delta_0)=(170^\circ,40^\circ)$.  
Monte Carlo rotation trials ($N_\mathrm{rot}=1000$) yield a global significance  
of $p_\mathrm{global}\approx0.019$, which is *marginal but intriguing*.  

The result motivates deeper analyses using **CHIME Cat 2**, **ASKAP**, and future multi-tracer catalogs  
to confirm or refute possible CPT–symmetric hemispheric asymmetries.

---

## 🧩 Repository Structure

FRB_Siamese_Universes_v1.3_finalCPT/
│
├── FRB_Siamese_Universes_v1.3_finalCPT.pdf # Full research paper
├── frb_qso_heatmap_axis.py # FRB–QSO cross-correlation map
├── frb_qso_autocorr_sigma.py # σ-map and autocorrelation
├── frb_sigma_radial_density.py # Radial Z–density analysis
├── frb_axis_rotation_test.py # Random-axis (look-elsewhere) correction
├── frb_qso_mirror_overlay.py # Mirror-sky reflection mapping
├── data/
│ ├── chimefrbcat1.csv
│ ├── desi_qso_dr1.csv
│ └── frb_dipole_results.json
├── figures/
│ ├── frb_sky_scatter_filtered.png
│ ├── frb_qso_mirror_overlay.png
│ ├── frb_sigma_axis_map.png
│ ├── Figure_1A.png ... Figure_1F.png
│ └── ...
└── README.md
---

## 🧠 Methodology Overview

The analysis integrates:
- **Isotropy test:** Dipole regression of FRB dispersion measures (DM).  
- **Mirror mapping:** Reflection of coordinates about a fixed CPT axis.  
- **σ–map:** Angular field of normalized fluctuations $Z_i = (N_i - \langle N \rangle)/\sqrt{\langle N \rangle}$.  
- **Cluster statistics:** Supra-threshold detection at $|Z| \ge 2.63$.  
- **Axis rotation test:** Monte Carlo randomization ($N_\mathrm{rot}=1000$) for global $p$-value.  
- **Visualization:** Full-sky and directional significance maps.

The angular coherence scale ($\theta_c\approx6^\circ$) corresponds to $\sim0.17$ Gpc at $z\!\approx\!1$  
under $\Lambda$CDM cosmology ($D_A \simeq 1.6$ Gpc).

---

## 📊 Key Results

| Metric | Value | Notes |
|:--------------------------|:-----------------|:-------------------------------|
| $R^2$ (isotropy) | $1\times10^{-4}$ | Null dipole result |
| $p$ (isotropy test) | 0.995 | Global isotropy |
| $\mathrm{DM}_+$ median | 1089.6 pc cm⁻³ | Northern CPT hemisphere |
| $\mathrm{DM}_-$ median | 1231.2 pc cm⁻³ | Southern CPT hemisphere |
| $N_{\rm clusters}$ | 35 | |Z| ≥ 2.63 |
| Excess near axis | $Z = 7.0 ± 1.1$ | Bootstrap |
| $p_{\rm global}$ | 0.019 | 1000 random axes |

---

## 🔬 Dependencies

Developed and tested in **Python 3.11** with:
```bash
pip install numpy pandas scipy matplotlib astropy tqdm

📂 Data and Code Availability

All scripts, figures, and processed CSV/JSON files are included in this repository.
Raw data sources:

CHIME/FRB Catalog 1: CHIME Collaboration (2021), ApJS, 257, 59

DESI DR1 QSOs: DESI Collaboration (2025), arXiv:2503.14745

Reproduction of figures requires running the Python scripts sequentially from the main directory.

🪞 Acknowledgments

This project was co-developed with ChatGPT (OpenAI) as an analytical collaborator and editor.
Dedicated “to the pursuit of symmetry beyond the visible cosmos.”
