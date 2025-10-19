#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Autocorrelación angular y mapa de significancia σ
para el campo de asimetría FRB–QSO siamés.
Autores: CosmicThinker & ChatGPT (Toko)
Versión: 1.0 — Octubre 2025
"""

import numpy as np
import matplotlib.pyplot as plt
import json, math
from scipy.fft import fft2, ifft2, fftshift
from scipy.ndimage import uniform_filter1d
import pandas as pd

# --- Cargar rejilla ---
with open("frb_qso_heatmap_grid.json", "r") as f:
    data = json.load(f)

arr = np.array(data["heatmap_norm"], dtype=np.float32)
arr -= np.nanmean(arr)
arr /= np.nanstd(arr) + 1e-8

h, w = arr.shape
px_per_deg = w / 360.0

# --- Autocorrelación 2D por FFT ---
F = fft2(arr)
AC = fftshift(np.real(ifft2(np.abs(F) ** 2)))
AC /= AC.max()

# --- Perfil radial ---
y, x = np.indices(arr.shape)
r = np.hypot(y - h / 2, x - w / 2)
bins = np.arange(0, int(np.min([h, w]) / 2))
radial = np.zeros_like(bins, dtype=np.float64)
for i in range(len(bins) - 1):
    mask = (r >= bins[i]) & (r < bins[i + 1])
    if np.any(mask):
        radial[i] = np.nanmean(AC[mask])

deg = bins / px_per_deg
corr = uniform_filter1d(radial, size=5)

def crossing_angle(deg, corr, target):
    idx = np.where(corr <= target)[0]
    if len(idx) == 0:
        return np.nan
    i = idx[0]
    if i > 0:
        x0, x1 = deg[i - 1], deg[i]
        y0, y1 = corr[i - 1], corr[i]
        return x0 + (target - y0) * (x1 - x0) / (y1 - y0)
    return deg[i]

theta_1e = crossing_angle(deg, corr, 1 / math.e)
theta_half = crossing_angle(deg, corr, 0.5)

# --- Parámetros físicos ---
f_sky = 0.66
def neff(theta_deg): return 4.0 * f_sky / (np.deg2rad(theta_deg) ** 2)
def L_gly(theta_deg, D_gly=11.0): return np.deg2rad(theta_deg) * D_gly

N1, N2 = neff(theta_1e), neff(theta_half)
L1, L2 = L_gly(theta_1e), L_gly(theta_half)

# --- σ necesario para 1 discrepancia ---
from scipy.stats import norm
def sigma_for_one(N):
    p = 1.0 / N
    z = norm.isf(p)  # unilateral
    return float(z)

sigma_1 = sigma_for_one(N1)
sigma_2 = sigma_for_one(N2)

# --- Resultados ---
results = {
    "theta_c_1e_deg": float(theta_1e),
    "theta_c_half_deg": float(theta_half),
    "N_eff_1e": float(N1),
    "N_eff_half": float(N2),
    "L_phys_1e_Gly": float(L1),
    "L_phys_half_Gly": float(L2),
    "sigma_for_one_discrepancy_1e": sigma_1,
    "sigma_for_one_discrepancy_half": sigma_2
}
with open("frb_autocorr_summary.json", "w") as f:
    json.dump(results, f, indent=2)

# --- Gráfica ---
plt.figure(figsize=(6, 4))
plt.plot(deg, corr, color="royalblue")
plt.axhline(1 / math.e, color="gray", ls="--", lw=1)
plt.axhline(0.5, color="gray", ls=":", lw=1)
plt.xlabel("Angular separation (degrees)")
plt.ylabel("Autocorrelation C(θ)")
plt.title("Autocorrelation of Siamese Asymmetry Field")
plt.grid(True)
plt.tight_layout()
plt.savefig("frb_autocorr_profile.png", dpi=200)
plt.show()

print("\n=== Resultados ===")
for k, v in results.items():
    print(f"{k:35s}: {v:.4f}")
print("\nGráficas y resumen guardados.")
