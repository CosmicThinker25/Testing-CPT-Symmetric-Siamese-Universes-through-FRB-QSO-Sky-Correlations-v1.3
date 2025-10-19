#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FRB–QSO Siamese — Directional Axis Significance Map (v1.2)
Autores: CosmicThinker & ChatGPT (Toko)
Versión: 1.2 — Octubre 2025

Genera un mapa de densidad de clústeres σ en todo el cielo
al contar cuántos picos caen dentro de 15° de cada dirección.
Así se visualiza dónde el cielo muestra mayor anisotropía.

Entradas:
  - frb_sigma_clusters.csv  (con centroid_ra_deg, centroid_dec_deg)
Salidas:
  - frb_sigma_axis_map.png  (Mollweide)
  - frb_sigma_axis_map.json  (valores numéricos)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json, math

# ---------------- Parámetros ----------------
CSV = "frb_sigma_clusters.csv"
THETA_MAX = 15.0   # radio angular en grados
F_SKY = 0.66
GRID_RA = 72       # nº de pasos en RA (360° / GRID_RA = resolución)
GRID_DEC = 36      # nº de pasos en Dec (-90° a +90°)
RA0, DEC0 = 170.0, 40.0  # eje CPT
# --------------------------------------------

df = pd.read_csv(CSV)
ra = np.deg2rad(df["centroid_ra_deg"].values)
dec = np.deg2rad(df["centroid_dec_deg"].values)
N = len(df)

# --- función de distancia angular entre dos direcciones ---
def ang_sep(ra1, dec1, ra2, dec2):
    return np.degrees(np.arccos(
        np.clip(np.sin(dec1)*np.sin(dec2) +
                np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2), -1, 1)
    ))

# --- mapa de densidades ---
ras = np.linspace(0, 360, GRID_RA)
decs = np.linspace(-90, 90, GRID_DEC)
density = np.zeros((GRID_DEC, GRID_RA))

for i, d in enumerate(decs):
    for j, r in enumerate(ras):
        theta = ang_sep(ra, dec, np.deg2rad(r), np.deg2rad(d))
        density[i, j] = np.sum(theta < THETA_MAX)

# --- normalización ---
DEG2_PER_SR = (180.0 / math.pi) ** 2
A_tot_deg2 = 4.0 * math.pi * DEG2_PER_SR * F_SKY
A_in_deg2 = 2.0 * math.pi * (1 - math.cos(np.deg2rad(THETA_MAX))) * DEG2_PER_SR * F_SKY
mu_in = N * (A_in_deg2 / A_tot_deg2)
sigma_in = np.sqrt(mu_in)
Z_map = (density - mu_in) / sigma_in

# --- guardar datos ---
data_out = {
    "theta_max_deg": THETA_MAX,
    "mu_in": float(mu_in),
    "sigma_in": float(sigma_in),
    "Z_map": Z_map.tolist(),
    "ras": ras.tolist(),
    "decs": decs.tolist()
}
with open("frb_sigma_axis_map.json", "w") as f:
    json.dump(data_out, f, indent=2)

# --- figura ---
plt.figure(figsize=(10, 5))
ax = plt.subplot(111, projection="mollweide")

# conversión a radianes centrada en RA=0
RA_plot = np.radians(ras - 180)
DEC_plot = np.radians(decs)
Z_plot = np.flipud(Z_map)  # Dec va de +90 a -90

im = ax.pcolormesh(RA_plot, DEC_plot, Z_plot, cmap="coolwarm", shading="auto", vmin=-4, vmax=4)
plt.colorbar(im, orientation="horizontal", pad=0.1, label="Directional significance Z")

# marcar el eje CPT
ax.plot(np.radians(RA0-180), np.radians(DEC0), 'yo', markersize=8, label="CPT (+)")
ax.plot(np.radians((RA0+180)-180), np.radians(-DEC0), 'ko', markersize=6, label="CPT (−)")

ax.grid(True, color="gray", alpha=0.5)
ax.set_title(f"FRB–QSO Siamese Directional Significance Map (θ≤{THETA_MAX}°)")
ax.legend(loc="lower left")
plt.tight_layout()
plt.savefig("frb_sigma_axis_map.png", dpi=250)
plt.show()
