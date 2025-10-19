#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FRB–QSO Siamese Mirror Overlay (v1.3)
Autores: CosmicThinker & ChatGPT (Toko)

Qué hace:
- Carga FRBs y QSOs.
- Aplica cortes básicos (DM >= 800, |b| > bcut, z < 3).
- Refleja los FRBs con respecto a un eje CPT (ra0, dec0).
- Dibuja:
  * Fondo de densidad QSO (pcolormesh) en Mollweide.
  * FRBs reales (blancos) y reflejados (verdes).
  * Eje CPT (+) y su antipoda (–).
- Guarda: frb_qso_mirror_overlay.png
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------- utilidades ----------
def reflect_cpt(ra_deg, dec_deg, ra0_deg, dec0_deg):
    """
    Reflexión siamés respecto a (ra0, dec0):
      ra'  = 2*ra0 - ra   (mod 360)
      dec' = 2*dec0 - dec (clamp [-90,90])
    """
    ra_ref = (2*ra0_deg - ra_deg) % 360.0
    dec_ref = 2*dec0_deg - dec_deg
    dec_ref = np.clip(dec_ref, -90.0, 90.0)
    return ra_ref, dec_ref

def radec_to_mollweide(ra_deg, dec_deg):
    """Convierte a radianes y centra RA en 0° → eje central de Mollweide."""
    return np.radians(ra_deg - 180.0), np.radians(dec_deg)

def galactic_lat_from_equatorial(ra_deg, dec_deg):
    """Aproximación rápida de |b| usando rotación estándar J2000 (lo suficiente para cortes)."""
    ra = np.radians(ra_deg); dec = np.radians(dec_deg)
    # Parámetros IAU 1958 / J2000
    ra_ngp = np.radians(192.85948)
    dec_ngp = np.radians(27.12825)
    l_omega = np.radians(32.93192)
    sinb = (np.sin(dec)*np.sin(dec_ngp) +
            np.cos(dec)*np.cos(dec_ngp)*np.cos(ra - ra_ngp))
    b = np.arcsin(np.clip(sinb, -1, 1))
    return np.degrees(b)

# ---------- main ----------
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="FRB–QSO Siamese Mirror Overlay")
    ap.add_argument("--frb", required=True, help="CSV de FRBs (p.ej., chimefrbcat1.csv)")
    ap.add_argument("--qso", required=True, help="CSV de QSOs (p.ej., quasar_catalog.csv)")
    ap.add_argument("--ra0", type=float, default=170.0, help="RA del eje CPT (deg)")
    ap.add_argument("--dec0", type=float, default=40.0, help="Dec del eje CPT (deg)")
    ap.add_argument("--bcut", type=float, default=20.0, help="Corte en |b| (deg)")
    ap.add_argument("--qso_bins_ra", type=int, default=144, help="Bins RA para densidad QSO")
    ap.add_argument("--qso_bins_dec", type=int, default=72, help="Bins Dec para densidad QSO")
    ap.add_argument("--dpi", type=int, default=220, help="DPI de la figura")
    args = ap.parse_args()

    # --- carga datos ---
    frb = pd.read_csv(args.frb)
    qso = pd.read_csv(args.qso)
    frb.columns = frb.columns.str.lower()
    qso.columns = qso.columns.str.lower()

    # --- cortes básicos ---
    # FRB: DM >= 800, |b| > bcut
    b_frb = galactic_lat_from_equatorial(frb["ra"].values, frb["dec"].values)
    frb_sel = frb[(frb["dm_fitb"] >= 800) & (np.abs(b_frb) >= args.bcut)].copy()

    # QSO: z < 3, |b| > bcut
    b_qso = galactic_lat_from_equatorial(qso["ra"].values, qso["dec"].values)
    qso_sel = qso[(qso["z"] < 3.0) & (np.abs(b_qso) >= args.bcut)].copy()

    print(f"FRBs usados: {len(frb_sel)} | QSOs usados: {len(qso_sel)} | |b|≥{args.bcut}°")

    # --- reflejo siamés de FRBs ---
    ra_ref, dec_ref = reflect_cpt(frb_sel["ra"].values, frb_sel["dec"].values,
                                  args.ra0, args.dec0)

    # --- densidad QSO en rejilla ---
    ra_bins = np.linspace(0, 360, args.qso_bins_ra + 1)
    dec_bins = np.linspace(-90, 90, args.qso_bins_dec + 1)
    H, yedges, xedges = np.histogram2d(qso_sel["dec"].values, qso_sel["ra"].values,
                                       bins=[dec_bins, ra_bins])
    # normalización log-suave
    H = np.log10(1 + H)

    # --- preparar proyección Mollweide ---
    ra_cent = 0.5*(ra_bins[:-1] + ra_bins[1:])
    dec_cent = 0.5*(dec_bins[:-1] + dec_bins[1:])
    RA_grid, DEC_grid = np.meshgrid(np.radians(ra_cent - 180.0),
                                    np.radians(dec_cent))
    H_plot = H  # ya con forma (dec, ra)

    # --- figura ---
    plt.figure(figsize=(10, 5))
    ax = plt.subplot(111, projection="mollweide")

    # Fondo: densidad QSO
    im = ax.pcolormesh(RA_grid, DEC_grid, H_plot,
                       cmap="magma", shading="auto", alpha=0.55)

    # FRBs reales (blancos) y reflejados (verdes)
    ax.scatter(*radec_to_mollweide(frb_sel["ra"].values, frb_sel["dec"].values),
               s=12, c="white", alpha=0.85, label="Real FRBs", linewidths=0)
    ax.scatter(*radec_to_mollweide(ra_ref, dec_ref),
               s=10, c="lime", alpha=0.85, label="Mirrored FRBs", linewidths=0)

    # Eje CPT (+) y antipoda (–)
    ax.plot(np.radians(args.ra0 - 180.0), np.radians(args.dec0),
            "o", color="gold", markersize=9, label="CPT (+)")
    ax.plot(np.radians((args.ra0 + 180.0) % 360 - 180.0), np.radians(-args.dec0),
            "o", color="black", markersize=7, label="CPT (−)")

    # Estética
    ax.grid(True, color="gray", alpha=0.5)
    ax.set_title("FRB–QSO Siamese Mirror Overlay", fontsize=13)
    cb = plt.colorbar(im, orientation="horizontal", pad=0.08, fraction=0.05)
    cb.set_label("QSO surface density (log scale)")

    ax.legend(loc="lower left", fontsize=8)
    plt.tight_layout()
    out = "frb_qso_mirror_overlay.png"
    plt.savefig(out, dpi=args.dpi)
    print(f"✅ Guardado: {out}")
