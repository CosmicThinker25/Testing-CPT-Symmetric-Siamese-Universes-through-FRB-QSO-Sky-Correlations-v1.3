# -*- coding: utf-8 -*-
"""
frb_cpt_dm_histogram.py
Genera un histograma comparativo de DM por hemisferio respecto al eje CPT.
"""

import argparse, json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True)
    ap.add_argument('--colmap', default='{"ra":"ra","dec":"dec","dm":"dm_fitb"}')
    ap.add_argument('--bmin', type=float, default=20.0)
    ap.add_argument('--dmmin', type=float, default=800.0)
    ap.add_argument('--ra0', type=float, default=170.0)
    ap.add_argument('--dec0', type=float, default=40.0)
    return ap.parse_args()

def load_filtered(csv_path, colmap, bmin, dmmin):
    df = pd.read_csv(csv_path)
    df = df[[colmap['ra'], colmap['dec'], colmap['dm']]].copy()
    df.columns = ['ra', 'dec', 'dm']
    df['ra'] = pd.to_numeric(df['ra'], errors='coerce')
    df['dec'] = pd.to_numeric(df['dec'], errors='coerce')
    df['dm'] = pd.to_numeric(df['dm'], errors='coerce')
    df = df.dropna(subset=['ra','dec','dm'])
    sc = SkyCoord(ra=df['ra'].values*u.deg, dec=df['dec'].values*u.deg, frame='icrs')
    b = sc.galactic.b.deg
    mask = (np.abs(b)>=bmin) & (df['dm']>=dmmin)
    return df.loc[mask,'ra'].values, df.loc[mask,'dec'].values, df.loc[mask,'dm'].values

def unit_vectors(ra,dec):
    ra,dec=np.deg2rad(ra),np.deg2rad(dec)
    x=np.cos(dec)*np.cos(ra); y=np.cos(dec)*np.sin(ra); z=np.sin(dec)
    return np.vstack([x,y,z]).T

def hemispheres(nvecs,dm,n0):
    sgn=np.sign(nvecs@n0)
    return dm[sgn>=0], dm[sgn<0]

def main():
    args=parse_args()
    colmap=json.loads(args.colmap)
    ra,dec,dm=load_filtered(args.csv,colmap,args.bmin,args.dmmin)
    nvecs=unit_vectors(ra,dec)
    n0=unit_vectors(np.array([args.ra0]),np.array([args.dec0]))[0]
    dm_pos,dm_neg=hemispheres(nvecs,dm,n0)

    plt.figure(figsize=(8,5))
    bins=np.linspace(min(dm),max(dm),25)
    plt.hist(dm_pos,bins=bins,color='royalblue',alpha=0.6,label=f'Hem. + ({len(dm_pos)})')
    plt.hist(dm_neg,bins=bins,color='indianred',alpha=0.6,label=f'Hem. – ({len(dm_neg)})')

    plt.axvline(np.median(dm_pos),color='blue',linestyle='--',lw=1.5)
    plt.axvline(np.median(dm_neg),color='red',linestyle='--',lw=1.5)

    plt.xlabel('Dispersion Measure (pc cm⁻³)')
    plt.ylabel('Número de FRBs')
    plt.title(f'Comparativa de DMs por hemisferio (RA₀={args.ra0}°, Dec₀={args.dec0}°)')
    plt.legend()
    plt.tight_layout()
    plt.savefig('frb_dm_histogram.png',dpi=200)
    plt.show()

    print(f"\nMediana hemisferio +: {np.median(dm_pos):.2f}")
    print(f"Mediana hemisferio –: {np.median(dm_neg):.2f}")
    print("→ Figura guardada como frb_dm_histogram.png")

if __name__=="__main__":
    main()
