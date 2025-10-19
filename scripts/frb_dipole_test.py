
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def sph_to_cart(ra_deg, dec_deg):
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    x = np.cos(dec) * np.cos(ra)
    y = np.cos(dec) * np.sin(ra)
    z = np.sin(dec)
    return np.vstack([x, y, z]).T

def fit_dipole(nvecs, y, w=None):
    # Fit y = a + b * (n dot d) with |d|=1
    N = len(y)
    if w is None:
        w = np.ones(N)
    w = np.asarray(w)

    # Fibonacci sphere directions
    num_dirs = 4000
    k = np.arange(num_dirs)
    golden = (np.sqrt(5)-1)/2 * 2*np.pi
    z = 1 - 2*(k+0.5)/num_dirs
    r = np.sqrt(1 - z*z)
    lon = k * golden
    d_dirs = np.vstack([r*np.cos(lon), r*np.sin(lon), z]).T

    y = np.asarray(y)
    y_mean = np.average(y, weights=w)
    yc = y - y_mean

    best = None
    best_ss = np.inf

    for d in d_dirs:
        ndot = nvecs @ d
        X = np.column_stack([np.ones(N), ndot])
        # Weighted least squares (diagonal weights)
        WX = X * w[:, None]
        Wy = yc * w
        beta, *_ = np.linalg.lstsq(WX, Wy, rcond=None)
        a_c, b = beta
        a = a_c + y_mean
        yhat = a + b * ndot
        ss = np.sum(w * (y - yhat)**2)
        if ss < best_ss:
            best_ss = ss
            best = (a, b, d, yhat)

    a, b, d_hat, yhat = best
    ss_tot = np.sum(w * (y - np.average(y, weights=w))**2)
    r2 = 1 - best_ss / ss_tot if ss_tot > 0 else np.nan
    ra = (np.rad2deg(np.arctan2(d_hat[1], d_hat[0])) + 360.0) % 360.0
    dec = np.rad2deg(np.arcsin(d_hat[2]))
    return dict(a=a, b=b, ra=ra, dec=dec, r2=r2), yhat

def permutation_pval(nvecs, y, w, obs_b, nperm=200, seed=123):
    rng = np.random.default_rng(seed)
    cnt = 0
    for _ in range(nperm):
        idx = rng.permutation(len(y))
        res, _ = fit_dipole(nvecs[idx], y, w)
        if abs(res["b"]) >= abs(obs_b):
            cnt += 1
    return (cnt + 1) / (nperm + 1)

def main():
    ap = argparse.ArgumentParser(description="FRB Dipole fit for RM and DM")
    ap.add_argument("--csv", required=True, help="Input CSV with FRB data")
    ap.add_argument("--colmap", default=None, help="JSON mapping for column names")
    ap.add_argument("--min_dm", type=float, default=None)
    ap.add_argument("--max_dm", type=float, default=None)
    ap.add_argument("--abs_rm_lt", type=float, default=None, help="Keep |RM| < this")
    ap.add_argument("--nperm", type=int, default=200)
    ap.add_argument("--seed", type=int, default=123)
    args = ap.parse_args()

    df = pd.read_csv(args.csv)

    cols = {
        "ra": "ra", "dec": "dec",
        "rm": "rm", "rm_err": "rm_err",
        "dm": "dm", "dm_err": "dm_err"
    }
    if args.colmap:
        cols.update(json.loads(args.colmap))

    if args.min_dm is not None:
        df = df[df[cols["dm"]] >= args.min_dm]
    if args.max_dm is not None:
        df = df[df[cols["dm"]] <= args.max_dm]
    if args.abs_rm_lt is not None:
        df = df[df[cols["rm"]].abs() < args.abs_rm_lt]

    df = df.dropna(subset=[cols["ra"], cols["dec"]]).copy()

    nvecs = sph_to_cart(df[cols["ra"]].values, df[cols["dec"]].values)

    results = {}
    for key in ["rm", "dm"]:
        if key not in cols or cols[key] not in df.columns:
            continue
        y = df[cols[key]].values.astype(float)
        w = None
        wcol = cols.get(key + "_err")
        if wcol in df.columns:
            err = df[wcol].values.astype(float)
            err[~np.isfinite(err) | (err <= 0)] = np.nan
            if np.isfinite(err).any():
                w = 1.0 / np.nan_to_num(err, nan=np.nanmedian(err))**2
        if w is None:
            w = np.ones_like(y)

        res, yhat = fit_dipole(nvecs, y, w)
        pval = permutation_pval(nvecs, y, w, obs_b=res["b"], nperm=args.nperm, seed=args.seed)
        res["p_perm"] = pval
        results[key] = res

    # Save results
    Path("frb_dipole_results.json").write_text(json.dumps(results, indent=2))

    # Plot sky scatter colored by RM if available, else DM
    fig = plt.figure(figsize=(8,6))
    colorcol = cols["rm"] if cols["rm"] in df.columns else (cols["dm"] if cols["dm"] in df.columns else None)
    if colorcol:
        sc = plt.scatter(df[cols["ra"]], df[cols["dec"]], c=df[colorcol], s=14)
        plt.colorbar(sc, label=colorcol.upper())
    else:
        plt.scatter(df[cols["ra"]], df[cols["dec"]], s=14)
    plt.xlabel("RA [deg]"); plt.ylabel("Dec [deg]")
    plt.title("FRB Sky")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig("frb_sky_scatter.png", dpi=160)
    plt.close(fig)

    print("=== Dipole fits ===")
    for k, res in results.items():
        print(f"[{k.upper()}] a={res['a']:.3f}, b={res['b']:.3f}, "
              f"dir=(RA,Dec)=({res['ra']:.1f}°, {res['dec']:.1f}°), "
              f"R2={res['r2']:.4f}, p_perm={res['p_perm']:.4f}")
    print("Saved: frb_dipole_results.json, frb_sky_scatter.png")

if __name__ == "__main__":
    main()
