#!/usr/bin/env python3
"""Growth rates of several axial harmonics from one broadband pair of runs (scripts/analyse_broadband.sh).

Reads
  * the level-set twin-difference export (export_ls: perturbation.bin/txt, times.txt), which is decomposed per axial harmonic ("harmonic-resolved DMD"): the base state is invariant under translation along the axis, so the linearized operator commutes with it and its eigenmodes are Fourier harmonics in z; the export is transformed along z and an exact DMD is run on the complex coefficients c_n(r, t) of each harmonic separately. This removes the ill-conditioning of a single DMD on the whole field, in which the time series of harmonics with close growth rates are nearly collinear over the window and are lost in the rank truncation;
  * the DMD archive of the narrow-band level-set twin difference (eigenvalues.txt, modes.bin, modes.txt): each growing, near-stationary mode is assigned the axial harmonic n of the box that carries most of its energy (FFT of the mode along z, rows of the band);
  * optionally the interface records fitted per harmonic (records/hNN.csv from run_lsa record);
and tabulates, per harmonic n with k_n R0 = 2 pi n R0 / L < 1: the DMD rate, the interface rate, Chandrasekhar's exact viscous rate and the two-fluid viscous rate (Tomotika 1935), with a figure of the measured rates on the theoretical curves.

Usage:
    python3 plotting/broadband_spectrum.py --archive runs/analysis/bb72_08/broadband/dmd_ls --records runs/analysis/bb72_08/broadband/records --zlen 72 --nharm 16 --outdir runs/analysis/bb72_08/broadband
"""
import argparse
import csv
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dispersion  # noqa: E402
from plot_mode_shape_archive import growth_rate_two_fluid, read_keyvals  # noqa: E402


def dmd_modes(archive, nharm, freq_rel=0.1):
    """(n, rate, amplitude) for the growing near-stationary DMD modes, the harmonic n from the axial FFT of the mode."""
    info = read_keyvals(os.path.join(archive, "modes.txt"))
    nsp, nm = int(info["n_space"]), int(info["n_modes"])
    nx, ny, nc = (int(v) for v in info["grid_nx_ny_ncomp"].split())
    m = np.fromfile(os.path.join(archive, "modes.bin"), "<f8").reshape(nm, 2, nsp)
    ev = np.loadtxt(os.path.join(archive, "eigenvalues.txt"), comments="#", ndmin=2)
    out = []
    rows_per_z = nsp // ny            # the band keeps the same radial cells at every z
    for j in range(nm):
        re, im, amp = ev[j, 1], ev[j, 2], ev[j, 5]
        if not (re > 0 and abs(im) <= freq_rel * re):
            continue
        phi = (m[j, 0] + 1j * m[j, 1]).reshape(ny, rows_per_z)
        power = np.sum(np.abs(np.fft.fft(phi, axis=0)) ** 2, axis=1)[1:nharm + 1]
        n = int(np.argmax(power)) + 1
        purity = power[n - 1] / np.sum(power)
        out.append((n, re, amp, purity))
    return out


def harmonic_dmd(export, nharm, zlen, t_max=None, skip_frac=0.1, band=2.0, r0=1.0, tol=1e-10, freq_rel=0.1):
    """Harmonic-resolved DMD of a twin-difference export: for each axial harmonic n, the dominant growing near-stationary rate of the complex Fourier coefficient field c_n(r, t), restricted to the band |r - r0| < band * dr."""
    info = {}
    with open(os.path.join(export, "perturbation.txt")) as fh:
        for line in fh:
            if "=" in line:
                k, v = line.split("=", 1)
                info[k.strip()] = v.strip()
    nsp, nsn = int(info["n_space"]), int(info["n_snap"])
    nx, ny, nc = (int(v) for v in info["grid_nx_ny_ncomp"].split())
    lo_r, lo_z, hi_r, hi_z = (float(v) for v in info["domain_lo_r_lo_z_hi_r_hi_z"].split())
    X = np.fromfile(os.path.join(export, "perturbation.bin"), "<f8").reshape((nsp, nsn), order="F")
    t = np.loadtxt(os.path.join(export, "times.txt"))
    X = X[:nx * ny].reshape(ny, nx, nsn)                 # first component (the level set): [z, r, t]
    dr = (hi_r - lo_r) / nx
    rc = lo_r + (np.arange(nx) + 0.5) * dr
    keep = np.abs(rc - r0) < band * dr if band > 0 else np.ones(nx, bool)
    zc = lo_z + (np.arange(ny) + 0.5) * (hi_z - lo_z) / ny
    sel = np.ones(nsn, bool) if t_max is None else t <= t_max + 1e-9
    idx = np.where(sel)[0]
    idx = idx[int(skip_frac * len(idx)):]
    dt = t[idx[1]] - t[idx[0]]
    out = {}
    for n in range(1, nharm + 1):
        k = 2 * np.pi * n / zlen
        C = np.tensordot(np.exp(-1j * k * zc), X[:, keep][:, :, idx], axes=(0, 0)) * (2.0 / ny)
        C = C * np.sqrt(rc[keep])[:, None]
        A, B = C[:, :-1], C[:, 1:]
        U, sv, Vh = np.linalg.svd(A, full_matrices=False)
        r = max(1, int(np.sum(sv > tol * sv[0]))) if sv[0] > 0 else 0
        if r == 0:
            continue
        U, sv, Vh = U[:, :r], sv[:r], Vh[:r]
        St = U.conj().T @ B @ Vh.conj().T / sv
        mu, W = np.linalg.eig(St)
        Phi = B @ Vh.conj().T / sv @ W
        b = np.abs(np.linalg.lstsq(Phi, C[:, 0], rcond=None)[0])
        om = np.log(mu.astype(complex)) / dt
        ok = (om.real > 0) & (np.abs(om.imag) <= freq_rel * om.real)
        if ok.any():
            j = np.argmax(np.where(ok, b, -1))
            out[n] = (om[j].real, r)
    return out


def interface_rates(records_dir):
    rates = {}
    if not records_dir:
        return rates
    for f in sorted(glob.glob(os.path.join(records_dir, "h*.csv"))):
        n = int(os.path.basename(f)[1:3])
        with open(f) as fh:
            rows = list(csv.DictReader(fh))
        if rows:
            rates[n] = (float(rows[0]["rate"]), float(rows[0]["efolds"]))
    return rates


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--archive", required=True)
    ap.add_argument("--records", default=None)
    ap.add_argument("--zlen", type=float, required=True, help="box length L (= ns.yblob of the broadband deck)")
    ap.add_argument("--nharm", type=int, required=True)
    ap.add_argument("--mu-l", type=float, default=0.02, dest="mu")
    ap.add_argument("--rho-g", type=float, default=1.225e-3, dest="rho_g")
    ap.add_argument("--mu-g", type=float, default=2.6e-4, dest="mu_g")
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--export", default=None, help="level-set twin-difference export for the harmonic-resolved DMD")
    ap.add_argument("--t-max", type=float, default=None, dest="t_max", help="end of the linear range of the sum of harmonics")
    a = ap.parse_args()
    hdmd = harmonic_dmd(a.export, a.nharm, a.zlen, a.t_max) if a.export and os.path.isdir(a.export) else {}

    modes = dmd_modes(a.archive, a.nharm) if os.path.isdir(a.archive) else []
    best = {}
    for n, rate, amp, pur in modes:
        if n not in best or amp > best[n][1]:
            best[n] = (rate, amp, pur)
    iface = interface_rates(a.records)

    rows = []
    print(f"{'n':>3} {'k R0':>7} {'harm.DMD':>10} {'DMD':>10} {'purity':>7} {'interface':>10} {'e-folds':>8} {'Chandrasekhar':>14} {'two-fluid':>10} {'hDMD vs 2f':>10} {'intf vs 2f':>10}")
    for n in range(1, a.nharm + 1):
        x = 2 * np.pi * n / a.zlen
        if x >= 1:
            continue
        ch = float(dispersion.growth_rate_viscous_exact(x, mu=a.mu))
        tf = growth_rate_two_fluid(x, a.mu, a.rho_g, a.mu_g / a.mu) if a.mu > 0 else np.nan
        d = best.get(n, (np.nan, np.nan, np.nan))
        f = iface.get(n, (np.nan, np.nan))
        h = hdmd.get(n, (np.nan, 0))[0]
        rows.append((n, x, d[0], d[2], f[0], f[1], ch, tf, h))
        dd = 100 * (h / tf - 1) if np.isfinite(h) else np.nan
        ff = 100 * (f[0] / tf - 1) if np.isfinite(f[0]) else np.nan
        print(f"{n:3d} {x:7.4f} {h:10.6f} {d[0]:10.6f} {d[2]:7.3f} {f[0]:10.6f} {f[1]:8.2f} {ch:14.6f} {tf:10.6f} {dd:+10.3f} {ff:+10.3f}")

    if a.outdir:
        os.makedirs(a.outdir, exist_ok=True)
        with open(os.path.join(a.outdir, "spectrum.csv"), "w") as fh:
            fh.write("n,kR0,rate_dmd,dmd_purity,rate_interface,efolds_interface,chandrasekhar,two_fluid_viscous,rate_harmonic_dmd\n")
            for r in rows:
                fh.write(",".join(f"{v:.8g}" if isinstance(v, float) else str(v) for v in r) + "\n")
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from pub_style import apply_style, cm_both, savefig_all
            apply_style()
            xs = np.linspace(0.01, 0.999, 300)
            fig, ax = plt.subplots(figsize=(5.4, 3.4))
            ax.plot(xs, dispersion.growth_rate_viscous_exact(xs, mu=a.mu), "k-", lw=1.0, label="Chandrasekhar (viscous)")
            r = np.array(rows, dtype=float)
            ax.plot(r[:, 1], r[:, 7], "k--", lw=0.8, label="two-fluid viscous")
            ax.plot(r[:, 1], r[:, 4], "o", ms=4, mfc="none", label="interface fit")
            ax.plot(r[:, 1], r[:, 8], "s", ms=4, mfc="none", label="harmonic-resolved DMD")
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 0.36)
            ax.set_xlabel(r"$k R_0$")
            ax.set_ylabel(r"growth rate $s\, t_\sigma$")
            ax.legend(frameon=False, loc="lower center")
            cm_both(ax)
            fig.tight_layout()
            savefig_all(fig, "broadband_spectrum", outdir=a.outdir, close=True)
        except Exception as exc:  # the table is the primary output
            print("figure not written:", exc)
        print(f"written to {a.outdir}/")


if __name__ == "__main__":
    main()
