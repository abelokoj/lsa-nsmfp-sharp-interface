"""Superposition test of Paper 1 (campaign_p1_superposition.txt).

For each stored time, compare the twin difference of the summed shape, D_C, with the sum of the twin differences of the two single shapes,
D_A + D_B, translated along the axis by delta = g L / (2 pi), the common shift of both harmonics of the broadband perturbation.
Reports the relative error ||D_C - S(D_A + D_B)|| / ||D_C|| for the velocity (whole domain) and for the level set in a band around r = R0,
against the interface amplitude of C (max |phi'| in the band, which equals the displacement there).

Usage (from $LSA, or with the export folders copied to the workstation):
  python3 superposition_check.py runs/analysis/sp035_16 runs/analysis/sp070_16 runs/analysis/spsum_16 [out.csv]
Each folder must hold export/ and export_ls/ as written by the post job (analyse_campaign.sh).
"""
import sys
import numpy as np

G = 2.399963229728653  # golden angle, the phase step of the broadband perturbation (solver_changes/broadband)


def load(d):
    meta = {}
    for line in open(f"{d}/perturbation.txt"):
        if "=" in line:
            k, v = line.split("=", 1)
            meta[k.strip()] = v.strip()
    nsp, nsn = int(meta["n_space"]), int(meta["n_snap"])
    nx, ny, nc = (int(x) for x in meta["grid_nx_ny_ncomp"].split())
    lo_r, lo_z, hi_r, hi_z = (float(x) for x in meta["domain_lo_r_lo_z_hi_r_hi_z"].split())
    X = np.fromfile(f"{d}/perturbation.bin", "<f8").reshape((nsp, nsn), order="F")
    t = np.loadtxt(f"{d}/times.txt")
    return X, t, (nx, ny, nc), (lo_r, hi_r, lo_z, hi_z)


def fields(X, grid):
    nx, ny, nc = grid
    return X.reshape((nc, ny, nx, -1))  # column k of X is component-major, z, then r fastest (C order of (nc, ny, nx))


def shift_z(F, delta, length):
    """F(z + delta) on a periodic axis (axis 1 of (nc, ny, nx, nt)), spectrally exact."""
    ny = F.shape[1]
    kz = 2 * np.pi * np.fft.fftfreq(ny, d=length / ny)
    return np.real(np.fft.ifft(np.fft.fft(F, axis=1) * np.exp(1j * kz * delta)[None, :, None, None], axis=1))


def self_check():
    ny, L = 64, 3.0
    z = (np.arange(ny) + 0.5) * L / ny
    f = np.sin(2 * np.pi * z / L)[None, :, None, None]
    d = 0.37
    assert np.allclose(shift_z(f, d, L)[0, :, 0, 0], np.sin(2 * np.pi * (z + d) / L), atol=1e-12)


def main(a, b, c, out):
    self_check()
    res = {}
    for kind, sub in (("vel", "export"), ("ls", "export_ls")):
        XA, tA, gA, dom = load(f"{a}/{sub}")
        XB, tB, gB, domB = load(f"{b}/{sub}")
        XC, tC, gC, domC = load(f"{c}/{sub}")
        assert gA == gB == gC and np.allclose(dom, domB) and np.allclose(dom, domC), "grids or domains differ"
        n = min(len(tA), len(tB), len(tC))
        assert np.allclose(tA[:n], tB[:n]) and np.allclose(tA[:n], tC[:n]), "snapshot times differ"
        FA, FB, FC = (fields(X[:, :n], gA) for X in (XA, XB, XC))
        lo_r, hi_r, lo_z, hi_z = dom
        L = hi_z - lo_z
        S = shift_z(FA + FB, G * L / (2 * np.pi), L)
        if kind == "ls":
            nx = gA[0]
            dr = (hi_r - lo_r) / nx
            rc = lo_r + (np.arange(nx) + 0.5) * dr
            band = np.abs(rc - 1.0) <= 2.0 * dr  # two cells either side of the unperturbed interface r = R0
            FC, S = FC[:, :, band, :], S[:, :, band, :]
            res["amp"] = np.abs(FC).max(axis=(0, 1, 2))
            # C must carry harmonics 1 and 2 of the box at equal amplitude; a solver built without the broadband change gives harmonic 1 only
            c = np.abs(np.fft.rfft(FC[0, :, :, 0].mean(axis=1)))
            print(f"initial shape of C: |harmonic 2| / |harmonic 1| = {c[2] / c[1]:.3f} (expected 1; about 0 means the broadband build was not used)")
        num = np.sqrt(((FC - S) ** 2).sum(axis=(0, 1, 2)))
        den = np.sqrt((FC ** 2).sum(axis=(0, 1, 2)))
        res[kind] = np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)
        res["t"] = tA[:n]
    ok = res["amp"] <= 0.1  # the admissible window of Paper 1 closes at a = 0.1 R0
    with open(out, "w") as f:
        f.write("t,amp_C,relerr_vel,relerr_ls_band\n")
        for row in zip(res["t"], res["amp"], res["vel"], res["ls"]):
            f.write("%.6e,%.6e,%.6e,%.6e\n" % row)
    for kind in ("vel", "ls"):
        e = res[kind][ok & np.isfinite(res[kind])]
        print(f"{kind}: max relative error {100*e.max():.3f}% over {e.size} snapshots with amp_C <= 0.1 R0 (median {100*np.median(e):.3f}%)")
    # One power law over the whole window hides a U shape: a fixed error floor dominates at small amplitude (exponent near -1) and the
    # interaction of the harmonics at large amplitude (exponent near +1 for a smooth map), so each side of the minimum is fitted separately.
    for kind in ("vel", "ls"):
        e, a = res[kind], res["amp"]
        m = ok & np.isfinite(e) & (a > 0) & (e > 0)
        if m.sum() < 3:
            continue
        i = np.flatnonzero(m)[np.argmin(e[m])]
        line = f"{kind}: minimum {100*e[i]:.3f}% at amp_C = {a[i]:.4f} R0 (t = {res['t'][i]:.3f})"
        for side, sel in (("below", m & (a <= a[i])), ("above", m & (a >= a[i]))):
            if sel.sum() >= 3:
                line += f"; {side} it error ~ amp^{np.polyfit(np.log(a[sel]), np.log(e[sel]), 1)[0]:.2f}"
        print(line)
    print(f"written {out}")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    main(*sys.argv[1:4], sys.argv[4] if len(sys.argv) > 4 else "superposition.csv")
