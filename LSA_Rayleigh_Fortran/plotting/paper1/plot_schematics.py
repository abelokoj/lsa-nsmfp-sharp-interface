"""
plot_schematics.py
==================
The two schematic figures of Paper 1, in the house style of pub_style.py (PNG, PDF and SVG):

    figures/nsmfp_variants.*         the continuous-twin and restarting variants of NS-MFP (from the thesis,
                                     figures/make_ch4_figures.py, Fig. 3)
    figures/ch3_column_schematic.*   the liquid column in its periodic domain with the no-slip wall (from the thesis
                                     chapter 3 schematic)

Both are drawings without axes, so the tick and grid rules of the style do not apply; fonts, sizes and colours do.

    python plot_schematics.py
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pub_style import apply_style, savefig_all, subcaption, PALETTE

apply_style()
FIGS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'figures')
BASE, PERT, DIFF, GREY = 'k', PALETTE[0], PALETTE[1], '0.35'


def variants():
    fig, axs = plt.subplots(2, 1, figsize=(6.4, 3.4), sharex=True)
    t = np.linspace(0, 6, 400)
    ax = axs[0]
    ax.plot(t, 0 * t, color=BASE, lw=1.2, label=r'twin run $S_b(t)$')
    ax.plot(t, 0.15 * np.exp(0.42 * t), color=PERT, lw=1.5, label=r'perturbed run $S_p(t)$')
    for k in range(7):
        ax.plot([k, k], [0, 0.15 * np.exp(0.42 * k)], color=DIFF, lw=1.0)
    ax.text(6.05, 0.15 * np.exp(0.42 * 6) * 0.9, r'$\delta S^{(k)} = S_p - S_b$', color=DIFF)
    subcaption(ax, r'(a) continuous twin: two runs, differences sampled every $\tau$', loc='left')
    ax.legend(loc='upper left', frameon=False)
    ax.set_xlim(-0.05, 6.05); ax.axis('off')
    ax = axs[1]
    ax.plot(t, 0 * t, color=BASE, lw=1.2)
    # the solver's driver (commit dddb827) rescales only the first member, so each cycle starts at the amplitude the previous one reached
    for k in range(6):
        a0 = 0.15 * np.exp(0.42 * k); tt = np.linspace(k, k + 1, 50)
        ax.plot(tt, a0 * np.exp(0.42 * (tt - k)), color=PERT, lw=1.5)
        ax.plot([k, k], [0, a0], ':', color=BASE, lw=0.8)
        ax.plot([k + 1, k + 1], [0, a0 * np.exp(0.42)], color=DIFF, lw=1.0)
    subcaption(ax, r'(b) restarting: each member restarts from $\bar S$ plus the previous difference (only the first is rescaled)',
               loc='left')
    ax.text(5.9, -0.35, r'$t$')
    ax.set_xlim(-0.05, 6.05); ax.axis('off')
    fig.tight_layout()
    savefig_all(fig, 'nsmfp_variants', outdir=FIGS, close=True)


def column():
    fig, ax = plt.subplots(figsize=(6.4, 2.9))
    z = np.linspace(0, 2 * np.pi, 400)
    R0, eps, Rw = 1.0, 0.18, 1.9
    r = R0 + eps * np.cos(z)
    ax.fill_between(z, -r, r, color=PERT, alpha=0.15, lw=0)
    ax.plot(z, r, color=PERT, lw=1.5); ax.plot(z, -r, color=PERT, lw=1.5)
    ax.plot(z, 0 * z, color=GREY, lw=0.8, ls='-.')
    for s in (+1, -1):
        ax.plot(z, s * Rw + 0 * z, color='k', lw=1.4)
        ax.fill_between(z, s * Rw, s * (Rw + 0.12), color='k', alpha=0.15, hatch='////', lw=0)
    arrow = dict(arrowstyle='<->', color='k', lw=0.8)
    ax.annotate('', xy=(np.pi, R0 - eps), xytext=(np.pi, 0), arrowprops=arrow)
    ax.text(np.pi + 0.08, 0.35, r'$R_0 + \eta$')
    ax.annotate('', xy=(0.5, Rw), xytext=(0.5, 0), arrowprops=arrow)
    ax.text(0.6, 1.25, r'$R_w$')
    ax.annotate('', xy=(2 * np.pi, -1.68), xytext=(0, -1.68), arrowprops=dict(arrow, color=GREY))
    ax.text(np.pi, -1.62, r'$\lambda = 2\pi/k$ (periodic)', ha='center', va='bottom', color=GREY,
            bbox=dict(fc='white', ec='none', pad=1))
    ax.text(4.4, -0.45, r'liquid $\rho_1,\ \mu_1$', color=PERT, va='center')
    ax.text(4.5, 1.45, r'gas $\rho_2,\ \mu_2$', color=GREY)
    ax.text(1.8, 2.12, 'no-slip wall')
    ax.text(-0.15, 0.0, 'axis', color=GREY, ha='right', va='center')
    ax.set_xlim(-0.6, 2 * np.pi + 0.1); ax.set_ylim(-2.15, 2.4)
    ax.axis('off')
    fig.tight_layout()
    savefig_all(fig, 'ch3_column_schematic', outdir=FIGS, close=True)


if __name__ == '__main__':
    variants()
    column()
