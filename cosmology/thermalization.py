#!/usr/bin/env python3
"""Thermalization of the dark sector: Fig. 1 of arXiv:2506.04718.

Plots, against the SM temperature T, the inverse-decay rate Gamma_ID
(Eqs. 20 and 22), the 2 -> 2 scattering rate Gamma_I (Eq. 24) and the
Hubble rate H(T), for the four parameter sets of the figure, and prints
the temperatures T0 and T1 between which Gamma_ID > H (Eq. 25).

Two normalizations of the O(1) factors that the paper writes as "~":

  paper    the factors read off the published figure (default). With
           g_Y = 0.357 the paper's curves are, to 0.1%,
             Gamma_ID = 0.4636 Gamma_0 M/T                   (T > M)
             Gamma_ID = 1.076  Gamma_0 (M/T)^{3/2} e^{-M/T}  (T < M)
             Gamma_I  = 0.9657 a_X a_Y eps^2 Q_X^2 s T^3 / (12 pi |D(s)|^2),
                        s = 4 <p>^2, <p> = 3.15 T
  derived  the coefficients derived in the appendices of the dark_radiation
           draft: Gamma_ID from the full thermal average of Eq. (20),
           Bose-Einstein for Z', valid for any T, and
           Gamma_I = (4 pi/3) a_X a_Y eps^2 Q_X^2 sbar/|D(sbar)|^2 n_f,eq,
           sbar = 2 <p>^2.

Gamma_0 = a_Y Q_Y^2 eps^2 M/6 (Eq. 21) with Q_Y = 1, H uses g_* = 106.75
and the reduced Planck mass. The Z' width in |D(s)|^2 is that of its
decays into the massless dark Weyl fermions (charges 7 and 8); the height
of the resonance spike depends on the temperature grid.

Usage:  python3 thermalization.py [--normalization paper|derived] [--outdir DIR]
"""
import argparse
from pathlib import Path

import numpy as np
from scipy import integrate
from scipy.special import zeta
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

MP = 2.435e18            # reduced Planck mass [GeV]
GSTAR = 106.75
GY = 0.357               # hypercharge coupling
ALPHA_Y = GY**2 / (4 * np.pi)
QY = 1.0
QX = 8                   # reference dark charge in Gamma_I (paper)
QX_MASSLESS = (7, 8)     # massless dark Weyl fermions, for the Z' width
P_MEAN = 7 * np.pi**4 / (180 * zeta(3))   # <p>/T for fermions, 3.15
N_F = 3 * zeta(3) / (2 * np.pi**2)        # n_f,eq/T^3 with g_f = 2, 0.18

PANELS = [  # (M_Z' [GeV], g_X, eps), in the order of Fig. 1
    (30.0, 1e-2, 1e-3),
    (30.0, 1e-3, 1e-3),
    (5.0, 1e-3, 1e-3),
    (10.0, 1e-2, 1e-6),
]


def hubble(T):
    return np.sqrt(np.pi**2 * GSTAR / 90) * T**2 / MP


def gamma0(M, eps):
    """Rest-frame width Z' -> f fbar, Eq. (21)."""
    return ALPHA_Y * QY**2 * eps**2 * M / 6


def width_zp(M, gX):
    """Z' width into the massless dark Weyl fermions."""
    return sum((gX * q)**2 for q in QX_MASSLESS) * M / (24 * np.pi)


def gamma_id_paper(T, M, eps):
    x = M / T
    return np.where(T > M, 0.4636 * gamma0(M, eps) * x,
                    1.076 * gamma0(M, eps) * x**1.5 * np.exp(-np.minimum(x, 700)))


def gamma_id_derived(T, M, eps):
    """Eq. (20), second line, with Bose-Einstein f(k) and g_Z' = 3."""
    out = []
    for t in np.atleast_1d(T):
        if M / t > 200:
            out.append(0.0)
            continue

        def f(k):
            return 1 / np.expm1(np.sqrt(k * k + M * M) / t)
        kmax = M + 60 * t
        num = integrate.quad(lambda k: k * k * f(k) * M / np.sqrt(k * k + M * M), 0, kmax, limit=200)[0]
        den = integrate.quad(lambda k: k * k * f(k), 0, kmax, limit=200)[0]
        n_zp = 3 * den / (2 * np.pi**2)
        out.append(gamma0(M, eps) * num / den * n_zp / (N_F * t**3))
    return np.array(out)


def gamma_i(T, M, gX, eps, normalization):
    aX = gX**2 / (4 * np.pi)
    G = width_zp(M, gX)
    base = aX * ALPHA_Y * eps**2 * QX**2
    if normalization == "paper":
        s = 4 * (3.15 * T)**2
        return 0.9657 * base * s * T**3 / (12 * np.pi * ((s - M**2)**2 + (M * G)**2))
    s = 2 * (P_MEAN * T)**2
    return 4 * np.pi / 3 * QY**2 * base * s / ((s - M**2)**2 + (M * G)**2) * N_F * T**3


def window(T, rate, H):
    """(T1, T0): rate > H for T1 < T < T0 on the grid, by log interpolation.
    None for an edge of the grid that is still inside the window."""
    d = np.log(np.maximum(rate, 1e-300)) - np.log(H)
    above = d > 0
    if not above.any():
        return None

    def cross(i):
        return 10**np.interp(0, sorted([d[i], d[i + 1]]),
                             np.log10(T[i:i + 2]) if d[i] < d[i + 1] else np.log10(T[i:i + 2])[::-1])
    i0, i1 = np.argmax(above), len(above) - 1 - np.argmax(above[::-1])
    t1 = cross(i0 - 1) if i0 > 0 else None
    t0 = cross(i1) if i1 < len(T) - 1 else None
    return t1, t0


def panel(ax, M, gX, eps, normalization, T):
    gid = gamma_id_paper(T, M, eps) if normalization == "paper" else gamma_id_derived(T, M, eps)
    gi = gamma_i(T, M, gX, eps, normalization)
    H = hubble(T)
    ax.loglog(T, gi, "k-", lw=1.2, label=r"$\Gamma_I$")
    ax.loglog(T, gid, "r:", lw=1.5, label=r"$\Gamma_{ID}$")
    ax.loglog(T, H, "b--", lw=1.5, label="Hubble rate")
    ax.set_xlim(0.1, 1e3)
    ax.set_ylim(1e-22, 1e-7)
    ax.set_xlabel(r"$T$ (GeV)")
    ax.set_ylabel(r"(GeV)")
    g = f"{gX:g}" if gX >= 0.01 else rf"10^{{{int(round(np.log10(gX)))}}}"
    ax.set_title(rf"$M_{{Z'}}={M:g}\,\mathrm{{GeV}},\ g_X={g},\ \epsilon=10^{{{int(round(np.log10(eps)))}}}$",
                 fontsize=10)
    ax.legend(fontsize=8, loc="best")
    return window(T, gid, H)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--normalization", choices=["paper", "derived"], default="paper")
    p.add_argument("--outdir", default="figures")
    p.add_argument("--points", type=int, default=1000, help="temperatures on a log grid from 0.1 to 1000 GeV")
    a = p.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    T = np.logspace(-1, 3, a.points)
    suffix = "" if a.normalization == "paper" else "_derived"

    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
    print(f"normalization: {a.normalization}")
    print(f"{'M_Zp':>6} {'g_X':>7} {'eps':>7}   T1 [GeV]   T0 [GeV]   (Gamma_ID > H for T1 < T < T0)")
    for n, (M, gX, eps) in enumerate(PANELS, 1):
        ax = axes.flat[n - 1]
        cr = panel(ax, M, gX, eps, a.normalization, T)
        f1, a1 = plt.subplots(figsize=(4.6, 3.5))
        panel(a1, M, gX, eps, a.normalization, T)
        f1.tight_layout()
        f1.savefig(out / f"thermalize{n}{suffix}.pdf")
        plt.close(f1)
        if cr is None:
            print(f"{M:6g} {gX:7g} {eps:7g}   never thermalized by inverse decay")
            continue
        t1, t0 = (f"{t:8.3g}" if t is not None else (">1000" if k else "<0.1").rjust(8) for k, t in enumerate(cr))
        print(f"{M:6g} {gX:7g} {eps:7g}   {t1}   {t0}")
    fig.tight_layout()
    fig.savefig(out / f"fig1{suffix}.pdf")
    print(f"figures written to {out}/")


if __name__ == "__main__":
    main()
