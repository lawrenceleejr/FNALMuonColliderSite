#!/usr/bin/env python3
"""Draw the chain-timing figures from static/geo/chain_maps.npz.

One panel per file; PDF for the paper, SVG (mesh rasterized) for the site.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
STYLE = "/root/.claude/skills/synced/tufte-plots/assets/tufte.mplstyle"
if os.path.exists(STYLE):
    plt.style.use(STYLE)
plt.rcParams["figure.dpi"] = 130

d = np.load(os.path.join(ROOT, "static", "geo", "chain_maps.npz"))
T, E, R = d["T_EDGES"], d["E_EDGES"], d["R_EDGES"]
MAP_TE, MAP_ER, MAP_ER_NEAR = d["MAP_TE"], d["MAP_ER"], d["MAP_ER_NEAR"]

FIGS = os.path.join(ROOT, "static", "figs")
os.makedirs(FIGS, exist_ok=True)

def save(fig, name):
    for ext in ("pdf", "svg"):
        fig.savefig(os.path.join(FIGS, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)

CMAP = "magma"

def logdens(m, xedges, yedges):
    """Per-bin counts -> density per (decade x . decade y)."""
    dx = np.diff(np.log10(xedges))
    dy = np.diff(np.log10(yedges))
    return m / (dy[:, None] * dx[None, :])

# ----------------------------------------------------------------------
# Figure 1: (time in cycle, E_nu) at the deep hall
# ----------------------------------------------------------------------
z = logdens(MAP_TE, T, E)
fig, ax = plt.subplots(figsize=(6.6, 4.0))
vmax = z.max()
pc = ax.pcolormesh(T * 1e3, E, z, norm=LogNorm(vmin=vmax / 1e5, vmax=vmax),
                   cmap=CMAP, rasterized=True)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(2e-2, 200)
ax.set_ylim(2, 5200)
ax.set_xlabel("time in the 5 Hz cycle  [ms]")
ax.set_ylabel(r"$E_\nu$  [GeV]")
cb = fig.colorbar(pc, ax=ax, pad=0.015)
cb.set_label(r"$\nu$ / cm$^2$ / cycle per (decade $E$ $\cdot$ decade $t$)")
cb.outline.set_visible(False)
# direct stage labels: dark ink on the white above each chirp segment
lab = dict(color="0.25", ha="center", fontsize=8.5)
ax.text(0.075, 700, "RCS1\n63$\\to$314 GeV", **lab)
ax.text(0.52, 1650, "RCS2\n$\\to$750 GeV", **lab)
ax.text(1.45, 3200, "RCS3\n$\\to$1.5 TeV", **lab)
ax.text(4.2, 6.5, "RCS4 $\\to$5 TeV", color="w", ha="center", fontsize=8.5)
ax.text(42, 2500, "collider store, 5 TeV", color="0.15", ha="center",
        fontsize=9)
ax.text(42, 1250, "fades as $e^{-t/104\\,\\mathrm{ms}}$", color="0.15",
        ha="center", fontsize=8)
ax.text(0.98, 0.03,
        "each column is a train of ~ps pulses:\n"
        "21.0 / 49.1 / 36.7 $\\mu$s spacing (RCS1-2 / RCS3-4 / collider);\n"
        "single turns resolved left of ~0.5 ms",
        transform=ax.transAxes, ha="right", va="bottom",
        fontsize=7.5, color="0.35")
ax.set_title("One machine cycle at the deep hall (9.25 km): the chain chirp,"
             " then the store", fontsize=10, loc="left")
save(fig, "chain_timing_te")

# ----------------------------------------------------------------------
# Figure 2: (r_perp, E_nu) at the deep hall (all plumes converged)
# ----------------------------------------------------------------------
z = logdens(MAP_ER, R, E)
fig, ax = plt.subplots(figsize=(6.6, 4.0))
vmax = z.max()
pc = ax.pcolormesh(R, E, z, norm=LogNorm(vmin=vmax / 3e6, vmax=vmax),
                   cmap=CMAP, rasterized=True)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(0.02, 80)
ax.set_ylim(2, 5200)
ax.set_xlabel(r"$r_\perp$ from the convergence axis  [m]")
ax.set_ylabel(r"$E_\nu$  [GeV]")
cb = fig.colorbar(pc, ax=ax, pad=0.015)
cb.set_label(r"$\nu$ / cm$^2$ / yr per decade of $E_\nu$")
cb.outline.set_visible(False)
ax.text(0.085, 2500, "RCS4 pencil chirp\n$r\\sim L/\\gamma$", color="0.15",
        ha="center", fontsize=8.5)
ax.text(1.6, 1300, "collider store (smeared,\n$\\sigma_\\theta$ = 0.15 mrad: no prism)",
        color="0.15", ha="center", fontsize=8.5)
ax.text(3.6, 45, "RCS2$-$3", color="w", ha="center", fontsize=8.5)
ax.text(16, 8.5, "RCS1", color="w", ha="center", fontsize=8.5)
ax.set_title("Fluence at the deep hall: every stage converged on one axis",
             fontsize=10, loc="left")
save(fig, "chain_flux_er")

# ----------------------------------------------------------------------
# Figure 3: (r_perp, E_nu) at the civic-envelope near hall (1.3 km)
# ----------------------------------------------------------------------
z = logdens(MAP_ER_NEAR, R, E)
fig, ax = plt.subplots(figsize=(6.6, 4.0))
vmax = z.max()
pc = ax.pcolormesh(R, E, z, norm=LogNorm(vmin=vmax / 3e7, vmax=vmax),
                   cmap=CMAP, rasterized=True)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(0.02, 80)
ax.set_ylim(2, 5200)
ax.set_xlabel(r"$r_\perp$ from the collider axis  [m]")
ax.set_ylabel(r"$E_\nu$  [GeV]")
cb = fig.colorbar(pc, ax=ax, pad=0.015)
cb.set_label(r"$\nu$ / cm$^2$ / yr per decade of $E_\nu$")
cb.outline.set_visible(False)
ax.text(0.15, 1100, "collider store\n($\\sigma_\\theta$ = 0.15 mrad)",
        color="0.15", ha="center", fontsize=8.5)
ax.text(11, 900, "RCS pencils cross the plane\n50$-$68 m below the hall\n(ring-averaged here)",
        color="0.25", ha="center", fontsize=8)
ax.set_title("Fluence at the near hall (1.3 km, civic envelope): "
             "collider only", fontsize=10, loc="left")
save(fig, "chain_flux_er_near")
