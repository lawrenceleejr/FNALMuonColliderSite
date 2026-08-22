#!/usr/bin/env python3
"""Figures for the co-tilted chain variant (tools/uiuc_chain.py).

Standalone: rebuilds terrain and beam lines from the solved parameters in
static/geo/uiuc_chain.json, and the far-hall payoff series from the .npz.
"""
import json
import math
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
STYLE = "/root/.claude/skills/synced/tufte-plots/assets/tufte.mplstyle"
if os.path.exists(STYLE):
    plt.style.use(STYLE)
plt.rcParams["figure.dpi"] = 130

import corridor_layout as cl      # terrain north of 41.79 (fast import)

R_E = 6371000.0
IP_LAT = 41.8443
LAT = lambda y: IP_LAT + y / cl.MLAT
Y = lambda lat: (lat - IP_LAT) * cl.MLAT

_tool = open(os.path.join(ROOT, "static", "tool", "index.html")).read()
SOUTH_PROF = json.loads(re.search(r"const PROFILE_M = (\[\[.*?\]\]);", _tool).group(1))

def terrain(lat):
    if lat >= 41.79:
        return cl.elev_at_lat(lat)
    for (la1, e1), (la2, e2) in zip(SOUTH_PROF, SOUTH_PROF[1:]):
        if la1 <= lat <= la2:
            return e1 + (lat - la1) / (la2 - la1) * (e2 - e1)
    return SOUTH_PROF[0][1]

def beam_z(z0, th, y0, y):
    d = y - y0
    return z0 + th * d + d * d / (2 * R_E)

D = json.load(open(os.path.join(ROOT, "static", "geo", "uiuc_chain.json")))
npz = np.load(os.path.join(ROOT, "static", "geo", "uiuc_chain.npz"))
FIGS = os.path.join(ROOT, "static", "figs")

def save(fig, name):
    for ext in ("pdf", "svg"):
        fig.savefig(os.path.join(FIGS, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)

# ----------------------------------------------------------------------
# Figure A: the north side -- where every beam emerges, with terrain
# ----------------------------------------------------------------------
BEAMS = [
    ("collider", 0.0, D["collider_baseline"], "#0b6e4f"),
    ("RCS1/2", 0.0, D["rings"]["RCS1/2"], "#6d3580"),
    ("RCS3/4", -390.0, D["rings"]["RCS3/4"], "#b5541c"),
]
ys = np.linspace(-600, 4600, 800)
fig, ax = plt.subplots(figsize=(6.8, 3.9))
terr = np.array([terrain(LAT(y)) for y in ys])
ax.fill_between(ys / 1e3, terr, 150, color="0.90", lw=0)
ax.plot(ys / 1e3, terr, color="0.45", lw=1.2)
LABEL_X = {"collider": 1.35, "RCS1/2": 2.05, "RCS3/4": 3.15}
for name, y0, s, col in BEAMS:
    z0, th = s["z0_masl"], s["theta_mrad"] * 1e-3
    zb = beam_z(z0, th, y0, ys)
    ax.plot(ys / 1e3, np.where(zb < terr + 60, zb, np.nan), color=col, lw=1.8)
    y_em = Y(s["north_emergence_lat"])
    ax.plot(y_em / 1e3, terrain(s["north_emergence_lat"]), "o", color=col,
            ms=5, zorder=5)
    ax.annotate(name + "\n" + "%.4f" % s["north_emergence_lat"] +
                "$^{\\circ}$N\n" + "%.0f m inside" % s["inside_fence_m"],
                (y_em / 1e3, terrain(s["north_emergence_lat"])),
                xytext=(LABEL_X[name], 262), fontsize=7.5, color=col,
                ha="center", va="bottom",
                arrowprops=dict(arrowstyle="-", color=col, lw=0.7,
                                shrinkB=3))
y_f = Y(41.8699) / 1e3
ax.axvline(y_f, color="0.2", lw=1.0, ls=(0, (4, 3)))
ax.text(y_f + 0.06, 163, "Fermilab fence\n41.8699$^{\\circ}$N", fontsize=7.5, color="0.2")
ax.plot(1.3, terrain(41.8560) - 15, "s", color="0.25", ms=5)
ax.text(1.3, terrain(41.8560) - 34, "near hall\n(15 m deep)", fontsize=7.5,
        color="0.25", ha="center")
ax.text(-0.28, 160, "IP straights\n(side-by-side, 33$-$45 m deep)", fontsize=7.5,
        color="0.3", ha="left")
ax.text(3.7, 236, "terrain", fontsize=8, color="0.45")
ax.set_xlim(-0.6, 4.6)
ax.set_ylim(152, 292)
ax.set_xlabel("km north of the IP")
ax.set_ylabel("elevation  [m ASL]")
ax.set_title("Co-tilted chain: every north beam emerges on site",
             fontsize=10, loc="left")
save(fig, "uiuc_chain_profile")

# ----------------------------------------------------------------------
# Figure B: nutau arrivals vs time at the UIUC far hall, baseline vs variant
# ----------------------------------------------------------------------
T = npz["T_EDGES"]
T_CTR = np.sqrt(T[:-1] * T[1:])
fig, ax = plt.subplots(figsize=(6.6, 3.8))
for key, col, lab, dy, ls in (
        ("VAR_NUTAU", "#b5541c", "co-tilted chain (chirp + store)", 1.9, "-"),
        ("BASE_NUTAU", "0.4", "baseline (store only)", 0.42, (0, (3, 2)))):
    y = npz[key]
    ax.plot(T_CTR * 1e3, np.ma.masked_where(~(y > 0), y),
            drawstyle="steps-mid", color=col, lw=1.5, ls=ls)
    i = np.where(y > 0)[0][-1]
    ax.text(210, y[i] * dy, lab, color=col, fontsize=8.5, va="center")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(2e-2, 320)
ax.set_ylim(4e-6, 0.3)
ax.set_xlabel("time in the 5 Hz cycle  [ms]  (bunch-centric)")
ax.set_ylabel(r"oscillated $\bar\nu_\tau$ / cm$^2$ / ms (per cycle)")
ax.text(0.5, 3e-2, "chirp at 198 km:\n(L/E)$^2$ reaches P$\\sim$10$^{-3}$$-$10$^{-2}$\n"
        "in the soft turns", fontsize=7.8, color="#b5541c", ha="center")
ax.text(0.98, 0.96, "per tonne-year on axis: 0.14 $\\to$ 0.29 $\\nu_\\tau$ CC "
        "($\\tau$-threshold included);\nwhole-plane $\\bar\\nu_\\tau$ (E $>$ 5 GeV): "
        "1.0 $\\times$ 10$^{14}$ $\\to$ 9.7 $\\times$ 10$^{14}$ /yr;\n"
        "store portions coincide by construction",
        transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color="0.3")
ax.set_title(r"UIUC far hall: what aiming the chain south buys in $\bar\nu_\tau$",
             fontsize=10, loc="left")
save(fig, "uiuc_chain_nutau")
