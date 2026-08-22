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
fv = np.load(os.path.join(ROOT, "static", "geo", "chain_flavor.npz"))
E_CTR = np.sqrt(E[:-1] * E[1:])
T_CTR = np.sqrt(T[:-1] * T[1:])

def mean_profile(m):
    tot = m.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(tot > 0, (m * E_CTR[:, None]).sum(axis=0) / tot, np.nan)

import matplotlib.patheffects as pe
MEAN_LINE = dict(color="w", lw=1.8,
                 path_effects=[pe.withStroke(linewidth=3.2, foreground="#2a2a2a")])

FIGS = os.path.join(ROOT, "static", "figs")
os.makedirs(FIGS, exist_ok=True)

def save(fig, name):
    for ext in ("pdf", "svg"):
        fig.savefig(os.path.join(FIGS, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)

CMAP = "magma"
CMAP_OBJ = plt.get_cmap("magma").copy()
CMAP_OBJ.set_bad(CMAP_OBJ(0.0))     # exactly-empty bins print as floor, not white

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
                   cmap=CMAP_OBJ, rasterized=True)
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
prof = mean_profile(MAP_TE)
ax.plot(T_CTR * 1e3, prof, **MEAN_LINE)
ax.text(0.55, 62, "fluence-weighted $\\langle E_\\nu\\rangle(t)$",
        color="w", fontsize=8.5, rotation=27,
        path_effects=[pe.withStroke(linewidth=2.6, foreground="#2a2a2a")])
ax.set_title("One machine cycle at the on-site convergence hall (9.25 km): chirp,"
             " then store", fontsize=10, loc="left")
save(fig, "chain_timing_te")

# ----------------------------------------------------------------------
# Figure 1b: three wall-clock cycles, linear time, with the physical
# interleave (each store runs the full 0.2 s; the next bunch's chirp
# overlaps its last 7 ms), and the mean-energy profile on top.
# ----------------------------------------------------------------------
TM, MAPM, MEANM = fv["T_MULTI"], fv["MAP_MULTI"], fv["MEAN_MULTI"]
zm = MAPM / (np.diff(np.log10(E))[:, None] * (np.diff(TM) * 1e3)[None, :])
fig, ax = plt.subplots(figsize=(6.6, 4.0))
vmax = zm.max()
pc = ax.pcolormesh(TM * 1e3, E, zm, norm=LogNorm(vmin=vmax / 3e4, vmax=vmax),
                   cmap=CMAP_OBJ, rasterized=True)
ax.set_yscale("log")
ax.set_xlim(0, 600)
ax.set_ylim(5, 5200)
ax.set_xlabel("wall-clock time  [ms]")
ax.set_ylabel(r"$E_\nu$  [GeV]")
cb = fig.colorbar(pc, ax=ax, pad=0.015)
cb.set_label(r"$\nu$ / cm$^2$ per (decade $E$ $\cdot$ ms)")
cb.outline.set_visible(False)
tm_ctr = 0.5 * (TM[:-1] + TM[1:]) * 1e3
ax.plot(tm_ctr, MEANM, **MEAN_LINE)
ax.text(52, 1150, r"fluence-weighted $\langle E_\nu\rangle(t)$", color="0.15",
        fontsize=8.5)
for k in range(4):
    ax.axvline(200 * k, color="0.35", lw=0.8, ls=(0, (2, 3)))
ax.text(8, 3600, "inject / dump every 200 ms", fontsize=7.5, color="w",
        path_effects=[pe.withStroke(linewidth=2.2, foreground="#2a2a2a")])
ax.annotate("chirp of the next bunch:\n7 ms, 63 GeV$\\to$5 TeV,\n"
            "pulls $\\langle E\\rangle$ down through\nthe store's tail",
            (196, 620), xytext=(228, 60), fontsize=7.5, color="w",
            path_effects=[pe.withStroke(linewidth=2.2, foreground="#2a2a2a")],
            arrowprops=dict(arrowstyle="->", color="w", lw=0.9))
ax.text(100, 22, "store fades $e^{-t/104\\,\\mathrm{ms}}$; 85% decayed at dump",
        fontsize=7.5, color="w", ha="center",
        path_effects=[pe.withStroke(linewidth=2.2, foreground="#2a2a2a")])
ax.set_title("Three 5 Hz cycles, wall-clock: sawtooth store, chirp in the "
             "overlap window", fontsize=10, loc="left")
save(fig, "chain_timing_multirep")

# ----------------------------------------------------------------------
# Figure 2: (r_perp, E_nu) at the deep hall (all plumes converged)
# ----------------------------------------------------------------------
z = logdens(MAP_ER, R, E)
fig, ax = plt.subplots(figsize=(6.6, 4.0))
vmax = z.max()
pc = ax.pcolormesh(R, E, z, norm=LogNorm(vmin=vmax / 3e6, vmax=vmax),
                   cmap=CMAP_OBJ, rasterized=True)
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
ax.set_title("Fluence at the convergence hall (9.25 km): every stage on one axis",
             fontsize=10, loc="left")
save(fig, "chain_flux_er")

# ----------------------------------------------------------------------
# Figure 3: (r_perp, E_nu) at the civic-envelope near hall (1.3 km)
# ----------------------------------------------------------------------
z = logdens(MAP_ER_NEAR, R, E)
fig, ax = plt.subplots(figsize=(6.6, 4.0))
vmax = z.max()
pc = ax.pcolormesh(R, E, z, norm=LogNorm(vmin=vmax / 3e7, vmax=vmax),
                   cmap=CMAP_OBJ, rasterized=True)
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

# ----------------------------------------------------------------------
# Flavor vs time: per-species interaction rate through one (bunch-centric)
# cycle at each hall, and the oscillation flavor evolution.
# ----------------------------------------------------------------------
C_NUMU, C_NUBE, C_NC = "#0b6e4f", "#b5541c", "0.45"

def mstep(ax, y, col, lw=1.6):
    ym = np.ma.masked_where(~(y > 0), y)
    ax.plot(T_CTR * 1e3, ym, drawstyle="steps-mid", color=col, lw=lw)

def flavor_fig(prefix, title, name):
    numu, nube, nc = fv[prefix + "_NUMU_CC"], fv[prefix + "_NUBE_CC"], fv[prefix + "_NC"]
    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    for y, col, lab, dy in ((numu, C_NUMU, r"$\nu_\mu$ CC", 1.35),
                            (nube, C_NUBE, r"$\bar\nu_e$ CC", 0.42),
                            (nc, C_NC, r"NC ($\nu+\bar\nu$)", 1.0)):
        m = y > 0
        mstep(ax, y, col)
        if m.any():
            i = np.where(m)[0][-1]
            ax.text(210, y[i] * dy, lab, color=col, fontsize=9, va="center")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(2e-2, 320)
    ax.set_xlabel("time in the 5 Hz cycle  [ms]  (bunch-centric: chain first, then store)")
    ax.set_ylabel(r"interactions / (t$\cdot$yr) per ms of cycle")
    ax.set_title(title, fontsize=10, loc="left")
    return fig, ax

fig, ax = flavor_fig("DEEP", "Flavor composition vs time, on-site convergence hall "
                     "(9.25 km): shares barely move", "chain_flavor_deep")
ax.set_ylim(3e2, 3e7)
ax.text(0.02, 0.05, "cycle-averaged shares 48% : 24% : 29% "
        r"($\nu_\mu$CC : $\bar\nu_e$CC : NC)" + "\nand nearly time-invariant "
        r"($\sigma \propto E$ for both species);" + "\nthe moving part is the "
        "oscillated flavor, next figure",
        transform=ax.transAxes, fontsize=7.8, color="0.35", va="bottom")
save(fig, "chain_flavor_deep")

fig, ax = flavor_fig("NEAR", "Flavor composition vs time at the near hall "
                     "(1.3 km): silent until the store", "chain_flavor_near")
ax.set_ylim(3e4, 3e9)
ax.text(0.02, 0.05, "collider only: the chirp misses this hall (RCS pencils\n"
        "cross 50–68 m below), so the first 7 ms are empty;\n"
        "shares 44% : 27% : 29%, constant",
        transform=ax.transAxes, fontsize=7.8, color="0.35", va="bottom")
save(fig, "chain_flavor_near")

# oscillation flavor evolution
fig, ax = plt.subplots(figsize=(6.6, 3.6))
for key, col, lab, dy in (("FAR_POSC", "#b5541c", "UIUC far hall (197.4 km)", 1.5),
                      ("DEEP_POSC", "#6d3580", "convergence hall (9.25 km)", 1.4),
                      ("NEAR_POSC", "0.45", "near hall (1.3 km)", 0.55)):
    y = fv[key]
    ym = np.ma.masked_invalid(y)
    ax.plot(T_CTR * 1e3, ym, drawstyle="steps-mid", color=col, lw=1.7)
    m = np.isfinite(y)
    i = np.where(m)[0][-1]
    ax.text(210, y[i] * dy, lab, color=col, fontsize=9, va="center")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(2e-2, 320)
ax.set_ylim(2e-10, 3e-5)
ax.set_xlabel("time in the 5 Hz cycle  [ms]  (bunch-centric)")
ax.set_ylabel(r"$\langle P(\nu_\mu\to\nu_\tau)\rangle$ of the arriving fluence")
ax.text(0.13, 9e-7, "RCS1: 47$\\times$ the store's oscillated\nfraction, 560$\\times$ the RCS4-end dip",
        fontsize=8, color="#6d3580", ha="center")
ax.text(30, 2.1e-9, "5 TeV store", fontsize=8, color="0.35")
ax.text(2.6, 8.6e-6, "far hall: south beam, $\\bar\\nu_\\mu\\to\\bar\\nu_\\tau$ (CP mirror)",
        fontsize=7.5, color="#b5541c", ha="center")
ax.text(0.02, 0.04, "vacuum $\\Delta m^2_{31}=2.5\\times10^{-3}$ eV$^2$, amp. 0.95;  "
        "$\\nu_e\\to\\nu_\\tau$ adds ~5% with the same shape;\n"
        "means are dominated by the soft tail of the accepted spectrum",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=7.5, color="0.35")
ax.set_title("Flavor evolution through the cycle: the chirp is the corridor's "
             "oscillation window", fontsize=10, loc="left")
save(fig, "chain_flavor_osc")

# ----------------------------------------------------------------------
# Total on-axis flux vs time, three sites
# ----------------------------------------------------------------------
SITES = (("NEAR", "0.45", "near hall (1.3 km)"),
         ("DEEP", "#6d3580", "convergence hall (9.25 km)"),
         ("FAR", "#b5541c", "UIUC far hall (197.4 km)"))

fig, ax = plt.subplots(figsize=(6.6, 3.8))
for key, col, lab in SITES:
    y = fv[key + "_FLUX"]
    ax.plot(T_CTR * 1e3, np.ma.masked_where(~(y > 0), y),
            drawstyle="steps-mid", color=col, lw=1.6)
    i = np.where(y > 0)[0][-1]
    ax.text(210, y[i], lab, color=col, fontsize=9, va="center")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(2e-2, 320)
ax.set_ylim(3, 3e6)
ax.set_xlabel("time in the 5 Hz cycle  [ms]  (bunch-centric)")
ax.set_ylabel(r"$\nu$ / cm$^2$ / ms on axis (per cycle)")
ax.text(0.3, 5e5, "chirp: only the convergence hall\nis on the RCS pencils' axis",
        fontsize=7.8, color="#6d3580", ha="center")
ax.text(30, 1.1e4, r"5 TeV store: $\propto 1/L^2$", fontsize=8, color="0.35")
ax.set_title("Total on-axis flux vs time at the three sites", fontsize=10,
             loc="left")
save(fig, "chain_flux_total_t")

# ----------------------------------------------------------------------
# Oscillated nutau arrivals vs time, three sites: the (L/E)^2 cancellation
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.6, 3.8))
for key, col, lab in SITES:
    y = fv[key + "_NUTAU"]
    ax.plot(T_CTR * 1e3, np.ma.masked_where(~(y > 0), y),
            drawstyle="steps-mid", color=col, lw=1.6)
dy = {"NEAR": 2.6, "DEEP": 1.0, "FAR": 0.38}
for key, col, lab in SITES:
    y = fv[key + "_NUTAU"]
    i = np.where(y > 0)[0][-1]
    ax.text(210, y[i] * dy[key], lab, color=col, fontsize=9, va="center")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(2e-2, 320)
ax.set_ylim(3e-6, 0.3)
ax.set_xlabel("time in the 5 Hz cycle  [ms]  (bunch-centric)")
ax.set_ylabel(r"oscillation-made $\nu_\tau$ / cm$^2$ / ms (per cycle)")
ax.text(0.35, 6e-2, "chirp: (L/E) blows up as E falls\n"
        "$-$ 30% of the hall's $\\nu_\\tau$ in 3.5% of the cycle",
        fontsize=7.8, color="#6d3580", ha="center")
ax.text(35, 2.5e-3, "store: flux $\\propto 1/L^2$ and P $\\propto L^2$ cancel $-$\n"
        "all three sites receive the same $\\nu_\\tau$ areal density\n"
        "(0.006$-$0.011 /cm$^2$/cycle; residual spread is\nthe soft-tail acceptance)",
        fontsize=7.8, color="0.35", ha="center", va="center")
ax.text(0.98, 0.97, "oscillation only ($\\nu_\\mu\\to\\nu_\\tau$ vacuum; "
        "$\\nu_e$ channel ~5%);\nrock-produced $\\nu_\\tau$ not included "
        "(beam-size study sec. 5)",
        transform=ax.transAxes, ha="right", va="top", fontsize=7.2, color="0.35")
ax.set_title(r"$\nu_\tau$ arrivals vs time: the (L/E)$^2$ cancellation, live",
             fontsize=10, loc="left")
save(fig, "chain_nutau_t")
