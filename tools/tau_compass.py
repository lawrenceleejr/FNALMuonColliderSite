#!/usr/bin/env python3
"""Tau appearance around the compass, one polar panel per stage, plus the
on-axis nu_tau CC rate per tonne-year against DUNE.

Colour at (bearing, range L) is the flux-weighted oscillation probability a
detector ON THE AXIS of that stage's beam would see at that range,
    Pbar(L) = < 0.95 sin^2(1.267 dm^2 L / E_nu) >  over the stage's point spectrum,
dm^2 = 2.5e-3 eV^2.  The RCS straights are 1/gamma-sharp pencils, so a point on
their axis sees the hard forward spectrum f(x) = 2x^2(3-2x) (x = E_nu/E_mu,
mean 0.7); the collider's plume is divergence-smeared (sigma_theta = 7/gamma),
so a point on its axis sees the whole-plume spectrum g(y) = 5/3 - 3y^2 + 4y^3/3
(mean 0.35) cut off below y = 1/(1+(gamma sigma)^2) = 1/50 by the Gaussian tail,
and its on-axis density is diluted by DIV_SUPP = 2 gamma^2 sigma^2 ~ 101.  Each
ramp is averaged with decays per turn ~ 1/E (uniform in ln E).  Pbar depends on
range only, so each panel is a set of rings; the geography says what sits at
each ring.

Rates: nu_tau + nubar_tau CC per tonne-year for a point detector on the axis,
    R(L) = sum_E w_E N_E gamma_E^2/(pi L^2) * int spec(x) P(L, xE) [sigma_nu + sigma_nubar](xE) xi_tau(xE) dx * N_A/tonne,
with N_E the north-aimed decays per year per sign from static/geo/chain_timing.json,
CSMS cross-sections (CC fraction 0.71 of the CC+NC totals in corridor_layout.py)
and the tau-threshold suppression table of uiuc_chain.py.  Because the on-axis
flux falls as 1/L^2 while P rises as L^2, R is baseline-independent inside the
quadratic regime -- the (L/E)^2 cancellation the timing study noted.  DUNE's
expected nu_tau CC yields (40 kt fiducial, 1300 km) are drawn as per-tonne lines.

Output: static/figs/tau_compass.{svg,pdf}; static/geo/tau_compass.json.
"""
import json, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import compass_common as cc
STYLE = "/root/.claude/skills/synced/tufte-plots/assets/tufte.mplstyle"
if os.path.exists(STYLE):
    plt.style.use(STYLE)
plt.rcParams["figure.dpi"] = 130

# ---------------------------------------------------------------- beam and physics inputs
CT = json.load(open(os.path.join(ROOT, "static", "geo", "chain_timing.json")))
DEC = {s["stage"]: s["decays_north_per_cycle"] for s in CT["stages"]}       # per cycle, per sign
SEC_YR = 1.2e7                                                               # s/yr, as corridor_layout.py
CYCLES_YR = CT["rep_rate_Hz"] * SEC_YR
STAGES = [("RCS1", 63.0, 314.0, DEC["RCS1"], False), ("RCS2", 314.0, 750.0, DEC["RCS2"], False),
          ("RCS3", 750.0, 1500.0, DEC["RCS3"], False), ("RCS4", 1500.0, 5000.0, DEC["RCS4"], False),
          ("Collider store", 5000.0, 5000.0, DEC["collider"], True)]
M_MU = 0.1056584                       # GeV
SIG_THETA = 0.15e-3                    # rad, MINT smeared divergence
GAMMA_COLL = 5000.0 / M_MU
DIV_SUPP = 2 * GAMMA_COLL ** 2 * SIG_THETA ** 2
N_A_T = 6.022e29                       # nucleons per tonne
DM2, P_AMP = 2.5e-3, 0.95
L_MAX_PER_GEV = math.pi / 2 / (1.267 * DM2)          # km per GeV: first oscillation maximum
EARTH_DIAM = 2 * cc.R_E
# DUNE reference (TDR vol. II, arXiv:2002.03005, sec. 4.1.1.3): nu_tau CC interactions per year at the
# 40 kt fiducial far detector, 1300 km, 1.2 MW, before detector efficiency: ~130 in the standard
# CP-optimised configuration, ~1000 with two NuMI-like parabolic horns (tau-optimised).
DUNE_MASS_T = 40000.0
DUNE_YIELDS = {"DUNE, CP-optimised beam: ~130 $\\nu_\\tau$ CC / yr in 40 kt": 130.0,
               "DUNE, $\\tau$-optimised horns: ~1000 / yr in 40 kt": 1000.0}

def posc(E, L):
    return P_AMP * np.sin(1.267 * DM2 * L / E) ** 2
def g_plume(y):
    return 5.0 / 3.0 - 3.0 * y ** 2 + 4.0 * y ** 3 / 3.0
def f_axis(x):
    return 2.0 * x ** 2 * (3.0 - 2.0 * x)
Y_MIN_COLL = 1.0 / (1.0 + (GAMMA_COLL * SIG_THETA) ** 2)
def g_point_collider(y):
    return np.where(y >= Y_MIN_COLL, g_plume(y), 0.0)
POINT = {nm: (g_point_collider if smeared else f_axis) for nm, _, _, _, smeared in STAGES}

# CSMS (arXiv:1106.3723) CC+NC totals per nucleon, isoscalar, from corridor_layout.py; CC ~ 0.71 of the total
_SIG_NU = [(1e1, 8.2e-38), (1e3, 8.2e-36), (2e3, 15.8e-36), (5e3, 35.6e-36), (1e4, 62e-36)]
_SIG_NUB = [(1e1, 4.8e-38), (1e3, 4.8e-36), (2e3, 9.4e-36), (5e3, 22.8e-36), (1e4, 42e-36)]
CC_FRAC = 0.71
def sig_cc(tab, E):
    xs = np.log([e for e, _ in tab]); ys = np.log([s for _, s in tab])
    E = np.asarray(E, float)
    out = np.exp(np.interp(np.log(np.maximum(E, 10.0)), xs, ys)) * CC_FRAC
    return np.where(E < 10.0, out * E / 10.0, out)
_TAU_R = [(3.6, 0.0), (5, 0.05), (10, 0.18), (20, 0.34), (50, 0.52), (100, 0.62), (300, 0.74), (1000, 0.85), (5000, 0.95)]
def tau_thresh(E):
    xs = np.log([e for e, _ in _TAU_R]); ys = [r for _, r in _TAU_R]
    return np.interp(np.log(np.maximum(np.asarray(E, float), 3.61)), xs, ys, left=0.0, right=0.95)

Y = np.linspace(0.0005, 0.9995, 1000)
def _egrid(E1, E2):
    E = np.geomspace(E1, E2, 160) if E2 > E1 else np.array([E1])
    return E, np.ones(len(E)) / len(E)                     # decays uniform in ln E over a linear ramp
def pbar(L, E1, E2, spec, e_cut=None):
    E, w = _egrid(E1, E2)
    W = w[:, None] * spec(Y)[None, :]
    Enu = Y[None, :] * E[:, None]
    if e_cut is not None:
        W = np.where(Enu >= e_cut, W, 0.0)
    W = W / W.sum()
    return np.array([(W * posc(Enu, Lk)).sum() for Lk in np.atleast_1d(L)])
def mean_Enu(E1, E2, spec):
    E, w = _egrid(E1, E2)
    W = w[:, None] * spec(Y)[None, :]
    return float((W * (Y[None, :] * E[:, None])).sum() / W.sum())
def rate_per_t_yr(L, nm):
    """nu_tau + nubar_tau CC per tonne-year for a point detector on the axis at range L (km)."""
    _, E1, E2, n_cycle, smeared = next(s for s in STAGES if s[0] == nm)
    E, w = _egrid(E1, E2)
    S = POINT[nm](Y); S = S / np.trapezoid(S, Y)
    Enu = Y[None, :] * E[:, None]
    N_sign = n_cycle * CYCLES_YR
    out = []
    for Lk in np.atleast_1d(L):
        dens = N_sign * (E / M_MU) ** 2 / (math.pi * (Lk * 1e5) ** 2)      # on-axis decays / cm^2 per sign
        if smeared:
            dens = dens / DIV_SUPP
        per_pair = np.trapezoid(S[None, :] * posc(Enu, Lk) * (sig_cc(_SIG_NU, Enu) + sig_cc(_SIG_NUB, Enu))
                            * tau_thresh(Enu), Y, axis=1)                   # cm^2 per (mu+, mu-) decay pair
        out.append(float((w * dens * per_pair).sum() * N_A_T))
    return np.array(out)

L = np.geomspace(10.0, 1500.0, 400)
CURVES = {}
for nm, E1, E2, _, _ in STAGES:
    CURVES[nm] = dict(point=pbar(L, E1, E2, POINT[nm]), plane=pbar(L, E1, E2, g_plume, e_cut=3.5),
                      rate=rate_per_t_yr(L, nm), E1=E1, E2=E2, Enu=mean_Enu(E1, E2, POINT[nm]))
    CURVES[nm]["Lmax"] = L_MAX_PER_GEV * CURVES[nm]["Enu"]

# ---------------------------------------------------------------- figure
D_MIN, D_MAX = 10.0, 1500.0
RINGS = (20, 50, 100, 200, 500, 1000)
LO, HI = -9.0, -2.0
CMAP = plt.get_cmap("viridis")
NORM = Normalize(LO, HI)
fig = plt.figure(figsize=(13.6, 10.6))
gs = fig.add_gridspec(2, 3, left=0.025, right=0.985, top=0.875, bottom=0.155, wspace=0.14, hspace=0.36)
lakes = cc.lakes_polar(D_MAX, min_vertices=40, decimate_to=150)
th_e = np.linspace(0, 2 * np.pi, 721)
LE = np.geomspace(D_MIN, D_MAX, 601)

LABEL = {"UIUC South Farms": ("UIUC", (6, -9), "left"), "Purdue": ("Purdue", (6, 2), "left"),
         "Soudan mine (MINOS far)": ("Soudan", (-6, -7), "right"), "Ash River (NOvA far)": ("Ash River", (-6, 5), "right"),
         "SURF (DUNE far)": ("SURF", (4, -9), "left"), "Chicago": ("Chicago", (5, 0), "left")}
def fmt_e(v):
    e = int(math.floor(math.log10(v))); m = v / 10 ** e
    return r"%.1f\times10^{%d}" % (m, e)
axes = []
for k, (nm, E1, E2, _, smeared) in enumerate(STAGES):
    ax = fig.add_subplot(gs[k // 3, k % 3], projection="polar")
    axes.append(ax)
    rmap = cc.polar_axes(ax, D_MIN, D_MAX, RINGS, tilt_labels=False, fs=5.8, label_color="0.15",
                         ring_color="w", rlabel_pos=250)
    ax.set_yticklabels([])
    P = pbar(0.5 * (LE[1:] + LE[:-1]), E1, E2, POINT[nm])
    TH, RR = np.meshgrid(th_e, rmap(LE), indexing="ij")
    C = np.tile(np.log10(np.maximum(P, 1e-12))[None, :], (len(th_e) - 1, 1))
    ax.pcolormesh(TH, RR, C, cmap=CMAP, norm=NORM, shading="flat", zorder=0, rasterized=True)
    # contours at 1, 2, 5 x 10^k; decade rings labelled
    Pc = CURVES[nm]["point"]
    for dec in range(-10, 0):
        for m, lw, a in ((1.0, 0.75, 0.9), (2.0, 0.35, 0.7), (5.0, 0.35, 0.7)):
            lev = m * 10.0 ** dec
            if Pc.min() < lev < Pc.max():
                Ld = float(np.interp(math.log10(lev), np.log10(Pc), L))
                ax.plot(th_e, np.full_like(th_e, rmap(Ld)), color="w", lw=lw, alpha=a, zorder=3)
                if m == 1.0:
                    ax.text(math.radians(300), rmap(Ld), "$10^{%d}$" % dec, fontsize=6.2, color="w", ha="center",
                            va="bottom", zorder=30)
    for lnm, great, az, Ds in lakes:
        if not great:
            continue
        ax.fill(az, rmap(Ds), color="w", alpha=.18, lw=0.6, ec="w", zorder=2)
    ax.plot(th_e, np.full_like(th_e, rmap(198.4)), color="#ff9d3a", lw=1.4, ls=(0, (4, 2.5)), zorder=6)

    for d in RINGS:                                             # ring labels, drawn on top
        ax.text(math.radians(250), rmap(d), "%g km" % d, fontsize=5.8, color="0.15", ha="center", va="center",
                zorder=30, bbox=dict(boxstyle="round,pad=0.15", fc="w", ec="none", alpha=.8))
    for lm, (la, lo, cat) in cc.LANDMARKS.items():
        if lm not in LABEL:
            continue
        D, a = cc.inv(la, lo)
        mk = {"university": "*", "mine": "D", "public": "^", "city": "o", "federal": "s"}[cat]
        ax.plot([math.radians(a)], [rmap(D)], mk, ms=8 if mk == "*" else 5, color="w", mec="0.1", mew=0.6, zorder=5)
        if k == 0:
            txt, off, ha = LABEL[lm]
            ax.annotate(txt, (math.radians(a), rmap(D)), xytext=off, textcoords="offset points", fontsize=6.4,
                        color="0.1", ha=ha, va="center", zorder=4.5,
                        bbox=dict(boxstyle="round,pad=0.12", fc="w", ec="none", alpha=.75))
    if k == 0:
        ax.text(math.radians(118), rmap(198.4) + 0.04, "198 km: UIUC /\ncollider circle", fontsize=6.2,
                color="#b5541c", ha="center", va="bottom", zorder=30, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.15", fc="w", ec="none", alpha=.8))
    c = CURVES[nm]
    ttl = ("%s   %g$\\rightarrow$%g GeV" % (nm, E1, E2)) if E2 > E1 else ("%s   %g GeV" % (nm, E1))
    ax.set_title(ttl + "   $\\langle E_\\nu\\rangle_{axis}\\approx$ %.0f GeV" % c["Enu"], fontsize=8.8, pad=9)
    rate_uiuc = float(rate_per_t_yr(198.4, nm)[0])
    ax.text(0.5, -0.07,
            "1st oscillation maximum at $L=%s$ km (%.0f Earth diameters): this whole map is $\\bar P\\propto L^2$\n"
            "UIUC 198 km: $\\bar P=%s$   $\\cdot$   Soudan 736 km: $%s$\n"
            "on-axis $\\nu_\\tau+\\bar\\nu_\\tau$ CC: %.3g per tonne$\\cdot$yr, at any range"
            % (fmt_e(c["Lmax"]), c["Lmax"] / EARTH_DIAM, fmt_e(float(pbar(198.4, E1, E2, POINT[nm])[0])),
               fmt_e(float(pbar(736.0, E1, E2, POINT[nm])[0])), rate_uiuc),
            transform=ax.transAxes, fontsize=6.3, color="0.25", ha="center", va="top", linespacing=1.35)

# sixth cell: Pbar(L) on top, rate per tonne-year below, DUNE lines
sub = gs[1, 2].subgridspec(2, 1, hspace=0.55)
ax6 = fig.add_subplot(sub[0]); ax7 = fig.add_subplot(sub[1])
cols = ["#7b3294", "#c2a5cf", "#5aae61", "#1b7837", "#b5541c"]
for (nm, E1, E2, _, _), col in zip(STAGES, cols):
    c = CURVES[nm]
    ax6.plot(L, c["point"], color=col, lw=1.8, label=nm)
    ax6.plot(L, c["plane"], color=col, lw=0.9, ls="--", alpha=.8)
    ax7.plot(L, c["rate"], color=col, lw=1.8, label=nm)
tot = sum(CURVES[nm]["rate"] for nm, *_ in STAGES)
ax7.plot(L, tot, color="0.15", lw=1.2, ls=":", label="all five, co-tilted chain")
for Lm, txt, y in ((198.4, "UIUC", 1.6e-10), (655, "L. Superior", 1.6e-10), (736, "Soudan", 4e-9),
                   (811, "Ash River", 1.6e-10), (1289, "SURF", 1.6e-10)):
    for a in (ax6, ax7):
        a.axvline(Lm, color="0.75", lw=0.6, ls=":", zorder=0)
    ax6.text(Lm * 0.97, y, txt, rotation=90, fontsize=6.2, color="0.4", ha="right", va="bottom")
for (lab, n), ls in zip(DUNE_YIELDS.items(), ("--", "-.")):
    r = n / DUNE_MASS_T
    ax7.axhline(r, color="0.35", lw=1.0, ls=ls, zorder=1)
    ax7.text(1450, r * 1.2, lab, fontsize=6.2, color="0.3", ha="right", va="bottom",
             bbox=dict(boxstyle="round,pad=0.1", fc="w", ec="none", alpha=.85), zorder=20)
ax6.set_xscale("log"); ax6.set_yscale("log"); ax6.set_xlim(10, 1500); ax6.set_ylim(1e-10, 5e-2)
ax6.set_ylabel(r"$\bar P(\bar\nu_\mu\to\bar\nu_\tau)$ on axis", fontsize=7.5)
ax6.tick_params(labelsize=6.5); ax6.set_xticklabels([])
ax6.legend(fontsize=6.3, frameon=False, loc="upper left", ncol=2, columnspacing=0.9,
           title="solid: on the beam axis   dashed: whole plane, $E_\\nu>$3.5 GeV", title_fontsize=6.0, alignment="left")
ax6.set_title("the same curves, all stages", fontsize=8.8, pad=8)
ax7.set_xscale("log"); ax7.set_yscale("log"); ax7.set_xlim(10, 1500)
ax7.set_ylim(1e-4, 30.0)
ax7.set_xlabel("range L (km)", fontsize=8)
ax7.set_ylabel(r"$\nu_\tau+\bar\nu_\tau$ CC per tonne$\cdot$yr, on axis", fontsize=7.5)
ax7.tick_params(labelsize=6.5)
lg = ax7.legend(fontsize=6.0, frameon=True, loc="upper left", ncol=3, columnspacing=0.9, handlelength=1.8,
                framealpha=0.9, edgecolor="none", borderpad=0.4)
lg.set_zorder(20)
ax7.set_title("rate for a point detector on the axis: baseline-independent", fontsize=8.8, pad=8)

sm = plt.cm.ScalarMappable(norm=NORM, cmap=CMAP)
cax = fig.add_axes([0.075, 0.06, 0.545, 0.016])
cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
cb.set_ticks(range(int(LO), int(HI) + 1))
cb.set_ticklabels(["$10^{%d}$" % k for k in range(int(LO), int(HI) + 1)])
cb.ax.tick_params(labelsize=7)
cb.set_label(r"$\bar\nu_\tau$ appearance fraction $\bar P$ of the stage's $\bar\nu_\mu$ flux, on the beam axis, "
             r"at that range ($\Delta m^2=2.5\times10^{-3}$ eV$^2$, amplitude 0.95); contours at 1, 2, 5 $\times10^k$", fontsize=7.6)
fig.suptitle("Tau appearance around the compass, stage by stage", fontsize=11.5, y=0.975)
fig.text(0.5, 0.937, "colour = fraction of the stage's $\\bar\\nu_\\mu$ flux that has become $\\bar\\nu_\\tau$ by the range where "
         "the beam surfaces at that bearing (it depends on range only: rings). White: lakes. Orange dashed: the 198 km "
         "circle the collider and the co-tilted chain exit on.", ha="center", fontsize=8.6, color="0.3")
fig.text(0.665, 0.083,
         "Point spectra: the RCS pencils are $1/\\gamma$-sharp, so a point on the axis sees the hard forward spectrum "
         "($\\langle E_\\nu\\rangle = 0.7\\,E_\\mu$);\nthe collider plume is smeared to $\\sigma_\\theta = 7/\\gamma$, so its axis sees the "
         "whole-plume spectrum ($\\langle E_\\nu\\rangle = 0.35\\,E_\\mu$, cut at $E_\\nu/E_\\mu = 1/50$)\nat $1/101$ of the pencil density. "
         "Rates: north-aimed decays per cycle from chain_timing.json, 5 Hz, $1.2\\times10^7$ s/yr, both signs,\n"
         "CSMS CC cross-sections, $\\tau$-threshold suppression; on-axis flux $\\propto 1/L^2$ and $\\bar P\\propto L^2$ cancel, so the "
         "rate is flat in $L$.\nDUNE lines: TDR vol. II $\\S$4.1.1.3 $\\nu_\\tau$ CC interactions per year before detector efficiency (130 CP-optimised,\n"
         "~1000 $\\tau$-optimised) in the 40 kt fiducial far detector at 1300 km, per tonne.",
         ha="left", va="top", fontsize=6.6, color="0.35")
for ext in ("pdf", "svg"):
    fig.savefig(os.path.join(ROOT, "static", "figs", "tau_compass." + ext), bbox_inches="tight")

# ---------------------------------------------------------------- report
KEY = [(198.4, "UIUC"), (320, "25 mrad"), (447, "35 mrad"), (655, "Lake Superior"), (736, "Soudan"), (811, "Ash River"), (1289, "SURF")]
out = dict(dm2_eV2=DM2, amplitude=P_AMP, sec_per_year=SEC_YR, cycles_per_year=CYCLES_YR, div_supp=DIV_SUPP,
           spectrum="point: RCS f(x)=2x^2(3-2x) on axis; collider g(y)=5/3-3y^2+4y^3/3 for y>1/50 (sigma_theta=7/gamma); "
                    "plane: g(y) with E_nu>3.5 GeV; decays uniform in ln E over each ramp",
           rate_model="nu_tau+nubar_tau CC per tonne-year, point detector on axis; CSMS CC (0.71 x CC+NC totals), "
                      "tau threshold table from uiuc_chain.py; both muon signs",
           dune=dict(source="DUNE TDR vol. II arXiv:2002.03005 sec. 4.1.1.3; interactions before detector efficiency, 1.2 MW",
                     fiducial_mass_t=DUNE_MASS_T, baseline_km=1300, nu_tau_cc_per_year={"CP-optimised": 130.0, "tau-optimised horns": 1000.0},
                     per_tonne_year={"CP-optimised": 130.0 / DUNE_MASS_T, "tau-optimised horns": 1000.0 / DUNE_MASS_T}),
           L_km=[round(x, 1) for x in L[::20]], stages={})
for nm, E1, E2, n_cycle, smeared in STAGES:
    c = CURVES[nm]
    out["stages"][nm] = dict(E_GeV=[E1, E2], decays_north_per_cycle_per_sign=n_cycle, smeared=smeared,
                             mean_Enu_axis_GeV=round(c["Enu"], 1), first_osc_max_km=round(c["Lmax"]),
                             first_osc_max_earth_diameters=round(c["Lmax"] / EARTH_DIAM, 1),
                             pbar_point=[float("%.3e" % v) for v in c["point"][::20]],
                             pbar_plane_above_3p5GeV=[float("%.3e" % v) for v in c["plane"][::20]],
                             rate_per_t_yr=[float("%.3e" % v) for v in c["rate"][::20]],
                             at={lab: float("%.3e" % pbar(Lk, E1, E2, POINT[nm])[0]) for Lk, lab in KEY},
                             rate_at={lab: float("%.3e" % rate_per_t_yr(Lk, nm)[0]) for Lk, lab in KEY})
    print("%-15s <Enu> %6.0f GeV  Lmax %.2e km (%.0f Earth diam) | rate/t/yr at UIUC %.3g, Soudan %.3g | Pbar UIUC %.1e Soudan %.1e" % (
        nm, c["Enu"], c["Lmax"], c["Lmax"] / EARTH_DIAM, out["stages"][nm]["rate_at"]["UIUC"], out["stages"][nm]["rate_at"]["Soudan"],
        out["stages"][nm]["at"]["UIUC"], out["stages"][nm]["at"]["Soudan"]))
tot_uiuc = sum(out["stages"][nm]["rate_at"]["UIUC"] for nm, *_ in STAGES)
print("all five at UIUC: %.3g per t-yr; DUNE per t-yr: %s" % (tot_uiuc, {k: "%.2e" % (v / DUNE_MASS_T) for k, v in DUNE_YIELDS.items()}))
out["total_rate_per_t_yr_all_stages_on_axis"] = tot_uiuc
json.dump(out, open(os.path.join(ROOT, "static", "geo", "tau_compass.json"), "w"), indent=1)
