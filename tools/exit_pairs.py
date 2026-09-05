#!/usr/bin/env python3
"""Paired exit curves: for one straight at 35 m depth, where its two ends
surface for each tilt, in every bearing.

Each colour is one tilt.  The dotted curve is the exit of the end that goes
UP (the near emergence, s ~ d0/theta); the solid curve is the exit of the end
that goes DOWN and is brought back to the surface by Earth curvature (the far
exit, s ~ 2 R_E theta).  At zero tilt the two coincide.  Faint circles are the
smooth-sphere values; the bold curves use GEBCO 2020 terrain sampled along 48
bearings (data/gebco_bearings.json, tools/fetch_gebco_bearings.py), with
exits over the Great Lakes taken at the water surface.

Output: static/figs/exit_pairs.{svg,pdf}; static/geo/exit_pairs.json.
"""
import json, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import compass_common as cc
STYLE = "/root/.claude/skills/synced/tufte-plots/assets/tufte.mplstyle"
if os.path.exists(STYLE):
    plt.style.use(STYLE)
plt.rcParams["figure.dpi"] = 130

D0 = 35.0                                        # m, straight depth
TILTS = [0.0, 2.0, 5.0, 15.4, 25.0, 35.0, 50.0, 65.0]      # mrad
_cm = plt.get_cmap("plasma")
COL = [_cm(0.02 + 0.86 * i / (len(TILTS) - 1)) for i in range(len(TILTS))]

# ---------------------------------------------------------------- terrain
G = json.load(open(os.path.join(ROOT, "data", "gebco_bearings.json")))
RAD = np.array(G["radii_km"], float)             # km
BRG = [float(b) for b in G["bearings_deg"]]
Z0 = float(G["ip_elev_m"])
LAKE_SURFACE = {"Lake Superior": 183.0, "Lake Michigan": 176.0, "Lake Huron": 176.0,
                "Lake Erie": 174.0, "Lake Ontario": 75.0}
_lk = json.load(open(os.path.join(ROOT, "data", "ne_10m_lakes.geojson")))
GL = []
for f in _lk["features"]:
    nm = f["properties"].get("name")
    if nm in LAKE_SURFACE:
        for r in cc._rings(f["geometry"]):
            GL.append((nm, LAKE_SURFACE[nm], [(c[0], c[1]) for c in r[::max(1, len(r) // 600)]]))

def surface_elev(az, D, e):
    """GEBCO gives the lake BED; a beam 'surfaces' at the water surface."""
    if D < 40 or e is None:
        return e
    la, lo = cc.fwd(az, D)
    for nm, zs, poly in GL:
        if e < zs and cc.pin(lo, la, poly):
            return zs
    return e

PROF = {}                                        # bearing -> ground height rel. to IP tangent plane (m)
for az in BRG:
    el = G["elev"]["%.1f" % az]
    el = np.array([surface_elev(az, D, e) for D, e in zip(RAD, el)], float)
    PROF[az] = el - Z0 - (RAD * 1000.0) ** 2 / (2 * cc.R_E * 1000.0)

def first_exit(theta_end, hg):
    """First range (km) where the beam end (signed slope theta_end, rad) is at or above ground."""
    diff = -D0 + theta_end * RAD * 1000.0 - hg
    idx = np.where(diff >= 0)[0]
    if len(idx) == 0:
        return None
    i = idx[0]
    if i == 0:
        return RAD[0]
    d1, d2 = diff[i - 1], diff[i]
    f = -d1 / (d2 - d1) if d2 != d1 else 0.0
    return float(RAD[i - 1] + f * (RAD[i] - RAD[i - 1]))

RES = []
for th_mrad in TILTS:
    th = th_mrad * 1e-3
    sn, sf = cc.exits_smooth(th, D0)
    near = [first_exit(+th, PROF[az]) for az in BRG]
    far = [first_exit(-th, PROF[az]) for az in BRG]
    RES.append(dict(tilt_mrad=th_mrad, smooth_near_km=round(sn, 3), smooth_far_km=round(sf, 2),
                    near_km=[None if v is None else round(v, 3) for v in near],
                    far_km=[None if v is None else round(v, 2) for v in far]))

def stats(v):
    v = [x for x in v if x is not None]
    return (min(v), float(np.median(v)), max(v)) if v else (None, None, None)

# ---------------------------------------------------------------- figure
D_MIN, D_MAX = 0.3, 1500.0
fig = plt.figure(figsize=(9.8, 10.6))
fig.subplots_adjust(top=0.915, bottom=0.115)
ax = fig.add_subplot(111, projection="polar")
rmap = cc.polar_axes(ax, D_MIN, D_MAX, rings=(0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000),
                     tilt_labels=False, fs=6.6, rlabel_pos=197)
ax.set_yticklabels([])
for d in (0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000):        # ring labels, drawn on top
    ax.text(math.radians(197), rmap(d), "%g km" % d, fontsize=6.6, color="0.35", ha="center", va="center",
            zorder=30, bbox=dict(boxstyle="round,pad=0.15", fc="w", ec="none", alpha=.85))

# DOE land, lakes
faz, fD = cc.fence_polar()
ax.fill(faz, rmap(fD), color="0.80", lw=0.8, ec="0.45", zorder=1)
cc.draw_lakes(ax, rmap, cc.lakes_polar(D_MAX))
for nm, la, lo in cc.LAKE_LABELS:
    D, a = cc.inv(la, lo)
    ax.text(math.radians(a), rmap(D), nm, fontsize=7.2, color="#1f5a7a", ha="center", va="center",
            style="italic", zorder=3, bbox=dict(boxstyle="round,pad=0.1", fc="w", ec="none", alpha=.5))

# the pairs
AZ = np.radians(BRG + [BRG[0]])
circ = np.radians(np.arange(0, 361, 2))
for i, (r, col) in enumerate(zip(RES, COL)):
    sn, sf = r["smooth_near_km"], r["smooth_far_km"]
    near = np.array([np.nan if v is None else v for v in r["near_km"] + r["near_km"][:1]])
    far = np.array([np.nan if v is None else v for v in r["far_km"] + r["far_km"][:1]])
    if r["tilt_mrad"] == 0:
        ax.plot(circ, np.full_like(circ, rmap(sf)), color=col, lw=0.8, alpha=.35, zorder=3)
        ax.plot(AZ, rmap(far), color=col, lw=1.9, zorder=5)
    else:
        ax.plot(circ, np.full_like(circ, rmap(sn)), color=col, lw=0.8, alpha=.35, ls=(0, (1.2, 1.6)), zorder=3)
        ax.plot(circ, np.full_like(circ, rmap(sf)), color=col, lw=0.8, alpha=.35, zorder=3)
        ax.plot(AZ, rmap(near), color=col, lw=1.9, ls=(0, (1.2, 1.6)), zorder=5)
        ax.plot(AZ, rmap(far), color=col, lw=1.9, zorder=5)
    lab = "0 mrad (both ends)" if r["tilt_mrad"] == 0 else "%g mrad" % r["tilt_mrad"]
    az_lab = (215, 160)[i % 2]
    ax.text(math.radians(az_lab), rmap(sf) + 0.03, lab, fontsize=7.2, color=col, fontweight="bold",
            ha="center", va="bottom", zorder=25, bbox=dict(boxstyle="round,pad=0.12", fc="w", ec="none", alpha=.9))

# two worked pairs: the baseline (15.4 mrad, up-north) and an RCS aimed at Lake Huron (50 mrad, up-SW)
def pair(th_mrad, az_far, label_far, label_near, col, dr_far, dr_near, ha_far, ha_near, da_near=0.0):
    i = BRG.index(min(BRG, key=lambda b: abs(((b - az_far) + 180) % 360 - 180)))
    j = BRG.index(min(BRG, key=lambda b: abs(((b - (az_far + 180) % 360) + 180) % 360 - 180)))
    r = next(x for x in RES if x["tilt_mrad"] == th_mrad)
    Df, Dn = r["far_km"][i], r["near_km"][j]
    tf, tn = math.radians(az_far), math.radians((az_far + 180) % 360)
    ax.plot([tf, tf], [0, rmap(Df)], color=col, lw=3.2, alpha=.35, zorder=4, solid_capstyle="butt")
    ax.plot([tn, tn], [0, rmap(Dn)], color=col, lw=3.2, alpha=.35, zorder=4, solid_capstyle="butt")
    ax.plot([tf], [rmap(Df)], "o", ms=9, color=col, mec="w", mew=1.2, zorder=22)
    ax.plot([tn], [rmap(Dn)], "o", ms=9, color=col, mec="w", mew=1.2, mfc="w", zorder=22)
    ax.plot([tn], [rmap(Dn)], "o", ms=4, color=col, zorder=23)
    ax.annotate(label_far % Df, (tf, rmap(Df)), xytext=(tf, rmap(Df) + dr_far), fontsize=7.2, color=col,
                fontweight="bold", ha=ha_far, va="center", zorder=26,
                bbox=dict(boxstyle="round,pad=0.15", fc="w", ec=col, lw=0.6, alpha=.92))
    ax.annotate(label_near % Dn, (tn, rmap(Dn)), xytext=(tn + math.radians(da_near), rmap(Dn) + dr_near), fontsize=7.2, color=col,
                fontweight="bold", ha=ha_near, va="center", zorder=26,
                bbox=dict(boxstyle="round,pad=0.15", fc="w", ec=col, lw=0.6, alpha=.92))
    return Df, Dn

c154 = COL[TILTS.index(15.4)]
c50 = COL[TILTS.index(50.0)]
Df_b, Dn_b = pair(15.4, 180.0, "baseline, down-going end:\nUIUC %.0f km", "baseline, up-going end:\n%.1f km, inside the fence",
                  c154, 0.16, 0.42, "center", "left", da_near=22)
Df_h, Dn_h = pair(50.0, 52.5, "RCS3/4 aimed at Lake Huron:\ndown-going end %.0f km", "its up-going end:\n%.2f km, inside the fence",
                  c50, 0.18, 0.62, "center", "center", da_near=22)

# landmarks
LM = {"UIUC South Farms": ("*", 11, "#0b6e4f"), "Purdue": ("*", 10, "#0b6e4f"),
      "Kettle Moraine SF": ("^", 7, "#4d7a2f"), "Soudan mine (MINOS far)": ("D", 6, "#7a2e21"),
      "Ash River (NOvA far)": ("D", 6, "#7a2e21"), "Chicago": ("o", 4.5, "0.4")}
for nm, (mk, ms, col) in LM.items():
    la, lo, _ = cc.LANDMARKS[nm]
    D, a = cc.inv(la, lo)
    ax.plot([math.radians(a)], [rmap(D)], mk, ms=ms, color=col, mec="w", mew=0.7, zorder=4.5)
def lab(nm, dr, da, ha, col, fs=7):
    la, lo, _ = cc.LANDMARKS[nm]
    D, a = cc.inv(la, lo)
    ax.annotate(nm, (math.radians(a), rmap(D)), xytext=(math.radians(a + da), rmap(D) + dr), fontsize=fs,
                color=col, ha=ha, va="center", zorder=4,
                bbox=dict(boxstyle="round,pad=0.1", fc="w", ec="none", alpha=.7))
lab("Purdue", 0.02, -5, "left", "#0b6e4f")
lab("Kettle Moraine SF", 0.0, 8, "left", "#4d7a2f")
lab("Soudan mine (MINOS far)", -0.08, -6, "right", "#7a2e21")
lab("Ash River (NOvA far)", 0.07, 5, "left", "#7a2e21")
lab("Chicago", 0.0, 7, "left", "0.4")

H = [Line2D([], [], color="0.3", lw=1.9, ls=(0, (1.2, 1.6)), label="up-going end surfaces here (near emergence)"),
     Line2D([], [], color="0.3", lw=1.9, label="down-going end surfaces here (far exit)"),
     Line2D([], [], color="0.3", lw=0.8, alpha=.4, label="same, smooth sphere"),
     Patch(color="0.80", label="Fermilab site"), Patch(color="#3d7ea6", alpha=.55, label="lakes (Natural Earth 10m)"),
     Line2D([], [], color=c154, lw=3.2, alpha=.35, label="one straight: its two ends are opposite each other")]
fig.legend(handles=H, loc="lower center", ncol=3, fontsize=7.2, frameon=False, bbox_to_anchor=(0.5, 0.055),
           handlelength=2.6, columnspacing=1.4)
ax.set_title("Where the two ends of one straight surface, tilt by tilt\n"
             "35 m-deep straight; one colour per tilt $-$ dotted: the end that goes up, solid: the end that goes down; "
             "at 0 mrad they are the same curve",
             fontsize=9.6, pad=22)
fig.text(0.5, 0.015,
         "Bold curves: GEBCO 2020 terrain along 48 bearings (0.25$-$1000 km), exits over the Great Lakes taken at the water surface. "
         "Faint circles: smooth sphere,\n"
         "$s_{near}=R_E(-\\theta+\\sqrt{\\theta^2+2d_0/R_E})\\approx d_0/\\theta$ and "
         "$s_{far}=R_E(\\theta+\\sqrt{\\theta^2+2d_0/R_E})\\approx 2R_E\\theta$. "
         "Log range; rings in km. The fence is 0.9$-$4.9 km out, so the dotted curves show which tilts emerge on site.",
         ha="center", va="center", fontsize=6.9, color="0.35")
for ext in ("pdf", "svg"):
    fig.savefig(os.path.join(ROOT, "static", "figs", "exit_pairs." + ext), bbox_inches="tight")

# ---------------------------------------------------------------- report
for r in RES:
    n = stats(r["near_km"]); f = stats(r["far_km"])
    iN, iS = BRG.index(0.0), BRG.index(180.0)
    r["near_terrain_min_med_max_km"] = [round(x, 2) for x in n]
    r["far_terrain_min_med_max_km"] = [round(x, 1) for x in f]
    r["meridian"] = dict(up_north_end_km=r["near_km"][iN], down_south_end_km=r["far_km"][iS],
                         up_south_end_km=r["near_km"][iS], down_north_end_km=r["far_km"][iN])
    print("%5.1f mrad  smooth near %6.2f far %6.1f | terrain near %5.2f-%5.2f (med %5.2f)  far %6.1f-%6.1f (med %6.1f) | "
          "meridian upN %s downS %s upS %s downN %s" % (r["tilt_mrad"], r["smooth_near_km"], r["smooth_far_km"],
          n[0], n[2], n[1], f[0], f[2], f[1], r["meridian"]["up_north_end_km"], r["meridian"]["down_south_end_km"],
          r["meridian"]["up_south_end_km"], r["meridian"]["down_north_end_km"]))
fd = {("%.1f" % b): round(cc.fence_distance(b)) for b in BRG}
json.dump(dict(straight_depth_m=D0, ip_elev_m=Z0, bearings_deg=BRG, fence_m=fd, tilts=RES,
               worked_pairs=dict(baseline_15p4_up_north=dict(far_south_km=Df_b, near_north_km=Dn_b),
                                 rcs_50_to_lake_huron_az52=dict(far_km=Df_h, near_km=Dn_h))),
          open(os.path.join(ROOT, "static", "geo", "exit_pairs.json"), "w"), indent=1)
print("baseline pair: far S %.1f km, near N %.2f km; Huron pair: far %.1f, near %.2f" % (Df_b, Dn_b, Df_h, Dn_h))
