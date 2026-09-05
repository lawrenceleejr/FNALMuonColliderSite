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
D_CAP = 200.0                                    # m, the deepest any ring here may be
TILTS = [0.0, 2.0, 5.0, 15.4, 25.0, 35.0, 50.0, 65.0]      # mrad
_cm = plt.get_cmap("plasma")
COL = [_cm(0.02 + 0.86 * i / (len(TILTS) - 1)) for i in range(len(TILTS))]

# ---------------------------------------------------------------- terrain
G = json.load(open(os.path.join(ROOT, "data", "gebco_bearings.json")))
RAD = np.array(G["radii_km"], float)             # km
BRG = sorted(float(k) for k in G["elev"].keys())        # every bearing the cache holds (2.5 deg steps)
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

# Near field at 10 m: USGS NED (tools/fetch_ned_near.py -> data/ned_near.json) replaces GEBCO inside
# NED_MAX km; the two grids are tied together by their mean difference over the overlap (1-5 km).
NED_PATH = os.path.join(ROOT, "data", "ned_near.json")
NED = json.load(open(NED_PATH)) if os.path.exists(NED_PATH) else None
NED_MAX, SEAM_OFFSET_M, NED_USED = 6.0, 0.0, False
if NED and NED.get("ip_elev_m") is not None and all(v is not None for k in NED["elev"] for v in NED["elev"][k]):
    NRAD = np.array(NED["radii_km"], float)
    diffs = []
    for az in BRG:
        key = "%.1f" % az
        if key not in NED["elev"]:
            continue
        gi = {float(D): e for D, e in zip(RAD, G["elev"][key]) if 1.0 <= D <= 5.0}
        ni = {float(D): e for D, e in zip(NRAD, NED["elev"][key])}
        diffs += [gi[D] - ni[D] for D in gi if D in ni and gi[D] is not None and ni[D] is not None]
    SEAM_OFFSET_M = float(np.mean(diffs)) if diffs else 0.0          # GEBCO minus NED where both exist
    NED_USED = all("%.1f" % az in NED["elev"] for az in BRG)
    Z0 = float(NED["ip_elev_m"])                                      # reference the profiles to the 10 m grade at the IP
    print("NED near field: %d bearings, seam offset GEBCO-NED = %+.2f m (sd %.2f m over %d pairs), IP grade %.1f m" % (
        len(NED["elev"]), SEAM_OFFSET_M, float(np.std(diffs)) if diffs else 0.0, len(diffs), Z0))

def merged_profile(az):
    """(ranges km, elevations m) along a bearing: NED to NED_MAX km, GEBCO (shifted onto NED) beyond."""
    key = "%.1f" % az
    if NED_USED:
        r_n, e_n = list(NRAD), list(NED["elev"][key])
        far = [(D, e - SEAM_OFFSET_M if e is not None else None) for D, e in zip(RAD, G["elev"][key]) if D > NED_MAX]
        return np.array(r_n + [d for d, _ in far], float), [e for e in e_n] + [e for _, e in far]
    return RAD, list(G["elev"][key])

PROF, RADS = {}, {}                              # bearing -> (ground height rel. to IP tangent plane (m), ranges km)
for az in BRG:
    rad, el = merged_profile(az)
    el = np.array([surface_elev(az, D, e) for D, e in zip(rad, el)], float)
    PROF[az] = el - Z0 - (rad * 1000.0) ** 2 / (2 * cc.R_E * 1000.0)
    RADS[az] = rad

def first_exit(theta_end, hg, rad, d0=D0):
    """First range (km) where the beam end (signed slope theta_end, rad) is at or above ground."""
    diff = -d0 + theta_end * rad * 1000.0 - hg
    idx = np.where(diff >= 0)[0]
    if len(idx) == 0:
        return None
    i = idx[0]
    if i == 0:
        return float(rad[0])
    d1, d2 = diff[i - 1], diff[i]
    f = -d1 / (d2 - d1) if d2 != d1 else 0.0
    return float(rad[i - 1] + f * (rad[i] - rad[i - 1]))

RES = []
for th_mrad in TILTS:
    th = th_mrad * 1e-3
    sn, sf = cc.exits_smooth(th, D0)
    near = [first_exit(+th, PROF[az], RADS[az]) for az in BRG]
    far = [first_exit(-th, PROF[az], RADS[az]) for az in BRG]
    RES.append(dict(tilt_mrad=th_mrad, smooth_near_km=round(sn, 3), smooth_far_km=round(sf, 2),
                    near_km=[None if v is None else round(v, 3) for v in near],
                    far_km=[None if v is None else round(v, 2) for v in far]))

RES_CAP = []                                     # the same, for a straight at the 200 m ceiling
for th_mrad in TILTS:
    th = th_mrad * 1e-3
    sn, sf = cc.exits_smooth(th, D_CAP)
    near = [first_exit(+th, PROF[az], RADS[az], D_CAP) for az in BRG]
    far = [first_exit(-th, PROF[az], RADS[az], D_CAP) for az in BRG]
    RES_CAP.append(dict(tilt_mrad=th_mrad, smooth_near_km=round(sn, 3), smooth_far_km=round(sf, 2),
                        near_km=[None if v is None else round(v, 3) for v in near],
                        far_km=[None if v is None else round(v, 2) for v in far]))

def stats(v):
    v = [x for x in v if x is not None]
    return (min(v), float(np.median(v)), max(v)) if v else (None, None, None)

# ---------------------------------------------------------------- figure
D_MIN, D_MAX = 0.2, 1000.0
fig = plt.figure(figsize=(9.8, 10.6))
fig.subplots_adjust(top=0.915, bottom=0.135)
ax = fig.add_subplot(111, projection="polar")
rmap = cc.polar_axes(ax, D_MIN, D_MAX, rings=(0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500),
                     tilt_labels=False, fs=6.6, rlabel_pos=197)
ax.set_yticklabels([])
for d in (0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500):        # ring labels, drawn on top
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

# the pairs: periodic Catmull-Rom (cubic Hermite) interpolation through the sampled bearings, in log range,
# evaluated every half degree -- local, C1-smooth, no ringing at lake edges
AZ_FINE = np.radians(np.arange(0, 360.5, 0.5))
def smooth(vals):
    v = np.array([np.nan if x is None else x for x in vals], float)
    if np.isnan(v).any():
        return np.full_like(AZ_FINE, np.nan)
    y = np.log(v); n = len(y); step = 360.0 / n
    out = np.empty_like(AZ_FINE)
    for k, a in enumerate(np.degrees(AZ_FINE)):
        u = (a % 360.0) / step; i = int(np.floor(u)); t = u - i
        y0, y1, y2, y3 = y[(i - 1) % n], y[i % n], y[(i + 1) % n], y[(i + 2) % n]
        out[k] = 0.5 * ((2 * y1) + (-y0 + y2) * t + (2 * y0 - 5 * y1 + 4 * y2 - y3) * t * t + (-y0 + 3 * y1 - 3 * y2 + y3) * t ** 3)
    return np.exp(out)
circ = np.radians(np.arange(0, 360.5, 0.5))
for i, (r, col) in enumerate(zip(RES, COL)):
    sn, sf = r["smooth_near_km"], r["smooth_far_km"]
    near, far = smooth(r["near_km"]), smooth(r["far_km"])
    if r["tilt_mrad"] == 0:
        ax.plot(circ, np.full_like(circ, rmap(sf)), color=col, lw=0.8, alpha=.35, zorder=3)
        ax.plot(AZ_FINE, rmap(far), color=col, lw=1.9, zorder=5)
    else:
        ax.plot(circ, np.full_like(circ, rmap(sn)), color=col, lw=0.8, alpha=.35, ls=(0, (1.2, 1.6)), zorder=3)
        ax.plot(circ, np.full_like(circ, rmap(sf)), color=col, lw=0.8, alpha=.35, zorder=3)
        ax.plot(AZ_FINE, rmap(near), color=col, lw=1.9, ls=(0, (1.2, 1.6)), zorder=5)
        ax.plot(AZ_FINE, rmap(far), color=col, lw=1.9, zorder=5)
    lab = "0 mrad (both ends)" if r["tilt_mrad"] == 0 else "%g mrad" % r["tilt_mrad"]
    az_lab = (215, 160)[i % 2]
    ax.text(math.radians(az_lab), rmap(sf) + 0.03, lab, fontsize=7.2, color=col, fontweight="bold",
            ha="center", va="bottom", zorder=25, bbox=dict(boxstyle="round,pad=0.12", fc="w", ec="none", alpha=.9))

# the accelerator rings themselves, in plan: racetracks with the east straight on the meridian, body to the west
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    import corridor_layout as cl
RINGS = [("collider ring, C 11.0 km", cl.y_coll - cl.IP_Y, cl.LS_COLLIDER, cl.R_COLL, "#0b6e4f", 10.0),
         ("RCS3/4, C 14.7 km", cl.y_rcs34 - cl.IP_Y, cl.LS_RCS34, cl.R_RCS34, "#7a2e21", -12.0),
         ("RCS1/2, C 6.3 km", cl.y_rcs12 - cl.IP_Y, cl.LS_RCS12, cl.R_RCS12, "#1c5a96", -8.0)]
for nm, yc, Ls, Rr, colr, da in RINGS:
    pts = cl.racetrack(yc, Ls, Rr, n_arc=90)                      # (x east, y north) metres from the IP
    az_r = np.unwrap(np.array([math.atan2(x, y) for x, y in pts]))
    D_r = np.array([math.hypot(x, y) / 1000.0 for x, y in pts])
    ax.plot(az_r, rmap(D_r), color=colr, lw=1.3, alpha=.9, zorder=6)
    if nm.startswith("RCS3/4"):
        xw, yw = -Rr, yc + Ls / 2 + Rr                                # label at the north arc apex
    else:
        xw, yw = -2 * Rr, yc                                            # label near the west straight's midpoint
    ax.text(math.atan2(xw, yw) + math.radians(da if not nm.startswith("RCS3/4") else 0.0),
            rmap(math.hypot(xw, yw) / 1000.0) + 0.04, nm, fontsize=6.4, color=colr,
            ha="center", va="bottom", zorder=24, bbox=dict(boxstyle="round,pad=0.12", fc="w", ec="none", alpha=.85))

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

# the depth ceiling: a level straight at 200 m surfaces at sqrt(2 R d0) = 50 km both ways, and no up-going end of any
# tilt can surface beyond that curve; a straight that is at once at UIUC and farther north than 12.8 km would have to
# be deeper than the ceiling (s_near s_far = 2 R d0)
lev_cap = smooth(RES_CAP[0]["far_km"])
ax.plot(AZ_FINE, rmap(lev_cap), color="0.3", lw=1.3, ls=(0, (4, 2)), zorder=6)
ax.text(math.radians(252), rmap(RES_CAP[0]["smooth_far_km"]) + 0.05,
        "level straight at the 200 m ceiling: %.0f km\nno up-going end surfaces beyond this" % RES_CAP[0]["smooth_far_km"],
        fontsize=6.8, color="0.3", ha="center", va="bottom", zorder=25, bbox=dict(boxstyle="round,pad=0.15", fc="w", ec="none", alpha=.9))

# how far a CO-TILTED ring can reach before it breaks the ceiling.  If the whole ring lies in the tilted plane, its
# depth swings by theta * (Ls + 2R) end to end; with 35 m of cover at the shallow end that must stay under 200 m.
# Tilting only the straights (level arcs + vertical achromats) costs theta * Ls instead, a few tens of metres, so it
# is never the binding constraint -- the ceiling limits the CO-TILTED option, not the beam.
MIN_COVER = 35.0
COTILT = []
for nm, yc, Ls, Rr, colr, da in RINGS:
    ext = Ls + 2 * Rr
    th = (D_CAP - MIN_COVER) / ext                                   # rad
    COTILT.append(dict(ring=nm.split(",")[0], extent_m=ext, max_cotilt_mrad=th * 1e3,
                       max_far_exit_km=cc.exits_smooth(th, D0)[1], swing_at_50mrad_m=50e-3 * ext,
                       swing_at_baseline_m=15.4e-3 * ext, straight_only_swing_at_50mrad_m=50e-3 * Ls))
BIND = min(COTILT, key=lambda c: c["max_far_exit_km"])
ax.plot(AZ_FINE, rmap(np.full_like(AZ_FINE, BIND["max_far_exit_km"])), color="#7a2e21", lw=1.2, ls=(0, (7, 2, 1.5, 2)), zorder=6)
ax.text(math.radians(288), rmap(BIND["max_far_exit_km"]) + 0.045,
        "%s co-tilted hits the ceiling here: %.0f km (%.0f mrad)\nfarther needs level arcs and vertical achromats"
        % (BIND["ring"], BIND["max_far_exit_km"], BIND["max_cotilt_mrad"]),
        fontsize=6.8, color="#7a2e21", ha="center", va="bottom", zorder=25,
        bbox=dict(boxstyle="round,pad=0.15", fc="w", ec="none", alpha=.9))

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
     Line2D([], [], color=c154, lw=3.2, alpha=.35, label="one straight: its two ends are opposite each other"),
     Line2D([], [], color="#0b6e4f", lw=1.3, label="the rings in plan (collider, RCS3/4, RCS1/2)"),
     Line2D([], [], color="0.3", lw=1.3, ls=(0, (4, 2)), label="level straight at the 200 m ceiling: outer limit of every up-going end"),
     Line2D([], [], color="#7a2e21", lw=1.2, ls=(0, (7, 2, 1.5, 2)), label="farthest a co-tilted ring reaches inside the ceiling")]
fig.legend(handles=H, loc="lower center", ncol=3, fontsize=6.8, frameon=False, bbox_to_anchor=(0.5, 0.044),
           handlelength=2.6, columnspacing=1.4)
ax.set_title("Where the two ends of one straight surface, tilt by tilt\n"
             "35 m-deep straights, no ring deeper than 200 m; one colour per tilt $-$ dotted: the end that goes up, solid: the end that goes down; "
             "at 0 mrad they coincide",
             fontsize=9.6, pad=22)
fig.text(0.5, 0.012,
         ("Bold curves: USGS NED 10 m terrain inside 6 km and GEBCO 2020 beyond, along %d bearings (0.1$-$1000 km), Catmull-Rom smoothed; exits over the Great Lakes taken at the water surface. " % len(BRG) if NED_USED else
          "Bold curves: GEBCO 2020 terrain along %d bearings (0.25$-$1000 km), Catmull-Rom smoothed; exits over the Great Lakes taken at the water surface. " % len(BRG)) + 
         "Faint circles: smooth sphere,\n"
         "$s_{near}=R_E(-\\theta+\\sqrt{\\theta^2+2d_0/R_E})\\approx d_0/\\theta$ and "
         "$s_{far}=R_E(\\theta+\\sqrt{\\theta^2+2d_0/R_E})\\approx 2R_E\\theta$. "
         "Log range, 0.2$-$1000 km; rings in km. The fence is 0.9$-$4.9 km out, so the dotted curves show which tilts emerge on site; "
         "the near-field zoom shows how depth up to the 200 m ceiling moves them.",
         ha="center", va="center", fontsize=6.9, color="0.35")
for ext in ("pdf", "svg"):
    fig.savefig(os.path.join(ROOT, "static", "figs", "exit_pairs." + ext), bbox_inches="tight")

# ---------------------------------------------------------------- near field, zoomed: 35 m -> the 200 m ceiling
D2_MIN, D2_MAX = 0.3, 70.0
RINGS2 = (0.5, 1, 2, 5, 10, 20, 50)
LAB_AZ = {0.0: 215, 2.0: 152, 5.0: 232, 15.4: 128, 25.0: 75, 35.0: 108, 50.0: 300, 65.0: 240}
fig2 = plt.figure(figsize=(9.8, 10.3))
fig2.subplots_adjust(top=0.915, bottom=0.13)
ax2 = fig2.add_subplot(111, projection="polar")
rmap2 = cc.polar_axes(ax2, D2_MIN, D2_MAX, rings=RINGS2, tilt_labels=False, fs=6.6, rlabel_pos=180)
ax2.set_yticklabels([])
cc.draw_lakes(ax2, rmap2, cc.lakes_polar(D2_MAX), z=1)
ax2.fill(faz, rmap2(fD), color="0.80", lw=0.8, ec="0.45", zorder=2)
for nm, yc, Ls, Rr, colr, da in RINGS:
    pts = cl.racetrack(yc, Ls, Rr, n_arc=90)
    az_r = np.unwrap(np.array([math.atan2(x, y) for x, y in pts]))
    D_r = np.array([math.hypot(x, y) / 1000.0 for x, y in pts])
    ax2.plot(az_r, rmap2(D_r), color=colr, lw=1.3, alpha=.9, zorder=7)
TICK_AZ = np.radians(np.arange(0, 360, 60) + 20)
for i, (r, rc, col) in enumerate(zip(RES, RES_CAP, COL)):
    lvl = r["tilt_mrad"] == 0
    n35 = smooth(r["far_km"] if lvl else r["near_km"])
    n200 = smooth(rc["far_km"] if lvl else rc["near_km"])
    ls = "-" if lvl else (0, (1.2, 1.6))
    if r["tilt_mrad"] == 15.4:                              # only the baseline's band is filled, or the plot washes out
        ax2.fill_between(AZ_FINE, rmap2(n35), rmap2(n200), color=col, alpha=.11, lw=0, zorder=3)
    for t in TICK_AZ:                                        # a rung of the ladder: same tilt, 35 m -> the ceiling
        j = int(np.argmin(np.abs(AZ_FINE - t)))
        ax2.plot([t, t], [rmap2(n35[j]), rmap2(n200[j])], color=col, lw=0.7, alpha=.33, zorder=4)
    ax2.plot(AZ_FINE, rmap2(n35), color=col, lw=2.0, ls=ls, zorder=6)
    ax2.plot(AZ_FINE, rmap2(n200), color=col, lw=0.9, ls=ls, zorder=6)
    lab2 = "0 mrad, both ends" if lvl else "%g mrad" % r["tilt_mrad"]
    a2 = math.radians(LAB_AZ[r["tilt_mrad"]])
    j = int(np.argmin(np.abs(AZ_FINE - (a2 if a2 <= AZ_FINE[-1] else a2 + 2 * math.pi))))
    ax2.text(a2, 0.5 * (rmap2(n35[j]) + rmap2(n200[j])), lab2, fontsize=7.0, color=col, fontweight="bold",
             ha="center", va="center", zorder=25,
             bbox=dict(boxstyle="round,pad=0.13", fc="w", ec="none", alpha=.88))
for d in RINGS2:                                             # ring labels last, on a bearing the curves cross steeply
    ax2.text(math.radians(186), rmap2(d), "%g km" % d, fontsize=6.6, color="0.35", ha="center", va="center",
             zorder=30, bbox=dict(boxstyle="round,pad=0.15", fc="w", ec="none", alpha=.88))
# the UIUC straight due north: how far depth pushes its up-going end, and the depth at which it reaches the fence
iN = BRG.index(0.0); i15 = TILTS.index(15.4)
fence_N = cc.fence_distance(0.0) / 1000.0
def up_north_at(d):
    return first_exit(15.4e-3, PROF[0.0], RADS[0.0], d)
lo, hi = 35.0, 200.0
d_fence = None
if up_north_at(lo) < fence_N < up_north_at(hi):
    for _ in range(40):
        m = 0.5 * (lo + hi)
        if up_north_at(m) < fence_N:
            lo = m
        else:
            hi = m
    d_fence = 0.5 * (lo + hi)
# the same bisection for every tilt, on the terrain profile: the depth at which the up-going end leaves the fence,
# and how many bearings still emerge on site once the straight sits at the ceiling
def d_at_fence(th_mrad, az):
    b = min(BRG, key=lambda x: abs(((x - az) + 180) % 360 - 180))
    sf = cc.fence_distance(b) / 1000.0
    f = lambda d: first_exit(th_mrad * 1e-3, PROF[b], RADS[b], d)
    if not (f(1.0) < sf < f(D_CAP)):
        return None
    a, c = 1.0, D_CAP
    for _ in range(45):
        m = 0.5 * (a + c)
        a, c = (m, c) if f(m) < sf else (a, m)
    return 0.5 * (a + c)
CEIL_TABLE = []
for r, rc in zip(RES, RES_CAP):
    on35 = sum(1 for k, b in enumerate(BRG) if r["near_km"][k] < cc.fence_distance(b) / 1000.0)
    on200 = sum(1 for k, b in enumerate(BRG) if rc["near_km"][k] < cc.fence_distance(b) / 1000.0)
    CEIL_TABLE.append(dict(tilt_mrad=r["tilt_mrad"], up_north_35m_km=r["near_km"][iN], up_north_ceiling_km=rc["near_km"][iN],
                           depth_at_fence_north_m=d_at_fence(r["tilt_mrad"], 0.0),
                           bearings_on_site_at_35m=on35, bearings_on_site_at_ceiling=on200, n_bearings=len(BRG)))
    print("  %5.1f mrad: up N %6.2f -> %6.2f km, leaves the fence at d0 = %s, on site %3d/%d at 35 m, %3d/%d at the ceiling"
          % (r["tilt_mrad"], r["near_km"][iN], rc["near_km"][iN],
             ("%.0f m" % CEIL_TABLE[-1]["depth_at_fence_north_m"]) if CEIL_TABLE[-1]["depth_at_fence_north_m"] else "n/a",
             on35, len(BRG), on200, len(BRG)))
n35N, n200N = RES[i15]["near_km"][iN], RES_CAP[i15]["near_km"][iN]
ax2.plot([0, 0], [rmap2(n35N), rmap2(n200N)], color=c154, lw=5, alpha=.35, solid_capstyle="butt", zorder=8)
ax2.plot([0], [rmap2(n35N)], "o", ms=8, color=c154, mec="w", mew=1.0, zorder=22)
ax2.plot([0], [rmap2(n200N)], "o", ms=8, color=c154, mec="w", mew=1.0, mfc="w", zorder=22)
ax2.plot([0], [rmap2(fence_N)], "_", ms=14, color="0.2", mew=1.5, zorder=23)
ax2.annotate("the UIUC straight's up-going end, due north:\n%.1f km at 35 m $\\rightarrow$ %.1f km at the 200 m ceiling\n"
             "it leaves the fence (%.2f km) at d$_0$ $\\approx$ %.0f m" % (n35N, n200N, fence_N, d_fence if d_fence else float("nan")),
             (0, rmap2(n200N)), xytext=(math.radians(38), rmap2(D2_MAX) - 0.06), fontsize=7.0, color=c154, fontweight="bold",
             ha="left", va="center", zorder=26, bbox=dict(boxstyle="round,pad=0.15", fc="w", ec=c154, lw=0.6, alpha=.93),
             arrowprops=dict(arrowstyle="-", color=c154, lw=0.6))
H2 = [Line2D([], [], color="0.3", lw=2.0, ls=(0, (1.2, 1.6)), label="up-going end, straight at 35 m"),
      Line2D([], [], color="0.3", lw=0.9, ls=(0, (1.2, 1.6)), label="the same end at the 200 m ceiling"),
      Line2D([], [], color="0.3", lw=0.7, alpha=.33, label="rungs: the depth freedom in between"),
      Line2D([], [], color="0.3", lw=2.0, label="0 mrad: both ends, 35 m and at the ceiling"),
      Patch(color=c154, alpha=.22, label="the baseline tilt's band, filled"),
      Patch(color="0.80", label="Fermilab site"), Line2D([], [], color="#0b6e4f", lw=1.3, label="the rings in plan"),
      Line2D([], [], color=c154, lw=5, alpha=.35, label="the UIUC straight due north, 35 $\\rightarrow$ 200 m")]
fig2.legend(handles=H2, loc="lower center", ncol=3, fontsize=6.8, frameon=False, bbox_to_anchor=(0.5, 0.05),
            handlelength=2.6, columnspacing=1.4)
ax2.set_title("Near field, zoomed: where the up-going end surfaces for straights between 35 m and the 200 m ceiling\n"
              "one colour per tilt; bold at 35 m, thin at the ceiling $-$ $s_{near}\\approx d_0/\\theta$, so extra depth pushes every up-going exit outward",
              fontsize=9.6, pad=22)
fig2.text(0.5, 0.012,
          "Deeper is farther: at the ceiling a level straight surfaces at %.0f km both ways and a 65 mrad straight's up end at %.1f km. "
          "The fence is 0.9$-$4.9 km out, so on-site emergence needs\n$d_0 \\lesssim \\theta\\, s_{fence}$: about %.0f m for the UIUC straight due north; "
          "at 50 mrad it is %.0f m on the tightest bearing and the full 200 m everywhere the fence is past %.1f km. "
          "Terrain: USGS NED 10 m inside 6 km, GEBCO beyond; 144 bearings, Catmull-Rom smoothed."
          % (RES_CAP[0]["smooth_far_km"], RES_CAP[-1]["smooth_near_km"], d_fence if d_fence else float("nan"),
             1e3 * 50e-3 * min(cc.fence_distance(b) / 1000.0 for b in BRG), 0.200 / 50e-3),
          ha="center", va="center", fontsize=6.9, color="0.35")
for ext in ("pdf", "svg"):
    fig2.savefig(os.path.join(ROOT, "static", "figs", "exit_pairs_near." + ext), bbox_inches="tight")

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
for rc in RES_CAP:
    iN_, iS_ = BRG.index(0.0), BRG.index(180.0)
    rc["meridian"] = dict(up_north_end_km=rc["near_km"][iN_], down_south_end_km=rc["far_km"][iS_],
                          up_south_end_km=rc["near_km"][iS_], down_north_end_km=rc["far_km"][iN_])
lev = stats(RES_CAP[0]["far_km"])
for c in COTILT:
    print("co-tilted %-9s extent %5.0f m: ceiling allows %5.1f mrad = far exit %5.0f km; at 15.4 mrad it swings %3.0f m, "
          "at 50 mrad %3.0f m (straights-only at 50 mrad: %.0f m)"
          % (c["ring"], c["extent_m"], c["max_cotilt_mrad"], c["max_far_exit_km"], c["swing_at_baseline_m"],
             c["swing_at_50mrad_m"], c["straight_only_swing_at_50mrad_m"]))
print("ceiling %.0f m: level straight surfaces at %.1f km smooth, %.1f-%.1f km with terrain; UIUC straight's up end %.2f km at 35 m -> %.2f km at the ceiling; "
      "reaches the fence (%.2f km) at d0 = %s m" % (D_CAP, RES_CAP[0]["smooth_far_km"], lev[0], lev[2], n35N, n200N, fence_N, "%.0f" % d_fence if d_fence else "n/a"))
json.dump(dict(straight_depth_m=D0, depth_ceiling_m=D_CAP, ip_elev_m=Z0, bearings_deg=BRG, fence_m=fd, tilts=RES, tilts_at_ceiling=RES_CAP,
               ceiling=dict(level_exit_smooth_km=RES_CAP[0]["smooth_far_km"], level_exit_terrain_min_max_km=[lev[0], lev[2]],
                            uiuc_straight_up_end_km_at_35m=n35N, uiuc_straight_up_end_km_at_ceiling=n200N, fence_north_km=round(fence_N, 3),
                            uiuc_straight_depth_at_fence_m=round(d_fence, 1) if d_fence else None,
                            min_cover_m=MIN_COVER, co_tilted_rings=COTILT, per_tilt=CEIL_TABLE,
                            note="a co-tilted ring's depth swings by theta*(Ls+2R); tilting only the straights costs theta*Ls"),
               terrain=dict(near_field="USGS NED 10 m to %.0f km" % NED_MAX if NED_USED else "GEBCO 2020 only",
                            far_field="GEBCO 2020", seam_offset_gebco_minus_ned_m=round(SEAM_OFFSET_M, 2)),
               worked_pairs=dict(baseline_15p4_up_north=dict(far_south_km=Df_b, near_north_km=Dn_b),
                                 rcs_50_to_lake_huron_az52=dict(far_km=Df_h, near_km=Dn_h))),
          open(os.path.join(ROOT, "static", "geo", "exit_pairs.json"), "w"), indent=1)
print("baseline pair: far S %.1f km, near N %.2f km; Huron pair: far %.1f, near %.2f" % (Df_b, Dn_b, Df_h, Dn_h))
