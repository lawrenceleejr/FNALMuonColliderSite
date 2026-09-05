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

def first_exit(theta_end, hg, rad):
    """First range (km) where the beam end (signed slope theta_end, rad) is at or above ground."""
    diff = -D0 + theta_end * rad * 1000.0 - hg
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

def stats(v):
    v = [x for x in v if x is not None]
    return (min(v), float(np.median(v)), max(v)) if v else (None, None, None)

# ---------------------------------------------------------------- figure
D_MIN, D_MAX = 0.1, 1500.0
fig = plt.figure(figsize=(9.8, 10.6))
fig.subplots_adjust(top=0.915, bottom=0.135)
ax = fig.add_subplot(111, projection="polar")
rmap = cc.polar_axes(ax, D_MIN, D_MAX, rings=(0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000),
                     tilt_labels=False, fs=6.6, rlabel_pos=197)
ax.set_yticklabels([])
for d in (0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000):        # ring labels, drawn on top
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

# one straight with BOTH ends far: UIUC (S) and Lake Superior (N) -- a chord 10.2 km below the IP
S1, S2 = 198.4, 655.0
d_chord = S1 * S2 / (2 * cc.R_E)
gcol = "0.35"
ax.plot([math.pi, math.pi], [0, rmap(S1)], color=gcol, lw=1.2, ls=(0, (3, 2)), alpha=.8, zorder=4)
ax.plot([0, 0], [0, rmap(S2)], color=gcol, lw=1.2, ls=(0, (3, 2)), alpha=.8, zorder=4)
ax.plot([0], [rmap(S2)], "s", ms=7, color=gcol, mec="w", mew=1.0, zorder=22)
ax.annotate("one straight, BOTH ends far:\nUIUC 198 km S and Lake Superior 655 km N\nneeds the straight %.1f km below the IP\n"
            "($d_0 = s_1 s_2 / 2R_E$; deepest point 14 km)" % d_chord, (0, rmap(S2)), xytext=(math.radians(24), rmap(1050.0)),
            fontsize=6.8, color=gcol, ha="left", va="center", zorder=26,
            bbox=dict(boxstyle="round,pad=0.15", fc="w", ec=gcol, lw=0.6, alpha=.92),
            arrowprops=dict(arrowstyle="-", color=gcol, lw=0.6))

# the chords through UIUC and a northern lake (two_ends.py): the true great-circle chord from the South Farms (dotted, misses
# the IP), and for Green Bay also the same two ranges as a straight THROUGH the IP (dashed spokes)
TEJ = json.load(open(os.path.join(ROOT, "static", "geo", "two_ends.json")))["chord_for_plots"]
for k, CH in enumerate(TEJ["chords"][:2]):
    ccol = CH["color"]; pt = CH["point"]; tcx = CH["true_chord"]; tip = CH["through_ip"]
    az_c, D_c, dep_c = cc.great_circle_polar(TEJ["uiuc"], (pt["lat"], pt["lon"]))
    ax.plot(az_c, rmap(D_c), color=ccol, lw=1.6, ls=(0, (1.5, 1.5)), zorder=7)
    ax.plot([az_c[-1]], [rmap(D_c[-1])], "s", ms=7, color=ccol, mec="w", mew=1.0, zorder=22)
    i_min = int(np.argmin(D_c))
    ax.plot([az_c[i_min]], [rmap(D_c[i_min])], "o", ms=5, color=ccol, mec="w", zorder=22)
    if CH["name"] == "Green Bay":
        tg, tu = math.radians(pt["bearing_deg"]), math.radians(tip["uiuc_end_bearing_deg"])
        ax.plot([tu, tu], [0, rmap(S1)], color=ccol, lw=1.2, ls=(0, (3, 2)), alpha=.85, zorder=4)
        ax.plot([tg, tg], [0, rmap(pt["range_km"])], color=ccol, lw=1.2, ls=(0, (3, 2)), alpha=.85, zorder=4)
        ax.plot([tg], [rmap(pt["range_km"])], "s", ms=7, color=ccol, mec="w", mew=1.0, zorder=22)
        ax.plot([tu], [rmap(S1)], "s", ms=7, color=ccol, mec="w", mew=1.0, zorder=22)
        txt = ("chord through UIUC and Green Bay\n"
               "dotted: the true chord, South Farms to mid-bay (%.0f km at %.0f$^\\circ$):\n"
               "it passes %.0f km east of the IP, %.1f km deep there, %.1f km at midpoint\n"
               "dashed: the same two ranges as a straight through the IP:\n"
               "%.1f km deep, UIUC end %.0f km west of the South Farms"
               % (pt["range_km"], pt["bearing_deg"], tcx["closest_approach_to_ip_km"], tcx["depth_at_closest_approach_km"],
                  tcx["midpoint_depth_km"], tip["depth_at_ip_km"], tip["uiuc_end_offset_km"]))
        xy, ha = (0.985, 0.162), "right"
    else:
        txt = ("chord through UIUC and Lake Winnebago (dotted)\n"
               "South Farms to the lake's water (%.0f km at %.0f$^\\circ$): it crosses the\n"
               "Fermilab site %.1f km WEST of the IP, %.1f km deep there (%.1f km at midpoint)\n"
               "as a straight through the IP it is %.1f km deep, tilt %.1f mrad,\n"
               "with its UIUC end %.0f km east of the South Farms"
               % (pt["range_km"], pt["bearing_deg"], tcx["closest_approach_to_ip_km"], tcx["depth_at_closest_approach_km"],
                  tcx["midpoint_depth_km"], tip["depth_at_ip_km"], tip["tilt_mrad"], abs(tip["uiuc_end_offset_km"])))
        xy, ha = (0.015, 0.162), "left"
    ax.annotate(txt, (az_c[i_min], rmap(D_c[i_min])), xytext=xy, textcoords="figure fraction", fontsize=6.5, color=ccol,
                ha=ha, va="center", zorder=26, bbox=dict(boxstyle="round,pad=0.2", fc="w", ec=ccol, lw=0.6, alpha=.93),
                arrowprops=dict(arrowstyle="-", color=ccol, lw=0.6, alpha=.7))

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
     Line2D([], [], color="0.35", lw=1.2, ls=(0, (3, 2)), label="a chord with both ends far (10 km deep)"),
     Line2D([], [], color="#2a7f9e", lw=1.6, ls=(0, (1.5, 1.5)), label="the UIUC$-$Green Bay chord (true, misses the IP)"),
     Line2D([], [], color="#2a7f9e", lw=1.2, ls=(0, (3, 2)), label="same ranges as a straight through the IP"),
     Line2D([], [], color="#6b8e23", lw=1.6, ls=(0, (1.5, 1.5)), label="the UIUC$-$Lake Winnebago chord (true, crosses the site)")]
fig.legend(handles=H, loc="lower center", ncol=3, fontsize=6.8, frameon=False, bbox_to_anchor=(0.5, 0.044),
           handlelength=2.6, columnspacing=1.4)
ax.set_title("Where the two ends of one straight surface, tilt by tilt\n"
             "35 m-deep straight; one colour per tilt $-$ dotted: the end that goes up, solid: the end that goes down; "
             "at 0 mrad they are the same curve",
             fontsize=9.6, pad=22)
fig.text(0.5, 0.012,
         ("Bold curves: USGS NED 10 m terrain inside 6 km and GEBCO 2020 beyond, along %d bearings (0.1$-$1000 km), Catmull-Rom smoothed; exits over the Great Lakes taken at the water surface. " % len(BRG) if NED_USED else
          "Bold curves: GEBCO 2020 terrain along %d bearings (0.25$-$1000 km), Catmull-Rom smoothed; exits over the Great Lakes taken at the water surface. " % len(BRG)) + 
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
               terrain=dict(near_field="USGS NED 10 m to %.0f km" % NED_MAX if NED_USED else "GEBCO 2020 only",
                            far_field="GEBCO 2020", seam_offset_gebco_minus_ned_m=round(SEAM_OFFSET_M, 2)),
               worked_pairs=dict(baseline_15p4_up_north=dict(far_south_km=Df_b, near_north_km=Dn_b),
                                 rcs_50_to_lake_huron_az52=dict(far_km=Df_h, near_km=Dn_h))),
          open(os.path.join(ROOT, "static", "geo", "exit_pairs.json"), "w"), indent=1)
print("baseline pair: far S %.1f km, near N %.2f km; Huron pair: far %.1f, near %.2f" % (Df_b, Dn_b, Df_h, Dn_h))
