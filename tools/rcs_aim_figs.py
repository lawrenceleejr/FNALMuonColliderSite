#!/usr/bin/env python3
"""Figures for the RCS aiming study (tools/rcs_aim.py)."""
import json, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
STYLE = "/root/.claude/skills/synced/tufte-plots/assets/tufte.mplstyle"
if os.path.exists(STYLE): plt.style.use(STYLE)
plt.rcParams["figure.dpi"] = 130
FIGS = os.path.join(ROOT, "static", "figs")

def save(fig, name):
    for ext in ("pdf", "svg"):
        fig.savefig(os.path.join(FIGS, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig); print("wrote", name)

R_E = 6371.0            # km
IP_LAT = 41.8443
MLAT = 111.132          # km/deg
G = json.load(open(os.path.join(ROOT, "data", "gebco_meridian.json")))["profile"]
D = json.load(open(os.path.join(ROOT, "static", "geo", "rcs_aim.json")))
lat = np.array([a for a, _ in G]); el = np.array([e for _, e in G], dtype=float)
y = (lat - IP_LAT) * MLAT                       # km north of the IP

# ---------------------------------------------------------------- Figure A
# The two chords, on the real meridian section: UIUC vs Lake Superior.
fig, ax = plt.subplots(figsize=(7.0, 4.1))
ax.fill_between(y, el / 1000.0, -12, color="0.90", lw=0)
ax.plot(y, el / 1000.0, color="0.45", lw=0.9)
ax.fill_between(y, 0, el / 1000.0, where=el < 0, color="#3d7ea6", alpha=.55, lw=0)
ax.axhline(0, color="0.6", lw=0.6, ls=(0, (4, 3)))

def chord(depth_m, th_mrad, col, lab, yl):
    th = th_mrad * 1e-3
    z0 = (np.interp(IP_LAT, lat, el) - depth_m) / 1000.0
    zz = z0 + th * yl + yl ** 2 / (2 * R_E)
    return zz

for s, col in ((D["scenarios"][0], "#0b6e4f"), (D["scenarios"][2], "#b5541c")):
    far = s[s["far_side"]]
    y0 = 0.0
    yf = (far["exit_lat"] - IP_LAT) * MLAT
    yl = np.linspace(min(0, yf), max(0, yf), 900)
    zz = chord(s["depth_m"], s["theta_mrad"], col, s["label"], yl)
    ax.plot(yl, zz, color=col, lw=1.9)
    ax.plot([yf], [np.interp(far["exit_lat"], lat, el) / 1000.0], "o",
            color=col, ms=6, zorder=6)
    per = far.get("perigee")
    if per:
        yp = (per["lat"] - IP_LAT) * MLAT
        ax.annotate("perigee %.1f km" % per["depth_km"],
                    (yp, chord(s["depth_m"], s["theta_mrad"], col, "", np.array([yp]))[0]),
                    xytext=(yp, -11.1), fontsize=7.5, color=col, ha="center",
                    arrowprops=dict(arrowstyle="-", color=col, lw=0.7))
ax.annotate("UIUC South Farms\n198 km, 15.4 mrad\n(land: fenced farm strip)",
            (-197.7, 0.22), xytext=(-330, -3.4), fontsize=8, color="#0b6e4f",
            ha="center", arrowprops=dict(arrowstyle="->", color="#0b6e4f", lw=0.8))
ax.annotate("Lake Superior open water\n656 km, 51.8 mrad\n"
            "~300 m water column,\n30+ km offshore",
            (655.7, 0.05), xytext=(690, -4.6), fontsize=8, color="#b5541c",
            ha="center", arrowprops=dict(arrowstyle="->", color="#b5541c", lw=0.8))
ax.text(505, 0.95, "Keweenaw Peninsula\n(+416 m: NOT water)", fontsize=7.2,
        color="0.25", ha="center")
ax.annotate("", (600, 0.45), xytext=(530, 0.78),
            arrowprops=dict(arrowstyle="->", color="0.45", lw=0.7))
ax.text(-1230, -0.9, "Gulf of Mexico\n1280 km needs 101 mrad", fontsize=7.2,
        color="0.35", ha="center")
ax.set_xlim(-1330, 830)
ax.set_ylim(-11.8, 1.9)
ax.set_xlabel("km along the corridor meridian  (south $\\leftarrow$ | $\\rightarrow$ north)")
ax.set_ylabel("elevation  [km]")
ax.set_title("Where the last RCS can point: the two viable targets on one meridian",
             fontsize=10, loc="left")
save(fig, "rcs_aim_chords")

# ---------------------------------------------------------------- Figure B
# depth and angle are independent knobs
fig, ax = plt.subplots(figsize=(6.6, 3.8))
d0 = np.linspace(10, 460, 200)
th = -51.80e-3
far = np.array([656.0 + 0.0 * v for v in d0])       # far range: flat in depth
near = d0 / abs(th) / 1000.0
ax.plot(d0, near, color="#b5541c", lw=1.9)
ax.plot(d0, [655.4 + (662.9 - 655.4) * (v - 15) / (440 - 15) for v in d0],
        color="0.45", lw=1.9)
ax.set_yscale("log")
ax.set_xlabel("straight depth below grade  [m]   (the other knob)")
ax.set_ylabel("surfacing distance  [km]")
ax.text(235, 170, "FAR exit (north, into the lake):\n655 $\\to$ 663 km across the whole "
        "depth range\n$-$ set by the ANGLE alone", fontsize=8, color="0.3",
        ha="center")
ax.text(255, 1.05, "NEAR exit (south): $s \\simeq d_0/|\\theta|$\n$-$ set by the DEPTH alone",
        fontsize=8, color="#b5541c", ha="center")
ax.axhspan(0.1, 2.632, color="#0b6e4f", alpha=.10, lw=0)
ax.text(462, 0.30, "on Fermilab land", fontsize=7.5, color="#0b6e4f", ha="right")
for v, lb in ((15, "15 m"), (120, "120 m"), (440, "440 m")):
    ax.plot([v], [v / abs(th) / 1000.0], "o", color="#b5541c", ms=5)
    ax.annotate(lb, (v, v / abs(th) / 1000.0), xytext=(v, v / abs(th) / 1000.0 * 1.9),
                fontsize=7.2, color="#b5541c", ha="center")
ax.set_xlim(0, 470); ax.set_ylim(0.15, 1500)
ax.set_title("Depth and angle decouple: angle picks the target, depth picks the "
             "near emergence", fontsize=10, loc="left")
save(fig, "rcs_aim_knobs")
