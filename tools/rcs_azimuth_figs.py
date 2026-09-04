#!/usr/bin/env python3
"""Polar opportunity map for the azimuth-freedom study (tools/rcs_azimuth.py).

Azimuth around, required tilt radial.  Each water body occupies an arc-band;
the site's own geometry sets theta_free(bearing) -- the tilt above which the
NEAR beam is in navigable airspace before it leaves the fence.  A target is
"free" (needs no corridor, no easement, no overflight) when its band lies
outside that curve.
"""
import json, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
STYLE = "/root/.claude/skills/synced/tufte-plots/assets/tufte.mplstyle"
if os.path.exists(STYLE): plt.style.use(STYLE)
plt.rcParams["figure.dpi"] = 130
import rcs_azimuth as ra

D = json.load(open(os.path.join(ROOT, "static", "geo", "rcs_azimuth.json")))
OPP, CAND = D["opportunity_map"], D["verified_candidates"]
VER = {c["body"]: c for c in CAND if c["is_water"]}

fig = plt.figure(figsize=(7.2, 6.4))
ax = fig.add_subplot(111, projection="polar")
ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)

# theta_free(far azimuth): uses the fence distance on the OPPOSITE bearing
az = np.arange(0, 360, 2.0)
tf = np.array([(15.0 + 500 * 0.3048) / ra.fence_distance((a + 180) % 360) * 1e3
               for a in az])
th = np.radians(az)
ax.plot(np.append(th, th[0]), np.append(tf, tf[0]), color="#0b6e4f", lw=1.8, zorder=5)
ax.fill_between(np.append(th, th[0]), np.append(tf, tf[0]), 120,
                color="#0b6e4f", alpha=.09, zorder=1)

COL = {"Lake Huron": "#b5541c", "Lake Superior": "#6d3580",
       "Lake Ontario": "#1c5a96", "Lake Erie": "#8a6d1f",
       "Lake Michigan": "#3d7ea6"}
for nm, w in OPP.items():
    if nm not in COL:
        continue
    v = VER.get(nm)
    a0 = math.radians(w["az_deg"])
    t0, t1 = w["tilt_mrad"]
    ax.plot([a0, a0], [t0, t1], color=COL[nm], lw=6, alpha=.85, zorder=4,
            solid_capstyle="butt")
    if v:
        ax.plot([math.radians(v["azimuth_deg"])], [v["tilt_mrad"]], "o",
                color=COL[nm], ms=8, mec="w", mew=1.2, zorder=6)

lab = [("Lake Huron", 50, 49.5, 93, 62,
        "LAKE HURON $-$ best overall\n630 km, 49.5 mrad\n179 m water, FREE"),
       ("Lake Ontario", 74, 67.2, 104, 92,
        "Lake Ontario\n856 km, 67.2 mrad, 176 m\nFREE (steeper achromat)"),
       ("Lake Erie", 82, 49.3, 34, 108,
        "Lake Erie: free,\nbut only 44 m of water"),
       ("Lake Superior", 8, 51.4, 97, 340,
        "Lake Superior\n654 km, 319 m water\n$-$ deepest, but NOT free\n(0.58 km off-site)"),
       ("Lake Michigan", 24, 24.4, 30, 300,
        "Lake Michigan: too close\n310 km $\\Rightarrow$ 4$\\times$ dose,\n3.8 km overflight")]
for nm, a, tt, rr, aa, txt in lab:
    ax.annotate(txt, (math.radians(a), tt), xytext=(math.radians(aa), rr),
                fontsize=7.4, color=COL[nm], ha="center", va="center",
                arrowprops=dict(arrowstyle="-", color=COL[nm], lw=0.7, alpha=.75))
ax.text(math.radians(183), 24, "corridor-bound", fontsize=8.5, color="0.45",
        ha="center", style="italic")
ax.text(math.radians(120), 110, "azimuth FREE", fontsize=8.5, color="#0b6e4f",
        ha="center", style="italic")
ax.set_ylim(0, 118)
ax.set_yticks([20, 40, 60, 80, 100])
ax.set_yticklabels(["20", "40", "60", "80", "100 mrad"], fontsize=7.5, color="0.4")
ax.set_xticks(np.radians(np.arange(0, 360, 45)))
ax.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"], fontsize=9)
ax.grid(color="0.85", lw=0.6)
ax.set_title("Aiming the last RCS: required tilt vs bearing, against the site's own\n"
             "$\\theta_{free}$ — inside the green curve you need a corridor, outside you don't",
             fontsize=9.5, pad=22)
fig.text(0.5, 0.045,
    "$\\theta_{free}$ is set by the site's own geometry: above it the near beam reaches navigable airspace "
    "BEFORE leaving the fence\n$-$ no corridor, no easement, no overflight. It dips to ~35 mrad toward ENE "
    "because the IP sits on the site's EAST edge:\n4.7 km of Fermilab land lies WSW, only 0.9 km east. "
    "Bands are GEBCO-verified submerged runs; dots are the verified aim points.",
    ha="center", va="center", fontsize=7.3, color="0.35")
for ext in ("pdf", "svg"):
    fig.savefig(os.path.join(ROOT, "static", "figs", "rcs_azimuth_map." + ext),
                bbox_inches="tight")
print("wrote rcs_azimuth_map")
