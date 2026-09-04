#!/usr/bin/env python3
"""Polar landmark map for the collider in its default civic-envelope aim.

Same frame as the RCS opportunity map (bearing around, range radial) but for
the ring whose tilt is FIXED: 15.40 mrad up-north puts the south beam at
UIUC (~198 km) and sends the north beam climbing over the ComEd corridor.
A fixed tilt sweeps a fixed-range circle around the compass, so the plot
shows what else sits on -- or near -- that circle, and everything a
neutrino-beam siting study wants to know the bearing and range of: the
Great Lakes (real footprints, projected), DOE and university partners,
mines and quarries (candidate far halls), public land, airports, and the
population centres the beam must not graze.

Radius is log(range) so that 10 km neighbours and 1300 km mines share one
frame; every range ring is also labelled with the tilt that reaches it,
theta = D / (2 R_E).  The site's own theta_free(bearing) curve is drawn
too: the collider's 198 km circle lies entirely inside it, which is exactly
why this ring is corridor-bound in every bearing and why the ComEd
right-of-way was chosen.

Output: static/figs/collider_map.{svg,pdf}; landmarks in static/geo/collider_map.json.
"""
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
STYLE = "/root/.claude/skills/synced/tufte-plots/assets/tufte.mplstyle"
if os.path.exists(STYLE):
    plt.style.use(STYLE)
plt.rcParams["figure.dpi"] = 130

R_E = 6371.0088
IP_LAT, IP_LON = 41.8443, -88.222972
TH_COLL = 15.40e-3
D0_KM = 0.035                         # km, straight depth
D_COLL = R_E * (TH_COLL + math.sqrt(TH_COLL**2 + 2 * D0_KM / R_E))   # km, far-exit circle (198.4)

def inv(lat, lon):
    p1, p2 = math.radians(IP_LAT), math.radians(lat)
    dl = math.radians(lon - IP_LON)
    h = math.sin((p2 - p1) / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    D = 2 * R_E * math.asin(min(1.0, math.sqrt(h)))
    az = math.degrees(math.atan2(math.sin(dl) * math.cos(p2),
                                 math.cos(p1) * math.sin(p2)
                                 - math.sin(p1) * math.cos(p2) * math.cos(dl))) % 360
    return D, az

# ----------------------------------------------------------------------
# Site boundary -> theta_free(bearing)
# ----------------------------------------------------------------------
_b = json.load(open(os.path.join(ROOT, "data", "fnal_boundary.geojson")))
MLAT = 111.132
MLON = 111.320 * math.cos(math.radians(IP_LAT))
BXY = [((lo - IP_LON) * MLON * 1000, (la - IP_LAT) * MLAT * 1000)
       for lo, la in _b["features"][0]["geometry"]["coordinates"][0]]

def _pin(x, y, poly):
    ins = False
    for (x1, y1), (x2, y2) in zip(poly, poly[1:] + [poly[0]]):
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1:
            ins = not ins
    return ins

def fence_distance(az_deg):
    a = math.radians(az_deg)
    ux, uy = math.sin(a), math.cos(a)
    prev, d = 0.0, 25.0
    while d <= 25000:
        if not _pin(ux * d, uy * d, BXY):
            lo, hi = prev, d
            for _ in range(50):
                m = 0.5 * (lo + hi)
                if _pin(ux * m, uy * m, BXY):
                    lo = m
                else:
                    hi = m
            return 0.5 * (lo + hi)
        prev, d = d, d + 25.0
    return 25000.0

D0 = 35.0                              # m, the collider straight's depth
FT = 0.3048

# ----------------------------------------------------------------------
# Landmarks (lat, lon). Coordinates are gazetteer-grade (~0.01 deg), which
# is far finer than this log-range plot can show.
# ----------------------------------------------------------------------
LANDMARKS = {
    "university": [
        ("UIUC South Farms (the aim)", 40.0602, -88.2230),
        ("Purdue", 40.4237, -86.9212),
        ("Notre Dame", 41.7056, -86.2353),
        ("UW-Madison", 43.0766, -89.4125),
        ("U. Iowa", 41.6611, -91.5302),
        ("Michigan State", 42.7018, -84.4822),
        ("Indiana U.", 39.1754, -86.5127),
    ],
    "federal": [
        ("Argonne NL", 41.7183, -87.9786),
        ("Midewin Natl Tallgrass Prairie", 41.36, -88.13),
    ],
    "mine": [
        ("Soudan mine (MINOS far)", 47.8232, -92.2373),
        ("Ash River (NOvA far)", 48.3794, -92.8313),
        ("SURF / Homestake (DUNE far)", 44.3525, -103.7510),
        ("Thornton Quarry", 41.573, -87.611),
        ("Ottawa silica mines", 41.36, -88.85),
        ("S. Illinois coal basin", 37.75, -88.75),
        ("Iron Range (Hibbing)", 47.427, -92.938),
        ("Iron Mountain (Menominee Range)", 45.82, -88.07),
        ("Marquette Range (Empire/Tilden)", 46.45, -87.60),
    ],
    "storage": [
        ("Manlove gas storage (Mt. Simon)", 40.40, -88.20),
        ("Herscher gas storage (Mt. Simon)", 41.05, -88.10),
    ],
    "public": [
        ("Kettle Moraine SF (N)", 43.60, -88.20),
        ("Shawnee NF", 37.45, -88.65),
        ("Chequamegon-Nicolet NF", 45.75, -88.85),
        ("Ottawa NF", 46.35, -89.30),
        ("Hiawatha NF", 46.05, -86.75),
        ("Starved Rock SP", 41.319, -88.994),
        ("Moraine Hills SP", 42.30, -88.23),
        ("Chain O'Lakes SP", 42.46, -88.20),
        ("Illinois Beach SP", 42.43, -87.81),
        ("Indiana Dunes NP", 41.65, -87.05),
    ],
    "city": [
        ("Chicago Loop", 41.8781, -87.6298),
        ("Milwaukee", 43.0389, -87.9065),
        ("Rockford", 42.2711, -89.0940),
        ("Peoria", 40.6936, -89.5890),
        ("Springfield IL", 39.7817, -89.6501),
        ("Indianapolis", 39.7684, -86.1581),
        ("Bloomington-Normal", 40.4842, -88.9937),
        ("St. Louis", 38.6270, -90.1994),
        ("Aurora", 41.7606, -88.3201),
        ("Naperville", 41.7508, -88.1535),
        ("Elgin", 42.0354, -88.2826),
    ],
    "airport": [
        ("O'Hare", 41.9742, -87.9073),
        ("Midway", 41.7868, -87.7522),
        ("Willard (UIUC, IP2 aim)", 40.0392, -88.2781),
    ],
    "corridor": [
        ("Smith Rd / ComEd ROW\n(houses north of the lab)", 41.9356, -88.2230),
        ("Wayne (384-500 ft overflight)", 41.953, -88.2230),
        ("Cary / Fox River Grove\n(untilted N exit, 38 km)", 42.1884, -88.2230),
        ("Wheatland Twp (untilted S exit, 28 km)", 41.5930, -88.2230),
    ],
}
STYLE_BY = {
    "university": dict(marker="*", ms=11, color="#0b6e4f", z=7),
    "federal":    dict(marker="s", ms=6.5, color="#1c5a96", z=6),
    "mine":       dict(marker="D", ms=6, color="#7a2e21", z=6),
    "storage":    dict(marker="h", ms=7, color="#7a2e21", z=6),
    "public":     dict(marker="^", ms=7, color="#4d7a2f", z=5),
    "city":       dict(marker="o", ms=4.5, color="0.45", z=4),
    "airport":    dict(marker="x", ms=6, color="0.35", z=4),
    "corridor":   dict(marker="|", ms=9, color="#b5541c", z=6),
}

# ----------------------------------------------------------------------
# Lakes: real footprints projected into (bearing, range)
# ----------------------------------------------------------------------
_lk = json.load(open(os.path.join(ROOT, "data", "ne_10m_lakes.geojson")))
GREAT = {"Lake Superior", "Lake Michigan", "Lake Huron", "Lake Erie",
         "Lake Ontario", "Lake Winnebago", "Lake Winnipeg", "Lake Nipigon",
         "Lake of the Woods", "Lake Simcoe", "Lake St. Clair", "Lake Geneva"}

def rings(g):
    if g["type"] == "Polygon":
        return [g["coordinates"][0]]
    if g["type"] == "MultiPolygon":
        return [p[0] for p in g["coordinates"]]
    return []


# ----------------------------------------------------------------------
# What the 198 km circle lands on, per degree of bearing, and what each
# bearing would cost on the near side (sub-500 ft band beyond the fence).
# ----------------------------------------------------------------------
def fwd(az_deg, D):
    p1, l1 = math.radians(IP_LAT), math.radians(IP_LON)
    a, dl = math.radians(az_deg), D / R_E
    p2 = math.asin(math.sin(p1) * math.cos(dl) + math.cos(p1) * math.sin(dl) * math.cos(a))
    l2 = l1 + math.atan2(math.sin(a) * math.sin(dl) * math.cos(p1),
                         math.cos(dl) - math.sin(p1) * math.sin(p2))
    return math.degrees(p2), math.degrees(l2)

_near = []
for f in _lk["features"]:
    nm = f["properties"].get("name") or "unnamed lake"
    for r in rings(f["geometry"]):
        if min(inv(c[1], c[0])[0] for c in r[::5]) < D_COLL + 150:
            _near.append((nm, [(c[0], c[1]) for c in r]))

def lake_at(lat, lon):
    for nm, poly in _near:
        if _pin(lon, lat, poly):
            return nm
    return None

CIRCLE = []
for a in range(360):
    la, lo = fwd(a, D_COLL)
    sf = fence_distance((a + 180) % 360)
    CIRCLE.append(dict(az=a, lat=round(la, 4), lon=round(lo, 4), water=lake_at(la, lo),
                       near_fence_m=round(sf), offsite_sub500ft_km=round(((D0 + 500 * FT) / TH_COLL - sf) / 1000, 2)))
arcs, cur = [], None
for c in CIRCLE + [dict(az=360, water=None)]:
    if c["water"] and cur and cur["water"] == c["water"] and c["az"] == cur["end"] + 1:
        cur["end"] = c["az"]
    else:
        if cur:
            arcs.append(cur)
        cur = dict(water=c["water"], start=c["az"], end=c["az"]) if c["water"] else None
WATER_ARCS = [dict(body=a["water"], az_from=a["start"], az_to=a["end"], deg=a["end"] - a["start"] + 1) for a in arcs]

# ----------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------
D_MIN, D_MAX = 8.0, 1500.0
rmap = lambda D: np.log10(np.maximum(np.asarray(D, float), D_MIN) / D_MIN)
R_MAX = float(rmap(D_MAX))

fig = plt.figure(figsize=(9.4, 10.2))
fig.subplots_adjust(top=0.93, bottom=0.115)
ax = fig.add_subplot(111, projection="polar")
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.set_ylim(0, R_MAX)

# theta_free(bearing) -> D_free; shade the free region (outside)
az = np.arange(0, 361, 2.0)
tfree = np.array([(D0 + 500 * FT) / fence_distance((a + 180) % 360) for a in az])
Dfree = 2 * R_E * tfree
th = np.radians(az)
ax.fill_between(th, rmap(Dfree), R_MAX, color="#0b6e4f", alpha=.06, lw=0, zorder=0)
ax.plot(th, rmap(Dfree), color="#0b6e4f", lw=1.3, zorder=3)

# lakes
for f in _lk["features"]:
    nm = f["properties"].get("name")
    for r in rings(f["geometry"]):
        r = r[::max(1, len(r) // 400)]          # decimate: invisible at this scale, keeps the SVG small
        pts = [inv(c[1], c[0]) for c in r]
        Ds = np.array([p[0] for p in pts])
        if Ds.min() > D_MAX or len(pts) < 12:
            continue
        azs = np.unwrap(np.radians([p[1] for p in pts]))
        col = "#3d7ea6" if nm in GREAT else "#7fb0cc"
        ax.fill(azs, rmap(Ds), color=col, alpha=.55 if nm in GREAT else .4,
                lw=0, zorder=2)

# the collider's fixed-tilt circle
ax.plot(np.radians(np.arange(0, 361)), np.full(361, rmap(D_COLL)),
        color="#0b6e4f", lw=2.0, ls=(0, (5, 3)), zorder=6)

# the two beams: south chord underground to UIUC; north near-beam's sub-500 ft band
ax.plot([math.pi, math.pi], [0, rmap(D_COLL)], color="#0b6e4f", lw=1.1, alpha=.7, zorder=5)
ax.plot([0, 0], [0, rmap((D0 + 500 * FT) / TH_COLL / 1000)], color="#b5541c", lw=3.5,
        solid_capstyle="butt", zorder=5)

# landmarks
OUT = []
for cat, items in LANDMARKS.items():
    st = STYLE_BY[cat]
    for name, la, lo in items:
        D, a = inv(la, lo)
        if D > D_MAX:
            continue
        hollow = cat == "storage"
        ax.plot([math.radians(a)], [rmap(D)], st["marker"], ms=st["ms"],
                color="w" if hollow else st["color"],
                mec=st["color"] if (hollow or cat == "corridor") else "w",
                mew=1.2 if hollow else 0.8, zorder=st["z"])
        OUT.append(dict(name=name.replace("\n", " "), category=cat,
                        lat=la, lon=lo, range_km=round(D, 1),
                        bearing_deg=round(a, 1),
                        tilt_to_reach_mrad=round((D / (2 * R_E) - D0_KM / D) * 1e3, 2),
                        on_collider_circle=bool(abs(D - D_COLL) < 15)))

# label placement: (name, az, D, dr (log units), da (deg), ha)
def label(name, la, lo, dr=0.0, da=0.0, ha="left", col="0.25", fs=7.0, bold=False):
    D, a = inv(la, lo)
    ax.annotate(name, (math.radians(a), rmap(D)),
                xytext=(math.radians(a + da), rmap(D) + dr),
                fontsize=fs, color=col, ha=ha, va="center",
                fontweight="bold" if bold else "normal",
                bbox=dict(boxstyle="round,pad=0.12", fc="w", ec="none", alpha=0.72),
                arrowprops=dict(arrowstyle="-", color=col, lw=0.5, alpha=.6)
                if (abs(dr) > 0.04 or abs(da) > 3) else None, zorder=9)

G, F, M, P, C, A, K = ("#0b6e4f", "#1c5a96", "#7a2e21", "#4d7a2f", "0.45", "0.35", "#b5541c")
# south sector
label("UIUC South Farms\nthe collider's aim (198 km)", 40.0602, -88.2230, dr=0.14, da=9, col=G, fs=7.6, bold=True)
label("Purdue (192 km: on the circle)", 40.4237, -86.9212, dr=0.02, da=-4, col=G)
label("Indiana U.", 39.1754, -86.5127, dr=0.05, col=G)
label("Indianapolis", 39.7684, -86.1581, dr=0.05, col=C)
label("Bloomington-Normal", 40.4842, -88.9937, dr=-0.06, da=-4, ha="right", col=C)
label("Manlove gas storage", 40.40, -88.20, dr=-0.04, da=-6, col=M, fs=6.5)
label("Herscher gas storage", 41.05, -88.10, dr=0.0, da=8, ha="right", col=M, fs=6.5)
label("Peoria", 40.6936, -89.5890, dr=0.05, ha="right", col=C)
label("Springfield", 39.7817, -89.6501, dr=0.05, ha="right", col=C)
label("St. Louis", 38.6270, -90.1994, dr=0.05, ha="right", col=C)
label("S. Illinois coal basin", 37.75, -88.75, dr=0.0, da=6, ha="right", col=M)
label("Shawnee NF", 37.45, -88.65, dr=0.02, da=-5, col=P)
label("Midewin (federal)", 41.36, -88.13, dr=0.04, da=-8, col=F)
label("untilted S exit, 28 km", 41.5930, -88.2230, dr=0.0, da=6, ha="right", col=K, fs=6.5)
label("Aurora", 41.7606, -88.3201, dr=0.22, da=25, ha="right", col=C, fs=6.5)
label("Naperville", 41.7508, -88.1535, dr=0.02, da=-10, col=C, fs=6.5)
label("Starved Rock SP", 41.319, -88.994, dr=-0.10, da=-4, ha="right", col=P)
label("Ottawa silica mines", 41.36, -88.85, dr=0.10, da=-2, ha="right", col=M)
# east sector
label("Argonne NL", 41.7183, -87.9786, dr=0.16, da=6, col=F)
label("Thornton Quarry", 41.573, -87.611, dr=0.06, da=6, col=M)
label("Indiana Dunes NP", 41.65, -87.05, dr=0.05, da=6, col=P)
label("Notre Dame", 41.7056, -86.2353, dr=0.05, da=4, col=G)
label("Chicago", 41.8781, -87.6298, dr=0.06, da=8, col=C)
label("O'Hare", 41.9742, -87.9073, dr=0.06, da=8, col=A, fs=6.5)
label("Michigan State", 42.7018, -84.4822, dr=-0.02, da=-3, ha="right", col=G)
# north sector
label("near beam: north over the ComEd ROW\nSmith Rd 10 km (houses), Wayne 12 km", 41.9356, -88.2230,
      dr=0.18, da=-25, ha="right", col=K, fs=7.2, bold=True)
label("Elgin", 42.0354, -88.2826, dr=0.06, da=8, col=C, fs=6.5)
label("untilted N exit, 38 km", 42.1884, -88.2230, dr=0.02, da=12, col=K, fs=6.5)
label("Moraine Hills SP", 42.30, -88.23, dr=0.09, da=10, col=P)
label("Chain O'Lakes SP", 42.46, -88.20, dr=0.0, da=-10, ha="right", col=P)
label("Rockford", 42.2711, -89.0940, dr=0.04, da=8, col=C)
label("Milwaukee", 43.0389, -87.9065, dr=0.05, da=7, col=C)
label("Kettle Moraine SF\n(195 km due N: the mirror aim)", 43.60, -88.20, dr=-0.06, da=-8, ha="right", col=P)
label("UW-Madison", 43.0766, -89.4125, dr=0.06, da=-4, ha="right", col=G)
label("U. Iowa", 41.6611, -91.5302, dr=0.05, ha="right", col=G)
label("Chequamegon-Nicolet & Ottawa NFs", 45.75, -88.85, dr=-0.04, da=-8, ha="right", col=P)
label("Menominee & Marquette iron ranges", 45.82, -88.07, dr=0.09, da=7, col=M)
label("Hiawatha NF", 46.05, -86.75, dr=-0.04, da=10, col=P)
label("Iron Range", 47.427, -92.938, dr=-0.20, da=-13, ha="right", col=M)
label("Soudan mine (MINOS far, 735 km)", 47.8232, -92.2373, dr=-0.26, da=-4, ha="right", col=M)
label("Ash River (NOvA far)", 48.3794, -92.8313, dr=0.12, da=6, col=M)
label("SURF / Homestake\n(DUNE far, 1290 km)", 44.3525, -103.7510, dr=-0.08, da=8, col=M)

# lake names
for nm, la, lo in (("Lake Michigan", 44.0, -86.5), ("Lake Superior", 48.3, -87.5),
                   ("Lake Huron", 45.0, -82.4), ("Lake Erie", 42.2, -81.2),
                   ("L. Ontario", 43.7, -77.9), ("L. Winnebago", 44.15, -88.42)):
    D, a = inv(la, lo)
    ax.text(math.radians(a), rmap(D), nm, fontsize=7.5, color="#1f5a7a",
            ha="center", va="center", style="italic", zorder=8,
            bbox=dict(boxstyle="round,pad=0.1", fc="w", ec="none", alpha=0.55))

# range rings, labelled in km AND the tilt that reaches them
RINGS = (20, 50, 100, 200, 500, 1000)
ax.set_yticks([float(rmap(d)) for d in RINGS])
ax.set_yticklabels(["%g km\n%.1f mrad" % (d, d / (2 * R_E) * 1e3) for d in RINGS],
                   fontsize=6.6, color="0.4")
ax.set_rlabel_position(260)
ax.set_xticks(np.radians(np.arange(0, 360, 45)))
ax.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"], fontsize=9)
ax.grid(color="0.86", lw=0.6)
ax.text(math.radians(290), rmap(260), "collider circle: 15.40 mrad from 35 m\nsurfaces %.0f km out in every bearing" % D_COLL,
        fontsize=7.4, color="#0b6e4f", ha="center", va="center",
        bbox=dict(boxstyle="round,pad=0.15", fc="w", ec="none", alpha=0.8), zorder=9)
ax.set_title("The collider in its default aim, and everything around it\n"
             "bearing around, log range radial (each ring also gives the tilt that reaches it); "
             "dashed = the 15.40 mrad exit circle",
             fontsize=9.5, pad=24)

from matplotlib.lines import Line2D
from matplotlib.patches import Patch
H = [Line2D([], [], marker="*", ms=10, color=G, ls="", mec="w", label="university partner"),
     Line2D([], [], marker="s", ms=6, color=F, ls="", mec="w", label="DOE / federal land"),
     Line2D([], [], marker="D", ms=5.5, color=M, ls="", mec="w", label="mine, quarry, far-detector site"),
     Line2D([], [], marker="h", ms=7, color="w", mec=M, mew=1.2, ls="", label="Mt. Simon gas storage"),
     Line2D([], [], marker="^", ms=7, color=P, ls="", mec="w", label="public land"),
     Line2D([], [], marker="o", ms=4.5, color=C, ls="", mec="w", label="population centre"),
     Line2D([], [], marker="x", ms=6, color=A, ls="", label="airport"),
     Line2D([], [], marker="|", ms=9, color=K, ls="", label="corridor mark on the meridian"),
     Line2D([], [], color="#0b6e4f", lw=2, ls=(0, (5, 3)), label="15.40 mrad exit circle"),
     Line2D([], [], color="#0b6e4f", lw=1.3, label=r"$\theta_{free}$(bearing); shaded = azimuth free"),
     Patch(color="#3d7ea6", alpha=.55, label="lake (Natural Earth 10m footprint)"),
     Line2D([], [], color=K, lw=3.5, label="near beam below 500 ft (north)")]
fig.legend(handles=H, loc="lower center", ncol=4, fontsize=7.0, frameon=False,
           bbox_to_anchor=(0.5, 0.052), handletextpad=0.5, columnspacing=1.2)
fig.text(0.5, 0.016,
         "Lakes are real footprints projected into this frame; every landmark is placed by its true bearing and range from the IP. "
         "The collider's whole 198 km circle lies inside $\\theta_{free}$,\n"
         "so at 15.4 mrad the near beam is corridor-bound in every bearing $-$ which is why the study rides the ComEd "
         "right-of-way north. Ranges and bearings: static/geo/collider_map.json.",
         ha="center", va="center", fontsize=7.0, color="0.35")

for ext in ("pdf", "svg"):
    fig.savefig(os.path.join(ROOT, "static", "figs", "collider_map." + ext),
                bbox_inches="tight")
json.dump(dict(collider_tilt_mrad=15.40, straight_depth_m=D0, exit_circle_km=round(D_COLL, 1),
               note="tilt_to_reach_mrad = D/(2 R_E) - d0/D for a 35 m straight; circle[] gives where the far "
                    "exit lands per degree of bearing and the near-side sub-500 ft band beyond the fence",
               water_arcs_on_circle=WATER_ARCS, circle=CIRCLE,
               landmarks=sorted(OUT, key=lambda d: d["range_km"])),
          open(os.path.join(ROOT, "static", "geo", "collider_map.json"), "w"), indent=1)
print("wrote collider_map; %d landmarks; circle %.0f km" % (len(OUT), D_COLL))
print("water arcs on the circle:", WATER_ARCS)
for a in (0, 24, 45, 90, 145, 180, 225, 270, 315):
    c = CIRCLE[a]
    print("  az %3d -> %.3f, %.3f  water=%s  near fence %4d m  off-site band %.1f km" % (a, c["lat"], c["lon"], c["water"], c["near_fence_m"], c["offsite_sub500ft_km"]))
print("on the collider circle (+-15 km):",
      [d["name"] for d in OUT if d["on_collider_circle"]])
