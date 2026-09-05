#!/usr/bin/env python3
"""One straight, both ends far: what it costs in depth.

A straight whose two ends surface at great-circle ranges s1 and s2 on
opposite bearings is a chord of the Earth.  Its depth below the IP and its
tilt follow from the two-exit conditions s_near s_far = 2 R_E d0 and
s_far - s_near = 2 R_E theta:

    d0 = s1 s2 / (2 R_E),     theta = (s2 - s1) / (2 R_E),
    deepest point ((s1+s2)/2)^2 / (2 R_E), located (s2-s1)/2 toward the far end.

With s1 fixed at the UIUC South Farms (198.4 km) this gives the depth needed
for the other end to come out at any northern target -- Green Bay, Lake
Superior -- and compares it with real excavation depths.  It also sizes the
feasible alternative: the ring's two straights aimed separately (east
straight up-north to UIUC, west straight down-north to the water), which
needs inclined arcs and vertical achromats but no extra depth.

Output: static/figs/two_ends.{svg,pdf}; static/geo/two_ends.json.
"""
import json, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import compass_common as cc
STYLE = "/root/.claude/skills/synced/tufte-plots/assets/tufte.mplstyle"
if os.path.exists(STYLE):
    plt.style.use(STYLE)
plt.rcParams["figure.dpi"] = 130

R = cc.R_E
S_UIUC = cc.inv(40.0602, -88.2230)[0]                    # 198.4 km at bearing 180.0

def chord(s1, s2):
    """Depth at the IP (km), tilt (rad), deepest point (km) and its offset (km) for exits at s1, s2."""
    return dict(d0_km=s1 * s2 / (2 * R), tilt_mrad=(s2 - s1) / (2 * R) * 1e3,
                perigee_km=((s1 + s2) / 2) ** 2 / (2 * R), perigee_offset_km=(s2 - s1) / 2)

# ---------------------------------------------------------------- where is Green Bay from here?
_lk = json.load(open(os.path.join(ROOT, "data", "ne_10m_lakes.geojson")))
LM = max((r for f in _lk["features"] if f["properties"].get("name") == "Lake Michigan"
          for r in cc._rings(f["geometry"])), key=len)
LM = [(c[0], c[1]) for c in LM]
def in_michigan(az, D):
    la, lo = cc.fwd(az, D)
    return cc.pin(lo, la, LM)
GB = []                                                 # (az, first water km, last water km) for the bay, 250-400 km
for az10 in range(0, 200):                              # 0 .. 19.9 deg in 0.1 steps
    az = az10 / 10.0
    wet = [D for D in np.arange(250, 400, 1.0) if in_michigan(az, D)]
    if wet:
        GB.append((az, float(min(wet)), float(max(wet))))
gb_az, gb_d1, gb_d2 = GB[0] if GB else (None, None, None)
gb_offset_uiuc = S_UIUC * math.sin(math.radians(gb_az)) if gb_az is not None else None
gb_south_end = cc.fwd(180 + gb_az, S_UIUC) if gb_az is not None else None
print("Green Bay water first reachable at bearing %.1f deg, %.0f-%.0f km; the UIUC end then lands %.1f km west of the "
      "South Farms at %.4f N %.4f W" % (gb_az, gb_d1, gb_d2, gb_offset_uiuc, gb_south_end[0], -gb_south_end[1]))

# ---------------------------------------------------------------- targets north on the UIUC line (bearing 0)
gb_miss = gb_d1 * math.sin(math.radians(gb_az))           # km between the UIUC line and the bay's water at that range
TARGETS = [("north fence, on site", 2.85), ("Zone D level-ring exit (Cary / Fox River Grove)", 38.2),
           ("Kettle Moraine SF", 195.2), ("Shawano County, WI (Green Bay's latitude, %.0f km west of the bay)" % gb_miss, 298.0),
           ("Green Bay water (needs bearing %.1f$^\\circ$)" % gb_az, gb_d1),
           ("Lake Superior, south shore", 626.0), ("Lake Superior, mid-crossing", 655.0), ("Lake Superior, north shore", 695.0)]
REF = [("collider straight", 0.035), ("deep reference hall", 0.107), ("Gotthard base tunnel, max overburden", 2.3),
       ("SNOLAB", 2.07), ("Mponeng mine, deepest excavation", 4.0), ("Kola superdeep borehole", 12.26)]
rows = []
for nm, s2 in TARGETS:
    c = chord(S_UIUC, s2); c.update(target=nm, north_exit_km=s2); rows.append(c)
    print("%-62s s2 %6.1f km -> depth at IP %7.3f km, deepest %7.3f km, tilt %6.2f mrad" % (nm, s2, c["d0_km"], c["perigee_km"], c["tilt_mrad"]))

# ---------------------------------------------------------------- the feasible alternative: two straights, one ring
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    import corridor_layout as cl
LS, RA, C = cl.LS_RCS34, cl.R_RCS34, cl.C_RCS34            # m
P_TOP = 5000.0                                              # GeV, RCS4 top energy
def two_straights(s_far_west_km, depth_m):
    d = depth_m / 1000.0
    thA = S_UIUC / (2 * R) - d / S_UIUC                     # east straight, up-north, down-going end at UIUC
    thB = s_far_west_km / (2 * R) - d / s_far_west_km       # west straight, down-north, down-going end in the water
    rise_per_lap = (thA + thB) * LS                         # m: both straights rise in the direction of travel
    slope = rise_per_lap / (2 * math.pi * RA)               # each arc descends half of it over pi R
    bend_total = 2 * (thA + thB) + 4 * slope                # rad of vertical bending per lap
    BL = bend_total * P_TOP / 0.2998                        # T m at 5 TeV
    return dict(depth_m=depth_m, far_west_km=s_far_west_km, tilt_east_mrad=thA * 1e3, tilt_west_mrad=thB * 1e3,
                near_east_km=d / thA, near_west_km=d / thB, arc_slope_mrad=slope * 1e3,
                vertical_bend_per_lap_mrad=bend_total * 1e3, BL_Tm_at_5TeV=BL, m_of_8T=BL / 8, fraction_of_ring=BL / 8 / C)
ALT = [two_straights(s, d) for s in (655.0, 300.0) for d in (35.0, 80.0)]
for a in ALT:
    print("two straights: west end %.0f km, depth %.0f m -> tilts %.1f / %.1f mrad, near ends %.2f / %.2f km, arc slope %.1f mrad, "
          "%.0f mrad of vertical bend per lap = %.0f T m = %.0f m of 8 T (%.1f%% of the ring)" % (
              a["far_west_km"], a["depth_m"], a["tilt_east_mrad"], a["tilt_west_mrad"], a["near_east_km"], a["near_west_km"],
              a["arc_slope_mrad"], a["vertical_bend_per_lap_mrad"], a["BL_Tm_at_5TeV"], a["m_of_8T"], 100 * a["fraction_of_ring"]))

# ---------------------------------------------------------------- the chords the polar plots draw
UIUC = (40.0602, -88.2230)

def lake_rings(name):
    return [[(c[0], c[1]) for c in r] for f in _lk["features"] if f["properties"].get("name") == name
            for r in cc._rings(f["geometry"])]
WIN = max(lake_rings("Lake Winnebago"), key=len)
def in_poly(poly, lat, lon):
    return cc.pin(lon, lat, poly)
def interior_points(poly, margin_km=2.0, step_deg=0.01):
    """Grid points inside the polygon whose 4 neighbours at +-margin are also inside (i.e. >= margin from shore)."""
    lats = [c[1] for c in poly]; lons = [c[0] for c in poly]
    out = []
    la = min(lats)
    while la <= max(lats):
        lo = min(lons)
        while lo <= max(lons):
            dla = margin_km / 111.132; dlo = margin_km / (111.32 * math.cos(math.radians(la)))
            if all(in_poly(poly, la + a, lo + b) for a, b in ((0, 0), (dla, 0), (-dla, 0), (0, dlo), (0, -dlo))):
                out.append((la, lo))
            lo += step_deg
        la += step_deg
    return out
def chord_record(name, short, color, point, note):
    az_gc, D_gc, dep_gc = cc.great_circle_polar(UIUC, point)
    i = int(np.argmin(D_gc))
    D_pt, az_pt = cc.inv(*point)
    via = chord(S_UIUC, D_pt)
    la1, lo1, la2, lo2 = map(math.radians, (UIUC[0], UIUC[1], point[0], point[1]))
    arc = cc.R_E * math.acos(math.sin(la1) * math.sin(la2) + math.cos(la1) * math.cos(la2) * math.cos(lo2 - lo1))
    return dict(name=name, short=short, color=color,
                point=dict(lat=round(point[0], 4), lon=round(point[1], 4), bearing_deg=round(az_pt, 1), range_km=round(D_pt, 1), note=note),
                true_chord=dict(closest_approach_to_ip_km=round(float(D_gc[i]), 2), closest_approach_bearing_deg=round(float(np.degrees(az_gc[i])) % 360, 1),
                                depth_at_closest_approach_km=round(float(dep_gc[i]), 2), midpoint_depth_km=round(float(dep_gc.max()), 2), arc_km=round(arc, 1),
                                inside_fence=bool(D_gc[i] * 1000 < cc.fence_distance(float(np.degrees(az_gc[i])) % 360))),
                through_ip=dict(uiuc_end_bearing_deg=round((180 + az_pt) % 360, 1), uiuc_end_offset_km=round(S_UIUC * math.sin(math.radians(az_pt)), 1),
                                depth_at_ip_km=round(via["d0_km"], 2), tilt_mrad=round(via["tilt_mrad"], 2), deepest_km=round(via["perigee_km"], 2)))
# Lake Winnebago: the lake point (>= 2 km from shore) whose chord to UIUC passes closest to the IP, and the lake centre
WPTS = interior_points(WIN, 2.0)
def miss(pt):
    az_gc, D_gc, _ = cc.great_circle_polar(UIUC, pt, n=300)
    return float(D_gc.min())
w_best = min(WPTS, key=miss)
w_centre = (float(np.mean([p[0] for p in WPTS])), float(np.mean([p[1] for p in WPTS])))
print("Lake Winnebago: %d interior points; chord-to-UIUC closest to the IP from %.4f N %.4f W (miss %.2f km); lake centre %.4f N %.4f W (miss %.2f km)"
      % (len(WPTS), w_best[0], -w_best[1], miss(w_best), w_centre[0], -w_centre[1], miss(w_centre)))
# the bay proper lies west of the Door Peninsula (lon < -87.75); pick the bearing with the longest run of BAY water
def bay_run(az):
    wet = [D for D in np.arange(250, 400, 1.0) if in_michigan(az, D) and cc.fwd(az, D)[1] < -87.75]
    return (az, float(min(wet)), float(max(wet))) if wet else None
runs = [r for r in (bay_run(a / 10.0) for a in range(20, 120)) if r]
best = max(runs, key=lambda r: r[2] - r[1])
gb_mid_az, gb_mid_D = best[0], 0.5 * (best[1] + best[2])
GBP = cc.fwd(gb_mid_az, gb_mid_D)
az_gc, D_gc, dep_gc = cc.great_circle_polar(UIUC, GBP)
i_min = int(np.argmin(D_gc))
via_ip = chord(S_UIUC, gb_mid_D)                               # the same two ranges, but as a straight through the IP
CHORD = dict(uiuc=UIUC, green_bay_point=dict(lat=round(GBP[0], 4), lon=round(GBP[1], 4), bearing_deg=gb_mid_az, range_km=gb_mid_D,
                                              water_run_km=[best[1], best[2]]),
             true_chord=dict(closest_approach_to_ip_km=round(float(D_gc[i_min]), 2), closest_approach_bearing_deg=round(float(np.degrees(az_gc[i_min])) % 360, 1),
                             depth_at_closest_approach_km=round(float(dep_gc[i_min]), 2), midpoint_depth_km=round(float(dep_gc.max()), 2),
                             arc_km=round(float(cc.R_E * 2 * math.asin(min(1.0, math.sqrt(0) + 0) or 0) + 0), 1)),
             through_ip=dict(uiuc_end_bearing_deg=round((180 + gb_mid_az) % 360, 1), uiuc_end_offset_km=round(S_UIUC * math.sin(math.radians(gb_mid_az)), 1),
                             depth_at_ip_km=round(via_ip["d0_km"], 2), tilt_mrad=round(via_ip["tilt_mrad"], 2), deepest_km=round(via_ip["perigee_km"], 2)))
la1, lo1, la2, lo2 = map(math.radians, (UIUC[0], UIUC[1], GBP[0], GBP[1]))
CHORD["true_chord"]["arc_km"] = round(cc.R_E * math.acos(math.sin(la1) * math.sin(la2) + math.cos(la1) * math.cos(la2) * math.cos(lo2 - lo1)), 1)
print("chord for the plots: Green Bay mid-bay point %.4f N %.4f W (bearing %.1f, %.0f km, %.0f km of water); true UIUC-Green Bay chord passes "
      "%.1f km from the IP at bearing %.0f, %.2f km deep there, %.2f km at midpoint, arc %.0f km; through-the-IP version: UIUC end %.0f km W of the "
      "South Farms, %.2f km deep at the IP, tilt %.1f mrad" % (GBP[0], -GBP[1], gb_mid_az, gb_mid_D, best[2] - best[1],
      CHORD["true_chord"]["closest_approach_to_ip_km"], CHORD["true_chord"]["closest_approach_bearing_deg"], CHORD["true_chord"]["depth_at_closest_approach_km"],
      CHORD["true_chord"]["midpoint_depth_km"], CHORD["true_chord"]["arc_km"], CHORD["through_ip"]["uiuc_end_offset_km"], CHORD["through_ip"]["depth_at_ip_km"],
      CHORD["through_ip"]["tilt_mrad"]))

CHORDS = [chord_record("Green Bay", "Green Bay", "#2a7f9e", GBP, "mid-bay water, the bearing with the longest run of bay water"),
          chord_record("Lake Winnebago", "Winnebago", "#6b8e23", w_best, "the lake point (>= 2 km from shore) whose chord to UIUC passes closest to the IP"),
          chord_record("Lake Winnebago (centre)", "Winnebago centre", "#6b8e23", w_centre, "lake centre, for reference; not drawn")]
CHORD["chords"] = CHORDS
for c in CHORDS:
    print("chord UIUC-%s: point %.4f N %.4f W (bearing %.1f, %.0f km); true chord misses the IP by %.2f km at bearing %.0f (%s the fence), %.2f km deep there, "
          "%.2f at midpoint; through the IP: UIUC end %.1f km off the South Farms, %.2f km deep, tilt %.1f mrad" % (
          c["name"], c["point"]["lat"], -c["point"]["lon"], c["point"]["bearing_deg"], c["point"]["range_km"], c["true_chord"]["closest_approach_to_ip_km"],
          c["true_chord"]["closest_approach_bearing_deg"], "inside" if c["true_chord"]["inside_fence"] else "outside",
          c["true_chord"]["depth_at_closest_approach_km"], c["true_chord"]["midpoint_depth_km"], c["through_ip"]["uiuc_end_offset_km"],
          c["through_ip"]["depth_at_ip_km"], c["through_ip"]["tilt_mrad"]))
TARGETS_W = chord(S_UIUC, CHORDS[1]["point"]["range_km"]); TARGETS_W.update(target="Lake Winnebago (needs bearing %.1f$^\\circ$)" % CHORDS[1]["point"]["bearing_deg"],
                                                                            north_exit_km=CHORDS[1]["point"]["range_km"])
rows.append(TARGETS_W)


# ---------------------------------------------------------------- figure
s2 = np.geomspace(2.0, 1000.0, 400)
fig, ax = plt.subplots(figsize=(8.4, 5.6))
fig.subplots_adjust(left=0.1, right=0.97, top=0.86, bottom=0.2)
ax.plot(s2, S_UIUC * s2 / (2 * R), color="#0b6e4f", lw=2.0, label="depth of the straight at the IP, $d_0 = s_1 s_2 / 2R_E$")
ax.plot(s2, ((S_UIUC + s2) / 2) ** 2 / (2 * R), color="#0b6e4f", lw=1.0, ls="--", label="deepest point of the chord")
for nm, d in REF:
    ax.axhline(d, color="0.6", lw=0.6, ls=":")
    if nm.startswith("Gotthard"):
        ax.text(950, d * 1.08, nm, fontsize=6.6, color="0.4", va="bottom", ha="right")
    else:
        ax.text(2.1, d * 1.08, nm, fontsize=6.6, color="0.4", va="bottom")
ax.axhspan(0.03, 0.15, color="#0b6e4f", alpha=.07, lw=0)
ax.text(2.1, 0.062, "this study's tunnels", fontsize=6.6, color="#0b6e4f", va="center")
for r in rows:
    if r["north_exit_km"] is None:
        continue
    ax.plot([r["north_exit_km"]], [r["d0_km"]], "o", ms=5, color="#b5541c", mec="w", zorder=5)
short = {"north fence, on site": "on-site fence", "Zone D level-ring exit (Cary / Fox River Grove)": "Cary / Fox River Grove\n(level ring's exit)",
         "Kettle Moraine SF": "Kettle Moraine", "Lake Superior, south shore": "", "Lake Superior, north shore": "",
         "Lake Superior, mid-crossing": "Lake Superior"}
for r in rows:
    nm = r["target"]; lab = short.get(nm, nm)
    if nm.startswith("Shawano"):
        lab = "Green Bay's latitude on the UIUC line\n(Shawano Co., %.0f km W of the bay)" % gb_miss
    if nm.startswith("Green Bay water"):
        lab = "Green Bay water\n(bearing %.1f$^\\circ$, not on the UIUC line)" % gb_az
    OFF = {"north fence, on site": ((28, -18), "left", "top"), "Kettle Moraine SF": ((-10, 14), "right", "bottom"),
           "Lake Superior, mid-crossing": ((0, -26), "center", "top")}
    off, ha, va = OFF.get(nm, ((0, -26), "center", "top"))
    if nm.startswith("Shawano"):
        off, ha, va = (-4, -44), "center", "top"
    if nm.startswith("Green Bay water"):
        off, ha, va = (16, 18), "left", "bottom"
    if nm.startswith("Lake Winnebago"):
        lab = "Lake Winnebago\n(bearing %.0f$^\\circ$; the true chord\ncrosses the site 3.9 km W of the IP)" % CHORDS[1]["point"]["bearing_deg"]
        off, ha, va = (-14, 26), "right", "bottom"
    if lab:
        ax.annotate("%s\n%.2f km deep" % (lab, r["d0_km"]) if r["d0_km"] >= 1 else "%s\n%.0f m deep" % (lab, r["d0_km"] * 1000),
                    (r["north_exit_km"], r["d0_km"]), xytext=off, textcoords="offset points", fontsize=6.6,
                    color="#b5541c", ha=ha, va=va, arrowprops=dict(arrowstyle="-", color="#b5541c", lw=0.5))
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(2, 1000); ax.set_ylim(0.01, 30)
ax.set_xlabel("range of the OTHER end, north on the UIUC line, $s_2$ (km)", fontsize=8.5)
ax.set_ylabel("depth (km)", fontsize=8.5)
ax.tick_params(labelsize=7.5)
ax.legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(0.30, 1.0))
ax.set_title("One straight with BOTH ends far: south end at UIUC (198 km), north end at $s_2$\n"
             "the straight must sit $d_0 = s_1 s_2/2R_E$ below the IP $-$ kilometres, for any northern water", fontsize=9.5, pad=10)
fig.text(0.5, 0.035,
         "A straight is a chord: its two exits obey $s_{near}\\,s_{far} = 2R_E d_0$. Pushing the near exit from 2 km to 300 km (Green Bay's latitude) "
         "takes the straight from 35 m to 4.7 km;\nLake Superior needs 10 km, deeper than any mine. The feasible route is the ring's two straights aimed "
         "separately at ordinary depth (see text): inclined arcs and ~%d mrad of vertical bending per lap." % round(ALT[0]["vertical_bend_per_lap_mrad"]),
         ha="center", va="center", fontsize=6.8, color="0.35")
for ext in ("pdf", "svg"):
    fig.savefig(os.path.join(ROOT, "static", "figs", "two_ends." + ext), bbox_inches="tight")

json.dump(dict(s1_uiuc_km=round(S_UIUC, 2), formula="d0 = s1 s2 / (2 R_E); tilt = (s2 - s1)/(2 R_E); deepest = ((s1+s2)/2)^2/(2 R_E)",
               chord_for_plots=CHORD,
               green_bay=dict(first_water_bearing_deg=gb_az, water_km=[gb_d1, gb_d2], uiuc_end_offset_km=gb_offset_uiuc,
                              miss_from_uiuc_line_km=gb_miss,
                              south_end_latlon=gb_south_end, all_bearings=GB[:40]),
               targets=rows, reference_depths_km=dict(REF), two_straights_alternative=ALT,
               rcs34=dict(Ls_m=LS, R_m=RA, C_m=C)),
          open(os.path.join(ROOT, "static", "geo", "two_ends.json"), "w"), indent=1)
print("wrote two_ends")
