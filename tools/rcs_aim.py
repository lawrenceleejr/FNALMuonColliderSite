#!/usr/bin/env python3
"""Independent (depth, angle) aiming for the last RCS -- including water targets.

The co-tilted chain variant (tools/uiuc_chain.py) locks each ring's straight
depth to its tilt through the UIUC exit constraint.  This tool inverts the
parameterisation: **depth and angle are free inputs**, and the outputs are
where each end of the straight surfaces.  That is the right frame for asking
"can the last RCS point somewhere else, further away, in water?"

Master relation.  A chord that re-surfaces a great-circle distance D from a
straight at depth d0 needs a tilt

    theta = d0/D - D/(2 R_E)   (positive = up to north)

so |theta| ~ D/(2 R_E) for any D >> d0: the tilt IS the range, ~1.57 mrad per
100 km, independent of depth.  Depth only sets where the NEAR end surfaces,
s_near ~ d0/|theta|.  Hence the two knobs are cleanly separable:

    angle  ->  how far the far exit lands
    depth  ->  how soon the near exit lands (and nothing else)

Land/water classification comes from a cached GEBCO 2020 profile along the
corridor meridian (data/gebco_meridian.json, 0.02 deg = 2.2 km sampling) plus
Natural Earth 10m lake polygons (data/ne_10m_lakes.geojson).  GEBCO is the
authority here: the Natural Earth polygon is coarse enough to swallow the
Keweenaw Peninsula (46.95-47.46 N, up to +416 m) into Lake Superior, which
would put a "water" exit on dry land.

Outputs static/geo/rcs_aim.json.
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

R_E = 6371000.0
MER = -88.222972
IP_LAT = 41.8443
MLAT = 111132.92 - 559.82 * math.cos(2 * math.radians(IP_LAT)) \
    + 1.175 * math.cos(4 * math.radians(IP_LAT))
FENCE_N, FENCE_S = 41.8699, 41.82065
SUPERIOR_SURFACE = 183.0        # m MSL
FT = 0.3048

# ----------------------------------------------------------------------
# Terrain / bathymetry along the meridian
# ----------------------------------------------------------------------
_g = json.load(open(os.path.join(ROOT, "data", "gebco_meridian.json")))
GEB = _g["profile"]                                   # [[lat, elev_m], ...]

def gebco(lat):
    if lat <= GEB[0][0]:
        return GEB[0][1]
    if lat >= GEB[-1][0]:
        return GEB[-1][1]
    for (a, e1), (b, e2) in zip(GEB, GEB[1:]):
        if a <= lat <= b:
            return e1 + (lat - a) / (b - a) * (e2 - e1)
    return GEB[-1][1]

_lk = json.load(open(os.path.join(ROOT, "data", "ne_10m_lakes.geojson")))

def _rings(g):
    if g["type"] == "Polygon":
        return [g["coordinates"][0]]
    if g["type"] == "MultiPolygon":
        return [p[0] for p in g["coordinates"]]
    return []

def _pin(x, y, poly):
    ins = False
    for (x1, y1), (x2, y2) in zip(poly, poly[1:] + [poly[0]]):
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1:
            ins = not ins
    return ins

def lake_name(lat):
    for f in _lk["features"]:
        for r in _rings(f["geometry"]):
            lons = [c[0] for c in r]
            lats = [c[1] for c in r]
            if min(lons) - .2 <= MER <= max(lons) + .2 and \
               min(lats) - .2 <= lat <= max(lats) + .2:
                if _pin(MER, lat, [(c[0], c[1]) for c in r]):
                    return f["properties"].get("name") or "(unnamed lake)"
    return None

def local_relief(lat, half=0.06):
    vals = [e for la, e in GEB if lat - half <= la <= lat + half]
    return (max(vals) - min(vals)) if vals else 0.0

def classify(lat):
    """Land/water at a meridian latitude.

    GEBCO decides, the lake polygon only names the body.  Two traps this
    guards against, both real on this meridian:
      * NE's coarse Lake Superior polygon swallows the KEWEENAW PENINSULA
        (46.90-47.46 N, +30 to +416 m) -- dry land, 400 m of relief.
      * GEBCO's 15-arcsec grid cannot resolve a narrow inland reservoir:
        at Kentucky Lake it returns 100-128 m (28 m of relief), i.e. the
        river valley averaged with its bluffs.  Such a claim is reported
        as UNVERIFIED, never as a water target.
    """
    e = gebco(lat)
    nm = lake_name(lat)
    if e < 0:
        return dict(water=True, body=nm or "open water (below MSL)",
                    bed_masl=e,
                    water_column_m=round(
                        (SUPERIOR_SURFACE - e) if (nm or "").startswith("Lake Superior")
                        else -e, 0))
    if nm:
        flat = local_relief(lat) < 3.0
        return dict(water=bool(flat), body=nm, bed_masl=e, flat=bool(flat),
                    verdict=("flat water surface" if flat else
                             "polygon claims lake, GEBCO shows %.0f m of relief "
                             "in +-6.6 km: land or unresolved reservoir"
                             % local_relief(lat)))
    return dict(water=False, body=None, ground_masl=e)

# ----------------------------------------------------------------------
# Geometry: depth and angle as independent knobs
# ----------------------------------------------------------------------
def surf(lat):
    return gebco(lat)

def beam_z(z0, th, y):
    return z0 + th * y + y * y / (2 * R_E)

def crossing(z0, th, lo, hi):
    f = lambda y: beam_z(z0, th, y) - surf(IP_LAT + y / MLAT)
    a, b, fa = lo, hi, None
    fa = f(a)
    if fa * f(b) > 0:
        return None
    for _ in range(220):
        m = 0.5 * (a + b)
        if (f(m) > 0) == (fa > 0):
            a, fa = m, f(m)
        else:
            b = m
    return 0.5 * (a + b)

def theta_for_range(depth, D, north=True):
    """Tilt that puts the far exit D metres away (sign: + = up to north)."""
    t = depth / D - D / (2 * R_E)
    return t if north else -t

def aim(depth, theta_mrad, label=""):
    """Full report for a straight at `depth` below grade, tilt `theta_mrad`
    (positive = up to north)."""
    th = theta_mrad * 1e-3
    z0 = surf(IP_LAT) - depth
    out = dict(label=label, depth_m=depth, theta_mrad=theta_mrad)
    for side, sgn, lo, hi in (("north", +1, 100.0, 1.30e6),
                              ("south", -1, -1.30e6, -100.0)):
        y = crossing(z0, th, lo, hi) if sgn > 0 else crossing(z0, th, lo, hi)
        if y is None:
            out[side] = dict(exit=None)
            continue
        lat = IP_LAT + y / MLAT
        D = abs(y)
        cls = classify(lat)
        graze = abs(th + y / R_E)
        # chord perigee (below the geoid) and its latitude
        y_per = -th * R_E
        per = None
        if (y_per > 0) == (sgn > 0) and abs(y_per) < D:
            per = dict(lat=round(IP_LAT + y_per / MLAT, 3),
                       depth_km=round((surf(IP_LAT + y_per / MLAT)
                                       - beam_z(z0, th, y_per)) / 1000, 2))
        out[side] = dict(exit_lat=round(lat, 4), range_km=round(D / 1000, 1),
                         graze_mrad=round(graze * 1e3, 2), perigee=per,
                         **{k: (round(v, 1) if isinstance(v, float) else v)
                            for k, v in cls.items()})
    # civic envelope on the near side = whichever end surfaces closer in
    cands = [(s, out[s]["range_km"]) for s in ("north", "south")
             if out[s].get("exit_lat") is not None]
    near = min(cands, key=lambda kv: kv[1])[0]
    out["far_side"] = "north" if near == "south" else "south"
    ne = out[near]
    y_em = (ne["exit_lat"] - IP_LAT) * MLAT
    sgn = 1 if near == "north" else -1
    fence_y = (FENCE_N - IP_LAT) * MLAT if near == "north" \
        else (FENCE_S - IP_LAT) * MLAT
    h_fence = beam_z(z0, th, fence_y) - surf(IP_LAT + fence_y / MLAT)
    y500 = y_em + sgn * (500 * FT) / abs(th)
    out["near_side"] = dict(
        side=near, emerges_on_site=bool(abs(y_em) < abs(fence_y)),
        height_at_fence_m=round(h_fence, 0),
        reaches_500ft_at_km=round(abs(y500) / 1000, 2),
        sub500ft_offsite_km=round(max(0.0, (abs(y500) - abs(fence_y))) / 1000, 2))
    return out

# ----------------------------------------------------------------------
# Ring feasibility: planar span vs tilt, and the achromat alternative
# ----------------------------------------------------------------------
LS34, R34, C34 = 450.0, 2200.0, 14720.0

def ring_feasibility(theta_mrad, depth):
    th = theta_mrad * 1e-3
    span = LS34 / 2 + R34                     # m from straight centre to apex
    lift = span * abs(th)          # the shallow apex rises by this, either sign
    return dict(planar_apex_lift_m=round(lift, 1),
                planar_shallow_apex_depth_m=round(depth - lift, 1),
                planar_ok=bool(depth - lift >= 6.0),
                achromat_theta_mrad=theta_mrad,
                achromat_BL_Tm_at_5TeV=round(abs(th) * 5000 / 0.2998, 0),
                achromat_len_m_at_8T=round(abs(th) * 5000 / 0.2998 / 8, 0),
                achromat_frac_of_ring_pct=round(
                    2 * (abs(th) * 5000 / 0.2998 / 8) / C34 * 100, 2))

# ----------------------------------------------------------------------
# Dose (King eq.10 raw, bounding each stage at top energy)
# ----------------------------------------------------------------------
NDEC = {"RCS3": 1.86e17, "RCS4": 1.73e17}
ETOP = {"RCS3": 1.5, "RCS4": 5.0}

def king_raw_mSv(stage, L_m):
    return 1.1e-18 * NDEC[stage] * ETOP[stage] ** 4 / (L_m / 1000) ** 2 * 1e3

# ----------------------------------------------------------------------
# Candidate scan: every water crossing the meridian offers
# ----------------------------------------------------------------------
def water_windows():
    """CONTIGUOUS meridian bands of sub-MSL GEBCO elevation (no positive gap
    tolerated), with the tilt each end needs from a 20 m-deep straight."""
    out, cur = [], None
    for lat, e in GEB:
        if e < 0:
            body = lake_name(lat) or "open water (below MSL)"
            if cur is None:
                cur = [lat, lat, body, e, e]
            else:
                cur[1] = lat
                cur[3] = min(cur[3], e)
                cur[4] = max(cur[4], e)
        elif cur:
            out.append(cur)
            cur = None
    if cur:
        out.append(cur)
    res = []
    for a, b, body, emin, emax in out:
        if (b - a) * 111.13 < 8:            # skip sub-8 km slivers
            continue
        Da = abs(a - IP_LAT) * MLAT
        Db = abs(b - IP_LAT) * MLAT
        north = a > IP_LAT
        res.append(dict(
            body=body, lat_range=[round(a, 3), round(b, 3)],
            side="N" if north else "S",
            length_km=round((b - a) * 111.13, 1),
            range_km=[round(min(Da, Db) / 1000, 1), round(max(Da, Db) / 1000, 1)],
            tilt_window_mrad=[round(abs(theta_for_range(20.0, max(Da, Db))) * 1e3, 2),
                              round(abs(theta_for_range(20.0, min(Da, Db))) * 1e3, 2)],
            deepest_bed_masl=emin,
            max_water_column_m=round(SUPERIOR_SURFACE - emin, 0) if north and emin < 0 else None))
    return res

# ----------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------
WIN = water_windows()
# best target: the longest open-water window
best = max(WIN, key=lambda w: w["length_km"])
Dmid = 0.5 * (best["range_km"][0] + best["range_km"][1]) * 1000
th_best = theta_for_range(20.0, Dmid, north=(best["side"] == "N")) * 1e3

# th_best is the tilt for mid-target; sign already set by the window's side
SCEN = [
    aim(35.0, 15.40, "UIUC, south aim (co-tilted chain reference)"),
    aim(15.0, th_best, "Lake Superior, 15 m straight"),
    aim(35.0, th_best, "Lake Superior, 35 m straight"),
    aim(120.0, th_best, "Lake Superior, 120 m straight"),
    aim(440.0, th_best, "Lake Superior, 440 m straight (deep-ring option)"),
    aim(35.0, 100.61, "Gulf of Mexico, south aim (tilt disqualifier)"),
]
report = {
    "master_relation": "theta = d0/D - D/(2 R_E); |theta| ~ 1.57 mrad per 100 km "
                       "of range, independent of depth. Depth sets only the "
                       "near-side emergence distance s ~ d0/|theta|.",
    "data": {"bathymetry": _g["dataset"], "lake_polygons": "Natural Earth 10m",
             "caveat": "GEBCO decides land/water; the NE polygon only names the "
                       "body. NE's coarse Lake Superior polygon includes the "
                       "Keweenaw Peninsula (46.95-47.46 N, up to +416 m), which "
                       "is land."},
    "water_windows_on_the_meridian": WIN,
    "recommended_target": best,
    "scenarios": SCEN,
    "ring_feasibility": {
        s["label"]: ring_feasibility(s["theta_mrad"], s["depth_m"]) for s in SCEN},
    "exit_dose_raw_mSv_yr": {
        s["label"]: {
            side: {st: round(king_raw_mSv(st, s[side]["range_km"] * 1000), 4)
                   for st in NDEC}
            for side in ("north", "south") if s[side].get("exit_lat")}
        for s in SCEN},
}
os.makedirs(os.path.join(ROOT, "static", "geo"), exist_ok=True)
with open(os.path.join(ROOT, "static", "geo", "rcs_aim.json"), "w") as f:
    json.dump(report, f, indent=1)

print("=== water windows on the corridor meridian (GEBCO-verified) ===")
for w in WIN:
    print(" %-34s %-14s %5s  %6.1f km long  D=%.0f-%.0f km  tilt %.2f-%.2f mrad"
          % (w["body"], "%.3f-%.3f" % tuple(w["lat_range"]), w["side"],
             w["length_km"], w["range_km"][0], w["range_km"][1],
             *w["tilt_window_mrad"]))
print("\n=== scenarios (depth and angle independent) ===")
for s in SCEN:
    n, so = s["north"], s["south"]
    print(" %-42s d0=%4.0f m  th=%+7.2f mrad" % (s["label"], s["depth_m"], s["theta_mrad"]))
    for side, r in (("N", n), ("S", so)):
        if not r.get("exit_lat"):
            continue
        tag = ("WATER: " + str(r["body"])) if r.get("water") else "land"
        per = (" perigee %.2f km @ %.2f N" % (r["perigee"]["depth_km"], r["perigee"]["lat"])
               if r.get("perigee") else "")
        print("    %s exit %8.1f km  %.4f N  graze %5.2f mrad  %s%s"
              % (side, r["range_km"], r["exit_lat"], r["graze_mrad"], tag, per))
    ns = s["near_side"]
    print("    near side (%s): on site %s, %+.0f m at the fence, 500 ft at %.2f km, "
          "%.2f km sub-500 ft off site"
          % (ns["side"], ns["emerges_on_site"], ns["height_at_fence_m"],
             ns["reaches_500ft_at_km"], ns["sub500ft_offsite_km"]))
    rf = report["ring_feasibility"][s["label"]]
    print("    ring: planar apex lift %.0f m -> %s; achromat %.0f T.m at 5 TeV "
          "(%.0f m of 8 T x2 = %.2f%% of the ring)"
          % (rf["planar_apex_lift_m"],
             "OK" if rf["planar_ok"] else "IMPOSSIBLE (apex above grade)",
             rf["achromat_BL_Tm_at_5TeV"], rf["achromat_len_m_at_8T"],
             rf["achromat_frac_of_ring_pct"]))
    dz = report["exit_dose_raw_mSv_yr"][s["label"]]
    print("    King-raw dose at exits: " + "; ".join(
        "%s %s" % (side, ", ".join("%s %.3g mSv/yr" % (st, v) for st, v in d.items()))
        for side, d in dz.items()))
print("\nWrote static/geo/rcs_aim.json")
