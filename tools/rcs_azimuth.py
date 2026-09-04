#!/usr/bin/env python3
"""Azimuth freedom: what a steep-angled straight can point at, in any bearing.

The corridor's meridian (-88.222972) was not chosen for physics.  It was
chosen because the collider's shallow 15.40 mrad north beam needs ~9.5 km of
institutional land to climb through the private-airspace band, and the ComEd
Aurora-Wayne right-of-way is what happens to run north from the site.  That
constraint is a function of the ANGLE:

    emergence      s_em  = d0 / theta
    clears 500 ft  s_500 = (d0 + 152.4 m) / theta
    off-site band  = max(0, s_500 - s_fence)          ~ 1/theta
    far target     D     = 2 R_E theta

so a steep beam clears navigable airspace almost immediately and needs no
corridor in any direction.  Past

    theta_free = (d0 + 152.4) / s_fence      (~60 mrad for d0 = 15 m)

the near beam is above 500 ft BEFORE it leaves the fence: the bearing becomes
free, and each ring can be aimed independently of every other.

This tool maps that freedom.  Great-lake and ocean polygons are transformed
into (azimuth, required-tilt) space -- the opportunity map -- and candidate
aim points are then verified against GEBCO 2020 point queries (cached in
data/gebco_targets.json), because the polygons are only good enough to
GENERATE candidates: on the corridor meridian, Natural Earth's Lake Superior
polygon includes 400 m of the Keweenaw Peninsula.

Ring fit is checked too: rotating a 14.72 km racetrack inside the real site
boundary is not free, and that -- not the beam -- is what actually limits
azimuth.

Outputs static/geo/rcs_azimuth.json.
"""
import json
import math
import os
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

R_E = 6371.0088                       # km
IP_LAT, IP_LON = 41.8443, -88.222972
FT = 0.3048
D0_DEFAULT = 15.0                     # m, shallow straight (see rcs-aim study)

# Published surface elevations (m MSL).  These are ESSENTIAL: a plain
# "elevation < 0" water test silently rejects Huron, Erie and most of
# Michigan, whose BEDS sit above sea level (Erie's surface is 174 m and it is
# only ~64 m deep).  Water means "below the local water surface", not "below
# sea level".
LAKE_SURFACE = {
    "Lake Superior": 183.0, "Lake Michigan": 176.0, "Lake Huron": 176.0,
    "Lake Erie": 174.0, "Lake Ontario": 75.0, "Lake St. Clair": 175.0,
    "Lake Winnipeg": 217.0, "Lake Nipigon": 260.0, "Lake Simcoe": 219.0,
    "Lake Nipissing": 196.0, "Lake of the Woods": 323.0, "Rainy Lake": 337.0,
    "ocean": 0.0,
}

def water_at(nm, elev, relief=None):
    """Is a GEBCO sample submerged in body `nm`?  Two signatures:
    (a) elev well below the surface -> real bathymetry, column = surface-elev;
    (b) elev within ~3 m of the surface AND locally flat -> GEBCO carries no
        bathymetry there and is returning the water surface itself."""
    if nm not in LAKE_SURFACE or elev is None:
        return None
    s = LAKE_SURFACE[nm]
    if elev < s - 1.0:
        return dict(water=True, column_m=round(s - elev, 0), source="bathymetry")
    if abs(elev - s) <= 3.0 and (relief is None or relief <= 4.0):
        return dict(water=True, column_m=0.0, source="surface elevation "
                    "(no bathymetry in GEBCO here; depth unknown)")
    return dict(water=False, column_m=None,
                source="%.0f m, i.e. %.0f m ABOVE the %s surface" % (elev, elev - s, nm))

# ----------------------------------------------------------------------
# Spherical geodesy (adequate at these scales; the study's other pages use
# the same class of approximation, and a real aim needs a geodetic solution)
# ----------------------------------------------------------------------
def fwd(lat, lon, az_deg, D_km):
    p1, l1 = math.radians(lat), math.radians(lon)
    a, d = math.radians(az_deg), D_km / R_E
    p2 = math.asin(math.sin(p1) * math.cos(d) + math.cos(p1) * math.sin(d) * math.cos(a))
    l2 = l1 + math.atan2(math.sin(a) * math.sin(d) * math.cos(p1),
                         math.cos(d) - math.sin(p1) * math.sin(p2))
    return math.degrees(p2), (math.degrees(l2) + 540) % 360 - 180

def inv(lat, lon):
    """distance (km) and initial bearing (deg) from the IP to (lat, lon)."""
    p1, p2 = math.radians(IP_LAT), math.radians(lat)
    dl = math.radians(lon - IP_LON)
    dp = p2 - p1
    h = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    D = 2 * R_E * math.asin(min(1.0, math.sqrt(h)))
    az = math.degrees(math.atan2(math.sin(dl) * math.cos(p2),
                                 math.cos(p1) * math.sin(p2)
                                 - math.sin(p1) * math.cos(p2) * math.cos(dl))) % 360
    return D, az

def tilt_for(D_km, d0_m=D0_DEFAULT):
    """Required tilt (mrad) to surface at D_km, from depth d0_m."""
    D = D_km * 1000.0
    return (d0_m / D + D / (2 * R_E * 1000.0)) * 1e3

def near_side(theta_mrad, d0_m, s_fence_m):
    th = theta_mrad * 1e-3
    s_em = d0_m / th
    s_500 = (d0_m + 500 * FT) / th
    return dict(emergence_km=round(s_em / 1000, 3),
                clears_500ft_km=round(s_500 / 1000, 2),
                offsite_sub500ft_km=round(max(0.0, s_500 - s_fence_m) / 1000, 2),
                on_site=bool(s_em < s_fence_m))

# ----------------------------------------------------------------------
# Site boundary: fence distance vs azimuth, and racetrack fit vs azimuth
# ----------------------------------------------------------------------
_b = json.load(open(os.path.join(ROOT, "data", "fnal_boundary.geojson")))
BND = [(la, lo) for lo, la in _b["features"][0]["geometry"]["coordinates"][0]]
MLAT = 111.132
MLON = 111.320 * math.cos(math.radians(IP_LAT))
BXY = [((lo - IP_LON) * MLON * 1000, (la - IP_LAT) * MLAT * 1000) for la, lo in BND]

def _pin(x, y, poly):
    ins = False
    for (x1, y1), (x2, y2) in zip(poly, poly[1:] + [poly[0]]):
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1:
            ins = not ins
    return ins

def fence_distance(az_deg):
    """Metres from the IP to the site boundary along a bearing."""
    a = math.radians(az_deg)
    ux, uy = math.sin(a), math.cos(a)
    lo, hi = 0.0, 20000.0
    if not _pin(0, 0, BXY):
        return 0.0
    for _ in range(60):
        m = 0.5 * (lo + hi)
        if _pin(ux * m, uy * m, BXY):
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)

# --- ring fit: does ANY placement of the ring fit at this azimuth? ---------
# Ring position on site is a free design variable (the baseline nests the
# three rings at different centres), so the question is not "does the ring
# centred on the IP fit" but "does there EXIST a placement".  Answer it by
# eroding the site polygon by the rotated ring outline: build a
# distance-to-boundary grid once (positive inside), then for each azimuth take
# the element-wise min of the grid shifted by every outline point.  Any cell
# >= margin is a legal ring centre.
import numpy as _np

_GRID = 25.0                                     # m
_xs = _np.arange(min(x for x, _ in BXY) - 200, max(x for x, _ in BXY) + 200, _GRID)
_ys = _np.arange(min(y for _, y in BXY) - 200, max(y for _, y in BXY) + 200, _GRID)
_XX, _YY = _np.meshgrid(_xs, _ys)

def _build_dist():
    ins = _np.zeros(_XX.shape, dtype=bool)
    for (x1, y1), (x2, y2) in zip(BXY, BXY[1:] + [BXY[0]]):
        if (y1 > _YY) is None:
            continue
        cond = ((y1 > _YY) != (y2 > _YY)) & \
               (_XX < (x2 - x1) * (_YY - y1) / (y2 - y1 + 1e-12) + x1)
        ins ^= cond
    d = _np.full(_XX.shape, 1e9)
    for (x1, y1), (x2, y2) in zip(BXY, BXY[1:] + [BXY[0]]):
        dx, dy = x2 - x1, y2 - y1
        L2 = dx * dx + dy * dy
        if L2 == 0:
            continue
        s = _np.clip(((_XX - x1) * dx + (_YY - y1) * dy) / L2, 0.0, 1.0)
        d = _np.minimum(d, _np.hypot(_XX - x1 - s * dx, _YY - y1 - s * dy))
    return _np.where(ins, d, -d)

_DIST = _build_dist()

def _outline(Ls, R_arc, n=40):
    pts = []
    for i in range(n + 1):
        th = math.pi * i / n
        pts.append((-R_arc + R_arc * math.cos(th), Ls / 2 + R_arc * math.sin(th)))
    for i in range(n + 1):
        th = math.pi * i / n
        pts.append((-R_arc - R_arc * math.cos(th), -Ls / 2 - R_arc * math.sin(th)))
    return pts

def fit_at(az_deg, Ls, R_arc, margin=120.0):
    """Best achievable clearance at this azimuth over all ring placements."""
    a = math.radians(az_deg)
    ca, sa = math.cos(a), math.sin(a)
    best = _DIST.copy()
    for u, v in _outline(Ls, R_arc):
        x = u * ca + v * sa
        y = -u * sa + v * ca
        ix, iy = int(round(x / _GRID)), int(round(y / _GRID))
        sh = _np.full_like(_DIST, -1e9)
        h, w = _DIST.shape
        y0s, y1s = max(0, iy), min(h, h + iy)
        x0s, x1s = max(0, ix), min(w, w + ix)
        if y0s >= y1s or x0s >= x1s:
            return -1e9
        sh[y0s - iy:y1s - iy, x0s - ix:x1s - ix] = _DIST[y0s:y1s, x0s:x1s]
        best = _np.minimum(best, sh)
    return float(best.max())

def fit_window(Ls, R_arc, margin=120.0, step=2.0):
    return [a for a in [i * step for i in range(int(180 / step))]
            if fit_at(a, Ls, R_arc) >= margin]

# ----------------------------------------------------------------------
# Water polygons -> (azimuth, tilt) opportunity space
# ----------------------------------------------------------------------
def _rings(g):
    if g["type"] == "Polygon":
        return [g["coordinates"][0]]
    if g["type"] == "MultiPolygon":
        return [p[0] for p in g["coordinates"]]
    return []

BODIES = {}
_lk = json.load(open(os.path.join(ROOT, "data", "ne_10m_lakes.geojson")))
for f in _lk["features"]:
    nm = f["properties"].get("name")
    if not nm:
        continue
    for r in _rings(f["geometry"]):
        BODIES.setdefault(nm, []).extend([(c[1], c[0]) for c in r])
_oc = json.load(open(os.path.join(ROOT, "data", "ne_50m_ocean.geojson")))
for f in _oc["features"]:
    for r in _rings(f["geometry"]):
        BODIES.setdefault("ocean", []).extend([(c[1], c[0]) for c in r])

AZ_BIN = 2.0
D_MAX = 1400.0
# Only bodies GEBCO can actually adjudicate: the Great Lakes, the ocean, and
# the large Canadian lakes.  Small reservoirs are excluded on principle --
# GEBCO's 15" grid cannot resolve them (see the rcs-aim study on Kentucky
# Lake), so a "water" claim there would be unverifiable.
MAJOR = ["Lake Superior", "Lake Michigan", "Lake Huron", "Lake Erie",
         "Lake Ontario", "Lake Winnipeg", "Lake Nipigon", "Lake of the Woods",
         "Rainy Lake", "Lake Nipissing", "Lake Simcoe", "Lake St. Clair"]
# NOTE: "ocean" is deliberately excluded. Natural Earth's ocean polygon has
# the world coastline as its EXTERIOR ring, so a point-in-polygon test calls
# the whole continent ocean. The nearest real ocean water (Atlantic ~1200 km,
# Gulf ~1280 km) needs 94-101 mrad anyway -- past the useful range.
POLY = {}
for nm in MAJOR:
    if nm in BODIES:
        POLY[nm] = BODIES[nm]

def _poly_rings(nm):
    """Rebuild per-ring lon/lat polygons for a named body."""
    out = []
    src = _oc if nm == "ocean" else _lk
    for f in src["features"]:
        if nm != "ocean" and f["properties"].get("name") != nm:
            continue
        for r in _rings(f["geometry"]):
            out.append([(c[0], c[1]) for c in r])
    return out

RINGS = {nm: _poly_rings(nm) for nm in POLY}

def inside(nm, lat, lon):
    for r in RINGS[nm]:
        lons = [c[0] for c in r]
        lats = [c[1] for c in r]
        if min(lons) <= lon <= max(lons) and min(lats) <= lat <= max(lats):
            if _pin(lon, lat, r):
                return True
    return False

def runs_along(nm, az_deg, step=4.0):
    """Contiguous radial intervals (km) actually INSIDE the body along a
    bearing -- ray-cast, so crescents and multi-lobed lakes are handled."""
    out, cur = [], None
    D = 100.0
    while D <= D_MAX:
        la, lo = fwd(IP_LAT, IP_LON, az_deg, D)
        if inside(nm, la, lo):
            cur = [D, D] if cur is None else [cur[0], D]
        elif cur:
            out.append(cur)
            cur = None
        D += step
    if cur:
        out.append(cur)
    return [r for r in out if r[1] - r[0] >= step]

def opportunity_map(d0=D0_DEFAULT):
    """For each major body: the best (azimuth, radial run) it offers."""
    out = {}
    for nm in POLY:
        best = None
        for i in range(int(360 / AZ_BIN)):
            az = i * AZ_BIN
            for r in runs_along(nm, az):
                depth = r[1] - r[0]
                if best is None or depth > best["run_km"]:
                    best = dict(az_deg=az, D_km=[round(r[0], 1), round(r[1], 1)],
                                run_km=round(depth, 1),
                                tilt_mrad=[round(tilt_for(r[0], d0), 2),
                                           round(tilt_for(r[1], d0), 2)])
        if best:
            out[nm] = best
    return out

# ----------------------------------------------------------------------
# GEBCO verification of specific aim points (cached)
# ----------------------------------------------------------------------
CACHE_FN = os.path.join(ROOT, "data", "gebco_targets.json")
CACHE = json.load(open(CACHE_FN)) if os.path.exists(CACHE_FN) else {}

def gebco_batch(pts):
    """Batched GEBCO lookup (the API takes 100 locations per call).  Single-
    point calls get throttled into None-storms; batch + pace + retry."""
    import time
    need = [(la, lo) for la, lo in pts
            if "%.4f,%.4f" % (la, lo) not in CACHE]
    for i in range(0, len(need), 95):
        ch = need[i:i + 95]
        loc = "|".join("%.4f,%.4f" % p for p in ch)
        u = "https://api.opentopodata.org/v1/gebco2020?locations=" + loc
        vals = None
        for a in range(5):
            try:
                r = json.load(urllib.request.urlopen(u, timeout=90))
                vals = [x["elevation"] for x in r["results"]]
                break
            except Exception:
                time.sleep(2.5 * (a + 1))
        for j, p in enumerate(ch):
            CACHE["%.4f,%.4f" % p] = vals[j] if vals else None
        time.sleep(1.2)
    return [CACHE.get("%.4f,%.4f" % (la, lo)) for la, lo in pts]

def gebco(lat, lon):
    return gebco_batch([(lat, lon)])[0]

def verify_aim(az_deg, D_km, d0=D0_DEFAULT, label="", body=None, probe=16,
               step=8.0):
    """Sample GEBCO along the aim bearing; report the contiguous SUBMERGED run
    containing the aim point, judged against the body's surface elevation."""
    nm = body or label
    th = tilt_for(D_km, d0)
    lat, lon = fwd(IP_LAT, IP_LON, az_deg, D_km)
    pts = [fwd(IP_LAT, IP_LON, az_deg, D_km + k * step)
           for k in range(-probe, probe + 1)]
    els = gebco_batch([(round(a, 4), round(b, 4)) for a, b in pts])
    rel = (max(e for e in els if e is not None) -
           min(e for e in els if e is not None)) if any(e is not None for e in els) else None
    w = [water_at(nm, e) for e in els]
    c = probe
    ok = lambda i: w[i] is not None and w[i]["water"]
    lo_i = hi_i = c
    if ok(c):
        while lo_i > 0 and ok(lo_i - 1):
            lo_i -= 1
        while hi_i < len(w) - 1 and ok(hi_i + 1):
            hi_i += 1
    cent = water_at(nm, els[c], rel)
    cols = [w[i]["column_m"] for i in range(lo_i, hi_i + 1)
            if w[i] and w[i]["water"] and w[i]["column_m"] is not None]
    s_f = fence_distance((az_deg + 180) % 360)
    th_free = (d0 + 500 * FT) / s_f * 1e3
    return dict(label=label, body=nm, azimuth_deg=az_deg, D_km=round(D_km, 1),
                tilt_mrad=round(th, 2), aim_point=[round(lat, 4), round(lon, 4)],
                gebco_masl=els[c], is_water=bool(cent and cent["water"]),
                water_column_m=(cent or {}).get("column_m"),
                evidence=(cent or {}).get("source"),
                submerged_run_km=round((hi_i - lo_i) * step, 0) if ok(c) else 0.0,
                max_column_in_run_m=(max(cols) if cols else None),
                near_side_bearing_deg=round((az_deg + 180) % 360, 1),
                fence_distance_m=round(s_f, 0),
                near_side=near_side(th, d0, s_f),
                theta_free_mrad=round(th_free, 1),
                free_azimuth=bool(th >= th_free))

# ----------------------------------------------------------------------
# Assemble
# ----------------------------------------------------------------------
OPP = opportunity_map()
# prefetch all verification points in batches so the API is never hammered
_pre = []
for nm, w in OPP.items():
    Dm = 0.5 * (w["D_km"][0] + w["D_km"][1])
    _pre += [fwd(IP_LAT, IP_LON, w["az_deg"], Dm + k * 8.0) for k in range(-12, 13)]
_pre += [fwd(IP_LAT, IP_LON, 0.0, 656.0 + k * 8.0) for k in range(-12, 13)]
_pre += [fwd(IP_LAT, IP_LON, 180.0, 198.0 + k * 8.0) for k in range(-12, 13)]
gebco_batch([(round(a, 4), round(b, 4)) for a, b in _pre])

CANDS = []
for nm, w in sorted(OPP.items(), key=lambda kv: -kv[1]["run_km"]):
    Dmid = 0.5 * (w["D_km"][0] + w["D_km"][1])
    CANDS.append(verify_aim(w["az_deg"], Dmid, label=nm, body=nm))
CANDS.append(verify_aim(0.0, 656.0, label="Lake Superior, meridian (rcs-aim study)",
                        body="Lake Superior"))
CANDS.append(verify_aim(180.0, 198.0, label="UIUC South Farms (collider aim)",
                        body="UIUC"))

TRADE = [dict(tilt_mrad=t, far_D_km=round(2 * R_E * t * 1e-3, 0),
              **near_side(t, D0_DEFAULT, 2632.0)) for t in
         (15.4, 20, 30, 40, 51.8, 60, 63.6, 70, 80, 100)]

RINGSPEC = (("collider_C11.0km", 700.0, 1528.0),
            ("rcs12_C6.28km", 500.0, 841.0),
            ("rcs34_C14.72km", 450.0, 2200.0))
FITS = {}
for _nm, _Ls, _R in RINGSPEC:
    _cl = {a: round(fit_at(a, _Ls, _R), 0) for a in [i * 6.0 for i in range(30)]}
    FITS[_nm] = dict(best_clearance_m_by_azimuth=_cl,
                     max_clearance_m=max(_cl.values()),
                     azimuths_with_120m=[a for a, c in _cl.items() if c >= 120],
                     azimuths_with_0m=[a for a, c in _cl.items() if c >= 0])

report = dict(
    premise="the corridor meridian follows the ComEd ROW only because the "
            "collider's 15.40 mrad north beam needs ~9.5 km of institutional "
            "land to climb through the private-airspace band. Off-site "
            "sub-500 ft band ~ 1/theta, so a steep beam needs no corridor in "
            "any bearing and its azimuth is free.",
    theta_free_mrad_at_d0_15m=round((15.0 + 500 * FT) / 2632.0 * 1e3, 1),
    trade=TRADE,
    opportunity_map=OPP,
    verified_candidates=CANDS,
    ring_fit_vs_azimuth=FITS,
    caveats=[
        "spherical geodesy; a real aim needs a geodetic (WGS84) solution",
        "NE polygons GENERATE candidates only; GEBCO verifies. On the "
        "meridian NE's Lake Superior polygon includes 400 m of the Keweenaw",
        "ring fit uses a conservative rotated bounding rectangle of the "
        "racetrack against the real site boundary",
        "land use under the near-side sub-500 ft band is UNVERIFIED for any "
        "bearing other than the studied north (ComEd) and south corridors",
    ])
os.makedirs(os.path.join(ROOT, "static", "geo"), exist_ok=True)
json.dump(report, open(os.path.join(ROOT, "static", "geo", "rcs_azimuth.json"), "w"),
          indent=1)
json.dump(CACHE, open(CACHE_FN, "w"))

print("theta_free = %.1f mrad (d0 = 15 m, fence 2.63 km) -> D = %.0f km\n"
      % (report["theta_free_mrad_at_d0_15m"],
         2 * R_E * report["theta_free_mrad_at_d0_15m"] * 1e-3))
print("=== opportunity map: water bodies in (azimuth, tilt) space ===")
print("%-20s %-10s %-18s %-16s %s" % ("body", "best az", "D window km",
                                        "tilt mrad", "run"))
for nm, w in sorted(OPP.items(), key=lambda kv: -kv[1]["run_km"]):
    print("%-20s %-10s %-18s %-16s %.0f km"
          % (nm, "%.0f deg" % w["az_deg"],
             "%.0f - %.0f" % (w["D_km"][0], w["D_km"][1]),
             "%.1f - %.1f" % tuple(w["tilt_mrad"]), w["run_km"]))
print("\n=== GEBCO-verified candidate aims ===")
for c in CANDS:
    ns = c["near_side"]
    print(" %-42s az %5.1f  D %6.1f km  tilt %5.1f mrad" %
          (c["label"], c["azimuth_deg"], c["D_km"], c["tilt_mrad"]))
    print("    aim %.4f, %.4f  GEBCO %s  %s  submerged run %.0f km  max column %s"
          % (c["aim_point"][0], c["aim_point"][1],
             ("%+d m" % c["gebco_masl"]) if c["gebco_masl"] is not None else "NA",
             "WATER" if c["is_water"] else "land", c["submerged_run_km"],
             ("%.0f m" % c["max_column_in_run_m"])
             if c["max_column_in_run_m"] is not None else "n/a"))
    print("    evidence: %s" % c["evidence"])
    print("    near side on bearing %5.1f: fence %4.0f m, emerges %.2f km (%s), "
          "500 ft at %.2f km, off-site band %.2f km -> azimuth %s (theta_free %.1f)"
          % (c["near_side_bearing_deg"], c["fence_distance_m"], ns["emergence_km"],
             "on site" if ns["on_site"] else "OFF SITE", ns["clears_500ft_km"],
             ns["offsite_sub500ft_km"], "FREE" if c["free_azimuth"] else "corridor-bound",
             c["theta_free_mrad"]))
print("\n=== ring fit vs azimuth (best over ALL placements, mod 180 deg) ===")
print(" %-20s %-9s %s" % ("ring", "max clr", "clearance (m) at az 0,6,...,174"))
for k, v in FITS.items():
    row = " ".join("%4.0f" % v["best_clearance_m_by_azimuth"][a]
                   for a in [i * 6.0 for i in range(0, 30, 2)])
    print(" %-20s %6.0f m  %s" % (k, v["max_clearance_m"], row))
    ok = v["azimuths_with_120m"]
    print("   %-18s fits with 120 m margin at %d/30 azimuths%s"
          % ("", len(ok), (": %.0f-%.0f deg" % (min(ok), max(ok))) if ok else ""))
print("\nWrote static/geo/rcs_azimuth.json")
