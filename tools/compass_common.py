"""Shared pieces for the polar 'compass' figures: spherical geodesy about the
IP, the Fermilab fence, Natural Earth lake footprints projected into
(bearing, log-range), and a polar-axes helper with a log radial scale whose
rings are labelled in km and in the tilt that reaches that range."""
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
R_E = 6371.0088                      # km
DEPTH_CAP_M = 200.0                  # no ring in this study may sit deeper than this below grade
IP_LAT, IP_LON = 41.8443, -88.222972
FT = 0.3048

# ----------------------------------------------------------------------
def inv(lat, lon):
    """(range km, bearing deg) from the IP to (lat, lon)."""
    p1, p2 = math.radians(IP_LAT), math.radians(lat)
    dl = math.radians(lon - IP_LON)
    h = math.sin((p2 - p1) / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    D = 2 * R_E * math.asin(min(1.0, math.sqrt(h)))
    az = math.degrees(math.atan2(math.sin(dl) * math.cos(p2),
                                 math.cos(p1) * math.sin(p2)
                                 - math.sin(p1) * math.cos(p2) * math.cos(dl))) % 360
    return D, az

def fwd(az_deg, D):
    """(lat, lon) a range D km from the IP along bearing az."""
    p1, l1 = math.radians(IP_LAT), math.radians(IP_LON)
    a, dl = math.radians(az_deg), D / R_E
    p2 = math.asin(math.sin(p1) * math.cos(dl) + math.cos(p1) * math.sin(dl) * math.cos(a))
    l2 = l1 + math.atan2(math.sin(a) * math.sin(dl) * math.cos(p1),
                         math.cos(dl) - math.sin(p1) * math.sin(p2))
    return math.degrees(p2), math.degrees(l2)

# ----------------------------------------------------------------------
# Straight-chord geometry for a straight at depth d0 (m) tilted theta (rad)
# up-to-one-end, on a smooth sphere: the up-going end surfaces at s_near,
# the down-going end at s_far (both km).
def exits_smooth(theta, d0_m):
    d0 = d0_m / 1000.0
    q = math.sqrt(theta * theta + 2 * d0 / R_E)
    return R_E * (-theta + q), R_E * (theta + q)

def tilt_for_far(D_km, d0_m):
    """Tilt whose down-going end surfaces at D (km) from a d0-deep straight."""
    return D_km / (2 * R_E) - (d0_m / 1000.0) / D_km

# ----------------------------------------------------------------------
_b = json.load(open(os.path.join(ROOT, "data", "fnal_boundary.geojson")))
FENCE_LL = _b["features"][0]["geometry"]["coordinates"][0]          # (lon, lat)
MLAT = 111.132
MLON = 111.320 * math.cos(math.radians(IP_LAT))
FENCE_XY = [((lo - IP_LON) * MLON * 1000, (la - IP_LAT) * MLAT * 1000) for lo, la in FENCE_LL]

def pin(x, y, poly):
    ins = False
    for (x1, y1), (x2, y2) in zip(poly, poly[1:] + [poly[0]]):
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1:
            ins = not ins
    return ins

def fence_distance(az_deg):
    """Metres from the IP to the site boundary along a bearing."""
    a = math.radians(az_deg)
    ux, uy = math.sin(a), math.cos(a)
    prev, d = 0.0, 25.0
    while d <= 25000:
        if not pin(ux * d, uy * d, FENCE_XY):
            lo, hi = prev, d
            for _ in range(50):
                m = 0.5 * (lo + hi)
                if pin(ux * m, uy * m, FENCE_XY):
                    lo = m
                else:
                    hi = m
            return 0.5 * (lo + hi)
        prev, d = d, d + 25.0
    return 25000.0

def fence_polar():
    """The boundary as (bearing rad, range km) arrays, unwrapped for polar fill."""
    pts = [inv(la, lo) for lo, la in FENCE_LL]
    az = np.unwrap(np.radians([p[1] for p in pts]))
    return az, np.array([p[0] for p in pts])

# ----------------------------------------------------------------------
GREAT = {"Lake Superior", "Lake Michigan", "Lake Huron", "Lake Erie",
         "Lake Ontario", "Lake Winnebago", "Lake Winnipeg", "Lake Nipigon",
         "Lake of the Woods", "Lake Simcoe", "Lake St. Clair"}

def _rings(g):
    if g["type"] == "Polygon":
        return [g["coordinates"][0]]
    if g["type"] == "MultiPolygon":
        return [p[0] for p in g["coordinates"]]
    return []

def lakes_polar(d_max_km, min_vertices=12, decimate_to=400):
    """[(name, is_great, az_rad_unwrapped, D_km)] for every lake ring within range."""
    lk = json.load(open(os.path.join(ROOT, "data", "ne_10m_lakes.geojson")))
    out = []
    for f in lk["features"]:
        nm = f["properties"].get("name") or "lake"
        for r in _rings(f["geometry"]):
            r = r[::max(1, len(r) // decimate_to)]
            if len(r) < min_vertices:
                continue
            pts = [inv(c[1], c[0]) for c in r]
            Ds = np.array([p[0] for p in pts])
            if Ds.min() > d_max_km:
                continue
            out.append((nm, nm in GREAT, np.unwrap(np.radians([p[1] for p in pts])), Ds))
    return out

def draw_lakes(ax, rmap, lakes, face="#3d7ea6", alpha=.55, edge=None, lw=0.0, z=2):
    for nm, great, az, Ds in lakes:
        ax.fill(az, rmap(Ds), color=face, alpha=alpha if great else 0.7 * alpha,
                lw=lw, ec=edge, zorder=z)

# ----------------------------------------------------------------------
LANDMARKS = {                                   # name: (lat, lon, category)
    "UIUC South Farms": (40.0602, -88.2230, "university"),
    "Purdue": (40.4237, -86.9212, "university"),
    "Kettle Moraine SF": (43.60, -88.20, "public"),
    "Soudan mine (MINOS far)": (47.8232, -92.2373, "mine"),
    "Ash River (NOvA far)": (48.3794, -92.8313, "mine"),
    "SURF (DUNE far)": (44.3525, -103.7510, "mine"),
    "Chicago": (41.8781, -87.6298, "city"),
    "Argonne": (41.7183, -87.9786, "federal"),
}
LAKE_LABELS = (("Lake Michigan", 44.0, -86.5), ("Lake Superior", 48.3, -87.5),
               ("Lake Huron", 45.0, -82.4), ("Lake Erie", 42.2, -81.2),
               ("L. Ontario", 43.7, -77.9))

# ----------------------------------------------------------------------
def make_rmap(d_min):
    return lambda D: np.log10(np.maximum(np.asarray(D, float), d_min) / d_min)

def polar_axes(ax, d_min, d_max, rings, tilt_labels=True, d0_m=35.0, label_color="0.4",
               fs=6.6, ring_color="0.86", rlabel_pos=250):
    """N-up, clockwise, log radius. Returns rmap."""
    rmap = make_rmap(d_min)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_ylim(0, float(rmap(d_max)))
    ax.set_yticks([float(rmap(d)) for d in rings])
    if tilt_labels:
        labs = []
        for d in rings:
            th = tilt_for_far(d, d0_m) * 1e3
            labs.append("%g km\n%.1f mrad" % (d, th) if th > 0.05 else "%g km" % d)
    else:
        labs = ["%g km" % d for d in rings]
    ax.set_yticklabels(labs, fontsize=fs, color=label_color)
    ax.set_rlabel_position(rlabel_pos)
    ax.set_xticks(np.radians(np.arange(0, 360, 45)))
    ax.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"], fontsize=8.5)
    ax.grid(color=ring_color, lw=0.6)
    return rmap

# ----------------------------------------------------------------------
def slerp(p1, p2, t):
    """Point at fraction t along the great circle from p1 to p2 (lat, lon in degrees)."""
    la1, lo1, la2, lo2 = map(math.radians, (p1[0], p1[1], p2[0], p2[1]))
    v1 = (math.cos(la1) * math.cos(lo1), math.cos(la1) * math.sin(lo1), math.sin(la1))
    v2 = (math.cos(la2) * math.cos(lo2), math.cos(la2) * math.sin(lo2), math.sin(la2))
    om = math.acos(max(-1.0, min(1.0, sum(a * b for a, b in zip(v1, v2)))))
    s1, s2 = math.sin((1 - t) * om) / math.sin(om), math.sin(t * om) / math.sin(om)
    v = [s1 * a + s2 * b for a, b in zip(v1, v2)]
    return math.degrees(math.atan2(v[2], math.hypot(v[0], v[1]))), math.degrees(math.atan2(v[1], v[0]))

def great_circle_polar(p1, p2, n=600):
    """The surface trace of the chord p1-p2 as (bearing rad, unwrapped; range km) about the IP, plus the
    chord's depth below the surface at each point (km): depth(x) = R - sqrt(R^2 cos^2(a) + x^2), x measured
    along the chord from its midpoint, a = half the central angle."""
    pts = [slerp(p1, p2, i / n) for i in range(n + 1)]
    D = np.array([inv(la, lo)[0] for la, lo in pts])
    az = np.unwrap(np.radians([inv(la, lo)[1] for la, lo in pts]))
    la1, lo1, la2, lo2 = map(math.radians, (p1[0], p1[1], p2[0], p2[1]))
    ang = math.acos(max(-1.0, min(1.0, math.sin(la1) * math.sin(la2) + math.cos(la1) * math.cos(la2) * math.cos(lo2 - lo1))))
    half = R_E * math.sin(ang / 2)                                   # half the chord length
    x = np.linspace(-half, half, n + 1)                               # along-chord offset from the midpoint
    depth = R_E - np.sqrt((R_E * math.cos(ang / 2)) ** 2 + x ** 2)
    return az, D, depth
