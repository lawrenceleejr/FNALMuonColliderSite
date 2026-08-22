#!/usr/bin/env python3
"""
Corridor-aligned muon collider siting layout for Fermilab.

Places a 10 TeV muon collider ring and its RCS chain on the FNAL site such
that every long straight section lies on the meridian through a neutrino
detector sited on the ComEd Aurora-Wayne transmission right-of-way north of
the laboratory. All decay-neutrino plumes from the long straights are then
confined to a single N-S corridor which doubles as the physics beamline
(detector at ~8 km) and the controlled radiological corridor.

Computes:
  * racetrack ring placements fitted inside the real OSM site boundary
  * beam/plume trajectories including Earth curvature, ring pitch and real
    terrain (SRTM profile along the meridian)
  * detector hall depth, plume exit points north and south
  * neutrino fluxes and event rates at the detector
  * annual-dose estimates at the plume exit regions using the analytic
    model of B.J. King, arXiv:physics/9908017 (conservative "equilibrium
    approximation"), with explicitly-labelled correction factors

Pure stdlib. Outputs: static/geo/layout.geojson, static/geo/summary.json,
static/map/layout_data.js, static/figs/corridor_profile.svg, and a printed report.
"""

import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# ----------------------------------------------------------------------
# Fixed geographic inputs
# ----------------------------------------------------------------------
DETECTOR_LAT = 41.0 + 55/60 + 39.3/3600        # 41 55'39.3" N
DETECTOR_LON = -(88.0 + 13/60 + 22.7/3600)     # 88 13'22.7" W
MERIDIAN = DETECTOR_LON                        # beam azimuth = 0 deg (true N)

R_EARTH = 6.371e6                              # m

# ----------------------------------------------------------------------
# Machine / beam assumptions (IMCC-class 10 TeV design, tunable)
# ----------------------------------------------------------------------
E_MU = 5.0             # TeV per beam (sqrt(s) = 10 TeV)
GAMMA = E_MU * 1e12 / 105.658e6
N_MU_YEAR = 1.8e12 * 5 * 1.2e7      # muons/yr/sign injected (1.8e12/bunch, 5 Hz, 1.2e7 s)
CHAIN_TRANSMISSION = 0.90           # fraction of injected muons alive at collider
C_COLLIDER = 11000.0    # m  (11 km per the corridor study baseline)
LS_COLLIDER = 700.0     # m  IP insertion straight (FF + matching), also utility straight
C_RCS12 = 6283.0        # m  Tevatron-equal circumference, two stacked RCS (63->750 GeV)
LS_RCS12 = 500.0        # m
C_RCS34_TARGET = 16500.0  # m  FNAL "site filler" bound for the final RCS pair (-> 5 TeV)
DEPTH_COLLIDER = 100.0  # m below grade at the IP (NuMI/MINOS-hall depth class, in bedrock)
DEPTH_RCS34 = 80.0
DEPTH_RCS12 = 60.0
BOUNDARY_MARGIN = 150.0  # m tunnel setback from the site boundary

# Forward-plume model. MINT (arXiv:2608.02718), decaying muons along the real
# IMCC hybrid v0.6+v0.9 interaction-region lattice, finds the IP-straight plume
# is set by the muon beam divergence, sigma_theta ~ 0.1-0.2 mrad >> 1/gamma =
# 21 urad, and that the energy-angle ("prism") correlation is washed out. The
# plume is therefore modelled as two components: a divergence-smeared component
# (weight 1 - PENCIL_FRAC, effective divergence SIG_THETA) and an optional true
# 1/gamma pencil (weight PENCIL_FRAC) from a dedicated dispersion-free waisted
# drift (beta* ~ l gives sigma_theta = sqrt(eps_N/(gamma*l)) ~ 1 urad) that the
# current lattice does NOT have. PENCIL_FRAC = 0 is the honest baseline;
# PENCIL_FRAC_DED is the design ask (everything but the +-180 m FF/chicanes).
TAU_MU = 2.19698e-6     # s
SIG_THETA = 0.15e-3     # rad: effective divergence of the smeared component (MINT)
PENCIL_FRAC = 0.0       # baseline: current IMCC-class insertion, no dedicated drift
PENCIL_FRAC_DED = 340.0 / LS_COLLIDER   # dedicated-drift scenario
DIV_SUPP = 2 * GAMMA**2 * SIG_THETA**2  # on-axis density dilution of the smeared component (~101)
# Fraction of stored muons that decay within the 0.2 s store (the rest are
# dumped and make no corridor neutrinos); 1 - exp(-0.2/(gamma*tau)) = 0.854.
F_STORE = 1.0 - math.exp(-0.2 / (GAMMA * TAU_MU))
# sigma_nu propagator suppression vs the linear low-energy extrapolation,
# flux-weighted over the two species (CSMS 1106.3723). With the prism washed
# out the fluence-mean energies drop to 0.35/0.30 E_mu (half the on-axis
# values), where the suppression is milder: x0.92.
XSEC_FLATTEN = 0.92
SEG_SPREAD = 0.5e-3     # rad: optional further +-0.5 mrad vertical segmentation
WOBBLE = 1.0e-3         # rad: +-1 mrad IMCC mover system (arcs, utility and RCS straights)
SHOWER_W = 2.0          # m: transverse washout scale of hadronic/EM showers in soil

# ----------------------------------------------------------------------
# Geodesy helpers (local equirectangular frame centred on the meridian)
# ----------------------------------------------------------------------
LAT0 = 41.845

def m_per_deg_lat(lat):
    p = math.radians(lat)
    return 111132.92 - 559.82*math.cos(2*p) + 1.175*math.cos(4*p)

def m_per_deg_lon(lat):
    p = math.radians(lat)
    return 111412.84*math.cos(p) - 93.5*math.cos(3*p)

MLAT = m_per_deg_lat(LAT0)
MLON = m_per_deg_lon(LAT0)

def to_xy(lat, lon):
    """x east, y north, metres, origin at (LAT0, MERIDIAN)."""
    return (lon - MERIDIAN) * MLON, (lat - LAT0) * MLAT

def to_ll(x, y):
    return LAT0 + y / MLAT, MERIDIAN + x / MLON

# ----------------------------------------------------------------------
# Site boundary and terrain
# ----------------------------------------------------------------------
with open(os.path.join(ROOT, 'data', 'fnal_boundary.geojson')) as f:
    _b = json.load(f)
BOUNDARY_LL = _b['features'][0]['geometry']['coordinates'][0]   # [lon, lat]
BOUNDARY_XY = [to_xy(la, lo) for lo, la in BOUNDARY_LL]

with open(os.path.join(ROOT, 'data', 'elevation_profile_meridian.json')) as f:
    _e = json.load(f)
PROF = [(p['lat'], p['elev_m']) for p in _e['points']]

def elev_at_lat(lat):
    for (la1, e1), (la2, e2) in zip(PROF, PROF[1:]):
        if la1 <= lat <= la2:
            t = (lat - la1) / (la2 - la1)
            return e1 + t * (e2 - e1)
    return PROF[0][1] if lat < PROF[0][0] else PROF[-1][1]

def point_in_poly(x, y, poly):
    inside = False
    for (x1, y1), (x2, y2) in zip(poly, poly[1:]):
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1:
            inside = not inside
    return inside

def dist_to_boundary(x, y):
    best = 1e18
    for (x1, y1), (x2, y2) in zip(BOUNDARY_XY, BOUNDARY_XY[1:]):
        dx, dy = x2 - x1, y2 - y1
        t = max(0.0, min(1.0, ((x - x1)*dx + (y - y1)*dy) / (dx*dx + dy*dy)))
        best = min(best, math.hypot(x - x1 - t*dx, y - y1 - t*dy))
    return best

def fits(pts, margin):
    return all(point_in_poly(x, y, BOUNDARY_XY) and dist_to_boundary(x, y) >= margin
               for x, y in pts)

# ----------------------------------------------------------------------
# Racetrack geometry: east straight ON the meridian (x = 0), body west.
# ----------------------------------------------------------------------
def racetrack(y_c, Ls, R, n_arc=24):
    """Polyline (x, y) of a racetrack whose east straight is centred at
    (0, y_c), runs N-S, with arcs of radius R and the west straight at
    x = -2R.  Returns closed polyline."""
    y_n, y_s = y_c + Ls/2, y_c - Ls/2
    pts = [(0.0, y_s), (0.0, y_n)]
    for i in range(1, n_arc + 1):                      # north arc, E -> W
        a = math.pi * i / n_arc
        pts.append((-R + R*math.cos(a), y_n + R*math.sin(a)))
    pts.append((-2*R, y_s))                            # west straight, N -> S
    for i in range(1, n_arc + 1):                      # south arc, W -> E
        a = math.pi * i / n_arc
        pts.append((-R - R*math.cos(a), y_s - R*math.sin(a)))
    pts.append((0.0, y_s))
    return pts

def northmost_fit(Ls, R, margin, y_hi=3000.0, y_lo=-3500.0, step=10.0):
    y = y_hi
    while y > y_lo:
        if fits(racetrack(y, Ls, R), margin):
            return y
        y -= step
    return None

# Collider: fixed C, slide as far north as possible, then step 350 m south so
# BOTH straights' plume emergences (at the corridor-study tilt) fall on-site --
# the west boundary sits lower (~41.8656) than the NE corner.
R_COLL = (C_COLLIDER - 2*LS_COLLIDER) / (2*math.pi)
y_coll = northmost_fit(LS_COLLIDER, R_COLL, BOUNDARY_MARGIN)
assert y_coll is not None, "collider does not fit"
y_coll -= 350.0

# RCS3/4 "site filler": maximise C = 2*pi*R + 2*Ls subject to fit
best34 = None
R_test = 2350.0
while R_test >= 1700.0:
    Ls_test = 1300.0
    while Ls_test >= 100.0:
        y = northmost_fit(Ls_test, R_test, 100.0, step=25.0)
        if y is not None:
            C = 2*math.pi*R_test + 2*Ls_test
            if best34 is None or C > best34[0]:
                best34 = (C, R_test, Ls_test, y)
        Ls_test -= 50.0
    R_test -= 25.0
C_RCS34, R_RCS34, LS_RCS34, y_rcs34 = best34

# RCS1/2: fixed Tevatron-equal C, stack its straight over the collider straight
R_RCS12 = (C_RCS12 - 2*LS_RCS12) / (2*math.pi)
y_rcs12 = y_coll if fits(racetrack(y_coll, LS_RCS12, R_RCS12), 100.0) \
    else northmost_fit(LS_RCS12, R_RCS12, 100.0)

# Stack the RCS3/4 straight over the collider straight too, if it fits there
y_stacked = northmost_fit(LS_RCS34, R_RCS34, 100.0, y_hi=y_coll, step=10.0)
if y_stacked is not None and abs(y_stacked - y_coll) < 400:
    y_rcs34 = y_stacked

IP_Y = y_coll                       # IP at centre of the collider east straight
IP_LAT, _ = to_ll(0, IP_Y)
IP_ELEV = elev_at_lat(IP_LAT)
IP_TUNNEL_Z = IP_ELEV - DEPTH_COLLIDER    # m ASL

# ----------------------------------------------------------------------
# Plume trajectories along the meridian
# ----------------------------------------------------------------------
def plume_depth(x, d0, pitch_down_north, north=True):
    """Depth (m) of the plume centreline below local terrain at range x (m)
    from the source straight, following a straight chord on a spherical
    Earth with real terrain. pitch_down_north > 0 tilts the straight down
    toward the north."""
    s = 1.0 if north else -1.0
    lat = IP_LAT + s * x / MLAT
    dz_terrain = elev_at_lat(lat) - IP_ELEV
    return d0 + dz_terrain - x*x/(2*R_EARTH) + s * pitch_down_north * x

def exit_range(d0, pitch, north=True, x_max=150e3):
    """Range (m) at which the plume centreline reaches the surface."""
    x, step = 1000.0, 1000.0
    prev = plume_depth(x, d0, pitch, north)
    while x < x_max:
        x += step
        d = plume_depth(x, d0, pitch, north)
        if d <= 0 < prev:
            lo, hi = x - step, x
            for _ in range(40):
                mid = 0.5*(lo + hi)
                if plume_depth(mid, d0, pitch, north) > 0:
                    lo = mid
                else:
                    hi = mid
            return 0.5*(lo + hi)
        prev = d
    return None

DET_RANGE = (DETECTOR_LAT - IP_LAT) * m_per_deg_lat((DETECTOR_LAT + IP_LAT)/2)
DET_DEPTH = plume_depth(DET_RANGE, DEPTH_COLLIDER, 0.0, north=True)
DET_ELEV = elev_at_lat(DETECTOR_LAT)

# Convergence pitches so the stacked RCS straight plumes hit the same hall
def convergence_pitch(d_ring):
    return (DET_DEPTH - plume_depth(DET_RANGE, d_ring, 0.0, True)) / DET_RANGE

PITCH_RCS34 = convergence_pitch(DEPTH_RCS34)
PITCH_RCS12 = convergence_pitch(DEPTH_RCS12)

# ----------------------------------------------------------------------
# Neutrino flux and event rate at the detector
# ----------------------------------------------------------------------
N_COLL = N_MU_YEAR * CHAIN_TRANSMISSION * F_STORE  # muons/yr/sign decaying in the collider
N_DEC_NORTH = N_COLL * LS_COLLIDER / C_COLLIDER    # decays/yr aimed north (one sign)
L_CM = DET_RANGE * 100.0
FLUX_NAIVE = N_DEC_NORTH * GAMMA**2 / (math.pi * L_CM**2)  # per SPECIES: 1/gamma-pencil limit

def onaxis_factor(f):
    """On-axis density vs the 1/gamma-pencil limit for pencil fraction f."""
    return f + (1.0 - f) / DIV_SUPP

FLUX_CORE = FLUX_NAIVE * onaxis_factor(PENCIL_FRAC)        # baseline (f = 0)
FLUX_CORE_DED = FLUX_NAIVE * onaxis_factor(PENCIL_FRAC_DED)
# With the prism washed out, any point inside the smeared core sees the
# angle-integrated spectrum: <E> = 0.35 E_mu (numu) / 0.30 E_mu (nubar_e) --
# half the on-axis means a true pencil would deliver there.
E_NU_MEAN = 0.325 * E_MU * 1000.0                  # GeV, fluence-weighted blend
# CSMS (arXiv:1106.3723) CC+NC totals per nucleon, isoscalar, log-log interp.
# Each decay emits ONE numu (at <E> = 0.7 Emu on axis) AND ONE nubar_e (0.6 Emu),
# so the rate is Phi * [sigma_nu(0.7 Emu) + sigma_nubar(0.6 Emu)].
_SIG_NU  = [(1e1, 8.2e-38), (1e3, 8.2e-36), (2e3, 15.8e-36), (5e3, 35.6e-36), (1e4, 62e-36)]
_SIG_NUB = [(1e1, 4.8e-38), (1e3, 4.8e-36), (2e3, 9.4e-36), (5e3, 22.8e-36), (1e4, 42e-36)]
def sigma_csms(tab, E_gev):
    import bisect
    xs = [math.log(e) for e, s in tab]; ys = [math.log(s) for e, s in tab]
    x = math.log(E_gev)
    i = 0 if x <= xs[0] else len(xs)-2 if x >= xs[-1] else bisect.bisect(xs, x)-1
    t = (x - xs[i]) / (xs[i+1] - xs[i])
    return math.exp(ys[i] + t*(ys[i+1] - ys[i]))
SIGMA_SOFT = (sigma_csms(_SIG_NU, 0.35*E_MU*1e3) +
              sigma_csms(_SIG_NUB, 0.30*E_MU*1e3))  # cm^2 per decay, smeared spectrum
SIGMA_HARD = (sigma_csms(_SIG_NU, 0.7*E_MU*1e3) +
              sigma_csms(_SIG_NUB, 0.6*E_MU*1e3))   # cm^2 per decay, on-axis pencil spectrum

def rate_per_kg(f):
    return FLUX_NAIVE * 6.022e26 * (f * SIGMA_HARD + (1 - f) / DIV_SUPP * SIGMA_SOFT)

RATE_PER_KG = rate_per_kg(PENCIL_FRAC)             # interactions/kg/yr, both species
RATE_PER_KG_DED = rate_per_kg(PENCIL_FRAC_DED)
# Containment radii of the smeared core: r50 = 1.20 sigma_theta L and r99 =
# 3.0 sigma_theta L (Gaussian), with the kinematic 1/gamma and 9.95/gamma
# profile added in quadrature (matches a full convolution MC to ~5%).
R50_DET = DET_RANGE * math.hypot(1.0/GAMMA, 1.20*SIG_THETA)
R99_DET = DET_RANGE * math.hypot(9.95/GAMMA, 3.0*SIG_THETA)
R50_KIN_DET = DET_RANGE / GAMMA                    # dedicated-pencil core, 50% flux

# ----------------------------------------------------------------------
# Dose model (King, physics/9908017)
#   D_ss [Sv/yr] = 1.1e-18 * (l/C) * Nmu * E[TeV]^4 / L[km]^2   (eq. 10)
#   D_ave[Sv/yr] = 3.7e-23 * Nmu * E^3 / L^2                    (eq. 7, arcs)
# ----------------------------------------------------------------------
def dose_ss_raw(l_straight, L_m):
    return 1.1e-18 * (l_straight / C_COLLIDER) * N_COLL * E_MU**4 / (L_m/1000.0)**2

def dilution(spread_halfangle, L_m):
    """Peak-dose reduction when a pencil plume is spread vertically over
    +-spread_halfangle, relative to the min(cone, shower) width."""
    w_pencil = max(2*L_m/GAMMA, SHOWER_W)
    w_band = 2*spread_halfangle*L_m + w_pencil
    return w_band / w_pencil

def exit_report(pitch):
    rows = []
    for north in (True, False):
        x = exit_range(DEPTH_COLLIDER, pitch, north)
        if x is None:
            rows.append(dict(side='N' if north else 'S', exit_km=None))
            continue
        lat = IP_LAT + (1 if north else -1) * x / MLAT
        graze = abs(-x/R_EARTH + (1 if north else -1)*(-pitch))   # rad, approx
        raw = dose_ss_raw(LS_COLLIDER, x)          # all 700 m as a 1/gamma pencil
        # baseline: every decay smeared into SIG_THETA (the MINT lattice result)
        base = raw * XSEC_FLATTEN / DIV_SUPP
        # optional +-SEG_SPREAD vertical segmentation on top of the smeared plume
        w_div = 2 * 1.20 * SIG_THETA * x + SHOWER_W
        seg = base / ((2*SEG_SPREAD*x + w_div) / w_div)
        # dedicated-drift scenario: PENCIL_FRAC_DED of the decays in a sharp pencil
        ded = XSEC_FLATTEN * (dose_ss_raw(LS_COLLIDER*PENCIL_FRAC_DED, x)
                              + dose_ss_raw(LS_COLLIDER*(1-PENCIL_FRAC_DED), x) / DIV_SUPP)
        rows.append(dict(side='N' if north else 'S', exit_km=x/1000, exit_lat=lat,
                         graze_mrad=graze*1000, footprint_len_km=(w_div/graze)/1000 if graze > 0 else None,
                         core_width_m=w_div,
                         dose_raw_Sv=raw, dose_mitigated_mSv=base*1000,
                         dose_segmented_mSv=seg*1000, dose_dedicated_mSv=ded*1000))
    return rows

PITCHES = [0.0, 1e-3, 2e-3, 5e-3]
EXITS = {p: exit_report(p) for p in PITCHES}

# Collider utility straight (west): a FODO straight with no final focus is a
# TRUE 1/gamma pencil (its divergence ~ sqrt(eps_N/(gamma*beta)) << 1/gamma),
# so unlike the IP straight it gets no free divergence smearing and NEEDS the
# in-straight vertical dogleg (+-SEG_SPREAD) plus the +-1 mrad mover.
x_u = exit_range(DEPTH_COLLIDER, 0.0, True)
DOSE_UTILITY_RAW = dose_ss_raw(LS_COLLIDER, x_u)
DOSE_UTILITY = DOSE_UTILITY_RAW * XSEC_FLATTEN / dilution(SEG_SPREAD + WOBBLE, x_u)
# RCS straights: bound decays at top energy; decay fractions per ring (assumption)
RCS_DECAY_FRac = {'RCS1': 0.09, 'RCS2': 0.09, 'RCS3': 0.07, 'RCS4': 0.07}
RCS_TOP_E = {'RCS1': 0.314, 'RCS2': 0.75, 'RCS3': 1.5, 'RCS4': 4.5}  # TeV (RCS4 tops at ~4.5 in the 14.72 km racetrack)
RCS_C = {'RCS1': C_RCS12, 'RCS2': C_RCS12, 'RCS3': C_RCS34, 'RCS4': C_RCS34}
RCS_LS = {'RCS1': LS_RCS12, 'RCS2': LS_RCS12, 'RCS3': LS_RCS34, 'RCS4': LS_RCS34}
def rcs_dose_bound(name):
    n_dec = N_MU_YEAR * RCS_DECAY_FRac[name] * RCS_LS[name] / RCS_C[name]
    x = exit_range(DEPTH_RCS34 if '3' in name or '4' in name else DEPTH_RCS12,
                   PITCH_RCS34 if '3' in name or '4' in name else PITCH_RCS12, True) or 40e3
    raw = 1.1e-18 * n_dec * RCS_TOP_E[name]**4 / (x/1000.0)**2 / 1.0
    return raw * XSEC_FLATTEN / dilution(WOBBLE, x)
RCS_DOSES = {k: rcs_dose_bound(k) for k in RCS_DECAY_FRac}

# Arc "disk" dose: between the ring and ~35 km the disk is deep underground, so
# the relevant number is the in-plane dose where the disk grazes the surface
# (an annulus at the exit range, all azimuths), King eq. 7, diluted by the mover.
X_ARC_EXIT = exit_range(DEPTH_COLLIDER, 0.0, True)
DOSE_ARC_RAW = 3.7e-23 * N_COLL * E_MU**3 / (X_ARC_EXIT/1000.0)**2
DOSE_ARC_EXIT = DOSE_ARC_RAW / dilution(WOBBLE, X_ARC_EXIT)

# ----------------------------------------------------------------------
# GeoJSON output
# ----------------------------------------------------------------------
def line_feature(name, pts_xy, props=None):
    coords = [[round(to_ll(x, y)[1], 6), round(to_ll(x, y)[0], 6)] for x, y in pts_xy]
    p = {"name": name}
    p.update(props or {})
    return {"type": "Feature", "properties": p,
            "geometry": {"type": "LineString", "coordinates": coords}}

def point_feature(name, lat, lon, props=None):
    p = {"name": name}
    p.update(props or {})
    return {"type": "Feature", "properties": p,
            "geometry": {"type": "Point", "coordinates": [round(lon, 6), round(lat, 6)]}}

features = [
    _b['features'][0],
    line_feature("Collider ring (10 TeV, C=%.1f km, depth %g m)" % (C_COLLIDER/1000, DEPTH_COLLIDER),
                 racetrack(y_coll, LS_COLLIDER, R_COLL, 48), {"ring": "collider"}),
    line_feature("RCS3/RCS4 tunnel (C=%.2f km, depth %g m)" % (C_RCS34/1000, DEPTH_RCS34),
                 racetrack(y_rcs34, LS_RCS34, R_RCS34, 48), {"ring": "rcs34"}),
    line_feature("RCS1/RCS2 tunnel (C=%.2f km, depth %g m)" % (C_RCS12/1000, DEPTH_RCS12),
                 racetrack(y_rcs12, LS_RCS12, R_RCS12, 48), {"ring": "rcs12"}),
    point_feature("Interaction point (IP)", IP_LAT, MERIDIAN,
                  {"depth_m": DEPTH_COLLIDER, "elev_grade_m": round(IP_ELEV, 1)}),
    point_feature("Neutrino detector hall", DETECTOR_LAT, DETECTOR_LON,
                  {"range_from_IP_km": round(DET_RANGE/1000, 3),
                   "hall_depth_m": round(DET_DEPTH, 1),
                   "on": "ComEd Aurora-Wayne ROW"}),
    line_feature("Neutrino corridor (IP -> detector)",
                 [(0, IP_Y), (0, (DETECTOR_LAT - LAT0) * MLAT)], {"role": "physics beam"}),
]
for p in PITCHES:
    for row in EXITS[p]:
        if row['exit_km'] is None:
            continue
        features.append(point_feature(
            "Plume exit (%s, pitch %.0f mrad down-N)" % (row['side'], p*1000),
            row['exit_lat'], MERIDIAN,
            {"exit_km": round(row['exit_km'], 1),
             "dose_mitigated_mSv_yr": round(row['dose_mitigated_mSv'], 2)}))

geojson = {"type": "FeatureCollection",
           "attribution": "Site boundary (c) OpenStreetMap contributors, ODbL",
           "features": features}
os.makedirs(os.path.join(ROOT, 'static', 'geo'), exist_ok=True)
with open(os.path.join(ROOT, 'static', 'geo', 'layout.geojson'), 'w') as f:
    json.dump(geojson, f, indent=1)
os.makedirs(os.path.join(ROOT, 'static', 'map'), exist_ok=True)
with open(os.path.join(ROOT, 'static', 'map', 'layout_data.js'), 'w') as f:
    f.write("// generated by tools/corridor_layout.py\nconst LAYOUT = ")
    json.dump(geojson, f)
    f.write(";\nconst COMED = ")
    with open(os.path.join(ROOT, 'data', 'comed_corridor.geojson')) as g:
        f.write(g.read())
    f.write(";\n")

# ----------------------------------------------------------------------
# Elevation-profile SVG (corridor cross-section, 45 km S to 55 km N)
# ----------------------------------------------------------------------
def profile_svg(path):
    xs = [x for x in range(-45000, 55001, 500)]
    W, H = 1000, 420
    x0, x1 = -45000.0, 55000.0
    z0, z1 = -80.0, 320.0     # m ASL window
    def X(x): return (x - x0) / (x1 - x0) * (W - 90) + 60
    def Z(z): return H - 50 - (z - z0) / (z1 - z0) * (H - 90)
    terr, beam = [], []
    for x in xs:
        lat = IP_LAT + x / MLAT
        terr.append((X(x), Z(elev_at_lat(lat))))
        beam.append((X(x), Z(IP_TUNNEL_Z - x*x/(2*R_EARTH))))
    def pl(pts, style):
        return '<polyline fill="none" %s points="%s"/>' % (
            style, ' '.join('%.1f,%.1f' % p for p in pts))
    svg = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" font-family="Georgia,serif">' % (W, H),
           '<rect width="%d" height="%d" fill="white"/>' % (W, H),
           pl(terr, 'stroke="#7a5c3e" stroke-width="1.8"'),
           pl(beam, 'stroke="#0b6e4f" stroke-width="1.6" stroke-dasharray="7 4"'),
           ]
    # markers
    for xm, label, dy in [(0, 'IP straight (-%g m)' % DEPTH_COLLIDER, -8),
                          (DET_RANGE, 'detector hall (-%d m)' % round(DET_DEPTH), -8)]:
        lat = IP_LAT + xm / MLAT
        svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#888" stroke-width="0.8"/>' %
                   (X(xm), Z(elev_at_lat(lat)), X(xm), Z(IP_TUNNEL_Z - xm*xm/(2*R_EARTH))))
        svg.append('<text x="%.1f" y="%.1f" font-size="12" text-anchor="middle">%s</text>' %
                   (X(xm), Z(elev_at_lat(lat)) + dy, label))
    for p, col in [(0.0, '#b3261e'), (2e-3, '#c77b21')]:
        for north in (True, False):
            xe = exit_range(DEPTH_COLLIDER, p, north)
            if xe:
                xv = xe if north else -xe
                svg.append('<circle cx="%.1f" cy="%.1f" r="4" fill="%s"/>' %
                           (X(xv), Z(elev_at_lat(IP_LAT + xv/MLAT)), col))
    svg.append('<text x="60" y="24" font-size="15">Corridor elevation profile along %-.6f&#176; W '
               '(south &#8592; | &#8594; north), 10&#215; vertical exaggeration implicit</text>' % -MERIDIAN)
    svg.append('<text x="60" y="42" font-size="12" fill="#7a5c3e">terrain</text>')
    svg.append('<text x="130" y="42" font-size="12" fill="#0b6e4f">collider-straight plume centreline (pitch 0)</text>')
    svg.append('<text x="430" y="42" font-size="12" fill="#b3261e">exit points: pitch 0</text>')
    svg.append('<text x="560" y="42" font-size="12" fill="#c77b21">pitch 2 mrad down-N</text>')
    for xm in range(-40000, 50001, 10000):
        svg.append('<text x="%.1f" y="%d" font-size="11" text-anchor="middle">%+d km</text>' %
                   (X(xm), H - 28, xm//1000))
    svg.append('</svg>')
    with open(path, 'w') as f:
        f.write('\n'.join(svg))

os.makedirs(os.path.join(ROOT, 'static', 'figs'), exist_ok=True)
profile_svg(os.path.join(ROOT, 'static', 'figs', 'corridor_profile.svg'))

# ----------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------
def ll(y):
    return to_ll(0, y)[0]

summary = {
  "corridor": {
    "meridian_lon": MERIDIAN,
    "detector": {"lat": DETECTOR_LAT, "lon": DETECTOR_LON,
                 "grade_elev_m": round(DET_ELEV, 1),
                 "hall_depth_m": round(DET_DEPTH, 1),
                 "range_from_IP_km": round(DET_RANGE/1000, 3)},
    "onsite_meridian_extent_lat": [41.82065, 41.86990],
  },
  "collider": {"C_m": C_COLLIDER, "Ls_m": LS_COLLIDER, "R_arc_m": round(R_COLL, 1),
               "depth_m": DEPTH_COLLIDER,
               "IP_lat": round(IP_LAT, 6),
               "east_straight_lat": [round(ll(y_coll - LS_COLLIDER/2), 6), round(ll(y_coll + LS_COLLIDER/2), 6)],
               "west_straight_lon_offset_m": -2*round(R_COLL, 1)},
  "rcs34": {"C_m": round(C_RCS34, 1), "Ls_m": LS_RCS34, "R_arc_m": R_RCS34,
            "depth_m": DEPTH_RCS34, "straight_center_lat": round(ll(y_rcs34), 6),
            "C_target_m": C_RCS34_TARGET,
            "C_shortfall_pct": round(100*(1 - C_RCS34/C_RCS34_TARGET), 1),
            "convergence_pitch_mrad_downN": round(PITCH_RCS34*1000, 2)},
  "rcs12": {"C_m": C_RCS12, "Ls_m": LS_RCS12, "R_arc_m": round(R_RCS12, 1),
            "depth_m": DEPTH_RCS12, "straight_center_lat": round(ll(y_rcs12), 6),
            "convergence_pitch_mrad_downN": round(PITCH_RCS12*1000, 2)},
  "beam_physics": {
    "gamma": round(GAMMA), "one_over_gamma_urad": round(1e6/GAMMA, 1),
    "muons_per_year_per_sign": N_MU_YEAR,
    "store_decay_fraction": round(F_STORE, 3),
    "decays_aimed_north_per_year": N_DEC_NORTH,
    "plume_model": "two-component per MINT arXiv:2608.02718",
    "sigma_theta_mrad": SIG_THETA*1e3,
    "pencil_frac": [PENCIL_FRAC, PENCIL_FRAC_DED],
    "onaxis_density_dilution": round(DIV_SUPP, 1),
    "core_flux_nu_cm2_yr": [FLUX_CORE, FLUX_CORE_DED],
    "core_r50_r99_at_detector_m": [round(R50_DET, 2), round(R99_DET, 2)],
    "dedicated_pencil_r50_at_detector_m": round(R50_KIN_DET, 2),
    "mean_fluence_Enu_GeV": E_NU_MEAN,
    "events_per_kg_per_year": [RATE_PER_KG, RATE_PER_KG_DED],
    "events_per_tonne_per_year": [RATE_PER_KG*1000, RATE_PER_KG_DED*1000],
  },
  "dose": {
    "model": "King physics/9908017 eq.10/eq.7 equilibrium approximation; "
             "IP-straight plume smeared by SIG_THETA per arXiv:2608.02718",
    "plume": {"sigma_theta_mrad": SIG_THETA*1e3, "onaxis_dilution": round(DIV_SUPP, 1),
              "pencil_frac_scenarios": [PENCIL_FRAC, PENCIL_FRAC_DED]},
    "exits_by_pitch": {("%.0f_mrad" % (p*1000)): EXITS[p] for p in PITCHES},
    "utility_straight_raw_Sv": round(DOSE_UTILITY_RAW, 2),
    "utility_straight_dogleg_wobbled_mSv": round(DOSE_UTILITY*1000, 3),
    "rcs_straights_wobbled_mSv": {k: round(v*1000, 4) for k, v in RCS_DOSES.items()},
    "arc_exit_annulus_raw_mSv": round(DOSE_ARC_RAW*1000, 4),
    "arc_exit_annulus_wobbled_mSv": round(DOSE_ARC_EXIT*1000, 4),
    "limits_mSv_yr": {"DOE_public": 1.0, "FNAL_design_goal": 0.1},
  },
}
with open(os.path.join(ROOT, 'static', 'geo', 'summary.json'), 'w') as f:
    json.dump(summary, f, indent=1)

# ----------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------
print("=" * 74)
print("CORRIDOR-ALIGNED MUON COLLIDER LAYOUT -- FNAL")
print("=" * 74)
print("Corridor meridian (detector longitude): %.6f" % MERIDIAN)
print("Detector: %.6f N, %.6f W  grade %.0f m ASL" % (DETECTOR_LAT, -DETECTOR_LON, DET_ELEV))
print()
print("Collider  C=%.1f km  R_arc=%.0f m  straights %g m  depth %g m" %
      (C_COLLIDER/1000, R_COLL, LS_COLLIDER, DEPTH_COLLIDER))
print("  IP at lat %.6f  (east straight %.6f - %.6f)" %
      (IP_LAT, ll(y_coll - LS_COLLIDER/2), ll(y_coll + LS_COLLIDER/2)))
print("RCS3/4    C=%.2f km (target %.1f km, shortfall %.1f%%)  R=%.0f m  Ls=%.0f m"
      % (C_RCS34/1000, C_RCS34_TARGET/1000, 100*(1 - C_RCS34/C_RCS34_TARGET), R_RCS34, LS_RCS34))
print("  straight centre lat %.6f, depth %g m, convergence pitch %.2f mrad down-N"
      % (ll(y_rcs34), DEPTH_RCS34, PITCH_RCS34*1000))
print("RCS1/2    C=%.2f km  R=%.0f m  Ls=%.0f m  depth %g m  pitch %.2f mrad down-N"
      % (C_RCS12/1000, R_RCS12, LS_RCS12, DEPTH_RCS12, PITCH_RCS12*1000))
print()
print("Detector: range %.2f km from IP, hall depth %.0f m below grade" %
      (DET_RANGE/1000, DET_DEPTH))
print("  plume: sigma_theta %.2f mrad (arXiv:2608.02718), on-axis dilution /%.0f, store factor %.3f"
      % (SIG_THETA*1e3, DIV_SUPP, F_STORE))
print("  core r50/r99 %.2f/%.2f m; core flux %.2e nu/cm2/yr/species; <Enu> ~ %.0f GeV" %
      (R50_DET, R99_DET, FLUX_CORE, E_NU_MEAN))
print("  interaction rate %.2e /kg/yr  (%.2e per tonne-year)" %
      (RATE_PER_KG, RATE_PER_KG*1000))
print("  dedicated-drift scenario (f=%.2f): core r50 %.2f m, flux %.2e, %.2e per tonne-year" %
      (PENCIL_FRAC_DED, R50_KIN_DET, FLUX_CORE_DED, RATE_PER_KG_DED*1000))
print()
print("PLUME EXITS AND PEAK ANNUAL DOSE (King eq.10; corrections as labelled)")
for p in PITCHES:
    print(" pitch %.0f mrad down-N:" % (p*1000))
    for row in EXITS[p]:
        if row['exit_km'] is None:
            print("   %s: no exit within 150 km" % row['side'])
        else:
            print("   %s: exit %6.1f km (lat %.4f)  graze %.2f mrad  core %4.0f m wide x %4.1f km"
                  "  raw %7.2f Sv/yr -> smeared %7.2f mSv/yr (+seg %6.2f; dedicated-drift %7.1f)"
                  % (row['side'], row['exit_km'], row['exit_lat'], row['graze_mrad'],
                     row['core_width_m'], row['footprint_len_km'] or -1,
                     row['dose_raw_Sv'], row['dose_mitigated_mSv'],
                     row['dose_segmented_mSv'], row['dose_dedicated_mSv']))
print()
print("Wobbled sources (+-1 mrad mover):")
print("  collider utility straight (700 m, raw %.1f Sv/yr): dogleg+wobble -> %.1f mSv/yr"
      % (DOSE_UTILITY_RAW, DOSE_UTILITY*1000))
for k, v in RCS_DOSES.items():
    print("  %s straights: %.4f mSv/yr" % (k, v*1000))
print("  arc disk at surface-grazing annulus (%.0f km, all azimuths): raw %.3f -> wobbled %.4f mSv/yr"
      % (X_ARC_EXIT/1000, DOSE_ARC_RAW*1000, DOSE_ARC_EXIT*1000))
print()
print("Wrote static/geo/layout.geojson, static/geo/summary.json, static/map/layout_data.js, static/figs/corridor_profile.svg")
