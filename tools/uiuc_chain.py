#!/usr/bin/env python3
"""Design variant: the co-tilted chain -- every straight plume to UIUC.

The baseline pitches the RCS planes DOWN-to-north to converge their pencils
on the deep-reference hall at 9.25 km.  This variant co-tilts every ring
UP-to-north, like the civic-envelope collider, so that every east-straight
SOUTH pencil exits at the collider's own UIUC South Farms exit -- putting
the full 63 GeV - 5 TeV chirp, with its (L/E)^2-enhanced oscillation
probability, onto the far hall.

Each ring is a planar racetrack: two N-S straights joined by semicircular
arcs, so the plane spans y_c +- (Ls/2 + R) north-south and a tilt theta
raises the NORTH arc apex by (Ls/2 + R) * theta above the straight.  Per
ring, (z0, theta) must satisfy three constraints simultaneously:

  A. the south pencil passes through the common UIUC exit point;
  B. the north pencil emerges INSIDE the Fermilab fence (41.8699 N);
  C. the north arc apex keeps a minimum rock/drift cover.

A gives a one-parameter family z0(theta); B caps z0 from above (shallow
straights emerge sooner); C floors it (shallow straights lift the apex
out of the ground).  The window closes tightly for the big RCS3/4 ring
and comfortably for RCS1/2.  Straights sit SIDE-BY-SIDE with the IP
straight (10-25 m E-W offsets), not stacked, so tunnels never intersect
and every emergence lands in one on-site zone.

Outputs: static/geo/uiuc_chain.json, and the far-hall payoff series in
static/geo/uiuc_chain.npz (figures via tools/uiuc_chain_figs.py).
"""
import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import chain_timing as ct           # per-turn bursts + samplers (runs its MC)
cl = ct.cl

R_E = 6371000.0
MMU = ct.MMU
RNG = np.random.default_rng(20260824)

# ----------------------------------------------------------------------
# Terrain along the meridian: corridor_layout's profile north of 41.79,
# the interactive tool's SOUTH profile (40.16-41.78) below it.
# ----------------------------------------------------------------------
_tool = open(os.path.join(ROOT, "static", "tool", "index.html")).read()
_m = re.search(r"const PROFILE_M = (\[\[.*?\]\]);", _tool)
SOUTH_PROF = json.loads(_m.group(1))          # [lat, elev_m], 40.16..41.78

def terrain(lat):
    if lat >= 41.79:
        return cl.elev_at_lat(lat)
    for (la1, e1), (la2, e2) in zip(SOUTH_PROF, SOUTH_PROF[1:]):
        if la1 <= lat <= la2:
            t = (lat - la1) / (la2 - la1)
            return e1 + t * (e2 - e1)
    return SOUTH_PROF[0][1]

IP_LAT = 41.8443
Y = lambda lat: (lat - IP_LAT) * cl.MLAT      # metres north of the IP
LAT = lambda y: IP_LAT + y / cl.MLAT
FENCE_N = 41.8699
Z_IP = 191.0                                  # m ASL, collider straightaway
TH_COLL = 15.40e-3                            # baseline tilt, up-north

def beam_z(z0, th, y0, y):
    """Beam elevation (MSL-equivalent): straight line in the local tangent
    frame at y0, plus the sphere correction (tangent plane sits (y-y0)^2/2R
    above the sphere)."""
    d = y - y0
    return z0 + th * d + d * d / (2 * R_E)

def surface_crossing(z0, th, y0, y_lo, y_hi):
    f = lambda y: beam_z(z0, th, y0, y) - terrain(LAT(y))
    a, b = y_lo, y_hi
    fa = f(a)
    for _ in range(200):
        m = 0.5 * (a + b)
        if (f(m) > 0) == (fa > 0):
            a = m
            fa = f(m)
        else:
            b = m
    return 0.5 * (a + b)

# ----------------------------------------------------------------------
# 1. The common exit: the baseline collider's own south emergence
# ----------------------------------------------------------------------
y_exit = surface_crossing(Z_IP, TH_COLL, 0.0, -215000, -150000)
z_exit = terrain(LAT(y_exit))
y_em_coll = surface_crossing(Z_IP, TH_COLL, 0.0, 500, 12000)

# ----------------------------------------------------------------------
# 2. Per-ring solve.  Family A: the ring's straight centre (y_i, z0) with
# tilt th must put its beam line through (y_exit, z_exit):
#     z0(th) = z_exit - th*(y_exit - y_i) - (y_exit - y_i)^2/(2 R_E)
# Then scan th for the window where B and C both hold.
# ----------------------------------------------------------------------
COVER_MIN = 6.0        # m of cover required over the north-arc apex
RINGS = [
    # name,     y_i,   Ls,    R_arc
    ("RCS1/2",  0.0,   cl.LS_RCS12, cl.R_RCS12),
    ("RCS3/4", -390.0, cl.LS_RCS34, cl.R_RCS34),
]

def solve_ring(y_i, Ls, R_arc, th):
    d = y_exit - y_i
    z0 = z_exit - th * d - d * d / (2 * R_E)
    y_apex = y_i + Ls / 2 + R_arc
    z_apex = beam_z(z0, th, y_i, y_apex)          # plane elevation at apex
    cover = terrain(LAT(y_apex)) - z_apex
    y_em = surface_crossing(z0, th, y_i, y_i + 200, y_i + 15000)
    depth_ip = terrain(LAT(y_i)) - z0
    # near-hall crossing (detector lat 41.8560)
    y_nh = Y(41.8560)
    d_nh = terrain(41.8560) - beam_z(z0, th, y_i, y_nh)
    return dict(theta_mrad=th * 1e3, z0_masl=z0, depth_at_straight_m=depth_ip,
                north_emergence_lat=LAT(y_em),
                inside_fence_m=(Y(FENCE_N) - y_em),
                apex_lat=LAT(y_apex), apex_cover_m=cover,
                near_hall_depth_m=d_nh)

solutions, windows = {}, {}
for name, y_i, Ls, R_arc in RINGS:
    ths = np.arange(14.80e-3, 15.60e-3, 0.005e-3)
    ok = []
    for th in ths:
        s = solve_ring(y_i, Ls, R_arc, th)
        if s["inside_fence_m"] > 0 and s["apex_cover_m"] >= COVER_MIN:
            ok.append((th, s))
    windows[name] = (ok[0][0] * 1e3, ok[-1][0] * 1e3) if ok else None
    if ok:
        # pick the window midpoint (balanced margins)
        th_mid = 0.5 * (ok[0][0] + ok[-1][0])
        solutions[name] = solve_ring(y_i, Ls, R_arc, th_mid)

coll = dict(theta_mrad=TH_COLL * 1e3, z0_masl=Z_IP,
            depth_at_straight_m=terrain(IP_LAT) - Z_IP,
            north_emergence_lat=LAT(y_em_coll),
            inside_fence_m=Y(FENCE_N) - y_em_coll,
            apex_lat=LAT(0 + cl.LS_COLLIDER / 2 + cl.R_COLL),
            apex_cover_m=terrain(LAT(cl.LS_COLLIDER / 2 + cl.R_COLL))
            - beam_z(Z_IP, TH_COLL, 0, cl.LS_COLLIDER / 2 + cl.R_COLL),
            near_hall_depth_m=terrain(41.8560)
            - beam_z(Z_IP, TH_COLL, 0, Y(41.8560)))

# ----------------------------------------------------------------------
# 3. Far-hall payoff: baseline (store only) vs co-tilted chain
# (chirp + store), on the far-hall time axis. Includes a tau-threshold
# factor for nutau CC rates (approximate DIS suppression).
# ----------------------------------------------------------------------
L_FAR = abs(y_exit)
_TAU_R = [(3.6, 0.0), (5, 0.05), (10, 0.18), (20, 0.34), (50, 0.52),
          (100, 0.62), (300, 0.74), (1000, 0.85), (5000, 0.95)]

def tau_thresh(E):
    E = np.asarray(E, dtype=float)
    xs = np.log([e for e, r in _TAU_R])
    ys = [r for e, r in _TAU_R]
    return np.interp(np.log(np.maximum(E, 3.61)), xs, ys, left=0.0, right=0.95)

def posc(E, Lkm):
    return 0.95 * np.sin(1.267 * 2.5e-3 * Lkm / np.maximum(E, 1e-3)) ** 2

T_EDGES = ct.T_EDGES
NB = len(T_EDGES) - 1
R_PROBE_FAR = 2.0
AREA = math.pi * (R_PROBE_FAR * 100) ** 2

def far_series(components):
    flux = np.zeros(NB); nutau = np.zeros(NB)
    ntcc = 0.0; numucc = 0.0; pancake = 0.0
    cache = {}
    for blist, L, sig in components:
        for tt, E, nd in blist:
            i = int(np.searchsorted(T_EDGES, tt)) - 1
            key = (round(E, 1), round(sig * 1e6), round(L))
            if key not in cache:
                n = 1500000 if E >= 4999 else 150000
                g = E / MMU
                acc_f = acc_p = acc_ptot = acc_tcc = acc_mcc = 0.0
                for kind in ("numu", "nubare"):
                    x = ct.sample_michel(n, kind)
                    cth = 2 * RNG.random(n) - 1
                    En = g * (x * MMU / 2) * (1 + cth)
                    th = np.arctan2(np.sqrt(1 - cth ** 2), g * (1 + cth))
                    phi = 2 * np.pi * RNG.random(n)
                    tx = th * np.cos(phi) + RNG.normal(0, sig, n)
                    ty = th * np.sin(phi) + RNG.normal(0, sig, n)
                    ins = np.hypot(tx, ty) < R_PROBE_FAR / L
                    Ei = En[ins]
                    acc_f += ins.sum() / n
                    if kind == "numu":
                        P = posc(Ei, L / 1000)
                        acc_p += P.sum() / n
                        # whole-plane count, restricted to tau-CC-able E>5 GeV
                        hi = En > 5.0
                        acc_ptot += posc(En[hi], L / 1000).sum() / n
                        # nutau CC in a thin on-axis target (CSMS nu * thresh)
                        s_t = np.array([cl.sigma_csms(cl._SIG_NU, max(e, 10.0))
                                        for e in Ei]) * tau_thresh(Ei) * P
                        acc_tcc += s_t.sum() / n
                        s_m = np.array([cl.sigma_csms(cl._SIG_NU, max(e, 10.0))
                                        for e in Ei])
                        acc_mcc += s_m.sum() / n * 0.72
                cache[key] = (acc_f, acc_p, acc_ptot, acc_tcc, acc_mcc)
            acc_f, acc_p, acc_ptot, acc_tcc, acc_mcc = cache[key]
            if 0 <= i < NB:
                flux[i] += nd * acc_f
                nutau[i] += nd * acc_p
            ntcc += nd * acc_tcc
            numucc += nd * acc_mcc
            pancake += nd * acc_ptot
    w_ms = np.diff(T_EDGES) * 1e3
    return dict(flux=flux / AREA / w_ms, nutau=nutau / AREA / w_ms,
                nutau_cc_per_t_yr=ntcc / AREA * ct.CYCLES_YR * 6.022e29,
                numu_cc_per_t_yr=numucc / AREA * ct.CYCLES_YR * 6.022e29,
                pancake_per_yr=pancake * ct.CYCLES_YR)

print("far-hall series: baseline (store only) ...", flush=True)
store = [(ct.bursts["collider"], L_FAR, cl.SIG_THETA)]
BASE = far_series(store)
print("far-hall series: co-tilted chain (chirp + store) ...", flush=True)
chain = [(ct.bursts[nm], L_FAR + abs(yi), ct.SIG_DIV_RCS)
         for nm, yi in (("RCS1", 0.0), ("RCS2", 0.0),
                        ("RCS3", 390.0), ("RCS4", 390.0))]
VAR = far_series(store + chain)

# ----------------------------------------------------------------------
# 4. Exit-strip dose budget at UIUC and at the on-site north zone
# ----------------------------------------------------------------------
N_DEC_STAGE = {"RCS1": 7.4e17, "RCS2": 6.8e17, "RCS3": 1.86e17, "RCS4": 1.73e17}
E_TOP = {"RCS1": 0.314, "RCS2": 0.75, "RCS3": 1.5, "RCS4": 5.0}

def king_raw(ndec, E_tev, L_m):
    return 1.1e-18 * ndec * E_tev ** 4 / (L_m / 1000) ** 2

dose_uiuc = {nm: king_raw(N_DEC_STAGE[nm], E_TOP[nm], L_FAR) * 1e3
             for nm in N_DEC_STAGE}                       # mSv/yr, raw pencil
L_north = {"RCS1": 2200, "RCS2": 2200, "RCS3": 2800, "RCS4": 2800}
dose_north = {nm: king_raw(N_DEC_STAGE[nm], E_TOP[nm], L_north[nm])
              for nm in N_DEC_STAGE}                      # Sv/yr, on-site

out = {
    "concept": "co-tilted chain: every east-straight south pencil exits at "
               "the collider's UIUC South Farms point; straights side-by-side "
               "with the IP straight (10-25 m E-W offsets)",
    "common_exit": {"lat": round(LAT(y_exit), 4), "range_km": round(-y_exit / 1e3, 1),
                    "elev_masl": round(z_exit, 1)},
    "collider_baseline": {k: (round(v, 4) if isinstance(v, float) else v)
                          for k, v in coll.items()},
    "rings": {nm: {k: (round(v, 4) if isinstance(v, float) else v)
                   for k, v in s.items()} for nm, s in solutions.items()},
    "feasible_tilt_windows_mrad": {nm: ([round(w[0], 3), round(w[1], 3)] if w else None)
                                   for nm, w in windows.items()},
    "apex_cover_min_m": COVER_MIN,
    "far_hall_payoff": {
        "flux_on_axis_ratio_variant_over_baseline":
            round(float(np.nansum(VAR["flux"] * np.diff(T_EDGES)) /
                        np.nansum(BASE["flux"] * np.diff(T_EDGES))), 3),
        "nutau_per_cm2_cycle": {"baseline": float(np.nansum(BASE["nutau"] * np.diff(T_EDGES) * 1e3)),
                                 "variant": float(np.nansum(VAR["nutau"] * np.diff(T_EDGES) * 1e3))},
        "nutau_cc_per_t_yr": {"baseline": BASE["nutau_cc_per_t_yr"],
                               "variant": VAR["nutau_cc_per_t_yr"]},
        "numu_cc_per_t_yr": {"baseline": BASE["numu_cc_per_t_yr"],
                              "variant": VAR["numu_cc_per_t_yr"]},
        "whole_plane_oscillated_nutau_per_yr_Egt5GeV": {
            "baseline": BASE["pancake_per_yr"], "variant": VAR["pancake_per_yr"]},
        "tau_threshold_model": _TAU_R,
    },
    "uiuc_exit_dose_raw_mSv_yr": {k: round(v, 3) for k, v in dose_uiuc.items()},
    "north_onsite_strip_dose_raw_Sv_yr": {k: round(v, 2) for k, v in dose_north.items()},
    "notes": [
        "west straights carry no physics: dogleg+wobble them as the baseline "
        "does the collider utility straight (their far bands dilute ~x40)",
        "Zone C RCS exit strips of the baseline are eliminated: every chain "
        "exit is on-site (north) or inside the UIUC strip (south)",
        "open item inherited from the civic-envelope baseline and tripled "
        "here: the tilted-ring ARC disks graze at azimuth-dependent ranges "
        "(near-north azimuths exit close); needs its own study with movers",
        "exit aim sensitivity: +-0.1 mrad of plane tilt moves a south exit "
        "+-1.3 km; per-ring trim is a survey deliverable",
    ],
}
os.makedirs(os.path.join(ROOT, "static", "geo"), exist_ok=True)
with open(os.path.join(ROOT, "static", "geo", "uiuc_chain.json"), "w") as f:
    json.dump(out, f, indent=1)
np.savez_compressed(os.path.join(ROOT, "static", "geo", "uiuc_chain.npz"),
                    T_EDGES=T_EDGES, BASE_FLUX=BASE["flux"], BASE_NUTAU=BASE["nutau"],
                    VAR_FLUX=VAR["flux"], VAR_NUTAU=VAR["nutau"],
                    RING_SOLUTIONS=json.dumps(out["rings"]),
                    COLL=json.dumps(out["collider_baseline"]),
                    EXIT=json.dumps(out["common_exit"]))
print(json.dumps(out, indent=1))
