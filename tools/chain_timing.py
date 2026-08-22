#!/usr/bin/env python3
"""Timing study of the corridor neutrino beam through one 5 Hz machine cycle.

Every long straight of the acceleration chain is co-aligned on the corridor
meridian, and the RCS ring planes are pitched (corridor_layout's convergence
pitches, 4.32 / 2.16 mrad down-N) so their straight plumes converge on the
same deep detector hall as the collider plume (9.25 km from the IP, 107 m
below grade).  This script models what that hall sees through one cycle:

    RCS1 ramp -> RCS2 ramp -> RCS3 ramp -> RCS4 ramp -> collider store

with per-turn neutrino bursts (the bunch passes the aligned straight once per
turn per sign; each pass is a ~5 ps pulse at the hall), muon decay along the
ramps, boosted Michel spectra per species, and each component's transverse
profile at the hall plane:

  * RCS straights are FODO straights with no final focus, so their plumes are
    TRUE 1/gamma pencils: sigma_theta ~ sqrt(eps_N/(gamma*beta)) is 3-30 urad,
    below 1/gamma at every chain energy;
  * the collider IP straight is divergence-smeared, sigma_theta = SIG_THETA
    per MINT (arXiv:2608.02718), with the energy-angle correlation washed out.

Ramp durations solve  f_dec = T ln(gf/gi) / (tau (gf-gi))  for the per-stage
decay fractions the dose model already assumes (RCS1/2: 9%, RCS3/4: 7%),
which lands them in the same class as the IMCC interim-report RCS chain
(0.3-6 ms).  The collider stores for the rest of the 0.2 s cycle.

Outputs:
  static/geo/chain_timing.json  -- per-stage timing/flux table
  static/geo/chain_maps.npz     -- the three 2D fluence maps
(figures are drawn from the .npz by tools/chain_figs.py)
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import corridor_layout as cl   # runs the layout; provides geometry + constants

TAU = cl.TAU_MU
MMU = 105.658e-3               # GeV
C_LIGHT = 299792458.0
N0 = 1.8e12                    # muons/sign injected into RCS1
CYCLE = 0.2                    # s (5 Hz)
CYCLES_YR = 5 * 1.2e7
RNG = np.random.default_rng(20260822)

# ----------------------------------------------------------------------
# Chain definition (energies from content/design.md section 2)
# ----------------------------------------------------------------------
def ramp_time(Ei, Ef, f_dec):
    gi, gf = Ei / MMU, Ef / MMU
    return f_dec * TAU * (gf - gi) / math.log(gf / gi)

STAGES = [
    # name,   ring,    C [m],      Ls [m],       E_i,    E_f [GeV], f_dec
    ("RCS1", "rcs12", cl.C_RCS12, cl.LS_RCS12,    63.0,  314.0, 0.09),
    ("RCS2", "rcs12", cl.C_RCS12, cl.LS_RCS12,   314.0,  750.0, 0.09),
    ("RCS3", "rcs34", cl.C_RCS34, cl.LS_RCS34,   750.0, 1500.0, 0.07),
    ("RCS4", "rcs34", cl.C_RCS34, cl.LS_RCS34,  1500.0, 5000.0, 0.07),
]
SIG_DIV_RCS = 20e-6            # rad: FODO-straight divergence bound (<< 1/gamma)

# Distance from each straight's centre to the deep hall, along the corridor.
Y_HALL = (cl.DETECTOR_LAT - cl.LAT0) * cl.MLAT
L_STAGE = {"rcs12": Y_HALL - cl.y_rcs12, "rcs34": Y_HALL - cl.y_rcs34,
           "coll": cl.DET_RANGE}

# ----------------------------------------------------------------------
# Per-turn burst lists: (t_in_cycle, E_mu [GeV], decays aimed north)
# ----------------------------------------------------------------------
def stage_turns(t0, N_in, C, Ls, Ei, Ef, T_ramp):
    T_rev = C / C_LIGHT
    n_turns = max(1, int(round(T_ramp / T_rev)))
    out, N = [], N_in
    for k in range(n_turns):
        E = Ei + (Ef - Ei) * (k + 0.5) / n_turns   # linear-in-time ramp
        g = E / MMU
        dec_turn = N * C / (g * C_LIGHT * TAU)
        out.append((t0 + (k + 1) * T_rev, E, dec_turn * Ls / C))
        N -= dec_turn
    return out, N

bursts, summary_rows = {}, []
t_cursor, N = 0.0, N0
for name, ring, C, Ls, Ei, Ef, f in STAGES:
    T = ramp_time(Ei, Ef, f)
    blist, N_out = stage_turns(t_cursor, N, C, Ls, Ei, Ef, T)
    bursts[name] = blist
    summary_rows.append(dict(
        stage=name, E_GeV=[Ei, Ef], ramp_ms=round(T * 1e3, 3),
        turns=len(blist), T_rev_us=round(C / C_LIGHT * 1e6, 2),
        decays_north_per_cycle=sum(b[2] for b in blist),
        survival=round(N_out / N, 4)))
    t_cursor, N = t_cursor + T, N_out

T_CHAIN = t_cursor
g5 = cl.E_MU * 1e3 / MMU
T_rev_c = cl.C_COLLIDER / C_LIGHT
store, Nc, t = [], N, T_CHAIN
while t < CYCLE:
    dec_turn = Nc * cl.C_COLLIDER / (g5 * C_LIGHT * TAU)
    store.append((t, cl.E_MU * 1e3, dec_turn * cl.LS_COLLIDER / cl.C_COLLIDER))
    Nc -= dec_turn
    t += T_rev_c
bursts["collider"] = store
summary_rows.append(dict(
    stage="collider", E_GeV=[5000.0, 5000.0],
    ramp_ms=round((CYCLE - T_CHAIN) * 1e3, 1), turns=len(store),
    T_rev_us=round(T_rev_c * 1e6, 2),
    decays_north_per_cycle=sum(b[2] for b in store),
    survival=round(Nc / N, 4)))

# ----------------------------------------------------------------------
# Decay sampling: boosted Michel spectra (unpolarised), massless daughters.
#   numu:    f(x) = 2 x^2 (3 - 2x),  nubar_e: f(x) = 12 x^2 (1 - x),
# x = 2 E*/m_mu; lab E = g E*(1 + cos th*), tan th = sin th*/(g(1 + cos th*)).
# ----------------------------------------------------------------------
def sample_michel(n, kind):
    out = np.empty(0)
    while len(out) < n:
        m = 2 * (n - len(out)) + 64
        x = RNG.random(m) ** (1.0 / 3.0)          # propose ~ x^2
        acc = (3 - 2 * x) / 3.0 if kind == "numu" else (1 - x)
        x = x[RNG.random(m) < acc]
        out = np.concatenate([out, x])
    return out[:n]

def sample_decays(E_mu, n, sig_div):
    """(E_nu, theta_x, theta_y) for n decays x 2 species; parent directions
    carry Gaussian divergence sig_div per transverse plane."""
    g = E_mu / MMU
    E_all, tx_all, ty_all = [], [], []
    for kind in ("numu", "nubare"):
        x = sample_michel(n, kind)
        cth = 2 * RNG.random(n) - 1
        Estar = x * (MMU / 2.0)
        E = g * Estar * (1 + cth)
        th = np.arctan2(np.sqrt(1 - cth ** 2), g * (1 + cth))
        phi = 2 * np.pi * RNG.random(n)
        tx = th * np.cos(phi) + RNG.normal(0, sig_div, n)
        ty = th * np.sin(phi) + RNG.normal(0, sig_div, n)
        E_all.append(E); tx_all.append(tx); ty_all.append(ty)
    return np.concatenate(E_all), np.concatenate(tx_all), np.concatenate(ty_all)

def sigma_of(E, first_half_numu):
    """CSMS CC+NC per nucleon; array E, boolean species mask."""
    s = np.array([cl.sigma_csms(cl._SIG_NU, max(e, 10.0)) for e in E])
    sb = np.array([cl.sigma_csms(cl._SIG_NUB, max(e, 10.0)) for e in E])
    return np.where(first_half_numu, s, sb)

# ----------------------------------------------------------------------
# Map 1: (time in cycle, E_nu) -- fluence through a 0.5 m-radius on-axis
# disc at the deep hall, per cycle.
# ----------------------------------------------------------------------
T_EDGES = np.logspace(np.log10(2e-5), np.log10(0.2), 241)
E_EDGES = np.logspace(np.log10(2.0), np.log10(5200.0), 161)
MAP_TE = np.zeros((len(E_EDGES) - 1, len(T_EDGES) - 1))
R_PROBE = 0.5                                       # m
N_MC_TE = 150000

def accumulate_te(blist, L, sig_div):
    area = math.pi * (R_PROBE * 100.0) ** 2         # cm^2
    groups = {}                                     # (t-bin, E rounded) -> decays
    for tt, E, nd in blist:
        i = int(np.searchsorted(T_EDGES, tt)) - 1
        if 0 <= i < len(T_EDGES) - 1:
            key = (i, round(E, 1))
            groups[key] = groups.get(key, 0.0) + nd
    for (i, E), nd in groups.items():
        Enu, tx, ty = sample_decays(E, N_MC_TE, sig_div)
        inside = np.hypot(tx, ty) < (R_PROBE / L)
        if inside.any():
            h, _ = np.histogram(Enu[inside], bins=E_EDGES)
            MAP_TE[:, i] += h * (nd / N_MC_TE) / area

print("building (t, E) map ...", flush=True)
for name, ring, C, Ls, Ei, Ef, f in STAGES:
    accumulate_te(bursts[name], L_STAGE[ring], SIG_DIV_RCS)
accumulate_te(bursts["collider"], L_STAGE["coll"], cl.SIG_THETA)

# ----------------------------------------------------------------------
# Maps 2/3: (r_perp, E_nu) fluence per year at (a) the deep hall, where all
# components converge on axis, and (b) the civic-envelope near hall (1.3 km,
# 15 m deep), where the RCS plumes -- still aimed at the deep hall -- pass
# ~50-70 m below and only the collider component is on axis.
# ----------------------------------------------------------------------
R_EDGES = np.logspace(np.log10(0.02), np.log10(80.0), 161)
MAP_ER = np.zeros((len(E_EDGES) - 1, len(R_EDGES) - 1))
MAP_ER_NEAR = np.zeros_like(MAP_ER)
RING_AREAS = math.pi * (R_EDGES[1:] ** 2 - R_EDGES[:-1] ** 2) * 1e4   # cm^2

def accumulate_er(blist, L, sig_div, the_map, offset=0.0,
                  n_groups=24, n_mc=400000):
    tot = sum(b[2] for b in blist)
    if tot <= 0:
        return
    Es = np.array([b[1] for b in blist])
    Ns = np.array([b[2] for b in blist])
    order = np.argsort(Es)
    Es, Ns = Es[order], Ns[order]
    cums = np.cumsum(Ns) / tot
    qs = np.linspace(0, 1, n_groups + 1)
    for a, b in zip(qs[:-1], qs[1:]):
        sel = (cums > a) & (cums <= b)
        if not sel.any():
            continue
        nd = Ns[sel].sum()
        E_rep = float(np.average(Es[sel], weights=Ns[sel]))
        Enu, tx, ty = sample_decays(E_rep, n_mc, sig_div)
        r = np.hypot(tx * L, ty * L + offset)
        h, _, _ = np.histogram2d(Enu, r, bins=[E_EDGES, R_EDGES])
        the_map += h * (nd / n_mc) * CYCLES_YR / RING_AREAS[None, :]

print("building (E, r) map at the deep hall ...", flush=True)
for name, ring, C, Ls, Ei, Ef, f in STAGES:
    accumulate_er(bursts[name], L_STAGE[ring], SIG_DIV_RCS, MAP_ER)
accumulate_er(bursts["collider"], L_STAGE["coll"], cl.SIG_THETA, MAP_ER)

L_NEAR = 1300.0        # m north of the IP; hall 15 m deep (civic envelope)

def near_offset(ring):
    """Vertical distance of an RCS plume axis below the near hall."""
    Lr = L_NEAR + (cl.y_coll - (cl.y_rcs12 if ring == "rcs12" else cl.y_rcs34))
    d0 = cl.DEPTH_RCS12 if ring == "rcs12" else cl.DEPTH_RCS34
    pitch = cl.PITCH_RCS12 if ring == "rcs12" else cl.PITCH_RCS34
    depth = d0 + pitch * Lr - Lr ** 2 / (2 * cl.R_EARTH)
    return depth - 15.0

print("building (E, r) map at the near hall ...", flush=True)
for name, ring, C, Ls, Ei, Ef, f in STAGES:
    Lr = L_NEAR + (cl.y_coll - (cl.y_rcs12 if ring == "rcs12" else cl.y_rcs34))
    accumulate_er(bursts[name], Lr, SIG_DIV_RCS, MAP_ER_NEAR,
                  offset=near_offset(ring))
accumulate_er(bursts["collider"], L_NEAR, cl.SIG_THETA, MAP_ER_NEAR)

# ----------------------------------------------------------------------
# Per-stage fluence and interaction rate through the 0.5 m on-axis probe
# at the deep hall (thin-target rate per tonne-year, CSMS sigma).
# ----------------------------------------------------------------------
def probe_numbers(blist, L, sig_div, n_mc=200000):
    tot = sum(b[2] for b in blist)
    if tot <= 0:
        return dict(fluence=0, rate_per_t_yr=0, mean_E_GeV=0)
    Es = np.array([b[1] for b in blist])
    Ns = np.array([b[2] for b in blist])
    E_rep = float(np.average(Es, weights=Ns))
    Enu, tx, ty = sample_decays(E_rep, n_mc, sig_div)
    inside = np.hypot(tx, ty) < (R_PROBE / L)
    if not inside.any():
        return dict(fluence=0, rate_per_t_yr=0, mean_E_GeV=0)
    area = math.pi * (R_PROBE * 100.0) ** 2
    w = tot / n_mc * CYCLES_YR / area              # fluence per accepted sample
    species = np.arange(2 * n_mc) < n_mc
    sig = sigma_of(Enu[inside], species[inside])
    return dict(fluence=float(inside.sum() * w),
                rate_per_t_yr=float(sig.sum() * w * 6.022e29),
                mean_E_GeV=float(Enu[inside].mean()))

print("probe rates ...", flush=True)
for row in summary_rows:
    ring = dict(RCS1="rcs12", RCS2="rcs12", RCS3="rcs34",
                RCS4="rcs34", collider="coll")[row["stage"]]
    sig_div = cl.SIG_THETA if row["stage"] == "collider" else SIG_DIV_RCS
    row["at_deep_hall_r50cm"] = probe_numbers(
        bursts[row["stage"]], L_STAGE[ring], sig_div)

# ----------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------
out = {
    "cycle_s": CYCLE, "rep_rate_Hz": 5,
    "chain_time_ms": round(T_CHAIN * 1e3, 2),
    "store_time_ms": round((CYCLE - T_CHAIN) * 1e3, 1),
    "chain_survival_to_collider": round(N / N0, 3),
    "store_decay_fraction": round(1 - Nc / N, 3),
    "stages": summary_rows,
    "hall": {"range_from_IP_km": round(cl.DET_RANGE / 1000, 3),
             "depth_m": round(cl.DET_DEPTH, 1),
             "stage_distances_km": {k: round(v / 1000, 3)
                                    for k, v in L_STAGE.items()}},
    "near_hall_rcs_offsets_m": {r: round(near_offset(r), 1)
                                for r in ("rcs12", "rcs34")},
    "pulse_trains_us": {"rcs12": round(cl.C_RCS12 / C_LIGHT * 1e6, 2),
                        "rcs34": round(cl.C_RCS34 / C_LIGHT * 1e6, 2),
                        "collider": round(T_rev_c * 1e6, 2)},
    "notes": ["RCS plumes are 1/gamma pencils (FODO straights); collider "
              "plume smeared by sigma_theta per arXiv:2608.02718",
              "RCS decay fractions match the dose model (9/9/7/7%); the "
              "implied chain survival 0.72 is less generous than the "
              "CHAIN_TRANSMISSION=0.90 used for collider normalisation",
              "wobble/segmentation mitigations OFF in these maps; they "
              "spread the RCS bands vertically at the hall (see safety)"],
}
with open(os.path.join(ROOT, "static", "geo", "chain_timing.json"), "w") as f:
    json.dump(out, f, indent=1)
np.savez_compressed(os.path.join(ROOT, "static", "geo", "chain_maps.npz"),
                    T_EDGES=T_EDGES, E_EDGES=E_EDGES, R_EDGES=R_EDGES,
                    MAP_TE=MAP_TE, MAP_ER=MAP_ER, MAP_ER_NEAR=MAP_ER_NEAR)
print(json.dumps(out, indent=1))
