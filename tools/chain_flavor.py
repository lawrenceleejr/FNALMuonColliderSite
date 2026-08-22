#!/usr/bin/env python3
"""Flavor-vs-time and multi-cycle extensions of the chain timing model.

Produces (static/geo/chain_flavor.npz):
  * MAP_TE_MULTI  -- (E, t) fluence map over three wall-clock 0.2 s cycles,
    with the physical interleave: the bunch injected at t = k*0.2 s stores for
    the full cycle (dumped at 85% decayed when the next bunch arrives), and
    the NEXT bunch's 7 ms chain chirp runs concurrently during the last 7 ms
    of each store.  Single-cycle figures elsewhere follow ONE bunch
    (chain first, then store) -- same physics, bunch-centric phase.
  * Per-flavor interaction-rate time series at the deep hall (9.25 km, all
    five components) and the near hall (1.3 km, collider only; the RCS
    pencils cross that plane 50-68 m below the detector):
    numu CC, nubar_e CC, NC(nu+nubar), per tonne-year per ms of cycle.
  * Oscillation flavor evolution: fluence-weighted <P(numu->nutau)> vs time
    in the cycle at both halls (vacuum, dm2_31 = 2.5e-3 eV^2, amplitude
    0.95; the nue->nutau channel adds ~5% and tracks the same shape).

North-aimed flux is one muon sign: exactly one numu and one nubar_e per
decay (CP mirror selectable by circulation direction).  The interacting
flavor SHARES are nearly time-invariant (sigma ~ E for both species), which
these series quantify; the oscillated component is the piece that moves.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import chain_timing as ct          # runs the per-turn model (and corridor_layout)
cl = ct.cl

RNG = np.random.default_rng(20260823)
E_EDGES = ct.E_EDGES
T_EDGES = ct.T_EDGES               # bunch-centric single-cycle bins (log)
R_PROBE = ct.R_PROBE               # 0.5 m on-axis probe
AREA = math.pi * (R_PROBE * 100.0) ** 2      # cm^2
N_MC = 200000
CC_FRAC_NU, CC_FRAC_NUB = 0.72, 0.70         # CSMS CC/(CC+NC), ~flat 0.1-5 TeV
DM2 = 2.5e-3                                  # eV^2
P_AMP = 0.95

def posc(E_GeV, L_km):
    return P_AMP * np.sin(1.267 * DM2 * L_km / np.maximum(E_GeV, 1e-3)) ** 2

# ----------------------------------------------------------------------
# Per-(E, sig_div, L) decay templates, cached: everything downstream is a
# linear rescaling by the number of decays, so each unique parent energy is
# simulated once.
# ----------------------------------------------------------------------
_cache = {}

def template(E_mu, sig_div, L, r_probe=R_PROBE):
    key = (round(E_mu, 1), round(sig_div * 1e6, 1), round(L), round(r_probe, 2))
    if key in _cache:
        return _cache[key]
    g = E_mu / ct.MMU
    # the 5 TeV store template repeats in every time bin of the multi-cycle
    # map, so its MC noise would print as horizontal stripes; oversample it
    n_mc = 2000000 if E_mu >= 4999 else N_MC
    out = {}
    for kind, tab, ccf in (("nu", cl._SIG_NU, CC_FRAC_NU),
                           ("nub", cl._SIG_NUB, CC_FRAC_NUB)):
        x = ct.sample_michel(n_mc, "numu" if kind == "nu" else "nubare")
        cth = 2 * RNG.random(n_mc) - 1
        E = g * (x * ct.MMU / 2.0) * (1 + cth)
        th = np.arctan2(np.sqrt(1 - cth ** 2), g * (1 + cth))
        phi = 2 * np.pi * RNG.random(n_mc)
        tx = th * np.cos(phi) + RNG.normal(0, sig_div, n_mc)
        ty = th * np.sin(phi) + RNG.normal(0, sig_div, n_mc)
        inside = np.hypot(tx, ty) < r_probe / L
        Ein = E[inside]
        sig = np.array([cl.sigma_csms(tab, max(e, 10.0)) for e in Ein])
        out[kind] = dict(
            frac=inside.mean(),
            hist=np.histogram(Ein, bins=E_EDGES)[0] / n_mc,
            sig_mean=sig.sum() / n_mc,          # <sigma * acceptance> per decay
            cc=ccf,
            p_mean=(posc(Ein, L / 1000.0).sum() / n_mc),  # <P * acceptance>
        )
    _cache[key] = out
    return out

# ----------------------------------------------------------------------
# Flavor time series on the bunch-centric single-cycle axis
# ----------------------------------------------------------------------
def flavor_series(components, t_edges, r_probe=R_PROBE):
    nb = len(t_edges) - 1
    area = math.pi * (r_probe * 100.0) ** 2          # cm^2
    r_numu_cc = np.zeros(nb); r_nube_cc = np.zeros(nb); r_nc = np.zeros(nb)
    flu_nu = np.zeros(nb); flu_all = np.zeros(nb); p_num = np.zeros(nb)
    for blist, L, sig in components:
        for tt, E, nd in blist:
            i = int(np.searchsorted(t_edges, tt)) - 1
            if not (0 <= i < nb):
                continue
            tp = template(E, sig, L, r_probe)
            # events per tonne-year in this time bin (thin target on axis)
            w = nd * ct.CYCLES_YR / area * 6.022e29
            r_numu_cc[i] += w * tp["nu"]["sig_mean"] * tp["nu"]["cc"]
            r_nube_cc[i] += w * tp["nub"]["sig_mean"] * tp["nub"]["cc"]
            r_nc[i] += w * (tp["nu"]["sig_mean"] * (1 - tp["nu"]["cc"]) +
                            tp["nub"]["sig_mean"] * (1 - tp["nub"]["cc"]))
            flu_nu[i] += nd * tp["nu"]["frac"]
            flu_all[i] += nd * (tp["nu"]["frac"] + tp["nub"]["frac"])
            p_num[i] += nd * tp["nu"]["p_mean"]      # numu(bar) -> nutau(bar)
    widths_ms = np.diff(t_edges) * 1e3
    with np.errstate(invalid="ignore", divide="ignore"):
        p_mean = np.where(flu_nu > 0, p_num / np.maximum(flu_nu, 1e-30), np.nan)
    return dict(numu_cc=r_numu_cc / widths_ms, nube_cc=r_nube_cc / widths_ms,
                nc=r_nc / widths_ms, p_osc=p_mean,
                # on-axis fluence rate and oscillated-nutau arrival rate,
                # per cm^2 per ms of cycle time (per cycle)
                flux=flu_all / area / widths_ms,
                nutau=p_num / area / widths_ms)

print("flavor series: deep hall ...", flush=True)
deep_components = [(ct.bursts[name], ct.L_STAGE[ring], ct.SIG_DIV_RCS)
                   for name, ring, C, Ls, Ei, Ef, f in ct.STAGES]
deep_components.append((ct.bursts["collider"], ct.L_STAGE["coll"], cl.SIG_THETA))
FL_DEEP = flavor_series(deep_components, T_EDGES)

print("flavor series: near hall ...", flush=True)
FL_NEAR = flavor_series([(ct.bursts["collider"], ct.L_NEAR, cl.SIG_THETA)],
                        T_EDGES)

# UIUC far hall, 197.4 km down the SOUTH beam (CP-mirror flavours:
# nubar_mu + nu_e). Collider-only: the south-going RCS plumes, pitched
# down-to-north, leave the ground 11.5 / 21 km south and pass 3.4-3.9 km
# ABOVE the surface at the exit. Wider probe (2 m << the 36 m core) for
# MC statistics; all series are per-cm^2, so probe size divides out.
L_FAR = 197400.0
print("flavor series: UIUC far hall ...", flush=True)
FL_FAR = flavor_series([(ct.bursts["collider"], L_FAR, cl.SIG_THETA)],
                       T_EDGES, r_probe=2.0)

# ----------------------------------------------------------------------
# Multi-cycle wall-clock map: three 0.2 s cycles, linear time.
# Store: bunch injected at t = k*T fills [kT, (k+1)T); chain of the NEXT
# bunch runs during the last T_CHAIN of each cycle.
# ----------------------------------------------------------------------
print("multi-cycle map ...", flush=True)
N_CYC = 3
T_CYC = ct.CYCLE
T_MULTI = np.linspace(0, N_CYC * T_CYC, 721)
MAP_MULTI = np.zeros((len(E_EDGES) - 1, len(T_MULTI) - 1))
E_CTR = np.sqrt(E_EDGES[:-1] * E_EDGES[1:])

g5 = cl.E_MU * 1e3 / ct.MMU
T_rev_c = cl.C_COLLIDER / ct.C_LIGHT
N_at_inject = ct.N0 * ct.out["chain_survival_to_collider"] \
    if hasattr(ct, "out") else ct.N0 * 0.725
store_wall = []
t, Nc = 0.0, N_at_inject
while t < T_CYC:
    dec = Nc * cl.C_COLLIDER / (g5 * ct.C_LIGHT * ct.TAU)
    store_wall.append((t, cl.E_MU * 1e3, dec * cl.LS_COLLIDER / cl.C_COLLIDER))
    Nc -= dec
    t += T_rev_c
chain_wall = [(T_CYC - ct.T_CHAIN + tt, E, nd)
              for name, ring, C, Ls, Ei, Ef, f in ct.STAGES
              for (tt, E, nd) in ct.bursts[name]]
chain_L = {}
for name, ring, C, Ls, Ei, Ef, f in ct.STAGES:
    for (tt, E, nd) in ct.bursts[name]:
        chain_L[round(E, 1)] = ct.L_STAGE[ring]

def add_multi(blist, L_of, sig):
    nb = len(T_MULTI) - 1
    for k in range(N_CYC):
        for tt, E, nd in blist:
            i = int(np.searchsorted(T_MULTI, tt + k * T_CYC)) - 1
            if not (0 <= i < nb):
                continue
            L = L_of if isinstance(L_of, float) else L_of[round(E, 1)]
            tp = template(E, sig, L)
            MAP_MULTI[:, i] += (tp["nu"]["hist"] + tp["nub"]["hist"]) * nd / AREA

add_multi(store_wall, ct.L_STAGE["coll"], cl.SIG_THETA)
add_multi(chain_wall, chain_L, ct.SIG_DIV_RCS)

# ----------------------------------------------------------------------
# Fluence-weighted mean energy profiles for overlay
# ----------------------------------------------------------------------
def mean_profile(m):
    tot = m.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(tot > 0, (m * E_CTR[:, None]).sum(axis=0) / tot, np.nan)

MEAN_MULTI = mean_profile(MAP_MULTI)

np.savez_compressed(
    os.path.join(ROOT, "static", "geo", "chain_flavor.npz"),
    T_EDGES=T_EDGES, E_EDGES=E_EDGES, T_MULTI=T_MULTI,
    MAP_MULTI=MAP_MULTI, MEAN_MULTI=MEAN_MULTI,
    DEEP_NUMU_CC=FL_DEEP["numu_cc"], DEEP_NUBE_CC=FL_DEEP["nube_cc"],
    DEEP_NC=FL_DEEP["nc"], DEEP_POSC=FL_DEEP["p_osc"],
    DEEP_FLUX=FL_DEEP["flux"], DEEP_NUTAU=FL_DEEP["nutau"],
    NEAR_NUMU_CC=FL_NEAR["numu_cc"], NEAR_NUBE_CC=FL_NEAR["nube_cc"],
    NEAR_NC=FL_NEAR["nc"], NEAR_POSC=FL_NEAR["p_osc"],
    NEAR_FLUX=FL_NEAR["flux"], NEAR_NUTAU=FL_NEAR["nutau"],
    FAR_NUMU_CC=FL_FAR["numu_cc"], FAR_NUBE_CC=FL_FAR["nube_cc"],
    FAR_NC=FL_FAR["nc"], FAR_POSC=FL_FAR["p_osc"],
    FAR_FLUX=FL_FAR["flux"], FAR_NUTAU=FL_FAR["nutau"])

tot_deep = (FL_DEEP["numu_cc"] + FL_DEEP["nube_cc"] + FL_DEEP["nc"])
mask = tot_deep > 0
shares = {k: float((FL_DEEP[k][mask] * 1).sum() / tot_deep[mask].sum())
          for k in ("numu_cc", "nube_cc", "nc")}
widths_ms = np.diff(T_EDGES) * 1e3
def per_cycle(series):
    return float(np.nansum(series * widths_ms))
print(json.dumps({
    "nutau_per_cm2_per_cycle": {"deep": per_cycle(FL_DEEP["nutau"]),
                                "near": per_cycle(FL_NEAR["nutau"]),
                                "far_UIUC": per_cycle(FL_FAR["nutau"])},
    "nutau_chirp_fraction_deep": float(
        np.nansum((FL_DEEP["nutau"] * widths_ms)[T_EDGES[:-1] < ct.T_CHAIN]) /
        per_cycle(FL_DEEP["nutau"])),
    "posc_store_far": float(np.nanmedian(FL_FAR["p_osc"][T_EDGES[:-1] > 0.02])),
    "interaction_shares_deep_cycleavg": {k: round(v, 3) for k, v in shares.items()},
    "posc_range_deep": [float(np.nanmin(FL_DEEP["p_osc"])),
                        float(np.nanmax(FL_DEEP["p_osc"]))],
    "posc_store_near": float(np.nanmedian(FL_NEAR["p_osc"])),
    "multi_map_sum_nu_cm2_3cycles": float(MAP_MULTI.sum()),
}, indent=1))
