#!/usr/bin/env python3
"""Three-way far-hall comparison for the RCS-to-water reassignment.

Merges into static/geo/rcs_aim.json:
  baseline (collider store only) | store + full co-tilted chain |
  store + RCS1/2 only, i.e. RCS3/4 reassigned to the Lake Superior aim.

Reuses uiuc_chain.far_series (which includes the tau-threshold factor), so
the numbers are the same machinery as the co-tilted-chain page.
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import uiuc_chain as uc
ct, cl, L = uc.ct, uc.cl, uc.L_FAR

store = [(ct.bursts["collider"], L, cl.SIG_THETA)]
def ch(ns): return [(ct.bursts[n], L + abs(y), ct.SIG_DIV_RCS) for n, y in ns]
SETS = {
    "baseline_store_only": store,
    "store_plus_full_chain": store + ch([("RCS1", 0), ("RCS2", 0),
                                         ("RCS3", 390), ("RCS4", 390)]),
    "store_plus_RCS12_only": store + ch([("RCS1", 0), ("RCS2", 0)]),
}
fx = lambda r: float((r["flux"] * np.diff(uc.T_EDGES)).sum())
res = {}
for k, v in SETS.items():
    r = uc.far_series(v)
    res[k] = dict(flux_integral=fx(r), nutau_cc_per_t_yr=r["nutau_cc_per_t_yr"],
                  numu_cc_per_t_yr=r["numu_cc_per_t_yr"],
                  plane_nutau_Egt5_per_yr=r["pancake_per_yr"])
b, f, p = (res[k] for k in SETS)
frac = lambda a, key: (p[key] - b[key]) / (f[key] - b[key])
res["reassignment_verdict"] = {
    "nutau_cc_gain_retained_pct": round(100 * frac(p, "nutau_cc_per_t_yr"), 1),
    "plane_nutau_gain_retained_pct": round(100 * frac(p, "plane_nutau_Egt5_per_yr"), 1),
    "onaxis_flux_gain_given_up_pct": round(100 * (1 - frac(p, "flux_integral")), 1),
    "reading": "the oscillated-nutau payoff lives in the LOW-energy turns "
               "(RCS1/2), where (L/E)^2 is largest; RCS3/4 -- the high-energy "
               "ring and the variant's only new dose problem at UIUC "
               "(3.0 mSv/yr raw) -- is the one UIUC needs least.",
}
fn = os.path.join(ROOT, "static", "geo", "rcs_aim.json")
d = json.load(open(fn)); d["uiuc_chain_split"] = res
json.dump(d, open(fn, "w"), indent=1)
print(json.dumps(res, indent=1))
