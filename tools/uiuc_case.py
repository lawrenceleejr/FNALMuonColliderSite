#!/usr/bin/env python3
"""The case for a detector at UIUC.

physics_case.py optimised one figure of merit -- tau appearance significance per
tonne -- and found it scales as L^2, which sends the tau detector 700-800 km out.
That is true and this file does not unsay it.  But it is one figure of merit, and
a detector is a fixed object: a given frontal area, a given length, a given mass.
Ask instead what ONE such detector sees at each hall on the corridor, and the
answer is different:

  * Close in, the beam is a pencil.  At 2 km the top-energy pencil is 7 cm across;
    at the 9.25 km hall, 34 cm.  A kilotonne detector there is dark for every
    stage but the lowest -- its mass is not in the beam.
  * At UIUC (198 km) the smallest spot is 7 m and every stage lights the whole
    detector.  It is the FIRST hall on the corridor where that is true.
  * Beyond UIUC the tau yield of the detector is flat in L (the study's
    1/L^2 x L^2 cancellation) while everything else falls as 1/L^2.  At 830 km
    the same detector collects the same tau events and 17x fewer of everything
    else.

So UIUC is the unique hall that is both fully lit and still in the 10^10-a-year
regime, and everything that is rate-limited rather than S/B-limited -- the
charged-current cross-sections, the essentially unmeasured high-energy nu_e, the
neutrino-electron scattering that gives a leptonic sin^2 theta_W, charm, exotics
-- is an order of magnitude better there than at any far hall.  The tau programme
at UIUC is then a detector requirement, stated as a number: 1e-6 fakes per
interaction, against OPERA's 1e-4 at one-sixth the tau flight path.

Output: static/figs/uiuc_case.{svg,pdf}; static/geo/uiuc_case.json.
"""
import contextlib, io, json, math, os, sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
with contextlib.redirect_stdout(io.StringIO()):
    import physics_case as pc                     # Q per site, literature inputs, OPERA benchmarks
ns, tc = pc.ns, pc.tc
STYLE = "/root/.claude/skills/synced/tufte-plots/assets/tufte.mplstyle"
if os.path.exists(STYLE):
    plt.style.use(STYLE)
plt.rcParams["figure.dpi"] = 130

# ---------------------------------------------------------------- one detector, carried to every hall
A_DET = 100.0                                 # m^2 frontal area (10 m x 10 m)
LEN_DET = 15.0                                # m along the beam
RHO_FE = 7.8                                  # t/m^3
M_DET = A_DET * LEN_DET * RHO_FE              # 11 700 t
M_BIG = 40000.0                               # a DUNE-mass alternative, same frontal area assumed lit
N_E_PER_T = 6.022e23 * 26 / 55.85 * 1e6       # electrons per tonne of iron
SEC_YR = 3.1536e7
# neutrino-electron elastic scattering, sigma = k E: the textbook values
K_NUE = {"nubar_mu": 1.33e-42, "nu_e": 9.5e-42, "nu_mu": 1.55e-42, "nubar_e": 3.97e-42}   # cm^2 / GeV
# literature, with the two numbers fetched from INSPIRE / arXiv in this session
CHARM2 = dict(events=2677 + 2752, sin2thw=0.2324, err=0.0083,
              source="CHARM II, Phys. Lett. B 335 (1994) 246: 2677 nu_mu-e and 2752 nubar_mu-e events, sin^2 theta_W = 0.2324 +- 0.0083")
NUTEV = dict(sin2thw=0.2277, stat=0.0013, syst=0.0009,
             source="NuTeV, Phys. Rev. Lett. 88 (2002) 091802: 0.2277 +- 0.0013 +- 0.0009, three standard deviations above the SM")

def stage_meanE(nm):
    _, E1, E2, _, _ = next(s for s in tc.STAGES if s[0] == nm)
    E, w = tc._egrid(E1, E2)
    return float((w * E).sum() / w.sum())                      # mean stored-muon energy over the ramp

def spot_r50(nm, E_mu, L_km):
    """Radius holding half the flux: 1/gamma cone for the pencils, the MINT divergence for the store."""
    if nm == "Collider store":
        return 1.177 * tc.SIG_THETA * L_km * 1e3
    return L_km * 1e3 * tc.M_MU / E_mu

def a_eff(nm, E_mu, L_km):
    """Effective area of the plume: peak density x a_eff = the whole plume."""
    r = spot_r50(nm, E_mu, L_km)
    return 2 * math.pi * (r / 1.177) ** 2 if nm == "Collider store" else math.pi * r * r

def lit(nm, L_km):
    return min(1.0, a_eff(nm, stage_meanE(nm), L_km) / A_DET)

def stage_yields(nm, L, side):
    """Per tonne-year on axis for one stage: CC by flavour, tau-flavour CC, nu-e elastic."""
    _, E1, E2, n_cycle, smeared = next(s for s in tc.STAGES if s[0] == nm)
    E, w = tc._egrid(E1, E2)
    Sm = tc.POINT[nm](tc.Y); Sm = Sm / np.trapezoid(Sm, tc.Y)
    Se = (ns.g_e_point_collider if smeared else ns.f_e_axis)(tc.Y); Se = Se / np.trapezoid(Se, tc.Y)
    Enu = tc.Y[None, :] * E[:, None]
    dens = n_cycle * tc.CYCLES_YR * (E / tc.M_MU) ** 2 / (math.pi * (L * 1e5) ** 2)
    if smeared:
        dens = dens / tc.DIV_SUPP
    mu_sp, e_sp = ("nubar_mu", "nu_e") if side == "nubar" else ("nu_mu", "nubar_e")
    I = lambda S, k: float((w * dens * np.trapezoid(S[None, :] * k * Enu, tc.Y, axis=1)).sum())
    r = ns.rates(nm, L, side)
    return dict(cc_mu=r["bg_main"], cc_e=r["bg_e"], tau=r["main"] + r["minor"],
                nue_el=(I(Sm, K_NUE[mu_sp]) + I(Se, K_NUE[e_sp])) * N_E_PER_T,
                nue_el_from_e=I(Se, K_NUE[e_sp]) * N_E_PER_T)

def detector_at(L, side, mass=M_DET):
    """What a fixed detector collects per year at range L: every stage scaled by the mass it lights."""
    tot = dict(cc=0.0, cc_e=0.0, cc_mu=0.0, tau=0.0, nue_el=0.0)
    per = {}
    for nm, *_ in tc.STAGES:
        y = stage_yields(nm, L, side); f = lit(nm, L); m = mass * f
        per[nm] = dict(lit_fraction=f, cc=(y["cc_mu"] + y["cc_e"]) * m, tau=y["tau"] * m,
                       cc_e=y["cc_e"] * m, cc_mu=y["cc_mu"] * m, nue_el=y["nue_el"] * m,
                       spot_r50_m=spot_r50(nm, stage_meanE(nm), L))
        for k in tot:
            tot[k] += per[nm][k]
    tot["rate_hz"] = tot["cc"] / SEC_YR
    return tot, per

HALLS = [("on-site hall", 2.0, "nu"), ("deep hall", 9.25, "nu"), ("UIUC", 198.4, "nubar"), ("830 km", 830.0, "nubar")]
_AT = {lab: detector_at(L, side) for lab, L, side in HALLS}
AT = {lab: v[0] for lab, v in _AT.items()}
PER = {lab: v[1] for lab, v in _AT.items()}

# the same detector carried continuously along the south line
LGRID = np.geomspace(1.0, 1500.0, 160)
REF = {nm: stage_yields(nm, 198.4, "nubar") for nm, *_ in tc.STAGES}
CC_L = np.array([sum((REF[nm]["cc_mu"] + REF[nm]["cc_e"]) * (198.4 / L) ** 2 * M_DET * lit(nm, L)
                     for nm, *_ in tc.STAGES) for L in LGRID])
TAU_L = np.array([sum(REF[nm]["tau"] * M_DET * lit(nm, L) for nm, *_ in tc.STAGES) for L in LGRID])
L_FULL = float(LGRID[np.argmax(np.all([[lit(nm, L) >= 0.999 for nm, *_ in tc.STAGES] for L in LGRID], axis=1))])

# tau appearance significance against the fake rate, 5 years
F = np.geomspace(1e-7, 1e-3, 200)
Q_UIUC = pc.CHIRP["UIUC, 198 km S"]["Q_per_t_yr"]
Q_830 = pc.CHIRP["830 km S (NE Mississippi)"]["Q_per_t_yr"]
Z = {"UIUC, 40 kt": pc.EPS_TAU * np.sqrt(M_BIG * 5 * Q_UIUC / F),
     "UIUC, 11.7 kt": pc.EPS_TAU * np.sqrt(M_DET * 5 * Q_UIUC / F),
     "830 km, 10 kt": pc.EPS_TAU * np.sqrt(1e4 * 5 * Q_830 / F)}
f_for = lambda Q, M, z=5.0: (pc.EPS_TAU / z) ** 2 * M * 5 * Q          # fake rate that gives z sigma in 5 yr
REQ = {"UIUC, 40 kt": f_for(Q_UIUC, M_BIG), "UIUC, 11.7 kt": f_for(Q_UIUC, M_DET), "830 km, 10 kt": f_for(Q_830, 1e4)}

# sin^2 theta_W from nu-e scattering: per-event sensitivity |d ln sigma / d s|
s = 0.231
SENS = {"nubar_mu": abs((2 * (-0.5 + s) / 3 + 2 * s) / ((-0.5 + s) ** 2 / 3 + s ** 2)),
        "nu_e": abs((2 * (0.5 + s) + 2 * s / 3) / ((0.5 + s) ** 2 + s ** 2 / 3))}
N_NUE = AT["UIUC"]["nue_el"]
SENS_AVG = 0.86 * SENS["nu_e"] + 0.14 * SENS["nubar_mu"]
DSTAT_1YR = 1.0 / (SENS_AVG * math.sqrt(N_NUE))
DSTAT_5YR = DSTAT_1YR / math.sqrt(5)
Q2 = {nm: 2 * 0.511e-3 * pc.STAGE_E[nm] for nm, *_ in tc.STAGES}       # Q^2 <= 2 m_e E_nu

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(11.0, 8.2))
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], hspace=0.40, wspace=0.26,
                      left=0.075, right=0.975, top=0.885, bottom=0.10)
C_CC, C_TAU, C_E, C_NUE = "#1f6f8b", "#7a2e21", "#0b6e4f", "#b5541c"

# (a) one detector, carried along the corridor
axa = fig.add_subplot(gs[0, :])
axa.plot(LGRID, CC_L, color=C_CC, lw=2.0, label="all charged-current interactions per year")
axa.plot(LGRID, TAU_L, color=C_TAU, lw=2.0, label="$\\tau$-flavour CC per year")
axa.axvspan(1.0, L_FULL, color="0.5", alpha=.10, lw=0)
axa.text(math.sqrt(1.0 * L_FULL), 3e6, "the beam is a pencil:\nmost of the detector is dark", fontsize=6.8,
         color="0.35", ha="center", va="center")
axa.axvline(L_FULL, color="0.45", lw=0.8, ls=(0, (4, 2)))
axa.text(L_FULL * 0.94, 2e6, "every stage lights\nthe whole detector\nfrom %.0f km" % L_FULL, fontsize=6.6,
         color="0.4", ha="right", va="center")
for lab, L, side in HALLS:
    t = AT[lab]
    axa.plot([L, L], [t["tau"], t["cc"]], color="0.6", lw=0.6, zorder=3)
    axa.plot([L], [t["cc"]], "o", ms=6, color=C_CC, mec="w", mew=0.8, zorder=5)
    axa.plot([L], [t["tau"]], "o", ms=6, color=C_TAU, mec="w", mew=0.8, zorder=5)
    above = L > 100
    axa.annotate("%s\n%.1e / yr  =  %s" % (lab, t["cc"], ("%.0f Hz" % t["rate_hz"]) if t["rate_hz"] < 2000 else "%.0f kHz" % (t["rate_hz"] / 1e3)),
                 (L, t["cc"]), xytext=(0, 12 if above else -14), textcoords="offset points", fontsize=6.4, color=C_CC,
                 ha="center", va="bottom" if above else "top", zorder=6,
                 bbox=dict(boxstyle="round,pad=0.12", fc="w", ec="none", alpha=.85))
    axa.annotate("%.0f $\\tau$ / yr" % t["tau"], (L, t["tau"]), xytext=(0, -11), textcoords="offset points",
                 fontsize=6.4, color=C_TAU, ha="center", va="top", zorder=6)
axa.set_xscale("log"); axa.set_yscale("log")
axa.set_xlim(1.0, 1500.0); axa.set_ylim(3, 1e13)
axa.set_xlabel("range from the IP (km)", fontsize=8.2)
axa.set_ylabel("events per year in one %.1f kt detector" % (M_DET / 1e3), fontsize=7.8)
axa.set_title("(a)  the same %.0f m $\\times$ %.0f m $\\times$ %.0f m iron detector, carried along the corridor" % (10, 10, LEN_DET),
              fontsize=8.6, loc="left", pad=5)
axa.legend(fontsize=6.8, frameon=False, loc="upper right")
axa.text(0.99, 0.52, "flat: the $\\tau$ yield of a lit detector does not depend on range\n"
                     "$1/L^2$: everything else does", transform=axa.transAxes, fontsize=6.4,
         color="0.35", ha="right", va="center", style="italic")
axa.tick_params(labelsize=7.2)

# (b) what that detector sees per year, hall by hall
axb = fig.add_subplot(gs[1, 0])
keys = [("cc", "all CC", C_CC), ("cc_e", "e-flavour CC", C_E), ("nue_el", "$\\nu$$-$e elastic", C_NUE), ("tau", "$\\tau$-flavour CC", C_TAU)]
x = np.arange(len(HALLS)); wdt = 0.19
for k, (key, lab, col) in enumerate(keys):
    vals = [AT[h[0]][key] for h in HALLS]
    axb.bar(x + (k - 1.5) * wdt, vals, width=wdt, color=col, lw=0, label=lab)
axb.set_yscale("log"); axb.set_ylim(5, 5e12)
axb.set_xticks(x); axb.set_xticklabels(["%s\n%g km" % (h[0], h[1]) for h in HALLS], fontsize=6.8)
axb.set_ylabel("events per year, %.1f kt" % (M_DET / 1e3), fontsize=7.8)
axb.set_title("(b)  what it collects, hall by hall", fontsize=8.6, loc="left", pad=5)
axb.legend(fontsize=6.4, frameon=False, loc="upper right", ncol=2)
axb.axhline(CHARM2["events"], color=C_NUE, lw=0.7, ls=":")
axb.text(1.5, CHARM2["events"] * 1.35, "CHARM II, whole run", fontsize=5.9, color=C_NUE, ha="center",
         bbox=dict(boxstyle="round,pad=0.12", fc="w", ec="none", alpha=.9))
axb.axhline(19, color=C_TAU, lw=0.7, ls=":")
axb.text(1.5, 19 * 1.35, "every $\\nu_\\tau$ ever identified", fontsize=5.9, color=C_TAU, ha="center",
         bbox=dict(boxstyle="round,pad=0.12", fc="w", ec="none", alpha=.9))
axb.tick_params(labelsize=7.0)

# (c) the tau requirement, stated as a number
axc = fig.add_subplot(gs[1, 1])
for (lab, z), col, ls in zip(Z.items(), (C_TAU, C_TAU, "0.45"), ("-", (0, (4, 2)), "-")):
    axc.plot(F, z, color=col, lw=2.0 if "40" in lab or "830" in lab else 1.2, ls=ls, label=lab)
axc.axhline(5, color="0.5", lw=0.8, ls=":")
axc.text(9e-4, 5.4, "5$\\sigma$", fontsize=6.6, color="0.45", ha="right")
axc.axvline(pc.OPERA_FAKE_RATE, color="0.35", lw=0.9)
axc.text(pc.OPERA_FAKE_RATE * 0.88, 0.42, "OPERA achieved:\n17 GeV, 0.5 mm $\\tau$ flight", fontsize=6.3, color="0.35",
         va="bottom", ha="right")
for lab, col, dy in (("UIUC, 40 kt", C_TAU, 1.0), ("830 km, 10 kt", "0.45", 1.0)):
    axc.plot([REQ[lab]], [5], "o", ms=6, color=col, mec="w", mew=0.8, zorder=6)
    axc.annotate("%.0e" % REQ[lab], (REQ[lab], 5), xytext=(0, -12), textcoords="offset points", fontsize=6.4,
                 color=col, ha="center", va="top")
axc.set_xscale("log"); axc.set_yscale("log")
axc.set_xlim(1e-7, 1e-3); axc.set_ylim(0.3, 100)
axc.set_xlabel("fake $\\tau$ per interaction, at 9 % $\\tau$ efficiency", fontsize=8.2)
axc.set_ylabel("$\\bar\\nu_\\tau$ appearance significance, 5 yr", fontsize=7.8)
axc.set_title("(c)  the $\\tau$ programme at UIUC is a detector number: $10^{-6}$", fontsize=8.6, loc="left", pad=5)
axc.legend(fontsize=6.6, frameon=False, loc="lower left")
axc.tick_params(labelsize=7.2)

fig.suptitle("The case for a detector at UIUC", fontsize=12, y=0.965)
fig.text(0.5, 0.928, "closer in the beam is a pencil and a kilotonne detector is dark; farther out the $\\tau$ yield is the same "
                     "and everything else is 17$\\times$ rarer. 198 km is the first fully-lit hall and the last high-rate one.",
         ha="center", va="center", fontsize=8.2, color="0.35")
fig.text(0.5, 0.012,
         "Per-tonne on-axis rates from tools/nutau_sites.py, scaled by the fraction of the detector the plume lights "
         "($\\pi r_{50}^2$ for the 1/$\\gamma$ pencils, $2\\pi\\sigma^2$ for the divergence-smeared store, at each stage's mean energy). "
         "$\\nu$$-$e elastic: $\\sigma = k E$ with the textbook $k$ per species, on iron.\n"
         "Significance uses the chirp as a weighted fit with OPERA's 9 % $\\tau$ efficiency (physics_case.py). "
         "CHARM II PLB 335 (1994) 246: 5 429 $\\nu$$-$e events, sin$^2\\theta_W$ = 0.2324 $\\pm$ 0.0083. "
         "Near halls are on the north line ($\\mu^-$ decays): CC rates identical, $\\tau$ rates 1.6$\\times$ the south values shown.",
         ha="center", va="center", fontsize=6.6, color="0.35")
for ext in ("pdf", "svg"):
    fig.savefig(os.path.join(ROOT, "static", "figs", "uiuc_case." + ext), bbox_inches="tight")

# ---------------------------------------------------------------- report and JSON
def _c(o):
    if isinstance(o, np.ndarray): return [float(v) for v in o]
    if isinstance(o, dict): return {k: _c(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [_c(v) for v in o]
    if isinstance(o, (np.floating, np.integer)): return float(o)
    return o

print("one detector: %.0f m^2 frontal, %.0f m long, %.0f t of iron" % (A_DET, LEN_DET, M_DET))
print("every stage lights the whole detector from %.0f km" % L_FULL)
print("\n%-14s %10s %10s %10s %10s %10s %10s" % ("hall", "CC / yr", "rate", "e-flav CC", "nu-e el", "tau / yr", "RCS4 lit"))
for lab, L, side in HALLS:
    t = AT[lab]; p = PER[lab]
    print("%-14s %10.2e %7.0f Hz %10.2e %10.2e %10.0f %9.0f %%" % (lab, t["cc"], t["rate_hz"], t["cc_e"], t["nue_el"], t["tau"], 100 * p["RCS4"]["lit_fraction"]))
print("\nspot radius r50 per stage (m):")
for nm, *_ in tc.STAGES:
    print("  %-14s " % nm + "  ".join("%s %8.2f" % (h[0][:6], PER[h[0]][nm]["spot_r50_m"]) for h in HALLS))
u = AT["UIUC"]
print("\nUIUC, %.1f kt, per year: %.2e CC (%.0f Hz); e-flavour CC %.2e; nu-e elastic %.2e (CHARM II's %d every %.1f h); tau-flavour %.0f"
      % (M_DET / 1e3, u["cc"], u["rate_hz"], u["cc_e"], u["nue_el"], CHARM2["events"], CHARM2["events"] / (u["nue_el"] / 8766), u["tau"]))
print("sin^2 theta_W from nu-e: per-event sensitivity %.2f -> statistical %.1e per year, %.1e in 5 yr (CHARM II total %.4f; NuTeV %.4f)"
      % (SENS_AVG, DSTAT_1YR, DSTAT_5YR, CHARM2["err"], math.hypot(NUTEV["stat"], NUTEV["syst"])))
print("Q^2 <= 2 m_e E_nu per stage: " + ", ".join("%s %.2f GeV^2" % (nm.split()[0], Q2[nm]) for nm in Q2))
print("\ntau appearance, 5 yr, fake rate for 5 sigma: " + "; ".join("%s %.1e" % (k, v) for k, v in REQ.items()))
print("  = %.0fx / %.0fx / %.0fx better than OPERA's %.1e" % tuple([pc.OPERA_FAKE_RATE / REQ[k] for k in REQ] + [pc.OPERA_FAKE_RATE]))
print("  UIUC 40 kt: Z = %.1f at 1e-5, %.1f at 1e-6; produces %.0f nubar_tau+nu_tau per year"
      % (pc.EPS_TAU * math.sqrt(M_BIG * 5 * Q_UIUC / 1e-5), pc.EPS_TAU * math.sqrt(M_BIG * 5 * Q_UIUC / 1e-6),
         sum(REF[nm]["tau"] for nm, *_ in tc.STAGES) * M_BIG))

json.dump(_c(dict(
    detector=dict(frontal_m2=A_DET, length_m=LEN_DET, mass_t=M_DET, big_mass_t=M_BIG),
    fully_lit_from_km=L_FULL,
    halls={lab: dict(range_km=L, side=side, totals=AT[lab], per_stage=PER[lab]) for lab, L, side in HALLS},
    curve=dict(L_km=LGRID, cc_per_yr=CC_L, tau_per_yr=TAU_L),
    tau_requirement=dict(fake_rate_for_5sigma_5yr=REQ, opera_fake_rate=pc.OPERA_FAKE_RATE, Q_uiuc=Q_UIUC, Q_830=Q_830),
    sin2thw=dict(nue_events_per_yr_uiuc=N_NUE, per_event_sensitivity=SENS_AVG, stat_1yr=DSTAT_1YR, stat_5yr=DSTAT_5YR,
                 Q2_max_GeV2=Q2, charm2=CHARM2, nutev=NUTEV),
)), open(os.path.join(ROOT, "static", "geo", "uiuc_case.json"), "w"), indent=1)
print("\nwrote static/figs/uiuc_case.{svg,pdf}, static/geo/uiuc_case.json")
