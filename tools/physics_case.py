#!/usr/bin/env python3
"""The physics case: what the corridor beam is worth, against the world's samples.

Everything here is derived from the study's own machinery (tau_compass.py for the
oscillated on-axis rates, nutau_sites.py for the sign-resolved split, chain_timing.py
for the per-stage decay budget) except the literature comparison numbers, which are
collected in LIT below with their sources.

Three things the numbers have to establish, and one they have to admit:

  * SIZE.  The per-tonne on-axis rate is baseline-independent (1/L^2 flux times L^2
    appearance probability), so the sample is set by detector mass alone, and the
    unoscillated charged-current rate -- the DIS and electroweak sample -- is larger
    than every accelerator neutrino experiment ever run, by orders of magnitude.
  * SIGN.  mu+ and mu- counter-rotate, so each direction of a straight carries one
    sign: the south beam is 91 % nubar_tau, the north 96 % nu_tau.  No one has ever
    identified a nubar_tau.
  * PROVENANCE.  The flux is muon decay: the Michel spectrum, pure QED, no hadron
    production anywhere in the chain.

  * The admission: signal-to-background.  Everything terrestrial sits deep in the
    (L/E)^2 regime, so the tau appearance fraction is 1e-7 to 1e-4 and the tau sample
    hides under an unoscillated charged-current flux 10^4 to 10^7 times larger.  S/B
    scales as L^2 and as 1/E^2, so the two levers are a long baseline and a time cut
    on the low-energy turns of the chirp.  This file computes where those two put us
    relative to what OPERA actually achieved.

Output: static/figs/physics_case.{svg,pdf}; static/geo/physics_case.json.
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
    import nutau_sites as ns
tc = ns.tc
import compass_common as cc
STYLE = "/root/.claude/skills/synced/tufte-plots/assets/tufte.mplstyle"
if os.path.exists(STYLE):
    plt.style.use(STYLE)
plt.rcParams["figure.dpi"] = 130

M_TAU = 1.77686                       # GeV
CTAU_TAU = 87.11e-6                   # m
INELASTICITY = 0.6                    # <E_tau> / <E_nu> in a CC event, DIS-averaged

# ---------------------------------------------------------------- literature comparison
LIT = dict(
    DONUT=dict(label="DONUT (1997 run)", nu_tau_identified=9, interactions=None,
               charge_resolved=False, mean_E_nu_GeV=100.0,
               sigma_nutau_uncertainty=0.47,
               source="Phys. Rev. D 78, 052002 (2008): 9 nu_tau CC interactions; "
                      "sigma = (0.39 +- 0.13 +- 0.13) x 10^-38 E cm^2/GeV"),
    OPERA=dict(label="OPERA (2008-12)", nu_tau_identified=10, interactions=19505,
               background=2.0, tau_efficiency=0.09, charge_resolved=False,
               mean_E_nu_GeV=17.0,
               source="Phys. Rev. Lett. 120, 211801 (2018): 10 candidates, 2.0 expected "
                      "background, 19 505 located interactions"),
    DUNE_cpv=dict(label="DUNE, 7 yr, CP-optimised", nu_tau_identified=25 * 7,
                  mean_E_nu_GeV=8.0, charge_resolved=False,
                  source="DUNE TDR vol. II arXiv:2002.03005 sec. 8.9.1: ~25 selected/yr"),
    DUNE_tau=dict(label="DUNE, 7 yr, $\\tau$-optimised", nu_tau_identified=150 * 7,
                  mean_E_nu_GeV=20.0, charge_resolved=False,
                  source="DUNE TDR vol. II sec. 8.9.1: ~150 selected/yr with a high-energy beam"),
    NuTeV=dict(label="NuTeV", cc_events=1.62e6,
               source="Phys. Rev. Lett. 88, 091802 (2002): 1.62e6 CC events; "
                      "sin^2 theta_W = 0.2277 +- 0.0016, ~3 sigma from the SM"),
)
EPS_TAU = LIT["OPERA"]["tau_efficiency"]                       # take OPERA's own tau efficiency
OPERA_FAKE_RATE = LIT["OPERA"]["background"] / LIT["OPERA"]["interactions"]   # 1.0e-4 per interaction

# ---------------------------------------------------------------- per-stage, per-site
SITES = [("on-site hall, 2 km N", 2.0, "nu"), ("deep hall, 9.25 km N", 9.25, "nu"),
         ("UIUC, 198 km S", 198.4, "nubar"), ("400 km S", 400.0, "nubar"),
         ("655 km S (W. Tennessee)", 655.0, "nubar"), ("830 km S (NE Mississippi)", 830.0, "nubar"),
         ("655 km N (Lake Superior)", 655.0, "nu")]

def stage_mean_E(nm):
    _, E1, E2, _, _ = next(s for s in tc.STAGES if s[0] == nm)
    E, w = tc._egrid(E1, E2)
    S = tc.POINT[nm](tc.Y); S = S / np.trapezoid(S, tc.Y)
    return float((w * np.trapezoid(S[None, :] * (tc.Y[None, :] * E[:, None]), tc.Y, axis=1)).sum() / w.sum())

STAGE_E = {nm: stage_mean_E(nm) for nm, *_ in tc.STAGES}
TAU_FLIGHT = {nm: INELASTICITY * STAGE_E[nm] / M_TAU * CTAU_TAU for nm in STAGE_E}   # m

def stage_rows(L, side):
    out = {}
    for nm, *_ in tc.STAGES:
        r = ns.rates(nm, L, side)
        sig, bg = r["main"] + r["minor"], r["bg_main"] + r["bg_e"]
        nut, nub = (r["minor"], r["main"]) if side == "nubar" else (r["main"], r["minor"])
        out[nm] = dict(nutau_per_t_yr=nut, nubartau_per_t_yr=nub, tau_total_per_t_yr=sig,
                       cc_per_t_yr=bg, S_over_B=sig / bg,
                       fake_rate_needed=sig / bg * EPS_TAU,
                       vs_opera=(sig / bg * EPS_TAU) / OPERA_FAKE_RATE,
                       tau_flight_mm=TAU_FLIGHT[nm] * 1e3)
    tot_s = sum(v["tau_total_per_t_yr"] for v in out.values())
    tot_b = sum(v["cc_per_t_yr"] for v in out.values())
    out["FULL CHAIN"] = dict(tau_total_per_t_yr=tot_s, cc_per_t_yr=tot_b, S_over_B=tot_s / tot_b,
                             nutau_per_t_yr=sum(v["nutau_per_t_yr"] for v in out.values()),
                             nubartau_per_t_yr=sum(v["nubartau_per_t_yr"] for v in out.values()),
                             fake_rate_needed=tot_s / tot_b * EPS_TAU,
                             vs_opera=(tot_s / tot_b * EPS_TAU) / OPERA_FAKE_RATE, tau_flight_mm=None)
    # the low-energy window: RCS1 + RCS2 only, the turns the chirp puts first
    ws = sum(out[n]["tau_total_per_t_yr"] for n in ("RCS1", "RCS2"))
    wb = sum(out[n]["cc_per_t_yr"] for n in ("RCS1", "RCS2"))
    out["RCS1+RCS2 window"] = dict(tau_total_per_t_yr=ws, cc_per_t_yr=wb, S_over_B=ws / wb,
                                   nutau_per_t_yr=sum(out[n]["nutau_per_t_yr"] for n in ("RCS1", "RCS2")),
                                   nubartau_per_t_yr=sum(out[n]["nubartau_per_t_yr"] for n in ("RCS1", "RCS2")),
                                   fake_rate_needed=ws / wb * EPS_TAU,
                                   vs_opera=(ws / wb * EPS_TAU) / OPERA_FAKE_RATE, tau_flight_mm=None)
    return out

SITE_ROWS = {lab: dict(range_km=L, sign_side=side, stages=stage_rows(L, side)) for lab, L, side in SITES}

# ---------------------------------------------------------------- gating on the chirp
# The chain is a chirp: the 5 Hz cycle steps the stored energy from 63 GeV to 5 TeV in a known
# sequence, so a detector knows the beam energy at every instant to within a turn.  Because the
# appearance probability goes as (L/E)^2, the earliest, lowest-energy turns carry a hundred times
# the tau fraction of the store.  Gating on them is free signal-to-background.
def gate(L, side, e_max):
    """Signal, unoscillated CC and S/B per tonne-year for the turns with ramp energy below e_max."""
    sig = bg = 0.0
    for nm, E1, E2, n_cycle, smeared in tc.STAGES:
        if E1 > e_max:
            continue
        E, w = tc._egrid(E1, E2)
        m = E <= e_max
        if not m.any():
            continue
        Sm = tc.POINT[nm](tc.Y); Sm = Sm / np.trapezoid(Sm, tc.Y)
        Se = (ns.g_e_point_collider if smeared else ns.f_e_axis)(tc.Y); Se = Se / np.trapezoid(Se, tc.Y)
        Enu = tc.Y[None, :] * E[m][:, None]
        dens = n_cycle * tc.CYCLES_YR * (E[m] / tc.M_MU) ** 2 / (math.pi * (L * 1e5) ** 2)
        if smeared:
            dens = dens / tc.DIV_SUPP
        P_mu = tc.posc(Enu, L); P_e = P_mu * (ns.AMP_E / ns.AMP_MU)
        xi = tc.tau_thresh(Enu)
        s_nu, s_nub = tc.sig_cc(tc._SIG_NU, Enu), tc.sig_cc(tc._SIG_NUB, Enu)
        sig_main, sig_e = (s_nub, s_nu) if side == "nubar" else (s_nu, s_nub)
        I = lambda S, P, sg, thr=True: float(
            (w[m] * dens * np.trapezoid(S[None, :] * P * sg * (xi if thr else 1.0), tc.Y, axis=1)).sum() * tc.N_A_T)
        sig += I(Sm, P_mu, sig_main) + I(Se, P_e, sig_e)
        bg += I(Sm, np.ones_like(P_mu), sig_main, False) + I(Se, np.ones_like(P_mu), sig_e, False)
    return sig, bg

E_CUTS = [80, 100, 125, 160, 200, 250, 314, 500, 750, 1500, 5000]
GATES = {}
for lab, L, side in SITES:
    rows = []
    for ec in E_CUTS:
        s, b = gate(L, side, float(ec))
        if b <= 0:
            continue
        rows.append(dict(e_max_GeV=ec, tau_per_t_yr=s, cc_per_t_yr=b, S_over_B=s / b,
                         fake_rate_needed=s / b * EPS_TAU,
                         better_than_opera_by=OPERA_FAKE_RATE / (s / b * EPS_TAU)))
    GATES[lab] = rows

# ---------------------------------------------------------------- the chirp as a fit, not a cut
# A hard time cut throws away the high-energy turns entirely.  The chirp is better used as a
# weighted fit: every turn is a separate measurement with its own S/B, and the significance of the
# whole adds in quadrature, sum_i s_i^2 / b_i.  With a tau efficiency eps and a fake rate f per
# interaction that gives Z^2 = (eps^2 / f) sum_i S_i^2 / B_i, so the beam's figure of merit is the
# baseline- and detector-independent quantity Q = sum_i S_i^2 / B_i per tonne-year.
def chirp_Q(L, side, nbin=40):
    """Q = sum over turns of S^2/B, per tonne-year; and the turn that dominates it."""
    Q, best = 0.0, None
    for nm, E1, E2, n_cycle, smeared in tc.STAGES:
        E, w = tc._egrid(E1, E2)
        edges = np.linspace(E.min(), E.max(), nbin + 1) if E.max() > E.min() else np.array([E.min(), E.max() + 1])
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (E >= lo) & (E < hi) if hi < edges[-1] else (E >= lo) & (E <= hi)
            if not m.any():
                continue
            Sm = tc.POINT[nm](tc.Y); Sm = Sm / np.trapezoid(Sm, tc.Y)
            Se = (ns.g_e_point_collider if smeared else ns.f_e_axis)(tc.Y); Se = Se / np.trapezoid(Se, tc.Y)
            Enu = tc.Y[None, :] * E[m][:, None]
            dens = n_cycle * tc.CYCLES_YR * (E[m] / tc.M_MU) ** 2 / (math.pi * (L * 1e5) ** 2)
            if smeared:
                dens = dens / tc.DIV_SUPP
            P_mu = tc.posc(Enu, L); P_e = P_mu * (ns.AMP_E / ns.AMP_MU)
            xi = tc.tau_thresh(Enu)
            s_nu, s_nub = tc.sig_cc(tc._SIG_NU, Enu), tc.sig_cc(tc._SIG_NUB, Enu)
            sg_main, sg_e = (s_nub, s_nu) if side == "nubar" else (s_nu, s_nub)
            I = lambda S, P, sg, thr=True: float(
                (w[m] * dens * np.trapezoid(S[None, :] * P * sg * (xi if thr else 1.0), tc.Y, axis=1)).sum() * tc.N_A_T)
            Si = I(Sm, P_mu, sg_main) + I(Se, P_e, sg_e)
            Bi = I(Sm, np.ones_like(P_mu), sg_main, False) + I(Se, np.ones_like(P_mu), sg_e, False)
            if Bi <= 0:
                continue
            q = Si * Si / Bi
            Q += q
            if best is None or q > best[1]:
                best = (0.5 * (lo + hi), q, Si, Bi)
    return Q, best

def significance(Q, mass_t, years, f, eps=EPS_TAU):
    """Z for a background-limited counting fit over the chirp."""
    return eps * math.sqrt(Q * mass_t * years / f)

CHIRP = {}
for lab, L, side in SITES:
    Q, best = chirp_Q(L, side)
    CHIRP[lab] = dict(range_km=L, sign_side=side, Q_per_t_yr=Q,
                      dominant_turn_GeV=best[0], dominant_share=best[1] / Q,
                      Z_10kt_5yr_opera_rejection=significance(Q, 1e4, 5, OPERA_FAKE_RATE),
                      Z_10kt_5yr_10x_opera=significance(Q, 1e4, 5, OPERA_FAKE_RATE / 10),
                      mass_t_for_5sigma_5yr_opera=25.0 ** 1 * (5.0 / EPS_TAU) ** 2 * OPERA_FAKE_RATE / (Q * 5) / 25.0)

# existing underground laboratories, already on the compass, at the baselines the tau programme wants
LABS = [("Soudan (MINOS far hall)", 736.1, 335.9), ("Ash River (NOvA far hall)", 811.3, 335.2),
        ("SURF (DUNE far site)", 1289.1, 287.7)]
LAB_REACH = []
for nm, D, az in LABS:
    Q, _ = chirp_Q(D, "nu")
    LAB_REACH.append(dict(lab=nm, range_km=D, bearing_deg=az, tilt_mrad=D / (2 * cc.R_E) * 1e3,
                          up_going_end_km=cc.exits_smooth(D / (2 * cc.R_E), 35.0)[0], Q_per_t_yr=Q,
                          kt_for_5sigma_5yr=(5.0 / EPS_TAU) ** 2 * OPERA_FAKE_RATE / (Q * 5) / 1e3,
                          Z_40kt_5yr=significance(Q, 4e4, 5, OPERA_FAKE_RATE),
                          within_study_tilt_range=bool(D / (2 * cc.R_E) * 1e3 <= 67.0)))

# ---------------------------------------------------------------- figure
Q198 = CHIRP["UIUC, 198 km S"]["Q_per_t_yr"]
Q198N = CHIRP["655 km N (Lake Superior)"]["Q_per_t_yr"] / (655.0 / 198.4) ** 2
def mass_for_5sigma(L, Qref, f, years=5.0):
    Q = Qref * (np.asarray(L, float) / 198.4) ** 2
    return (5.0 / EPS_TAU) ** 2 * f / (Q * years) / 1e3            # kt

fig = plt.figure(figsize=(11.0, 8.0))
gs = fig.add_gridspec(2, 2, hspace=0.36, wspace=0.24, left=0.075, right=0.975, top=0.885, bottom=0.115)

# (a) how much detector, how far away
axa = fig.add_subplot(gs[0, 0])
LL = np.geomspace(120, 1500, 300)
for Qref, col, lab in ((Q198N, "#1f6f8b", "$\\nu_\\tau$, north beam ($\\mu^-$ decays)"),
                       (Q198, "#7a2e21", "$\\bar\\nu_\\tau$, south beam ($\\mu^+$ decays)")):
    axa.plot(LL, mass_for_5sigma(LL, Qref, OPERA_FAKE_RATE), color=col, lw=2.0, label=lab)
    axa.plot(LL, mass_for_5sigma(LL, Qref, OPERA_FAKE_RATE / 10), color=col, lw=1.1, ls=(0, (4, 2)))
axa.axhspan(0, 40, color="#0b6e4f", alpha=.07, lw=0)
axa.axhline(40, color="#0b6e4f", lw=0.9, ls=":")
axa.text(126, 44, "DUNE far-detector mass", fontsize=6.4, color="#0b6e4f", ha="left", va="bottom")
axa.axhline(10, color="0.5", lw=0.8, ls=":")
axa.text(126, 10.6, "10 kt", fontsize=6.4, color="0.45", ha="left", va="bottom")
MARK = [("UIUC", 198.4, "#7a2e21", (0, 10), "center"), ("W. Tennessee", 655.0, "#7a2e21", (-6, 11), "center"),
        ("NE Mississippi", 830.0, "#7a2e21", (16, 9), "left"), ("Soudan", 736.1, "#1f6f8b", (-13, -11), "right"),
        ("Ash River", 811.3, "#1f6f8b", (13, -5), "left"), ("SURF", 1289.1, "#1f6f8b", (0, 11), "center")]
for nm, L, col, off, ha in MARK:
    Qr = Q198 if col == "#7a2e21" else Q198N
    m = float(mass_for_5sigma(L, Qr, OPERA_FAKE_RATE))
    axa.plot([L], [m], "o", ms=5.2, color=col, mec="w", mew=0.8, zorder=6)
    axa.annotate("%s\n%.0f kt" % (nm, m), (L, m), xytext=off, textcoords="offset points",
                 fontsize=6.2, color=col, ha=ha, va="bottom" if off[1] > 0 else "top", zorder=7,
                 bbox=dict(boxstyle="round,pad=0.10", fc="w", ec="none", alpha=.8))
axa.set_xscale("log"); axa.set_yscale("log")
axa.set_xlim(120, 1500); axa.set_ylim(2, 4000)
axa.set_xlabel("baseline (km)", fontsize=8.2)
axa.set_ylabel("detector mass for a 5$\\sigma$ appearance signal in 5 yr (kt)", fontsize=7.8)
axa.set_title("(a)  the $\\tau$ programme wants a long baseline, not the 198 km hall", fontsize=8.6, loc="left", pad=5)
axa.legend(fontsize=6.6, frameon=False, loc="upper right")
axa.text(0.03, 0.03, "solid: OPERA's demonstrated rejection\n($10^{-4}$ fakes per interaction, 9 % $\\tau$ efficiency)\n"
                     "dashed: 10$\\times$ better",
         transform=axa.transAxes, fontsize=6.2, color="0.35", va="bottom", linespacing=1.35)
axa.tick_params(labelsize=7.2)

# (b) the tau flies far enough to see
axb = fig.add_subplot(gs[0, 1])
axb.axhspan(0.05, 1.0, color="#b3261e", alpha=.08, lw=0)
axb.text(3000, 0.14, "emulsion territory:\n$\\tau$ flight under a millimetre", fontsize=6.4, color="#b3261e", va="center", ha="right")
axb.axhspan(2.0, 300, color="#0b6e4f", alpha=.07, lw=0)
axb.text(3.4, 150, "resolvable by a silicon vertex tracker", fontsize=6.4, color="#0b6e4f", va="center")
Es = [STAGE_E[nm] for nm, *_ in tc.STAGES]; Ls_ = [TAU_FLIGHT[nm] * 1e3 for nm, *_ in tc.STAGES]
axb.plot(Es, Ls_, "-", color="#7a2e21", lw=1.4, zorder=4)
axb.plot(Es, Ls_, "o", ms=7, color="#7a2e21", mec="w", mew=1.0, zorder=5)
BOFF = {"RCS1": (8, -7, "left"), "RCS2": (8, -7, "left"), "RCS3": (8, -7, "left"),
        "RCS4": (-6, 9, "right"), "Collider store": (8, -9, "left")}
for (nm, *_), E, Lf in zip(tc.STAGES, Es, Ls_):
    dx, dy, ha = BOFF[nm]
    axb.annotate(nm, (E, Lf), xytext=(dx, dy), textcoords="offset points", fontsize=6.3,
                 color="#7a2e21", fontweight="bold", ha=ha)
for key, dy, dx, ha in (("OPERA", -13, -6, "right"), ("DONUT", 10, 0, "center"),
                        ("DUNE_cpv", -13, 0, "center"), ("DUNE_tau", 11, 6, "left")):
    d = LIT[key]; E = d["mean_E_nu_GeV"]; Lf = INELASTICITY * E / M_TAU * CTAU_TAU * 1e3
    axb.plot([E], [Lf], "s", ms=5.5, color="0.4", mec="w", mew=0.8, zorder=5)
    axb.annotate("%s\n%.2f mm" % (d["label"].split(" (")[0].split(",")[0], Lf), (E, Lf),
                 xytext=(dx, dy), textcoords="offset points", fontsize=6.2, color="0.35",
                 ha=ha, va="bottom" if dy > 0 else "top")
axb.set_xscale("log"); axb.set_yscale("log")
axb.set_xlim(3, 4000); axb.set_ylim(0.05, 300)
axb.set_xlabel("mean neutrino energy on axis (GeV)", fontsize=8.2)
axb.set_ylabel("$\\tau$ flight path, $\\gamma c\\tau$ (mm)", fontsize=7.8)
axb.set_title("(b)  at these energies the $\\tau$ decay is a tracking problem", fontsize=8.6, loc="left", pad=5)
axb.tick_params(labelsize=7.2)

# (c) the size of the sample
axc = fig.add_subplot(gs[1, 0])
BARS = [("OPERA, whole run (2008$-$12)", LIT["OPERA"]["interactions"], "0.55"),
        ("NuTeV, whole run (1996$-$97)", LIT["NuTeV"]["cc_events"], "0.55"),
        ("830 km hall, per tonne-year", SITE_ROWS["830 km S (NE Mississippi)"]["stages"]["FULL CHAIN"]["cc_per_t_yr"], "#7a2e21"),
        ("UIUC 198 km, per tonne-year", SITE_ROWS["UIUC, 198 km S"]["stages"]["FULL CHAIN"]["cc_per_t_yr"], "#7a2e21"),
        ("deep hall 9.25 km, per tonne-year", SITE_ROWS["deep hall, 9.25 km N"]["stages"]["FULL CHAIN"]["cc_per_t_yr"], "#b5541c"),
        ("on-site hall 2 km, per tonne-year", SITE_ROWS["on-site hall, 2 km N"]["stages"]["FULL CHAIN"]["cc_per_t_yr"], "#b5541c")]
yy = np.arange(len(BARS))[::-1]
for (lab, v, col), y in zip(BARS, yy):
    axc.barh(y, v, height=0.6, color=col, alpha=.85, lw=0)
    axc.text(v * 1.5, y, "%s" % ("%.1e" % v).replace("e+0", "$\\times10^{").replace("e+", "$\\times10^{") + "}$",
             fontsize=6.6, color=col if col != "0.55" else "0.4", va="center", ha="left", fontweight="bold")
axc.set_yticks(yy); axc.set_yticklabels([b[0] for b in BARS], fontsize=6.8)
axc.set_xscale("log"); axc.set_xlim(1e4, 1e13)
axc.set_xlabel("charged-current interactions", fontsize=8.2)
axc.set_title("(c)  two whole experiments, and one tonne-year of this beam",
              fontsize=8.6, loc="left", pad=5)
axc.tick_params(labelsize=7.2)
nutev_min = LIT["NuTeV"]["cc_events"] / (BARS[5][1] / 3.1536e7) / 60.0
axc.text(0.985, 0.70, "one tonne in the on-site hall collects\nNuTeV's entire dataset every %.0f minutes" % nutev_min,
         transform=axc.transAxes, fontsize=6.8, color="#b5541c", ha="right", va="center", fontweight="bold",
         bbox=dict(boxstyle="round,pad=0.25", fc="w", ec="#b5541c", lw=0.5, alpha=.92))

# (d) the sign
axd = fig.add_subplot(gs[1, 1])
NS = SITE_ROWS["655 km N (Lake Superior)"]["stages"]["FULL CHAIN"]
SS = SITE_ROWS["830 km S (NE Mississippi)"]["stages"]["FULL CHAIN"]
grp = [("north beam\n($\\mu^-$ decays)", NS["nutau_per_t_yr"], NS["nubartau_per_t_yr"]),
       ("south beam\n($\\mu^+$ decays)", SS["nutau_per_t_yr"], SS["nubartau_per_t_yr"])]
x = np.arange(2)
for i, (lab, nt, nb) in enumerate(grp):
    axd.bar(i - 0.19, nt, width=0.34, color="#1f6f8b", lw=0)
    axd.bar(i + 0.19, nb, width=0.34, color="#7a2e21", lw=0)
    axd.text(i - 0.19, nt + 0.006, "%.3f\n(%.0f %%)" % (nt, 100 * nt / (nt + nb)), fontsize=6.5,
             color="#1f6f8b", ha="center", va="bottom", fontweight="bold")
    axd.text(i + 0.19, nb + 0.006, "%.3f\n(%.0f %%)" % (nb, 100 * nb / (nt + nb)), fontsize=6.5,
             color="#7a2e21", ha="center", va="bottom", fontweight="bold")
axd.set_xticks(x); axd.set_xticklabels([g[0] for g in grp], fontsize=7.4)
axd.set_ylim(0, 0.33)
axd.set_ylabel("$\\tau$-flavour CC per tonne-year, on axis", fontsize=7.8)
axd.set_title("(d)  one straight, two sign-tagged beams, fired at once", fontsize=8.6, loc="left", pad=5)
axd.legend(handles=[Patch(color="#1f6f8b", label="$\\nu_\\tau$"), Patch(color="#7a2e21", label="$\\bar\\nu_\\tau$")],
           fontsize=7.2, frameon=False, loc="upper center", ncol=2)
axd.text(0.5, 0.70, "the world has identified %d $\\nu_\\tau$ interactions\n(DONUT 9, OPERA 10)\nand has never charge-tagged one"
         % (LIT["DONUT"]["nu_tau_identified"] + LIT["OPERA"]["nu_tau_identified"]),
         transform=axd.transAxes, fontsize=7.2, color="0.3", ha="center", va="center", fontweight="bold",
         bbox=dict(boxstyle="round,pad=0.3", fc="w", ec="0.7", lw=0.6))
axd.tick_params(labelsize=7.2)

fig.suptitle("What the corridor beam is worth", fontsize=12, y=0.965)
fig.text(0.5, 0.928, "$\\mu^+$ and $\\mu^-$ counter-rotate, so one straight fires a sign-tagged $\\bar\\nu_\\tau$ beam south and a "
                     "$\\nu_\\tau$ beam north, from a flux that is pure muon decay",
         ha="center", va="center", fontsize=8.4, color="0.35")
fig.text(0.5, 0.012,
         "Rates from tools/tau_compass.py and tools/nutau_sites.py (on-axis point detector, CSMS cross-sections, "
         "$\\tau$-threshold suppression). Reach uses the chirp as a weighted fit,\n"
         "$Z = \\epsilon\\,\\sqrt{(M T/f)\\sum_i S_i^2/B_i}$, with OPERA's own $\\tau$ efficiency and fake rate. "
         "Comparison numbers: DONUT PRD 78 052002, OPERA PRL 120 211801, NuTeV PRL 88 091802, DUNE TDR arXiv:2002.03005.",
         ha="center", va="center", fontsize=6.7, color="0.35")
for ext in ("pdf", "svg"):
    fig.savefig(os.path.join(ROOT, "static", "figs", "physics_case." + ext), bbox_inches="tight")

# ---------------------------------------------------------------- report and machine-readable output
def _clean(o):
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return float(o)
    return o

print("tau flight path (gamma c tau, at <E_tau> = %.1f <E_nu>):" % INELASTICITY)
for nm, *_ in tc.STAGES:
    print("  %-14s <E_nu> %7.1f GeV -> %7.2f mm" % (nm, STAGE_E[nm], TAU_FLIGHT[nm] * 1e3))
for k in ("DONUT", "OPERA", "DUNE_cpv", "DUNE_tau"):
    E = LIT[k]["mean_E_nu_GeV"]
    print("  %-14s <E_nu> %7.1f GeV -> %7.2f mm   [%s]" % (LIT[k]["label"].split(" (")[0], E,
          INELASTICITY * E / M_TAU * CTAU_TAU * 1e3, k))

print("\nper-site totals, on axis, per tonne-year (full chain):")
for lab, S in SITE_ROWS.items():
    f = S["stages"]["FULL CHAIN"]
    print("  %-28s L %6.1f km %-5s | nu_tau %.4f  nubar_tau %.4f  CC %.3e  S/B %.2e"
          % (lab, S["range_km"], S["sign_side"], f["nutau_per_t_yr"], f["nubartau_per_t_yr"],
             f["cc_per_t_yr"], f["S_over_B"]))

print("\nsignal-to-background bought by gating on the low-energy turns of the chirp:")
for lab in ("UIUC, 198 km S", "655 km S (W. Tennessee)", "830 km S (NE Mississippi)"):
    print("  %s" % lab)
    for r in GATES[lab]:
        if r["e_max_GeV"] in (100, 314, 5000):
            print("    turns below %5.0f GeV: %7.2f tau/kt-yr, %10.3e CC/kt-yr, S/B %8.2e, "
                  "needs a fake rate %5.1fx better than OPERA's"
                  % (r["e_max_GeV"], r["tau_per_t_yr"] * 1e3, r["cc_per_t_yr"] * 1e3, r["S_over_B"],
                     r["better_than_opera_by"]))

print("\nreach, using the chirp as a weighted fit (eps = %.2f, OPERA's fake rate %.1e per interaction):"
      % (EPS_TAU, OPERA_FAKE_RATE))
for lab, c in CHIRP.items():
    need = (5.0 / EPS_TAU) ** 2 * OPERA_FAKE_RATE / (c["Q_per_t_yr"] * 5) / 1e3
    print("  %-28s Q %.2e /t-yr | %8.1f kt for 5 sigma in 5 yr | Z(10 kt, 5 yr) %5.2f, "
          "%5.2f at 10x better rejection" % (lab, c["Q_per_t_yr"], need,
          c["Z_10kt_5yr_opera_rejection"], c["Z_10kt_5yr_10x_opera"]))
print("  existing underground laboratories, already on the compass:")
for r in LAB_REACH:
    print("    %-26s %7.1f km, bearing %5.1f deg, tilt %5.1f mrad%s | %5.1f kt for 5 sigma, Z(40 kt, 5 yr) %5.1f"
          % (r["lab"], r["range_km"], r["bearing_deg"], r["tilt_mrad"],
             "" if r["within_study_tilt_range"] else " (beyond this study's verified tilt range)",
             r["kt_for_5sigma_5yr"], r["Z_40kt_5yr"]))

cc_onsite = SITE_ROWS["on-site hall, 2 km N"]["stages"]["FULL CHAIN"]["cc_per_t_yr"]
print("\nthe non-tau dataset: %.2e CC interactions per tonne-year in the on-site hall (%.0f Hz per tonne)."
      % (cc_onsite, cc_onsite / 3.1536e7))
print("  NuTeV's entire 1996-97 sample (%.2e events) arrives every %.0f minutes per tonne."
      % (LIT["NuTeV"]["cc_events"], LIT["NuTeV"]["cc_events"] / (cc_onsite / 3.1536e7) / 60))

json.dump(_clean(dict(
    note="physics case: sample size, sign, tau flight path and appearance reach for the corridor beam",
    inputs=dict(tau_efficiency=EPS_TAU, opera_fake_rate_per_interaction=OPERA_FAKE_RATE,
                inelasticity_Etau_over_Enu=INELASTICITY, literature=LIT),
    stage_energies_GeV=STAGE_E, tau_flight_m=TAU_FLIGHT,
    sites=SITE_ROWS, chirp_gates=GATES, reach=CHIRP, underground_labs=LAB_REACH,
    on_site_cc_per_t_yr=cc_onsite,
    nutev_equivalent_minutes_per_tonne=LIT["NuTeV"]["cc_events"] / (cc_onsite / 3.1536e7) / 60),
), open(os.path.join(ROOT, "static", "geo", "physics_case.json"), "w"), indent=1)
print("\nwrote static/figs/physics_case.{svg,pdf}, static/geo/physics_case.json")
