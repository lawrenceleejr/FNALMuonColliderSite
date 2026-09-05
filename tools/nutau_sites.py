#!/usr/bin/env python3
"""nu_tau versus nubar_tau at four sites, per tonne-year, against DUNE.

mu+ and mu- counter-rotate, so each direction of a straight carries ONE sign:
the beam going one way is mu+ decays (nubar_mu + nu_e), the other way mu-
decays (nu_mu + nubar_e).  A site therefore sees mostly nubar_tau (from
nubar_mu -> nubar_tau, amplitude 0.95) with a nu_tau minority from
nu_e -> nu_tau (amplitude sin^2 2theta13 sin^2 theta23 ~ 0.050), or the
mirror image.  Which sign a site gets is set by the circulation sense; the
numbers here take mu+ southbound on the east straights (UIUC gets nubar) and
report both orientations.

Rates use the on-axis machinery of tau_compass.py (point detector on the
axis; CSMS CC cross-sections; tau-threshold suppression; north-aimed decays
from chain_timing.json).  Because the on-axis flux falls as 1/L^2 and the
appearance probability rises as L^2, the per-tonne rate is the same at every
far site; what changes with baseline is the signal-to-background (~L^2), the
spot size (~L), and the whole-plane nu_tau count (~L^2).

Sites: FNAL on site (2.0 km north, the beam still 4.5 m underground), UIUC
(198.4 km S), Green Bay (300 km N; the bay's water is off the UIUC line by
2.8 deg, see two_ends.py), Lake Superior (655 km N).  DUNE reference: TDR
vol. II, 130 nu_tau CC / yr in 40 kt in the CP-optimised neutrino mode,
before efficiency; antineutrino mode scaled by the TDR's nubar/nu ratio of
selected nu_tau-CC backgrounds (27/46 and 20/32, ~0.6).

Output: static/figs/nutau_sites.{svg,pdf}; static/geo/nutau_sites.json.
"""
import io, contextlib, json, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
with contextlib.redirect_stdout(io.StringIO()):
    import tau_compass as tc                      # reuses (and regenerates) the compass machinery
STYLE = "/root/.claude/skills/synced/tufte-plots/assets/tufte.mplstyle"
if os.path.exists(STYLE):
    plt.style.use(STYLE)
plt.rcParams["figure.dpi"] = 130

AMP_MU = tc.P_AMP                                 # nu_mu -> nu_tau: cos^4(th13) sin^2(2 th23) ~ 0.95
AMP_E = 0.0875 * 0.57                             # nu_e -> nu_tau: sin^2(2 th13) sin^2(th23) ~ 0.050 (NuFIT-class values)
def f_e_axis(x):                                  # on-axis nu_e (nubar_e) spectrum, E_nu = x E_mu
    return 12.0 * x ** 2 * (1.0 - x)
def g_e_plume(y):                                 # whole-plume nu_e spectrum, mean 0.30
    return 2.0 - 6.0 * y ** 2 + 4.0 * y ** 3
def g_e_point_collider(y):
    return np.where(y >= tc.Y_MIN_COLL, g_e_plume(y), 0.0)

SITES = [("FNAL site (2.0 km N)", 2.0), ("UIUC (198 km S)", 198.4), ("Green Bay (300 km N)", 300.0), ("Lake Superior (655 km N)", 655.0)]
Y = tc.Y

def rates(nm, L, side):
    """Per tonne-year on axis for stage nm at range L.  side='nubar': mu+ decays arrive (nubar_mu, nu_e);
    side='nu': mu- decays arrive (nu_mu, nubar_e).  Returns dict of signal and background components."""
    _, E1, E2, n_cycle, smeared = next(s for s in tc.STAGES if s[0] == nm)
    E, w = tc._egrid(E1, E2)
    Sm = tc.POINT[nm](Y); Sm = Sm / np.trapezoid(Sm, Y)
    Se = (g_e_point_collider if smeared else f_e_axis)(Y); Se = Se / np.trapezoid(Se, Y)
    Enu = Y[None, :] * E[:, None]
    N_sign = n_cycle * tc.CYCLES_YR
    dens = N_sign * (E / tc.M_MU) ** 2 / (math.pi * (L * 1e5) ** 2)
    if smeared:
        dens = dens / tc.DIV_SUPP
    P_mu = tc.posc(Enu, L)                                        # amplitude 0.95 already inside
    P_e = P_mu * (AMP_E / AMP_MU)
    xi = tc.tau_thresh(Enu)
    s_nu, s_nub = tc.sig_cc(tc._SIG_NU, Enu), tc.sig_cc(tc._SIG_NUB, Enu)
    if side == "nubar":       # nubar_mu -> nubar_tau (sigma_nubar); nu_e -> nu_tau (sigma_nu)
        sig_main, sig_e = s_nub, s_nu
    else:                     # nu_mu -> nu_tau; nubar_e -> nubar_tau
        sig_main, sig_e = s_nu, s_nub
    def integ(S, P, sig, thr=True):
        return float((w * dens * np.trapezoid(S[None, :] * P * sig * (xi if thr else 1.0), Y, axis=1)).sum() * tc.N_A_T)
    main = integ(Sm, P_mu, sig_main)                              # tau from the mu-type neutrino
    minor = integ(Se, P_e, sig_e)                                 # tau from the e-type neutrino
    bg_main = integ(Sm, np.ones_like(P_mu), sig_main, thr=False)  # unoscillated mu-type CC, the tau background
    bg_e = integ(Se, np.ones_like(P_mu), sig_e, thr=False)
    return dict(main=main, minor=minor, bg_main=bg_main, bg_e=bg_e)

OUT = dict(sites={}, notes=[])
for site, L in SITES:
    side = "nu" if "N)" in site else "nubar"                      # mu+ southbound on the east straights
    rec = dict(range_km=L, sign_side=side, stages={})
    tot = dict(nutau=0.0, nubartau=0.0, bg=0.0, store_nutau=0.0, store_nubartau=0.0, store_bg=0.0)
    for nm, *_ in tc.STAGES:
        r = rates(nm, L, side)
        if side == "nubar":
            nub, nut = r["main"], r["minor"]
        else:
            nut, nub = r["main"], r["minor"]
        bg = r["bg_main"] + r["bg_e"]
        rec["stages"][nm] = dict(nutau_cc_per_t_yr=nut, nubartau_cc_per_t_yr=nub, mu_type_plus_e_type_cc_per_t_yr=bg)
        tot["nutau"] += nut; tot["nubartau"] += nub; tot["bg"] += bg
        if nm == "Collider store":
            tot["store_nutau"], tot["store_nubartau"], tot["store_bg"] = nut, nub, bg
    rec["totals"] = tot
    rec["S_over_B_full_chain"] = (tot["nutau"] + tot["nubartau"]) / tot["bg"]
    rec["S_over_B_store"] = (tot["store_nutau"] + tot["store_nubartau"]) / tot["store_bg"]
    # spot radii on axis: collider smeared r50 = 1.2 sigma L; RCS pencil 1/gamma L at the stage's lowest/highest energy
    rec["spot_r50_collider_m"] = 1.2 * tc.SIG_THETA * L * 1e3
    rec["spot_rcs4_top_m"] = L * 1e3 * tc.M_MU / 5000.0
    rec["spot_rcs1_bottom_m"] = L * 1e3 * tc.M_MU / 63.0
    OUT["sites"][site] = rec
    print("%-26s %6.1f km  side=%-5s | full chain: nu_tau %.3f  nubar_tau %.3f  per t-yr (bg %.2e, S/B %.1e) | store: %.3f / %.3f (S/B %.1e) | spots: coll %.1f m, RCS4 %.1f m, RCS1 %.0f m" % (
        site, L, side, tot["nutau"], tot["nubartau"], tot["bg"], rec["S_over_B_full_chain"], tot["store_nutau"], tot["store_nubartau"],
        rec["S_over_B_store"], rec["spot_r50_collider_m"], rec["spot_rcs4_top_m"], rec["spot_rcs1_bottom_m"]))

# whole-plane nu_tau counts per year (every neutrino crossing the exit plane, E > 3.5 GeV), which DO grow as L^2
plane = {}
for site, L in SITES:
    tot_plane = 0.0
    for nm, E1, E2, n_cycle, smeared in tc.STAGES:
        Pp = float(tc.pbar(L, E1, E2, tc.g_plume, e_cut=3.5)[0])
        tot_plane += n_cycle * tc.CYCLES_YR * Pp                   # per sign, per year
    plane[site] = tot_plane
    OUT["sites"][site]["plane_crossing_tau_per_year_one_sign"] = tot_plane
    print("   whole plane, one sign: %.2e tau-flavour neutrinos/yr crossing the exit plane at %s" % (tot_plane, site))

# DUNE reference
DUNE = dict(source="DUNE TDR vol. II arXiv:2002.03005, sec. 4.1.1.3 (130 nu_tau CC/yr, CP-optimised, 40 kt, before efficiency); "
                   "antineutrino mode scaled by the TDR's selected nu_tau-CC background ratios 27/46 and 20/32",
            nu_mode_cc_per_yr=130.0, nubar_mode_cc_per_yr=130.0 * 0.6, mass_t=40000.0)
DUNE["nu_mode_per_t_yr"] = DUNE["nu_mode_cc_per_yr"] / DUNE["mass_t"]
DUNE["nubar_mode_per_t_yr"] = DUNE["nubar_mode_cc_per_yr"] / DUNE["mass_t"]
DUNE["nominal_3p5_plus_3p5_avg_per_t_yr"] = 0.5 * (DUNE["nu_mode_per_t_yr"] + DUNE["nubar_mode_per_t_yr"])
DUNE["S_over_B_approx"] = 0.05
OUT["dune"] = DUNE
OUT["notes"] = ["on-axis point detector; per-tonne rates are baseline-independent (1/L^2 flux x L^2 probability)",
                "sign side: mu+ southbound on the east straights -> UIUC sees nubar_mu, nu_e; northbound beams carry mu- decays. Reversing the circulation swaps every site",
                "nu_e -> nu_tau amplitude 0.050 (sin^2 2th13 sin^2 th23); vacuum formulas: matter effects cancel at these L << 4E/A",
                "DUNE nubar-mode yield is a scaling, not a TDR number"]
json.dump(OUT, open(os.path.join(ROOT, "static", "geo", "nutau_sites.json"), "w"), indent=1)

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(12.0, 5.4))
gs = fig.add_gridspec(1, 2, left=0.07, right=0.985, top=0.83, bottom=0.24, wspace=0.28, width_ratios=[1.35, 1])
ax = fig.add_subplot(gs[0]); ax2 = fig.add_subplot(gs[1])
C_NU, C_NUB, C_D = "#1c5a96", "#b5541c", "0.45"
group_labels, group_x, sub_labels, sub_x = [], [], [], []
x = 0
for site, L in SITES:
    rec = OUT["sites"][site]; t = rec["totals"]
    x0 = x
    for cfg, (nut, nub) in (("store", (t["store_nutau"], t["store_nubartau"])), ("chain", (t["nutau"], t["nubartau"]))):
        ax.bar(x, nut, width=0.8, color=C_NU, edgecolor="w", lw=0.5)
        ax.bar(x, nub, width=0.8, bottom=nut, color=C_NUB, edgecolor="w", lw=0.5)
        ax.text(x, nut + nub, "%.2f" % (nut + nub), ha="center", va="bottom", fontsize=6.4, color="0.2")
        sub_labels.append(cfg); sub_x.append(x); x += 1
    group_labels.append(site.replace(" (", "\n(")); group_x.append(0.5 * (x0 + x - 1))
    x += 0.7
for mode, val, frac_nu in (("DUNE\n$\\nu$ mode", DUNE["nu_mode_per_t_yr"], 0.9), ("DUNE\n$\\bar\\nu$ mode", DUNE["nubar_mode_per_t_yr"], 0.2)):
    ax.bar(x, val * frac_nu, width=0.8, color=C_NU, edgecolor="w", lw=0.5, hatch="//")
    ax.bar(x, val * (1 - frac_nu), width=0.8, bottom=val * frac_nu, color=C_NUB, edgecolor="w", lw=0.5, hatch="//")
    ax.text(x, val, "%.4f" % val, ha="center", va="bottom", fontsize=6.4, color="0.2")
    group_labels.append(mode); group_x.append(x); x += 1
ax.set_xticks(sub_x); ax.set_xticklabels(sub_labels, fontsize=6.2, color="0.35")
for gx, gl in zip(group_x, group_labels):
    ax.text(gx, -0.11, gl, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=6.6, color="0.15")
ax.set_yscale("log"); ax.set_ylim(3e-4, 2.0)
ax.set_ylabel(r"$\nu_\tau$ / $\bar\nu_\tau$ CC per tonne$\cdot$yr, on axis", fontsize=8.5)
ax.tick_params(labelsize=7)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=C_NU, label=r"$\nu_\tau$ CC"), Patch(color=C_NUB, label=r"$\bar\nu_\tau$ CC"),
                   Patch(facecolor="w", edgecolor="0.4", hatch="//", label="DUNE, sign split illustrative")],
          fontsize=7, frameon=False, loc="upper right", ncol=3)
ax.set_title("per tonne on the axis: the same at every far site; the SIGN is set by which way $\\mu^+$ circulates\n"
             "(here $\\mu^+$ southbound on the east straights: UIUC gets $\\bar\\nu_\\tau$, the northern sites $\\nu_\\tau$; the minority sign is $\\nu_e\\to\\nu_\\tau$)",
             fontsize=8.6, pad=10)
# right: S/B and spot size vs baseline
Lg = np.geomspace(1.5, 1500, 200)
sb_chain = OUT["sites"]["UIUC (198 km S)"]["S_over_B_full_chain"] * (Lg / 198.4) ** 2
sb_store = OUT["sites"]["UIUC (198 km S)"]["S_over_B_store"] * (Lg / 198.4) ** 2
ax2.plot(Lg, sb_chain, color="0.15", lw=1.8, label="full chain")
ax2.plot(Lg, sb_store, color="#b5541c", lw=1.4, label="collider store")
for site, L in SITES:
    rec = OUT["sites"][site]
    ax2.plot([L], [rec["S_over_B_full_chain"]], "o", ms=6, color="0.15", mec="w", zorder=5)
    OFF = {"FNAL site (2.0 km N)": ((8, -4), "left", "top"), "UIUC (198 km S)": ((-8, 4), "right", "bottom"),
           "Green Bay (300 km N)": ((8, -6), "left", "top"), "Lake Superior (655 km N)": ((-6, 8), "right", "bottom")}
    off, ha, va = OFF[site]
    r50 = rec["spot_r50_collider_m"]
    ax2.annotate(site.split(" (")[0] + ("\n$r_{50}$ %.0f m" % r50 if r50 >= 10 else "\n$r_{50}$ %.2g m" % r50), (L, rec["S_over_B_full_chain"]),
                 xytext=off, textcoords="offset points", fontsize=6.3, ha=ha, va=va, color="0.25")
ax2.axhline(DUNE["S_over_B_approx"], color="0.5", lw=1.0, ls="--")
ax2.text(1.7, DUNE["S_over_B_approx"] * 1.3, "DUNE: $\\sim$130 $\\nu_\\tau$ CC against a few thousand $\\nu_\\mu$ CC per year, S/B $\\sim$ 5%", fontsize=6.6, color="0.4")
ax2.set_xscale("log"); ax2.set_yscale("log")
ax2.set_xlim(1.5, 1500); ax2.set_ylim(1e-12, 1)
ax2.set_xlabel("baseline L (km)", fontsize=8.5)
ax2.set_ylabel(r"$\tau$ CC / ($\nu_\mu$-type + $\nu_e$-type CC), on axis", fontsize=8.5)
ax2.tick_params(labelsize=7)
ax2.legend(fontsize=7, frameon=False, loc="lower right")
ax2.set_title("what baseline buys: signal-to-background $\\propto L^2$ (rate does not change); $r_{50}$ = collider spot radius", fontsize=8.6, pad=10)
fig.suptitle(r"$\nu_\tau$ versus $\bar\nu_\tau$ at four sites, per tonne-year, against DUNE", fontsize=11, y=0.975)
fig.text(0.5, 0.03,
         "Rates for a point detector on the beam axis (tau_compass.py machinery: north-aimed decays per cycle, 5 Hz, $1.2\\times10^7$ s/yr, CSMS CC, $\\tau$ threshold). "
         "$\\nu_e\\to\\nu_\\tau$ amplitude 0.050 vs 0.95 for $\\nu_\\mu\\to\\nu_\\tau$.\nAt 2 km the collider spot is 0.36 m and the pencils' 4 cm to 3 m: a tonne cannot sit on the axis, "
         "and S/B is $10^{-11}$. DUNE: TDR 130 $\\nu_\\tau$ CC/yr in 40 kt ($\\nu$ mode, before efficiency); $\\bar\\nu$ mode scaled by 0.6 from the TDR's background ratios.",
         ha="center", va="center", fontsize=6.6, color="0.35")
for ext in ("pdf", "svg"):
    fig.savefig(os.path.join(ROOT, "static", "figs", "nutau_sites." + ext), bbox_inches="tight")
print("wrote nutau_sites")
