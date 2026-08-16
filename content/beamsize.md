---
title: How Big Is the Beam? Spot Size from the Collider Optics
weight: 4
---

Everything else in this study treats the neutrino beam as a line. It is not
quite one, and the width matters in three separate places: it sets how small
a detector can be and still catch the whole beam, it sets how wide the
surface-grazing exit strips are, and it decides whether β* and the beam
emittance have any bearing on the siting argument at all. This page works
the width out from the collider's beam parameters. Every number is
recomputed live in the [interactive tool](../tool/), where ε<sub>N</sub> and
β* are sliders.

**Headline:** the far-field spot is set by the decay kinematics alone —
**r = L/γ**, about 2.7 cm at the on-site detector and 4.2 m at the UIUC
exit. The muon beam's own divergence matters only inside the final-focus
drift (±L\* ≈ ±6 m of the 700 m straight): decays there feed a **wide,
faint halo carrying ~1.7 % of the flux at ~2 × 10⁻⁵ of the core's surface
density**, and no plausible β\* changes the core. A one-tonne tungsten
cylinder 16 cm across intercepts 99 % of the interactions that a
perfectly-aimed beam delivers to it — about **7 × 10¹¹ per year** at
1.3 km — and the grazing exit that the safety case turns on is an
**8.3 m × 534 m** strip.

## 1. The profile is a single universal curve

Muon decay is isotropic in the muon rest frame (for an unpolarised beam,
per neutrino species). Boosting to the lab by γ = E<sub>μ</sub>/m<sub>μ</sub>
= 47 300 at E<sub>μ</sub> = 5 TeV gives the standard beamed distribution

<div style="text-align:center;font-size:1.05rem;margin:.7rem 0">
dN/dΩ = N γ² / [ π (1 + γ²θ²)² ]
</div>

Its on-axis value, N γ²/π, is exactly the peak flux
Φ = Nγ²/πL² used everywhere else in this study — the flux formula and this
profile are the same statement, and that is a useful internal check.

Writing **x = γθ = γr/L**, the profile depends on nothing but x. The beam is
therefore **self-similar**: the same picture at every distance, scaled by
L/γ. Integrating:

| quantity | formula | 
|---|---|
| enclosed flux | F(x) = x²/(1+x²) |
| mean neutrino energy | ⟨E⟩(x) = ⟨E⟩₀ / (1+x²) |
| enclosed interactions (σν ∝ E) | F(x) = 1 − (1+x²)⁻² |

So **exactly half the neutrinos land inside 1/γ** — a clean statement of
what the "1/γ cone" means. The tails are heavy: 90 % needs 3/γ and 99 %
needs 9.95/γ.

The interaction profile is *tighter* than the flux profile, and this is a
real effect rather than a detail. The boost correlates energy with angle,
so off-axis neutrinos are softer; since σν grows with E, the interaction
density goes as (1+x²)⁻³ instead of (1+x²)⁻². The result:

| containment | flux | interactions |
|---|---|---|
| 50 % | 1.00 L/γ | 0.64 L/γ |
| 90 % | 3.00 L/γ | 1.47 L/γ |
| 99 % | 9.95 L/γ | 3.00 L/γ |

**75 % of all interactions happen inside 1/γ**, and 99 % inside 3/γ. A
detector only has to cover 3/γ, not 10/γ.

## 2. Where β\* comes in — and why it drops out

The muon beam has its own angular spread, and getting its reach right
requires respecting a fact an earlier revision of this page fumbled: **a
drift conserves angles**. In the field-free drift around the IP the beam's
angular divergence is

<div style="text-align:center;font-size:1.05rem;margin:.7rem 0">
σ<sub>θ</sub>\* = √( ε<sub>N</sub> / (γ β\*) )
</div>

*everywhere between the final-focus quadrupoles* — the familiar
β(s) = β\* + s²/β\* describes the growth of the beam *envelope*, not a
shrinking angular spread. For the IMCC 10 TeV parameters
(ε<sub>N</sub> = 25 µm·rad, β\* = 1.5 mm):

* σ<sub>θ</sub>\* = **0.59 mrad** — 28 times the 1/γ = 21 µrad decay cone —
  for every decay inside the ±L\* ≈ ±6 m final-focus drift;
* beam waist σ<sub>x</sub>\* = 0.89 µm, utterly negligible at kilometre range;
* beyond the quads the beam is recollimated: in the matching sections the
  divergence is µrad-scale, far below 1/γ.

The right picture is therefore **two populations**: a *collinear core* from
the (700 − 2L\*) ≈ 688 m of recollimated straight, and a *final-focus halo*
from the ±L\* drift, carrying weight 2L\*/L<sub>s</sub> ≈ **1.7 %** of the
flux, spread over the 0.59 mrad divergence. At the detector the halo is a
77 cm-radius wash — ~800× the core's area — so its surface density is
**~2 × 10⁻⁵ of the core's**: irrelevant to rates and siting, but ~140×
brighter than this page previously claimed. Two same-order effects belong
in the same bucket and await a real lattice: the **beam–beam deflection**
of the outgoing beam at the IP (maximum kick ~0.5 mrad at the design
disruption — comparable to σ<sub>θ</sub>\* itself) and dispersive angles in
the chromatic-correction sections. This halo, not any property of the
core, is the "IP-adjacent fan" that the effective-pencil-length parameter
of the [safety assessment](../safety/) stands in for.

**The far-field core is kinematic.** β\* sets the halo's angular width
(σ<sub>θ</sub>\* ∝ 1/√β\*) but its *weight* is fixed by geometry
(2L\*/L<sub>s</sub>); the tool's β\* slider spans 0.1–20 mm — a factor of
200 — and the core containment radii at the fence, Smith Road, and
Champaign County do not move.

## 3. The spot at each location

Baseline configuration (IP at 41.8443° N, 191 m ASL, 15.40 mrad),
E<sub>μ</sub> = 5 TeV:

| location | L from IP | L/γ | 50 % flux | 99 % flux | 99 % of events | on the ground |
|---|---|---|---|---|---|---|
| detector hall | 1.30 km | 2.7 cm | **2.7 cm** | 27 cm | 8.2 cm | — |
| north beam leaves the ground | 2.18 km | 4.6 cm | **4.6 cm** | 46 cm | 14 cm | 9.2 cm × 5.9 m |
| Fermilab fence | 2.84 km | 6.0 cm | **6.0 cm** | 60 cm | 18 cm | — |
| Smith Road | 10.1 km | 21 cm | **21 cm** | 2.1 m | 64 cm | — |
| under Urbana | 193 km | 4.1 m | **4.1 m** | 40 m | 12 m | — |
| UIUC South Farms exit | 197 km | 4.2 m | **4.2 m** | 42 m | 13 m | 8.3 m × 534 m |

Two consequences are worth pulling out.

**The north emergence is smaller than a parking space.** Where the beam
breaks the surface inside the Fermilab fence, at a grazing angle of
15.7 mrad, the 50 %-flux core paints a streak **9 cm wide and 5.9 m long**.
The "surface-grazing exit" that dominates muon-collider siting discussions
is, at this range, a patch of lawn.

**The southern exit strip is long and thin.** At 197 km and 15.6 mrad the
same core becomes **8.3 m × 534 m** — the 1/sin(graze) = ×64 elongation is
the entire reason exit strips are a land-use question at all. Taking the
99 % contour instead gives roughly 83 m × 5.3 km, which is the scale that
matches the footprint lengths quoted in the
[safety assessment](../safety/) (that model additionally spreads the plume
with ±0.5 mrad segmentation, which widens the strip and lowers the peak
dose in the same proportion).

## 4. What this means for the detector

At the detector, 99 % of *interactions* fall within 8.2 cm of the axis. A
cylinder of tungsten 16 cm in diameter and 2.4 m long weighs one tonne and
swallows essentially the whole beam. The study's headline rate —
~7 × 10¹¹ interactions per year in one tonne — therefore describes **an object you
could carry on a truck**, not a cavern; the hall is sized by the
instrumentation and the access, not by the beam.

This is the practical difference between a muon-collider neutrino beam and
a conventional one. A horn-focused beam is metres wide by construction, so
"per tonne" means a kilotonne detector. Here the beam is delivered
pre-collimated by the boost itself, and the same physics reach comes from a
target small enough to sit inside a spectrometer.

The trade to keep in mind is alignment, not size: with a 2.7 cm core at
1.3 km, pointing the beam is a 20 µrad problem. That is well inside what the
IMCC's ±1 mrad mover system resolves, but it does mean the detector's
position is a survey deliverable, and that the beam should be steerable onto
it rather than the other way round.

## 5. ντ appearance: the corridor as an oscillation experiment

The beam contains **no ντ at production** — a muon decay makes exactly one
μ-type and one e-type neutrino — so every ντ at a detector is
oscillation-made, driven by Δm²₃₁ = 2.5 × 10⁻³ eV² with the near-maximal
atmospheric amplitude (P ≈ 0.95 sin²Δ₃₁ for νμ→ντ; the e-type channel adds
~5 % through sin²2θ₁₃). At these energies and baselines the phase
Δ₃₁ = 1.267 Δm²L/E is at most ~10⁻⁴, so P ∝ (L/E)² — which produces the
punchline of this section:

**The ντ event rate per tonne is the same at every detector on the line.**
The flux falls as 1/L² and the appearance probability grows as L², and they
cancel *exactly* in the small-phase regime. On axis, with the CSMS
cross-sections, the Michel spectrum, and the τ-mass threshold factor:

| location | L | ντ CC per tonne-yr | νμ CC per tonne-yr (background) | S/B | ντ crossing this plane per year |
|---|---|---|---|---|---|
| on-site detector | 1.3 km | **6.6** | 4.1 × 10¹² | 1.6 × 10⁻¹² | 1.7 × 10⁷ |
| north exit | 2.2 km | **6.6** | 1.4 × 10¹² | 4.6 × 10⁻¹² | 4.8 × 10⁷ |
| Smith Road | 10.1 km | **6.6** | 6.7 × 10¹⁰ | 1.0 × 10⁻¹⁰ | 1.0 × 10⁹ |
| UIUC South Farms exit | 197.4 km | **6.6** | 1.8 × 10⁸ | 3.8 × 10⁻⁸ | **4.0 × 10¹¹** |

**Production over a year:** by the time the south beam surfaces at the
South Farms, **~4 × 10¹¹ of its 6.2 × 10¹⁸ forward-beamed neutrinos
(6 × 10⁻⁸) have become ντ**; with the second IP's south beam to Willard,
the machine delivers **~8 × 10¹¹ oscillation-made ντ per year** to ground
level. (Each straight's two exits carry opposite muon signs — μ⁻ decays aim
one way, μ⁺ the other — so each exit plane is a single-sign beam: ντ south
of IP1, ν̄τ from the CP-mirror channel likewise.) The near detector, 150×
closer, intercepts a beam in which only 1.7 × 10⁷ ντ/yr yet exist.

Three physics notes worth the ink:

1. **What "same rate, 23 000× better S/B" buys.** Six-ish ντ CC per
   tonne-year is DONUT's and OPERA's entire careers in a tonne — a
   100-tonne far detector on the South Farms would make **~660 ντ CC/yr**,
   two orders of magnitude beyond the world sample. But it sits under
   1.8 × 10⁸ νμ CC/tonne-yr: even emulsion-grade τ identification
   (~10⁻⁵–10⁻⁶ mis-ID per CC, charm-dominated) leaves the signal a factor
   ~10–100 under the fakes. The honest statement is that the far site is
   where ντ appearance *could* be attacked — the near site cannot at all —
   and that it needs τ-ID progress, kinematic rejection, or sign-selected
   charm vetoes beyond demonstrated performance.
2. **The oscillated flux is not a pencil.** Off axis the energy drops as
   (1+γ²θ²)⁻¹, so P grows as (1+γ²θ²)² — exactly cancelling the beamed
   flux profile: the ντ surface density is nearly **flat in solid angle**
   out to several 1/γ rather than peaked. The per-tonne numbers above are
   on-axis σ-weighted rates and barely change with detector radius.
3. **A free sterile-neutrino lever arm.** L/E ≈ 0.06 km/GeV at the far
   site puts the *first oscillation maximum* at Δm² ≈ 20 eV² — the corridor
   is, incidentally, an eV²–10 eV²-scale sterile-search geometry with a
   known-flavour TeV beam, something no accelerator facility currently
   offers.

## 6. Caveats

* Unpolarised muons are assumed. Polarisation changes the rest-frame
  angular distribution of each species (and hence the flavour composition
  off-axis) without moving the 1/γ scale.
* The two flavours in each beam have different *energy* spectra (⟨E⟩ ≈
  0.7 E<sub>μ</sub> for ν<sub>μ</sub>, 0.6 for ν<sub>e</sub>), so their
  interaction profiles differ slightly through σν ∝ E.
* Decays are taken as uniform along the straight and the beam as parallel
  outside the final focus. A real lattice has β varying through the
  insertion, which redistributes the small halo but not the core.
* Cross-sections: CSMS (arXiv:1106.3723) CC+NC per nucleon, log-log
  interpolated, as elsewhere in this study; ντ CC includes a τ-mass
  threshold factor (negligible above ~1 TeV). Oscillations: vacuum,
  two-channel (νμ→ντ dominant, νe→ντ at the 5 % level); matter effects are
  irrelevant at Δ₃₁ ≲ 10⁻⁴.

*Live recomputation, with ε<sub>N</sub> and β\* as inputs and
publication-quality figure export: [the beam-geometry tool](../tool/).*
