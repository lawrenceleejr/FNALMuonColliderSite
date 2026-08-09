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
exit. β\* and the emittance control the muon beam's own divergence, which
beats the decay cone over only **±4 cm of the 700 m straight**; changing β\*
by a factor of 200 barely moves the spot. A one-tonne tungsten cylinder
16 cm across intercepts 99 % of the interactions at the detector, and the
grazing exit that the safety case turns on is an **8.3 m × 534 m** strip.

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

The muon beam has its own angular spread. With normalised emittance
ε<sub>N</sub> and betatron function β(s) = β\* + s²/β\* along the straight,
the RMS divergence is

<div style="text-align:center;font-size:1.05rem;margin:.7rem 0">
σ<sub>θ</sub>(s) = √( ε<sub>N</sub> / (γ β(s)) )
</div>

which is *largest at the interaction point*, where β is smallest. For the
IMCC 10 TeV parameters (ε<sub>N</sub> = 25 µm·rad, β\* = 1.5 mm):

* σ<sub>θ</sub>\* = **0.59 mrad** — 28 times the 1/γ = 21 µrad decay cone;
* beam waist σ<sub>x</sub>\* = 0.89 µm, utterly negligible at kilometre range;
* but β grows as s²/β\*, so σ<sub>θ</sub> falls below 1/γ within
  **s = ±4.2 cm of the IP** — 1.2 × 10⁻⁴ of the 700 m insertion.

Past the first final-focus quadrupole (s ≈ 6 m, β ≈ 24 km) the divergence is
0.15 µrad, less than 1 % of the decay cone. **The far-field spot is
kinematic.** The tool's β\* slider spans 0.1–20 mm — a factor of 200 — and
moves the divergence-dominated region only from ±1 cm to ±15 cm. No
plausible optics choice changes the beam size at the fence, at Smith Road,
or in Champaign County.

What the divergence *does* make is a faint wide halo from the final-focus
region: at the detector it is a 77 cm-radius wash carrying 1.2 × 10⁻⁴ of the
flux spread over ~800× the core area, so about **1.5 × 10⁻⁷ of the core's
surface density**. That is the "IP-adjacent fan" that the effective-pencil-length
parameter of the [safety assessment](../safety/) stands in for, and it is
the one place where a real lattice — not this analytic model — is needed.

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
~10¹⁰ interactions per tonne-year — therefore describes **an object you
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

## 5. Caveats

* Unpolarised muons are assumed. Polarisation changes the rest-frame
  angular distribution of each species (and hence the flavour composition
  off-axis) without moving the 1/γ scale.
* The two flavours in each beam have different *energy* spectra (⟨E⟩ ≈
  0.7 E<sub>μ</sub> for ν<sub>μ</sub>, 0.6 for ν<sub>e</sub>), so their
  interaction profiles differ slightly through σν ∝ E; the ~0.65 blend used
  here is adequate for sizing.
* Decays are taken as uniform along the straight and the beam as parallel
  outside the final focus. A real lattice has β varying through the
  insertion, which redistributes the small halo but not the core.
* σν ≈ 0.35 × 10⁻³⁸ cm²/GeV per nucleon, flattened above ~1 TeV, as
  elsewhere in this study.

*Live recomputation, with ε<sub>N</sub> and β\* as inputs and
publication-quality figure export: [the beam-geometry tool](../tool/).*
