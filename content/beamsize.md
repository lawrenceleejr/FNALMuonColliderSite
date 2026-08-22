---
title: How Big Is the Beam? Spot Size from the Collider Optics
weight: 4
---

Everything else in this study treats the neutrino beam as a line. It is not
one, and the width matters in three separate places: it sets how small a
detector can be and still catch the beam, it sets how wide the
surface-grazing exit strips are, and it decides whether the machine optics
have any bearing on the siting argument at all. This page works the width
out from the collider's beam parameters — and it has been rewritten, because
the first simulation of this exact beam on a real lattice,
**MINT ([arXiv:2608.02718](https://arxiv.org/abs/2608.02718), Choi, Hostert,
Li & Liu)**, contradicts the assumption an earlier revision of this page was
built on. That revision claimed the far-field spot was set by the decay
kinematics alone, r = L/γ, with the muon beam's divergence confined to a
1.7 %-weight halo from a ±6 m final-focus drift. MINT, decaying muons along
the IMCC hybrid v0.6+v0.9 interaction-region lattice — the same
ε<sub>N</sub> = 25 µm·rad and β\* = 1.5 mm this study assumes — finds the
opposite:

> "the resulting neutrino flux spot size is dominated by the muon beam
> divergence rather than the decay kinematics, θ<sub>μ</sub> ∼ O(0.1 mrad)
> ≫ 1/γ ∼ 0.02 mrad, washing out neutrino energy-angle correlations (the
> so-called prism effect)."

**Headline:** on the current IMCC-class lattice the spot is set by the
optics, not the boost: **σ<sub>θ</sub> ≈ 0.15 mrad ≈ 7 × (1/γ)**, so the
50 %-containment radius at the on-site detector is **24 cm, not 2.7 cm**,
the on-axis flux density is **~100× lower** than the pencil formula gives,
and the energy–radius "prism" is gone. The 1/γ pencil is still achievable —
but only from a **dedicated dispersion-free drift** built beyond the final
focus and chicanes, which becomes this study's one explicit machine-design
ask (§4). The siting geometry is untouched (it is centreline, not width),
and the [safety case comes out stronger](../safety/), because the smearing
that costs detector rate is exactly the plume dilution the dose model was
hoping for.

## 1. The decay kernel is still universal — it is just no longer the answer

Muon decay is isotropic in the rest frame (per species, unpolarised beam).
Boosting by γ = 47 300 gives the beamed profile dN/dΩ = Nγ²/[π(1+γ²θ²)²]:
half the neutrinos inside 1/γ = 21 µrad, 99 % inside 9.95/γ. That kernel is
exact and this study keeps using it — for the **RCS plumes**, which come
from FODO straights with no final focus and really are 1/γ pencils (see
[the timing study](../timing/)), and as the kernel under the collider
plume's convolution.

What changed is the convolution. Each decaying muon points not along the
axis but along its own trajectory, and near the IP the trajectory spread is
enormous compared to 1/γ.

## 2. Where the divergence actually comes from

The physics an earlier revision of this page got right: **a field-free drift
conserves angles** — σ<sub>θ</sub> = √(ε<sub>N</sub>/γβ\*) everywhere
between the final-focus quadrupoles, and β(s) = β\* + s²/β\* describes the
envelope, not a shrinking angular spread. With β\* = 1.5 mm that is
σ<sub>θ</sub>\* = 0.59 mrad = 28 × (1/γ).

What it got wrong, in two places, both settled by MINT's Appendix B:

* **The final-focus drift is ~150 m long, not ±6 m.** "In the current
  design, the angular divergence for the longest drift section, ∼150 m, is
  about 0.1 − 0.2 mrad, which is the more appropriate scale for the angular
  divergence of the neutrino beam." All 150 m of drift radiates at that
  divergence — angles are conserved, so there is no "±L\* only" bookkeeping
  to hide behind.
* **The matching sections are not µrad-quiet.** The interaction region is
  dispersion-free, but "in the matching sections … the dispersive
  contribution can dominate the local angular spread" — larger, not smaller.

Two independent checks pin the effective number. MINT's on-axis flux at
5 km (2 × 10¹⁴ ν/cm²/yr from 6.1 × 10¹⁸ forward neutrinos) implies a
Gaussian σ<sub>θ</sub> = 0.14 mrad; their 60 %-containment radius of 1.3 m
at 5 km implies 0.15–0.20 mrad. The pure-pencil formula Φ = Nγ²/πL²
evaluated for their own straight overshoots their simulated peak by 87× —
which is just 2γ²σ<sub>θ</sub>² for that σ. **This study now takes
σ<sub>θ</sub> = 0.15 mrad, an on-axis density dilution of 2γ²σ<sub>θ</sub>²
≈ 101, and no energy–radius correlation inside the smeared core.** Beyond
the core, MINT finds a chicane plateau out to ~3 mrad and an arc plateau
beyond, together carrying O(1 %) of the flux at ~10⁻⁴ of the core density —
at 197 km, 3 mrad is a 590 m-radius faint halo this study previously did
not know it had.

The prism's loss matters twice over. The spot no longer selects energy by
radius, and the *spectrum* at any point inside the core is the
angle-integrated one: fluence-mean energies **0.35 E<sub>μ</sub> (νμ) and
0.30 E<sub>μ</sub> (ν̄e)** — half the on-axis means the pencil would deliver
— which softens per-fluence interaction rates by a further ~1.8× on top of
the density dilution.

## 3. The spot at each location

Baseline configuration (IP at 41.8443° N, 191 m ASL, 15.40 mrad),
E<sub>μ</sub> = 5 TeV, σ<sub>θ</sub> = 0.15 mrad. The "dedicated drift"
column is the §4 scenario — what a purpose-built pencil would restore.

| location | L from IP | r₅₀ (current lattice) | r₉₉ | r₅₀ (dedicated drift) |
|---|---|---|---|---|
| detector hall | 1.30 km | **24 cm** | 65 cm | 2.7 cm |
| north beam leaves the ground | 2.18 km | **40 cm** | 1.08 m | 4.6 cm |
| Fermilab fence | 2.84 km | **52 cm** | 1.41 m | 6.0 cm |
| Smith Road | 10.1 km | **1.8 m** | 5.0 m | 21 cm |
| under Urbana | 193 km | **35 m** | 96 m | 4.1 m |
| UIUC South Farms exit | 197.4 km | **36 m** | 98 m | 4.2 m |

(r₅₀ = 1.20 σ<sub>θ</sub>L and r₉₉ = 3.0 σ<sub>θ</sub>L in quadrature with
the kinematic 1/γ and 9.95/γ — a convolution MC reproduces these to ~5 %.
With the prism washed out, interaction containment ≈ flux containment; the
old "99 % of events inside 3L/γ" tightening is gone with the correlation
that produced it.)

On the ground, at the two grazing emergences:

* **The north emergence is still a lawn-scale object.** At 2.18 km and
  15.7 mrad grazing, the 50 %-flux streak is **0.8 m × 50 m** (99 %:
  2.2 m × 140 m) — a mowing strip, not a parking space, but still 655 m
  inside the fence and entirely on DOE land.
* **The southern exit strip is now a field, not a sidewalk.** At 197 km the
  50 % core is **72 m × 4.6 km ≈ 0.33 km²** (99 %: 196 m × 12.6 km). The
  land action at the UIUC South Farms is still single-owner and fenceable,
  but the [safety assessment](../safety/) §4 now describes a strip of
  research farmland measured in tens of hectares, not a sidewalk — with the
  compensation that the same smearing has already diluted the peak dose by
  the 101 the old model had to assume mitigations for.

## 4. Flux, rates, and the dedicated-drift ask

On-axis at the 1.3 km detector hall, per species, with the 0.854 store-decay
factor now included:

| quantity | current lattice (f = 0) | dedicated drift (f = 0.49) | old page (pencil, f = 1) |
|---|---|---|---|
| on-axis fluence (ν/cm²/yr) | **2.2 × 10¹⁵** | 1.1 × 10¹⁷ | 2.6 × 10¹⁷ |
| interactions per tonne-year | **2.9 × 10¹⁰** | 2.6 × 10¹² | 7 × 10¹² |
| in the ⌀16 cm × 2.4 m, 1 t tungsten cylinder | **2.4 × 10¹⁰/yr** (8 % intercept) | 2.7 × 10¹¹/yr | 7 × 10¹¹/yr (89 % intercept) |
| 99 % of flux within | 65 cm (57 t of tungsten) | 27 cm | 27 cm |

The one-tonne-cylinder detector concept does not survive the current
lattice: the beam is no longer smaller than the target, and a cylinder that
does contain it weighs 57 t. The better model is MINT's own benchmark —
trade target mass for reconstruction: a ~3 t, 39 m gaseous-argon TPC with a
vertex tracker sized to the beam centre, a dipole for TeV charge ID, and a
front muon monitor for the **rock-muon halo this page previously ignored**
(MINT: ~2 penetrating, highly polarised ⟨E⟩ ≈ 1.2 TeV muons per bunch
crossing through a 1.3 m-radius face at 5 km — scaled to our 1.3 km hall,
~5.6 µ/m² per crossing, arriving within tens of ps of the neutrinos, so
they are tagged geometrically, not by timing).

**The dedicated-drift ask, stated once and explicitly.** A dispersion-free
drift of half-length ℓ with a waist β\* ≈ ℓ has σ<sub>θ</sub> = √(ε/ℓ) =
1.2 µrad for ℓ = 350 m — a true 1/γ pencil at sub-mm beam size. The current
insertion cannot do this (it is a final focus; 0.59 mrad at the waist is
what β\* = 1.5 mm costs), and MINT says the fix out loud: "dedicated
straight sections for neutrino and muon beam dump physics may be required."
The corridor's physics numbers therefore carry a **pencil fraction f** — the
fraction of straight decays in such a drift: f = 0 is today's lattice,
f = 340/700 ≈ 0.49 is everything that is not the ±180 m of final focus and
chicanes. Every f > 0 sharpens the exit strips back toward the raw King
doses (the [safety assessment](../safety/) §4 puts numbers on that trade) —
physics rate and exit dose are now the same knob, and this study no longer
pretends otherwise.

## 5. ντ appearance along the line

The structural results survive untouched, because they are ratios: the
oscillation phase Δ₃₁ ∝ L/E is ≲10⁻⁴ everywhere on the line, so
P ∝ (L/E)² cancels the 1/L² flux and **the ντ event rate per tonne is the
same at every detector on the line**; the νμ CC background falls as 1/L², so
the far site keeps its **23 000× signal-to-background advantage** over the
on-site hall (S/B 3.8 × 10⁻⁸ at the UIUC exit vs 1.6 × 10⁻¹² at 1.3 km);
and L/E ≈ 0.06 km/GeV at the far site still parks the first oscillation
maximum at Δm² ≈ 20 eV², the free sterile-neutrino lever arm.

The absolute rate does not survive. The old 6.6 ντ CC per tonne-year was an
on-axis pencil number; with the density ÷101, the store factor, and the
softened spectrum (P σ ∝ 1/E, so the soft tail partly compensates), the
baseline is **≈ 0.1–0.2 ντ CC per tonne-year (f = 0), or ≈ 2.9 with the
dedicated drift** — a 100 t far detector on the South Farms collects ~15 or
~290 ντ CC/yr respectively (the old page said 660). The total oscillation
inventory is barely touched (it is a whole-plane count): **~3.4 × 10¹¹
oscillation-made ντ cross the South Farms exit plane per year**, ~7 × 10¹¹
to ground level counting the second IP's beam.

Two additions from MINT that this section previously did not know about:

1. **Rock-made ντ are a real background.** Primary neutrinos produce
   ντ + ν̄τ in the rock upstream of any detector (CC charm → D<sub>(s)</sub>
   → τντ dominant, plus inverse τ decay, resonant D\*<sub>s</sub>, ℓτ
   tridents): 2.0 × 10⁹/yr through a 1.3 m face at 5 km in MINT's geometry.
   Scaled along the corridor, they outnumber oscillation ντ by ~10² at the
   on-site hall — the near site's ντ are essentially all rock-made — and
   contribute at the ~10 % level at the UIUC exit, softer than the beam
   (⟨E⟩ ~ 300 GeV vs ~1.7 TeV) and separable on energy. These are scalings
   of MINT's benchmark, not simulations; a MINT run at 197 km would settle
   them.
2. **Wrong-sign contamination is bounded at O(10⁻⁹)** of the primary rate,
   which is what the "each exit plane is a single-sign beam" claim needed.

The honest statement stands, sharpened: the far site is the only place ντ
appearance can be attacked at all — now for background reasons as well as
rate — and it needs τ-ID progress *and* the dedicated drift to be more than
a curiosity.

## 6. Caveats

* σ<sub>θ</sub> = 0.15 mrad is a one-parameter Gaussian stand-in for MINT's
  simulated lattice (their Figs. 3, 5, 7); the chicane/arc plateaus are
  carried as a note, not modelled. The corridor's 700 m straight is a
  *proposal* — MINT simulated the IMCC ±180 m insertion — so the smeared
  fraction of a real corridor lattice is a lattice-design output, not an
  input.
* Unpolarised muons; polarisation reshapes each species' spectrum without
  moving the angular scales.
* Cross-sections: CSMS (arXiv:1106.3723) CC+NC per nucleon, log-log
  interpolated, evaluated at the fluence-mean energies; oscillations vacuum,
  two-channel, matter effects irrelevant at these phases.
* Numbers regenerate from
  [`tools/corridor_layout.py`](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/tools/corridor_layout.py)
  (σ<sub>θ</sub> and f are explicit inputs); the RCS pencils and the
  time-structure of all five beams are in [the timing study](../timing/).

*Live recomputation, with σ<sub>θ</sub>, ε<sub>N</sub>, β\* and the pencil
fraction f as inputs: [the beam-geometry tool](../tool/).*
