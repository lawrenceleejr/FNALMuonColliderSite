---
title: "The Case for UIUC: One Detector, Carried Along the Corridor"
weight: 13
---

[The physics case](../case/) optimised one figure of merit — τ-appearance
significance per tonne — found that it scales as L², and sent the τ detector
700 km out. That is true, and this page does not unsay it. But it is one figure
of merit, and a detector is not a tonne: it is a fixed object with a frontal
area, a length and a mass. Ask instead **what one real detector sees at each
hall on the corridor**, and the corridor's own geometry answers differently.
Tool: [`uiuc_case.py`](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/tools/uiuc_case.py);
machine-readable [uiuc_case.json](../geo/uiuc_case.json).

<div style="margin:1.2rem 0"><img src="../figs/uiuc_case.svg" alt="One 11.7 kt detector carried along the corridor: charged-current and tau yields against range, what it collects hall by hall, and tau appearance significance against fake rate" style="max-width:100%"></div>

## 1. Closer in, the beam is a pencil

Take one detector — 10 m × 10 m frontal area, 15 m of iron along the beam,
**11.7 kt** — and carry it to every hall. Each stage of the chirp illuminates
a spot of radius L/γ (the 1/γ pencils) or 1.18 σ<sub>θ</sub>L (the
divergence-smeared store), and the detector only counts the mass that spot
covers:

| stage | spot r₅₀ at 2 km | at 9.25 km | **at UIUC, 198 km** | at 830 km |
|---|---|---|---|---|
| RCS1 (110 GeV) | 1.35 m | 6.3 m | 134 m | 560 m |
| RCS2 (350 GeV) | 0.42 m | 2.0 m | 42 m | 175 m |
| RCS3 (760 GeV) | 0.20 m | 0.9 m | 19 m | 81 m |
| RCS4 (2.0 TeV) | **7 cm** | **34 cm** | **7.2 m** | 30 m |
| collider store | 0.35 m | 1.6 m | 35 m | 147 m |

At the on-site hall the top-energy pencil is seven centimetres across; at the
deep reference hall, a third of a metre. A kilotonne detector there is
**dark** — for RCS4 it lights 0 % of its mass at either hall, and the whole
complex delivers the same 5 × 10¹¹ interactions a year to 11.7 kt as it would
to a few tonnes, because the count is set by the plume, not the detector. The
near halls are where a compact, dense target belongs.

**UIUC is the first hall on the corridor where every stage lights the whole
detector**: the RCS4 spot reaches 10 m at 158 km. And it is the last where the
rate is still large. The same detector collects:

| hall | range | CC interactions / yr | rate | τ-flavour CC / yr | RCS4 mass lit |
|---|---|---|---|---|---|
| on-site hall | 2 km | 5.1 × 10¹¹ | 16 kHz | 17 | 0 % |
| deep hall | 9.25 km | 5.1 × 10¹¹ | 16 kHz | 330 | 0 % |
| **UIUC** | **198 km** | **2.6 × 10¹⁰** | **813 Hz** | **1 763** | **100 %** |
| far hall | 830 km | 1.5 × 10⁹ | 46 Hz | 1 763 | 100 % |

Read the last two rows together. The τ yield is identical — that is the
1/L² × L² cancellation, now for a real detector — and everything else is
**17.5× larger at UIUC**. Beyond 158 km the τ curve is flat and the rest of
the physics falls as 1/L². Every measurement that is limited by rate rather
than by signal-to-background is therefore an order of magnitude better at UIUC
than anywhere farther, and unavailable to a large detector anywhere nearer.
813 Hz in 11.7 kt is 0.07 Hz per tonne: no pile-up, every event isolated, a
triggerless readout if you want one. At the deep hall the same detector would
run at 16 kHz and see almost nothing of the high-energy beam.

## 2. What is rate-limited, and therefore UIUC's

**High-energy ν<sub>e</sub> physics, from zero.** Half of every muon decay is
an electron-flavour neutrino, so the south beam is 50 % ν<sub>e</sub> by
number, and at UIUC that is **1.5 × 10¹⁰ ν<sub>e</sub> charged-current
interactions a year** in 11.7 kt. Every high-energy neutrino beam ever built
was a ν<sub>μ</sub> beam with a 1–2 % ν<sub>e</sub> contamination from kaon
decay; the world's sample of ν<sub>e</sub> CC interactions above a few tens of
GeV is of order 10⁴, every one of them carrying a 5–10 % flux uncertainty
inherited from the kaons. The ν<sub>e</sub> cross-sections, structure functions
and the ν<sub>e</sub>/ν<sub>μ</sub> universality test at Q² up to 10⁴ GeV² are
essentially unmeasured, and this hall measures them with 10¹⁰ events from a
spectrum that is analytic.

**A leptonic sin²θ<sub>W</sub>, at a new Q².** Neutrino–electron elastic
scattering is the clean electroweak probe — no nucleon, no parton
distribution, no nuclear correction; NuTeV's 3σ anomaly lives in exactly the
hadronic corrections it avoids. It is also rare, about 5 × 10⁻⁴ of the
charged-current rate, which is why CHARM II's entire programme
([PLB 335, 246](https://doi.org/10.1016/0370-2693(94)91421-8)) was 5 429 events
for sin²θ<sub>W</sub> = 0.2324 ± 0.0083. At UIUC the 11.7 kt detector collects
**1.4 × 10⁷ ν–e events a year — CHARM II's whole run every 3.4 hours** — 86 %
of them ν<sub>e</sub>e, the channel with the charged-current interference.
The statistical floor is 9 × 10⁻⁵ per year, 4 × 10⁻⁵ in five, a factor of
four below LEP and SLD's combined 1.6 × 10⁻⁴, and at Q² = 2m<sub>e</sub>E<sub>e</sub>
≲ 0.1–2 GeV² — two orders of magnitude below NuTeV's, four below the Z pole,
where no precision leptonic measurement exists.

That statistical floor is not the measurement; it is the challenge to the
flux claim. Every 0.1 % of flux normalisation is 3 × 10⁻⁴ in
sin²θ<sub>W</sub>, so beating CHARM II by a hundred needs the stored-muon
count and polarisation to 10⁻³ — which is the standard the muon-decay beam is
supposed to meet, and this is the measurement that would prove it. There is
also a route that needs no normalisation at all: the ν̄<sub>μ</sub> and
ν<sub>e</sub> in the same beam come from the same decays, so their flux
*ratio* is fixed by muon-decay kinematics alone, and the two ν–e channels'
relative rate and spectral shape carry sin²θ<sub>W</sub> without the muon count
entering. That analysis has to separate two indistinguishable single-electron
final states by spectrum, and it is owed; the point is that the beam offers it
and no other does.

**And the rest.** Charm at ~5 % of CC is 10⁹ events a year; the dimuon
sample that anchors the strange-quark distribution today is ~10⁴. Structure
functions to Q² ~ 10⁴ GeV² with a known flux. Neutrino tridents, heavy neutral
leptons, dark-sector searches in a beam whose energy is known turn by turn.
All of it scales with rate, all of it is 17× UIUC's over the far hall.

**The ν̄<sub>τ</sub> count.** 1 763 τ-flavour interactions a year in 11.7 kt,
91 % of them ν̄<sub>τ</sub>; **6 000 a year in 40 kt**. Every ν<sub>τ</sub>
interaction ever identified, every thirty hours. The count was never UIUC's
problem.

## 3. The τ programme at UIUC is a number, and the number is 10⁻⁶

Purity is. The same L² that makes UIUC's ν–e sample seventeen times larger
makes its τ sample seventeen times dirtier than the far hall's, and the
question is what background rejection buys it back. Using the chirp as a
weighted fit, with OPERA's own 9 % τ efficiency:

| configuration | fake τ per interaction for 5σ in 5 yr | vs OPERA's 10⁻⁴ | at 10⁻⁵ | at 10⁻⁶ |
|---|---|---|---|---|
| far hall, 830 km, 10 kt | 3.2 × 10⁻⁵ | 3× better | 8.8σ | 28σ |
| **UIUC, 40 kt** | **7.2 × 10⁻⁶** | **14× better** | **4.2σ** | **13σ** |
| UIUC, 11.7 kt | 2.1 × 10⁻⁶ | 49× better | 2.3σ | 7.3σ |

So the requirement is stated: **a 40 kt magnetised detector at UIUC needs a
fake rate of 10⁻⁶ per interaction for a decisive result, 7 × 10⁻⁶ for
discovery.** Against OPERA's 10⁻⁴ that is one to two orders of magnitude, and
three things make it the right target rather than a wish. OPERA's τ flew half
a millimetre at 17 GeV; here it flies 3 mm at the lowest turns and 60 mm at
the highest, so the decay vertex is a tracking measurement, not an emulsion
one. OPERA's backgrounds were charm with a missed primary lepton, hadronic
re-interactions faking a kink, and large-angle muon scatters; at 100 GeV to
2 TeV the primary lepton is hard and hard to miss, re-interaction kinks are
resolved over centimetres, and Coulomb scattering falls as 1/p². And the beam
energy is known at every instant, so the fit knows which turns to trust. The
efficiency is upside, not a given: Z scales linearly with ε, and 9 % was set
by emulsion scanning, not by physics.

Be exact about what 10⁻⁶ delivers with 40 kt over five years. The
**background-free core is the turns below 100 GeV: about 15 identified
ν̄<sub>τ</sub> against 6 fakes** — a DONUT-sized sample, every event
charge-tagged, the first ν̄<sub>τ</sub> ever seen. The turns below 314 GeV
give ~130 identified against ~250 fakes, a measurement with S/B near one. The
full chirp gives ~2 700 identified τ candidates under 4 × 10⁵ fakes: a
cross-section sample that must be extracted by statistical subtraction with
the fake rate measured in the ν<sub>e</sub>-enriched control channels. The 13σ
is the weighted sum of all of that, and it is the low-energy turns that carry
it.

That is a harder detector than the far hall needs. It is also the same
detector — the far hall wants 10⁻⁵ to be comfortable — and building it once
makes every site on the corridor work. The reward for building the hard
version is that it sits in the beam with seventeen times the rest of the
physics.

## 4. Two halls, one beam

None of this is either/or, and the strongest configuration uses UIUC as the
anchor of a pair. The on-site hall at 2 km and UIUC at 198 km are a factor of
a hundred apart in L, they see the same stored muons in the same cycle, and
the Standard Model predicts their ratios exactly: charged-current rates in the
ratio (198/2)² of pure geometry — a flux-monitor cross-check at the 10⁻³
level — and τ appearance in the ratio 10⁴, i.e. **zero** at the near hall.
Any τ-flavour excess in the near hall above the L² law is the non-unitarity
floor of §6 of [the case](../case/), and the near hall has 10⁹ interactions a
year per tonne to bound it. UIUC is the far arm of that pair; nothing farther
out can pair with the site as cleanly, because the RCS pencils there have
outgrown any detector.

## 5. Why here

The corridor points at UIUC because UIUC is the partner: the exit strip is
university land, the emergence is the South Farms, and the whole siting
concept is *aim at a campus*. A detector there is not an addition to the
corridor, it is the reason the corridor has a direction. It would be the
largest particle-physics instrument ever hosted on a university campus, with
a 198 km baseline that is a two-hour drive from the source rather than a
flight — staffed, upgradeable, taught in. Soudan and Ash River are superb
halls and the far-site case for the small ring stands; but they are mines at
the end of a road, and the physics that only a big detector can do is not
there.

## 6. What is owed

* **The background simulation**, as before — it decides the τ programme at
  every site, and at UIUC it must reach 10⁻⁶. Charm with a missed primary
  lepton is the term to model first.
* **The ν–e analysis**: forward-electron identification at TeV energies (the
  E<sub>e</sub>θ² < 2m<sub>e</sub> cut is microradians here), the
  ν<sub>e</sub> quasi-elastic background, and the flux normalisation to 10⁻³
  that turns a 4 × 10⁻⁵ statistical floor into a result.
* **A detector concept** for 40 kt at 813–2 800 Hz: magnetised, few-centimetre
  vertex resolution, iron-scintillator or LAr with a silicon vertex layer — the
  design space is open.
* **A re-solve of the aim** if the collider and RCS3/4 stay on UIUC while
  RCS1/2 goes to a far hall (the split the far-site case suggests): UIUC keeps
  96 % of its rate but loses the low-energy turns that carry the τ core, so
  that split trades UIUC's ν̄<sub>τ</sub> discovery for the far hall's. The
  choice is a real one and this page's numbers are for the whole chain at
  UIUC.
