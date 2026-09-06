---
title: "The Case: Why This Is Two Discoveries for the Price of One Machine"
weight: 12
---

Every muon collider design treats its decay-neutrino flux as a liability. It is
the thing that drives the tunnel deep, forces the ring to wobble, and gets
written into siting studies as a constraint to be survived. This study starts
from the opposite premise: **aim it on purpose.** What comes out is not a
mitigated hazard but the most intense, best-characterised, and only
sign-tagged neutrino beam anyone has proposed — and it costs the energy-frontier
programme almost nothing, because the flux exists whether you use it or not.

The numbers below come from the study's own machinery
([`physics_case.py`](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/tools/physics_case.py),
building on [`tau_compass.py`](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/tools/tau_compass.py)
and [`nutau_sites.py`](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/tools/nutau_sites.py));
machine-readable [physics_case.json](../geo/physics_case.json). Comparison
numbers are from the published literature and are cited where used.

<div style="margin:1.2rem 0"><img src="../figs/physics_case.svg" alt="Four panels: detector mass needed for a 5 sigma tau appearance signal versus baseline, tau flight path versus beam energy, total interaction counts against past experiments, and the sign composition of the north and south beams" style="max-width:100%"></div>

## 1. The energy frontier gets there anyway

The argument for a 10 TeV muon collider does not need the neutrinos. Leptons
deliver their full energy to the hard collision, so the standard comparison in
the muon collider literature puts a 10 TeV μ⁺μ⁻ machine's reach for heavy pair
production alongside a 100 TeV proton collider's. Muons radiate
(m<sub>μ</sub>/m<sub>e</sub>)⁴ ≈ 1.8 × 10⁹ times less than electrons, so unlike
a circular e⁺e⁻ machine the ring can be small and still reach the energy
frontier: this study fits an 11.0 km collider and a 14.7 km final RCS **inside
the existing Fermilab boundary**, with real terrain and the real fence
([the ring stack](../rings/)).

That is the whole capital argument. Everything below is what you get for
building the halls and the detectors — a marginal cost against a machine whose
justification is already the energy frontier.

## 2. The best-understood neutrino beam ever built

This is the least contestable claim, and it is the one that makes the rest
work.

A conventional neutrino beam starts with protons on a target. The flux is
whatever the hadroproduction gives you, and knowing it means measuring π and K
yields off a replica target and propagating them through a horn — 5–10 % before
a near-detector constraint, a few percent after, and it is the dominant
systematic in every long-baseline programme. The intrinsic ν<sub>e</sub>
contamination, about 1 %, is both the principal background for ν<sub>μ</sub> →
ν<sub>e</sub> appearance and one of the harder things to pin down.

A muon storage ring has no hadrons in it. The flux is **muon decay** — a pure
QED three-body process whose spectrum and radiative corrections are known at
the 10⁻⁴ level:

<div style="text-align:center;font-family:monospace;margin:.8rem 0">
μ⁺ → e⁺ ν<sub>e</sub> ν̄<sub>μ</sub> &nbsp;&nbsp;&nbsp; μ⁻ → e⁻ ν̄<sub>e</sub> ν<sub>μ</sub>
</div>

so the beam is **exactly 50 % ν̄<sub>μ</sub> and 50 % ν<sub>e</sub> by number**
(or the mirror), with analytic spectra. There is no wrong-sign contamination,
no K-decay component, no target replica campaign. The ν<sub>e</sub> is not a
contaminant to be subtracted — it is half the beam, with a spectrum you can
write down.

What replaces the hadroproduction uncertainty is beam instrumentation: the
number of stored muons, the polarisation (which tilts the Michel spectrum and
is measurable from the decay-electron precession), and the divergence. Those
are accelerator measurements, and they are the kind that get to sub-percent.
One nuance this study adds: the collider straight's plume is
divergence-dominated — σ<sub>θ</sub> = 7/γ, per
[MINT](https://arxiv.org/abs/2608.02718) — so its on-axis flux depends on a
measured divergence, while the **RCS pencils are 1/γ-sharp and need no such
correction**. As §5 shows, the τ programme rides on the RCS turns, so it uses
the cleaner half of the beam.

## 3. The largest neutrino dataset ever taken, by three orders of magnitude

The per-tonne on-axis rate is **baseline-independent** — the flux falls as
1/L², the appearance probability rises as L², and they cancel — so the sample
size is set by detector mass, not by where you put it. What does change with
distance is the raw interaction rate, and near the ring it is extraordinary:

| hall | range | CC interactions per tonne-year, on axis |
|---|---|---|
| on-site hall | 2.0 km | **2.23 × 10¹⁰** (706 Hz per tonne) |
| deep reference hall | 9.25 km | 1.04 × 10⁹ |
| UIUC far hall | 198 km | 2.19 × 10⁶ |
| far-south hall | 830 km | 1.25 × 10⁵ |

For scale: **NuTeV's entire 1996–97 dataset — 1.62 × 10⁶ charged-current
events, still the reference measurement for sin²θ<sub>W</sub> in neutrino
scattering — arrives in the on-site hall every 38 minutes, per tonne.** OPERA's
whole run located 19,505 interactions; a tonne on site collects that in 28
seconds.

That is not a rhetorical flourish, it is a different kind of experiment. A
few-tonne detector in the near hall is a deep-inelastic-scattering facility with
10⁹ events a year, in a beam of known composition, with ν and ν̄ available
simultaneously (§4). Charm production is ~5 % of that, so the charm sample runs
to 10⁷–10⁸ events a year — against the ~10⁴ dimuon events that currently anchor
the strange-quark PDF.

## 4. One straight, two sign-tagged beams, fired at once

μ⁺ and μ⁻ counter-rotate in the same ring, so **each direction along a straight
carries one sign**. The south beam is μ⁺ decays; the north beam is μ⁻ decays.
Both fire simultaneously, from the same stored muons, through the same
apertures.

| | ν<sub>τ</sub> per tonne-year | ν̄<sub>τ</sub> per tonne-year | purity |
|---|---|---|---|
| north beam (μ⁻ decays) | **0.231** | 0.008 | 97 % ν<sub>τ</sub> |
| south beam (μ⁺ decays) | 0.014 | **0.137** | 91 % ν̄<sub>τ</sub> |

The minority component is not contamination from the beamline — it is
ν<sub>e</sub> → ν<sub>τ</sub>, a known amplitude. Add a magnetic field and the
τ charge is measured event by event, and the separation is complete.

**Nobody has ever done this.** Two experiments have identified ν<sub>τ</sub>
charged-current interactions individually: DONUT, with 9
([PRD 78, 052002](https://doi.org/10.1103/PhysRevD.78.052002)), and OPERA, with
10 ([PRL 120, 211801](https://doi.org/10.1103/PhysRevLett.120.211801)). Nineteen
events, in twenty-five years. **Not one of them was charge-tagged.** DONUT's
sample was an inseparable mixture of ν<sub>τ</sub> and ν̄<sub>τ</sub> from
D<sub>s</sub> decay; OPERA ran a ν<sub>μ</sub> beam, so its ten were
ν<sub>τ</sub>. The antineutrino of the third generation — a particle in the
Standard Model since 1975 — has never been individually identified.

A sign-tagged 91 %-pure ν̄<sub>τ</sub> beam is the only way to change that, and
a muon storage ring is the only way to make one.

## 5. The τ decay stops being an emulsion problem

The reason ν<sub>τ</sub> physics has 19 events is that the τ has to be caught
decaying, and at conventional beam energies it does not fly far enough to see
without photographic emulsion. At OPERA's 17 GeV the τ travels 0.50 mm; in
DUNE's CP-optimised beam, 0.24 mm.

This beam is a chirp from 63 GeV to 5 TeV, and the τ flight scales with energy:

| stage | ⟨E<sub>ν</sub>⟩ on axis | τ flight path |
|---|---|---|
| RCS1 | 110 GeV | **3.2 mm** |
| RCS2 | 350 GeV | 10.3 mm |
| RCS3 | 757 GeV | 22.3 mm |
| RCS4 | 2.0 TeV | 59.9 mm |
| collider store | 1.8 TeV | 53.1 mm |

Centimetres. A silicon vertex tracker resolves that; emulsion is not required,
which is exactly what lets the detector be tens of kilotonnes instead of tens
of tonnes. That single fact — **macroscopic τ flight at accelerator-scale
statistics** — is what this beam has and no other proposed beam does.

## 6. What it would pin down

**The ν<sub>τ</sub> and ν̄<sub>τ</sub> cross-sections.** The only direct
measurement is DONUT's, σ = (0.39 ± 0.13 ± 0.13) × 10⁻³⁸ E cm²/GeV — **47 %**.
The ν̄<sub>τ</sub> cross-section has never been measured at all. Every
ν<sub>τ</sub> appearance result, every IceCube flavour-ratio analysis and every
test of τ-flavour universality in neutrino scattering currently rests on a
number known to a factor of two. With a beam of *ab initio* known flux, this
becomes a percent-level measurement — and separately for each sign.

**Unitarity of the PMNS matrix, third row.** The τ row is by far the
least-constrained, and it is where heavy neutral leptons would show up. The
structure of the test is unusually clean here. Expanding the appearance
probability at small L/E and using unitarity,

<div style="text-align:center;font-family:monospace;margin:.8rem 0">
P(ν<sub>μ</sub>→ν<sub>τ</sub>) → | Σ<sub>i</sub> U*<sub>μi</sub> U<sub>τi</sub> Δ<sub>i</sub> |², &nbsp; Δ<sub>i</sub> = 1.267 Δm²<sub>i</sub> L/E
</div>

which vanishes as L² **only because Σ<sub>i</sub> U*<sub>μi</sub>U<sub>τi</sub> = 0.**
A non-unitary mixing matrix leaves a *zero-distance* appearance floor that does
not switch off as L → 0. The corridor delivers the ideal test bed: the same
beam seen at 2 km, 9.25 km and 198 km, where the Standard Model predicts a
strict L² scaling over four decades in probability, and the near hall — useless
for appearance physics on its own — has 10⁹ interactions a year with which to
bound the floor.

**sin²θ<sub>W</sub>, and the NuTeV anomaly.** NuTeV's 3σ deviation
([PRL 88, 091802](https://doi.org/10.1103/PhysRevLett.88.091802)) has never been
resolved, and its dominant systematics were exactly the ones this beam does not
have: the relative ν/ν̄ flux normalisation and the ν<sub>e</sub> contamination.
The Paschos–Wolfenstein ratio needs a ν beam and a ν̄ beam with a common
normalisation. **The corridor is that by construction** — north hall and south
hall, ν and ν̄, same stored muons, same cycle, with the only flux ratio between
them being the exactly-known 1/L² geometry.

**Strange-quark PDFs, charm production, and the structure functions**, from
10⁷–10⁸ charm events a year in a beam whose energy is known turn by turn.

## 7. What has to be true — and where it puts the far detector

The honest problem is signal-to-background. Everything terrestrial sits deep in
the (L/E)² regime — the first oscillation maximum for these energies is 5 × 10⁴
to 10⁶ km, four Earth diameters and up — so the τ appearance fraction is 10⁻⁷ to
10⁻⁴, and the τ sample hides under an unoscillated charged-current flux 10⁴ to
10⁷ times larger. This is the make-or-break number, and it has two levers.

**Lever one: distance.** The per-tonne signal is flat in L and the background
falls as 1/L², so **S/B ∝ L²**. This is the single most important design
consequence in this study, and it says the τ far detector does *not* belong at
UIUC:

| far hall | range | detector mass for 5σ in 5 yr, at OPERA's own τ efficiency and fake rate |
|---|---|---|
| UIUC | 198 km | 569 kt — not a proposal |
| W. Tennessee, due S | 655 km | 52 kt |
| Soudan (MINOS far hall) | 736 km | **17 kt** |
| NE Mississippi, due S | 830 km | 33 kt |
| Ash River (NOvA far hall) | 811 km | **14 kt** |
| SURF (DUNE far site) | 1289 km | 5.5 kt *(needs 101 mrad, beyond this study's verified range)* |

**Lever two: the chirp.** The 5 Hz cycle steps the stored energy from 63 GeV to
5 TeV in a known sequence, so a detector knows the beam energy at every instant
to within a turn — the [timing study](../timing/) already maps this. Because
P ∝ 1/E², the earliest turns carry a hundred times the τ fraction of the store.
At 830 km, gating on the turns below 100 GeV gives **S/B = 5 × 10⁻⁴, needing a
background rejection only 2.3× better than OPERA achieved**; the whole chirp
used as a weighted fit is better still.

Put together: **a 14–17 kt magnetised detector at Ash River or Soudan reaches
5σ in five years using nothing better than OPERA's demonstrated τ
identification.** Ten times better rejection — which the 6× longer τ flight
ought to buy — brings that to 10 kt with room to spare. This is a DUNE-scale
detector at a DUNE-scale baseline, in a hall that already exists.

**And the ring it needs is the cheapest one.** Breaking the sensitivity down by
stage, **RCS1 and RCS2 — the 6.28 km ring, the smallest in the complex —
supply 89.6 % of it** (RCS1 alone 62 %). The collider store contributes 4 %.
So the τ programme does not compete with the collider's aim at all; it asks
that *RCS1/2* be pointed at a far laboratory. Checking that against the
geometry this study has already established:

| | Soudan | Ash River |
|---|---|---|
| range, bearing | 736 km, 336° | 811 km, 335° |
| tilt | 57.8 mrad | 63.7 mrad |
| straight depth needed (6 m of arc-apex cover) | 69 m | 76 m |
| deepest point of the ring | 132 m | 145 m — inside the 200 m ceiling |
| co-tilt ceiling for this ring | 88.9 mrad | 88.9 mrad — no achromats needed |
| near beam emerges at | 1.19 km | 1.18 km — on site |

It closes, with margin, on a ring small enough to co-tilt bodily. The
[azimuth study](../rcs-azimuth/) already showed that above ~50 mrad the bearing
is free, and the [exit-pair map](../collider-map/) already showed that at these
tilts the up-going end is inside the fence on **every** bearing — so the
long-baseline τ configuration is also the *safest* one on site.

What is still owed, and what the case genuinely depends on:

* **A background simulation.** Everything above uses OPERA's fake rate as a
  benchmark. The dominant irreducible background here is charm, which at these
  energies also flies centimetres — the discriminants are the primary lepton
  veto, the charge correlation, and the kinematics, and none of that is
  simulated yet. This is the single piece of work that decides the τ programme.
* **A re-solve of the RCS1/2 geometry along bearing 335°**, with the real
  terrain, replacing the closed-form check in the table above.
* **A detector concept.** Magnetised, few-centimetre vertex resolution, 10–20 kt,
  tolerating ~10⁵ interactions per kt-year in the gated window.
* **The near-hall rate.** 706 Hz per tonne is a detector design problem, not a
  physics one, but it is not nothing.

## 8. The shape of the argument

A 10 TeV muon collider at Fermilab is justified by the energy frontier alone.
Its neutrino flux is, in every other treatment, the reason the project is hard.
Point it deliberately and the same flux is:

* the **only** route to an identified ν̄<sub>τ</sub>, from a 91 %-pure sign-tagged
  beam, reachable with a 14 kt detector in an existing hall;
* the **largest** neutrino dataset ever taken — NuTeV's entire run every 38
  minutes per tonne — and the only one where ν and ν̄ arrive together with a
  common normalisation;
* the **best-characterised** beam ever built, because its flux is muon decay and
  not hadroproduction;
* and therefore the instrument that measures the ν<sub>τ</sub>/ν̄<sub>τ</sub>
  cross-sections, tests third-row unitarity against a strict L² prediction, and
  re-opens sin²θ<sub>W</sub> with the systematics that broke NuTeV removed.

None of it asks the collider to be a different machine. It asks for a 6.28 km
ring to be tilted at a laboratory that is already built, and for the flux
everyone else buries to be aimed instead.
