---
title: "The Collider's Compass: Paired Exits, Tau Appearance, and Landmarks"
weight: 10
---

Three views of one frame — bearing around the Fermilab IP, log range radial.
The [azimuth-freedom map](../rcs-azimuth/) asked which bearings a steep beam
could take; these pages ask, for the whole compass, **where each tilt puts the
two ends of a straight**, **how much tau appearance each stage's beam carries
to each range**, and **what is actually out there**, including for the
collider in its default aim. Tools:
[`exit_pairs.py`](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/tools/exit_pairs.py),
[`tau_compass.py`](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/tools/tau_compass.py),
[`collider_map.py`](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/tools/collider_map.py);
machine-readable [exit_pairs.json](../geo/exit_pairs.json),
[tau_compass.json](../geo/tau_compass.json), [collider_map.json](../geo/collider_map.json).

## 1. Where the two ends of one straight surface, tilt by tilt

A straight at depth *d₀* tilted by θ has an end that goes **up** and an end
that goes **down**. The up end surfaces at s<sub>near</sub> ≈ d₀/θ; the down
end dives, and Earth's curvature brings it back to the surface at
s<sub>far</sub> ≈ 2 R<sub>E</sub> θ — on the opposite bearing. So every tilt
is a **pair** of exit curves, and the figure draws each pair in one colour:
dotted for the up-going end, solid for the down-going end. At zero tilt the
two coincide on a single curve. The bold curves are terrain-corrected with
GEBCO 2020 sampled along 48 bearings from 0.25 to 1000 km
([`fetch_gebco_bearings.py`](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/tools/fetch_gebco_bearings.py),
[data/gebco_bearings.json](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/data/gebco_bearings.json));
the faint circles are the smooth-sphere values. Exits over the Great Lakes are
taken at the water surface, not the bed.

<div style="margin:1.2rem 0"><img src="../figs/exit_pairs.svg" alt="Paired exit curves per tilt: dotted up-going end, solid down-going end, terrain-corrected, with the Fermilab fence and lakes" style="max-width:100%"></div>

| tilt | up-going end, smooth | with terrain (min–max over bearings) | down-going end, smooth | with terrain | up end inside the fence |
|---|---|---|---|---|---|
| **0 mrad** | 21.1 km | **11.8–37.3 km** (29 km due N, 16 km due S) | same | same | 0 of 48 bearings |
| 2 mrad | 11.9 km | 8.2–18.7 km | 37.4 km | 24–50 km | 0 of 48 |
| 5 mrad | 6.4 km | 4.8–9.9 km | 70 km | 55–81 km | 0 of 48 |
| **15.4 mrad** (collider) | 2.25 km | 1.9–3.0 km | 198.5 km | 195–206 km | 31 of 48 (off-site when the up end points E–NE, fence 0.9–1.5 km) |
| 25 mrad | 1.39 km | 1.2–1.5 km | 320 km | 316–325 km | 36 of 48 |
| 35 mrad | 1.00 km | 0.87–1.05 km | 447 km | 444–455 km | 45 of 48 (off-site only due E) |
| 50 mrad | 0.70 km | 0.63–0.72 km | 638 km | 635–643 km | **48 of 48** |
| 65 mrad | 0.54 km | 0.48–0.54 km | 829 km | 826–845 km | **48 of 48** |

Four readings:

1. **"Modulo terrain" is a factor of three at zero tilt.** A level 35 m
   straight surfaces anywhere from 12 to 37 km out depending on bearing:
   29 km due north, where the ground climbs toward the moraines, 16 km due
   south. For tilts below ~5 mrad the near exit moves by tens of percent
   with bearing and the far exit by ±30 %. Above 15 mrad terrain is a
   ±0.5 km effect on the near end and ±4 % on the far end: the smooth
   formula is then good enough, which is why the [aiming](../rcs-aim/) and
   [azimuth](../rcs-azimuth/) studies could use it.
2. **The fence test is on the plot.** The grey site is 0.9–4.9 km across
   from the IP, so a dotted curve inside the grey means the up-going end
   emerges on DOE land. Below 5 mrad no bearing does; at the collider's
   15.4 mrad, 31 of 48 bearings do — the failures are the east and
   north-east, where the fence is under 1.5 km; from 50 mrad on, every
   bearing does. This is θ<sub>free</sub>'s emergence half, drawn rather
   than tabulated.
3. **Two worked pairs.** The baseline collider: up end due north at
   2.2 km inside the fence, down end at UIUC, 198 km — one straight, two dots
   on opposite bearings, joined through the centre. The Lake Huron RCS aim
   (50 mrad, [azimuth study](../rcs-azimuth/)): down end at 637 km into the
   lake at bearing 52°, up end at 0.70 km on site at bearing 232°. The
   plot makes the coupling visible: choosing the far exit fixes the bearing
   of the near one.
4. **Which solid circle reaches which water** is now a single glance: 15.4
   mrad grazes Lake Michigan for 45° of bearing; 25 mrad crosses it fully;
   50 mrad reaches Superior, Huron and Erie; 65 mrad reaches Ontario and
   passes Ash River.

All curves are for a 35 m straight. The [RCS aiming](../rcs-aim/) study's
15 m straight halves every up-going range and leaves the down-going ones
essentially unchanged (s<sub>far</sub> depends on depth only through the
√(θ² + 2d₀/R<sub>E</sub>) term).

## 2. Tau appearance around the compass, stage by stage

<div style="margin:1.2rem 0"><img src="../figs/tau_compass.svg" alt="Five polar colour maps of the tau appearance fraction per stage, the curves, and the on-axis nu_tau CC rate per tonne-year against DUNE" style="max-width:100%"></div>

For each acceleration stage and for the 5 TeV store, colour is the
flux-weighted oscillation probability a detector **on the beam axis** would
see at that range:

<div style="text-align:center;font-family:monospace;margin:.8rem 0">
P̄(L) = ⟨ 0.95 sin²(1.267 Δm² L / E<sub>ν</sub>) ⟩<sub>spectrum</sub>, &nbsp; Δm² = 2.5 × 10⁻³ eV²
</div>

The RCS straights are 1/γ-sharp pencils, so a point on their axis sees the
hard forward spectrum (⟨E<sub>ν</sub>⟩ = 0.7 E<sub>μ</sub>). The collider's
plume is divergence-smeared to σ<sub>θ</sub> = 7/γ
([arXiv:2608.02718](https://arxiv.org/abs/2608.02718)), so a point on *its*
axis sees the whole-plume spectrum (⟨E<sub>ν</sub>⟩ = 0.35 E<sub>μ</sub>,
Gaussian-cut at E<sub>ν</sub>/E<sub>μ</sub> = 1/50) at 1/101 of the pencil
density. Each ramp is averaged with decays ∝ 1/E per turn. P̄ depends on
range alone, so each panel is a set of rings — contoured at 1, 2, 5 × 10ᵏ —
and the geography under them says what a given fraction would be delivered
to. Under each panel: the **ν<sub>τ</sub> oscillation baseline** for that
stage's energy (the first maximum of sin², L<sub>max</sub> = 496 km ×
E<sub>ν</sub>[GeV]), the fractions at UIUC and Soudan, and the on-axis rate.

| stage | ⟨E<sub>ν</sub>⟩ on axis | 1st oscillation maximum | P̄ at UIUC 198 km | P̄ at Soudan 736 km | P̄ at SURF 1289 km | **ν<sub>τ</sub>+ν̄<sub>τ</sub> CC on axis, per tonne·yr** |
|---|---|---|---|---|---|---|
| RCS1 (63→314 GeV) | 110 GeV | 54,000 km (4 Earth diameters) | **1.1 × 10⁻⁴** | 1.5 × 10⁻³ | 4.5 × 10⁻³ | 0.017 |
| RCS2 (314→750) | 351 GeV | 174,000 km (14) | 7.2 × 10⁻⁶ | 9.8 × 10⁻⁵ | 3.0 × 10⁻⁴ | 0.062 |
| RCS3 (750→1500) | 758 GeV | 376,000 km (29) | 1.4 × 10⁻⁶ | 2.0 × 10⁻⁵ | 6.1 × 10⁻⁵ | 0.041 |
| RCS4 (1500→5000) | 2036 GeV | 1.0 × 10⁶ km (79) | 2.5 × 10⁻⁷ | 3.5 × 10⁻⁶ | 1.1 × 10⁻⁵ | 0.108 |
| Collider store (5000) | 1806 GeV | 0.9 × 10⁶ km (70) | 1.3 × 10⁻⁶ | 1.8 × 10⁻⁵ | 5.2 × 10⁻⁵ | 0.14 |
| **all five, co-tilted chain** | | | | | | **0.37** |
| DUNE, CP-optimised beam ([TDR vol. II §4.1.1.3](https://arxiv.org/abs/2002.03005)) | ~3 GeV | 1,300 km *is* the baseline | | | | 130 / yr in 40 kt = **0.0033** |
| DUNE, τ-optimised horns (same source) | higher | | | | | 1000 / yr in 40 kt = **0.025** |

The rate is ν<sub>τ</sub>+ν̄<sub>τ</sub> CC per tonne-year for a point detector
on the axis: north-aimed decays per cycle from the
[timing study](../timing/) (5 Hz, 1.2 × 10⁷ s/yr, both muon signs), CSMS
charged-current cross-sections, the τ-threshold suppression table of the
[co-tilted chain](../uiuc-chain/) page. It is **baseline-independent**:
the on-axis flux falls as 1/L² while P̄ rises as L², so the product is flat
for every stage across the whole map — the (L/E)² cancellation the timing
study noted, now drawn. Cross-check: the store-only value, 0.14 per
tonne-year, reproduces the 0.137–0.140 the co-tilted chain page obtained
with an independent Monte Carlo and a 2 m probe; the full-chain total here
(0.37) is 25 % above that page's 0.29–0.30, the difference being that page's
probe averaging over the pencils' small spots.

Readings:

1. **Every baseline on the map is deep in the quadratic regime.** The
   first oscillation maximum for the softest stage is 54,000 km away — four
   Earth diameters — and for the store 900,000 km. So the maps are
   concentric rings whose spacing is the same for every stage, and moving
   from UIUC to Soudan buys the same ×14 in P̄ for all of them. No terrestrial
   baseline reaches the ν<sub>τ</sub> oscillation maximum for a beam this
   energetic; DUNE reaches its own because its neutrinos are 3 GeV.
2. **The stage ordering of P̄ is the whole story of the fraction.** At any
   range RCS1 carries 400× the tau fraction of RCS4 and 100× the store's,
   which is why the [co-tilted chain](../uiuc-chain/) gains so much and why
   sending RCS3/4 [elsewhere](../rcs-aim/) costs so little.
3. **The ordering of the on-axis rate is different, and it is the one that
   counts events.** Rate ∝ N<sub>decays</sub> × γ² × P̄ × σ ∝ N × E, so the
   high-energy stages win on axis despite their tiny fractions: RCS4 gives
   0.11 per tonne-year and the store 0.14, RCS1 only 0.017. The low-energy
   stages' advantage lives off axis — their spots are 30–300 m wide at
   198 km and their ν̄<sub>τ</sub> are spread through them — which is the
   argument for a *large-area* far detector rather than a dense one.
4. **The collider's smeared plume beats the RCS4 pencil on its own axis** in
   both fraction (×5) and rate, although RCS4 spans the same top energy: the
   MINT divergence softens the point spectrum. Per neutrino, not per fluence.
5. **Against DUNE.** Per tonne on axis, the full chain's 0.37 ν<sub>τ</sub> CC
   per year is **15× DUNE's τ-optimised beam and 110× its CP-optimised
   beam**; the store alone is 5.6× the τ-optimised value. DUNE's advantage is
   mass: 40,000 t. To equal its 130 CP-optimised interactions per year a far
   detector here needs **350 t on axis**; to equal the τ-optimised 1000 per
   year, **2,700 t** — inside a spot whose 50 % radius is 36 m for the
   collider and 4–330 m for the pencils, so a kiloton-class LAr or emulsion
   volume fits geometrically. Every one of those events is turn-tagged and
   energy-tagged by the [chirp](../timing/); none of DUNE's are.
6. **Bearing enters only through what sits at each ring.** On the 198 km
   circle — UIUC, Purdue, Kettle Moraine, 45° of Lake Michigan — every
   stage's fraction and rate are fixed; the choice among them is land use,
   not physics. Range buys fraction (×14 to Soudan) but not on-axis rate.

## 3. The landmarks around the collider's default aim

The collider's tilt is not a knob: **15.40 mrad up-to-north from a 35 m
straight** is fixed by the UIUC exit, and §1 shows what that means — its
down-going end surfaces at **198.4 km in every bearing**. The third map
puts the real geography on that circle and around it: lake footprints, DOE
and university partners, mines and quarries and the MINOS, NOvA and DUNE far
sites, national forests and state parks, Mt. Simon gas-storage fields,
airports, and the population centres the beam must not graze.

<div style="margin:1.2rem 0"><img src="../figs/collider_map.svg" alt="Polar landmark map around the collider: bearing around, log range radial, with the 198 km exit circle, theta_free, lake footprints and categorised landmarks" style="max-width:100%"></div>

Radius is log range so that 10 km neighbours and 1300 km mines share one
frame; every range ring is also labelled with the tilt that reaches it. The
green curve is the site's θ<sub>free</sub>(bearing) from the
[azimuth study](../rcs-azimuth/), drawn at its equivalent range 430–1220 km.
Sanity checks against known baselines: the map places Soudan at 736 km
(MINOS: 735 km, and NuMI's 58 mrad downward pitch is exactly
736 km / 2 R<sub>E</sub>), Ash River at 811 km (NOvA: 810 km), and SURF at
1289 km (DUNE: ~1300 km).

### What is on the circle

Rotating the ring about the IP without changing its tilt moves the far exit
around the 198 km circle and the near beam around the opposite bearing.
Per degree of bearing, from the JSON:

| far exit bearing | lands at | near side | fence there | off-site sub-500 ft band |
|---|---|---|---|---|
| **180° (baseline)** | **UIUC South Farms, 198.4 km** (Willard Airport 200.8 km at 181°) | N — ComEd ROW, Wayne | 2845 m | **9.3 km** |
| 145° | **Purdue, 191.9 km** — a second land-grant campus on the circle | NNW — Elgin side | 2936 m | 9.2 km |
| 11°–55° | **Lake Michigan open water**, 45° of arc (e.g. az 24°: 43.47° N, 87.22° W, mid-lake) | SSW–SW — Aurora, Sugar Grove | 3057–4453 m | 9.1–7.7 km |
| 0.5° | **Kettle Moraine SF (north unit), 195.2 km** — Wisconsin state forest | S — Wheatland Twp, Naperville | 2628 m | 9.5 km |
| 90° | southern Michigan farmland (41.82° N, 85.83° W) | W — the site's long axis | 4657 m | 7.5 km |
| 270° | the Mississippi valley at Clinton, Iowa | E — 0.9 km of DOE land | 900 m | 11.3 km |

The 198 km circle lies **entirely inside θ<sub>free</sub>**: at 15.4 mrad the
near beam needs (d₀ + 152 m)/θ = 12.2 km to clear 500 ft and the fence is
0.9–4.7 km away, so 7.5–11.3 km of the climb is over someone else's land
*whatever the bearing*. The corridor is unavoidable; the design question is
only *whose* land it crosses, and the ComEd right-of-way is the only bearing
that offers a utility easement rather than houses. Rotating 11–55° east buys a
Lake Michigan exit at the same tilt and pays with 9.1 km over Aurora and Sugar
Grove — the trade the baseline made, shown from both sides for the first
time. And the circle is a university circle: UIUC by design, Purdue at 192 km
by geography, Notre Dame and UW–Madison a few mrad inside.

### Under the two beams

**North, airborne past 2.2 km.** Smith Road at 10.1 km and Wayne at 12 km
are the houses the near beam climbs over at 380–500 ft — the corridor marks
on the plot. Beyond them the bearing passes Moraine Hills State Park (51 km)
and Chain O'Lakes State Park (68 km) with the beam already 0.9–1.4 km up,
then Kettle Moraine SF at 195 km, the Menominee Range iron mines at Iron
Mountain (442 km, 34.7 mrad), the Marquette Range (515 km, 40.4 mrad), and
Lake Superior (626–695 km, 49–54 mrad — the [meridian water
target](../rcs-aim/)). The two tick marks at 38 km and 28 km are the
[safety page's](../safety/) untilted "Zone D" default exits.

**South, underground to 198 km.** The chord passes under Midewin National
Tallgrass Prairie at 54 km (638 m deep), reaches its 790 m perigee at 98 km —
10 km west of the **Herscher** gas-storage field (88 km, chord at 783 m) —
and climbs past the **Manlove** field at 161 km (chord at 484 m, 2 km west of
the field centre) before surfacing at the South Farms. Both fields store gas
in the Mt. Simon sandstone, the horizon the [aquifer page](../aquifer/) puts
the perigee in; whether the chord meets an actual reservoir is a question for
the field records that page already lists as owed.

### The rest of the compass

| landmark | bearing | range | tilt to reach | why it is on the map |
|---|---|---|---|---|
| Soudan mine | 336° | 736 km | 57.8 mrad | MINOS far site — the precedent for a Fermilab beam aimed at a mine |
| Ash River | 335° | 811 km | 63.7 mrad | NOvA far site |
| SURF / Homestake | 288° | 1289 km | 101 mrad | DUNE far site; same tilt class as the Gulf of Mexico disqualifier |
| Iron Mountain, Marquette Range | 1.5°, 5° | 442, 515 km | 35, 40 mrad | working and former iron mines *on the corridor meridian* |
| S. Illinois coal basin, Shawnee NF | 187°, 185° | 463, 485 km | 36, 38 mrad | deep mines and 1,100 km² of national forest on the meridian's southward extension |
| Argonne NL | 127° | 25 km | — | the other DOE laboratory; 2 mrad |
| Thornton Quarry, Ottawa silica mines | 121°, 223° | 60, 75 km | — | the two large open excavations within 100 km |
| Chequamegon-Nicolet, Ottawa, Hiawatha NFs | 353°, 350°, 14° | 443–502 km | 35–39 mrad | the largest public land within 500 km, all north |
| Chicago Loop, O'Hare, Midway | 85°, 61°, 98° | 49, 31, 40 km | — | the sector to keep every beam out of |

Everything worth aiming at from a far-hall point of view — deep mines, big
water, big public land — sits at 440–800 km and 35–65 mrad, exactly the
θ<sub>free</sub> class where §1's dotted curves are inside the fence on every
bearing. The collider at 15.4 mrad can have none of it; the RCS straights
can — and §2 says the low-energy ones are the ones to send.

## 4. Open items

* GEBCO at 15″ (~460 m) sets the terrain resolution; the near-exit curves
  for tilts ≥ 25 mrad (0.5–1.5 km out) are at the limit of what it resolves
  and want a local DEM.
* The tau rates are for a point detector on the axis; a detector wider than
  the spot collects the whole plume, whose ν̄<sub>τ</sub> count grows as L²
  instead of staying flat — the aperture is the missing design parameter.
  The DUNE reference values are TDR interaction counts before detector
  efficiency; the corridor numbers carry no efficiency either.
* Landmark coordinates are gazetteer-grade (~0.01°), fine for this frame
  but not for an aim point. "Public land" is national forests, state parks
  and Midewin; county forest preserves and state wildlife areas are not
  drawn. Gas-storage reservoir depths and extents relative to the chord are
  not modelled.
