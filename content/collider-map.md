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
two coincide on a single curve. The bold curves are terrain-corrected along
144 bearings (every 2.5°): USGS 10 m NED inside 6 km, where the up-going
exits live, and GEBCO 2020 from 6 to 1000 km, the two grids tied together by
their +1.6 m mean offset in the overlap, smoothed with a periodic Catmull-Rom
spline
([`fetch_ned_near.py`](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/tools/fetch_ned_near.py),
[`fetch_gebco_bearings.py`](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/tools/fetch_gebco_bearings.py);
[data/ned_near.json](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/data/ned_near.json),
[data/gebco_bearings.json](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/data/gebco_bearings.json));
the faint circles are the smooth-sphere values. Exits over the Great Lakes are
taken at the water surface, not the bed. The three rings are drawn in plan
(collider, RCS3/4, RCS1/2: east straights on the meridian, bodies west),
and three chords are drawn for the next subsection: grey dashed, one straight
through the IP with both ends far (UIUC and Lake Superior); teal, the chord
through UIUC and Green Bay — dotted the true one from the South Farms to
mid-bay, which passes 11 km east of the IP, dashed the same two ranges as a
straight through the IP; olive dotted, the true chord through UIUC and Lake
Winnebago, which crosses the Fermilab site 3.9 km west of the IP.

<div style="margin:1.2rem 0"><img src="../figs/exit_pairs.svg" alt="Paired exit curves per tilt: dotted up-going end, solid down-going end, terrain-corrected, with the Fermilab fence and lakes" style="max-width:100%"></div>

| tilt | up-going end, smooth | with terrain (min–max over bearings) | down-going end, smooth | with terrain | up end inside the fence |
|---|---|---|---|---|---|
| **0 mrad** | 21.1 km | **11.0–37.2 km** (29 km due N, 16 km due S) | same | same | 0 of 144 bearings |
| 2 mrad | 11.9 km | 7.9–18.7 km | 37.4 km | 24–50 km | 0 of 144 |
| 5 mrad | 6.4 km | 4.1–9.9 km | 70 km | 55–81 km | 0 of 144 |
| **15.4 mrad** (collider) | 2.25 km | 1.7–2.5 km | 198.5 km | 195–206 km | 90 of 144 (off-site when the up end points E–NE, fence 0.9–1.5 km) |
| 25 mrad | 1.39 km | 1.2–1.6 km | 320 km | 316–325 km | 106 of 144 |
| 35 mrad | 1.00 km | 0.85–1.12 km | 447 km | 444–455 km | 132 of 144 (off-site only around due E) |
| 50 mrad | 0.70 km | 0.60–0.78 km | 638 km | 635–643 km | **144 of 144** |
| 65 mrad | 0.54 km | 0.46–0.57 km | 829 km | 826–845 km | **144 of 144** |

Four readings:

1. **"Modulo terrain" is a factor of three at zero tilt.** A level 35 m
   straight surfaces anywhere from 11 to 37 km out depending on bearing:
   29 km due north, where the ground climbs toward the moraines, 16 km due
   south. For tilts below ~5 mrad the near exit moves by tens of percent
   with bearing and the far exit by ±30 %. Above 15 mrad terrain is a
   ±0.5 km effect on the near end and ±4 % on the far end: the smooth
   formula is then good enough, which is why the [aiming](../rcs-aim/) and
   [azimuth](../rcs-azimuth/) studies could use it.
2. **The fence test is on the plot.** The grey site is 0.9–4.9 km across
   from the IP, so a dotted curve inside the grey means the up-going end
   emerges on DOE land. Below 5 mrad no bearing does; at the collider's
   15.4 mrad, 90 of 144 bearings do — the failures are the east and
   north-east, where the fence is under 1.5 km; from 50 mrad on, every
   bearing does. With the 10 m terrain the baseline's own up-going end sits
   at 2.13 km due north, 0.7 km inside the fence. This is θ<sub>free</sub>'s emergence half, drawn rather
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

### Both ends far? The depth it costs

The paired curves suggest a question: could the last RCS be angled so that
one end comes out at UIUC and the other in Green Bay or Lake Superior — both
ends far? A straight is a chord of the Earth, and its two exits obey
s<sub>near</sub> s<sub>far</sub> = 2 R<sub>E</sub> d₀. Fixing the south end
at UIUC (s₁ = 198.4 km) and asking the north end to surface at s₂ gives

<div style="text-align:center;font-family:monospace;margin:.8rem 0">
d₀ = s₁ s₂ / 2R<sub>E</sub>, &nbsp; θ = (s₂ − s₁) / 2R<sub>E</sub>, &nbsp; deepest point ((s₁+s₂)/2)² / 2R<sub>E</sub>
</div>

<div style="margin:1.2rem 0"><img src="../figs/two_ends.svg" alt="Depth of a straight whose south end is at UIUC as a function of where its north end surfaces, against real excavation depths" style="max-width:100%"></div>

| north end at | s₂ | straight depth at the IP | deepest point of the chord | tilt |
|---|---|---|---|---|
| the north fence (the baseline) | 2.9 km | **44 m** | 0.79 km | 15.4 mrad up-N |
| Cary / Fox River Grove (the level ring's exit) | 38 km | 0.60 km | 1.10 km | 12.6 |
| Kettle Moraine SF | 195 km | 3.0 km | 3.0 km | level |
| Green Bay's latitude on the UIUC line (Shawano Co., 15 km W of the bay) | 298 km | 4.6 km | 4.8 km | 7.8 down-N |
| **Lake Winnebago** (bearing 358°, see below) | 258 km | **4.0 km** | 4.1 km | 4.7 |
| **Green Bay water** (needs bearing 2.8°, see below) | 304 km | **4.7 km** | 5.0 km | 8.3 |
| **Lake Superior** (south shore / mid-crossing / north shore) | 626 / 655 / 695 km | **9.7 / 10.2 / 10.8 km** | 13.3 / 14.3 / 15.7 km | 34–39 |

For scale: this study's tunnels are 35–107 m deep; SNOLAB is 2.07 km, the
Gotthard base tunnel's greatest overburden 2.3 km, the deepest mine
(Mponeng) 4.0 km, the deepest borehole (Kola) 12.3 km. **Both ends far is
not a tunnel.** Even Kettle Moraine — the mirror-image exit at the same
198 km — needs the straight 3 km down; Green Bay needs 4.7 km, below every
mine but one; Lake Superior needs 10 km, below every mine. The physics is
just the sagitta: to have the up-going end travel 300 km before it reaches
daylight, it has to start 4.7 km below it. The tilt, incidentally, is
*small* for these chords — 8 mrad for Green Bay, and exactly level for
Kettle Moraine — which is why depth, not angle, is the cost.

Two geographic notes. Green Bay's water is not on the UIUC line: the bay is
first reachable at bearing 2.8°, 304 km out, and its open water (36 km of it
along the ray) at bearing 5.2°, 328 km. The polar plots draw that chord two
ways. The **true chord** from the South Farms to mid-bay is a 526 km arc that
passes **11 km east of the IP** (near Wheaton) at 5.1 km depth, 5.4 km at its
midpoint — a straight through both places is not at Fermilab. A straight
**through the IP** with the same two ranges is 5.1 km deep, tilted 10 mrad,
and puts its south end **18 km west of the South Farms** (9.7 km for the bay's
southern tip at 2.8°). In the IP-centred polar frame the true chord is the
curve that swings out to the east and back; the through-IP version is a pair
of opposite spokes.

**Lake Winnebago is the one lake that lines up.** It is the nearest water on
a UIUC line (258 km at bearing 358°), and the chord from the South Farms into
its water — choosing the lake point at least 2 km from shore that brings the
chord nearest Fermilab — passes **3.9 km west of the IP, inside the fence**,
at **4.0 km depth** (4.1 km at its midpoint, 458 km arc). UIUC, the west side
of the site and the lake are nearly collinear: the through-IP version of the
same chord is 4.0 km deep with a 4.7 mrad tilt and puts its south end only
6.9 km east of the South Farms; the chord to the lake's centre passes 6.9 km
west of the IP, just outside the fence. It is still a 4 km-deep straight —
the deepest mine on Earth — but of every both-ends-far chord on the compass it
is the only one whose two exits are water and a partner campus *and* whose
middle lies under Fermilab land. On the paired-exit plot it is the olive
dotted curve threading the west side of the grey site polygon. And the level-straight case
is the one the [safety page](../safety/) already knows: at 0.60 km depth the
UIUC-bound straight's other end would surface at Cary, the Zone D default
exit.

**The feasible version is two straights, not one.** A racetrack has two long
straights, and with the [bent-ring](../uiuc-chain/) vertical achromats each
can carry its own tilt: the east straight up-north at 15.4 mrad (down-going
end at UIUC), the west straight down-north (down-going end in the water) —
both at ordinary depth. For RCS3/4 (L<sub>s</sub> = 450 m, R = 2.2 km,
C = 14.7 km) aimed at Lake Superior:

| | west end in Lake Superior (655 km) | west end at Green Bay's range (300 km) |
|---|---|---|
| tilts, east / west straight | 15.4 / 51.4 mrad | 15.4 / 23.4 mrad |
| up-going ends (35 m straights) | 2.27 km N, 0.68 km S — both on site | 2.27 km N, 1.49 km S — both on site |
| arcs inclined by | 2.2 mrad | 1.3 mrad |
| vertical bending per lap | 142 mrad = 2,370 T·m at 5 TeV = **296 m of 8 T (2.0 % of the ring)** | 83 mrad = 1,380 T·m = 172 m of 8 T (1.2 %) |
| depth | free: 35 m as shown; at 80 m the UIUC straight's up end moves to 5.3 km, off site, so the UIUC straight stays shallow while the water straight may be deeper | same |

Both straights rise in the direction of travel, so the arcs must descend to
close the loop — a 2 mrad grade — and every straight end needs an achromat:
2(θ<sub>A</sub> + θ<sub>B</sub>) plus the arc slopes, 142 mrad per lap. That
is 4.6× the vertical bending of the co-tilted chain's bent ring (which
needed 15.4 mrad per end), the price of using both straights for physics
instead of one. One consequence for the next subsection: the sign a site
receives is set by which way μ⁺ circulates, and a particle southbound on the
east straight is northbound on the west one — so UIUC and Lake Superior would
receive the **same** sign from this ring, and the on-site up-going ends the
other.

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
to (white dotted and dash-dot: the UIUC–Green Bay and UIUC–Lake Winnebago
chords of §1, with their ends). Under each panel: the **ν<sub>τ</sub> oscillation baseline** for that
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

### ν<sub>τ</sub> versus ν̄<sub>τ</sub> at four sites, per tonne-year, against DUNE

μ⁺ and μ⁻ counter-rotate, so each direction along a straight carries one
sign: the beam going one way is μ⁺ decays (ν̄<sub>μ</sub> + ν<sub>e</sub>),
the other way μ⁻ decays (ν<sub>μ</sub> + ν̄<sub>e</sub>). A site therefore
sees mostly ν̄<sub>τ</sub> — from ν̄<sub>μ</sub> → ν̄<sub>τ</sub>, amplitude
0.95 — with a ν<sub>τ</sub> minority from ν<sub>e</sub> → ν<sub>τ</sub>
(amplitude sin²2θ₁₃ sin²θ₂₃ ≈ 0.050), or the mirror image. Which one is a
choice of circulation sense, made once. Taking μ⁺ southbound on the east
straights (UIUC receives μ⁺ decays; northbound beams from a separate ring or
sense carry μ⁻ decays), with the same on-axis machinery as above
([`nutau_sites.py`](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/tools/nutau_sites.py),
[nutau_sites.json](../geo/nutau_sites.json)):

<div style="margin:1.2rem 0"><img src="../figs/nutau_sites.svg" alt="nu_tau and nubar_tau CC per tonne-year at four sites for the store and the full chain, against DUNE, and signal-to-background versus baseline" style="max-width:100%"></div>

| site | L | arriving sign | ν<sub>τ</sub> CC / t·yr | ν̄<sub>τ</sub> CC / t·yr | store only (ν<sub>τ</sub> / ν̄<sub>τ</sub>) | all-flavour CC background / t·yr | S/B | collider spot r₅₀ |
|---|---|---|---|---|---|---|---|---|
| FNAL site, 2.0 km N (beam 4.5 m underground) | 2 km | μ⁻ (ν) | 0.231 | 0.008 | 0.088 / 0.003 | 2.2 × 10¹⁰ | 1 × 10⁻¹¹ | 0.36 m |
| **UIUC**, South Farms | 198 km | μ⁺ (ν̄) | 0.014 | **0.137** | 0.005 / 0.052 | 2.2 × 10⁶ | 7 × 10⁻⁸ | 36 m |
| **Green Bay** (bay water at bearing 2.8°) | 300 km | μ⁻ (ν) | **0.231** | 0.008 | 0.088 / 0.003 | 9.9 × 10⁵ | 2.4 × 10⁻⁷ | 54 m |
| **Lake Superior** | 655 km | μ⁻ (ν) | **0.231** | 0.008 | 0.088 / 0.003 | 2.1 × 10⁵ | 1.2 × 10⁻⁶ | 118 m |
| DUNE, ν mode ([TDR](https://arxiv.org/abs/2002.03005)) | 1300 km | ν (with ν̄ contamination) | ≈ 0.003 (130 / yr in 40 kt, mostly ν<sub>τ</sub>) | | | | ≈ 0.05 | 40 kt fiducial |
| DUNE, ν̄ mode (scaled ×0.6 from the TDR's background ratios) | 1300 km | ν̄ | ≈ 0.002 (≈ 80 / yr, mostly ν̄<sub>τ</sub>) | | | | ≈ 0.05 | |

Whole-plane τ-flavour neutrinos per year (one sign, E<sub>ν</sub> > 3.5 GeV),
which unlike the on-axis rate grow as L²: 1.4 × 10¹¹ at 2 km, 1.4 × 10¹⁵ at
UIUC, 3.2 × 10¹⁵ at Green Bay, 1.5 × 10¹⁶ at Lake Superior.

Readings:

1. **Per tonne on the axis the rate is the same at every far site**, as §2
   showed; what the sign choice changes is the total, because
   σ<sub>ν</sub> ≈ 2 σ<sub>ν̄</sub>: a site fed by μ⁻ decays sees
   0.24 τ CC per tonne-year (0.231 ν<sub>τ</sub> + 0.008 ν̄<sub>τ</sub>), one
   fed by μ⁺ decays 0.15 (0.137 ν̄<sub>τ</sub> + 0.014 ν<sub>τ</sub>). UIUC
   can have either by reversing the circulation — a free choice, made once
   for the whole complex.
2. **The minority sign is 5–9 % of the sample and comes from the
   ν<sub>e</sub>**, whose ×20 smaller appearance amplitude is partly repaid
   by σ<sub>ν</sub>/σ<sub>ν̄</sub> and a softer spectrum. A far hall that
   can tell τ⁻ from τ⁺ (a magnetised or emulsion detector) measures both
   ν<sub>μ</sub>→ν<sub>τ</sub> and ν<sub>e</sub>→ν<sub>τ</sub> in one beam.
3. **What baseline buys is signal-to-background, not rate.** The
   unoscillated CC background is 2 × 10⁶ per tonne-year at UIUC and falls as
   1/L²; S/B is 7 × 10⁻⁸ at UIUC, 2.4 × 10⁻⁷ at Green Bay, 1.2 × 10⁻⁶ at Lake
   Superior — every one of them four to six orders below DUNE's ≈ 5 %. τ
   identification here cannot rely on statistics; it has to come from
   topology, and TeV-scale τ leptons oblige: a 500 GeV τ flies 25 mm before
   decaying, resolvable in emulsion or pixels, where DUNE's few-GeV τ leptons
   travel tens of microns.
4. **The FNAL site is not a τ site.** On axis it has the same per-tonne rate
   as everywhere else, but the collider spot is 0.36 m across and the RCS
   pencils 4 cm to 3 m, so no tonne fits on the axis, and S/B is 10⁻¹¹. It
   is the near detector for flux, cross-sections and the turn-tagged
   spectrum, not for appearance.
5. **Against nominal DUNE, per tonne**: 70–80× in either sign (0.24 vs 0.0033
   in the ν-type comparison; 0.15 vs 0.0019 in the ν̄-type). DUNE's plan
   alternates modes over seven years and averages ≈ 0.0026 per tonne-year;
   the corridor delivers both signs simultaneously to opposite ends of every
   straight, each turn-tagged, at 2 × 10⁶ background events per tonne-year
   per far hall for a 100 m-scale spot. Equal yearly *counts* need
   330–870 t on axis.

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
teal and olive dotted curves are the UIUC–Green Bay and UIUC–Lake Winnebago
chords of §1, ending at their lake markers. The green curve is the site's
θ<sub>free</sub>(bearing) from the
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

* Terrain is USGS NED at 10 m inside 6 km and GEBCO at 15″ (~460 m) beyond;
  the two disagree by +1.6 ± 2.2 m where they overlap, and GEBCO is shifted
  onto NED. Bearings are sampled every 2.5° and spline-smoothed between
  samples; the far exits over 100 km away are still GEBCO-limited to a few
  hundred metres.
* The two-straight ring's 142 mrad of vertical bending per lap, with closed
  vertical dispersion in a rapid-cycling ring, is unverified lattice work; so
  is the 2 mrad arc grade's effect on the ring's other systems.
* The DUNE antineutrino-mode ν<sub>τ</sub> yield is a scaling of the TDR's
  neutrino-mode number, not a TDR value; the sign split of DUNE's samples is
  illustrative.
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
