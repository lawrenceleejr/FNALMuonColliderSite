---
title: "Design Variant: the Co-Tilted Chain — Everything to UIUC"
weight: 7
---

The baseline pitches the RCS ring planes *down*-to-north (4.32 / 2.16 mrad)
to converge their pencils on the deep-reference hall at 9.25 km — a choice
made before the civic-envelope/UIUC scenario existed, and never revisited.
This page revisits it. The variant: **co-tilt every ring up-to-north like
the collider, so that every east-straight south pencil exits at the
collider's own UIUC South Farms point** — putting the full 63 GeV → 5 TeV
chirp, with its (L/E)²-enhanced oscillation probability, onto the far hall.
Solved by [`tools/uiuc_chain.py`](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/tools/uiuc_chain.py)
with the corridor terrain model
(machine-readable: [uiuc_chain.json](../geo/uiuc_chain.json)).

**Headline:** the geometry closes — every ring has a feasible (depth, tilt)
pair that exits south at UIUC, emerges north *inside the Fermilab fence*,
and keeps the whole ring underground — but for the big RCS3/4 ring the *planar* solution
closes by **25 µrad of tilt**, the tightest tolerance in the whole study —
a tolerance §1 shows dissolves entirely if the ring is allowed one pair of
vertical achromats (tilted straight, level arcs).
The payoff at the far hall: **on-axis flux ×2.0** and **ντ CC ×2.1 per
tonne-year** with slightly *improved* S/B, plus a ×9 jump in CC-able
oscillated ν̄τ crossing the plane. The costs: one new ~3 mSv/yr pencil
component in the UIUC strip, three more LEP-scale ring tilts, and a
25 µrad plane-aim requirement.

## 1. The geometry: one constraint family, two hard walls

Each ring is a planar racetrack — two N–S straights joined by semicircular
arcs — so its plane spans y<sub>c</sub> ± (L<sub>s</sub>/2 + R) north–south,
and a tilt θ lifts the **north arc apex** (L<sub>s</sub>/2 + R)·θ above the
straight. Per ring, (z₀, θ) must satisfy three things at once:

* **A — south exit at UIUC:** the beam line passes through the collider's
  own exit (40.0602° N, 198.2 km, 221 m ASL in this model's frame). This
  ties z₀ to θ: one-parameter family.
* **B — north emergence on site:** shallower straights emerge sooner; the
  fence (41.8699° N) caps how *deep* the straight can be.
* **C — ring underground:** the tilted plane's north arc apex must keep
  cover (≥ 6 m here); this floors how *shallow* the straight can be.

The straights sit **side-by-side with the IP straight** (10–25 m E–W
offsets), not stacked — tunnels never intersect, the rings stay nested in
plan, and the E–W offsets reappear as metres of lateral spread inside the
same exit strip. Solutions, window midpoints:

| ring | tilt (up-N) | straight depth | north emergence | inside fence | N-arc apex cover | crosses near-hall plane at | feasible tilt window |
|---|---|---|---|---|---|---|---|
| collider (baseline, for reference) | 15.40 mrad | 33.6 m | 41.8648° N | 566 m | **5.9 m** | 14.1 m deep | — |
| RCS1/2 | 15.41 mrad | 32.6 m | 41.8642° N | 633 m | 15.9 m | 13.1 m | 15.36–15.45 mrad |
| RCS3/4 | 15.31 mrad | 45.1 m | 41.8680° N | **211 m** | 8.9 m | 19.3 m | **15.300–15.325 mrad** |

<div style="margin:1.2rem 0"><img src="../figs/uiuc_chain_profile.svg" alt="North-side elevation profile: three co-tilted beams emerging on site" style="max-width:100%"></div>

**Where the north beams emerge** — the question this variant must answer —
is: *all on site, in one zone.* The collider and RCS1/2 pencils surface
together at 41.864° N (566/633 m inside the fence, at the baseline's
existing emergence strip); RCS3/4, forced deeper by its arc-apex
constraint, surfaces at 41.8680° N — **211 m inside the fence**, the
variant's tightest land margin. All three then climb the same ComEd
corridor the baseline's civic-envelope walk already audits, laterally
offset by tens of metres.

Two structural notes. First, RCS3/4's whole existence squeezes through a
**25 µrad tilt window**: 15.300 mrad puts its north arc at 6 m cover;
15.325 mrad puts its north emergence at the fence. Second, the collider's
own north arc apex sits at **5.9 m cover in the published baseline** — the
civic envelope already lives with a shallow north arc; this variant adds
two more of the same species, not a new species.

### Do the straights need different depths?

RCS1/2 doesn't: it already sits at collider depth (32.6 vs 33.6 m).
RCS3/4's ~11 m offset is pure planar-ring arithmetic — its arc-apex lift
(L<sub>s</sub>/2 + R)·θ = 37.3 m exceeds the collider's 29 m — and scanning
the exit family shows the offset is **irreducible under planarity**: the
straight depth trades 1 m for 1 m of apex cover, and even with the north
arc *at grade* (zero cover) the planar RCS3/4 straight cannot come
shallower than 36.4 m, still 2.8 m below the collider:

| required apex cover | planar RCS3/4 straight depth | offset vs collider | fence margin |
|---|---|---|---|
| 8 m | 44.3 m | +10.7 m | 269 m |
| 6 m | 42.3 m | +8.7 m | 404 m |
| 4 m | 40.3 m | +6.7 m | 539 m |
| 0 m (arc at grade) | 36.4 m | +2.8 m | 808 m |

The escape hatch is to drop planarity where nothing needs it: only the
**east straight** must lie on the UIUC line — the arcs radiate sideways
and can sit wherever civil engineering likes. A **bent ring** (tilted
straight, level arcs, joined by vertical achromats at the straight's two
ends) detaches the apex constraint entirely: the straight then sits at
**any point on the exit family** — e.g. 39.7 m / 15.34 mrad to cross the
near-hall plane at exactly the collider beam's depth (fence margin 582 m),
or ~34 m for a ~1 km margin — and the 25 µrad window opens into a
±40 µrad-plus comfort zone set only by siting preference. The hardware is
standard geometry-matching: 15.4 mrad of vertical bend per end,
B·L = 257 T·m at 5 TeV (~32 m of 8 T per bend, achromat pairs to close
vertical dispersion), ~120 m of vertical dipoles per ring = 0.8 %
(RCS3/4) / 1.9 % (RCS1/2) of the circumference. Decays inside the
achromats spray over the 15.4 mrad vertical fan — self-segmenting, so
dose-benign. The cost is real lattice work (ramped vertical bends with
closed dispersion in a rapid-cycling ring, where the IMCC baseline is
planar), but vertical doglegs are already contemplated for the mover and
segmentation systems; this is more of the same hardware doing geometry
instead of mitigation.

## 2. What the far hall gains

With the chain's south pencils on the UIUC axis, the far hall receives the
chirp — and at 198 km, the chirp's low-energy turns carry oscillation
probabilities the 5 TeV store cannot reach: P ~ 10⁻³–10⁻² in the soft
turns, against 5 × 10⁻⁶ for the store.

| at the UIUC far hall (on axis) | baseline (store only) | **co-tilted chain** |
|---|---|---|
| total flux | 1.4 × 10¹¹ ν/cm²/yr | **×2.04** |
| oscillated ν̄τ areal density | 0.0060 /cm²/cycle | 0.0082 |
| **ν̄τ CC per tonne-year** (τ-threshold included) | **0.140** | **0.293** |
| ν̄μ CC per tonne-year (background) | 4.6 × 10⁵ | 7.8 × 10⁵ |
| S/B | 3.1 × 10⁻⁷ | **3.8 × 10⁻⁷** |
| oscillated ν̄τ crossing the plane, E > 5 GeV | 1.0 × 10¹⁴ /yr | **9.7 × 10¹⁴ /yr** |

<div style="margin:1.2rem 0"><img src="../figs/uiuc_chain_nutau.svg" alt="Oscillated nutau arrival rate at the far hall: baseline vs co-tilted chain" style="max-width:100%"></div>

Three readings of that table:

1. **The chirp doubles the far hall's τ-appearance rate and does not
   dilute it** — S/B improves, because the chirp's ν̄τ ride on lower-energy
   ν̄μ whose backgrounds are σ-suppressed faster than the τ-threshold
   suppresses the signal.
2. **The pencil prism spreads the real bounty off-axis.** The ×9.3 in
   plane-crossing ν̄τ (E > 5 GeV) lands in a wide low-energy pancake
   (the high-P neutrinos sit at large radius of each pencil's cone), which
   argues for a *large-area, low-threshold* far detector — emulsion-class
   or a segmented array across the strip — rather than a single on-axis
   tonne.
3. **Every event stays turn-tagged.** The chirp arrives in its own time
   gates (658 µs of flight later, an offset), so the far hall becomes an
   (L/E)-scanned, time-stamped appearance experiment — the thing no other
   accelerator geometry offers.

## 3. What it costs

* **The UIUC exit strip gains a pencil component.** King-raw, unwobbled,
  bounding each stage at its top energy: RCS4 adds **3.0 mSv/yr in an
  ~11 m-wide strip** (RCS1–3 add < 0.04 combined; E⁴). This sits on top of
  the collider's 0.84 mSv/yr over 74 m — the same fenced-strip /
  occupancy / segmentation trade as the collider's dedicated drift, now
  with a second knob. A programmed mover sweep (±1 mrad) dilutes RCS4 to
  ~0.1 mSv/yr at the cost of smearing the chirp over ±200 m — the
  [timing study](../timing/) §2 trade, verbatim.
* **The on-site north zone works harder.** The RCS pencils' emergence
  streaks carry Sv-class raw doses (RCS4: ~15 Sv/yr in a decimetres-wide
  streak at 2.8 km) — inside the same fenced, occupancy-controlled
  emergence zone the collider's own streak already requires; a wider fence,
  not a new class of problem.
* **Three more LEP-scale tilts and a 25 µrad aim.** Ring-plane tilts of
  15.3–15.4 mrad (precedent: LEP 14, the collider itself 15.4) on every
  ring, with RCS3/4's plane aimed to 25 µrad — and ±0.1 mrad of tilt moves
  a south exit ±1.3 km, so per-ring trim is a survey deliverable of the
  same class as the collider's own targeting claim.
* **Zone C disappears — a gain, not a cost:** the baseline's RCS exit
  strips at 12–40 km over third-party land are eliminated; every chain
  exit is on-site (north) or inside the UIUC strip (south). The west
  straights carry no physics and keep the baseline's utility treatment
  (doglegs + movers; their far bands dilute ~×40).

## 4. Open items

* **The tilted-arc disks.** A tilted ring's arc "disk" grazes the surface
  at azimuth-dependent ranges — near-north azimuths exit close, east–west
  at ~38 km, south at ~198 km. This needs its own mover-diluted dose study
  — and honestly, the civic-envelope baseline already owes it for the
  collider's own tilted arcs; this variant triples the need.
* Injection/extraction lines between co-tilted rings; RF and impedance of
  15 mrad-tilted RCS planes (assumed benign, unverified).
* The τ-threshold factor here is an approximate DIS suppression table;
  the ν<sub>e</sub>→ν<sub>τ</sub> channel (~5 %) and rock-produced ν̄τ
  backgrounds ([beam-size study](../beamsize/) §5) are not included.
* All doses are King-model raw/bounded, pre-FLUKA, as everywhere in this
  study.

*Sequel:* [aiming the last RCS](../rcs-aim/) takes the depth and angle apart
as independent knobs and finds that RCS3/4 — this variant's tightest
tolerance *and* its only new dose problem — is also the ring UIUC needs
least, and can be sent 656 km into Lake Superior instead at a cost of 3 % of
the plane-crossing ν̄τ gain.

**Verdict.** If ντ appearance is the far hall's purpose — and it is the
one measurement this corridor can do that nothing else can — the co-tilted
chain is the better design: twice the CC rate, nine times the plane-level
oscillated sample, no new off-site land, one dose knob it already had, and
one genuinely hard tolerance (RCS3/4's 25 µrad planar window) as the
price of admission — or ~120 m of vertical achromat per ring to buy the
tolerance away and put every straight at the same depth.
