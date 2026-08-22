---
title: Corridor-Aligned Siting of a 10 TeV Muon Collider at Fermilab
weight: 1
---

**Concept:** edit the published FNAL muon-collider siting concepts so that the
collider's interaction-point (IP) straight *and* the long straights of the
RCS acceleration chain all lie on a single true-north meridian. The
muon-decay neutrino plumes from every long straight then stack into one
narrow N–S corridor that leaves the site at its northeast corner, runs up the
ComEd *Aurora–Wayne* transmission right-of-way, and feeds an underground
neutrino experiment at **41°55′39.3″ N, 88°13′22.7″ W** — the most intense
TeV-scale neutrino beam ever contemplated, essentially for free.

All numbers below are produced by [`tools/corridor_layout.py`](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/tools/corridor_layout.py)
(pure-stdlib; re-run it after changing any assumption) and are stored in
[`geo/summary.json`](../geo/summary.json). Radiological implications are in
[the safety assessment](../safety/).

---

## 1. The geodetic coincidence that makes this work

Three facts, all verified against surveyed data (OpenStreetMap, retrieved
2026-07-30; SRTM terrain):

1. The requested detector site sits **96 m from the ComEd Aurora–Wayne
   138/345 kV corridor**, which runs almost exactly due north–south at
   longitude ≈ −88.2235 between the lab and the site.
2. The **northeast corner of the Fermilab site boundary is the vertex
   (41.869907, −88.222983)** — on the *same meridian* as the detector to
   within 1 m of longitude. A due-north beam from on-site passes over the
   fence exactly at the corner and immediately enters the utility corridor.
3. Along that meridian (**−88.222972**, fixed by the detector longitude) the
   site has **5.48 km of north–south extent** (lat 41.82065 → 41.86990) —
   enough for the straight sections and the arcs of every ring.

Because a due-north geodesic keeps constant longitude, "final straightaways
pointing N–S such that the beam reaches the experiment" means exactly: *put
every long straight on the meridian −88.222972 at azimuth 0°/180°*.

## 2. Baseline layout

All rings are racetracks whose **east straight lies on the corridor
meridian**; ring bodies extend west into the site. Depths are below local
grade (~225 m ASL), in the dolomite bedrock — the same depth class as the
existing NuMI/MINOS underground areas, i.e. established FNAL construction
practice, not novel civil engineering.

| Ring | Role | C | Arc radius | Straights | Depth | East-straight location | Ring-plane pitch |
|---|---|---|---|---|---|---|---|
| **Collider** | 2 × 5 TeV (√s = 10 TeV) | 11.00 km | 1 528 m | 2 × 700 m | 100 m (reference; §5 for the shallow civic-envelope option) | lat 41.841129 – 41.847431, **IP at 41.844280** | 0 (horizontal; see §5 and Safety §2 for the pitch knob) |
| **RCS3 / RCS4** | 0.75 → 1.5 → ~4.5 TeV | 14.72 km | 2 200 m | 2 × 450 m | 80 m | centred lat 41.840768 | 2.16 mrad down-to-N |
| **RCS1 / RCS2** | 63 → 314 → 750 GeV | 6.28 km | 841 m | 2 × 500 m | 60 m | centred lat 41.844280 (stacked 40 m above the IP straight) | 4.32 mrad down-to-N |

The collider is placed 350 m south of the northmost fit so that **both** ends
of the IP straight — and both surface-grazing emergences of the tilted
civic-envelope configuration (§5) — stay comfortably inside the boundary.
The 11.0 km circumference (R_arc = 1 528 m, 2 × 700 m insertions) is the
"give or take" middle of the 10–11.5 km range the site supports on this
alignment; its **west straight falls on meridian ≈ −88.2598°**, which is
where the second interaction point lives (§5).

Notes:

* **Fit is against the real boundary polygon** (OSM way 31974155) with a
  150 m setback for the collider and 100 m for the RCS tunnels; the
  placement code slides each racetrack as far north as it fits and, for
  RCS3/4, maximises circumference over (R, L_s).
* The small **down-to-north pitches** of the RCS planes make their straight
  plumes *converge on the same detector hall* as the collider plume despite
  the different tunnel depths — and because the RCS straights are FODO (no
  final focus), theirs are the corridor's only true 1/γ pencils; the chain
  they deliver to the hall is quantified in [the timing study](../timing/).
  LEP/LHC's plane is tilted 14 mrad, so 2–5 mrad is well within precedent.
* Tunnels at different depths cross in plan freely (60 / 80 / 100 m); the
  stacked straight sections over the IP form a single "corridor vault"
  that can share access shafts.
* Injector complex (proton driver, target, cooling channel, low-energy
  recirculating linacs to 63 GeV) is unconstrained by the corridor and stays
  near the existing Main-Injector/PIP-II campus, with a transfer line to the
  RCS1/2 straight.

## 3. What was *edited* relative to the published FNAL siting concepts

The reference concepts (Muon Collider Forum report, arXiv:2209.01318; IMCC
interim report, arXiv:2407.12450; FNAL site-filler studies) fit rings to the
site for maximum circumference, with straight-section azimuths falling
wherever the fit puts them. The corridor design makes five deliberate edits:

1. **Collider ring rotated and translated** so the IP insertion is at
   azimuth 0° on meridian −88.222972, IP at 41.844280 N — instead of an
   orientation-agnostic centrally-placed ring. C = 11 km (within the
   published concepts' 10 km-class envelope).
2. **RCS1/2 moved out of the Tevatron tunnel** into a new 6.28 km racetrack
   (same circumference, so the published optics carry over) at 60 m depth
   with its straight on the corridor. *Rationale:* the Tevatron tunnel is
   ~7 m deep and its six straights point at fixed, unhelpful azimuths; its
   plumes would surface ~9.5 km out in six directions through Warrenville /
   Batavia / West Chicago at up to ~0.5 mSv/yr (wobbled). Tevatron reuse
   remains a documented fallback for a staged programme (see [Safety](../safety/) §6).
3. **The 16.5 km "site filler" for RCS3/4 becomes a 14.72 km racetrack**
   with N–S straights. This is the one real machine cost of alignment: the
   racetrack constraint plus boundary setbacks forfeit **10.8 %** of
   circumference. Options: (a) raise average bending field/packing ~11 %,
   (b) accept ~4.5 TeV per beam (√s ≈ 9 TeV), or (c) keep the
   boundary-hugging 16.5 km shape and give up straight alignment for RCS4
   only (its wobbled straights contribute ≤0.8 mSv/yr — see [Safety](../safety/) §2).
4. **Ring planes pitched** (0 / 2.16 / 4.32 mrad in the deep reference;
   15.40 mrad for the collider in the civic-envelope configuration, §5) to
   focus all straight plumes on one detector hall.
5. **Deep siting** (60–100 m, bedrock) replaces FNAL's traditional ~7–10 m
   cut-and-cover, which is what keeps every plume underground for tens of
   kilometres and pushes the surface-grazing exits far from the site.

## 4. The neutrino beamline and detector

**Geometry.** From the IP the plume centreline runs due north, ~100 m deep
at the fence in the deep reference layout, passing under the ComEd
right-of-way. Earth curvature plus the terrain rise toward the detector puts
the centreline **107 m below grade** at the experiment (grade 238.5 m ASL →
hall floor near 130 m ASL), **9.25 km** from the IP. See
[the corridor profile](../figs/corridor_profile.svg). (In the civic-envelope
configuration of §5 the same detector site can instead be served by a
near-surface hall — the [tool](../tool/) baseline puts a detector 15 m deep
at 41.8560° N, 1.3 km from the IP, still on DOE land.)

**What the detector sees** (per 1.2 × 10⁷ s Snowmass-year, IMCC-class beam:
1.8 × 10¹² μ/bunch/sign at 5 Hz, 90 % chain transmission):

| Quantity | Value (current lattice, f = 0) | Value (dedicated drift, f = 0.49 — [beam-size](../beamsize/) §4) |
|---|---|---|
| Decays aimed north in the IP straight | 5.3 × 10¹⁸ /yr — each emitting **one νμ and one ν̄e** (or the CP mirror, selectable by circulation direction) | same |
| Core radius at detector (50 % flux) | **1.7 m** (σθ = 0.15 mrad, [arXiv:2608.02718](https://arxiv.org/abs/2608.02718)) | **0.20 m** (1/γ pencil) |
| Core fluence | 4.4 × 10¹³ ν/cm²/yr **per species** | 2.2 × 10¹⁵ |
| Fluence-mean Eν | ≈ 1.6 TeV (blend; the energy–angle prism is washed out) | ≈ 3.2 TeV on axis |
| Interaction rate (CSMS σ, both species) | **5.5 × 10⁵ per kg-year** (5.5 × 10⁸ per tonne-year) | 5.2 × 10¹⁰ per tonne-year |

Even a ~10-tonne instrumented target collects ~5 × 10⁹ TeV-scale ν
interactions per year on the current lattice — six orders of magnitude
beyond FASERν, with both a νμ and νe component and known parent kinematics —
and returns to ~5 × 10¹¹ (eight orders) if the dedicated pencil drift of
[the beam-size study](../beamsize/) §4 is built. The 50 %-containment radius
is 1.7 m and 99 % of the flux falls within 4.6 m (MINT,
arXiv:2608.02718; the old 0.20 m pencil figure survives only in the
dedicated-drift scenario). The aligned RCS straights are not a perturbation:
their plumes are true 1/γ pencils converged on this same hall by the ring
pitches, and they add **~40 % of the collider's interaction rate** as a
63 GeV → 5 TeV energy chirp in the first 7 ms of each cycle — see
[the timing study](../timing/). With the RCS movers on, those bands sweep
±9 m. A hall of order 12 m (E–W) × 25 m (vertical) therefore captures
everything, with the instrumentation rather than the beam setting its size.
Physics case reference: "The Neutrino Slice at Muon Colliders",
arXiv:2412.14115.

**Civil.** The hall is a ~110 m-deep cavern with a surface shaft on a parcel
adjacent to the ROW (the given coordinate is 96 m from the line — outside
the conductors' fall zone; ComEd easement coordination is the same class of
negotiation as any deep utility crossing). The plume passes ~90–110 m under
Wayne / West Chicago parcels between the site and the hall — no surface
facility is needed anywhere along the corridor except at the detector.

## 5. The civic-envelope configuration and the second interaction point

The deep reference layout above keeps every plume ~100 m underground for
tens of kilometres, at the price of two far-away surface-grazing exit strips
(Safety §3). The **civic-envelope configuration** — the baseline of the
[interactive tool](../tool/), the [envelope note](../envelope/), and the
[aquifer audit](../aquifer/) — instead *chooses* both emergence points:

* **IP straight raised to 191 m ASL** (34 m below grade) and **tilted up
  15.40 mrad** (0.88°) toward north. This is the *minimum* tilt for which
  both beams pierce the ground only on institutional land:
  * the **north beam emerges at 41.8640° N — on the Fermilab site**, 655 m
    inside the fence, crosses the boundary 8 m up, climbs the ComEd ROW,
    passes Smith Road 116 m (380 ft) above grade — 350 ft above a 10 m
    roofline — and enters
    federally navigable airspace (500 ft AGL) 9.3 km past the fence;
  * the **south beam dives under Aurora** (74 m deep at the boundary,
    below 152 m/500 ft within 5.6 km), crosses under east-central Illinois
    at up to 789 m depth, and **surfaces 198 km away at 40.067° N on the
    UIUC South Farms** — university land, with the approach under Urbana
    66–119 m deep.
* Anything *less* tilted pushes the north emergence off-site into the ROW
  (still workable, but no longer on DOE land); anything more tilted is
  farther beyond LEP's 14 mrad plane-tilt precedent than necessary.
* **Ring options for the 15.40 mrad plane:** (a) a *planar tilted ring* —
  the natural choice; across the 3.06 km ring width the plane rises/falls
  ±23.5 m, so the north-arc apex needs ~6 m of engineered cover (berm —
  exactly how the original Main Ring was covered) while the south arc is
  63 m deep; or (b) a planar horizontal ring at depth with **vertical
  chicanes** (double-bend "doglegs") bringing only the IP insertion onto
  the tilted line — more optics work, standard geometry.

**The second interaction point** sits diametrically opposite, at the centre
of the **west straight, on meridian −88.2598°** (2.6 km west of the corridor
— under the west side of the site). In the planar tilted ring the IP2
straight is *automatically* tilted 15.40 mrad the same way (up-to-north). Its
beams do not follow the ComEd corridor; the audit of that meridian gives:

* **north beam:** emerges on-site or immediately north into
  Fermilab-adjacent open land west of the villages, climbing on the same
  schedule as IP1's beam (500 ft AGL within ~9 km); the meridian passes over
  the West Chicago/Wayne rural fringe rather than a subdivision line;
* **south beam:** re-surfaces ~198 km south at ≈ 40.045° N, −88.2598° —
  which is **Willard Airport (CMI), owned and operated by the University of
  Illinois**. A **+0.13 mrad trim** of the IP2 straight (well inside the
  ±1 mrad mover budget the IMCC already carries for exactly this purpose)
  centres that exit on airport land. An exit strip on a
  university-controlled airfield is the same land-use class as the South
  Farms solution: institutionally owned, fenced, and monitorable.

Both IPs therefore satisfy the same rule: **every point where a collider
plume is within 500 ft of the ground surface is on DOE, easement, or
University of Illinois land.** The RCS straights stay in the deep reference
configuration (their plumes exit at 20–80 km with the mover/wobble
mitigation of Safety §2) — only the collider ring is tilted.

## 6. Open items

* Full lattice check: racetrack optics with 700 m insertions at 11 km
  circumference; vertical-dogleg utility straights; combined-function arcs.
* MARS/FLUKA confirmation of everything in the [safety assessment](../safety/) (the analytic model is
  deliberately conservative).
* Geotechnical: Maquoketa/Galena–Platteville profile along the corridor and
  at the hall; groundwater at the NE corner; the Sandwich Fault Zone
  crossing at ~41.55–41.65° N (30–150 m throw, right where the south chord
  runs 300–500 m deep); well-log control on the dipping stratigraphic model
  (the south exit path is drift over Pennsylvanian shale, not dolomite —
  see the [aquifer audit](../aquifer/) §1).
* Land/easement engagement plan for the two collider-straight exit strips
  ([Safety](../safety/) §3) and the ROW crossing.
