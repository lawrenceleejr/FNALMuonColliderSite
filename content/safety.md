---
title: Neutrino Radiation Safety Along the North-South Corridor
weight: 2
---

This note quantifies the off-site radiological picture for the
corridor-aligned layout in [the design note](../design/), for the communities that
live above and beyond the beam. Numbers come from
`tools/corridor_layout.py`; the model and every correction factor are stated
explicitly so they can be checked and later replaced by full MARS/FLUKA
transport.

**Headline:** everywhere *along* the corridor the beam is 90–170 m
underground and delivers no measurable surface dose. The radiological
question is confined to two narrow, predictable **surface-grazing exit
strips** tens of kilometres away — a few metres wide — which must be brought
under easement/control, exactly as the international design (IMCC) solves the
same problem by aiming its insertion plumes at the Mediterranean and a fenced
mountainside. Everything else (arcs, RCS straights, utility straight) is
driven below or near the 0.1 mSv/yr design goal by standard mitigation.

## 1. Physics of the hazard

Muons decay in flight; each straight section emits a neutrino "pencil" along
its axis with half-angle 1/γ = **21 µrad** at 5 TeV. Neutrinos themselves are
essentially harmless to people they pass through; the dose arises where the
pencil **grazes the ground surface**, because neutrino interactions in the
last tens of metres of soil produce hadronic/EM showers that reach people at
the surface. Between the ring and that exit point the pencil is deep
underground and the showers it makes are absorbed by rock far below anyone's
feet: **surface dose above a buried plume is negligible; only the exit
regions matter.** Dose scales like E⁴ × (straight length) for straights and
E³ for arcs — which is why a 10 TeV machine's straights dominate its entire
radiological design.

Model: B.J. King, *Potential Hazards from Neutrino Radiation at Muon
Colliders*, arXiv:physics/9908017 — Eq. 10 (straights), Eq. 7 (arc disk),
"equilibrium approximation" (conservative: assumes full shower buildup).
Beam: 1.8 × 10¹² μ/bunch/sign, 5 Hz, 1.2 × 10⁷ s/yr, 90 % transmission →
0.97 × 10²⁰ decays/yr/sign in the collider.

Corrections applied (each labelled in the output):

* **σν flattening** above ~1 TeV vs the linear low-energy extrapolation: ×0.5.
* **Vertical segmentation** of insertion straights (deliberate ±0.5 mrad
  vertical doglegs splitting the non-IP part of each insertion): peak
  dilution = smeared band ÷ pencil/shower width at the exit range (~×20–40).
* **Mover ("wobbling") system**, ±1 mrad, per IMCC interim report
  (arXiv:2407.12450) — applied to arcs, utility and RCS straights, *not* to
  the IP straight (it would destroy luminosity).
* Effective collinear ("pencil") length of the IP insertion:
  **150 m baseline / 30 m optimistic** per direction — the 700 m insertion is
  mostly final-focus quads and chromatic-correction dipoles whose divergence
  and bending smear the plume; only truly field-free collinear drift makes
  the sharp core. This parameter dominates the result and is the first thing
  a real lattice + FLUKA study must pin down.

Regulatory anchors: DOE public dose limit **1 mSv/yr** (10 CFR 835 / DOE
O 458.1); Fermilab off-site *design goal* **0.1 mSv/yr**; natural background
in Illinois ≈ 3 mSv/yr.

## 2. Zone-by-zone assessment

### Zone A — above the corridor, fence to detector and beyond (0–20 km)
Communities: Fermilab's own north buffer, unincorporated West Chicago,
Wayne, valley of Norton Creek, then St. Charles/Elgin townships.
Plume depth: 100 m at the fence → 108 m at the detector → ~170 m at 20 km
(terrain rises northward faster than the chord drops).
**Surface dose: none measurable.** No shower particle range in rock
approaches these depths; there is no activation of soil or groundwater at
radiologically relevant levels (neutrino interaction density is ~10⁻¹⁷ per
nucleon per year even in the core).

### Zone B — arc "disk" (all azimuths, grazing annulus ~38 km out)
The arcs radiate a 21 µrad-thick disk in the ring plane in *all* directions;
it grazes the surface on an annulus ≈ 38 km away (Elburn/DeKalb farmland to
the west, Joliet area to the south, Schaumburg to the east, McHenry County
to the north). Raw King estimate at the annulus: **0.32 mSv/yr**; with the
±1 mrad mover system: **0.008 mSv/yr** — an order of magnitude below the
design goal, consistent with IMCC's "negligible with movers" conclusion.

### Zone C — RCS straights (on-corridor, wobbled)
Bounding each RCS's decays at its top energy:

| Source | Peak annual dose at exit (wobbled) |
|---|---|
| RCS1 straights (0.31 TeV) | < 0.0001 mSv |
| RCS2 straights (0.75 TeV) | 0.0005 mSv |
| RCS3 straights (1.5 TeV) | 0.005 mSv |
| RCS4 straights (4.5–5 TeV) | **0.64 mSv** → < 0.1 mSv with the same ±0.5 mrad segmentation used on the collider insertions |

These land inside the *same* corridor strips as the collider plume (that is
the point of the alignment) — they add no new exposure geography.

### Zone D — the two collider-straight exit strips (the real issue)
The IP straight cannot be wobbled. Its pencil surfaces at a predictable,
fixed location in each direction; with the baseline horizontal ring:

| Direction | Exit | Location (meridian −88.2230) | Hot-core size | Peak dose, King raw | Mitigated (150 m pencil) | Optimistic (30 m pencil) |
|---|---|---|---|---|---|---|
| North | 37.8 km | lat 42.19 — Fox River valley between Cary and Fox River Grove (McHenry Co.) | ~4 m wide (E–W) × ~7 km (N–S) | 0.70 Sv/yr | **18 mSv/yr** | 3.5 mSv/yr |
| South | 28.3 km | lat 41.59 — Wheatland/Na-Au-Say townships between Oswego and Plainfield | ~4 m wide × ~7 km | 1.25 Sv/yr | **41 mSv/yr** | 8.3 mSv/yr |

The collider's *west* (utility) straight produces the same pair of strips on
meridian −88.2560 (through Crystal Lake to the north, Yorkville farmland to
the south); with full in-straight doglegs plus wobble its peak is
~28 mSv/yr raw-mitigated and it should be designed shorter and more
aggressively segmented than the IP side — or accepted as two further
easement strips.

**Pitch as a siting knob.** Tilting the collider plane slides the exits
(down-to-north pushes the north exit away and pulls the south exit in):

| Ring pitch (down-N) | North exit → dose | South exit → dose |
|---|---|---|
| 0 mrad | 37.8 km (Cary/Fox River Grove) → 18 mSv | 28.3 km (rural Wheatland Twp) → 41 mSv |
| +1 mrad | 43.7 km (near Island Lake) → 12 mSv | 24.0 km (Wolf's Crossing/Oswego) → 67 mSv |
| +2 mrad | 51.0 km (**Moraine Hills State Park**) → 8 mSv | 20.9 km (south Oswego) → 101 mSv |
| +5 mrad | 81.3 km (near WI border, Richmond) → 3 mSv | 14.6 km (Normantown) → 284 mSv |

The south side is the binding constraint (terrain falls toward the Illinois
River valley, pulling exits closer). The baseline recommendation is
**pitch 0**, which keeps *both* strips over the two least-developed
available areas.

## 3. What "controlled strip" actually means

The >1 mSv/yr core at an exit is a strip roughly **as wide as a sidewalk**
(2L/γ + shower width ≈ 4 m E–W; broadened to ~40 m with the ±0.5 mrad
segmentation, which also cuts the peak ~20×) and a few km long where the
pencil skims the surface at 4–6 mrad. Bringing such a strip under control is
a *pipeline-easement-class* land action, not a town-scale exclusion:
fencing/monitoring of a strip comparable to the transmission ROW the beam
already follows. Options, in order of preference:

1. **Land**: acquire easements over the two strips (south strip is active
   farmland; agriculture can continue outside the fenced core — occupancy
   factors for a field are a further ×10–20 in realistic dose).
2. **Pitch tuning** (±1–2 mrad) to drop a strip onto public land (Moraine
   Hills SP, Fox River floodplain) under a DOE–IDNR agreement.
3. **Deeper collider** (200 m, LHC-class): exits recede to ~52/40 km, doses
   fall roughly ∝ 1/depth.
4. **Machine-side**: shorten the true collinear drift (every metre matters:
   dose ∝ pencil length), increase FF divergence credit, staged-luminosity
   early running (dose ∝ N_μ) with in-strip monitoring before full current.

This mirrors the IMCC approach exactly — CERN's own 10 TeV siting answer is
"insertion plumes exit into the Mediterranean and a fenceable mountainside".
Fermilab has no sea; it has flat, surveyable, largely rural land north and
south, which is the equivalent asset.

## 4. Scenario study: a deliberate southern exit at UIUC (≈199 km)

The corridor meridian (−88.222972°) happens to run straight down the state
between Champaign and Urbana — within a few hundred metres of the UIUC Main
Quad — and through the university's **South Farms** agricultural land south of
campus. That invites a designed variant: pitch the ring so the south plume
surfaces there deliberately, on institutional land, instead of at the ~28 km
default. (Load the *"South exit at UIUC South Farms"* preset in the
[beam-geometry tool](../tool/) to explore it live.)

*(Two flavours of this scenario ship as presets in the tool: the one below,
which keeps the north beam underground past Smith Road, and the site default —
straightaway at 185 m ASL, 15.42 mrad — which instead surfaces the north beam
at 41.8746° N, 0.5 km past the Fermilab fence at the southern edge of West
Chicago (Roosevelt Rd corridor), so that it crosses Smith Road a full 300 ft
above the 10 m rooflines; detector hall −14 m just inside the fence. The
default's trade: the first ~2 km after emergence pass low over West Chicago's
Roosevelt Road commercial/residential edge, which inherits the low-overflight
caveat below.)*

**Geometry.** Straightaway at 76 m ASL (148 m below grade at the IP), tilted
**14.87 mrad (0.85°)** up to the north — a hair above the 14 mrad LEP
precedent. Then:

| Quantity | Value |
|---|---|
| South exit | **40.0600° N, −88.2230°** — South Farms/Savoy area, ~5 km south of downtown Urbana |
| Chord through the Earth | 198.9 km, dipping to 0.85 km (bottoming in the Ironton–Galesville sandstone) |
| Exit grazing angle | 16.4 mrad (0.94°) |
| Beam depth under the UIUC Main Quad / downtown Urbana | **88 m / 93 m** |
| North side | detector hall stays at −25 m; north beam emerges at 41.943° N, ~800 m **north** of Smith Rd (houses there stay under a buried beam) |
| Targeting precision | ±0.1 mrad of tilt moves the exit ±1.3 km; ±10 m of terrain/geoid knowledge moves it ±0.6 km — trimming onto a chosen parcel is a survey problem, not a physics one |

**Dose to Urbana-area communities.** Same zone logic as §2. Everywhere north
of the exit — all of Champaign, Urbana, and campus — the beam is 44–900 m
underground and delivers no measurable surface dose. The only nonzero-dose
location is the grazing strip at the exit itself, and 199 km of 1/L² is the
whole story: **King-raw peak ~25 mSv/yr in an ~8 m-wide core** (vs 1,250 at
the 28 km default exit), falling to **~0.5 mSv/yr mitigated** (150 m
effective pencil, ±0.5 mrad segmentation) and **~0.1 mSv/yr** in the
optimistic case — i.e. at the public limit to the design goal *inside the
fence line*, on farm fields. Segmentation trades peak against footprint: the
unsegmented core is ~8 m × ~600 m; fully segmented it stretches toward
~12 km along the meridian at proportionally lower dose. Either way the land
action is a fenced strip through research farmland with a single
institutional owner — about as tractable as exit-zone control gets. South of
the exit the beam climbs at 0.94° over rural Champaign County (Willard
Airport is ~5 km west of the meridian, laterally clear of the pencil).

**Groundwater and the Mahomet Aquifer.** Champaign–Urbana drinks from
aquifer sands in the deep glacial drift, and the beam does pass through that
interval (50–95 m deep) under the metro area on its way up. Neutrinos do not
activate anything they pass through; the only activation pathway is the
hadronic/EM showers from the rare ν interactions in the ground itself.
Rates: in the on-axis core at exit, ~9 × 10⁴ interactions per litre-year in
water. With a generous spallation yield (~10⁻² tritium atoms per GeV of
shower), saturation tritium concentration in stagnant in-core water is
**~0.1 Bq/L — four orders of magnitude below the 740 Bq/L EPA drinking-water
limit** — over a core only ~8 m wide, before any dilution or flow. Under
campus, only the same thin pencil through the aquifer sees these rates at
all. For calibration: Fermilab's NuMI facility manages *actual* tritiated
sump water produced by a primary proton beam striking absorbers — a source
term orders of magnitude harsher than the diffuse ν-interaction wisp here —
within regulatory limits.

**The opportunity.** The exit fluence is still 1.2 × 10¹³ ν/cm²·yr on axis —
**~9 × 10⁷ interactions per tonne-year** at ⟨Eν⟩ ≈ 3.3 TeV. A far detector
on South Farms would sit in the most intense high-energy neutrino beam ever
delivered to a university campus, which reframes the exit strip from a
liability into UIUC's stake in the facility.

**Caveats specific to this scenario.** All §5 caveats apply, plus: the
flat-layer stratigraphy is anchored at Fermilab, and the bedrock section
genuinely changes over 199 km (the drift thickens over the Mahomet Bedrock
Valley; Paleozoic units dip and thicken southward into the Illinois Basin),
so the rock column and the depth-to-aquifer picture near the exit need ISGS
well-log data before anything firmer than order-of-magnitude; and the dose
and activation numbers are the same conservative analytic machinery as the
rest of this note — FLUKA/MARS with a real lattice before public claims.

## 5. Honest caveats

* King's equilibrium model is an **order-of-magnitude** tool, conservative
  by construction but with genuine ±(factor of a few) uncertainty; the
  pencil-length parameter moves results ×5 by itself. All Zone-D numbers
  must be redone with a real lattice + FLUKA/MARS before any siting claim
  is made publicly.
* Terrain from SRTM (~10 m vertical accuracy); geoid vs ellipsoid effects
  along 80 km are cm-to-m level and ignored.
* Muon-beam parameters are IMCC-interim-class and will evolve.
* No credit taken for building shielding, occupancy, or duty factor except
  where stated.

## 6. Fallbacks

* **Tevatron-tunnel RCS1/2** (shallow, misaligned): plumes exit ~9.5 km out
  in six azimuths at ≤ ~0.3 mSv/yr wobbled — acceptable for a staged
  low-energy programme but the aligned deep racetrack is the baseline.
* If Zone-D easements prove impossible, the machine can run at reduced
  μ-current (linear dose scaling) or √s = 3 TeV (E⁴: ×123 lower straight
  dose — every exit strip drops below 0.5 mSv/yr even raw) while land
  actions proceed.
