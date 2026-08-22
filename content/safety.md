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
strips** tens of kilometres away — now tens of metres wide, per the measured
plume of [arXiv:2608.02718](https://arxiv.org/abs/2608.02718) — which must
be brought under easement/control, exactly as the international design
(IMCC) solves the same problem by aiming its insertion plumes at the
Mediterranean and a fenced mountainside. Everything else (arcs, RCS straights, utility straight) is
driven below or near the 0.1 mSv/yr design goal by standard mitigation.

## 1. Physics of the hazard

Muons decay in flight; each straight section emits a neutrino plume along
its axis. For a FODO straight (the RCS straights, the utility straight) that
plume is a true "pencil" of half-angle 1/γ = **21 µrad** at 5 TeV; for the
IP straight it is not: MINT ([arXiv:2608.02718](https://arxiv.org/abs/2608.02718)),
decaying muons along the real IMCC interaction-region lattice, finds the
final-focus optics smear it to σ<sub>θ</sub> ≈ **0.15 mrad ≈ 7/γ**, diluting
the on-axis density by 2γ²σ<sub>θ</sub>² ≈ **101** before any deliberate
mitigation (see [the beam-size study](../beamsize/)). Neutrinos themselves are
essentially harmless to people they pass through; the dose arises where the
pencil **grazes the ground surface**, because neutrino interactions in the
last tens of metres of soil produce hadronic/EM showers that reach people at
the surface. Between the ring and that exit point the pencil is deep
underground and the showers it makes are absorbed by rock far below anyone's
feet: **surface dose above a buried plume is negligible; only the exit
regions matter.** Dose scales like E⁴ × (straight length) for straights and
E³ for arcs — which is why a 10 TeV machine's straights dominate its entire
radiological design.

The transverse size of the IP-straight plume, and hence the width of its
exit strips, is therefore set by the lattice, not the boost: the 50 %-flux
core on the ground at the far southern exit is **72 m × 4.6 km** (it was
8.3 m × 534 m in this study's superseded pencil model). Wider strip, but the
same physics divides the peak dose by ~100 — the mitigation this page used
to have to assume is now measured into the beam.

Model: B.J. King, *Potential Hazards from Neutrino Radiation at Muon
Colliders*, arXiv:physics/9908017 — Eq. 10 (straights), Eq. 7 (arc disk),
"equilibrium approximation" (conservative: assumes full shower buildup).
Beam: 1.8 × 10¹² μ/bunch/sign, 5 Hz, 1.2 × 10⁷ s/yr, 90 % transmission,
85.4 % of stored muons decaying within the 0.2 s store (the rest are dumped)
→ 0.83 × 10²⁰ decays/yr/sign in the collider.

Corrections applied (each labelled in the output):

* **IP-plume divergence dilution ÷101** (2γ²σ<sub>θ</sub>², σ<sub>θ</sub> =
  0.15 mrad): the dominant factor, and no longer an assumption — it is the
  measured property of the IMCC lattice (arXiv:2608.02718). It replaces this
  page's former "effective pencil length" parameter (150 m baseline / 30 m
  optimistic), which discarded 79 % of the decays to get a similar credit;
  every decay is now counted. The old parameter guessed the right length —
  MINT's longest drift is ~150 m — and the wrong divergence.
* **σν propagator suppression** vs the linear low-energy extrapolation:
  **×0.92**, flux-weighted over the two species at the *fluence-mean*
  energies ⟨Eν⟩ ≈ 1.5–1.8 TeV (the prism washout halves the mean energy at
  any point in the core; at the old on-axis 3–3.5 TeV the factor was 0.84).
* **Vertical segmentation** of insertion straights (deliberate ±0.5 mrad
  vertical doglegs): now an *optional further* peak dilution of ~×3.5–4 on
  top of the divergence, no longer load-bearing.
* **Mover ("wobbling") system**, ±1 mrad, per IMCC interim report
  (arXiv:2407.12450) — applied to arcs, utility and RCS straights, *not* to
  the IP straight (it would destroy luminosity). The utility straight, being
  FODO, really is a 1/γ pencil and genuinely needs its doglegs and movers;
  the RCS wobble trades against the chirp physics of
  [the timing study](../timing/) §2.
* **The pencil fraction f** ([beam-size study](../beamsize/) §4): a dedicated
  dispersion-free drift would put a fraction f of the straight's decays back
  into a true 1/γ pencil, restoring detector rate and exit-strip dose
  *together*. All doses below quote f = 0 (today's lattice) with the
  f = 0.49 scenario shown where it bites.

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
to the north). Raw King estimate at the annulus: **0.26 mSv/yr**; with the
±1 mrad mover system: **0.007 mSv/yr** — an order of magnitude below the
design goal, consistent with IMCC's "negligible with movers" conclusion.

### Zone C — RCS straights (on-corridor, wobbled)
Bounding each RCS's decays at its top energy:

| Source | Peak annual dose at exit (wobbled) |
|---|---|
| RCS1 straights (0.31 TeV) | < 0.0001 mSv |
| RCS2 straights (0.75 TeV) | 0.0005 mSv |
| RCS3 straights (1.5 TeV) | 0.005 mSv |
| RCS4 straights (4.5 TeV) | **0.83 mSv** → ~0.1 mSv with a further ±0.5 mrad segmentation |

These land inside the *same* corridor strips as the collider plume (that is
the point of the alignment) — they add no new exposure geography. Note the
trade documented in [the timing study](../timing/) §2: the RCS wobble that
buys these numbers sweeps the chirp pencils across the deep hall, so it
should be *programmed*, not random, if the chain's neutrinos are to be used.

### Zone D — the two collider-straight exit strips (the real issue)
The IP straight cannot be wobbled. Its pencil surfaces at a predictable,
fixed location in each direction; with the baseline horizontal ring:

| Direction | Exit | Location (meridian −88.2230) | Hot-core size | King raw (all 700 m as a pencil) | Smeared plume (f = 0, no other knobs) | + ±0.5 mrad segmentation | Dedicated drift (f = 0.49) |
|---|---|---|---|---|---|---|---|
| North | 38.2 km | lat 42.19 — Fox River valley between Cary and Fox River Grove (McHenry Co.) | ~16 m wide (E–W) × ~2.6 km (N–S) | 2.5 Sv/yr | **23 mSv/yr** | 6.6 mSv/yr | 1.12 Sv/yr |
| South | 27.9 km | lat 41.59 — Wheatland/Na-Au-Say townships between Oswego and Plainfield | ~12 m wide × ~2.8 km | 4.7 Sv/yr | **43 mSv/yr** | 13 mSv/yr | 2.1 Sv/yr |

The collider's *west* (utility) straight produces the same pair of strips on
meridian −88.2560 (through Crystal Lake to the north, Yorkville farmland to
the south). Being FODO, it is a true 1/γ pencil with no free divergence
credit: with full in-straight doglegs plus wobble its peak is ~39 mSv/yr,
and it should be designed shorter and more aggressively segmented than the
IP side — or accepted as two further easement strips.

The last column is the price of the [beam-size study](../beamsize/) §4
physics scenario, stated without flinching: a dedicated 1/γ drift aimed at
the *default* exits is a Sv-class strip and is simply not sitable there.
A sharp pencil is only compatible with an exit on controlled land — which
is what §4 below is for.

**Pitch as a siting knob.** Tilting the collider plane slides the exits
(down-to-north pushes the north exit away and pulls the south exit in):

| Ring pitch (down-N) | North exit → dose (f = 0) | South exit → dose (f = 0) |
|---|---|---|
| 0 mrad | 38.2 km (Cary/Fox River Grove) → 23 mSv | 27.9 km (rural Wheatland Twp) → 43 mSv |
| +1 mrad | 43.9 km (near Island Lake) → 17 mSv | 23.7 km (Wolf's Crossing/Oswego) → 59 mSv |
| +2 mrad | 51.2 km (**Moraine Hills State Park**) → 13 mSv | 20.7 km (south Oswego) → 78 mSv |
| +5 mrad | 80.8 km (near WI border, Richmond) → 5 mSv | 14.5 km (Normantown) → 158 mSv |

The south side is the binding constraint (terrain falls toward the Illinois
River valley, pulling exits closer). The baseline recommendation is
**pitch 0**, which keeps *both* strips over the two least-developed
available areas.

## 3. What "controlled strip" actually means

The >1 mSv/yr core at an exit is a strip roughly **as wide as a country
lane** (2 × 1.2 σ<sub>θ</sub>L + shower width ≈ 12–16 m E–W at the default
exits; ±0.5 mrad segmentation broadens it further and cuts the peak another
~3.5×) and 2–3 km long where the plume skims the surface at 4–6 mrad. Bringing such a strip under control is
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
which keeps the north beam underground past Smith Road, and the study
baseline — straightaway at 191 m ASL, 15.40 mrad — which instead surfaces
the north beam at 41.8640° N, **still 655 m inside the Fermilab fence**, so
that it crosses the boundary 27 ft up, climbs the ComEd right-of-way, and
passes Smith Road 380 ft above grade — 350 ft above the 10 m
rooflines; detector ~15 m deep on
site. The baseline's civic-envelope walk — including its one overflight
caveat, 1.9 km of Wayne at 384–500 ft AGL — is in the
[envelope note](../envelope/).)*

**Geometry.** Straightaway at 76 m ASL (148 m below grade at the IP), tilted
**14.87 mrad (0.85°)** up to the north — a hair above the 14 mrad LEP
precedent. Then:

| Quantity | Value |
|---|---|
| South exit | **40.0600° N, −88.2230°** — South Farms/Savoy area, ~5 km south of downtown Urbana |
| Chord through the Earth | 198.9 km, dipping to 0.85 km (bottoming in the brackish upper Mt. Simon sandstone — the region's gas-storage horizon, sealed beneath the Eau Claire aquitard) |
| Exit grazing angle | 16.4 mrad (0.94°) |
| Beam depth under the UIUC Main Quad / downtown Urbana | **88 m / 93 m** |
| North side | detector hall stays at −25 m; north beam emerges at 41.943° N, ~800 m **north** of Smith Rd (houses there stay under a buried beam) |
| Targeting precision | ±0.1 mrad of tilt moves the exit ±1.3 km; ±10 m of terrain/geoid knowledge moves it ±0.6 km — trimming onto a chosen parcel is a survey problem, not a physics one |

**Dose to Urbana-area communities.** Same zone logic as §2. Everywhere north
of the exit — all of Champaign, Urbana, and campus — the beam is 44–900 m
underground and delivers no measurable surface dose. The only nonzero-dose
location is the grazing strip at the exit itself, and 199 km of 1/L² is the
whole story: the all-700-m King-raw figure of ~92 mSv/yr lands at
**0.84 mSv/yr from the measured plume alone** (σθ = 0.15 mrad, f = 0 — no
segmentation, no discarded decays, the ×0.92 σν factor), in a core
**~74 m wide × ~4.5 km long**; a further ±0.5 mrad segmentation takes it to
**~0.23 mSv/yr**. That is under the 1 mSv/yr public limit on the lattice as
it exists and within ~2× of the 0.1 mSv/yr design goal with one deliberate
knob — *inside the fence line*, on farm fields. The counterweight is the
[beam-size study](../beamsize/) §4 physics scenario: a dedicated 1/γ drift
(f = 0.49) aimed here concentrates **~41 mSv/yr into a ~10 m × 640 m core**
— viable only as a fenced, monitored strip with occupancy control, and the
clearest statement in this study that detector rate and exit dose are now
one knob. Either way the land action is a fenced strip through research
farmland with a single institutional owner — about as tractable as
exit-zone control gets. South of
the exit the beam climbs at 0.94° over rural Champaign County (Willard
Airport is ~5 km west of the meridian, laterally clear of the pencil).

**Groundwater and the Mahomet Aquifer.** Champaign–Urbana drinks from
wells in the Mahomet aquifer system — sands in the deep glacial drift.
Where the corridor crosses the buried Mahomet Bedrock Valley's main trunk
(~40.25–40.5° N) the beam is 350–550 m *beneath* the valley fill; it meets
basal-drift sands of the same family only on the final climb-out, at
~66–100 m depth under northwest Urbana (the same picture, unit by unit, as
the [aquifer audit](../aquifer/) §5). Neutrinos do not
activate anything they pass through; the only activation pathway is the
hadronic/EM showers from the rare ν interactions in the ground itself.
Rates: in the on-axis core at exit, ~8 × 10² interactions per litre-year in
water (the divergence dilution ÷101 applies here too). With a generous
spallation yield (~10⁻² tritium atoms per GeV of shower), saturation tritium
concentration in stagnant in-core water is **~0.001 Bq/L — nearly six orders
of magnitude below the 740 Bq/L EPA drinking-water limit** — before any
dilution or flow (with the f = 0.49 dedicated drift the old ~0.1 Bq/L core
figure returns, still four orders under the limit). Under
campus, only the same thin pencil through the aquifer sees these rates at
all. For calibration: Fermilab's NuMI facility manages *actual* tritiated
sump water produced by a primary proton beam striking absorbers — a source
term orders of magnitude harsher than the diffuse ν-interaction wisp here —
within regulatory limits.

**The opportunity.** The exit fluence is ~9.6 × 10¹⁰ ν/cm²·yr per species
on axis on the current lattice (f = 0) — **~1.3 × 10⁶ interactions per
tonne-year** at ⟨Eν⟩ ≈ 1.6 TeV — rising to 4.7 × 10¹² and ~1.1 × 10⁸/t-yr
with the dedicated drift. Either is the most intense high-energy neutrino
beam ever delivered to a university campus, which reframes the exit strip
from a liability into UIUC's stake in the facility. The 197 km baseline also
makes it the line's only viable ντ-appearance site: **~3.4 × 10¹¹
oscillation-made ντ cross the exit plane per year**, giving ~0.1–0.2
(f = 0) to ~2.9 (dedicated drift) **ντ CC per tonne-year** against a νμ CC
background 23 000× smaller relative to any on-site location — see the
[beam-size study](../beamsize/) §5, including the rock-produced ντ
background that page now carries.

**Caveats specific to this scenario.** All §5 caveats apply, plus: the
stratigraphy is an interpolated dipping model (±75–150 m on mid-corridor
unit tops — see the [aquifer audit](../aquifer/) §1), and the bedrock section
genuinely changes over 199 km (the drift thickens over the Mahomet Bedrock
Valley; Paleozoic units dip and thicken southward into the Illinois Basin),
so the rock column and the depth-to-aquifer picture near the exit need ISGS
well-log data before anything firmer than order-of-magnitude; and the dose
and activation numbers are the same conservative analytic machinery as the
rest of this note — FLUKA/MARS with a real lattice before public claims.

## 5. Honest caveats

* King's equilibrium model is an **order-of-magnitude** tool, conservative
  by construction but with genuine ±(factor of a few) uncertainty. The
  former dominant unknown — the effective pencil length — is replaced by
  σθ = 0.15 mrad from arXiv:2608.02718's simulated lattice; what remains
  open is the pencil fraction f of a *corridor* lattice, which moves Zone-D
  results ×50 between its extremes and is a design choice, not an unknown.
  All Zone-D numbers must be redone with a real lattice + FLUKA/MARS before
  any siting claim is made publicly.
* Terrain from SRTM (~10 m vertical accuracy); geoid vs ellipsoid effects
  along 80 km are cm-to-m level and ignored.
* Muon-beam parameters are IMCC-interim-class and will evolve.
* No credit taken for building shielding, occupancy, or duty factor except
  where stated.

## 6. Fallbacks

* **Tevatron-tunnel RCS1/2** (shallow, misaligned): plumes exit ~9.5 km out
  in six azimuths at ≤ ~0.5 mSv/yr wobbled — acceptable for a staged
  low-energy programme but the aligned deep racetrack is the baseline.
* If Zone-D easements prove impossible, the machine can run at reduced
  μ-current (linear dose scaling) or √s = 3 TeV (E⁴: ×123 lower straight
  dose — every exit strip drops below 0.5 mSv/yr even raw) while land
  actions proceed.
