# FNALMuonColliderSite

Corridor-aligned siting study for a 10 TeV muon collider at Fermilab: an
11 km ring whose interaction-point straight — and the long straights of the
RCS acceleration chain — all lie on the meridian **−88.222972°**
(true-north azimuth), so their muon-decay neutrino plumes stack into a
single N–S corridor. In the civic-envelope baseline the straight is tilted
15.40 mrad: the north beam feeds an on-site experiment, leaves the ground
**inside the Fermilab fence**, and climbs the ComEd *Aurora–Wayne*
right-of-way into navigable airspace; the south beam dives under Aurora and
resurfaces 198 km away on the **UIUC South Farms** (the second IP's south
beam exits onto UIUC's Willard Airport). Wherever a beam is within 500 ft
of the surface, the land is DOE, utility-easement, or University of
Illinois property.

## Contents

| Path | What it is |
|---|---|
| [`content/design.md`](content/design.md) | Layout, alignment geometry, detector siting, and the explicit list of edits vs published FNAL siting concepts |
| [`content/safety.md`](content/safety.md) | Neutrino radiation assessment for downstream communities (zones, doses, exit strips, mitigation) |
| [`content/envelope.md`](content/envelope.md) | Airspace & land-ownership legal geometry: the ±500 ft civic envelope, Causby, 14 CFR 91.119, subsurface easements, segment-by-segment walk |
| [`content/beamsize.md`](content/beamsize.md) | Neutrino spot size from the collider optics: the 1/γ profile, containment radii, why β* drops out, exit-strip footprints, detector sizing |
| [`content/aquifer.md`](content/aquifer.md) | Aquifer radiation audit: activation of every water-bearing unit vs EPA drinking-water limits |
| [`paper/corridor.tex`](paper/corridor.tex) | Short PRL-format paper draft of the concept |
| [`tools/corridor_layout.py`](tools/corridor_layout.py) | Parametric generator: fits rings inside the real OSM site boundary, computes plume geometry over real terrain, fluxes, event rates, and King-model doses. Pure stdlib — `python3 tools/corridor_layout.py` |
| [`static/geo/layout.geojson`](static/geo/layout.geojson), [`static/geo/summary.json`](static/geo/summary.json) | Generated layout and numeric summary |
| [`static/tool/index.html`](static/tool/index.html) | Interactive beam-geometry tool: terrain + stratigraphy ray trace, constraint solver, emergence map, interaction-length accounting, Geant4 export (rock-column GDML + macro) for [G4TargetPractice](https://github.com/lawrenceleejr/G4TargetPractice) |
| [`static/map/index.html`](static/map/index.html) | Interactive Leaflet map of the layout (open locally in a browser) |
| [`static/figs/corridor_profile.svg`](static/figs/corridor_profile.svg) | Corridor elevation cross-section: terrain vs plume centreline, detector, exit points |
| `data/` | Inputs: FNAL boundary (OSM way 31974155, ODbL), ComEd corridor lines (OSM, ODbL), SRTM elevation profile along the meridian |

## Key numbers (civic-envelope baseline)

* Collider C = 11.0 km (R_arc = 1 528 m, 2 × 700 m straights), IP at
  41.84428 N on the corridor meridian; straightaway 191 m ASL (34 m below
  grade), ring plane tilted 15.40 mrad (one LEP). Second IP on the west
  straight, meridian −88.2598°.
* North beam: emerges on site at 41.8640° N, crosses the fence 27 ft up,
  crosses Smith Rd 380 ft above grade (350 ft above a 10 m roofline), and
  enters navigable airspace 9.3 km out.
* South beam: 74–152 m under Aurora, perigee 789 m in the brackish Mt. Simon, exits at 40.067° N on
  UIUC South Farms; Urbana passage 66–119 m deep. IP2's south exit lands on
  UIUC's Willard Airport with a +0.13 mrad trim.
* RCS3/4 racetrack C = 14.72 km (10.8 % below the 16.5 km site-filler bound —
  the machine cost of straight alignment), RCS1/2 C = 6.28 km stacked over
  the IP straight; RCS planes pitched 2.16 / 4.32 mrad so all plumes
  converge on the detector hall.
* Detector: ~5 × 10¹⁵ ν/cm²·yr core fluence, ⟨Eν⟩ ≈ 3.2 TeV,
  ≈ 1.2 × 10¹¹ interactions per tonne-year (deep-hall reference; CSMS σ, both species).
* Aquifers: every potable-aquifer crossing beyond the fence ≥ ~45× below the EPA tritium
  MCL in the stagnant worst case; whole-chord production ≈ 0.06 Ci/yr —
  see `content/aquifer.md`.
* Airspace/land: the ±500 ft civic envelope holds everywhere except four
  enumerated deep/high crossings — see `content/envelope.md` before quoting
  any number.

**Status: conceptual study.** Dose figures use a deliberately conservative
analytic model (King, arXiv:physics/9908017) and must be confirmed with
FLUKA/MARS and a real lattice before public siting claims.
