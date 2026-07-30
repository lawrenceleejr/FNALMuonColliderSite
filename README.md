# FNALMuonColliderSite

Corridor-aligned siting study for a 10 TeV muon collider at Fermilab, in
which the collider's interaction-point straight and the long straights of
the RCS acceleration chain are all placed on the meridian **−88.222972°**
(true-north azimuth), so their muon-decay neutrino plumes stack into a
single N–S corridor. The corridor leaves the site at its northeast corner,
follows the ComEd *Aurora–Wayne* transmission right-of-way, and delivers an
extremely intense TeV neutrino beam to an underground experiment at
**41°55′39.3″ N, 88°13′22.7″ W** (8.71 km from the IP, hall ~109 m deep).

## Contents

| Path | What it is |
|---|---|
| [`docs/DESIGN.md`](docs/DESIGN.md) | Layout, alignment geometry, detector siting, and the explicit list of edits vs published FNAL siting concepts |
| [`docs/SAFETY.md`](docs/SAFETY.md) | Neutrino radiation assessment for downstream communities (zones, doses, exit strips, mitigation) |
| [`tools/corridor_layout.py`](tools/corridor_layout.py) | Parametric generator: fits rings inside the real OSM site boundary, computes plume geometry over real terrain, fluxes, event rates, and King-model doses. Pure stdlib — `python3 tools/corridor_layout.py` |
| [`geo/layout.geojson`](geo/layout.geojson), [`geo/summary.json`](geo/summary.json) | Generated layout and numeric summary |
| [`map/index.html`](map/index.html) | Interactive Leaflet map of the layout (open locally in a browser) |
| [`figs/corridor_profile.svg`](figs/corridor_profile.svg) | Corridor elevation cross-section: terrain vs plume centreline, detector, exit points |
| `data/` | Inputs: FNAL boundary (OSM way 31974155, ODbL), ComEd corridor lines (OSM, ODbL), SRTM elevation profile along the meridian |

## Key numbers (baseline)

* Collider C = 10 km, IP at 41.849142 N on the corridor meridian, 100 m deep.
* RCS3/4 racetrack C = 14.72 km (10.8 % below the 16.5 km site-filler bound —
  the machine cost of straight alignment), RCS1/2 C = 6.28 km stacked over
  the IP straight; RCS planes pitched 2.3 / 4.6 mrad so all plumes converge
  on the detector hall.
* Detector: 6 × 10¹⁵ ν/cm²·yr core fluence, ⟨Eν⟩ ≈ 3.2 TeV,
  ≈ 4 × 10¹⁰ interactions per tonne-year.
* Communities above the corridor: no measurable dose (plume 90–170 m deep).
  The two surface-grazing exit strips of the IP straight (~4 m × ~7 km, at
  38 km N and 28 km S) are the real radiological land-use item — see
  `docs/SAFETY.md` before quoting any number.

**Status: conceptual study.** Dose figures use a deliberately conservative
analytic model (King, arXiv:physics/9908017) and must be confirmed with
FLUKA/MARS and a real lattice before public siting claims.
