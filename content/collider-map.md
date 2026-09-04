---
title: "The Collider's Compass: Landmarks Around the Default Aim"
weight: 10
---

The [azimuth-freedom map](../rcs-azimuth/) was drawn for a ring whose tilt is
a free knob. The collider's is not. **15.40 mrad up-to-north from a 35 m-deep
straight** is fixed by the UIUC exit, and a fixed tilt surfaces at a fixed
range — **198.4 km in every bearing** (2 R<sub>E</sub> θ = 196.2 km, plus
2.2 km from the depth term). So the collider's version of that plot is a
*circle*, and the questions are what sits on it, what sits under the two
beams, and what else the compass holds within reach of *some* tilt.
[`tools/collider_map.py`](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/tools/collider_map.py)
projects the real geography into the same (bearing, log-range) frame: Natural
Earth lake footprints, DOE and university partners, mines and quarries,
public land, airports, and the population centres the beam must not graze
(machine-readable: [collider_map.json](../geo/collider_map.json), including
where the circle lands for every degree of bearing).

**Headline:** the 198 km circle lies **entirely inside θ<sub>free</sub>** —
the collider is corridor-bound at *every* bearing, with 7.5–11.3 km of
sub-500 ft overflight beyond the fence wherever it points. Its meridian is
therefore chosen by the **near** side, not the far side: the ComEd
right-of-way is why the far exit is at UIUC rather than at any other point
of the circle. And the circle is well populated — **Purdue at 192 km (SE)**,
**Kettle Moraine State Forest at 195 km due north** (the mirror-image aim),
and **45° of open Lake Michigan** between bearings 11° and 55°.

## 1. The map

<div style="margin:1.2rem 0"><img src="../figs/collider_map.svg" alt="Polar landmark map around the collider: bearing around, log range radial, with the 198 km exit circle, theta_free, lake footprints and categorised landmarks" style="max-width:100%"></div>

Radius is log(range) so that 10 km neighbours and 1300 km mines share one
frame; every range ring is also labelled with the tilt that reaches it,
θ ≈ D / 2 R<sub>E</sub>. Lakes are real footprints projected into the frame,
not markers. The green curve is the site's own θ<sub>free</sub>(bearing) from
the [azimuth study](../rcs-azimuth/), drawn at its equivalent range
2 R<sub>E</sub> θ<sub>free</sub> = 430–1220 km. Sanity checks against known
baselines: the map places Soudan at 736 km (MINOS: 735 km, and NuMI's 58 mrad
downward pitch is exactly 736 km / 2 R<sub>E</sub>), Ash River at 811 km
(NOvA: 810 km), and SURF at 1289 km (DUNE: ~1300 km).

## 2. What is on the circle

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

Three readings:

1. **The band never closes.** At 15.4 mrad the near beam needs
   (d₀ + 152 m)/θ = 12.2 km to clear 500 ft, and the fence is 0.9–4.7 km
   away depending on bearing, so 7.5–11.3 km of the climb is over someone
   else's land *whatever the bearing*. The best case (far exit east, near
   beam over the site's 4.7 km west axis) saves 1.8 km against the baseline;
   the worst (far exit west, near beam east over 0.9 km of DOE land) costs
   2 km. The corridor is unavoidable; the design question is only *whose*
   land it crosses — and the ComEd right-of-way is the only bearing that
   offers a utility easement rather than houses.
2. **Same tilt, different water.** Rotating the ring 11–55° east puts the
   far exit into Lake Michigan at the *same* 15.4 mrad — an exit where
   nobody stands, the [RCS aiming](../rcs-aim/) ideal. The price is on the
   near side: at az 24° the beam climbs SSW over Aurora and Sugar Grove for
   9.1 km, the densest residential sector on the compass, against the ComEd
   easement's 9.3 km. That is the trade the baseline made, and this map is
   the first place the study has shown both sides of it.
3. **The circle is a university circle.** UIUC by design, Purdue at 192 km
   by accident of geography, Notre Dame (166 km) and UW–Madison (168 km) a
   few mrad inside, Michigan State and Indiana at 25 mrad. If the far-hall
   concept ever needs a second host or a comparison site, the candidates are
   already on the plot.

## 3. Under the two beams

The meridian itself carries landmarks the [envelope](../envelope/) and
[aquifer](../aquifer/) pages already audit; the map puts them on one bearing
line.

**North, airborne past 2.2 km.** Smith Road at 10.1 km and Wayne at 12 km
are the houses the near beam climbs over at 380–500 ft — the corridor marks
on the plot. Beyond them the bearing passes Moraine Hills State Park (51 km)
and Chain O'Lakes State Park (68 km) with the beam already 0.9–1.4 km up,
then Kettle Moraine SF at 195 km, the Menominee Range iron mines at Iron
Mountain (442 km, 34.7 mrad), the Marquette Range (515 km, 40.4 mrad), and
Lake Superior (626–695 km, 49–54 mrad — the [meridian water
target](../rcs-aim/)). The two tick marks at 38 km and 28 km are the
[safety page's](../safety/) untilted "Zone D" default exits (Cary/Fox River
Grove and Wheatland Township) — where the chord would surface if the ring
were level.

**South, underground to 198 km.** The chord passes under Midewin National
Tallgrass Prairie at 54 km (638 m deep), then reaches its 790 m perigee at
98 km — 10 km west of the **Herscher** gas-storage field (88 km, chord at
783 m) — and climbs past the **Manlove** field at 161 km (chord at 484 m,
2 km west of the field centre) before surfacing at the South Farms. Both
fields store gas in the Mt. Simon sandstone, the same brackish horizon the
aquifer page puts the perigee in. Whether the chord intersects an actual
storage reservoir is a question for the field records that page already
lists as owed; the map's contribution is to show that the two fields the
region's gas-storage operators use are the two nearest subsurface facilities
to the beam.

## 4. The rest of the compass

Everything else on the plot is out of the collider's reach at 15.4 mrad and
in reach of the tilts the RCS studies contemplate:

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

The pattern the [azimuth study](../rcs-azimuth/) found reappears here from
the landmark side: everything worth aiming at from a *far-hall* point of view
— deep mines, big water, big public land — sits at 440–800 km and 35–65 mrad,
which is exactly the θ<sub>free</sub> class where the near beam clears the
fence in the air. The collider at 15.4 mrad cannot have any of it; the RCS
straights can.

## 5. Open items

* Landmark coordinates are gazetteer-grade (~0.01°), fine for this frame but
  not for an aim point; the circle landings and water arc use spherical
  geodesy and the Natural Earth 10m polygons.
* "Public land" here is national forests, state parks and Midewin; county
  forest preserves, state wildlife areas and the Fox River corridor are not
  drawn. The near-side band's land use on any bearing other than north remains
  the unverified item the azimuth study flagged.
* The gas-storage reservoirs' depths and lateral extents relative to the
  chord are not modelled — only the fields' surface positions are plotted.
