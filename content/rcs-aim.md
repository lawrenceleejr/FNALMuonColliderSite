---
title: "Aiming the Last RCS: Depth, Angle, and a Target in Lake Superior"
weight: 8
---

The [co-tilted chain](../uiuc-chain/) locks each ring's straight depth to its
tilt through the UIUC exit constraint — one knob, not two. This page inverts
the parameterisation: **depth and angle are independent inputs**, and the
outputs are where each end of the straight surfaces. That is the frame for
the question the variant provokes — *can the last RCS point somewhere else,
further away, into water, the way IMCC aims its insertion plumes at the
Mediterranean?* Solved by
[`tools/rcs_aim.py`](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/tools/rcs_aim.py)
against a cached GEBCO 2020 bathymetry/topography profile of the corridor
meridian ([rcs_aim.json](../geo/rcs_aim.json)).

**Headline:** the corridor meridian offers exactly **one** verified deep-water
target — **Lake Superior open water, 47.48–48.10° N, 626–695 km out, needing
49.1–54.5 mrad** — and it is a remarkably good one: a 69 km-long aim window
(±2.7 mrad, *10⁵× looser* than the co-tilted chain's 25 µrad planar window),
a ~300 m water column, 30+ km from any shore, and an RCS4 exit dose of
**0.28 mSv/yr raw and unwobbled** — under the public limit before any
mitigation. And the reassignment is nearly free for the far hall:
sending RCS3/4 to the lake keeps **97 % of the co-tilted chain's
plane-crossing ν̄τ gain** while deleting its only new dose problem.

## 1. The two knobs, and why they separate

A chord that re-surfaces a great-circle distance *D* from a straight at depth
*d₀* needs a tilt

<div style="text-align:center;font-family:monospace;margin:.8rem 0">
θ = d₀/D − D/(2 R<sub>E</sub>)
</div>

For any *D* ≫ *d₀* the first term is negligible, so **|θ| ≈ D/(2 R<sub>E</sub>)
— 1.57 mrad per 100 km of range, independent of depth.** The angle *is* the
target. Depth then only decides where the *near* end surfaces,
s ≈ d₀/|θ|. Clean separation:

| knob | controls | does not control |
|---|---|---|
| **angle** | which target the far exit lands in | the near emergence (to <1 %) |
| **depth** | how far out the near end surfaces | the far range (655 → 663 km over 15–440 m of depth) |

<div style="margin:1.2rem 0"><img src="../figs/rcs_aim_knobs.svg" alt="Far range flat in depth; near emergence proportional to depth" style="max-width:100%"></div>

The practical consequence is the opposite of intuition: for a water-aimed
ring you want the straight **shallow**. At 51.8 mrad, a 15 m straight
surfaces its near (south) end 0.3 km out — on DOE land — and climbs through
500 ft by 3.2 km, leaving only **0.59 km of sub-500 ft overflight** beyond
the fence. A 440 m straight surfaces at 8.3 km, *off* site, with 8.6 km of
sub-500 ft band. Deep is worse, and cheaper is better.

## 2. What is actually out there

GEBCO at 0.02° (2.2 km) sampling along the meridian, with Natural Earth 10m
lake polygons naming the bodies. **GEBCO decides land vs water and the
polygon only names it** — a discipline this study needs, because the coarse
polygon gets it wrong in exactly the interesting place:

| candidate | verdict |
|---|---|
| **Lake Superior, 47.480–48.100° N** | **✓ open water.** 68.9 km along the meridian, bed −17 to −117 m MSL → water column up to **300 m** (surface at +183 m). D = 626–695 km, tilt **49.09–54.50 mrad**. |
| Keweenaw Peninsula, 46.90–47.46° N | ✗ **land, +30 to +416 m.** Natural Earth's Lake Superior polygon swallows the whole peninsula; a naive point-in-polygon test puts a "water" exit on a 400 m hill. The beam passes 3.7 km *beneath* it. |
| Kentucky Lake, ~36.96° N (543 km S) | ✗ **unverified.** GEBCO returns 100–128 m with 28 m of relief in ±6.6 km — its 15″ grid is averaging the Tennessee River valley with its bluffs, so it cannot resolve a ~2 km-wide reservoir. Also narrow (≈5 km of meridian) and shallow-recreational: rejected rather than built on. |
| Gulf of Mexico, ≤30.36° N (1280 km S) | ✗ **tilt disqualifier.** Needs **101 mrad** (7× LEP), and the near-shore water column is only 2–26 m. |
| Lake Winnebago, Lake Michigan, Green Bay | ✗ not crossed by this meridian at all. |

<div style="margin:1.2rem 0"><img src="../figs/rcs_aim_chords.svg" alt="The UIUC and Lake Superior chords on the real meridian section" style="max-width:100%"></div>

The chord to Lake Superior dives to a **perigee of 8.6 km** under central
Wisconsin (44.8° N) — an order of magnitude deeper than the UIUC chord's
0.79 km, and far beneath any aquifer, well, or mine on the route. Whatever
the corridor's groundwater story is, this beam has none.

## 3. The recommended aim

Straight at **15 m depth, tilted 51.80 mrad down-to-north**:

| | value |
|---|---|
| Far (north) exit | **47.745° N, 655.4 km** — Lake Superior open water, ~290 m column, grazing 51.1 mrad |
| Aim window | 49.09–54.50 mrad, i.e. **±2.7 mrad** — vs ±0.012 mrad for the planar UIUC aim |
| Exit dose, RCS4, King-raw unwobbled | **0.277 mSv/yr** (RCS3: 0.002) — *below* the 1 mSv/yr public limit with no mitigation at all |
| Near (south) exit | 41.842° N, **0.3 km — on DOE land**, in the emergence zone the collider's own streak already requires |
| Near-side civic envelope | 121 m above grade at the south fence; 500 ft by 3.2 km → **0.59 km** of sub-500 ft overflight over Aurora NE |
| Chord perigee | 8.6 km beneath 44.8° N |
| Ring geometry | planar is **impossible** (a 51.8 mrad plane lifts the shallow arc apex 126 m — above grade); the [bent ring](../uiuc-chain/) is mandatory: 864 T·m at 5 TeV = 108 m of 8 T per end, 2 ends = **1.5 % of the ring** |

Two things to weigh honestly. The 0.59 km of sub-500 ft overflight is over
Aurora NE/Eola residential rather than the ComEd right-of-way the baseline
climbs — worse land, but *smaller* than the 1.9 km of 384–500 ft Wayne
overflight the [envelope note](../envelope/) already accepts, and the steep
51.8 mrad climb is what makes it short. And the near-side emergence streak
carries Sv-class raw dose in a decimetres-wide line 0.3 km south of the
straight: on-site, fenced, occupancy-controlled — the same object as the
co-tilted chain's north streaks, not a new species.

## 4. What UIUC gives up — and why it barely matters

RCS3 and RCS4 share one tunnel, so aiming "the last RCS" at the lake removes
both from the UIUC chirp. Same machinery as the
[co-tilted chain](../uiuc-chain/) §2, τ-threshold included:

| at the UIUC far hall | store only | store + full chain | **store + RCS1/2 only** (RCS3/4 → lake) |
|---|---|---|---|
| on-axis flux (integral) | 1.00× | 2.12× | 1.20× |
| ν̄τ CC per tonne-year | 0.137 | 0.297 | **0.200** |
| ν̄μ CC per tonne-year | 4.35 × 10⁵ | 7.90 × 10⁵ | 4.57 × 10⁵ |
| plane-crossing ν̄τ, E > 5 GeV | 1.04 × 10¹⁴ | 9.66 × 10¹⁴ | **9.37 × 10¹⁴** |

Reassigning RCS3/4 to the lake **retains 96.6 % of the full chain's
plane-crossing ν̄τ gain** and 39 % of its on-axis ν̄τ CC gain, while giving up
82 % of its flux gain. The reason is the physics the timing study already
exposed: oscillated ν̄τ come from the **low-energy** turns, where (L/E)² is
largest — RCS1/2, not RCS3/4. So the ring you most want to send away is the
one UIUC needs least:

* it is the **only new dose problem** the co-tilted chain created at UIUC
  (RCS4: 3.04 mSv/yr raw in an ~11 m strip) — **deleted**, replaced by
  0.28 mSv/yr into 300 m of water;
* it is the ring whose planar tilt window was 25 µrad — **replaced** by
  ±2.7 mrad;
* and it costs the far hall's τ programme almost nothing.

Meanwhile the lake exit is not merely a dump. A 656 km baseline with a
63 GeV → 5 TeV turn-tagged chirp reaches L/E ≈ 10 km/GeV and P up to ~10⁻²
— the corridor's largest oscillation phase anywhere. Whether a detector
belongs in 300 m of Lake Superior water is far outside this study's scope,
but the flux is there and the exit is, by construction, somewhere nobody
stands.

*Sequel:* [azimuth freedom](../rcs-azimuth/) drops the assumption that the
aim must lie on the ComEd meridian at all. Rotating the bearing ENE puts the
near beam across 4.7 km of Fermilab land instead of 2.6 km, which removes even
the 0.59 km of overflight above — and opens **Lake Huron at 630 km** (179 m of
water, azimuth-free) alongside a strictly better Lake Superior aim at 8°
(319 m of column, a 248 km submerged run).

## 5. Open items

* **The 0.59 km overflight** needs the envelope treatment the north side
  already gets, over a residential rather than utility corridor.
* **Ramped vertical achromats** at 864 T·m in a rapid-cycling ring is the
  real engineering ask, unverified here (as is the 257 T·m version the
  co-tilted chain needs).
* **International and lake-jurisdiction questions** the study has not
  touched: the aim window's north end approaches the Canadian border
  (48.10° N is ~40 km south of it), and the exit sits in US waters of a
  Great Lake under an entirely different regulatory regime than an Illinois
  farm strip — plausibly *easier* (no occupancy) or *harder* (Great Lakes
  compact, treaty, IJC) than the UIUC strip. That determination is a legal
  question, not a physics one, and this page does not pretend to make it.
* **GEBCO at 2.2 km sampling** is adequate to find a 69 km target and to
  reject a 2 km reservoir, but the final aim point wants NOAA Lake Superior
  bathymetry and a shoreline check on the real exit ellipse.
* Arc-disk azimuthal doses, still owed by the baseline itself, get worse
  with every additional tilted ring.
