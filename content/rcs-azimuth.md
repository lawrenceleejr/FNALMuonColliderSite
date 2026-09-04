---
title: "Azimuth Freedom: What a Steep Beam Unlocks"
weight: 9
---

The corridor meridian **−88.222972°** was never chosen for physics. It was
chosen because the collider's shallow 15.40 mrad north beam needs ~9.5 km of
institutional land to climb through the private-airspace band, and the ComEd
Aurora–Wayne right-of-way is what happens to run north from the site. That
is a constraint on the *angle*, and it dissolves at steep angles:

<div style="text-align:center;font-family:monospace;margin:.8rem 0">
emergence s = d₀/θ &nbsp;·&nbsp; clears 500 ft at (d₀+152 m)/θ &nbsp;·&nbsp; off-site band ∝ 1/θ
</div>

| tilt | far target | emerges | clears 500 ft | off-site sub-500 ft |
|---|---|---|---|---|
| 15.4 mrad (UIUC) | 196 km | 0.97 km | 10.9 km | **8.2 km** |
| 30 mrad | 382 km | 0.50 km | 5.6 km | 2.9 km |
| 49.5 mrad (Lake Huron) | 630 km | 0.30 km | 3.4 km | 0.7 km |
| 63.6 mrad | 810 km | 0.24 km | 2.6 km | **0** |

Past **θ_free = (d₀ + 152 m)/s_fence** the near beam is in navigable
airspace *before it leaves the fence*: no corridor, no easement, no
overflight — and therefore **no reason to stay on the ComEd meridian**. Each
ring's bearing becomes an independent design variable. Mapped by
[`tools/rcs_azimuth.py`](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/tools/rcs_azimuth.py)
([rcs_azimuth.json](../geo/rcs_azimuth.json)).

**Headline:** the freedom's real value is not reaching better water — it is
choosing **how much of Fermilab's own land the near beam climbs over**. The
IP sits on the site's *eastern* edge (4.7 km of DOE land to the WSW, only
0.9 km east), so aiming the far beam **ENE** puts the near beam across the
site's long axis and drops θ_free from 63.6 to **35 mrad**. That single
rotation turns a corridor-bound aim into a self-contained one — and the best
target it opens is **Lake Huron at azimuth 50°, 630 km, 49.5 mrad**: 179 m
of water, a 224 km submerged run, an RCS4 exit dose of **0.30 mSv/yr raw**,
and **zero off-site overflight**.

## 1. The opportunity map

<div style="margin:1.2rem 0"><img src="../figs/rcs_azimuth_map.svg" alt="Polar map of required tilt vs bearing against theta_free" style="max-width:100%"></div>

Every major lake transformed into (azimuth, required-tilt) space, GEBCO-verified
at the aim point. θ_free(bearing) is the site's own curve — outside it the
azimuth is free.

| target | az | D | tilt | water column | submerged run | RCS4 dose* | near-side fence | verdict |
|---|---|---|---|---|---|---|---|---|
| **Lake Huron** | 50° | 630 km | 49.5 mrad | **179 m** | 224 km | 0.30 mSv/yr | 4749 m WSW | **free** |
| Lake Ontario | 74° | 856 km | 67.2 mrad | 176 m | 256 km | 0.16 mSv/yr | 4856 m | free, steeper achromat |
| Lake Erie | 82° | 628 km | 49.3 mrad | 44 m | 192 km | 0.30 mSv/yr | 4721 m | free, but shallow |
| Lake Superior | 8° | 654 km | 51.4 mrad | **319 m** | 248 km | 0.28 mSv/yr | 2678 m S | corridor-bound (0.58 km) |
| Lake Michigan | 24° | 310 km | 24.4 mrad | 240 m | 256 km | 1.24 mSv/yr | 3057 m | corridor-bound (3.8 km) |

\* King-raw, unwobbled, RCS4 bounded at top energy.

Four readings:

1. **Lake Huron is the best overall.** It is the only candidate that is
   simultaneously deep (179 m), far enough for a low dose (630 km), and at a
   bearing whose near side crosses 4.7 km of Fermilab land — so the beam is
   above 500 ft with 1.4 km of DOE land still to spare.
2. **Lake Superior keeps the deepest water but pays 0.58 km of overflight.**
   Its 51.4 mrad falls short of the 62.5 mrad that its due-south near side
   demands. Still, rotating even 8° off the meridian is a large gain over the
   [previous study's](../rcs-aim/) az-0 aim: 319 m of column and a 248 km
   submerged run, against 213 m and a shoreline-limited crossing.
3. **Closer water is worse, not better.** Lake Michigan is 310 km away with
   240 m of water, but 1/L² makes the exit dose **4× higher** and the shallow
   24.4 mrad leaves 3.8 km of overflight. Distance is the mitigation.
4. **Ring fit is not the constraint.** Allowing the ring centre to move — it
   is a free design variable; the baseline already nests the three rings at
   different centres — every ring fits at **every** azimuth with >120 m
   clearance. RCS3/4 is tightest, and *gains* room when rotated: 158 m
   clearance at az 0 versus **305 m at az 36°**.

## 2. Why the near side is what matters

The asymmetry that makes this work is a property of the existing layout, not
of the beam. The corridor meridian runs along the site's eastern boundary, so
the IP has:

* **0.9 km** of DOE land due east,
* **2.6–2.8 km** due north and south (the ComEd and UIUC bearings),
* **4.7–4.9 km** to the W, WSW and SW.

θ_free is inversely proportional to that distance, so it varies by a factor
of five around the compass — from ~96 mrad on the SSE bearing to ~34 mrad on
the WSW one. Pointing the *far* beam ENE is exactly what puts the *near* beam
on the generous WSW diagonal. The physics of the far target and the land use
of the near emergence are coupled only through this 180° flip, and it happens
to favour the aim we want.

Two caveats stated plainly. The near-side emergence streak still carries
Sv-class raw dose in a decimetres-wide line ~300 m from the straight — on
site, fenced, occupancy-controlled, the same object every variant in this
study has. And **land use under the near-side band is unverified for every
bearing except the studied north (ComEd) and south corridors**; "zero
off-site overflight" is a statement about the 500 ft envelope, and the WSW
sector of the site still needs the segment-by-segment walk the
[envelope note](../envelope/) gives the north side.

## 3. What this does to the design

Nothing about the collider changes: it stays on the meridian aimed at UIUC,
because that aim is the physics programme and its 15.4 mrad *needs* the
corridor. What changes is that **the water-aimed ring is no longer tied to
the rest of the complex**:

* aim the collider and RCS1/2 south to UIUC on the ComEd meridian (the
  [co-tilted chain](../uiuc-chain/), which keeps ~97 % of the ν̄τ physics —
  the τ payoff lives in the low-energy turns);
* aim **RCS3/4 at Lake Huron, azimuth 50°, 49.5 mrad**, with its own ring
  centre and its own bearing;
* the vertical achromat that the [bent ring](../uiuc-chain/) already needs
  grows with tilt: 49.5 mrad is 826 T·m at 5 TeV, ~103 m of 8 T per end,
  1.4 % of the ring — the same class as the 864 T·m the meridian Lake
  Superior aim needs.

## 4. Open items

* **Land use on the WSW near-side band**, per above — the one thing that
  could still force a corridor.
* **Great Lakes jurisdiction**, now for Lake Huron: a boundary water under
  the Boundary Waters Treaty and the IJC rather than a domestic lake. This is
  a legal question and this study does not adjudicate it; it is plausibly
  harder than Lake Superior's US waters, which is a genuine argument for
  keeping the meridian aim and accepting 0.58 km of overflight.
* **Spherical geodesy** throughout; a real aim needs a WGS84 geodetic
  solution and NOAA lake bathymetry on the exit ellipse.
* Natural Earth's ocean polygon is deliberately excluded — its exterior ring
  is the world coastline, so point-in-polygon calls the whole continent
  ocean. The nearest real ocean water (Atlantic ~1200 km, Gulf ~1280 km)
  needs 94–101 mrad regardless.
* Arc-disk azimuthal doses, still owed by the baseline, now with rings on
  several different bearings.
