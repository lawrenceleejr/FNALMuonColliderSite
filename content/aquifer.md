---
title: Aquifer Radiation Audit
weight: 5
---

Northern Illinois drinks from the rock this beam tunnels through. This page
audits what the corridor beams do to every water-bearing unit they cross:
how much activation product is made, where, in what volume of rock, and how
the resulting concentrations compare to drinking-water standards. The rock
path and geometry are exactly those of the
[interactive tool](../tool/)'s baseline configuration (IP at 41.8443° N,
straightaway 191 m ASL — 34 m below grade — tilted 15.40 mrad); aquifer
units are tinted blue in the tool's cross-section.

**Headline:** the only place any groundwater standard could be approached is
a centimetres-wide filament within **2.1 km of the collision point — land
that is Fermilab's in every direction** (the boundary is ≥ 2.6 km away along
the beam). Every aquifer crossing beneath private or municipal land is
two to four orders of magnitude below the EPA drinking-water limit *before*
counting groundwater motion, which dilutes the numbers by further orders of
magnitude. Total tritium production along the entire 198 km chord is
**≈ 0.06 Ci/yr** — of the order of what a single hospital nuclear-medicine
department handles, spread through two hundred kilometres of deep rock.

## 1. Where the beam meets water

The corridor stratigraphy alternates aquifers and aquitards
(depths below grade; ±30%-ish, dipping gently south into the Illinois Basin):

| unit | depth | role along the corridor |
|---|---|---|
| Glacial drift | 0–25 m | sand-and-gravel lenses; the **Mahomet aquifer** valley fill near Champaign |
| Silurian dolomite | 25–110 m | the **shallow private-well aquifer** of Kane/DuPage (wells typ. 50–120 m) |
| Maquoketa shale | 110–170 m | regional aquitard — seals the Silurian from below |
| Galena–Platteville dolomite | 170–340 m | tight; hosts Fermilab's deep NuMI/MINOS halls |
| St. Peter sandstone | 340–470 m | **deep municipal sandstone aquifer** (Aurora, Naperville wellfields) |
| Prairie du Chien | 470–640 m | leaky dolomite |
| Ironton–Galesville | 640–880 m | the region's **best deep municipal aquifer** |
| Mt. Simon | 880–1400 m | brackish at depth; storage/injection horizon |

The IP straight itself sits **in the Silurian dolomite** (34 m below grade).
The south beam dives through the entire sequence to a perigee of **789 m**
(inside the Ironton–Galesville; it never reaches the Mt. Simon) and climbs
back out to the UIUC South Farms exit, 198 km away. The north beam crosses
only 2.2 km of shallow rock — all inside the site — before surfacing.

## 2. What activation the beam actually makes

A neutrino beam deposits nothing along its path except through the rare
interactions of the neutrinos themselves. The number is set by the flux and
the cross-section:

* decays aimed each way per year: **6.2 × 10¹⁸** (baseline machine, 700 m
  straight in the 11 km ring),
* on-axis flux at distance *L*: Φ = N γ² / (π L²), with γ = 47 300 at
  Eμ = 5 TeV,
* σν ≈ 1.1 × 10⁻³⁵ cm² per nucleon at ⟨Eν⟩ ≈ 3.3 TeV.

Along the **whole 198 km chord** (5.1 × 10⁷ g/cm², or 3.5 × 10⁻⁴ interaction
lengths) that gives **2.2 × 10¹⁵ interactions per year** — in total, in all
of the rock, most of them within the first few kilometres of the IP. Each
interaction is a ~3 TeV hadronic/EM shower whose dominant water-relevant
product is tritium (³H), at ~10⁻² atoms per GeV → **~30 ³H atoms per
interaction**. (Other nuclides — ⁷Be, ²²Na — are produced at lower yield,
are far less mobile in groundwater, and are bounded by the same geometry;
tritium, which travels *as* water, is the honest worst case. This is the
same nuclide hierarchy that governs NuMI/LBNF groundwater reviews at
Fermilab today.)

Total: **7 × 10¹⁶ ³H atoms per year ≈ 2.2 GBq ≈ 0.06 Ci/yr at saturation**,
distributed along 198 km of chord with a 1/L² weighting toward the IP.

## 3. Unit-by-unit audit

The conservative screening quantity is the **stagnant in-core saturation
concentration**: assume a parcel of groundwater sits motionless inside the
beam's 1/γ core at the unit's *closest approach* to the IP, for many
tritium mean-lives, with zero exchange. Real water fails all three
assumptions, so real concentrations are far lower.

| aquifer unit | path in unit | closest approach to IP | core radius there | stagnant in-core saturation | vs EPA MCL (740 Bq/L) |
|---|---|---|---|---|---|
| Silurian dolomite | 10.7 km | 0 (IP sits in it) | — | see §4 | exceeds only < 2.1 km from IP, on DOE land |
| St. Peter sandstone | 23.8 km | 22.4 km | 47 cm | **6.2 Bq/L** | 120× below |
| Ironton–Galesville | 87.3 km | 54.5 km | 115 cm | **1.1 Bq/L** | 700× below |
| glacial drift / Mahomet region | 1.6 km | 197 km | 4.2 m | **0.08 Bq/L** | 9 000× below |
| (aquitards: Maquoketa, Galena–Platteville, Prairie du Chien) | 75 km | 5.1 km | 11 cm | 121 Bq/L max | not potable units |

(EPA MCL for tritium: 20 000 pCi/L = 740 Bq/L. The Mt. Simon is never
touched.)

Every municipal aquifer crossing sits **two to four orders of magnitude
below the MCL even in this deliberately absurd stagnant-core limit.** The
deep sandstones that Aurora and the Fox Valley pump are crossed at 22 km and
54 km from the IP, where the beam core has spread to half-metre scale and
the flux has fallen four orders of magnitude from its on-site values.

## 4. The Silurian filament — the one real number

The IP straight lives in the Silurian, so close to the IP the stagnant-core
concentration is genuinely high: the in-core saturation falls through the
MCL at **L = 2.1 km** from the collision point. Three facts keep this a
non-issue, in order of importance:

1. **It is all on Fermilab land.** Along the beam the site boundary is
   2.6 km south and 2.8 km north of the IP; the 2.1 km MCL contour, in both
   directions, ends inside the fence. At the *south boundary* the beam is
   74 m deep and its in-core ceiling is ~450 Bq/L — under the MCL before it
   leaves DOE property. (The north beam has left the ground entirely by
   2.2 km, still on site.)
2. **The "filament" is centimetres wide.** The 1/γ core radius is 2 cm at
   1 km. The volume of Silurian water inside the over-MCL region is a few
   cubic metres — not an aquifer, a thread. Fermilab's existing radiological
   groundwater program (NuMI, LBNF target complex) manages *larger*
   activation source terms in the same dolomite today, with decades of
   monitoring-well precedent.
3. **Water moves.** Silurian groundwater advects of order metres per year;
   any parcel crosses the centimetre core in days-to-weeks, integrating a
   tiny fraction of the stagnant-limit dose, and is then diluted into the
   surrounding aquifer volume. The stagnant number is a screening ceiling,
   not a prediction — the realistic concentrations are orders of magnitude
   lower still.

## 5. Under Urbana and the South Farms

On the climb-out the beam re-enters the shallow section 24 km before the
exit, passing under Urbana at 66–119 m depth and under the UIUC campus into
the South Farms. The flux there is set by L ≈ 197 km — down four orders of
magnitude from the on-site values, in a beam core that has spread to four
metres across — giving a stagnant-limit saturation concentration of
**0.08 Bq/L** in the drift, four orders of magnitude below the MCL. The Mahomet aquifer proper lies west of the corridor and
deeper than the crossing; its assessment is bounded by the same number.

## 6. What would still be done

This audit is a screening argument, deliberately conservative at every step.
A real project would do what NuMI did: FLUKA/MARS transport of the shower
products in the actual rock chemistry, a hydrogeological transport model
(Darcy velocities from site borings), and monitoring wells over the first
kilometres of both beams — with the numbers above saying the program will
measure background.

*All inputs and formulas are in
[`tools/corridor_layout.py`](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/tools/corridor_layout.py)
and reproduced live by the [interactive tool](../tool/); this page's table
was generated with the same geometry code. σν and yield assumptions as in
the [safety assessment](../safety/); King, arXiv:physics/9908017 for the
dose framework; EPA 40 CFR 141 for the MCL.*
