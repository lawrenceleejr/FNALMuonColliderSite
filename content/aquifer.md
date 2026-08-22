---
title: Aquifer Radiation Audit
weight: 5
---

Northern Illinois drinks from the rock this beam tunnels through. This page
audits what the corridor beams do to every water-bearing unit they cross:
how much activation product is made, where, in what volume of rock, and how
the resulting concentrations compare to drinking-water standards. The
geometry is the [interactive tool](../tool/)'s baseline (IP at 41.8443° N,
straightaway 191 m ASL — 34 m below grade — tilted 15.40 mrad), traced
through a **dipping, latitude-dependent stratigraphic model** of the
corridor (§1) — the strata fall away south-eastward into the Illinois
Basin, so the rock under Champaign County is not the rock under Fermilab,
and this audit no longer pretends it is.

**Headline:** the only place any groundwater standard could be approached
is a decimetres-wide filament within **~200 m of the collision point —
inside the accelerator complex itself** (the site boundary is ≥ 2.6 km away
along the beam). This audit now uses the divergence-smeared plume of
[arXiv:2608.02718](https://arxiv.org/abs/2608.02718) (see the
[beam-size study](../beamsize/)): the on-axis interaction density falls
÷118 relative to the superseded pencil model, and every margin below
improves by that factor. The chord bottoms not in a drinking-water aquifer
but in the **brackish upper Mt. Simon** — the same horizon northern
Illinois uses for natural-gas storage — sealed under the Eau Claire
aquitard. Every potable-aquifer crossing beyond the fence is **at least
~5 000× below the EPA drinking-water limit** in a deliberately absurd
stagnant worst case (most are 10⁴–10⁶× below), *before* counting
groundwater motion, which lowers all of these by further orders of
magnitude. Total tritium production along the 198 km chord is unchanged at
**≈ 0.05 Ci/yr** — of the order of what a single hospital nuclear-medicine
department handles, spread through two hundred kilometres of deep rock.
(The dedicated-pencil scenario of [beam-size](../beamsize/) §4 restores the
previously published concentrations — the MCL contour returns to ~1.4 km,
still well inside the fence, and every off-site margin returns to the
"≥ ~45×" the earlier revision of this page defended.)

## 1. The real stratigraphy along the chord

The Paleozoic section dips gently south-east into the Illinois Basin —
roughly 1–1.5 m/km near Fermilab, steepening to 3–6 m/km down the east
flank of the La Salle Anticlinorium south of ~41°N. Over the corridor's
198 km that amounts to a **400–650 m structural drop** in the deep-unit
tops, while younger cover (Pennsylvanian, then Devonian–Mississippian
shales and carbonates) appears above the truncated older section south of
~41.3° N. Two structures interrupt the smooth picture: the **Sandwich
Fault Zone**, crossed at ~41.55–41.65° N (throw plausibly 30–150 m here,
down-to-northeast), and the **La Salle Anticlinorium's** east flank, which
the southern half of the corridor rides. The glacial drift is equally
non-uniform: ~15–35 m at Fermilab, thin over the bedrock high near
41.2–41.4° N, and swelling to **80–150 m** where the corridor crosses the
buried **Mahomet Bedrock Valley** (~40.25–40.5° N) and into
Champaign–Urbana, where ~90–100 m of drift sits on **Pennsylvanian shale
bedrock** — the Silurian dolomite that hosts the IP straight lies ~320 m
deep down there.

What the south beam actually traverses, in order (distances from the IP;
the [tool](../tool/) recomputes this live and its rock table now uses the
same model):

| km from IP | unit entered | role |
|---|---|---|
| 0 | Silurian dolomite (IP at the drift/rock contact, 34 m) | shallow domestic-well aquifer of Kane/DuPage (wells typ. 30–90 m; karstic) |
| 1.4 | Maquoketa shale | the regional aquitard |
| 5.6 | Galena–Platteville dolomite | tight; hosts NuMI/MINOS |
| ~14 | St. Peter sandstone | **deep municipal sandstone aquifer** (Fox Valley wellfields) |
| ~19 | Prairie du Chien–Franconia | leaky dolomites/siltstones |
| ~28 | Ironton–Galesville | **the region's best deep municipal aquifer** (only ~40–60 m thick) |
| ~31 | Eau Claire shale | aquitard sealing the Mt. Simon |
| **~48–108** | **Mt. Simon sandstone** | **brackish; the gas-storage horizon (Herscher, Ancona, Manlove); perigee 789 m at ~40.94° N** |
| ~125 / ~142 / ~150 / ~161 | Ironton–Galesville / St. Peter / Galena–Platteville / Maquoketa again | climbing the basin's south flank (mostly non-potable at these depths/latitudes) |
| ~169 | Silurian again | deep and mineralized here, not a supply |
| ~176 / ~179 | Devonian carbonate / New Albany Shale | clipped once on the ascent |
| ~182 | Pennsylvanian (Bond–Mattoon shales) | Champaign County bedrock |
| ~192 | glacial drift, incl. basal (Mahomet-class) sands | **the beam's one drift-aquifer passage, 95 m under northwest Urbana** |
| 197.4 | surface — UIUC South Farms | exit path entirely in drift |

Where the model comes from: the Fermilab column is anchored to NuMI/MINOS
civil-construction geology (the MINOS hall sits in the uppermost
Galena–Platteville at ~100 m); the Champaign end to ISGS drift and bedrock
mapping and the Mahomet Bedrock Valley studies; the deep mid-corridor tops
to the gas-storage fields (Herscher, Ancona, Manlove) and the regional
Cambrian–Ordovician aquifer studies. Control points are interpolated at
seven latitudes; mid-corridor tops carry **±75–150 m** uncertainty, so the
crossing distances above are good to ~±30 % until real well logs along the
meridian replace the interpolation — a one-day ISGS-records task flagged in
§6. (The previous revision of this page hung Fermilab's column, flat, from
the local surface for all 198 km; that got the unit wrong for over half
the chord — including claiming the perigee was in the Ironton–Galesville
aquifer when it is in fact ~100 m into the brackish Mt. Simon below it.)

## 2. What activation the beam actually makes

A neutrino beam deposits nothing along its path except through the rare
interactions of the neutrinos themselves. The number is set by the flux and
the cross-section:

* decays aimed each way per year: **5.3 × 10¹⁸** (baseline machine, 700 m
  straight in the 11 km ring, 85.4 % store-decay factor),
* on-axis flux at distance *L*: Φ = N γ² / (π L²), with γ = 47 300 at
  Eμ = 5 TeV,
* σν ≈ 1.1 × 10⁻³⁵ cm² per nucleon per neutrino, averaged over the
  angle-integrated spectrum (⟨Eν⟩ ≈ 1.6 TeV — whole-plane totals never
  depended on the prism).

Along the whole 198 km chord (5 × 10⁷ g/cm² at the corrected mean density
of ~2.5 g/cm³, or 3.3 × 10⁻⁴ interaction lengths) that gives
**≈ 2 × 10¹⁵ interactions per year** — in total, in all of the rock, most
of them within the first few kilometres of the IP. Each interaction is a
~3 TeV hadronic/EM shower whose dominant water-relevant product is tritium
(³H), at ~10⁻² atoms per GeV → **~30 ³H atoms per interaction**. (Other
nuclides — ⁷Be, ²²Na — are produced at lower yield, are far less mobile in
groundwater, and are bounded by the same geometry; tritium, which travels
*as* water, is the honest worst case. This is the same nuclide hierarchy
that governs NuMI/LBNF groundwater reviews at Fermilab today.)

Total: **6 × 10¹⁶ ³H atoms per year ≈ 1.9 GBq ≈ 0.05 Ci/yr at saturation**,
distributed along 198 km of chord with a 1/L² weighting toward the IP.

## 3. Unit-by-unit audit

The conservative screening quantity is the **stagnant in-core saturation
concentration**: assume a parcel of groundwater sits motionless inside the
beam's 1/γ core at the unit's *closest approach* to the IP, for many
tritium mean-lives, with zero exchange. Real water fails all three
assumptions, so real concentrations are far lower. With the dipping
geology, every deep unit is crossed **closer to the IP** than the flat
model claimed — the margins below are honest and smaller than previously
published here, and they still hold.

| unit | path in unit | closest approach | stagnant in-core saturation | vs EPA MCL (740 Bq/L) |
|---|---|---|---|---|
| Silurian dolomite (potable reach) | ~1.4 km on-site (+ deep mineralized re-crossing at 169 km) | 0 (IP sits at its top) | see §4 | exceeds only < ~0.2 km from IP, inside the complex |
| Maquoketa shale (aquitard) | ~12 km | 1.4 km | ~13 Bq/L | ~60× below — and a shale aquitard, not a water supply |
| Galena–Platteville (tight) | ~20 km | 5.6 km | ~0.9 Bq/L | ~900× below; not a supply unit here |
| **St. Peter sandstone** | ~13 km | **~14 km** (Kendall Co., under the Fox Valley cone of depression) | **~0.14 Bq/L** | **~5 000× below** |
| **Ironton–Galesville** | ~8 km | **~28 km** | **~0.03 Bq/L** | **~2 × 10⁴× below** |
| Mt. Simon (brackish/storage) | ~60 km | ~48 km | ~0.01 Bq/L | ~7 × 10⁴× below — and not potable |
| basal drift / Mahomet-class sands | ~6 km | ~192 km (NW Urbana) | **~8 × 10⁻⁴ Bq/L** | **~10⁶× below** |

The two aquifers the Fox Valley municipalities actually pump — the
St. Peter and the Ironton–Galesville — are crossed at ~14 and ~28 km, not
the 22 and 54 km the flat model gave, and squarely beneath the
Aurora–Oswego–Yorkville–Joliet **cone of depression** (potentiometric
declines exceeding 200 m; regional models project partial St. Peter
desaturation by mid-century). With the measured plume the honest margin statement is
"**≥ ~5 000× below the MCL in the stagnant limit**" (the dedicated-pencil
scenario returns it to the "≥ ~45×" this page previously defended) — over a
decimetres-to-metres-wide core, before any water movement, in units already
managed as a declining resource for reasons that have nothing to do with
this beam.

## 4. The Silurian filament — the one real number

The IP straight sits at the drift/Silurian contact (Fermilab's drift is
15–35 m thick; borings, not this model, would fix which side of the
contact the hall floor is on). Close to the IP the stagnant-core
concentration is genuinely high: the in-core saturation falls through the
MCL at **L ≈ 0.2 km** from the collision point (2.1 km in the superseded
pencil model; ~1.4 km in the dedicated-pencil scenario). Three facts keep this a
non-issue, in order of importance:

1. **It is all on Fermilab land.** Along the beam the site boundary is
   2.6 km south and 2.8 km north of the IP; the MCL contour, in both
   directions and in every scenario up to the full dedicated pencil, ends
   inside the fence. At the *south boundary* the beam is 74 m deep —
   already through the Maquoketa's roof and below the base of most domestic
   Silurian wells (typ. 30–90 m) — and its in-core ceiling is ~4 Bq/L
   (~450 Bq/L in the dedicated-pencil scenario), under the MCL before it
   leaves DOE property either way. (The north beam
   has left the ground entirely by 2.2 km, still on site.)
2. **The "filament" is decimetres wide.** The smeared-core radius is
   1.2 σθL ≈ 18 cm at 1 km — and the over-MCL reach now ends at ~200 m,
   where the core is ~4 cm. The volume of Silurian water inside the
   over-MCL region is under a cubic metre — not an aquifer, a thread. Fermilab's existing radiological
   groundwater program (NuMI, LBNF target complex) manages *larger*
   activation source terms in this same geology today, with decades of
   monitoring-well precedent.
3. **Water moves — here, fast.** The Niagaran Silurian of northeastern
   Illinois is crevice-and-solution-porosity rock: groundwater travels in
   fractures at **metres per day**, not metres per year. Any parcel crosses
   the centimetre-wide core in minutes-to-hours, integrating a vanishing
   fraction of the stagnant-limit dose, and is then diluted into the
   fracture network. (Fracture flow also means monitoring wells must be
   sited on the fracture fabric — a design input for §6, not a hazard.)
   The stagnant number is a screening ceiling, not a prediction.

## 5. Under Urbana, the Mahomet question, and the South Farms

Champaign–Urbana's water supply comes from wells in the **Mahomet aquifer**
— the sand-filled buried bedrock valley system whose main trunk the
corridor crosses at ~40.25–40.5° N. There the beam is **350–550 m beneath
the valley fill**: no contact at all where the aquifer is thickest and most
pumped. The beam re-enters the drift section only at ~192 km, crossing
basal-drift sands of the same family at **~66–100 m depth under northwest
Urbana** — the same depth class as the region's wells (Champaign–Urbana
sits at the valley system's southern margin; whether the specific sands the
beam threads are mapped Mahomet-member or till-bounded margin sands is a
well-log question). The number is what matters: at L ≈ 192 km the
stagnant-limit saturation is **~0.09 Bq/L — about 8 000× below the MCL** —
in a beam core four metres across, before any flow or dilution. The final
5 km to the South Farms exit run through drift and weathered Pennsylvanian
shale (not dolomite, as the flat model had it — the Geant4 rock-column
export now reflects this).

## 6. What would still be done

This audit is a screening argument, deliberately conservative at every
step, on an interpolated geological model with stated (±75–150 m
mid-corridor) error bars. A real project would, in order of cheapness:
pull the ~15 ISGS/ISWS well logs nearest the meridian to pin the unit tops
(a day's records work that would collapse the crossing-distance
uncertainty from ±30 % to a few percent); run FLUKA/MARS transport of the
shower products in the actual rock chemistry; build a fracture-flow
transport model for the on-site Silurian (Darcy velocities from site
borings, oriented to the karst fabric); and set monitoring wells over the
first kilometres of both beams — which the numbers above say will measure
background.

*Geometry and activation formulas are those of
[`tools/corridor_layout.py`](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/tools/corridor_layout.py)
and the [interactive tool](../tool/), which traces the same dipping model
live. Stratigraphic control: ISGS Bulletin 95 (Willman et al. 1975), the
ISGS 1:500 000 bedrock map (Kolata 2005), ISGS Circular 479 (Sandwich
Fault Zone), ISGS Bulletin 100 (Nelson 1995), the ISGS/ISWS
Cambrian–Ordovician aquifer studies, GSA Special Paper 258 and ISWS
reports on the Mahomet Bedrock Valley, gas-storage field records
(Herscher/Ancona/Manlove), and Fermilab NuMI/MINOS construction geology.
σν and yield assumptions as in the [safety assessment](../safety/); King,
arXiv:physics/9908017 for the dose framework; EPA 40 CFR 141 for the MCL.*
