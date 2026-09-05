---
title: "The Ring Stack: Depth, Tilt, and the Transfer Lines Between"
weight: 11
---

Three rings share this site — RCS1/2 (C = 6.28 km), the collider (11.00 km)
and RCS3/4 (14.72 km) — and every one of them puts its **east straight on the
corridor meridian**. That single constraint decides all three depths, because
a racetrack is planar: tilt it and the whole north–south reach turns into
depth. This page draws that, and then designs the lines that carry the beam
from one ring to the next. Tool:
[`ring_stack.py`](https://github.com/lawrenceleejr/FNALMuonColliderSite/blob/main/tools/ring_stack.py);
machine-readable [ring_stack.json](../geo/ring_stack.json).

<div style="margin:1.2rem 0"><img src="../figs/ring_stack.svg" alt="Plan of the three nested racetracks, north-south sections of the co-tilted and corridor-baseline configurations, and the closest approach of each pair of tunnels" style="max-width:100%"></div>

## 1. A planar ring turns tilt into depth

Because each ring is planar and tilted only about an east–west axis, its
elevation depends on the north coordinate alone:

<div style="text-align:center;font-family:monospace;margin:.8rem 0">
z(y) = z₀ + θ (y − y<sub>c</sub>) + (y − y<sub>c</sub>)² / 2R<sub>E</sub>
</div>

so the ring projects onto a north–south section as a **single straight line**
— the east straight occupying y<sub>c</sub> ± L<sub>s</sub>/2, the two arcs
sweeping out to y<sub>c</sub> ± (L<sub>s</sub>/2 + R). One section panel
therefore carries depth, tilt, reach and rock cover for all three rings at
once. The quoted "depth" of a ring is the depth *at its straight*; what the
tilt does to the rest of it is the swing column.

| ring | C | straight | arc R | N–S reach | tilt | depth at the straight | shallowest | deepest | swing |
|---|---|---|---|---|---|---|---|---|---|
| **co-tilted chain** ||||||||||
| RCS1/2 | 6.28 km | 500 m | 841 m | 2 182 m | +15.405 mrad | 32.6 m | **15.9 m** (N apex) | 47.5 m (S apex) | 31.6 m |
| collider | 11.00 km | 700 m | 1 528 m | 3 756 m | +15.400 | 33.6 m | **5.9 m** (N apex) | 59.7 m | 53.8 m |
| RCS3/4 | 14.72 km | 450 m | 2 200 m | 4 850 m | +15.312 | 45.1 m | **8.9 m** (N apex) | 81.6 m | 72.7 m |
| **corridor baseline** ||||||||||
| RCS1/2 | | | | | −4.32 mrad | 60 m | 53.4 m | 64.9 m | 11.5 m |
| collider | | | | | 0 (level) | 100 m | 96.2 m | 101.3 m | 5.0 m |
| RCS3/4 | | | | | −2.16 mrad | 80 m | 73.4 m | 86.1 m | 12.7 m |

Three readings.

**The 200 m ceiling never binds a ring.** The deepest point in the whole
complex is RCS3/4's south arc at **81.6 m**, well inside the
[ceiling](../collider-map/#the-200-m-ceiling-and-what-it-does-and-does-not-limit).
Co-tilting is expensive in depth *swing* — 73 m across RCS3/4 against 13 m in
the corridor baseline — but it starts from a much shallower straight (45 m
against 80 m), and the two effects nearly cancel at the deep end.

**What binds is the shallow end.** Tilt up-to-north and the north arc apex
rises out of the ground: the collider's sits under **5.9 m** of cover, RCS3/4's
under 8.9 m. That is the constraint that closes RCS3/4's
[tilt window](../uiuc-chain/) to 29 µrad, and it is the reason the co-tilted
chain cannot simply be made shallower to buy transfer-line room. Note that
5.9 m is already below the 6 m minimum the variant sets itself — the collider,
at exactly 15.400 mrad, sits a hair outside its own window (15.359–15.399 mrad
by the same test).

**Co-tilting squeezes the three tunnels together.** This is the finding the
section panel makes obvious and that the tilt table hides. Because all three
tilts agree to within 0.1 mrad, the three ring *planes* are parallel, and they
end up inside a **6.2 m elevation band** that stays 6.2 m wide across the whole
5 km footprint — RCS1/2 and the collider are only **0.99 m** apart in
elevation. In the corridor baseline the same three planes are 20–40 m apart.

## 2. The tunnels cross, and there is almost nothing between them

Three nested racetracks that share a tangent line must cross each other in
plan: with the straights ordered west to east as below, the collider's south
arc crosses RCS3/4's east straight at y = −548 m and RCS1/2's north arc crosses
RCS3/4's at y = +70 m. Where two tunnels cross, the
horizontal offset between them is zero **by definition**, so the only thing
keeping them apart is the elevation band above.

With a 5.5 m bore and the east straights 13 and 26 m off the meridian:

| pair | corridor baseline | co-tilted chain |
|---|---|---|
| RCS1/2 – collider | 46.8 m | 26.0 m (nested, no crossing) |
| collider – RCS3/4 | 15.6 m | **5.2 m** (crossing) |
| RCS1/2 – RCS3/4 | 18.7 m | **6.2 m** (crossing) |

5.2 m between the axes of two 5.5 m bores means the bores touch. This is a
real civil constraint that the co-tilted chain has not previously been
costed against, and it is *tighter than the 5.9 m arc cover*. Three things
can be done about it, and the numbers for each are in the JSON:

1. **Use the tilt windows.** Each ring's tilt is free over a window
   (RCS1/2 96 µrad, collider 40, RCS3/4 29), and moving within it slides z₀ by
   ~198 m per mrad. Pushing RCS1/2 to the shallow edge and RCS3/4 to the deep
   edge opens the planes to **8.2 m** — about 1.5 bore diameters. But that
   spends every margin at once: RCS3/4's north pencil then emerges 2 m inside
   the fence. Demanding 8 m of cover and 150 m of fence margin instead leaves
   only 4.0 m between planes, which is worse than doing nothing.
2. **Spread the rings in longitude.** The collider nests cleanly inside
   RCS3/4 — no crossing at all — once their straights are about **250 m** apart
   east–west instead of 13; at 260 m the closest approach jumps from 5.3 m to
   13.5 m, and at 300 m to 37 m. That trades the crossings for a corridor
   250 m wide at UIUC rather than 26 m, which the exit-strip analysis can
   absorb but which is a different design.
3. **Engineer the crossings.** There are only two of them, both known to the
   metre, and a locally deepened section or a shared cavern is ordinary
   practice. This is probably the right answer, but it is a cost line, not a
   free pass.

**Ordering the straights matters, and it is free.** The published layout
specifies 10–25 m east–west offsets but not their order. Of the six
orderings, **RCS1/2 < RCS3/4 < collider** (west to east) is the one to take:
it is the *chain* order, so no transfer line has to cross a third ring's
straight on the way, and it makes RCS1/2 nest inside the collider — the pair
with only 0.99 m of elevation between them, and so the pair that can least
afford to cross. Putting the collider in the middle instead drives the
750 GeV transfer line to within **2.0 m** of it.

## 3. Transfer lines

<div style="margin:1.2rem 0"><img src="../figs/transfer_lines.svg" alt="Plan and section of the transfer-line region, showing the mu-plus and mu-minus lines between RCS1/2, RCS3/4 and the collider for the as-built and concentric ring placements" style="max-width:100%"></div>

The chain is RCS1 → RCS2 → RCS3 → RCS4 → collider. Two of those four
transfers stay inside their own tunnel: RCS1/RCS2 and RCS3/RCS4 are stacked
pairs sharing one bore, so the transfer is a plain vertical dogleg — 1.5 m
over 200 m, 15 mrad, **16 T·m at 314 GeV** and **75 T·m at 1.5 TeV**, or 2 and
9 m of 8 T dipole. Nothing to site.

The two inter-ring transfers are the interesting ones, and the corridor
geometry sets their character. **The east straights are the only place any two
rings run parallel and close** — 13 m apart there, more than 180 m apart
everywhere else — so every transfer is a near-pure **translation**, which is
the most expensive kind: an offset *d* over a line of length *L* costs 2*d*/*L*
of bending at both ends, and

<div style="text-align:center;font-family:monospace;margin:.8rem 0">
∫B·dl = (2d/L) · p / 0.29979 &nbsp;&nbsp; [T·m, p in GeV/c]
</div>

At 5 TeV, p/0.29979 = 16 700 T·m, so a single metre of offset taken over a
100 m line already costs **334 T·m** — 42 m of 8 T dipole. The lines want to
be long.

| transfer | sign | p | length | offset (E / elevation) | bend | ∫B·dl | at 8 T | at 1.8 T |
|---|---|---|---|---|---|---|---|---|
| RCS2 → RCS3 | μ⁺ | 750 GeV | 745 m | +13 m / −17.6 m | 52 mrad | 129 T·m | 16 m | 72 m |
| RCS2 → RCS3 | μ⁻ | 750 GeV | — | — | — | **no room** | | |
| RCS4 → collider | μ⁺ | 5 TeV | — | — | — | **no room** | | |
| RCS4 → collider | μ⁻ | 5 TeV | 845 m | +13 m / +18.2 m | 43 mrad | 719 T·m | 90 m | 399 m |

**Half the lines have nowhere to go, and the reason is 390 m of latitude.**
μ⁺ and μ⁻ counter-rotate in the same ring, so each transfer needs a *mirror
pair* of lines — one leaving the south end of a straight, one the north. Those
two are mirror images only if the two rings share a straight centre. They do
not: RCS3/4 is the site filler and sits as far north as it will go, which is
still **390 m south** of the collider and RCS1/2. Its straight
(y = −615 … −165 m) therefore overlaps the collider's (−350 … +350 m) by only
185 m, and for one sign of each transfer the beam would have to travel
*backwards* to reach its injection septum. The best that sign can be given is
a 66 m line, against the ~120 m a septum, dogleg and matching section need.

**The fix is to make the three rings concentric.** Slide the collider and
RCS1/2 390 m south onto RCS3/4's straight centre. Both are far smaller than
RCS3/4, which already fits the site with margin, so nothing leaves the
boundary. Then:

| transfer | sign | length | bend | ∫B·dl | at 8 T |
|---|---|---|---|---|---|
| RCS2 → RCS3 | μ⁺ and μ⁻ | 355 m each | 126 mrad | 314 T·m | 39 m |
| RCS4 → collider | μ⁺ and μ⁻ | 455 m each | 74 mrad | 1 239 T·m | 155 m |

All four lines close, as exact mirror pairs. The 750 GeV line's clearance to
the collider — the ring it passes but does not connect — improves from 5.5 m
to 13.6 m; the 5 TeV line's clearance to RCS1/2 slips from 18.4 m to 16.0 m.
The costs are
real but modest: the lines are shorter, so each carries roughly 2–2.5× the
bending of the single long line it replaces, and the 5 TeV pair needs 155 m of
8 T dipole apiece — 1.1 % of RCS3/4's circumference, the same class as the
[vertical achromats](../rcs-azimuth/) the co-tilted chain already needs.

Two side effects of the shift, both benign and one useful:

* the IP moves 390 m south, so the UIUC exit moves 390 m — 0.2 % of the
  198 km baseline, and the far hall follows it;
* the collider's north-arc cover improves from **5.9 m to 12.9 m**, which
  fixes the one place the co-tilted chain breaks its own 6 m rule. RCS1/2's
  and RCS3/4's north pencils still emerge inside the fence, at 835 m and
  201 m.

The one thing the concentric shift does *not* fix is the plane crossing: at
window midpoints the collider and RCS3/4 come within 3.9 m of each other, a
little worse than the 5.2 m as built. Crossings stay a civil problem to be
engineered, not a geometry problem to be solved.

## 4. What is not designed here

* **Optics.** These are geometric lines: length, bend and ∫B·dl. No quadrupole
  layout, no dispersion matching, no aperture from the RCS emittance, and no
  check that the 126 mrad doglegs fit inside the chromatic acceptance of a
  muon beam at 750 GeV.
* **Extraction and injection hardware.** A 40 m septum allowance upstream and
  80 m for the injection kicker are assumed and reserved, not specified. Fast
  kickers for a 5 Hz, few-turn muon machine are their own problem.
* **The front end.** The chain here starts at RCS1's 63 GeV injection; the
  muon source, cooling channel and linac that feed it are not sited in this
  study, so the injector → RCS1 line is not drawn.
* **Whether the lines can share bores.** None of them ever gets more than 8 m
  from every ring at once — the whole transfer region is inside a 26 m-wide,
  20 m-deep envelope — so the realistic picture is a small number of wide
  enclosures rather than four independent tunnels.
