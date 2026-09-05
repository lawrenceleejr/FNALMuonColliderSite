#!/usr/bin/env python3
"""Depth and tilt of the three rings, and the transfer lines between them.

Every ring in this study is a planar racetrack whose EAST straight lies on
the corridor meridian, so the ring plane is tilted only about an east-west
axis and its elevation depends on the north coordinate alone:

    z(y) = z0 + theta (y - y_c) + (y - y_c)^2 / 2 R_E

(the quadratic term is the tangent plane rising away from the sphere).  A
tilted ring therefore projects onto a north-south section as a single
straight line, with the east straight occupying y_c +- Ls/2 and the two arcs
sweeping out to y_c +- (Ls/2 + R).  That makes one elevation panel enough to
carry depth, tilt, arc reach and cover for all three rings at once.

Two configurations are drawn:

  * "corridor baseline" (tools/corridor_layout.py): the collider deep and
    level, the two RCS planes pitched DOWN-to-north (4.32 and 2.16 mrad) so
    their pencils converge on the deep reference hall at 9.25 km.
  * "co-tilted chain" (tools/uiuc_chain.py): every ring tilted UP-to-north
    by ~15.4 mrad so that every south pencil exits at the UIUC South Farms.

The co-tilted case is the interesting one: because all three tilts agree to
0.1 mrad, the three ring planes end up inside a 6.2 m elevation band, and
wherever two rings cross in plan the only separation left is that band.
This tool measures those crossings, asks how far the published tilt windows
can prise the planes apart, and then designs the transfer lines.

Transfer lines.  The east straights are the only place any two rings run
parallel and close (they are 13-26 m apart there and > 180 m apart
everywhere else), so every transfer is a near-pure TRANSLATION: two equal
and opposite bends, costing 2 d / L of bending for an offset d over a line
of length L.  Because mu+ and mu- counter-rotate, each transfer needs a
mirror pair of lines, and because RCS3/4's straight sits 390 m south of the
other two, the two members of a pair get very different lengths.  Both are
sized here.

Output: static/figs/ring_stack.{svg,pdf}, static/figs/transfer_lines.{svg,pdf},
static/geo/ring_stack.json.
"""
import contextlib, io, json, math, os, re, sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
with contextlib.redirect_stdout(io.StringIO()):
    import corridor_layout as cl
import compass_common as cc
STYLE = "/root/.claude/skills/synced/tufte-plots/assets/tufte.mplstyle"
if os.path.exists(STYLE):
    plt.style.use(STYLE)
plt.rcParams["figure.dpi"] = 130

R_E = 6371000.0
IP_LAT = 41.8443
MLAT = cl.MLAT
FENCE_LAT = 41.8699
FENCE_Y = (FENCE_LAT - IP_LAT) * MLAT
DEPTH_CAP = cc.DEPTH_CAP_M                       # 200 m: no ring deeper than this
COVER_MIN = 6.0                                  # m of rock over the north-arc apex
LAT = lambda y: IP_LAT + np.asarray(y, float) / MLAT

# ---------------------------------------------------------------- terrain
_tool = open(os.path.join(ROOT, "static", "tool", "index.html")).read()
SOUTH_PROF = json.loads(re.search(r"const PROFILE_M = (\[\[.*?\]\]);", _tool).group(1))

def terrain_lat(lat):
    if lat >= 41.79:
        return cl.elev_at_lat(lat)
    for (la1, e1), (la2, e2) in zip(SOUTH_PROF, SOUTH_PROF[1:]):
        if la1 <= lat <= la2:
            return e1 + (lat - la1) / (la2 - la1) * (e2 - e1)
    return SOUTH_PROF[0][1]

def TER(y):
    y = np.asarray(y, float)
    return np.array([terrain_lat(IP_LAT + v / MLAT) for v in y.ravel()]).reshape(y.shape) if y.ndim \
        else terrain_lat(IP_LAT + float(y) / MLAT)

def beam_z(z0, th, y0, y):
    d = np.asarray(y, float) - y0
    return z0 + th * d + d * d / (2 * R_E)

def surface_crossing(z0, th, y0, ylo, yhi):
    f = lambda y: beam_z(z0, th, y0, y) - TER(y)
    a, b, fa = ylo, yhi, None
    fa = f(a)
    for _ in range(200):
        m = 0.5 * (a + b)
        if (f(m) > 0) == (fa > 0):
            a, fa = m, f(m)
        else:
            b = m
    return 0.5 * (a + b)

# ---------------------------------------------------------------- rings
#   name: (straight centre y, straight length, arc radius, circumference, colour)
RING = {
    "RCS1/2":   dict(Ls=cl.LS_RCS12,    R=cl.R_RCS12,    C=cl.C_RCS12,    col="#1f6f8b"),
    "collider": dict(Ls=cl.LS_COLLIDER, R=cl.R_COLL,     C=cl.C_COLLIDER, col="#0b6e4f"),
    "RCS3/4":   dict(Ls=cl.LS_RCS34,    R=cl.R_RCS34,    C=cl.C_RCS34,    col="#b5541c"),
}
ORDER = ["RCS1/2", "collider", "RCS3/4"]
# Straight-centre latitudes, in metres north of the IP.  "as built" is corridor_layout's own
# site fit: RCS3/4 is the site filler and sits as far north as it will go, 390 m south of the
# other two.  "concentric" slides the collider and RCS1/2 south onto RCS3/4's centre, which
# costs nothing in site fit (both are far smaller) and buys mirror-symmetric transfer lines.
Y_C = {"as built":   {"RCS1/2": 0.0, "collider": 0.0, "RCS3/4": -390.0},
       "concentric": {"RCS1/2": -390.0, "collider": -390.0, "RCS3/4": -390.0}}
PLACE = "as built"

def racetrack(nm, x_off=0.0, n=3000, place=None):
    """Closed plan centreline, in beam order for mu+ (east straight north->south):
    east straight N->S, south arc E->W, west straight S->N, north arc W->E."""
    r = RING[nm]; yc, Ls, R = Y_C[place or PLACE][nm], r["Ls"], r["R"]
    yn, ys = yc + Ls / 2, yc - Ls / 2
    m = max(3, int(n * Ls / (2 * Ls + 2 * math.pi * R)))
    na = max(8, int(n * math.pi * R / (2 * Ls + 2 * math.pi * R)))
    t = np.linspace(0, 1, m); a = np.linspace(0, math.pi, na)
    seg = [(np.full(m, x_off),               yn - t * Ls),                 # east straight, N -> S
           (x_off - R + R * np.cos(a),       ys - R * np.sin(a)),          # south arc, E -> W
           (np.full(m, x_off - 2 * R),       ys + t * Ls),                 # west straight, S -> N
           (x_off - R - R * np.cos(a),       yn + R * np.sin(a))]          # north arc, W -> E
    return np.concatenate([s[0] for s in seg]), np.concatenate([s[1] for s in seg])

def ring3d(nm, th, z0, x_off, n=3000, place=None):
    x, y = racetrack(nm, x_off, n, place)
    return x, y, beam_z(z0, th, Y_C[place or PLACE][nm], y)

# ---------------------------------------------------------------- the two configurations
# A. corridor baseline: depths below grade at the straight centre, RCS planes pitched down-north
BASE = {"RCS1/2":   dict(theta=-cl.PITCH_RCS12, depth=cl.DEPTH_RCS12),
        "collider": dict(theta=0.0,             depth=cl.DEPTH_COLLIDER),
        "RCS3/4":   dict(theta=-cl.PITCH_RCS34, depth=cl.DEPTH_RCS34)}
for nm, d in BASE.items():
    d["yc"] = Y_C["as built"][nm]
    d["z0"] = TER(d["yc"]) - d["depth"]

# B. co-tilted chain: tilt fixes z0 through the common UIUC exit
Y_EXIT = surface_crossing(191.0, 15.40e-3, 0.0, -215000, -150000)
Z_EXIT = TER(Y_EXIT)

def z0_for(th, yc):
    d = Y_EXIT - yc
    return Z_EXIT - th * d - d * d / (2 * R_E)

def ring_state(nm, th, place=None):
    r = RING[nm]; yc = Y_C[place or PLACE][nm]
    z0 = z0_for(th, yc); yap = yc + r["Ls"] / 2 + r["R"]
    return dict(theta=th, z0=z0, yc=yc, depth=TER(yc) - z0, y_apex=yap,
                cover=TER(yap) - beam_z(z0, th, yc, yap),
                y_emerge=surface_crossing(z0, th, yc, yc + 200, yc + 15000))

def tilt_window(nm, lo=15.20e-3, hi=15.60e-3, step=0.0005e-3, place=None):
    ok = [th for th in np.arange(lo, hi, step)
          if ring_state(nm, th, place)["cover"] >= COVER_MIN
          and ring_state(nm, th, place)["y_emerge"] < FENCE_Y]
    return (ok[0], ok[-1]) if ok else None

WIN = {nm: tilt_window(nm) for nm in ORDER}
_U = json.load(open(os.path.join(ROOT, "static", "geo", "uiuc_chain.json")))


# ---------------------------------------------------------------- E-W offsets and clearance
# Tunnels need separating.  Take a 5.5 m bore and ask for 2 diameters axis-to-axis where we can
# get it.  The published layout puts the straights 10-25 m apart but does not say in which order;
# of the six orderings, RCS1/2 < RCS3/4 < collider is the one to pick, for two reasons.  It is the
# CHAIN order, so no transfer line has to cross a third ring's straight on its way (the worst
# ordering leaves the collider sitting between RCS1/2 and RCS3/4, and the 750 GeV line then passes
# 2.0 m from it).  And it makes RCS1/2 nest inside the collider, which is the pair with only 1 m of
# elevation between them and so the pair that cannot afford to cross.
BORE = 5.5
X_OFF = {"RCS1/2": -26.0, "RCS3/4": -13.0, "collider": 0.0}

def pair_min(cfg, a, b, n=4000, place=None):
    A = ring3d(a, cfg[a]["theta"], cfg[a]["z0"], X_OFF[a], n, place)
    B = ring3d(b, cfg[b]["theta"], cfg[b]["z0"], X_OFF[b], n, place)
    d = np.sqrt((A[0][:, None] - B[0][None, :]) ** 2 + (A[1][:, None] - B[1][None, :]) ** 2
                + (A[2][:, None] - B[2][None, :]) ** 2)
    i = np.unravel_index(np.argmin(d), d.shape)
    return dict(pair="%s-%s" % (a, b), min_sep_m=float(d[i]), at_y_m=float(A[1][i[0]]),
                at_x_m=float(A[0][i[0]]), dz_m=float(A[2][i[0]] - B[2][i[1]]))

PAIRS = [("RCS1/2", "collider"), ("collider", "RCS3/4"), ("RCS1/2", "RCS3/4")]

# how far can the published tilt windows prise the three planes apart?
def z_at0(nm, th, place=None):
    yc = Y_C[place or PLACE][nm]
    return beam_z(z0_for(th, yc), th, yc, 0.0)

def best_spread(cover_min, fence_min):
    grid = {}
    for nm in ORDER:
        lo, hi = WIN[nm]
        grid[nm] = [th for th in np.arange(lo, hi + 1e-12, 0.001e-3)
                    if ring_state(nm, th)["cover"] >= cover_min
                    and FENCE_Y - ring_state(nm, th)["y_emerge"] >= fence_min]
    best = None
    for t1 in grid["RCS1/2"]:
        for t2 in grid["collider"]:
            for t3 in grid["RCS3/4"]:
                z = sorted([z_at0("RCS1/2", t1), z_at0("collider", t2), z_at0("RCS3/4", t3)])
                g = min(z[1] - z[0], z[2] - z[1])
                if best is None or g > best[0]:
                    best = (g, t1, t2, t3)
    if best is None:
        return None
    g, t1, t2, t3 = best
    return dict(min_plane_gap_m=g, cover_min_m=cover_min, fence_min_m=fence_min,
                rings={nm: dict(theta_mrad=th * 1e3, **{k: v for k, v in ring_state(nm, th).items()
                                                        if k in ("z0", "depth", "cover")},
                                inside_fence_m=FENCE_Y - ring_state(nm, th)["y_emerge"])
                       for nm, th in (("RCS1/2", t1), ("collider", t2), ("RCS3/4", t3))})

SPREAD = [best_spread(6.0, 0.0), best_spread(8.0, 150.0)]

# ---------------------------------------------------------------- transfer lines
# mu+ runs north->south on every east straight (mu- the other way).  A transfer takes the
# beam from a point on ring A's east straight to a point on ring B's east straight, running
# in the beam's own direction; the line is a chord with one bend at each end.
C_LIGHT = 0.2997925                                   # GeV/(T m) per unit charge
STAGE_P = {"RCS2 to RCS3": 750.0, "RCS4 to collider": 5000.0, "RCS1 to RCS2": 314.0, "RCS3 to RCS4": 1500.0}
SEPTUM_CLR = 40.0                                     # m of straight kept upstream of a septum
KICKER_CLR = 80.0                                     # m of straight kept downstream of an injection septum
L_MIN = 120.0                                         # m, shortest line we will call a transfer line

def straight_span(cfg, nm):
    yc = cfg[nm]["yc"]; return yc - RING[nm]["Ls"] / 2, yc + RING[nm]["Ls"] / 2

def design_transfer(cfg, a, b, p_gev, sign, place=None):
    """sign=+1: mu+, which runs north->south on every east straight; sign=-1: mu-, the other way.
    Extraction septum in ring A's east straight, injection septum in ring B's, the line a chord
    between them with one bend at each end.  Slopes are dz/ds ALONG THE DIRECTION OF TRAVEL, so a
    ring tilted up-to-north by theta presents slope -theta to a south-going beam."""
    ysA, ynA = straight_span(cfg, a); ysB, ynB = straight_span(cfg, b)
    thA, thB = cfg[a]["theta"], cfg[b]["theta"]
    zA = lambda y: beam_z(cfg[a]["z0"], thA, cfg[a]["yc"], y)
    zB = lambda y: beam_z(cfg[b]["z0"], thB, cfg[b]["yc"], y)
    dirn = -1.0 if sign > 0 else 1.0                      # dy/ds along the beam
    exA = np.arange(ysA + SEPTUM_CLR, ynA - SEPTUM_CLR + 1e-9, 2.5)
    inB = np.arange(ysB + SEPTUM_CLR, ynB - SEPTUM_CLR + 1e-9, 2.5)
    best, longest = None, None
    for ya in exA:
        for yb in inB:
            if dirn * (yb - ya) <= 0:                     # the line must run downstream
                continue
            if (dirn < 0 and yb - ysB < KICKER_CLR) or (dirn > 0 and ynB - yb < KICKER_CLR):
                continue                                  # leave B's straight room for the kicker
            L = math.hypot(yb - ya, X_OFF[b] - X_OFF[a])
            dz = zB(yb) - zA(ya)
            ah = abs(math.atan2(X_OFF[b] - X_OFF[a], abs(yb - ya)))
            av = abs(dz / L - thA * dirn) + abs(thB * dirn - dz / L)
            bend = 2 * ah + av
            bl = bend * p_gev / C_LIGHT
            rec = dict(transfer="%s to %s" % (a, b), sign="mu+" if sign > 0 else "mu-",
                       p_GeV=p_gev, extract_y_m=float(ya), inject_y_m=float(yb),
                       length_m=float(L), dx_m=X_OFF[b] - X_OFF[a], dz_m=float(dz),
                       bend_horiz_mrad=2 * ah * 1e3, bend_vert_mrad=av * 1e3,
                       bend_total_mrad=bend * 1e3, BL_Tm=bl,
                       dipole_m_at_8T=bl / 8.0, dipole_m_at_1p8T=bl / 1.8)
            if longest is None or L > longest["length_m"]:
                longest = rec
            if L >= L_MIN and (best is None or bl < best["BL_Tm"]):
                best = rec
    if best is not None:
        best["feasible"] = True
        return best
    if longest is not None:
        longest["feasible"] = False
        longest["why"] = ("no room: the longest line this sign can be given is %.0f m, below the %.0f m "
                          "a septum, dogleg and matching section need" % (longest["length_m"], L_MIN))
        return longest
    return dict(transfer="%s to %s" % (a, b), sign="mu+" if sign > 0 else "mu-", p_GeV=p_gev,
                feasible=False, why="no downstream injection point exists on ring B's east straight")

def line_clearance(cfg, rec, place=None, want=8.0):
    """How soon the line can have its own bore, and how close it passes to the ring it does not
    connect.  A transfer line necessarily starts inside its own ring's tunnel: what matters is the
    distance at which it has separated by `want` metres, and whether it then fouls a third ring."""
    if not rec.get("length_m"):
        return None
    a, b = rec["transfer"].split(" to ")
    ya, yb = rec["extract_y_m"], rec["inject_y_m"]
    t = np.linspace(0, 1, 400)
    lx = X_OFF[a] + t * (X_OFF[b] - X_OFF[a])
    ly = ya + t * (yb - ya)
    za = beam_z(cfg[a]["z0"], cfg[a]["theta"], cfg[a]["yc"], ya)
    zb = beam_z(cfg[b]["z0"], cfg[b]["theta"], cfg[b]["yc"], yb)
    lz = za + t * (zb - za)
    dist = {}
    for nm in ORDER:
        rx, ry, rz = ring3d(nm, cfg[nm]["theta"], cfg[nm]["z0"], X_OFF[nm], 2500, place)
        dist[nm] = np.sqrt((lx[:, None] - rx[None, :]) ** 2 + (ly[:, None] - ry[None, :]) ** 2
                           + (lz[:, None] - rz[None, :]) ** 2).min(axis=1)
    sarc = t * rec["length_m"]
    allc = np.minimum.reduce([dist[nm] for nm in ORDER])
    free = np.where(allc >= want)[0]
    third = [nm for nm in ORDER if nm not in (a, b)]
    return dict(own_bore_after_m=float(sarc[free[0]]) if len(free) else None,
                own_bore_until_m=float(sarc[free[-1]]) if len(free) else None,
                min_to_third_ring_m=float(dist[third[0]].min()) if third else None,
                third_ring=third[0] if third else None)

def design_stacked(p_gev, dz=1.5, L=200.0):
    """RCS1/RCS2 and RCS3/RCS4 share a tunnel as a vertically stacked pair, so their transfer is
    a plain vertical dogleg of dz over L: two bends, no horizontal component."""
    bend = 2 * dz / L
    bl = bend * p_gev / C_LIGHT
    return dict(kind="stacked pair, same tunnel", p_GeV=p_gev, length_m=L, dz_m=dz,
                bend_total_mrad=bend * 1e3, BL_Tm=bl, dipole_m_at_8T=bl / 8.0,
                dipole_m_at_1p8T=bl / 1.8)

def solve_place(place):
    """Tilts, depths and clearances for one placement of the straight centres."""
    win = {nm: tilt_window(nm, place=place) for nm in ORDER}
    if place == "as built":
        th = {"RCS1/2": _U["rings"]["RCS1/2"]["theta_mrad"] * 1e-3, "collider": 15.40e-3,
              "RCS3/4": _U["rings"]["RCS3/4"]["theta_mrad"] * 1e-3}
    else:                                   # window midpoints, the same rule uiuc_chain.py uses
        th = {nm: 0.5 * (win[nm][0] + win[nm][1]) for nm in ORDER}
    cfg = {nm: ring_state(nm, th[nm], place) for nm in ORDER}
    clear = [pair_min(cfg, a, b, place=place) for a, b in PAIRS]
    tr = []
    for a, b, key in (("RCS1/2", "RCS3/4", "RCS2 to RCS3"), ("RCS3/4", "collider", "RCS4 to collider")):
        for sign in (+1, -1):
            rec = design_transfer(cfg, a, b, STAGE_P[key], sign, place)
            rec["stage"] = key
            rec["clearance_m"] = line_clearance(cfg, rec, place)
            tr.append(rec)
    return dict(place=place, windows={nm: [w[0] * 1e3, w[1] * 1e3] for nm, w in win.items()},
                rings=cfg, clearances=clear, transfers=tr)

def ring_depth_profile(cfg, nm, place=None):
    """Depth below grade along the ring's north-south reach, and the extremes."""
    yc = cfg[nm]["yc"]; r = RING[nm]
    y = np.linspace(yc - r["Ls"] / 2 - r["R"], yc + r["Ls"] / 2 + r["R"], 400)
    d = np.array([TER(v) - beam_z(cfg[nm]["z0"], cfg[nm]["theta"], yc, v) for v in y])
    return y, d

def depth_summary(cfg, place=None):
    out = {}
    for nm in ORDER:
        y, d = ring_depth_profile(cfg, nm, place)
        r = RING[nm]; reach = r["Ls"] + 2 * r["R"]
        out[nm] = dict(circumference_m=r["C"], straight_m=r["Ls"], arc_R_m=r["R"],
                       ns_reach_m=reach, tilt_mrad=cfg[nm]["theta"] * 1e3,
                       depth_at_straight_m=cfg[nm]["depth"],
                       deepest_m=float(d.max()), deepest_at_y_m=float(y[int(np.argmax(d))]),
                       shallowest_m=float(d.min()), shallowest_at_y_m=float(y[int(np.argmin(d))]),
                       depth_swing_m=float(d.max() - d.min()),
                       apex_cover_m=cfg[nm].get("cover", float(d.min())),
                       under_ceiling=bool(d.max() <= DEPTH_CAP))
    return out

SOLVED = {p: solve_place(p) for p in Y_C}
for p in SOLVED:
    SOLVED[p]["depths"] = depth_summary(SOLVED[p]["rings"], p)
BASE_CLEAR = [pair_min(BASE, a, b, place="as built") for a, b in PAIRS]
BASE_DEPTH = depth_summary(BASE, "as built")
STACKED = [dict(stage="RCS1 to RCS2", **design_stacked(STAGE_P["RCS1 to RCS2"])),
           dict(stage="RCS3 to RCS4", **design_stacked(STAGE_P["RCS3 to RCS4"]))]

# ---------------------------------------------------------------- figure 1: depth and tilt
Y_LO, Y_HI = -3150.0, 3050.0
YG = np.linspace(Y_LO, Y_HI, 700)
FENCE_XY = cc.FENCE_XY                                   # (x east, y north) metres from the IP
CT = SOLVED["as built"]["rings"]

def depth_panel(ax, cfg, dmax, show_emerge=True, lw_ring=1.1):
    ax.axhspan(-7, 0, color="#b08a55", alpha=.20, lw=0)
    ax.plot([Y_LO, Y_HI], [0, 0], color="#7a5c33", lw=1.0, zorder=3)
    for nm in ORDER:
        c = cfg[nm]; r = RING[nm]; yc = c["yc"]
        y, d = ring_depth_profile(cfg, nm)
        ax.plot(y, d, color=r["col"], lw=lw_ring, alpha=.85, zorder=4)
        ysr = np.linspace(yc - r["Ls"] / 2, yc + r["Ls"] / 2, 40)
        ds = np.array([TER(v) - beam_z(c["z0"], c["theta"], yc, v) for v in ysr])
        ax.plot(ysr, ds, color=r["col"], lw=3.6, solid_capstyle="butt", zorder=6)
        ax.plot([y[-1]], [d[-1]], "o", ms=4.2, color=r["col"], mec="w", mew=0.8, zorder=7)
        if show_emerge and Y_LO < c.get("y_emerge", 1e9) < Y_HI:
            ax.plot([c["y_emerge"]], [0], "v", ms=5.5, color=r["col"], mec="w", mew=0.7, zorder=8)
    ax.axvline(FENCE_Y, color="0.45", lw=0.8, ls=(0, (5, 3)), zorder=2)
    ax.set_ylim(dmax, -9); ax.set_xlim(Y_LO, Y_HI)
    ax.set_ylabel("depth below grade (m)", fontsize=7.8)
    ax.tick_params(labelsize=7.0)

fig = plt.figure(figsize=(11.0, 7.0))
gs = fig.add_gridspec(3, 2, width_ratios=[1.0, 1.30], height_ratios=[1.0, 1.0, 0.62],
                      wspace=0.20, hspace=0.34, left=0.055, right=0.978, top=0.885, bottom=0.135)

# (a) plan, north up
axp = fig.add_subplot(gs[:, 0])
axp.fill([p[0] for p in FENCE_XY], [p[1] for p in FENCE_XY], color="0.90", ec="0.55", lw=0.8, zorder=1)
LAB = {"RCS1/2": (-1780, 40, "right"), "collider": (-3140, -560, "right"), "RCS3/4": (-4500, -1450, "right")}
for nm in ORDER:
    x, y = racetrack(nm, X_OFF[nm], 1600, "as built")
    axp.plot(x, y, color=RING[nm]["col"], lw=1.8, zorder=4)
    r = RING[nm]; yc = CT[nm]["yc"]
    axp.plot([X_OFF[nm]] * 2, [yc - r["Ls"] / 2, yc + r["Ls"] / 2], color=RING[nm]["col"],
             lw=4.4, solid_capstyle="butt", zorder=6)
    lx, ly, ha = LAB[nm]
    axp.text(lx, ly, "%s\nC = %.2f km" % (nm, r["C"] / 1000), fontsize=7.2, color=RING[nm]["col"],
             ha=ha, va="center", fontweight="bold", zorder=9,
             bbox=dict(boxstyle="round,pad=0.14", fc="w", ec="none", alpha=.85))
axp.axhline(FENCE_Y, color="0.45", lw=0.8, ls=(0, (5, 3)), zorder=2)
axp.text(-4650, FENCE_Y + 90, "Fermilab north fence", fontsize=6.7, color="0.45", ha="left", va="bottom")
axp.plot([0], [0], "*", ms=12, color="#0b6e4f", mec="w", mew=0.8, zorder=10)
axp.text(150, 60, "IP", fontsize=7.6, color="#0b6e4f", fontweight="bold", zorder=10)
axp.annotate("", (330, 2500), xytext=(330, 2050), arrowprops=dict(arrowstyle="-|>", color="0.35", lw=1.0))
axp.text(330, 2600, "N", fontsize=8, color="0.35", ha="center", fontweight="bold")
axp.annotate("all three east straights lie\nwithin 26 m of the meridian", (0, -900),
             xytext=(-1900, -2450), fontsize=6.8, color="0.3", ha="center", va="center", zorder=11,
             bbox=dict(boxstyle="round,pad=0.16", fc="w", ec="0.6", lw=0.5, alpha=.92),
             arrowprops=dict(arrowstyle="-", color="0.5", lw=0.6))
axp.set_xlim(-4780, 700); axp.set_ylim(Y_LO, Y_HI)
axp.set_xlabel("metres east of the meridian", fontsize=7.8)
axp.set_ylabel("metres north of the IP", fontsize=7.8)
axp.set_title("(a)  plan: three nested racetracks", fontsize=8.6, loc="left", pad=5)
axp.tick_params(labelsize=7.0)
axp.set_aspect("equal", adjustable="box")

# (b) co-tilted section
axb = fig.add_subplot(gs[0, 1])
depth_panel(axb, CT, 95)
axb.set_title("(b)  co-tilted chain: every ring aimed at UIUC", fontsize=8.6, loc="left", pad=5)
ANN = {"RCS1/2": (-3050, 26, "left", "#1f6f8b"), "collider": (-3050, 63, "left", "#0b6e4f"),
       "RCS3/4": (-1450, 90, "left", "#b5541c")}
for nm in ORDER:
    d = SOLVED["as built"]["depths"][nm]; tx, ty, ha, _ = ANN[nm]
    axb.text(tx, ty, "%s  %.3f mrad up-N\nstraight %.1f m deep; %.0f m of reach $\\Rightarrow$ %.0f m of swing"
             % (nm, d["tilt_mrad"], d["depth_at_straight_m"], d["ns_reach_m"], d["depth_swing_m"]),
             fontsize=6.4, color=RING[nm]["col"], ha=ha, va="center", fontweight="bold", zorder=20,
             bbox=dict(boxstyle="round,pad=0.14", fc="w", ec="none", alpha=.88))
axb.annotate("north-arc cover only %.1f m" % SOLVED["as built"]["depths"]["collider"]["apex_cover_m"],
             (1878, 5.9), xytext=(560, 14), fontsize=6.5, color="#0b6e4f", ha="right", va="center",
             fontweight="bold", zorder=21, bbox=dict(boxstyle="round,pad=0.14", fc="w", ec="#0b6e4f", lw=0.5, alpha=.92),
             arrowprops=dict(arrowstyle="-", color="#0b6e4f", lw=0.6))
# inset: the three planes are inside a 6 m band
axi = axb.inset_axes([0.625, 0.30, 0.26, 0.40])
yz = np.linspace(-250, 250, 60)
for nm in ORDER:
    c = CT[nm]
    axi.plot(yz, [TER(v) - beam_z(c["z0"], c["theta"], c["yc"], v) for v in yz], color=RING[nm]["col"], lw=1.8)
axi.set_ylim(41, 29); axi.set_xlim(-250, 250)
axi.tick_params(labelsize=5.4, length=2, pad=1)
axi.text(0.03, 0.06, "$\\times$10 zoom at the IP", transform=axi.transAxes, fontsize=5.8, color="0.4", ha="left", va="bottom")
axi.annotate("", (170, CT["RCS1/2"]["depth"] - 0.05), xytext=(170, 33.6 + 5.21),
             arrowprops=dict(arrowstyle="<->", color="0.3", lw=0.7))
axi.text(160, 36.5, "6.2 m", fontsize=5.8, color="0.25", ha="right", va="center", fontweight="bold")
for sp in axi.spines.values():
    sp.set_linewidth(0.5); sp.set_color("0.6")

# (c) corridor baseline section
axc = fig.add_subplot(gs[1, 1])
depth_panel(axc, BASE, 118, show_emerge=False)
axc.set_title("(c)  corridor baseline: RCS planes pitched down-to-north onto the 9.25 km hall", fontsize=8.6, loc="left", pad=5)
for nm, ty in (("RCS1/2", 50), ("RCS3/4", 70), ("collider", 90)):
    axc.text(-3050, ty, "%s  %+.2f mrad, %.0f m deep" % (nm, BASE_DEPTH[nm]["tilt_mrad"], BASE[nm]["depth"]),
             fontsize=6.4, color=RING[nm]["col"], ha="left", va="center", fontweight="bold", zorder=20,
             bbox=dict(boxstyle="round,pad=0.14", fc="w", ec="none", alpha=.88))

# (d) how close the tunnels come
axd = fig.add_subplot(gs[2, 1])
labs = ["RCS1/2\n& collider", "collider\n& RCS3/4", "RCS1/2\n& RCS3/4"]
ypos = np.arange(3)[::-1]
axd.axvspan(0, BORE, color="#b3261e", alpha=.13, lw=0)
axd.axvline(2 * BORE, color="0.55", lw=0.8, ls=(0, (4, 2)))
for k, (rows, mk, col, lab) in enumerate(((SOLVED["as built"]["clearances"], "o", "#b5541c", "co-tilted chain"),
                                          (BASE_CLEAR, "s", "0.45", "corridor baseline"))):
    axd.plot([r["min_sep_m"] for r in rows], ypos + (0.16 if k == 0 else -0.16), mk, ms=6.5,
             color=col, mec="w", mew=0.8, ls="", label=lab, zorder=5)
    for r, yy in zip(rows, ypos + (0.16 if k == 0 else -0.16)):
        axd.text(r["min_sep_m"] + 1.2, yy, "%.1f m" % r["min_sep_m"], fontsize=6.3, color=col,
                 va="center", ha="left", fontweight="bold")
axd.set_yticks(ypos); axd.set_yticklabels(labs, fontsize=6.4)
axd.set_xlim(0, 56); axd.set_ylim(-0.55, 2.55)
axd.set_xlabel("closest approach of the two tunnel centrelines (m)", fontsize=7.8)
axd.tick_params(labelsize=7.0)
axd.text(BORE / 2, 2.50, "bores\ntouch", fontsize=5.9, color="#b3261e", ha="center", va="top")
axd.text(2 * BORE + 0.8, -0.42, "2 bore diameters", fontsize=5.9, color="0.45", ha="left", va="bottom")
axd.legend(fontsize=6.6, frameon=False, loc="lower right", handletextpad=0.4)
axd.set_title("(d)  co-tilting squeezes the three tunnels together", fontsize=8.6, loc="left", pad=5)

H = [Line2D([], [], color="0.3", lw=3.6, label="east straight (on the meridian, aimed at UIUC)"),
     Line2D([], [], color="0.3", lw=1.1, label="the rest of the ring, projected onto the section"),
     Line2D([], [], marker="o", ls="", ms=4.2, color="0.3", label="north-arc apex: the shallowest cover"),
     Line2D([], [], marker="v", ls="", ms=5.5, color="0.3", label="where the north pencil reaches the surface"),
     Patch(color="#b08a55", alpha=.20, label="grade"),
     Line2D([], [], color="0.45", lw=0.8, ls=(0, (5, 3)), label="Fermilab north fence")]
fig.legend(handles=H, loc="lower center", ncol=3, fontsize=6.8, frameon=False,
           bbox_to_anchor=(0.5, 0.035), handlelength=2.4, columnspacing=1.6)
fig.suptitle("Depth and tilt of the three rings", fontsize=11.5, y=0.968)
fig.text(0.5, 0.930, "each ring is a planar racetrack, so its tilt turns the whole north$-$south reach into depth: "
                     "$z(y) = z_0 + \\theta\\,(y-y_c) + (y-y_c)^2/2R_E$", ha="center", va="center",
         fontsize=8.2, color="0.35")
fig.text(0.5, 0.011,
         "Deepest points %.0f m (RCS3/4's south arc), %.0f m (collider), %.0f m (RCS1/2): the 200 m ceiling never binds a ring. "
         "What binds is the shallow end $-$ %.1f m of cover over the collider's north arc $-$ and, once every ring is co-tilted to the\n"
         "same 15.4 mrad, the vertical room between them: the three planes end up inside a %.1f m band, so where two rings cross "
         "in plan the tunnels pass within %.1f m. The corridor baseline, whose planes are 20$-$40 m apart, keeps %.0f$-$%.0f m."
         % (SOLVED["as built"]["depths"]["RCS3/4"]["deepest_m"], SOLVED["as built"]["depths"]["collider"]["deepest_m"],
            SOLVED["as built"]["depths"]["RCS1/2"]["deepest_m"], SOLVED["as built"]["depths"]["collider"]["apex_cover_m"],
            max(abs(r["dz_m"]) for r in SOLVED["as built"]["clearances"]),
            min(r["min_sep_m"] for r in SOLVED["as built"]["clearances"][1:]),
            min(r["min_sep_m"] for r in BASE_CLEAR), max(r["min_sep_m"] for r in BASE_CLEAR)),
         ha="center", va="center", fontsize=6.8, color="0.35")
for ext in ("pdf", "svg"):
    fig.savefig(os.path.join(ROOT, "static", "figs", "ring_stack." + ext), bbox_inches="tight")

# ---------------------------------------------------------------- figure 2: transfer lines
def line_xyz(cfg, rec, n=200):
    a, b = rec["transfer"].split(" to ")
    ya, yb = rec["extract_y_m"], rec["inject_y_m"]
    t = np.linspace(0, 1, n)
    za = beam_z(cfg[a]["z0"], cfg[a]["theta"], cfg[a]["yc"], ya)
    zb = beam_z(cfg[b]["z0"], cfg[b]["theta"], cfg[b]["yc"], yb)
    return (X_OFF[a] + t * (X_OFF[b] - X_OFF[a]), ya + t * (yb - ya), za + t * (zb - za))

TCOL = {"mu+": "#7a2e21", "mu-": "#1c4f7a"}
fig2 = plt.figure(figsize=(10.6, 7.6))
gs2 = fig2.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.30, wspace=0.16,
                        left=0.065, right=0.982, top=0.885, bottom=0.155)
YL, YH = -1150.0, 700.0

def straight_bits(ax, cfg, vertical, place):
    """Draw the three east straights and the arcs peeling off them, in plan or in section."""
    for nm in ORDER:
        c = cfg[nm]; r = RING[nm]; yc = c["yc"]
        x, y = racetrack(nm, X_OFF[nm], 4000, place)
        z = beam_z(c["z0"], c["theta"], yc, y)
        keep = (y > YL) & (y < YH)
        if vertical:
            ax.plot(y[keep], TER(y[keep]) - z[keep], color=r["col"], lw=1.0, alpha=.5, zorder=3)
        else:
            ax.plot(y[keep], x[keep], color=r["col"], lw=1.0, alpha=.5, zorder=3)
        ys = np.linspace(yc - r["Ls"] / 2, yc + r["Ls"] / 2, 60)
        if vertical:
            ax.plot(ys, TER(ys) - beam_z(c["z0"], c["theta"], yc, ys), color=r["col"], lw=3.4,
                    solid_capstyle="butt", zorder=5)
        else:
            ax.plot(ys, np.full_like(ys, X_OFF[nm]), color=r["col"], lw=3.4, solid_capstyle="butt", zorder=5)

for col, (place, ttl) in enumerate((("as built", "(a)  as built: RCS3/4's straight sits 390 m south of the other two"),
                                    ("concentric", "(b)  concentric: the collider and RCS1/2 slid 390 m south onto RCS3/4's centre"))):
    S = SOLVED[place]; cfg = S["rings"]
    axpl = fig2.add_subplot(gs2[0, col]); axse = fig2.add_subplot(gs2[1, col], sharex=axpl)
    straight_bits(axpl, cfg, False, place)
    straight_bits(axse, cfg, True, place)
    for rec in S["transfers"]:
        c = TCOL[rec["sign"]]
        if not rec.get("feasible"):
            continue
        lx, ly, lz = line_xyz(cfg, rec)
        axpl.plot(ly, lx, color="w", lw=3.2, zorder=7, solid_capstyle="round")
        axpl.plot(ly, lx, color=c, lw=1.7, zorder=8)
        axse.plot(ly, TER(ly) - lz, color="w", lw=3.2, zorder=7, solid_capstyle="round")
        axse.plot(ly, TER(ly) - lz, color=c, lw=1.7, zorder=8)
        for ax, vv in ((axpl, lx), (axse, TER(ly) - lz)):
            ax.plot([ly[0]], [vv[0]], "o", ms=4.6, color=c, mec="w", mew=0.8, zorder=10)
            ax.plot([ly[-1]], [vv[-1]], "s", ms=4.6, color=c, mec="w", mew=0.8, zorder=10)
    axpl.set_ylim(-46, 16); axpl.set_xlim(YL, YH)
    axpl.set_ylabel("metres east of the meridian", fontsize=7.6)
    axpl.set_title(ttl, fontsize=8.4, loc="left", pad=5)
    axse.set_ylim(78, 18); axse.set_xlim(YL, YH)
    axse.set_ylabel("depth below grade (m)", fontsize=7.6)
    axse.set_xlabel("metres north of the IP", fontsize=8.0)
    axse.set_title("(%s)  the same region in section" % "cd"[col], fontsize=8.4, loc="left", pad=5)
    for ax in (axpl, axse):
        ax.tick_params(labelsize=7.0)
    plt.setp(axpl.get_xticklabels(), visible=False)
    # per-line annotation, stacked in the section panel
    rows = []
    for rec in S["transfers"]:
        if rec.get("feasible"):
            rows.append("%s, %s: %.0f m line, %+.0f m E / %+.1f m in depth, %.0f mrad of bend, "
                        "%.0f T$\\cdot$m = %.0f m of 8 T" % (rec["stage"], rec["sign"], rec["length_m"],
                        rec["dx_m"], rec["dz_m"], rec["bend_total_mrad"], rec["BL_Tm"], rec["dipole_m_at_8T"]))
        else:
            rows.append("%s, %s: NO ROOM $-$ %s" % (rec["stage"], rec["sign"],
                        "the straights overlap by too little" if "no room" in rec.get("why", "")
                        else "no downstream injection point"))
    for i, r in enumerate(rows):
        sign = "mu+" if "mu+" in r else "mu-"
        axse.text(YL + 40, 63.0 + 3.6 * i, r, fontsize=6.0, color=TCOL[sign] if "NO ROOM" not in r else "0.55",
                  ha="left", va="center", zorder=20,
                  bbox=dict(boxstyle="round,pad=0.10", fc="w", ec="none", alpha=.85))

H2 = [Line2D([], [], color=TCOL["mu+"], lw=1.7, label="$\\mu^+$ transfer line (south-going on the east straights)"),
      Line2D([], [], color=TCOL["mu-"], lw=1.7, label="$\\mu^-$ transfer line (north-going)"),
      Line2D([], [], marker="o", ls="", ms=4.6, color="0.3", label="extraction septum"),
      Line2D([], [], marker="s", ls="", ms=4.6, color="0.3", label="injection septum"),
      Line2D([], [], color="0.3", lw=3.4, label="east straight"),
      Line2D([], [], color="0.3", lw=1.0, alpha=.5, label="the arcs peeling off it")]
fig2.legend(handles=H2, loc="lower center", ncol=3, fontsize=6.9, frameon=False,
            bbox_to_anchor=(0.5, 0.055), handlelength=2.4, columnspacing=1.6)
fig2.suptitle("Transfer lines between the rings", fontsize=11.5, y=0.968)
fig2.text(0.5, 0.930, "the east straights are the only place two rings run parallel and close, so every transfer is a "
                      "translation: two opposite bends, $2d/L$ of bending for an offset $d$ over a line of length $L$",
          ha="center", va="center", fontsize=8.2, color="0.35")
fig2.text(0.5, 0.013,
          "RCS1$\\to$RCS2 and RCS3$\\to$RCS4 stay inside their own tunnel (vertically stacked rings, a 1.5 m dogleg over 200 m: "
          "%.0f and %.0f T$\\cdot$m, 2 and 9 m of 8 T) and are not drawn.\n"
          "$\\mu^+$ and $\\mu^-$ counter-rotate, so each transfer needs a mirror pair of lines. As built they are not mirrors: "
          "RCS3/4's straight is 390 m south of the other two, so for one sign of each transfer the two straights barely overlap "
          "and there is no room for a line.\nSliding the collider and RCS1/2 onto RCS3/4's centre restores the symmetry $-$ all four "
          "lines close, at roughly twice the bending each, and the collider's north-arc cover improves from %.1f m to %.1f m."
          % (STACKED[0]["BL_Tm"], STACKED[1]["BL_Tm"],
             SOLVED["as built"]["rings"]["collider"]["cover"], SOLVED["concentric"]["rings"]["collider"]["cover"]),
          ha="center", va="center", fontsize=6.8, color="0.35")
for ext in ("pdf", "svg"):
    fig2.savefig(os.path.join(ROOT, "static", "figs", "transfer_lines." + ext), bbox_inches="tight")

# ---------------------------------------------------------------- report and machine-readable output
def _clean(o):
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return float(o)
    return o

print("ring geometry and tilt windows (cover >= %.0f m over the north arc, north pencil inside the fence):" % COVER_MIN)
for nm in ORDER:
    lo, hi = WIN[nm]
    print("  %-9s C = %6.2f km  Ls = %4.0f m  R = %6.1f m  N-S reach %5.0f m  |  tilt window %.4f - %.4f mrad (%.0f urad)"
          % (nm, RING[nm]["C"] / 1000, RING[nm]["Ls"], RING[nm]["R"],
             RING[nm]["Ls"] + 2 * RING[nm]["R"], lo * 1e3, hi * 1e3, (hi - lo) * 1e6))
for lab, dep in (("corridor baseline", BASE_DEPTH), ("co-tilted, as built", SOLVED["as built"]["depths"]),
                 ("co-tilted, concentric", SOLVED["concentric"]["depths"])):
    print("\n%s:" % lab)
    for nm in ORDER:
        d = dep[nm]
        print("  %-9s tilt %+7.3f mrad | straight %5.1f m | shallowest %5.1f m at y %+6.0f | deepest %5.1f m at y %+6.0f "
              "| swing %5.1f m | under the %.0f m ceiling: %s"
              % (nm, d["tilt_mrad"], d["depth_at_straight_m"], d["shallowest_m"], d["shallowest_at_y_m"],
                 d["deepest_m"], d["deepest_at_y_m"], d["depth_swing_m"], DEPTH_CAP, d["under_ceiling"]))
print("\nclosest approach of the tunnel centrelines (E-W offsets %s):"
      % ", ".join("%s %+.0f m" % (k, v) for k, v in X_OFF.items()))
for lab, rows in (("corridor baseline", BASE_CLEAR), ("co-tilted, as built", SOLVED["as built"]["clearances"]),
                  ("co-tilted, concentric", SOLVED["concentric"]["clearances"])):
    print("  %-22s %s" % (lab, "  ".join("%s %5.2f m" % (r["pair"], r["min_sep_m"]) for r in rows)))
print("\ntransfer lines:")
for place in ("as built", "concentric"):
    print("  placement: %s" % place)
    for t in SOLVED[place]["transfers"]:
        if t.get("feasible"):
            c = t["clearance_m"] or {}
            print("    %-18s %-4s %5.0f GeV | L %5.0f m | %+5.1f m E, %+6.1f m in elevation | bend %5.1f mrad "
                  "(%.1f h + %.1f v) | %6.0f T m = %5.1f m of 8 T or %5.0f m of 1.8 T | nearest %s %.1f m"
                  % (t["stage"], t["sign"], t["p_GeV"], t["length_m"], t["dx_m"], t["dz_m"], t["bend_total_mrad"],
                     t["bend_horiz_mrad"], t["bend_vert_mrad"], t["BL_Tm"], t["dipole_m_at_8T"],
                     t["dipole_m_at_1p8T"], c.get("third_ring"), c.get("min_to_third_ring_m") or float("nan")))
        else:
            print("    %-18s %-4s %5.0f GeV | %s" % (t["stage"], t["sign"], t["p_GeV"], t.get("why")))
for st in STACKED:
    print("    %-18s both %5.0f GeV | %s: %.1f m over %.0f m, %.1f mrad, %.0f T m = %.1f m of 8 T"
          % (st["stage"], st["p_GeV"], st["kind"], st["dz_m"], st["length_m"], st["bend_total_mrad"],
             st["BL_Tm"], st["dipole_m_at_8T"]))
print("\nhow far the tilt windows can prise the three planes apart:")
for sp in SPREAD:
    if sp:
        print("  cover >= %.0f m, >= %.0f m inside the fence -> %.2f m between planes (%s)"
              % (sp["cover_min_m"], sp["fence_min_m"], sp["min_plane_gap_m"],
                 ", ".join("%s %.4f mrad" % (k, v["theta_mrad"]) for k, v in sp["rings"].items())))

json.dump(_clean(dict(
    note="depth and tilt of the three rings, and the transfer lines between them",
    geometry={nm: dict(circumference_m=RING[nm]["C"], straight_m=RING[nm]["Ls"], arc_radius_m=RING[nm]["R"],
                       ns_reach_m=RING[nm]["Ls"] + 2 * RING[nm]["R"],
                       straight_centre_y_m={p: Y_C[p][nm] for p in Y_C},
                       east_west_offset_m=X_OFF[nm]) for nm in ORDER},
    east_west_offset_rationale="chain order west to east (RCS1/2, RCS3/4, collider) so that no transfer line "
                               "crosses a third ring's straight, and RCS1/2 nests inside the collider",
    depth_ceiling_m=DEPTH_CAP, cover_min_m=COVER_MIN, bore_m=BORE,
    tilt_windows_mrad={nm: [WIN[nm][0] * 1e3, WIN[nm][1] * 1e3] for nm in ORDER},
    corridor_baseline=dict(rings={nm: dict(tilt_mrad=BASE[nm]["theta"] * 1e3, depth_m=BASE[nm]["depth"],
                                           z0_masl=BASE[nm]["z0"]) for nm in ORDER},
                           depths=BASE_DEPTH, clearances=BASE_CLEAR),
    co_tilted={p: dict(rings={nm: dict(tilt_mrad=S["rings"][nm]["theta"] * 1e3, z0_masl=S["rings"][nm]["z0"],
                                       depth_at_straight_m=S["rings"][nm]["depth"],
                                       apex_cover_m=S["rings"][nm]["cover"],
                                       inside_fence_m=FENCE_Y - S["rings"][nm]["y_emerge"]) for nm in ORDER},
                       depths=S["depths"], clearances=S["clearances"], transfers=S["transfers"])
               for p, S in SOLVED.items()},
    stacked_pairs=STACKED, plane_spread_options=SPREAD),
), open(os.path.join(ROOT, "static", "geo", "ring_stack.json"), "w"), indent=1)
print("\nwrote static/figs/ring_stack.{svg,pdf}, static/figs/transfer_lines.{svg,pdf}, static/geo/ring_stack.json")
