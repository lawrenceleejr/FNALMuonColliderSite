#!/usr/bin/env python3
"""Near-field terrain at 10 m: USGS NED (opentopodata 'ned10m') along the same
144 bearings as the GEBCO cache, 0.1-6 km from the IP every 100 m, so the
up-going exits (0.5-5 km) are judged against real local relief rather than
GEBCO's 460 m grid.  Cached to data/ned_near.json, written after every chunk."""
import json, math, os, time, urllib.request
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
R_E = 6371.0088
IP_LAT, IP_LON = 41.8443, -88.222972
OUT = os.path.join(ROOT, "data", "ned_near.json")
def fwd(az_deg, D):
    p1, l1 = math.radians(IP_LAT), math.radians(IP_LON)
    a, dl = math.radians(az_deg), D / R_E
    p2 = math.asin(math.sin(p1) * math.cos(dl) + math.cos(p1) * math.sin(dl) * math.cos(a))
    l2 = l1 + math.atan2(math.sin(a) * math.sin(dl) * math.cos(p1), math.cos(dl) - math.sin(p1) * math.sin(p2))
    return round(math.degrees(p2), 5), round(math.degrees(l2), 5)
RADII = [round(0.1 * k, 1) for k in range(1, 61)]
BEARINGS = [2.5 * k for k in range(144)]
C = json.load(open(OUT)) if os.path.exists(OUT) else dict(dataset="USGS NED 10 m via api.opentopodata.org/v1/ned10m",
                                                          ip=[IP_LAT, IP_LON], radii_km=RADII, bearings_deg=BEARINGS, ip_elev_m=None, elev={})
need = []
if C["ip_elev_m"] is None:
    need.append(("ip", 0, IP_LAT, IP_LON))
for az in BEARINGS:
    row = C["elev"].setdefault("%.1f" % az, [None] * len(RADII))
    for i, D in enumerate(RADII):
        if row[i] is None:
            la, lo = fwd(az, D); need.append(("%.1f" % az, i, la, lo))
print("points needed:", len(need), flush=True)
def fetch(chunk):
    loc = "|".join("%.5f,%.5f" % (la, lo) for _, _, la, lo in chunk)
    u = "https://api.opentopodata.org/v1/ned10m?locations=" + loc
    for a in range(6):
        try:
            r = json.load(urllib.request.urlopen(u, timeout=90))
            return [x["elevation"] for x in r["results"]]
        except Exception as e:
            print("  retry", a, repr(e)[:100], flush=True); time.sleep(3.0 * (a + 1))
    return None
t0 = time.time()
for i in range(0, len(need), 95):
    ch = need[i:i + 95]; vals = fetch(ch)
    if vals is None:
        print("chunk %d FAILED" % (i // 95), flush=True); continue
    for (key, idx, la, lo), v in zip(ch, vals):
        if key == "ip": C["ip_elev_m"] = v
        else: C["elev"][key][idx] = v
    json.dump(C, open(OUT, "w"))
    print("chunk %d/%d done (%.0f s)" % (i // 95 + 1, (len(need) + 94) // 95, time.time() - t0), flush=True)
    time.sleep(1.2)
missing = sum(1 for k in C["elev"] for v in C["elev"][k] if v is None)
print("finished; missing:", missing, "ip elev:", C["ip_elev_m"], flush=True)
