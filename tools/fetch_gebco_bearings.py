#!/usr/bin/env python3
"""Sample GEBCO 2020 (via opentopodata) along 48 bearings from the IP, from
0.25 km to 1000 km, and cache to data/gebco_bearings.json.  Batched 95 per
call, paced, retried, and written after every chunk so progress survives."""
import json, math, os, sys, time, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
R_E = 6371.0088
IP_LAT, IP_LON = 41.8443, -88.222972
OUT = os.path.join(ROOT, "data", "gebco_bearings.json")

def fwd(az_deg, D):
    p1, l1 = math.radians(IP_LAT), math.radians(IP_LON)
    a, dl = math.radians(az_deg), D / R_E
    p2 = math.asin(math.sin(p1) * math.cos(dl) + math.cos(p1) * math.sin(dl) * math.cos(a))
    l2 = l1 + math.atan2(math.sin(a) * math.sin(dl) * math.cos(p1), math.cos(dl) - math.sin(p1) * math.sin(p2))
    return round(math.degrees(p2), 4), round(math.degrees(l2), 4)

RADII = ([0.25 * k for k in range(1, 21)] + [5 + 0.5 * k for k in range(1, 51)]
         + [30 + 2.5 * k for k in range(1, 29)] + [100 + 10 * k for k in range(1, 21)]
         + [300 + 25 * k for k in range(1, 29)])
BEARINGS = [2.5 * k for k in range(144)]

if os.path.exists(OUT):
    C = json.load(open(OUT))
else:
    C = dict(dataset="GEBCO 2020 via api.opentopodata.org", ip=[IP_LAT, IP_LON],
             radii_km=RADII, bearings_deg=BEARINGS, ip_elev_m=None, elev={})

need = []
if C["ip_elev_m"] is None:
    need.append(("ip", 0, IP_LAT, IP_LON))
for az in BEARINGS:
    key = "%.1f" % az
    row = C["elev"].setdefault(key, [None] * len(RADII))
    for i, D in enumerate(RADII):
        if row[i] is None:
            la, lo = fwd(az, D)
            need.append((key, i, la, lo))
print("points needed:", len(need), flush=True)

def fetch(chunk):
    loc = "|".join("%.4f,%.4f" % (la, lo) for _, _, la, lo in chunk)
    u = "https://api.opentopodata.org/v1/gebco2020?locations=" + loc
    for a in range(6):
        try:
            r = json.load(urllib.request.urlopen(u, timeout=90))
            return [x["elevation"] for x in r["results"]]
        except Exception as e:
            print("  retry", a, repr(e)[:80], flush=True)
            time.sleep(3.0 * (a + 1))
    return None

t0 = time.time()
for i in range(0, len(need), 95):
    ch = need[i:i + 95]
    vals = fetch(ch)
    if vals is None:
        print("chunk %d FAILED" % (i // 95), flush=True)
        continue
    for (key, idx, la, lo), v in zip(ch, vals):
        if key == "ip":
            C["ip_elev_m"] = v
        else:
            C["elev"][key][idx] = v
    json.dump(C, open(OUT, "w"))
    print("chunk %d/%d done (%.0f s)" % (i // 95 + 1, (len(need) + 94) // 95, time.time() - t0), flush=True)
    time.sleep(1.2)
C["bearings_deg"] = sorted(float(k) for k in C["elev"].keys())
json.dump(C, open(OUT, "w"))
missing = sum(1 for k in C["elev"] for v in C["elev"][k] if v is None)
print("finished; missing:", missing, "ip elev:", C["ip_elev_m"], flush=True)
