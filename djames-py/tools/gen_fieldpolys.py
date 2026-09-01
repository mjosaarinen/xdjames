"""Precompute the canonical field polynomial for every parameter set.

The search in ff.find_irreducible is deterministic, so this only caches work:
it writes djames/fieldpoly.json, which the package loads at import.  Delete the
file and everything still runs, just slower on first use.
"""
import json, os, sys, time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from djames.ff import find_irreducible, ExtField

# (q, n) over every James / D-James parameter set in Tables 3-6, the toy sets,
# and the n=48 reference/demo field.
PAIRS = [
    (2, 189), (2, 283), (2, 390), (2, 578),
    (4, 105), (4, 149), (4, 223), (4, 309),
    (5, 94), (5, 132), (5, 206), (5, 274),
    (13, 77), (13, 91), (13, 171), (13, 192),
    (23, 74), (23, 78), (23, 163), (23, 167),
    (2, 32), (2, 48), (4, 24), (5, 21), (13, 18), (23, 16),
]

def main():
    path = os.path.join(os.path.dirname(__file__), "..", "djames", "fieldpoly.json")
    out = {}
    if os.path.exists(path):
        out = json.load(open(path))
    for q, n in PAIRS:
        key = "%d,%d" % (q, n)
        if key in out:
            continue
        t = time.time()
        cs = find_irreducible(q, n)
        out[key] = cs
        json.dump(out, open(path, "w"), indent=0, sort_keys=True)
        print("q=%-3d n=%-4d  %6.1fs  nonzero terms: %s"
              % (q, n, time.time() - t,
                 [i for i, c in enumerate(cs) if c]), flush=True)
    # Re-verify everything we shipped.
    for k, cs in sorted(out.items()):
        q, n = (int(x) for x in k.split(","))
        assert ExtField(q, n, fpoly=cs, verify=True)
    print("all %d field polynomials verified irreducible" % len(out))

if __name__ == "__main__":
    main()
