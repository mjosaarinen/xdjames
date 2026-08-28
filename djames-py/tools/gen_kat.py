"""Generate known-answer tests.

  python3 tools/gen_kat.py toy       -> kat/toy.json     (all toy sets, seconds)
  python3 tools/gen_kat.py q2        -> kat/q2.json      (real F_2 sets, minutes)
  python3 tools/gen_kat.py keygen    -> kat/keygen.json  (pk digests, all sets)
  python3 tools/gen_kat.py <name>    -> kat/<name>.json
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from djames import kat, params  # noqa: E402

KATDIR = os.path.join(os.path.dirname(__file__), "..", "kat")

GROUPS = {
    # Everything, fast: exercises all five fields and both schemes.
    "toy": (params.names("toy-"), True, 4),
    # The paper's headline sets.  Signing over F_2 costs q^r = 4 root-findings.
    "q2": (["d-james-128-q2", "james-128-q2",
            "d-james-256-q2", "james-256-q2"], True, 2),
    # Key generation only, for the sets whose signing is impractical in pure
    # Python (q^r reaches 529 for q = 23).
    "keygen": ([n for n in params.names() if not n.startswith("toy-")], False, 0),
}


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "toy"
    if which in GROUPS:
        names, do_sign, nmsg = GROUPS[which]
    else:
        names, do_sign, nmsg = [which], True, 4
    out = []
    for name in names:
        t = time.time()
        out.append(kat.generate(name, nmsg=nmsg, sign_vectors=do_sign))
        print("%-18s %7.1fs  pk %9d B  %s"
              % (name, time.time() - t, out[-1]["pk_len"],
                 out[-1]["pk_sha3_256"][:16]), flush=True)
    path = os.path.join(KATDIR, "%s.json" % which)
    kat.save(path, out)
    print("wrote", os.path.relpath(path))


if __name__ == "__main__":
    main()
