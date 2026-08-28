"""Emit kat/*.rsp: the same vectors in a line-based format.

The Rust implementation reads these with `include_str!` and a three-line
parser, which keeps it free of a JSON dependency while both implementations
are pinned to one source of truth.
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
KATDIR = os.path.join(os.path.dirname(__file__), "..", "kat")


def emit(vectors):
    out = []
    for i, v in enumerate(vectors):
        out.append("count = %d" % i)
        out.append("name = %s" % v["name"])
        out.append("seed = %s" % v["seed"])
        out.append("pk_len = %d" % v["pk_len"])
        out.append("pk_sha3_256 = %s" % v["pk_sha3_256"])
        for s in v["signatures"]:
            out.append("msg = %s" % s["msg"])
            out.append("sig = %s" % s["sig"])
        out.append("")
    return "\n".join(out) + "\n"


def main():
    for fn in sorted(os.listdir(KATDIR)):
        if not fn.endswith(".json"):
            continue
        with open(os.path.join(KATDIR, fn)) as f:
            vectors = json.load(f)
        path = os.path.join(KATDIR, fn[:-5] + ".rsp")
        with open(path, "w") as f:
            f.write(emit(vectors))
        print("wrote %s (%d sets)" % (os.path.relpath(path), len(vectors)))


if __name__ == "__main__":
    main()
