#!/usr/bin/env python3
"""
Minimal Python wrapper for testing the tblite command line interface.

The wrapper will assume a specific order in the arguments rather than
providing a generic command line interface by itself since it is
supposed to be used by meson for testing purposes only.
"""

try:
    import json
    import os
    import subprocess
    import sys
except ImportError:
    exit(77)

if len(sys.argv) < 4:
    raise RuntimeError("Requires at least four arguments")


class approx:
    def __init__(self, value, abs):
        self.value = value
        self.abs = abs
        self.log = []

    def __eq__(self, other):
        def compare(a, b, ctx):
            if isinstance(a, list) and isinstance(b, list):
                return all(compare(x, y, f"{ctx}[{idx}]") for idx, (x, y) in enumerate(zip(a, b)))

            if isinstance(a, dict) and isinstance(b, dict):
                try:
                    return all(compare(a[k], b[k], f"{ctx}[{k}]") for k in a.keys())
                except KeyError as e:
                    self.log.append(f"{ctx}: Missing key {e} in dictionary")
                    return False

            if isinstance(a, (int, float)) and isinstance(b, (int, float)):
                stat = abs(a - b) < self.abs
                if not stat:
                    self.log.append(f"{ctx}: {a} != {b} (diff: {abs(a - b)})")
                return stat

            stat = a == b
            if not stat:
                self.log.append(f"{ctx}: {a} != {b}")
            return stat

        stat = compare(self.value, other, "")
        if not stat:
            print("\n".join(self.log))

        return stat


def read_molden(fname):
    checked = ("Cell", "Atoms", "GTO")
    out, sec, buf = [], None, []

    def val(tok):
        try:
            return float(tok.replace("D", "E").replace("d", "e"))
        except ValueError:
            return tok

    # Store the previously collected section in out
    def flush():
        if sec is None:
            return

        if sec in checked:
            # Keep all lines if the full section should be checked
            out.append([sec, len(buf), [[val(t) for t in line.split()] for line in buf]])

        elif sec == "MO":
            # In MO section skip the MO coefficients due to phase dependence
            lines = []
            for line in buf:
                text = line.strip()
                if (
                    text.startswith("Sym=")
                    or text.startswith("Ene=")
                    or text.startswith("Spin=")
                    or text.startswith("Occup=")
                ):
                    lines.append([val(t) for t in text.split()])

            out.append([sec, len(buf), len(lines), lines])

        else:
            # Keep only the section name
            out.append([sec, len(buf)])

    with open(fname) as f:
        for line in f.read().splitlines():
            # Detect Molden section headers
            m = re.match(r"^\s*\[([^\]]+)\]", line)

            if m:
                # Store the previous section before starting a new one
                flush()

                # Use the section name without modifiers
                name = m.group(1)
                sec = name if name == "Molden Format" else name.split()[0]
                buf = [line]

            elif sec is not None:
                # Normal line inside the current section
                buf.append(line)

    # Store the final section
    flush()

    return out

thr = 1.0e-6

prog = sys.argv[1]
outp = sys.argv[2]
resp = sys.argv[3]

molden_ref = sys.argv[4] if len(sys.argv) > 4 else None
molden_out = os.path.basename(molden_ref) if molden_ref else None

with open(resp) as fd:
    wdir = os.path.dirname(fd.name)

    replacements = {
        "$ORIGIN": wdir,
        "$JSON": os.path.basename(outp),
    }

    if molden_ref is not None:
        replacements["$MOLDEN"] = molden_out

    args = []
    for arg in fd.read().strip().split("\n"):
        for key, value in replacements.items():
            arg = arg.replace(key, value)
        args.append(arg)

stat = subprocess.call(
    [prog, *args, "--json", os.path.basename(outp)],
    shell=False,
    stdin=None,
    stderr=subprocess.STDOUT,
)

if stat != 0:
    raise RuntimeError("Calculation failed")

with open(outp) as f:
    ref = json.load(f)
    del ref["version"]

with open(os.path.basename(outp)) as f:
    res = json.load(f)

assert approx(ref, abs=thr) == res

if molden_ref is not None:
    ref_molden = read_molden(molden_ref)
    res_molden = read_molden(molden_out)

    assert approx(ref_molden, abs=thr) == res_molden

