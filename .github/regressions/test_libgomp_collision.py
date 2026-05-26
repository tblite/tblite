"""Regression guard for #340: libgomp collision with torch <= 2.6.

Run by the `regression-libgomp-collision` job in `.github/workflows/wheel.yml`
after the manylinux wheel is built and installed alongside a torch release
known to ship a colliding bundled libgomp.

A child process imports torch first, then runs five identical GFN1-xTB
singlepoints. Acceptable outcomes are:
  (a) bit-identical energies (the manylinux wheel does not bundle libgomp,
      so torch's runtime is used and the SCF is deterministic), or
  (b) a clean ImportError mentioning libgomp / GOMP (load order makes the
      tblite extension unimportable; loud failure is the contract).

Any other outcome -- TBLiteRuntimeError, LAPACK info != 0, non-deterministic
energies -- is the silent-corruption class this regression exists to prevent.
"""
import subprocess
import sys
import textwrap

CHILD = textwrap.dedent("""
    import torch  # noqa: F401
    import numpy as np
    from tblite.interface import Calculator
    H2O = np.array([[0.0, 0.0, 0.0],
                    [1.4309, 1.1074, 0.0],
                    [-1.4309, 1.1074, 0.0]])
    energies = []
    for _ in range(5):
        calc = Calculator("GFN1-xTB", np.array([8, 1, 1]), H2O,
                          charge=0.0, uhf=0)
        calc.set("verbosity", 0)
        energies.append(float(calc.singlepoint().get("energy")))
    assert len(set(energies)) == 1, f"non-deterministic: {energies}"
    assert -10.0 < energies[0] < 0.0, energies[0]
""")


def main() -> int:
    proc = subprocess.run(
        [sys.executable, "-c", CHILD], capture_output=True, text=True
    )
    err = proc.stderr

    if proc.returncode == 0:
        print("OK: deterministic energies, torch-first import works.")
        return 0
    if "ImportError" in err and ("GOMP" in err or "libgomp" in err):
        print("OK: torch-first import fails loudly with libgomp diagnostic.")
        return 0
    raise AssertionError(
        "Outcome is neither deterministic success nor clean ImportError. "
        "Possible silent-corruption regression of #340.\n\n"
        "Subprocess stderr:\n" + err
    )


if __name__ == "__main__":
    sys.exit(main())
