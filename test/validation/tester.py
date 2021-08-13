#!/usr/bin/env python3
# This file is part of tblite.
# SPDX-Identifier: LGPL-3.0-or-later
#
# tblite is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tblite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with tblite.  If not, see <https://www.gnu.org/licenses/>.
"""
Minimal Python wrapper for testing the dftd4 command line interface.

The wrapper will assume a specific order in the arguments rather than
providing a generic command line interface by itself since it is
supposed to be used by meson for testing purposes only.
"""

try:
    import subprocess, sys, json, os, pytest
except ImportError:
    exit(77)

if len(sys.argv) < 4:
    raise RuntimeError("Requires at least four arguments")

thr = 1.0e-9
prog = sys.argv[1]
reference_file = sys.argv[2]
args = sys.argv[3:]

method_name = os.path.basename(reference_file)
test_name = os.path.basename(os.path.dirname(reference_file))
output_file = test_name + "-" + method_name

stat = subprocess.call(
    [prog] + args + ["--json", output_file],
    shell=False,
    stdin=None,
    stderr=subprocess.STDOUT,
)
if stat != 0:
    raise RuntimeError("Calculation failed")

with open(reference_file) as f:
    ref = json.load(f)
    del ref["version"]

with open(output_file) as f:
    res = json.load(f)

for key in ref:
    if key not in res:
        raise RuntimeError("Missing '" + key + "' entry in results")
    assert pytest.approx(res[key], abs=thr) == ref[key]
