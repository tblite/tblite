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
    import sys
    from test import test_dir
except ImportError:
    exit(77)

# if len(sys.argv) < 4:
#    raise RuntimeError("Requires at least four arguments")

thr = 5.0e-7
prog = sys.argv[1]
folder = sys.argv[2]
ref_file = sys.argv[3]
args = sys.argv[4:]
# if "xyz" in folder:
test_dir(folder, ref_file, prog, args)
