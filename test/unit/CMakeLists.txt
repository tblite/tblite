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

# Unit testing
set(
  tests
  "cgto-ortho"
  "coulomb-charge"
  "coulomb-multipole"
  "coulomb-thirdorder"
  "double-dictionary"
  "fit"
  "gfn1-xtb"
  "gfn2-xtb"
  "hamiltonian"
  "halogen"
  "integral-overlap"
  "integral-multipole"
  "ipea1-xtb"
  "ncoord"
  "npy"
  "npz"
  "post-processing"
  "repulsion"
  "slater-expansion"
  "spin"
  "solvation-born"
  "solvation-cpcm"
  "solvation-surface"
  "solvation-cds"
  "solvation-shift"
  "tagged-io"
  "wavefunction-restart"
  "xtb-param"
  "xtb-external"
  "ceh"
  "xtbml"
)
set(
  test-srcs
  "main.f90"
)
foreach(t IN LISTS tests)
  string(MAKE_C_IDENTIFIER ${t} t) 
  list(APPEND test-srcs "test_${t}.f90")
endforeach()

add_executable(
  "${PROJECT_NAME}-tester"
  "${test-srcs}"
)
target_link_libraries(
  "${PROJECT_NAME}-tester"
  PRIVATE
  "${PROJECT_NAME}-lib"
  "mstore::mstore"
)

foreach(t IN LISTS tests)
  add_test("${PROJECT_NAME}/${t}" "${PROJECT_NAME}-tester" "${t}")
endforeach()
