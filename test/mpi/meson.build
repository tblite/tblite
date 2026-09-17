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

mpi_tester = executable(
  'mpi-tester',
  sources: files('main.f90'),
  dependencies: [tblite_dep, mstore_dep, mpi_dep],
  link_language: 'fortran',
)

mpirun = find_program('mpirun', 'mpiexec', required: false)
if mpirun.found()
  test(
    'mpi-singlepoint',
    mpirun,
    args: ['-n', '2', mpi_tester.full_path()],
    suite: 'mpi',
  )
endif
