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

tblite_exe = executable(
  meson.project_name(),
  sources: files(
    'argument.f90',
    'cli.f90',
    'cli_help.f90',
    'driver.f90',
    'driver_fit.f90',
    'driver_param.f90',
    'driver_run.f90',
    'driver_tagdiff.f90',
    'main.f90',
  ),
  dependencies: tblite_dep,
  install: install,
)

test('app-help', tblite_exe, args: '--help')
test('app-fit-help', tblite_exe, args: ['fit', '--help'])
test('app-param-help', tblite_exe, args: ['param', '--help'])
test('app-run-help', tblite_exe, args: ['run', '--help'])
test('app-tagdiff-help', tblite_exe, args: ['tagdiff', '--help'])
