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

exe_args = []
if fc_id == 'gcc'
  exe_args += '-fall-intrinsics'
endif

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
    'features.F90',
    'main.f90',
  ),
  dependencies: tblite_dep,
  install: install,
  fortran_args: exe_args,
  link_language: 'fortran',
)

test('help', tblite_exe, args: '--help', suite: 'app')
test('version', tblite_exe, args: '--version', suite: 'app')
test('noargs', tblite_exe, should_fail: true, suite: 'app')
test('fit-help', tblite_exe, args: ['fit', '--help'], suite: 'app')
test('param-help', tblite_exe, args: ['param', '--help'], suite: 'app')
test('run-help', tblite_exe, args: ['--color', 'auto', 'run', '--help'], suite: 'app')
test('tagdiff-help', tblite_exe, args: ['tagdiff', '--help'], suite: 'app')

test(
  'tagdiff-example',
  tblite_exe,
  args: [
    'tagdiff',
    files('test'/'results.tag'),
    files('test'/'reference.tag'),
  ],
  suite: 'app',
)

test(
  'fit-example',
  tblite_exe,
  args: [
    'fit', '--dry-run',
    files('test'/'simple-xtb.toml'),
    files('test'/'fit-input.toml'),
  ],
  suite: 'app',
)

foreach m : ['gfn2', 'gfn1', 'ipea1']
  test(
    'param-@0@'.format(m),
    tblite_exe,
    args: [
      'param',
      '--method', m,
      '--output', '@0@-xtb.toml'.format(m),
    ],
    suite: 'app',
  )
endforeach
