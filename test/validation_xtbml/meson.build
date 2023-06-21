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

validation_tester = find_program(files('tester.py'))

run_element = [
  'ethene_xyz',
  'oh-',
  'ch3_uhf_1'
]

foreach t : run_element
    test(t,
      validation_tester,
      args: [
        tblite_exe,
        t,
        files(t/'ml_feature.csv'),
        #'--method', m,
        #'--grad', '@0@-@1@.txt'.format(m, t),
        t.contains('+') ? ['--charge', '+1'] : [],
        t.contains('-') ? ['--charge', '-1'] : [],
        t.contains('uhf_1') ? ['--spin', '1'] : [],
        t.contains('uhf_2') ? ['--spin', '2'] : [],
      ],
      suite: 'validation',
    )
endforeach
