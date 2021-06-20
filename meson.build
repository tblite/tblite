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

project(
  'tblite',
  'fortran',
  version: '0.1.0',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.55',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=both',
  ],
)
install = not (meson.is_subproject() and get_option('default_library') == 'static')

# General configuration information
lib_deps = []
inc_dirs = []
subdir('config')

# Collect source of the project
srcs = []
subdir('src')

# Library target
tblite_lib = library(
  meson.project_name(),
  sources: srcs,
  version: meson.project_version(),
  dependencies: lib_deps,
  include_directories: inc_dirs,
  install: install,
)

# Export dependency for other projects and test suite
tblite_inc = [tblite_lib.private_dir_include(), include_directories('include')]
tblite_dep = declare_dependency(
  link_with: tblite_lib,
  include_directories: tblite_inc,
  dependencies: lib_deps,
)

# Add executable targets
subdir('app')

# Package the license files
tblite_lic = files(
  'COPYING',
  'COPYING.LESSER',
)

tblite_header = files(
  'include/tblite.h',
)

if install
  # Distribute the license files in share/licenses/<name>
  install_data(
    tblite_lic,
    install_dir: get_option('datadir')/'licenses'/meson.project_name()
  )

  install_headers(
    tblite_header,
  )
  install_subdir(
    'include/tblite',
    install_dir: get_option('includedir'),
  )

  module_id = meson.project_name() / fc_id + '-' + fc.version()
  meson.add_install_script(
    find_program(files('config'/'install-mod.py')),
    get_option('includedir') / module_id,
  )

  pkg = import('pkgconfig')
  pkg.generate(
    tblite_lib,
    description: 'Light-weight tight-binding framework',
    subdirs: ['', module_id],
  )

  asciidoc = find_program('asciidoctor', required: false)
  if asciidoc.found()
    install_man(
      configure_file(
        command: [asciidoc, '-b', 'manpage', '@INPUT@', '-o', '@OUTPUT@'],
        input: files('man/tblite.1.adoc'),
        output: '@BASENAME@',
      )
    )
  endif
endif

# add the testsuite
subdir('test')

if get_option('python')
  subdir('python/tblite')
endif
