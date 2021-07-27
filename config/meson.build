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

os = host_machine.system()
fc = meson.get_compiler('fortran')
fc_id = fc.get_id()

if fc_id == 'gcc'
  add_project_arguments(
    '-ffree-line-length-none',
    '-fbacktrace',
    language: 'fortran',
  )
elif fc_id == 'intel'
  add_project_arguments(
    '-traceback',
    language: 'fortran',
  )
elif fc_id == 'pgi' or fc_id == 'nvidia_hpc'
  add_project_arguments(
    '-Mbackslash',
    '-Mallocatable=03',
    '-traceback',
    language: 'fortran',
  )
endif

if get_option('openmp')
  omp_dep = dependency('openmp')
  lib_deps += omp_dep
endif

lapack_vendor = get_option('lapack')
if lapack_vendor == 'auto'
  if fc_id == 'intel'
    lapack_vendor = 'mkl'
  endif
endif

if lapack_vendor == 'mkl'
  mkl_dep = []
  if fc_id == 'intel'
    mkl_dep += fc.find_library('mkl_intel_lp64')
    if get_option('openmp')
      mkl_dep += fc.find_library('mkl_intel_thread')
    endif
  elif fc_id == 'gcc'
    mkl_dep += fc.find_library('mkl_gf_lp64')
    if get_option('openmp')
      mkl_dep += fc.find_library('mkl_gnu_thread')
    endif
  else
    error('MKL not supported for this compiler')
  endif
  if not get_option('openmp')
    mkl_dep += fc.find_library('mkl_tbb_thread')
  endif
  mkl_dep += fc.find_library('mkl_core')
  lib_deps += mkl_dep

elif lapack_vendor == 'mkl-rt'
  mkl_dep = fc.find_library('mkl_rt')
  lib_deps += mkl_dep

elif lapack_vendor == 'openblas'
  openblas_dep = dependency('openblas', required: false)
  if not openblas_dep.found()
    openblas_dep = fc.find_library('openblas')
  endif
  lib_deps += openblas_dep
  if not fc.links('external dsytrs; call dsytrs(); end', dependencies: openblas_dep)
    lapack_dep = dependency('lapack', required: false)
    if not lapack_dep.found()
      lapack_dep = fc.find_library('lapack')
    endif
    lib_deps += lapack_dep
  endif

else
  lapack_dep = dependency('lapack', required: false)
  if not lapack_dep.found()
    lapack_dep = fc.find_library('lapack')
  endif
  lib_deps += lapack_dep
  blas_dep = dependency('blas', required: false)
  if not blas_dep.found()
    blas_dep = fc.find_library('blas')
  endif
  lib_deps += blas_dep
endif

# Create the tool chain library as subproject
mctc_dep = dependency('mctc-lib', required: false)
if not mctc_dep.found()
  mctc_prj = subproject(
    'mctc-lib',
    version: '>=0.2',
    default_options: [
      'default_library=static',
    ],
  )
  mctc_dep = mctc_prj.get_variable('mctc_dep')

  if install
    install_data(
      mctc_prj.get_variable('mctc_lic'),
      install_dir: get_option('datadir')/'licenses'/meson.project_name()/'mctc-lib'
    )
  endif
endif
lib_deps += mctc_dep

# Create the tool chain library as subproject
multicharge_dep = dependency('multicharge', required: false)
if not multicharge_dep.found()
  multicharge_prj = subproject(
    'multicharge',
    version: '>=0.1',
    default_options: [
      'default_library=static',
    ],
  )
  multicharge_dep = multicharge_prj.get_variable('multicharge_dep')

  if install
    install_data(
      multicharge_prj.get_variable('multicharge_lic'),
      install_dir: get_option('datadir')/'licenses'/meson.project_name()/'multicharge'
    )
  endif
endif
lib_deps += multicharge_dep

# Create DFT-D4 library as subproject
dftd4_dep = dependency('dftd4', required: false)
if not dftd4_dep.found()
  dftd4_prj = subproject(
    'dftd4',
    version: '>=0.1',
    default_options: [
      'default_library=static',
    ],
  )
  dftd4_dep = dftd4_prj.get_variable('dftd4_dep')

  if install
    install_data(
      dftd4_prj.get_variable('dftd4_lic'),
      install_dir: get_option('datadir')/'licenses'/meson.project_name()/'dftd4'
    )
  endif
endif
lib_deps += dftd4_dep

# Create DFT-D3 library as subproject
sdftd3_dep = dependency('s-dftd3', required: false)
if not sdftd3_dep.found()
  sdftd3_prj = subproject(
    's-dftd3',
    version: '>=0.4',
    default_options: [
      'default_library=static',
    ],
  )
  sdftd3_dep = sdftd3_prj.get_variable('sdftd3_dep')

  if install
    install_data(
      sdftd3_prj.get_variable('sdftd3_lic'),
      install_dir: get_option('datadir')/'licenses'/meson.project_name()/'s-dftd3'
    )
  endif
endif
lib_deps += sdftd3_dep

# Create TOML Fortran as subproject
tomlf_dep = dependency('toml-f', required: false)
if not tomlf_dep.found()
  tomlf_prj = subproject(
    'toml-f',
    version: '>=0.2',
    default_options: [
      'default_library=static',
    ],
  )
  tomlf_dep = tomlf_prj.get_variable('tomlf_dep')

  if install
    install_data(
      tomlf_prj.get_variable('tomlf_lic'),
      install_dir: get_option('datadir')/'licenses'/meson.project_name()/'toml-f'
    )
  endif
endif
lib_deps += tomlf_dep
