# This file is part of tblite.
# SPDX-Identifier: LGPL-3.0-or-later
#
# tblite is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tblite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with tblite.  If not, see <https://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.18)

set(_tblite_python_standalone FALSE)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(_tblite_python_standalone TRUE)
  project(
    tblite-python
    LANGUAGES C
    VERSION 0.7.0
    DESCRIPTION "Python API of the tblite project"
  )
endif()

include(GNUInstallDirs)

find_package(Python3 3.7 REQUIRED COMPONENTS "Interpreter" "Development.Module")
execute_process(
  COMMAND "${Python3_EXECUTABLE}" -c "import cffi, setuptools"
  RESULT_VARIABLE _tblite_cffi_status
  ERROR_QUIET
)
if(NOT _tblite_cffi_status EQUAL 0)
  message(FATAL_ERROR
    "The Python interpreter does not provide the cffi and setuptools modules"
  )
endif()

set(_tblite_target "")
set(_tblite_include_dir "")
# An in-tree build already provides the native target. A standalone Python
# build first tries an installed CMake package and then the source bundled in
# the Python sdist.
if(TARGET tblite::tblite)
  set(_tblite_target "tblite::tblite")
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../include/tblite.h")
    set(_tblite_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/../include")
  endif()
elseif(TARGET tblite::tblite-lib)
  set(_tblite_target "tblite::tblite-lib")
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../include/tblite.h")
    set(_tblite_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/../include")
  endif()
elseif(_tblite_python_standalone)
  find_package(tblite "${PROJECT_VERSION}" CONFIG QUIET)
  if(TARGET tblite::tblite)
    set(_tblite_target "tblite::tblite")
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/subprojects/tblite/CMakeLists.txt")
    enable_language(Fortran)
    set(BUILD_SHARED_LIBS FALSE CACHE BOOL "Build shared libraries" FORCE)
    set(TBLITE_WITH_API TRUE CACHE BOOL "Build the tblite C API" FORCE)
    set(TBLITE_WITH_PYTHON FALSE CACHE BOOL "Build the tblite Python API" FORCE)
    set(TBLITE_WITH_TESTS FALSE CACHE BOOL "Build tblite tests" FORCE)
    add_subdirectory("subprojects/tblite" "subprojects/tblite")
    set(_tblite_target "tblite::tblite")
    set(_tblite_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/subprojects/tblite/include")
  endif()
endif()

if(NOT _tblite_target)
  message(FATAL_ERROR
    "Could not find tblite >= ${PROJECT_VERSION}; install its CMake package "
    "or build from a Python sdist containing subprojects/tblite"
  )
endif()

if(DEFINED TBLITE_WITH_API AND NOT TBLITE_WITH_API)
  message(FATAL_ERROR "The tblite Python extension requires the tblite C API")
endif()

if(NOT _tblite_include_dir)
  get_target_property(
    _tblite_interface_include_dirs
    "${_tblite_target}"
    INTERFACE_INCLUDE_DIRECTORIES
  )
  unset(_tblite_include_dir)
  unset(_tblite_include_dir CACHE)
  find_path(
    _tblite_include_dir
    NAMES "tblite.h"
    HINTS ${_tblite_interface_include_dirs}
  )
  if(NOT _tblite_include_dir)
    message(FATAL_ERROR "Could not locate tblite.h from ${_tblite_target}")
  endif()
endif()

set(_tblite_header "${CMAKE_CURRENT_SOURCE_DIR}/include/_tblite.h")
set(_tblite_preprocessed_header "${CMAKE_CURRENT_BINARY_DIR}/_libtblite.h")
set(_tblite_cffi_source "${CMAKE_CURRENT_BINARY_DIR}/_libtblite.c")
# CFFI cannot parse all preprocessor constructs in the public C header.
add_custom_command(
  OUTPUT "${_tblite_preprocessed_header}"
  COMMAND
    "${CMAKE_C_COMPILER}"
    "-I${_tblite_include_dir}"
    "-DTBLITE_CFFI"
    "-E"
    "${_tblite_header}"
    "-o"
    "${_tblite_preprocessed_header}"
  DEPENDS
    "${_tblite_header}"
    "${_tblite_include_dir}/tblite.h"
  VERBATIM
)

add_custom_command(
  OUTPUT "${_tblite_cffi_source}"
  COMMAND
    "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ffi-builder.py"
    "${_tblite_preprocessed_header}"
    "_libtblite"
  DEPENDS
    "${_tblite_preprocessed_header}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ffi-builder.py"
  VERBATIM
)

Python3_add_library("_libtblite" MODULE WITH_SOABI "${_tblite_cffi_source}")
target_include_directories("_libtblite" PRIVATE "${_tblite_include_dir}")
target_link_libraries("_libtblite" PRIVATE "${_tblite_target}")
if(CMAKE_Fortran_COMPILER_LOADED)
  set_property(TARGET "_libtblite" PROPERTY LINKER_LANGUAGE "Fortran")
endif()

execute_process(
  COMMAND
    "${Python3_EXECUTABLE}"
    -c
    "import sysconfig, sys; print(sysconfig.get_path('platlib', vars={'base': '', 'platbase': '', 'platlibdir': sys.argv[1]}).lstrip('/'))"
    "${CMAKE_INSTALL_LIBDIR}"
  OUTPUT_VARIABLE _tblite_python_install_root
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE _tblite_sysconfig_status
)
if(NOT _tblite_sysconfig_status EQUAL 0)
  message(FATAL_ERROR "Could not determine the Python platform-library directory")
endif()

set(
  TBLITE_PYTHON_INSTALL_DIR
  "${_tblite_python_install_root}/tblite"
  CACHE PATH
  "Directory for the tblite Python package"
)

install(
  TARGETS "_libtblite"
  LIBRARY DESTINATION "${TBLITE_PYTHON_INSTALL_DIR}" COMPONENT "python"
  RUNTIME DESTINATION "${TBLITE_PYTHON_INSTALL_DIR}" COMPONENT "python"
)

set(
  _tblite_python_sources
  "tblite/__init__.py"
  "tblite/ase.py"
  "tblite/exceptions.py"
  "tblite/interface.py"
  "tblite/library.py"
  "tblite/qcschema.py"
  "tblite/test_ase.py"
  "tblite/test_interface.py"
  "tblite/test_library.py"
  "tblite/test_qcschema.py"
  "tblite/py.typed"
)
install(
  FILES ${_tblite_python_sources}
  DESTINATION "${TBLITE_PYTHON_INSTALL_DIR}"
  COMPONENT "python"
)
