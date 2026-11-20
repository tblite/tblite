# SPDX-License-Identifier: LGPL-3.0-or-later
# Find ScaLAPACK matching LAPACK and provide SCALAPACK::SCALAPACK.
# SCALAPACK_LIBRARY overrides automatic selection. With MKL, BLACS may also be
# specified through SCALAPACK_BLACS_LIBRARY when MPI cannot be identified.

include(FindPackageHandleStandardArgs)

set(_scalapack_dependencies TRUE)
if(NOT TARGET MPI::MPI_Fortran)
  find_package(MPI QUIET COMPONENTS Fortran)
  if(NOT MPI_Fortran_FOUND)
    set(_scalapack_dependencies FALSE)
  endif()
endif()
if(NOT TARGET LAPACK::LAPACK)
  find_package(LAPACK QUIET)
  if(NOT LAPACK_FOUND)
    set(_scalapack_dependencies FALSE)
  endif()
endif()

set(_scalapack_names scalapack)
unset(_scalapack_blacs)
string(
  TOLOWER
    "${MPI_Fortran_LIBRARY_VERSION_STRING};${MPI_Fortran_LIBRARIES};${MPI_Fortran_INCLUDE_DIRS};${MPI_Fortran_COMPILER}"
    _scalapack_mpi)
if(_scalapack_mpi MATCHES "open[ _-]?mpi|(^|[/;])ompi")
  list(APPEND _scalapack_names scalapack-openmpi)
  set(_scalapack_blacs openmpi)
elseif(_scalapack_mpi MATCHES "mpich|hydra")
  list(APPEND _scalapack_names scalapack-mpich)
  set(_scalapack_blacs intelmpi)
elseif(_scalapack_mpi MATCHES "intel.*mpi|impi|mpiif")
  set(_scalapack_blacs intelmpi)
endif()

# Use the actual LAPACK result, including its directory and integer interface.
unset(_scalapack_mkl)
foreach(_library IN LISTS LAPACK_LIBRARIES)
  if(_library MATCHES "mkl_((intel|gf)_(ilp64|lp64)|rt)\\.")
    set(_scalapack_mkl "${_library}")
    set(_scalapack_interface "${CMAKE_MATCH_3}")
    break()
  endif()
endforeach()

set(_scalapack_required SCALAPACK_LIBRARIES _scalapack_dependencies)
if(_scalapack_mkl AND NOT SCALAPACK_LIBRARY)
  if(NOT _scalapack_interface)
    # mkl_rt defaults to LP64; its runtime interface cannot be inferred here.
    set(_scalapack_interface lp64)
    if(BLA_SIZEOF_INTEGER STREQUAL "8")
      set(_scalapack_dependencies FALSE)
    endif()
  endif()
  get_filename_component(_scalapack_mkl_dir "${_scalapack_mkl}" DIRECTORY)
  set(_scalapack_suffixes "${CMAKE_FIND_LIBRARY_SUFFIXES}")
  if(_scalapack_mkl MATCHES "\\.a$")
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
  endif()
  # Avoid retaining a previous MKL interface or installation across reconfigure.
  unset(_scalapack_mkl_library CACHE)
  find_library(
    _scalapack_mkl_library
    NAMES mkl_scalapack_${_scalapack_interface}
    PATHS "${_scalapack_mkl_dir}")
  unset(_scalapack_blacs_library CACHE)
  unset(_scalapack_blacs_library)
  if(SCALAPACK_BLACS_LIBRARY)
    set(_scalapack_blacs_library "${SCALAPACK_BLACS_LIBRARY}")
  elseif(_scalapack_blacs)
    find_library(
      _scalapack_blacs_library
      NAMES mkl_blacs_${_scalapack_blacs}_${_scalapack_interface}
      PATHS "${_scalapack_mkl_dir}")
  endif()
  set(CMAKE_FIND_LIBRARY_SUFFIXES "${_scalapack_suffixes}")
  set(SCALAPACK_LIBRARIES
      "${_scalapack_mkl_library};${_scalapack_blacs_library}")
  list(APPEND _scalapack_required _scalapack_mkl_library
       _scalapack_blacs_library)
elseif(SCALAPACK_LIBRARY)
  set(SCALAPACK_LIBRARIES "${SCALAPACK_LIBRARY}")
else()
  find_package(PkgConfig QUIET)
  if(PkgConfig_FOUND)
    pkg_search_module(PC_SCALAPACK QUIET IMPORTED_TARGET ${_scalapack_names})
  endif()
  if(TARGET PkgConfig::PC_SCALAPACK)
    set(SCALAPACK_LIBRARIES PkgConfig::PC_SCALAPACK)
  endif()
endif()

find_package_handle_standard_args(Scalapack
                                  REQUIRED_VARS ${_scalapack_required})

if(Scalapack_FOUND AND NOT TARGET SCALAPACK::SCALAPACK)
  add_library(SCALAPACK::SCALAPACK INTERFACE IMPORTED)
  set_target_properties(
    SCALAPACK::SCALAPACK
    PROPERTIES INTERFACE_LINK_LIBRARIES
               "${SCALAPACK_LIBRARIES};LAPACK::LAPACK;MPI::MPI_Fortran")
endif()

mark_as_advanced(SCALAPACK_LIBRARY)
