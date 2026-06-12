@PACKAGE_INIT@

set("TBLITE_WITH_API" @TBLITE_WITH_API@)
set("TBLITE_WITH_OpenMP" @TBLITE_WITH_OpenMP@)
set("TBLITE_WITH_DDX" @TBLITE_WITH_DDX@)
set("TBLITE_WITH_TREXIO" @TBLITE_WITH_TREXIO@)
set("TBLITE_WITH_HDF5" @TBLITE_WITH_HDF5@)
set("TBLITE_USE_MCTCLIB" @TBLITE_USE_MCTCLIB@)
set("TBLITE_USE_MULTICHARGE" @TBLITE_USE_MULTICHARGE@)
set("TBLITE_USE_TOMLF" @TBLITE_USE_TOMLF@)
set("TBLITE_USE_SDFTD3" @TBLITE_USE_SDFTD3@)
set("TBLITE_USE_DFTD4" @TBLITE_USE_DFTD4@)
set("TBLITE_USE_DDX" @TBLITE_USE_DDX@)

enable_language("Fortran")
enable_language("C")

if(NOT TARGET "@PROJECT_NAME@::@PROJECT_NAME@")
  include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")

  include(CMakeFindDependencyMacro)

  if(NOT TARGET "OpenMP::OpenMP_Fortran" AND TBLITE_WITH_OpenMP)
    find_dependency("OpenMP")
  endif()

  if(NOT TARGET "LAPACK::LAPACK")
    find_dependency("LAPACK")
  endif()

  if(NOT TARGET "ddx" AND TBLITE_USE_DDX)
    find_dependency("ddX")
  endif()

  if(NOT TARGET "HDF5::HDF5" AND TBLITE_WITH_HDF5)
     find_dependency("HDF5")
  endif()

  if(NOT TARGET "TREXIO::TREXIO" AND TBLITE_WITH_TREXIO)
    find_dependency("TREXIO")
  endif()

  if(NOT TARGET "toml-f::toml-f" AND TBLITE_USE_TOMLF)
    find_dependency("toml-f")
  endif()

  if(NOT TARGET "mctc-lib::mctc-lib" AND TBLITE_USE_MCTCLIB)
    find_dependency("mctc-lib")
  endif()

  if(NOT TARGET "multicharge::multicharge" AND TBLITE_USE_MULTICHARGE)
    find_dependency("multicharge")
  endif()

  if(NOT TARGET "dftd4::dftd4" AND TBLITE_USE_DFTD4)
    find_dependency("dftd4")
  endif()

  if(NOT TARGET "s-dftd3::s-dftd3" AND TBLITE_USE_SDFTD3)
    find_dependency("s-dftd3")
  endif()
endif()
