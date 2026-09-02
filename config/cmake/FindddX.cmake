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

set(_pkg "ddX")
set(_url "https://github.com/ddsolvation/ddX")
set(_rev "c2adde3")

if(NOT DEFINED "${_pkg}_FIND_METHOD")
  if(DEFINED "${PROJECT_NAME}-dependency-method")
    set("${_pkg}_FIND_METHOD" "${${PROJECT_NAME}-dependency-method}")
  else()
    set("${_pkg}_FIND_METHOD" "cmake" "subproject" "fetch")
  endif()
  set("_${_pkg}_FIND_METHOD")
endif()

foreach(method IN ITEMS ${${_pkg}_FIND_METHOD})
  if(TARGET "ddx::ddx")
    break()
  endif()

  if("${method}" STREQUAL "cmake")
    message(STATUS "ddX: Find installed package")
    find_package("ddX" CONFIG QUIET)
    if(ddX_FOUND AND TARGET "ddx::ddx")
      message(STATUS "ddX: Found installed package")
      break()
    endif()
  endif()

  if("${method}" STREQUAL "subproject")
    set(DDX_SOURCE_DIR "${PROJECT_SOURCE_DIR}/subprojects/ddx")
    set(DDX_BINARY_DIR "${PROJECT_BINARY_DIR}/subprojects/ddx")
    if(EXISTS "${DDX_SOURCE_DIR}/CMakeLists.txt")
      message(STATUS "Include ddX from subprojects/ddx")
      set(DDX_BUILD_EXAMPLES OFF)
      if(NOT TBLITE_WITH_TESTS)
        set(BUILD_TESTING OFF)
      endif()
      add_subdirectory("${DDX_SOURCE_DIR}" "${DDX_BINARY_DIR}")
      unset(DDX_BUILD_EXAMPLES)
      unset(BUILD_TESTING)
      break()
    endif()
  endif()

  if("${method}" STREQUAL "fetch")
    message(STATUS "Retrieving ddX from ${_url}")
    include(FetchContent)
    FetchContent_Declare(
      "ddx"
      GIT_REPOSITORY "${_url}"
      GIT_TAG "${_rev}"
    )
    set(DDX_BUILD_EXAMPLES OFF)
    if(NOT TBLITE_WITH_TESTS)
      set(BUILD_TESTING OFF)
    endif()
    FetchContent_MakeAvailable("ddx")
    unset(DDX_BUILD_EXAMPLES)
    unset(BUILD_TESTING)
    break()
  endif()
endforeach()

if(NOT TARGET "ddx::ddx")
  message(FATAL_ERROR "Could not find dependency ddX")
endif()
set(ddX_FOUND TRUE)

if(DEFINED "_${_pkg}_FIND_METHOD")
  unset("${_pkg}_FIND_METHOD")
  unset("_${_pkg}_FIND_METHOD")
endif()
unset(_pkg)
unset(_url)
unset(_rev)
