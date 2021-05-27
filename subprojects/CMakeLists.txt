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

# Include the mctc-lib project
if(NOT TARGET mctc-lib)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mctc-lib/CMakeLists.txt)
    add_subdirectory("mctc-lib")
  else()
    set("mctc-lib-url" "https://github.com/grimme-lab/mctc-lib")
    message(STATUS "Retrieving mctc-lib from ${mctc-lib-url}")
    include(FetchContent)
    FetchContent_Declare(
      "mctc-lib"
      GIT_REPOSITORY "${mctc-lib-url}"
      GIT_TAG "HEAD"
    )
    FetchContent_MakeAvailable("mctc-lib")
  endif()
endif()

list(
  APPEND lib-deps
  "mctc-lib"
)

# Include the mstore project
if(NOT TARGET mstore)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mstore/CMakeLists.txt)
    add_subdirectory("mstore")
  else()
    set("mstore-url" "https://github.com/grimme-lab/mstore")
    message(STATUS "Retrieving mstore from ${mstore-url}")
    include(FetchContent)
    FetchContent_Declare(
      "mstore"
      GIT_REPOSITORY "${mstore-url}"
      GIT_TAG "HEAD"
    )
    FetchContent_MakeAvailable("mstore")
  endif()
endif()

# Include the multicharge project
if(NOT TARGET multicharge)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/multicharge/CMakeLists.txt)
    add_subdirectory("multicharge")
  else()
    set("multicharge-url" "https://github.com/grimme-lab/multicharge")
    message(STATUS "Retrieving multicharge from ${multicharge-url}")
    include(FetchContent)
    FetchContent_Declare(
      "multicharge"
      GIT_REPOSITORY "${multicharge-url}"
      GIT_TAG "HEAD"
    )
    FetchContent_MakeAvailable("multicharge")
  endif()
endif()

list(
  APPEND lib-deps
  "multicharge"
)

# Include the dftd4 project
if(NOT TARGET dftd4)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/dftd4/CMakeLists.txt)
    add_subdirectory("dftd4")
  else()
    set("dftd4-url" "https://github.com/dftd4/dftd4")
    message(STATUS "Retrieving dftd4 from ${dftd4-url}")
    include(FetchContent)
    FetchContent_Declare(
      "dftd4"
      GIT_REPOSITORY "${dftd4-url}"
      GIT_TAG "HEAD"
    )
    FetchContent_MakeAvailable("dftd4")
  endif()
endif()

list(
  APPEND lib-deps
  "dftd4"
)

# Include the s-dftd3 project
if(NOT TARGET s-dftd3)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/s-dftd3/CMakeLists.txt)
    add_subdirectory("s-dftd3")
  else()
    set("s-dftd3-url" "https://github.com/awvwgk/simple-dftd3")
    message(STATUS "Retrieving s-dftd3 from ${s-dftd3-url}")
    include(FetchContent)
    FetchContent_Declare(
      "s-dftd3"
      GIT_REPOSITORY "${s-dftd3-url}"
      GIT_TAG "HEAD"
    )
    FetchContent_MakeAvailable("s-dftd3")
  endif()
endif()

list(
  APPEND lib-deps
  "s-dftd3"
)

set(lib-deps "${lib-deps}" PARENT_SCOPE)
