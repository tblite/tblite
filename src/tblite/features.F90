! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

#ifndef TBLITE_HAS_HDF5
#define TBLITE_HAS_HDF5 0
#endif
#ifndef TBLITE_HAS_TREXIO
#define TBLITE_HAS_TREXIO 0
#endif

!> @file tblite/features.f90
!> Provides version and feature information

!> Interfaces to query the version information of the library.
module tblite_features
   implicit none
   private

   public :: get_tblite_feature
   public :: tblite_use_hdf5
   public :: tblite_use_trexio


   !> Logical flag indicating if HDF5 support is available
   logical, parameter :: tblite_use_hdf5 = TBLITE_HAS_HDF5 /= 0
   !> Logical flag indicating if TREXIO support is available
   logical, parameter :: tblite_use_trexio = TBLITE_HAS_TREXIO /= 0


contains


!> Getter function to retrieve tblite version
pure function get_tblite_feature(flag) result(use_feature)
   character(len=*), intent(in) :: flag
   logical :: use_feature

   select case(flag)
   case("hdf5")
      use_feature = tblite_use_hdf5
   case("trexio")
      use_feature = tblite_use_trexio
   case default
      use_feature = .false.
   end select

end function get_tblite_feature


end module tblite_features
