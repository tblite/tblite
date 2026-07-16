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

!> @file tblite/api/features.f90
!> Provides a stable feature query with #tblite_get_feature.

!> API export for the feature information of the library
module tblite_api_features
   use, intrinsic :: iso_c_binding
   use tblite_api_utils, only : c_f_character
   use tblite_features, only : get_tblite_feature
   implicit none
   private

   public :: get_feature_api, namespace

   character(len=*), parameter :: namespace = "tblite_"

contains


!> Query whether a specific feature is available
function get_feature_api(flag) result(available) &
      & bind(C, name=namespace//"get_feature")
   character(kind=c_char), intent(in) :: flag(*)
   integer(c_int) :: available

   character(len=:), allocatable :: flag_str
   logical :: has_flag

   call c_f_character(flag, flag_str)

   has_flag = get_tblite_feature(trim(flag_str))
   available = merge(1_c_int, 0_c_int, has_flag)
end function get_feature_api


end module tblite_api_features
