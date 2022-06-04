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

!> @file tblite/api/container.f90
!> Provides API exports for the #tblite_container handle.

!> API export for managing interaction containers
module tblite_api_container
   use, intrinsic :: iso_c_binding
   use tblite_api_version, only : namespace
   use tblite_container, only : container_type
   use tblite_external_field, only : electric_field
   implicit none
   private

   public :: vp_container, delete_container_api


   !> Void pointer to a container instance
   type :: vp_container
      !> Actual container
      class(container_type), allocatable :: ptr
   end type vp_container

   logical, parameter :: debug = .false.

contains


function new_electric_field_api(efield) result(vcont) &
      & bind(C, name=namespace//"new_electric_field")
   real(c_double), intent(in) :: efield(3)
   type(c_ptr) :: vcont
   type(vp_container), pointer :: cont

   allocate(cont)
   cont%ptr = electric_field(efield)
   vcont = c_loc(cont)
end function new_electric_field_api


subroutine delete_container_api(vcont) &
      & bind(C, name=namespace//"delete_container")
   type(c_ptr), intent(inout) :: vcont
   type(vp_container), pointer :: cont

   if (debug) print '("[Info]", 1x, a)', "delete_container"

   if (c_associated(vcont)) then
      call c_f_pointer(vcont, cont)

      deallocate(cont)
      vcont = c_null_ptr
   end if
end subroutine delete_container_api


end module tblite_api_container
