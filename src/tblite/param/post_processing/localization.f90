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

!> @file tblite/param/post_processing/localization.f90
!> Provides record for orbital localization post-processing
module tblite_param_post_processing_localization
   use mctc_env, only : error_type, fatal_error
   use tblite_param_post_processing_type, only : post_processing_record
   use tblite_toml, only : toml_table, get_value, set_value, add_table
   implicit none
   private

   public :: orbital_localization_record

   character(len=*), parameter :: k_localization = "localization", k_method = "method"

   !> Default orbital localization method, used if none is specified
   character(len=*), parameter :: default_method = "foster-boys"

   !> Record specifying the orbital localization post-processing
   type, public, extends(post_processing_record) :: orbital_localization_record
      !> Name of the requested localization method
      character(len=:), allocatable :: method
   contains
      !> Read parametrization data from TOML data structure
      procedure :: load_from_toml
      !> Write parametrization data to TOML data structure
      procedure :: dump_to_toml
      !> Populate parametrization record with default values
      procedure :: populate_default_param
   end type orbital_localization_record

contains

!> Populate parametrization record with default values
subroutine populate_default_param(param, method)
   !> Orbital localization post-processing record
   class(orbital_localization_record), intent(inout) :: param
   !> Name of the localization method, defaults to Foster-Boys
   character(len=*), intent(in), optional :: method

   param%method = default_method
   if (present(method)) param%method = method

end subroutine populate_default_param


!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   !> Instance of the parametrization data
   class(orbital_localization_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child
   integer :: stat

   self%method = default_method

   call get_value(table, k_localization, child, requested=.false.)
   if (.not.associated(child)) return

   call get_value(child, k_method, self%method, default_method, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot read entry for localization method")
      return
   end if

end subroutine load_from_toml


!> Write parametrization data to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(orbital_localization_record), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child

   call add_table(table, k_localization, child)
   call set_value(child, k_method, self%method)

end subroutine dump_to_toml

end module tblite_param_post_processing_localization
