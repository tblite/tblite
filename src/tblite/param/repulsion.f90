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

!> Definition of the repulsion interactions
module tblite_param_repulsion
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_param_serde, only : serde_record
   use tblite_toml, only : toml_table, get_value, set_value, add_table
   implicit none
   private

   public :: repulsion_record


   character(len=*), parameter :: k_effective = "effective", k_kexp = "kexp", &
      & k_klight = "klight"

   !> Parametrization records describing the repulsion interactions
   type, extends(serde_record) :: repulsion_record
      real(wp) :: kexp
      real(wp) :: klight
   contains
      !> Read parametrization data from TOML data structure
      procedure :: load_from_toml
      !> Write parametrization data to TOML data structure
      procedure :: dump_to_toml
   end type


contains


!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   !> Instance of the parametrization data
   class(repulsion_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child
   integer :: stat

   call get_value(table, k_effective, child, requested=.false.)
   if (.not.associated(child)) then
      call fatal_error(error, "No entry for effective repulsion found")
      return
   end if

   call get_value(child, k_kexp, self%kexp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for repulsion exponent")
      return
   end if

   call get_value(child, k_klight, self%klight, self%kexp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for light-atom repulsion exponent")
      return
   end if
end subroutine load_from_toml


!> Write parametrization data to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(repulsion_record), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child

   call add_table(table, k_effective, child)
   call set_value(child, k_kexp, self%kexp)
   if (abs(self%kexp - self%klight) > epsilon(self%kexp)) then
      call set_value(child, k_klight, self%klight)
   end if

end subroutine dump_to_toml


end module tblite_param_repulsion
