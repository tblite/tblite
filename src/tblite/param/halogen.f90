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

module tblite_param_halogen
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_param_serde, only : serde_record
   use tblite_toml, only : toml_table, get_value, set_value, add_table
   implicit none(type, external)
   private

   public :: halogen_record

   character(len=*), parameter :: k_classical = "classical", k_damping = "damping", &
      & k_rscale = "rscale"

   type, extends(serde_record) :: halogen_record
      real(wp) :: damping
      real(wp) :: rscale
   contains
      procedure :: load_from_toml
      procedure :: dump_to_toml
   end type


contains


!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   !> Instance of the parametrization data
   class(halogen_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child
   integer :: stat

   call get_value(table, k_classical, child, requested=.false.)
   if (.not.associated(child)) then
      call fatal_error(error, "No entry for classical halogen bonding correction found")
      return
   end if

   call get_value(child, k_damping, self%damping, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for halogen bonding damping parameter")
      return
   end if

   call get_value(child, k_rscale, self%rscale, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for halogen bonding radii scaling")
      return
   end if
end subroutine load_from_toml

!> Write parametrization data to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(halogen_record), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child

   call add_table(table, k_classical, child)
   call set_value(child, k_damping, self%damping)
   call set_value(child, k_rscale, self%rscale)
end subroutine dump_to_toml

end module tblite_param_halogen
