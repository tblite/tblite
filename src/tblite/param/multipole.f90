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

!> Definition of the anisotropic second-order electrostatic contributions
module tblite_param_multipole
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_param_serde, only : serde_record
   use tblite_toml, only : toml_table, get_value, set_value, add_table
   implicit none
   private

   public :: multipole_record


   character(len=*), parameter :: k_damped = "damped", k_dmp3 = "dmp3", k_dmp5 = "dmp5", &
      & k_kexp = "kexp", k_shift = "shift", k_rmax = "rmax"

   !> Representation of the multipolar electrostatics
   type, extends(serde_record) :: multipole_record
      !> Damping exponent for quadratic terms
      real(wp) :: dmp3
      !> Damping exponent for cubic terms
      real(wp) :: dmp5
      !> Exponent for multipole radii
      real(wp) :: kexp
      !> Shift for valence CN
      real(wp) :: shift
      !> Maximum multipole radus
      real(wp) :: rmax
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
   class(multipole_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child
   integer :: stat

   call get_value(table, k_damped, child, requested=.false.)
   if (.not.associated(child)) then
      call fatal_error(error, "No entry for damped multipole electrostatic found")
      return
   end if

   call get_value(child, k_dmp3, self%dmp3, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for multipole range separation parameter")
      return
   end if
   call get_value(child, k_dmp5, self%dmp5, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for multipole range separation parameter")
      return
   end if
   call get_value(child, k_kexp, self%kexp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for multipole damping function exponent")
      return
   end if
   call get_value(child, k_shift, self%shift, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for CN multipole shift")
      return
   end if
   call get_value(child, k_rmax, self%rmax, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for maximum multipole radius")
      return
   end if
end subroutine load_from_toml


!> Write parametrization data to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(multipole_record), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child

   call add_table(table, k_damped, child)
   call set_value(child, k_dmp3, self%dmp3)
   call set_value(child, k_dmp5, self%dmp5)
   call set_value(child, k_kexp, self%kexp)
   call set_value(child, k_shift, self%shift)
   call set_value(child, k_rmax, self%rmax)
end subroutine dump_to_toml

end module tblite_param_multipole
