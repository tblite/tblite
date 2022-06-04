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

!> @file tblite/param/dispersion.f90
!> Provides model for the dispersion corrections

!> Definition of the dispersion corrections
module tblite_param_dispersion
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io_symbols, only : to_number
   use tblite_param_serde, only : serde_record
   use tblite_toml, only : toml_table, get_value, set_value, add_table
   implicit none
   private

   public :: count


   character(len=*), parameter :: k_d3 = "d3", k_d4 = "d4", k_sc = "sc", &
      & k_s6 = "s6", k_s8 = "s8", k_s9 = "s9", k_a1 = "a1", k_a2 = "a2"

   !> Parametrization record specifying the dispersion model
   type, public, extends(serde_record) :: dispersion_record
      !> Scaling for dipole-dipole (C6) interactions
      real(wp) :: s6
      !> Scaling for dipole-quadrupole (C6) interactions
      real(wp) :: s8
      !> Scaling for critical radius in rational damping function
      real(wp) :: a1
      !> Offset for critical radius in rational damping function
      real(wp) :: a2
      !> Scaling for triple-dipole (C9) interactions
      real(wp) :: s9
      !> Use DFT-D3 type dispersion
      logical :: d3
      !> Use selfconsistent DFT-D4 type dispersion
      logical :: sc
   contains
      generic :: load => load_from_array
      generic :: dump => dump_to_array
      !> Read parametrization data from TOML data structure
      procedure :: load_from_toml
      !> Write parametrization data to TOML data structure
      procedure :: dump_to_toml
      !> Read parametrization data from parameter array
      procedure, private :: load_from_array
      !> Write parametrization data to parameter array
      procedure, private :: dump_to_array
   end type


   !> Masking for the dispersion model
   type, public :: dispersion_mask
   end type dispersion_mask


   interface count
      module procedure :: count_mask
   end interface count


contains


!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   !> Instance of the parametrization data
   class(dispersion_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child
   integer :: stat

   if (.not.any([table%has_key(k_d3), table%has_key(k_d4)])) then
      call fatal_error(error, "Dispersion model not provided in dispersion table")
      return
   end if

   call get_value(table, k_d3, child, requested=.false., stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot read D3 dispersion table")
      return
   end if
   self%d3 = associated(child)
   if (.not.associated(child)) then
      call get_value(table, k_d4, child, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot read D4 dispersion table")
         return
      end if
      call get_value(child, k_sc, self%sc, .true., stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot read self-consistency for D4 dispersion")
         return
      end if
   end if

   call get_value(child, k_s6, self%s6, 1.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot read s6 parameter for dispersion")
      return
   end if
   call get_value(child, k_s8, self%s8, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot read s8 parameter for dispersion")
      return
   end if
   call get_value(child, k_a1, self%a1, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot read a1 parameter for dispersion")
      return
   end if
   call get_value(child, k_a2, self%a2, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot read a2 parameter for dispersion")
      return
   end if
   call get_value(child, k_s9, self%s9, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot read s9 parameter for dispersion")
      return
   end if
end subroutine load_from_toml


!> Write parametrization data to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(dispersion_record), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child

   call add_table(table, merge(k_d3, k_d4, self%d3), child)
   if (.not.self%d3) call set_value(child, k_sc, self%sc)

   call set_value(child, k_s6, self%s6)
   call set_value(child, k_s8, self%s8)
   call set_value(child, k_a1, self%a1)
   call set_value(child, k_a2, self%a2)
   call set_value(child, k_s9, self%s9)
end subroutine dump_to_toml


!> Read parametrization data from parameter array
subroutine load_from_array(self, array, offset, base, mask, error)
   class(dispersion_record), intent(inout) :: self
   real(wp), intent(in) :: array(:)
   integer, intent(inout) :: offset
   type(dispersion_record), intent(in) :: base
   type(dispersion_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

   select type(self)
   type is (dispersion_record)
      self = base
   end select

end subroutine load_from_array

!> Write parametrization data to parameter array
subroutine dump_to_array(self, array, offset, mask, error)
   class(dispersion_record), intent(in) :: self
   real(wp), intent(inout) :: array(:)
   integer, intent(inout) :: offset
   type(dispersion_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

end subroutine dump_to_array

elemental function count_mask(mask) result(ncount)
   type(dispersion_mask), intent(in) :: mask
   integer :: ncount
   ncount = 0
end function count_mask


end module tblite_param_dispersion
