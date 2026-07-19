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

!> @file tblite/param/post_processing/molmom.f90
!> Provides record for molecular multipole moments post-processing
module tblite_param_post_processing_molmom
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_param_post_processing_type, only : post_processing_record
   use tblite_toml, only : toml_table, get_value, set_value, add_table
   implicit none
   private

   public :: count

   character(len=*), parameter :: k_dipm = "dipole", k_qp = "quadrupole", &
      & k_molmom = "molecular-multipole"

   !> Record specifying the molecular multipole moments post-processing
   type, public, extends(post_processing_record) :: molmom_record
      !> Compute dipole moments
      logical :: moldipm = .false.
      !> Compute quadrupole moments
      logical :: molqp = .false.
   contains
      !> Read parametrization data from TOML data structure
      procedure :: load_from_toml
      !> Write parametrization data to TOML data structure
      procedure :: dump_to_toml
      !> Populate parametrization record with default values
      procedure :: populate_default_param
   end type molmom_record

   !> Masking for the molecular multipole moments post-processing
   type, public :: molmom_mask
   end type molmom_mask

   interface count
      module procedure :: count_mask
   end interface count

contains

subroutine populate_default_param(param)
   !> Multipole moment post-processing record
   class(molmom_record), intent(inout) :: param

   param%moldipm = .true.
   param%molqp= .true.

end subroutine populate_default_param

!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   !> Instance of the parametrization data
   class(molmom_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child
   integer :: stat

   call get_value(table, k_molmom, child, requested=.false.)
   if (.not.associated(child)) return
   call get_value(child, k_dipm, self%moldipm, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot read entry for molecular dipole, boolean expected")
      return
   end if

   call get_value(child, k_qp, self%molqp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot read entry for molecular quadrupole, boolean expected")
      return
   end if

end subroutine load_from_toml


!> Write parametrization data to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(molmom_record), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child

   call add_table(table, k_molmom, child)

   call set_value(child, k_dipm, self%moldipm)
   call set_value(child, k_qp, self%molqp)

end subroutine dump_to_toml


elemental function count_mask(mask) result(ncount)
   type(molmom_mask), intent(in) :: mask
   integer :: ncount
   ncount = 0
end function count_mask

end module tblite_param_post_processing_molmom
