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
module tblite_param_ml_features
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io_symbols, only : to_number
   use tblite_param_serde, only : serde_record
   use tblite_toml, only : toml_table, get_value, set_value, add_table, toml_array
   use tblite_ml_features_methods, only : ml_features_method
   implicit none
   private

   public :: count

   character(len=*), parameter :: k_xtbml = "xtbml", k_xtbmlgeometry = "geometry", k_xtbmldensity = "density", &
      & k_xtbmlorbital = "orbital", k_xtbmlenergy = "energy", k_xtbmlconvolution = "convolution", k_xtbmla = "a", k_tensor = "tensorial-output"

   !> Parametrization record specifying the dispersion model
   type, public, extends(serde_record) :: ml_features_record
      !> Which ml features should be computed, 1 := xtbml, other fetaures not yet implemented
      integer :: ml_features
      !> Compute geometry-based xtbml features
      logical :: xtbml_geometry
      !> Compute density-based xtbml features
      logical :: xtbml_density
      !> Return vectorial information additional to norm of the corresponding multipole moments
      logical :: xtbml_tensor
      !> Compute orbital energy based xtbml features
      logical :: xtbml_orbital_energy
      !> Compute energy based features, necessary for partitioning weights
      logical :: xtbml_energy
      !> Compute extended feature i.e. do CN weigthed real space convolution
      logical :: xtbml_convolution
      !> Scaling for logistic function, convolution over an array of values is supported
      real(wp), allocatable :: xtbml_a(:)
      
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
   type, public :: ml_features_mask
   end type ml_features_mask


   interface count
      module procedure :: count_mask
   end interface count


contains


!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   use tblite_toml, only : len
   !> Instance of the parametrization data
   class(ml_features_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(toml_array), pointer :: array
   type(toml_table), pointer :: child
   integer :: stat, i

   if (.not.(table%has_key(k_xtbml))) then
      call fatal_error(error, "Currently only xtbml fetaures are implemented.")
      return
   end if

   call get_value(table, k_xtbml, child, requested=.false., stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot read xtbml feaure table")
      return
   end if
   self%ml_features = ml_features_method%xtbml
   if (self%ml_features == ml_features_method%xtbml) then
    call get_value(child, k_xtbmlgeometry, self%xtbml_geometry, .false., stat=stat)
    if (stat /= 0) then
        call fatal_error(error, "Cannot read entry for xtbml geometry based features, boolean expected")
        return
    end if

    call get_value(child, k_xtbmldensity, self%xtbml_density, .false., stat=stat)
    if (stat /= 0) then
        call fatal_error(error, "Cannot read entry for xtbml density based features, boolean expected")
        return
    end if
    if (self%xtbml_convolution) then 
        call get_value(child, k_tensor, self%xtbml_tensor, .false., stat=stat)
        if (stat /= 0) then
            call fatal_error(error, "Cannot read entry for xtbml tensorial-output, boolean expected")
            return
        end if 
    end if

    call get_value(child, k_xtbmlorbital, self%xtbml_orbital_energy, .false., stat=stat)
    if (stat /= 0) then
      call fatal_error(error, "Cannot read entry for xtbml orbital energy based features, boolean expected")
    return
    end if

    call get_value(child, k_xtbmlenergy, self%xtbml_energy, .false., stat=stat)
    if (stat /= 0) then
      call fatal_error(error, "Cannot read entry for xtbml geometry based features, boolean expected")
    return
    end if

    call get_value(child, k_xtbmlconvolution, self%xtbml_convolution, .false., stat=stat)
    if (stat /= 0) then
      call fatal_error(error, "Cannot read entry for xtbml convolution, boolean expected")
      return
    end if
    if (self%xtbml_convolution) then
        call get_value(child, k_xtbmla, array, stat=stat)
        if (stat /= 0) then
            call fatal_error(error, "Cannot read entry for xtbml convolution scale array, array of real values expected")
            return
        end if
        allocate(self%xtbml_a(len(array)))
        do i=1, size(self%xtbml_a)
         call get_value(array, i, self%xtbml_a(i))
      end do
   endif
   end if
   
   
end subroutine load_from_toml


!> Write parametrization data to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(ml_features_record), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child
   if (self%ml_features == ml_features_method%xtbml) then 

      call add_table(table, k_xtbml, child)

      call set_value(child, k_xtbmlgeometry, self%xtbml_geometry)
      call set_value(child, k_xtbmldensity, self%xtbml_density)
      if (self%xtbml_density) call set_value(child, k_tensor, self%xtbml_tensor)
      call set_value(child, k_xtbmlorbital, self%xtbml_orbital_energy)
      call set_value(child, k_xtbmlenergy, self%xtbml_energy)
      call set_value(child, k_xtbmlconvolution, self%xtbml_convolution)
      !if (self%xtbml_convolution) call set_value(child, k_xtbmla, self%xtbml_a)
   end if
end subroutine dump_to_toml


!> Read parametrization data from parameter array
subroutine load_from_array(self, array, offset, base, mask, error)
   class(ml_features_record), intent(inout) :: self
   real(wp), intent(in) :: array(:)
   integer, intent(inout) :: offset
   type(ml_features_record), intent(in) :: base
   type(ml_features_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

   select type(self)
   type is (ml_features_record)
      self = base
   end select

end subroutine load_from_array

!> Write parametrization data to parameter array
subroutine dump_to_array(self, array, offset, mask, error)
   class(ml_features_record), intent(in) :: self
   real(wp), intent(inout) :: array(:)
   integer, intent(inout) :: offset
   type(ml_features_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

end subroutine dump_to_array

elemental function count_mask(mask) result(ncount)
   type(ml_features_mask), intent(in) :: mask
   integer :: ncount
   ncount = 0
end function count_mask


end module tblite_param_ml_features
