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

!> @dir tblite/solvation
!> Contains the implementation of the implicit solvation models.

!> @file tblite/solvation.f90
!> Provides reexports of the implict solvation model related implementations.

!> Proxy module for implicit solvation models.
module tblite_solvation
   use mctc_env, only : error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_solvation_alpb, only : alpb_solvation, new_alpb, alpb_input
   use tblite_solvation_cpcm, only : cpcm_solvation, new_cpcm, cpcm_input
   use tblite_solvation_data, only : solvent_data, get_solvent_data
   use tblite_solvation_input, only : solvation_input
   use tblite_solvation_type, only : solvation_type
   implicit none
   private

   public :: alpb_solvation, new_alpb, alpb_input
   public :: cpcm_solvation, new_cpcm, cpcm_input
   public :: solvent_data, get_solvent_data
   public :: solvation_input, new_solvation, solvation_type


contains

!> Create new solvation model from input data
subroutine new_solvation(solv, mol, input, error)
   !> Instance of the solvation model
   class(solvation_type), allocatable, intent(out) :: solv
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input data
   type(solvation_input), intent(in) :: input
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   if (allocated(input%alpb)) then
      solv = alpb_solvation(mol, input%alpb)
      return
   end if

   if (allocated(input%cpcm)) then
      solv = cpcm_solvation(mol, input%cpcm)
      return
   end if

   call fatal_error(error, "Unknown solvation model")
end subroutine new_solvation

end module tblite_solvation
