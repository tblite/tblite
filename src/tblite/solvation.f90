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
   use tblite_solvation_cds, only : cds_solvation, new_cds, cds_input
   use tblite_solvation_data, only : solvent_data, get_solvent_data
   use tblite_solvation_input, only : solvation_input
   use tblite_solvation_type, only : solvation_type
   use tblite_data_alpb, only: get_alpb_param
   use tblite_data_cds, only: get_cds_param
   implicit none
   private

   public :: alpb_solvation, new_alpb, alpb_input
   public :: cpcm_solvation, new_cpcm, cpcm_input
   public :: cds_solvation, new_cds, cds_input
   public :: solvent_data, get_solvent_data
   public :: solvation_input, new_solvation, solvation_type
   public :: new_solvation_cds

contains

!> Create new solvation model from input data
subroutine new_solvation(solv, mol, input, error, method)
   !> Instance of the solvation model
   class(solvation_type), allocatable, intent(out) :: solv
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input data
   type(solvation_input), intent(in) :: input
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Method for parameter selection
   character(len=*), optional, intent(in) :: method
   !> scratch input
   type(alpb_input), allocatable :: scratch_input

   if (allocated(input%alpb)) then
      !> xTB like ALPB/GBSA with empirical parameters 
      if (input%alpb%xtb) then
         if ( .not. present(method)) then
            call fatal_error(error, "Unkown method for solvation model parameter selection")
            return
         end if
         scratch_input = input%alpb
         scratch_input%method = method
         call get_alpb_param(scratch_input, mol, error)
         solv = alpb_solvation(mol, scratch_input)
         return
      !> ALPB/GBSA without empirical parameters
      else
         solv = alpb_solvation(mol, input%alpb)
         return
      end if
   end if

   if (allocated(input%cpcm)) then
      solv = cpcm_solvation(mol, input%cpcm)
      return
   end if

   call fatal_error(error, "Unknown solvation model")
end subroutine new_solvation

!> Create new solvation model from input data
subroutine new_solvation_cds(solv, mol, input, error, method)
   !> Instance of the solvation model
   class(solvation_type), allocatable, intent(out) :: solv
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input data
   type(solvation_input), intent(in) :: input
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Method for parameter selection
   character(len=*), optional, intent(in) :: method
   !> scratch input
   type(cds_input), allocatable :: scratch_input

   if (allocated(input%cds).and.present(method)) then
      scratch_input = input%cds
      scratch_input%method = method
      call get_cds_param(scratch_input, mol, error)
      solv = cds_solvation(mol, scratch_input)
      return
    end if

   call fatal_error(error, "Unknown cds model")
end subroutine new_solvation_cds

end module tblite_solvation
