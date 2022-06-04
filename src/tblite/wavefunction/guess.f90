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

!> @file tblite/wavefunction/guess.f90
!> Provides guesses for the wavefunction

!> Implementation of the guess wavefunctions
module tblite_wavefunction_guess
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_disp_d4, only : get_eeq_charges
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_h0, only : get_occupation
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: sad_guess, eeq_guess

contains

subroutine sad_guess(mol, calc, wfn)
   type(structure_type), intent(in) :: mol
   type(xtb_calculator), intent(in) :: calc
   type(wavefunction_type), intent(inout) :: wfn

   wfn%qat(:, :) = 0.0_wp
   wfn%qat(:, 1) = mol%charge / mol%nat
   call shell_partition(mol, calc, wfn)
end subroutine sad_guess

subroutine eeq_guess(mol, calc, wfn)
   type(structure_type), intent(in) :: mol
   type(xtb_calculator), intent(in) :: calc
   type(wavefunction_type), intent(inout) :: wfn

   wfn%qat(:, :) = 0.0_wp
   call get_eeq_charges(mol, wfn%qat(:, 1))
   call shell_partition(mol, calc, wfn)
end subroutine eeq_guess

subroutine shell_partition(mol, calc, wfn)
   type(structure_type), intent(in) :: mol
   type(xtb_calculator), intent(in) :: calc
   type(wavefunction_type), intent(inout) :: wfn

   integer :: iat, ii, ish, spin

   call get_occupation(mol, calc%bas, calc%h0, wfn%nocc, wfn%n0at, wfn%n0sh)
   do spin = 1, size(wfn%qat, 2)
      do iat = 1, size(wfn%qat, 1)
         ii = calc%bas%ish_at(iat)
         do ish = 1, calc%bas%nsh_at(iat)
            wfn%qsh(ii+ish, spin) = (wfn%n0sh(ii+ish) / wfn%n0at(iat)) * wfn%qat(iat, spin)
         end do
      end do
   end do
end subroutine shell_partition

end module tblite_wavefunction_guess
