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

module test_solvation_shift
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_container, only : container_cache
   use tblite_scf_potential, only : potential_type
   use tblite_solvation_alpb
   use tblite_solvation_shift
   use tblite_data_shift
   use tblite_solvation_data
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: collect_solvation_shift


   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains


!> Collect all exported unit tests
subroutine collect_solvation_shift(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("solvation-shift", test_e_shift) &
      ]

end subroutine collect_solvation_shift


subroutine test_e(error, mol, input, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Solvation model input
   type(shift_input), intent(in) :: input

   !> Reference energy
   real(wp), intent(in) :: ref

   type(shift_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp) :: energy(mol%nat)
   type(shift_input) :: tmpinput

   energy = 0.0_wp

   tmpinput = input
   call get_shift_param(tmpinput, error)
   if(allocated(error))then
     call test_failed(error, "Failed to get solvation shift parameters")
   endif
   solv = shift_solvation(tmpinput)

   call solv%update(mol, cache)
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)
   call solv%get_engrad(mol, cache, energy)

   if (abs(sum(energy) - ref) > thr) then
      call test_failed(error, "Energy does not match reference")
      print *, sum(energy),'reference:',ref
   end if
end subroutine test_e

subroutine test_e_shift(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "04")
   call test_e(error, mol, shift_input(alpb=.true., solvent='water', method='gfn2'),          1.0807596985111987d-003)
   call test_e(error, mol, shift_input(alpb=.false., solvent='water', method='gfn2'),         1.8574431274405378d-003)
   call test_e(error, mol, shift_input(alpb=.true., solvent='water', method='gfn1'),          4.2456066722865648d-003)
   call test_e(error, mol, shift_input(alpb=.true., solvent='water', state=2, method='gfn2'), 7.9051551874469422d-003)
   call test_e(error, mol, shift_input(alpb=.true., solvent='water', state=3, method='gfn2'), 4.1120060345887655d-003)
   if(allocated(error)) return

end subroutine test_e_shift

end module test_solvation_shift
