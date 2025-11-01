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

module test_compressed_restart_calculator
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check
   use mctc_io, only : structure_type, new
   use tblite_context_type, only : context_type
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction, shell_partition
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   use tblite_scf_iterator, only : get_qat_from_qsh
   implicit none
   private

   public :: collect_compressed_restart_calculator

   real(wp), parameter :: acc = 0.01_wp
   real(wp), parameter :: thr = 1e-6_wp
   real(wp), parameter :: kt = 300.0_wp * 3.166808578545117e-06_wp

contains


!> Collect all exported unit tests
subroutine collect_compressed_restart_calculator(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("restart-compressed-closed-shell", test_restart_compressed_closed), &
      new_unittest("restart-compressed-spin-polarized", test_restart_compressed_spin) &
      ]

end subroutine collect_compressed_restart_calculator


subroutine test_restart_compressed_closed(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn_full, wfn_restart
   real(wp) :: energy_full, energy_restart
   real(wp), allocatable :: xyz(:, :)
   real(wp), allocatable :: qsh_save(:, :), dpat_save(:, :, :), qmat_save(:, :, :)
   real(wp), allocatable :: qat_from_qsh(:, :)
   integer, allocatable :: numbers(:)

   ! Test system: simple organic molecule
   ! C2O2N2H6 - a conjugated system
   allocate(numbers(12))
   numbers = [6, 6, 7, 7, 1, 1, 1, 1, 1, 1, 8, 8]

   allocate(xyz(3, 12))
   xyz(:, 1) = [-3.81469488143921_wp, +0.09993441402912_wp, 0.00000000000000_wp]
   xyz(:, 2) = [+3.81469488143921_wp, -0.09993441402912_wp, 0.00000000000000_wp]
   xyz(:, 3) = [-2.66030049324036_wp, -2.15898251533508_wp, 0.00000000000000_wp]
   xyz(:, 4) = [+2.66030049324036_wp, +2.15898251533508_wp, 0.00000000000000_wp]
   xyz(:, 5) = [-0.73178529739380_wp, -2.28237795829773_wp, 0.00000000000000_wp]
   xyz(:, 6) = [-5.89039325714111_wp, -0.02589114569128_wp, 0.00000000000000_wp]
   xyz(:, 7) = [-3.71254944801331_wp, -3.73605775833130_wp, 0.00000000000000_wp]
   xyz(:, 8) = [+3.71254944801331_wp, +3.73605775833130_wp, 0.00000000000000_wp]
   xyz(:, 9) = [+0.73178529739380_wp, +2.28237795829773_wp, 0.00000000000000_wp]
   xyz(:, 10) = [+5.89039325714111_wp, +0.02589114569128_wp, 0.00000000000000_wp]
   xyz(:, 11) = [-2.74426102638245_wp, +2.16115570068359_wp, 0.00000000000000_wp]
   xyz(:, 12) = [+2.74426102638245_wp, -2.16115570068359_wp, 0.00000000000000_wp]

   call new(mol, numbers, xyz)

   ! Create calculator
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return

   ! Full calculation
   call new_wavefunction(wfn_full, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   energy_full = 0.0_wp
   call xtb_singlepoint(ctx, mol, calc, wfn_full, acc, energy_full, verbosity=0)
   if (ctx%failed()) then
      call ctx%get_error(error)
      return
   end if

   ! Save compressed restart data
   allocate(qsh_save(calc%bas%nsh, 1))
   allocate(dpat_save(3, mol%nat, 1))
   allocate(qmat_save(6, mol%nat, 1))
   allocate(qat_from_qsh(mol%nat, 1))

   qsh_save(:, 1) = wfn_full%qsh(:, 1)
   dpat_save(:, :, 1) = wfn_full%dpat(:, :, 1)
   qmat_save(:, :, 1) = wfn_full%qpat(:, :, 1)

   ! Restart calculation with compressed data
   call new_wavefunction(wfn_restart, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)

   ! Set initial guess from compressed data
   wfn_restart%qsh(:, 1) = qsh_save(:, 1)
   wfn_restart%dpat(:, :, 1) = dpat_save(:, :, 1)
   wfn_restart%qpat(:, :, 1) = qmat_save(:, :, 1)

   ! Convert shell charges to atom charges
   call get_qat_from_qsh(calc%bas, wfn_restart%qsh, wfn_restart%qat)

   energy_restart = 0.0_wp
   call xtb_singlepoint(ctx, mol, calc, wfn_restart, acc, energy_restart, verbosity=0)
   if (ctx%failed()) then
      call ctx%get_error(error)
      return
   end if

   ! Check that energies match
   call check(error, energy_restart, energy_full, thr=thr, &
      & message="Energies do not match between full and compressed restart calculations")

end subroutine test_restart_compressed_closed


subroutine test_restart_compressed_spin(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn_full, wfn_restart
   real(wp) :: energy_full, energy_restart
   real(wp), allocatable :: xyz(:, :)
   real(wp), allocatable :: qsh_save(:, :), dpat_save(:, :, :), qmat_save(:, :, :)
   real(wp), allocatable :: qat_from_qsh(:, :)
   integer, allocatable :: numbers(:)

   ! Test system: small open-shell system (NO radical)
   allocate(numbers(2))
   numbers = [7, 8]

   allocate(xyz(3, 2))
   xyz(:, 1) = [0.0_wp, 0.0_wp, -1.1_wp]
   xyz(:, 2) = [0.0_wp, 0.0_wp, +1.1_wp]

   call new(mol, numbers, xyz, uhf=1)

   ! Create calculator
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return

      ! Full calculation (spin-polarized)
   call new_wavefunction(wfn_full, mol%nat, calc%bas%nsh, calc%bas%nao, 2, kt)
   energy_full = 0.0_wp
   call xtb_singlepoint(ctx, mol, calc, wfn_full, acc, energy_full, verbosity=0)
   if (ctx%failed()) then
      call ctx%get_error(error)
      return
   end if

   ! Save compressed restart data (spin-resolved)
   allocate(qsh_save(calc%bas%nsh, 2))
   allocate(dpat_save(3, mol%nat, 2))
   allocate(qmat_save(6, mol%nat, 2))
   allocate(qat_from_qsh(mol%nat, 2))

   qsh_save(:, :) = wfn_full%qsh(:, :)
   dpat_save(:, :, :) = wfn_full%dpat(:, :, :)
   qmat_save(:, :, :) = wfn_full%qpat(:, :, :)

   ! Restart calculation with compressed spin-resolved data
   call new_wavefunction(wfn_restart, mol%nat, calc%bas%nsh, calc%bas%nao, 2, kt)

   ! Set initial guess from compressed spin-resolved data
   wfn_restart%qsh(:, :) = qsh_save(:, :)
   wfn_restart%dpat(:, :, :) = dpat_save(:, :, :)
   wfn_restart%qpat(:, :, :) = qmat_save(:, :, :)

   ! Convert shell charges to atom charges
   call get_qat_from_qsh(calc%bas, wfn_restart%qsh, wfn_restart%qat)

   energy_restart = 0.0_wp
   call xtb_singlepoint(ctx, mol, calc, wfn_restart, acc, energy_restart, verbosity=0)
   if (ctx%failed()) then
      call ctx%get_error(error)
      return
   end if

   ! Check that energies match
   call check(error, energy_restart, energy_full, thr=thr, &
      & message="Energies do not match between full and compressed spin-polarized restart calculations")

end subroutine test_restart_compressed_spin

end module test_compressed_restart_calculator
