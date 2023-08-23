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

module test_ceh
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_context_type, only : context_type
   use tblite_lapack_solver, only : lapack_solver, lapack_algorithm
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction, &
   & new_wavefunction_derivative, wavefunction_derivative_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   use tblite_ceh_calculator, only : ceh_calculator
   use tblite_ceh_ceh, only: ceh_guess, new_ceh_calculator
   use tblite_blas, only: gemv
   implicit none
   private

   public :: collect_ceh

   real(wp), parameter :: acc = 0.01_wp
   real(wp), parameter :: thr = 100_wp*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: kt = 300.0_wp * 3.166808578545117e-06_wp

contains


!> Collect all exported unit tests
subroutine collect_ceh(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("dipmom-mol", test_d_mb01), &
      new_unittest("charge-mol", test_q_mb01) &
      ]

end subroutine collect_ceh


! subroutine numdiff_grad(ctx, mol, calc, wfn, numgrad)
!    type(context_type), intent(inout) :: ctx
!    type(structure_type), intent(in) :: mol
!    type(xtb_calculator), intent(in) :: calc
!    type(wavefunction_type), intent(in) :: wfn
!    real(wp), intent(out) :: numgrad(:, :)
! 
!    integer :: iat, ic
!    real(wp) :: er, el
!    type(structure_type) :: moli
!    type(wavefunction_type) :: wfni
!    real(wp), parameter :: step = 1.0e-9_wp
! 
!    do iat = 1, mol%nat
!       do ic = 1, 3
!          moli = mol
!          wfni = wfn
!          moli%xyz(ic, iat) = mol%xyz(ic, iat) + step
!          call xtb_singlepoint(ctx, moli, calc, wfni, acc, er, verbosity=0)
! 
!          moli = mol
!          wfni = wfn
!          moli%xyz(ic, iat) = mol%xyz(ic, iat) - step
!          call xtb_singlepoint(ctx, moli, calc, wfni, acc, el, verbosity=0)
! 
!          numgrad(ic, iat) = 0.5_wp*(er - el)/step
!       end do
!    end do
! end subroutine numdiff_grad


subroutine test_d_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(ceh_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(wavefunction_derivative_type) :: dwfn 
   real(wp) :: dipole(3), tmp(3)
   real(wp), parameter :: ref(3) = reshape([ &
         0.584361099036660_wp, &
         -1.47304239280996_wp, &
         -2.25861915370679_wp], shape(ref))

   call get_structure(mol, "MB16-43", "01")

   call new_ceh_calculator(calc, mol)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   call ceh_guess(ctx, calc, mol, error, wfn, dwfn)
   tmp = 0.0_wp
   dipole = 0.0_wp
   call gemv(mol%xyz, wfn%qat(:, 1), tmp)
   dipole(:) = tmp + sum(wfn%dpat(:, :, 1), 2)
   call check(error, sum(dipole), sum(ref), thr=1e-5_wp)

end subroutine test_d_mb01

subroutine test_q_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(ceh_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(wavefunction_derivative_type) :: dwfn 
   real(wp) :: dipole(3)
   real(wp), parameter :: ref = 0.5340202097_wp

   call get_structure(mol, "MB16-43", "01")

   call new_ceh_calculator(calc, mol)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   call ceh_guess(ctx, calc, mol, error, wfn, dwfn)
   call check(error, wfn%qat(mol%nat,1), ref, thr=1e-6_wp)

end subroutine test_q_mb01

end module test_ceh
