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

module test_xtb_param
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_context_type, only : context_type
   use tblite_param, only : param_record
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator, new_xtb_calculator
   use tblite_xtb_gfn2, only : export_gfn2_param, new_gfn2_calculator
   use tblite_xtb_gfn1, only : export_gfn1_param, new_gfn1_calculator
   use tblite_xtb_ipea1, only : export_ipea1_param, new_ipea1_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: collect_xtb_param

   real(wp), parameter :: acc = 1.0_wp
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: kt = 300.0_wp * 3.166808578545117e-06_wp

contains


!> Collect all exported unit tests
subroutine collect_xtb_param(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("gfn2-xtb-a", test_gfn2_mb02), &
      new_unittest("gfn2-xtb-b", test_gfn2_mb05), &
      new_unittest("gfn1-xtb-a", test_gfn1_mb01), &
      new_unittest("gfn1-xtb-b", test_gfn1_mb04), &
      new_unittest("ipea1-xtb-a", test_ipea1_mb03), &
      new_unittest("ipea1-xtb-b", test_ipea1_mb06) &
      ]

end subroutine collect_xtb_param


subroutine export_gen_param(method, param)
   character(len=*), intent(in) :: method
   type(param_record), intent(out) :: param
   select case(method)
   case("gfn1")
      call export_gfn1_param(param)
   case("gfn2")
      call export_gfn2_param(param)
   case("ipea1")
      call export_ipea1_param(param)
   end select
end subroutine export_gen_param


subroutine new_gen_calculator(calc, method, mol)
   type(xtb_calculator), intent(out) :: calc
   character(len=*), intent(in) :: method
   type(structure_type), intent(in) :: mol
   select case(method)
   case("gfn1")
      call new_gfn1_calculator(calc, mol)
   case("gfn2")
      call new_gfn2_calculator(calc, mol)
   case("ipea1")
      call new_ipea1_calculator(calc, mol)
   end select
end subroutine new_gen_calculator


subroutine test_gen(mol, method, error)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Method name
   character(len=*), intent(in) :: method
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(xtb_calculator) :: calc
   type(param_record) :: param
   type(wavefunction_type) :: wfn
   real(wp) :: energy1, energy2

   call export_gen_param(method, param)
   call new_xtb_calculator(calc, mol, param, error)
   if (allocated(error)) return

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy1, verbosity=0)

   call new_gen_calculator(calc, method, mol)

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy2, verbosity=0)

   call check(error, energy2, energy1, thr=thr)

end subroutine test_gen


subroutine test_gfn1_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call test_gen(mol, "gfn1", error)

end subroutine test_gfn1_mb01


subroutine test_gfn2_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "02")
   call test_gen(mol, "gfn2", error)

end subroutine test_gfn2_mb02


subroutine test_ipea1_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "03")
   call test_gen(mol, "ipea1", error)

end subroutine test_ipea1_mb03


subroutine test_gfn1_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "04")
   call test_gen(mol, "gfn1", error)

end subroutine test_gfn1_mb04


subroutine test_gfn2_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_gen(mol, "gfn2", error)

end subroutine test_gfn2_mb05


subroutine test_ipea1_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "06")
   call test_gen(mol, "ipea1", error)

end subroutine test_ipea1_mb06


end module test_xtb_param
