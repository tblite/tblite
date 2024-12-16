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
   use tblite_param, only : param_record, charge_record, dispersion_record, element_record, &
      & halogen_record, hamiltonian_record, multipole_record, repulsion_record, &
      & thirdorder_record, param_mask, count
   use tblite_param_molecular_moments, only:  molecular_multipole_record
   use tblite_toml, only : toml_table
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
      new_unittest("param-empty", test_param_empty, should_fail=.true.), &
      new_unittest("param-minimal", test_param_minimal), &
      new_unittest("param-invalid", test_param_invalid, should_fail=.true.), &
      new_unittest("element-empty", test_element_empty, should_fail=.true.), &
      new_unittest("charge-empty", test_charge_empty, should_fail=.true.), &
      new_unittest("dispersion-empty", test_dispersion_empty, should_fail=.true.), &
      new_unittest("halogen-empty", test_halogen_empty, should_fail=.true.), &
      new_unittest("hamiltonian-empty", test_hamiltonian_empty, should_fail=.true.), &
      new_unittest("multipole-empty", test_multipole_empty, should_fail=.true.), &
      new_unittest("repulsion-empty", test_repulsion_empty, should_fail=.true.), &
      new_unittest("thirdorder-empty", test_thirdorder_empty), &
      new_unittest("mol-multipole-empty", test_mol_multipole_empty), &
      new_unittest("mask-a", test_mask_gfn2), &
      new_unittest("mask-b", test_mask_gfn1), &
      new_unittest("gfn2-xtb-a", test_gfn2_mb02), &
      new_unittest("gfn2-xtb-b", test_gfn2_mb05), &
      new_unittest("gfn2-xtb-c", test_gfn2_round_trip), &
      new_unittest("gfn1-xtb-a", test_gfn1_mb01), &
      new_unittest("gfn1-xtb-b", test_gfn1_mb04), &
      new_unittest("gfn1-xtb-c", test_gfn1_round_trip), &
      new_unittest("ipea1-xtb-a", test_ipea1_mb03), &
      new_unittest("ipea1-xtb-b", test_ipea1_mb06), &
      new_unittest("ipea1-xtb-c", test_ipea1_round_trip) &
      ]

end subroutine collect_xtb_param


subroutine test_param_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(param_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_param_empty


subroutine test_param_minimal(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: io
   type(param_record) :: param

   open(newunit=io, status="scratch")
   write(io, '(a)') &
      "[hamiltonian.xtb]", &
      "wexp = 5.0000000000000000E-01", &
      "kpol = 2.0000000000000000E+00", &
      "enscale = 2.0000000000000000E-02", &
      "cn = ""gfn""", &
      "shell = {ss=1.85, pp=2.23, dd=2.23, sd=2.00, pd=2.00}", &
      "[element.H]", &
      "shells = [ ""1s"" ]", &
      "levels = [ -1.0707210999999999E+01 ]", &
      "slater = [ 1.2300000000000000E+00 ]", &
      "ngauss = [ 3 ]", &
      "refocc = [ 1.0000000000000000E+00 ]", &
      "shpoly = [ -9.5361800000000000E-03 ]", &
      "kcn = [ -5.0000000000000003E-02 ]", &
      "gam = 4.0577099999999999E-01", &
      "lgam = [ 1.0000000000000000E+00 ]", &
      "gam3 = 8.0000000000000016E-02", &
      "zeff = 1.1053880000000000E+00", &
      "arep = 2.2137169999999999E+00", &
      "xbond = 0.0000000000000000E+00", &
      "dkernel = 5.5638889999999996E-02", &
      "qkernel = 2.7430999999999999E-04", &
      "mprad = 1.3999999999999999E+00", &
      "mpvcn = 1.0000000000000000E+00", &
      "[element.C]", &
      "shells = [ ""2s"", ""2p"" ]", &
      "levels = [ -1.3970922000000002E+01, -1.0063292000000001E+01 ]", &
      "slater = [ 2.0964320000000001E+00, 1.8000000000000000E+00 ]", &
      "ngauss = [ 4, 4 ]", &
      "refocc = [ 1.0000000000000000E+00, 3.0000000000000000E+00 ]", &
      "shpoly = [ -2.2943210000000002E-02, -2.7110200000000002E-03 ]", &
      "kcn = [ -1.0214400000000000E-02, 1.6165700000000002E-02 ]", &
      "gam = 5.3801500000000002E-01", &
      "lgam = [ 1.0000000000000000E+00, 1.1056357999999999E+00 ]", &
      "gam3 = 1.5000000000000002E-01", &
      "zeff = 4.2310780000000001E+00", &
      "arep = 1.2476550000000000E+00", &
      "xbond = 0.0000000000000000E+00", &
      "dkernel = -4.1167399999999998E-03", &
      "qkernel = 2.1358300000000000E-03", &
      "mprad = 3.0000000000000000E+00", &
      "mpvcn = 3.0000000000000000E+00", &
      ""
   rewind io

   call param%load(io, error)
   close(io)
end subroutine test_param_minimal


subroutine test_param_invalid(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: io
   type(param_record) :: param

   open(newunit=io, status="scratch")
   write(io, '(a)') &
      "[hamiltonian.xtb]", &
      "wexp = 5.0000000000000000E-01", &
      "kpol = 2.0000000000000000E+00", &
      "enscale = 2.0000000000000000E-02", &
      "cn = ""gfn""", &
      "shell = {ss=1.85, pp=2.23, dd=2.23, sd=2.00, pd=2.00"
   rewind io

   call param%load(io, error)
   close(io)
end subroutine test_param_invalid


subroutine test_element_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(element_record) :: param

   table = toml_table()
   table%key = "Te"
   call param%load(table, error)
end subroutine test_element_empty

subroutine test_mol_multipole_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(molecular_multipole_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_mol_multipole_empty

subroutine test_charge_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(charge_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_charge_empty


subroutine test_dispersion_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(dispersion_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_dispersion_empty


subroutine test_halogen_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(halogen_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_halogen_empty


subroutine test_hamiltonian_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(hamiltonian_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_hamiltonian_empty


subroutine test_multipole_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(multipole_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_multipole_empty


subroutine test_repulsion_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(repulsion_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_repulsion_empty


subroutine test_thirdorder_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(thirdorder_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_thirdorder_empty


subroutine test_mask_gfn2(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(param_record), target :: param, base
   type(param_mask) :: mask1, mask2
   real(wp), allocatable :: array(:)
   type(toml_table) :: table
   integer :: io

   call export_gen_param("gfn2", base)

   mask1%ref => base%record
   open(newunit=io, status="scratch")
   write(io, '(a)') &
      "hamiltonian = {}", &
      "dispersion = {}", &
      "repulsion = {}", &
      "charge = {}", &
      "thirdorder = {}", &
      "multipole = {}", &
      "[element]", &
      "H = {slater=[false]}", &
      "C = {}", &
      ""
   rewind io

   call mask1%load(io, error)
   close(io)
   if (allocated(error)) return

   call check(error, count(mask1), 26)
   if (allocated(error)) return

   table = toml_table()
   call mask1%dump(table, error)
   if (allocated(error)) return

   mask2%ref => base%record
   call mask2%load(table, error)
   if (allocated(error)) return

   call check(error, count(mask2), 26)
   if (allocated(error)) return

   allocate(array(count(mask2)))
   call base%dump(array, mask2, error)
   if (allocated(error)) return

   call param%load(array, base, mask2, error)
   if (allocated(error)) return

end subroutine test_mask_gfn2


subroutine test_mask_gfn1(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(param_record), target :: param, base
   type(param_mask) :: mask
   real(wp), allocatable :: array(:)
   integer :: io

   call export_gen_param("gfn1", base)

   mask%ref => base%record
   open(newunit=io, status="scratch")
   write(io, '(a)') &
      "hamiltonian = {}", &
      "dispersion = {}", &
      "repulsion = {}", &
      "charge = {}", &
      "thirdorder = {}", &
      "halogen = {}", &
      "[element]", &
      "H = {slater=[false, true], gam3=false}", &
      "C = {}", &
      ""
   rewind io

   call mask%load(io, error)
   close(io)
   if (allocated(error)) return

   call check(error, count(mask), 28)
   if (allocated(error)) return

   allocate(array(count(mask)))
   call base%dump(array, mask, error)
   if (allocated(error)) return

   call param%load(array, base, mask, error)
   if (allocated(error)) return

end subroutine test_mask_gfn1


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


subroutine new_gen_calculator(calc, method, mol, error)
   type(xtb_calculator), intent(out) :: calc
   character(len=*), intent(in) :: method
   type(structure_type), intent(in) :: mol
   type(error_type), allocatable, intent(out) :: error
   select case(method)
   case("gfn1")
      call new_gfn1_calculator(calc, mol, error)
   case("gfn2")
      call new_gfn2_calculator(calc, mol, error)
   case("ipea1")
      call new_ipea1_calculator(calc, mol, error)
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

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy1, verbosity=0)

   call new_gen_calculator(calc, method, mol, error)
   if (allocated(error)) return

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy2, verbosity=0)

   call check(error, energy2, energy1, thr=thr)

end subroutine test_gen


subroutine test_round_trip(mol, method, error)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Method name
   character(len=*), intent(in) :: method
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(xtb_calculator) :: calc
   type(param_record) :: param1, param2
   type(wavefunction_type) :: wfn
   real(wp) :: energy1, energy2
   character(len=:), allocatable :: io

   call export_gen_param(method, param1)

   io = "." // method // "-" // get_name() // ".toml"
   call param1%dump(io, error)
   if (.not.allocated(error)) then
      call param2%load(io, error)
   end if
   call delete_file(io)
   if (allocated(error)) return

   call new_xtb_calculator(calc, mol, param1, error)
   if (allocated(error)) return

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy1, verbosity=0)

   call new_xtb_calculator(calc, mol, param2, error)
   if (allocated(error)) return

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy2, verbosity=0)

   call check(error, energy2, energy1, thr=thr)

end subroutine test_round_trip


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


subroutine test_gfn2_round_trip(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "07")
   call test_round_trip(mol, "gfn2", error)

end subroutine test_gfn2_round_trip


subroutine test_gfn1_round_trip(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "08")
   call test_round_trip(mol, "gfn1", error)

end subroutine test_gfn1_round_trip


subroutine test_ipea1_round_trip(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "09")
   call test_round_trip(mol, "ipea1", error)

end subroutine test_ipea1_round_trip


function get_name() result(name)
   character(len=20) :: name
   real :: val

   call random_number(val)
   write(name, '(a, z8.8)') "tblite-test-", int(val*1.0e9)
end function get_name


subroutine delete_file(file)
   character(len=*), intent(in) :: file
   integer :: unit
   logical :: exist
   inquire(file=file, exist=exist)
   if (exist) then
      open(newunit=unit, file=file)
      close(unit, status="delete")
   end if
end subroutine delete_file


end module test_xtb_param
