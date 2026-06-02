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

module test_molden
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_env_testing, only : new_unittest, &
      & unittest_type, check
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_basis_type, only : basis_type, cgto_type, &
      & new_basis
   use tblite_context_type, only : context_type
   use tblite_io_molden, only : load_molden, save_molden
   use tblite_wavefunction, only : wavefunction_type, &
      & new_wavefunction, eeq_guess
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: collect_molden

   real(wp), parameter :: kt = 300.0_wp * 3.166808578545117e-06_wp
   real(wp), parameter :: acc = 0.01_wp

contains

subroutine collect_molden(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("roundtrip", test_roundtrip), &
      new_unittest("reordered-sections", test_reordered_sections), &
      new_unittest("restart", test_restart_from_molden), &
      new_unittest("restart-unrestriced", test_restart_uhf_from_molden), &
      new_unittest("num-prim-fail", test_num_prim_fail, should_fail=.true.), &
      new_unittest("ang-mom-fail", test_angmom_fail, should_fail=.true.), &
      new_unittest("missing-atoms-fail", test_missing_atoms_fail, should_fail=.true.), &
      new_unittest("missing-gto-fail", test_missing_gto_fail, should_fail=.true.), &
      new_unittest("missing-mo-fail", test_missing_mo_fail, should_fail=.true.) &
      ]

end subroutine collect_molden

!> Test Molden roundtrip with standard section order
subroutine test_roundtrip(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check_roundtrip(".molden-roundtrip.molden", .false., error)
end subroutine test_roundtrip

!> Test Molden roundtrip with reordered sections
subroutine test_reordered_sections(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check_roundtrip(".molden-reordered.molden", .true., error)
end subroutine test_reordered_sections

subroutine check_roundtrip(filename, reorder, error)
   !> Molden file name
   character(len=*), intent(in) :: filename
   !> Whether to reorder sections before loading
   logical, intent(in) :: reorder
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=:), allocatable :: input
   type(structure_type) :: mol, mol_loaded
   type(basis_type) :: bas, bas_loaded
   type(wavefunction_type) :: wfn, wfn_loaded

   call make_roundtrip_data(mol, bas, wfn)

   call remove_file(filename)
   call save_molden(filename, mol, bas, wfn, error)
   if (allocated(error)) return

   if (reorder) then
      input = filename//".input"
      call remove_file(input)
      call rewrite_reordered(filename, input)
   else
      input = filename
   end if

   call load_molden(input, mol_loaded, bas_loaded, wfn_loaded, error)
   if (allocated(error)) return
   call remove_file(filename)
   call remove_file(input)

   call check_structure(error, mol_loaded, mol)
   if (allocated(error)) return
   call check_basis(error, bas_loaded, bas)
   if (allocated(error)) return
   call check_wavefunction(error, wfn_loaded, wfn)
   if (allocated(error)) return
end subroutine check_roundtrip

subroutine check_structure(error, actual, expected)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Read structure data
   type(structure_type), intent(in) :: actual
   !> Reference structure data
   type(structure_type), intent(in) :: expected

   integer :: i

   call check(error, actual%nat, expected%nat)
   if (allocated(error)) return
   call check(error, all(actual%num(actual%id) == expected%num(expected%id)), &
      & "Nuclear charges changed in Molden round trip")
   if (allocated(error)) return
   call check(error, all(abs(actual%xyz - expected%xyz) <= epsilon(1.0_wp)), &
      & "Nuclear coordinates changed in Molden round trip")
   if (allocated(error)) return
   call check(error, size(actual%sym) == size(expected%sym), &
      & "Nuclear symbol table size changed in Molden round trip")
   if (allocated(error)) return
   do i = 1, size(expected%sym)
      call check(error, trim(actual%sym(i)) == trim(expected%sym(i)), &
         & "Nuclear symbol " // trim(expected%sym(i)) // " changed in Molden round trip")
      if (allocated(error)) return
   end do
   call check(error, abs(actual%charge - expected%charge) <= epsilon(1.0_wp), &
      & "Molecular charge changed in Molden round trip")
   if (allocated(error)) return
   call check(error, actual%uhf, expected%uhf, &
      & "Molecular spin multiplicity changed in Molden round trip")
   if (allocated(error)) return
   if (allocated(expected%periodic)) then
      if (.not. allocated(actual%periodic)) then
         call fatal_error(error, "Periodicity missing in Molden round trip")
         return
      end if
      call check(error, all(actual%periodic .eqv. expected%periodic), &
         & "Periodicity changed in Molden round trip")
      if (allocated(error)) return
   else
      call check(error, .not. allocated(actual%periodic), &
         & "Unexpected periodicity in Molden round trip")
      if (allocated(error)) return
   end if
   if (allocated(expected%lattice)) then
      if (.not. allocated(actual%lattice)) then
         call fatal_error(error, "Lattice vectors missing in Molden round trip")
         return
      end if
      call check(error, all(abs(actual%lattice - expected%lattice) <= epsilon(1.0_wp)), &
         & "Lattice vectors changed in Molden round trip")
   else
      call check(error, .not. allocated(actual%lattice), &
         & "Unexpected lattice vectors in Molden round trip")
   end if
end subroutine check_structure

subroutine check_basis(error, actual, expected)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Read basis set information
   type(basis_type), intent(in) :: actual
   !> Reference basis set information
   type(basis_type), intent(in) :: expected

   integer :: isp, ish, nprim

   call check(error, actual%nsh, expected%nsh, &
      & "Basis shell count changed in Molden round trip")
   if (allocated(error)) return
   call check(error, actual%nao, expected%nao, &
      & "Basis AO count changed in Molden round trip")
   if (allocated(error)) return
   call check(error, actual%nao_cart, expected%nao_cart, &
      & "Basis Cartesian AO count changed in Molden round trip")
   if (allocated(error)) return
   call check(error, actual%maxl, expected%maxl, &
      & "Basis maximum angular momentum changed in Molden round trip")
   if (allocated(error)) return
   call check(error, all(actual%sh2at == expected%sh2at), &
      & "Basis nucleus map changed in Molden round trip")
   if (allocated(error)) return
   call check(error, all(actual%nao_sh == expected%nao_sh), &
      & "Basis shell AO counts changed in Molden round trip")
   if (allocated(error)) return
   call check(error, all(actual%iao_sh == expected%iao_sh), &
      & "Basis shell AO offsets changed in Molden round trip")
   if (allocated(error)) return
   call check(error, all(actual%ao2sh == expected%ao2sh), &
      & "Basis AO shell map changed in Molden round trip")
   if (allocated(error)) return
   call check(error, all(actual%ao2at == expected%ao2at), &
      & "Basis AO atom map changed in Molden round trip")
   if (allocated(error)) return

   ! Check CGTO primitives
   do isp = 1, size(expected%cgto, 2)
      do ish = 1, expected%nsh_id(isp)
         call check(error, actual%cgto(ish, isp)%ang, expected%cgto(ish, isp)%ang, &
            & "Shell angular momentum changed in Molden round trip")
         if (allocated(error)) return
         call check(error, actual%cgto(ish, isp)%nprim, expected%cgto(ish, isp)%nprim, &
            & "Shell primitive count changed in Molden round trip")
         nprim = expected%cgto(ish, isp)%nprim
         if (allocated(error)) return
         call check(error, all(abs(actual%cgto(ish, isp)%alpha(:nprim) - &
            & expected%cgto(ish, isp)%alpha(:nprim)) <= epsilon(1.0_wp)), &
            & "Primitive exponents changed in Molden round trip")
         if (allocated(error)) return
         call check(error, all(abs(actual%cgto(ish, isp)%coeff(:nprim) - &
            & expected%cgto(ish, isp)%coeff(:nprim)) <= epsilon(1.0_wp)), &
            & "Primitive coefficients changed in Molden round trip")
         if (allocated(error)) return
      end do
   end do
end subroutine check_basis

subroutine check_wavefunction(error, actual, expected)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Read wavefunction data
   type(wavefunction_type), intent(in) :: actual
   !> Reference wavefunction data
   type(wavefunction_type), intent(in) :: expected

   call check(error, actual%nspin, expected%nspin)
   if (allocated(error)) return
   call check(error, abs(actual%nocc - expected%nocc) <= epsilon(1.0_wp), &
      & "Wavefunction electron count changed in Molden round trip")
   if (allocated(error)) return
   call check(error, abs(actual%nuhf - expected%nuhf) <= epsilon(1.0_wp), &
      & "Wavefunction unpaired electron count changed in Molden round trip")
   if (allocated(error)) return
   call check(error, all(abs(actual%nel - expected%nel) <= epsilon(1.0_wp)), &
      & "Wavefunction spin electron counts changed in Molden round trip")
   if (allocated(error)) return
   call check(error, all(abs(actual%emo - expected%emo) <= epsilon(1.0_wp)), &
      & "Wavefunction MO energies changed in Molden round trip")
   if (allocated(error)) return
   call check(error, all(abs(actual%n0at - expected%n0at) <= epsilon(1.0_wp)), &
      & "Wavefunction atomic reference occupation changed in Molden round trip")
   if (allocated(error)) return
   call check(error, all(abs(actual%focc - expected%focc) <= epsilon(1.0_wp)), &
      & "Wavefunction MO occupations changed in Molden round trip")
   if (allocated(error)) return
   call check(error, all(abs(actual%coeff - expected%coeff) <= epsilon(1.0_wp)), &
      & "Wavefunction MO coefficients changed in Molden round trip")
end subroutine check_wavefunction

subroutine make_roundtrip_data(mol, bas, wfn)
   !> Molecular structure data
   type(structure_type), intent(out) :: mol
   !> Basis set information
   type(basis_type), intent(out) :: bas
   !> Wavefunction data
   type(wavefunction_type), intent(out) :: wfn

   integer :: iao
   integer, allocatable :: nshell(:)
   real(wp) :: lattice(3, 3)
   type(cgto_type), allocatable :: cgto(:, :)

   lattice = reshape([ &
      & 8.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 9.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 10.0_wp], [3, 3])

   call new(mol, [1, 4], reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 1.4_wp], [3, 2]), &
      & charge=0.0_wp, uhf=1, lattice=lattice)

   allocate(nshell(mol%nid), cgto(1, mol%nid))
   nshell(:) = 1
   cgto(:, :) = cgto_type()

   cgto(1, 1)%ang = 0
   cgto(1, 1)%nprim = 1
   cgto(1, 1)%alpha(1) = 1.2_wp
   cgto(1, 1)%coeff(1) = 0.8_wp

   cgto(1, 2)%ang = 1
   cgto(1, 2)%nprim = 2
   cgto(1, 2)%alpha(:2) = [0.9_wp, 0.4_wp]
   cgto(1, 2)%coeff(:2) = [0.7_wp, 0.2_wp]

   call new_basis(bas, mol, nshell, cgto, 1.0_wp)

   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, 2, 0.0_wp)
   wfn%nocc = 3.0_wp
   wfn%nuhf = 1.0_wp
   wfn%nel = [2.0_wp, 1.0_wp]
   wfn%n0at = [1.0_wp, 2.0_wp]
   wfn%n0sh = [1.0_wp, 2.0_wp]
   wfn%emo(:, 1) = [-0.7_wp, -0.2_wp, 0.1_wp, 0.3_wp]
   wfn%emo(:, 2) = [-0.6_wp, -0.1_wp, 0.2_wp, 0.4_wp]
   wfn%focc(:, 1) = [1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp]
   wfn%focc(:, 2) = [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]

   wfn%coeff(:, :, :) = 0.0_wp
   do iao = 1, bas%nao
      wfn%coeff(iao, iao, 1) = 1.0_wp
      wfn%coeff(iao, iao, 2) = 1.0_wp
   end do
   wfn%coeff(2, 3, 1) = 0.25_wp
   wfn%coeff(3, 2, 1) = -0.15_wp
   wfn%coeff(4, 4, 1) = 0.85_wp
   wfn%coeff(2, 4, 2) = -0.20_wp
   wfn%coeff(4, 2, 2) = 0.30_wp
end subroutine make_roundtrip_data

subroutine rewrite_reordered(input, output)
   !> Input file name
   character(len=*), intent(in) :: input
   !> Output file name
   character(len=*), intent(in) :: output

   character(len=512) :: lines(512)
   integer :: iline, ios, nline, unit

   nline = 0
   open(newunit=unit, file=input, status="old", action="read")
   do
      read(unit, '(a)', iostat=ios) lines(nline + 1)
      if (ios /= 0) exit
      nline = nline + 1
   end do
   close(unit)

   open(newunit=unit, file=output, &
      & status="replace", action="write")
   call write_section(unit, lines(:nline), "molden format")
   call write_section(unit, lines(:nline), "title")
   call write_section(unit, lines(:nline), "atoms")
   call write_section(unit, lines(:nline), "gto")
   call write_section(unit, lines(:nline), "6d")
   call write_section(unit, lines(:nline), "10f")
   call write_section(unit, lines(:nline), "15g")
   call write_section(unit, lines(:nline), "cell")
   call write_section(unit, lines(:nline), "core")
   call write_section(unit, lines(:nline), "pseudo")
   call write_section(unit, lines(:nline), "nval")
   call write_section(unit, lines(:nline), "mo")
   close(unit)

end subroutine rewrite_reordered

subroutine write_section(unit, lines, name)
   !> Output file unit
   integer, intent(in) :: unit
   !> All lines from the original file
   character(len=*), intent(in) :: lines(:)
   !> Section name to write
   character(len=*), intent(in) :: name

   integer :: first, iline, last

   first = 0
   do iline = 1, size(lines)
      if (section_id(lines(iline)) == name) then
         first = iline
         exit
      end if
   end do
   if (first == 0) return

   last = size(lines)
   do iline = first + 1, size(lines)
      if (len_trim(section_id(lines(iline))) > 0) then
         last = iline - 1
         exit
      end if
   end do

   do iline = first, last
      write(unit, '(a)') trim(lines(iline))
   end do
end subroutine write_section

function section_id(line) result(id)
   !> Input line
   character(len=*), intent(in) :: line
   character(len=32) :: id

   character(len=len(line)) :: work
   integer :: right, i, ic

   id = ""
   work = adjustl(line)
   if (len_trim(work) < 3) return
   if (work(1:1) /= "[") return

   right = index(work, "]")
   if (right <= 2) return
   id = adjustl(work(2:right-1))
   ! Lower case
   do i = 1, len_trim(id)
      ic = iachar(id(i:i))
      if (ic >= iachar("A") .and. ic <= iachar("Z")) then
         id(i:i) = achar(ic + iachar("a") - iachar("A"))
      end if
   end do
end function section_id

subroutine remove_file(filename)
   !> Molden file name to remove
   character(len=*), intent(in) :: filename

   integer :: io, stat
   logical :: exist

   inquire(file=filename, exist=exist)
   if (exist) then
      open(newunit=io, file=filename, status="old", action="write", iostat=stat)
      if (stat == 0) close(io, status="delete")
   end if
end subroutine remove_file

subroutine test_restart_from_molden(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check_restart_from_molden(".molden-restart.molden", 1, error)
end subroutine test_restart_from_molden

subroutine test_restart_uhf_from_molden(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check_restart_from_molden(".molden-restart-uhf.molden", 2, error)
end subroutine test_restart_uhf_from_molden

subroutine check_restart_from_molden(filename, nspin, error)
   !> Molden file name
   character(len=*), intent(in) :: filename
   !> Number of spin channels in the wavefunction to test
   integer, intent(in) :: nspin
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol, mol_loaded
   type(basis_type) :: bas, bas_loaded
   type(wavefunction_type) :: wfn, wfn_loaded
   type(xtb_calculator) :: calc
   real(wp) :: energy

   ! Perform GFN2-xTB calculation to save in trexio
   call make_restart_data(nspin, mol, bas, wfn, energy)

   call remove_file(filename)
   call save_molden(filename, mol, bas, wfn, error)
   if (allocated(error)) return

   call load_molden(filename, mol_loaded, bas_loaded, wfn_loaded, error)
   if (allocated(error)) return
   call remove_file(filename)

   ! Setup a new calculator
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return

   !> Check that loaded structure and basis set
   call check_structure(error, mol_loaded, mol)
   if (allocated(error)) return
   call check_basis(error, bas_loaded, calc%bas)
   if (allocated(error)) return

   ! Copy all necessary restart wavefunction information
   wfn%coeff = wfn_loaded%coeff
   wfn%focc = wfn_loaded%focc
   wfn%emo = wfn_loaded%emo
   wfn%nocc = wfn_loaded%nocc
   wfn%nel = wfn_loaded%nel
   ! Copy also the basis set to check for correct normalization and mappings
   calc%bas = bas_loaded

   ! Check for immediate convergence with Molden guess
   calc%max_iter = 2
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)
   call check(error, .not.ctx%failed(), &
      & "Calculation did not converge in < 3 iterations with Molden guess")
end subroutine check_restart_from_molden

subroutine make_restart_data(nspin, mol, bas, wfn, energy)
   !> Number of spin channels
   integer, intent(in) :: nspin
   !> Molecular structure data
   type(structure_type), intent(out) :: mol
   !> Basis set information
   type(basis_type), intent(out) :: bas
   !> Wavefunction data
   type(wavefunction_type), intent(out) :: wfn
   !> Total energy
   real(wp), intent(out) :: energy

   type(context_type) :: ctx
   type(xtb_calculator) :: calc
   type(error_type), allocatable :: error

   call get_structure(mol, "X23", "oxacb")

   call new_gfn2_calculator(calc, mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, nspin, kt)
   call eeq_guess(mol, calc, wfn, error)

   energy = 0.0_wp
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)

   bas = calc%bas
end subroutine make_restart_data

subroutine test_num_prim_fail(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: filename = ".molden-num-prim-fail.molden"
   type(structure_type) :: mol, mol_loaded
   type(basis_type) :: bas, bas_loaded
   type(wavefunction_type) :: wfn, wfn_loaded

   call make_roundtrip_data(mol, bas, wfn)
   bas%cgto(1, 2)%nprim = 0

   call remove_file(filename)
   call save_molden(filename, mol, bas, wfn, error)
   if (allocated(error)) return

   call load_molden(filename, mol_loaded, bas_loaded, wfn_loaded, error)
   call remove_file(filename)
end subroutine test_num_prim_fail

subroutine test_angmom_fail(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: filename = ".molden-angmom-fail.molden"
   type(structure_type) :: mol, mol_loaded
   type(basis_type) :: bas, bas_loaded
   type(wavefunction_type) :: wfn, wfn_loaded

   call make_roundtrip_data(mol, bas, wfn)
   bas%cgto(1, 2)%ang = 5

   call remove_file(filename)
   call save_molden(filename, mol, bas, wfn, error)
   call remove_file(filename)
end subroutine test_angmom_fail

subroutine test_missing_atoms_fail(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check_missing_section(".molden-missing-atoms.molden", "atoms", error)
end subroutine test_missing_atoms_fail

subroutine test_missing_gto_fail(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check_missing_section(".molden-missing-gto.molden", "gto", error)
end subroutine test_missing_gto_fail

subroutine test_missing_mo_fail(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check_missing_section(".molden-missing-mo.molden", "mo", error)
end subroutine test_missing_mo_fail

subroutine check_missing_section(filename, missing, error)
   !> Molden file name
   character(len=*), intent(in) :: filename
   !> Section to remove
   character(len=*), intent(in) :: missing
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=:), allocatable :: input
   type(structure_type) :: mol, mol_loaded
   type(basis_type) :: bas, bas_loaded
   type(wavefunction_type) :: wfn, wfn_loaded
   character(len=512) :: lines(1024)
   integer :: ios, nline, unit

   call make_roundtrip_data(mol, bas, wfn)

   call remove_file(filename)
   call save_molden(filename, mol, bas, wfn, error)
   if (allocated(error)) return

   input = filename//".input"
   call remove_file(input)

   ! Read Molden file into memory
   nline = 0
   open(newunit=unit, file=filename, status="old", action="read")
   do
      read(unit, '(a)', iostat=ios) lines(nline + 1)
      if (ios /= 0) exit
      nline = nline + 1
   end do
   close(unit)
   ! Write new file with the specified section missing
   open(newunit=unit, file=input, status="replace", action="write")
   if (trim(missing) /= "atoms") &
      call write_section(unit, lines(:nline), "atoms")
   if (trim(missing) /= "gto") &
      call write_section(unit, lines(:nline), "gto")
   if (trim(missing) /= "mo") &
      call write_section(unit, lines(:nline), "mo")
   close(unit)

   call load_molden(input, mol_loaded, bas_loaded, wfn_loaded, error)
   call remove_file(filename)
   call remove_file(input)
end subroutine check_missing_section

end module test_molden
