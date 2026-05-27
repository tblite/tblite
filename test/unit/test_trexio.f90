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

module test_trexio
   use mctc_env, only : error_type, wp
   use mctc_env_testing, only : new_unittest, unittest_type, check
   use mctc_io, only : structure_type, new
   use tblite_basis_type, only : basis_type, cgto_type, new_basis
   use tblite_features, only : get_tblite_feature
   use tblite_io_trexio, only : load_trexio, save_trexio
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   implicit none
   private

   public :: collect_trexio

contains

subroutine collect_trexio(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   if (get_tblite_feature("trexio")) then
      if (get_tblite_feature("hdf5")) then
         testsuite = [ &
            new_unittest("roundtrip-text", test_roundtrip_text), &
            new_unittest("roundtrip-hdf5", test_roundtrip_hdf5) &
         ]
      else
         testsuite = [new_unittest("roundtrip-text", test_roundtrip_text)]
      end if
   else
      testsuite = [new_unittest("trexio-disabled", test_trexio_disabled)]
   end if
end subroutine collect_trexio

subroutine test_trexio_disabled(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(error_type), allocatable :: actual
   type(structure_type) :: mol
   type(basis_type) :: bas
   type(wavefunction_type) :: wfn
   real(wp) :: energy

   call load_trexio(".trexio-disabled.trexio", mol, bas, wfn, energy, actual)
   call check(error, allocated(actual), "TREXIO read unexpectedly succeeded without TREXIO support")
   if (allocated(actual)) deallocate(actual)
end subroutine test_trexio_disabled

subroutine test_roundtrip_text(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check_roundtrip(".trexio-roundtrip.trexio", error)
end subroutine test_roundtrip_text

subroutine test_roundtrip_hdf5(error)
   type(error_type), allocatable, intent(out) :: error

   call check_roundtrip(".trexio-roundtrip.h5", error)
end subroutine test_roundtrip_hdf5

subroutine check_roundtrip(filename, error)
   !> TREXIO file name
   character(len=*), intent(in) :: filename
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: energy = -12.345_wp
   type(structure_type) :: mol
   type(structure_type) :: mol_loaded
   type(basis_type) :: bas
   type(basis_type) :: bas_loaded
   type(wavefunction_type) :: wfn
   type(wavefunction_type) :: wfn_loaded
   real(wp) :: energy_loaded

   call make_roundtrip_data(mol, bas, wfn)

   call remove_trexio_output(filename)
   call save_trexio(filename, mol, bas, wfn, energy, error)
   if (allocated(error)) return

   call load_trexio(filename, mol_loaded, bas_loaded, wfn_loaded, energy_loaded, error)
   if (allocated(error)) return

   call check_structure(error, mol_loaded, mol)
   if (allocated(error)) return
   call check_basis(error, bas_loaded, bas)
   if (allocated(error)) return
   call check_wavefunction(error, wfn_loaded, wfn)
   if (allocated(error)) return
   call check(error, abs(energy_loaded - energy) <= epsilon(1.0_wp), &
      & "State energy changed in TREXIO round trip")
   if (allocated(error)) return

   call remove_trexio_output(filename)
end subroutine check_roundtrip

subroutine check_structure(error, actual, expected)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Read structure data
   type(structure_type), intent(in) :: actual
   !> Reference structure data
   type(structure_type), intent(in) :: expected

   call check(error, actual%nat, expected%nat)
   if (allocated(error)) return
   call check(error, all(actual%num(actual%id) == expected%num(expected%id)), &
      & "Nuclear charges changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, all(abs(actual%xyz - expected%xyz) <= epsilon(1.0_wp)), &
      & "Nuclear coordinates changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, trim(actual%sym(1)) == "H" .and. trim(actual%sym(2)) == "O", &
      & "Nuclear labels changed in TREXIO round trip")
   call check(error, abs(actual%charge - expected%charge) <= epsilon(1.0_wp), &
      & "Molecular charge changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, actual%uhf, expected%uhf, &
      & "Molecular spin multiplicity changed in TREXIO round trip")
end subroutine check_structure

subroutine check_basis(error, actual, expected)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Read basis set information
   type(basis_type), intent(in) :: actual
   !> Reference basis set information
   type(basis_type), intent(in) :: expected

   call check(error, actual%nsh, expected%nsh, &
      & "Basis shell count changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, actual%nao, expected%nao, &
      & "Basis AO count changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, actual%nao_cart, expected%nao_cart, &
      & "Basis Cartesian AO count changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, actual%maxl, expected%maxl, &
      & "Basis maximum angular momentum changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, all(actual%sh2at == expected%sh2at), &
      & "Basis nucleus map changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, all(actual%nao_sh == expected%nao_sh), &
      & "Basis shell AO counts changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, all(actual%iao_sh == expected%iao_sh), &
      & "Basis shell AO offsets changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, all(actual%ao2sh == expected%ao2sh), &
      & "Basis AO shell map changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, all(actual%ao2at == expected%ao2at), &
      & "Basis AO atom map changed in TREXIO round trip")
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
      & "Wavefunction electron count changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, abs(actual%nuhf - expected%nuhf) <= epsilon(1.0_wp), &
      & "Wavefunction unpaired electron count changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, all(abs(actual%nel - expected%nel) <= epsilon(1.0_wp)), &
      & "Wavefunction spin electron counts changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, all(abs(actual%emo - expected%emo) <= epsilon(1.0_wp)), &
      & "Wavefunction MO energies changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, all(abs(actual%focc - expected%focc) <= epsilon(1.0_wp)), &
      & "Wavefunction MO occupations changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, all(abs(actual%coeff - expected%coeff) <= epsilon(1.0_wp)), &
      & "Wavefunction MO coefficients changed in TREXIO round trip")
end subroutine check_wavefunction

subroutine make_roundtrip_data(mol, bas, wfn)
   !> Molecular structure data
   type(structure_type), intent(out) :: mol
   !> Basis set information
   type(basis_type), intent(out) :: bas
   !> Wavefunction data
   type(wavefunction_type), intent(out) :: wfn

   integer, allocatable :: nshell(:)
   type(cgto_type), allocatable :: cgto(:, :)

   call new(mol, [1, 8], reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 1.4_wp], [3, 2]), charge=6.0_wp, uhf=1)

   allocate(nshell(mol%nid), cgto(1, mol%nid))
   nshell(:) = 1
   cgto(:, :) = cgto_type()

   cgto(1, 1)%ang = 0
   cgto(1, 1)%nprim = 1
   cgto(1, 1)%alpha(1) = 1.2_wp
   cgto(1, 1)%coeff(1) = 0.8_wp

   cgto(1, 2)%ang = 0
   cgto(1, 2)%nprim = 2
   cgto(1, 2)%alpha(:2) = [0.9_wp, 0.4_wp]
   cgto(1, 2)%coeff(:2) = [0.7_wp, 0.2_wp]

   call new_basis(bas, mol, nshell, cgto, 1.0_wp)

   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, 2, 0.0_wp)
   wfn%nocc = 3.0_wp
   wfn%nuhf = 1.0_wp
   wfn%nel = [2.0_wp, 1.0_wp]
   wfn%n0at = [1.0_wp, 8.0_wp]
   wfn%n0sh = [1.0_wp, 8.0_wp]
   wfn%emo(:, 1) = [-0.5_wp, 0.2_wp]
   wfn%emo(:, 2) = [-0.4_wp, 0.3_wp]
   wfn%focc(:, 1) = [1.5_wp, 0.5_wp]
   wfn%focc(:, 2) = [1.0_wp, 0.0_wp]
   wfn%coeff(:, :, 1) = reshape([1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [2, 2])
   wfn%coeff(:, :, 2) = reshape([0.6_wp, 0.8_wp, 0.8_wp, -0.6_wp], [2, 2])
end subroutine make_roundtrip_data

subroutine remove_trexio_output(filename)
   character(len=*), intent(in) :: filename

   integer :: io, stat
   logical :: exist

   inquire(file=filename, exist=exist)
   if (exist) then
      open(newunit=io, file=filename, status="old", action="write", iostat=stat)
      if (stat == 0) close(io, status="delete")
   end if
end subroutine remove_trexio_output

end module test_trexio
