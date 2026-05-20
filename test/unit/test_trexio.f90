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
   use tblite_basis_type, only : basis_type
   use tblite_features, only : get_tblite_feature
   use tblite_io_trexio, only : load_trexio, save_trexio, trexio_data_type
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
   type(error_type), allocatable, intent(out) :: error

   type(error_type), allocatable :: actual
   type(trexio_data_type) :: data

   call load_trexio(".trexio-disabled.trexio", data, actual)
   call check(error, allocated(actual), "TREXIO read unexpectedly succeeded without TREXIO support")
   if (allocated(actual)) deallocate(actual)
end subroutine test_trexio_disabled

subroutine test_roundtrip_text(error)
   type(error_type), allocatable, intent(out) :: error

   call check_roundtrip(".trexio-roundtrip.trexio", error)
end subroutine test_roundtrip_text

subroutine test_roundtrip_hdf5(error)
   type(error_type), allocatable, intent(out) :: error

   call check_roundtrip(".trexio-roundtrip.h5", error)
end subroutine test_roundtrip_hdf5

subroutine check_roundtrip(filename, error)
   character(len=*), intent(in) :: filename
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: energy = -12.345_wp
   type(structure_type) :: mol
   type(basis_type) :: bas
   type(wavefunction_type) :: wfn
   type(trexio_data_type) :: data

   call make_roundtrip_data(mol, bas, wfn)

   call remove_trexio_output(filename)
   call save_trexio(filename, mol, bas, wfn, energy, error)
   if (allocated(error)) return

   call load_trexio(filename, data, error)
   if (allocated(error)) return

   call check(error, data%nat, mol%nat)
   if (allocated(error)) return
   call check(error, all(abs(data%nuclear_charge - [1.0_wp, 8.0_wp]) <= epsilon(1.0_wp)), &
      & "Nuclear charges changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, all(abs(data%nuclear_coord - mol%xyz) <= epsilon(1.0_wp)), &
      & "Nuclear coordinates changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, trim(data%nuclear_label(1)) == "H" .and. trim(data%nuclear_label(2)) == "O", &
      & "Nuclear labels changed in TREXIO round trip")
   if (allocated(error)) return

   call check(error, data%electron_num, 3)
   if (allocated(error)) return
   call check(error, data%electron_up_num, 2)
   if (allocated(error)) return
   call check(error, data%electron_dn_num, 1)
   if (allocated(error)) return
   call check(error, abs(data%energy - energy) <= epsilon(1.0_wp), &
      & "State energy changed in TREXIO round trip")
   if (allocated(error)) return

   call check(error, data%nsh, bas%nsh)
   if (allocated(error)) return
   call check(error, data%nao, bas%nao)
   if (allocated(error)) return
   call check(error, all(data%basis_nucleus_index == bas%sh2at), &
      & "Basis nucleus map changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, all(data%basis_shell_ang_mom == [0, 0]), &
      & "Basis angular momenta changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, all(data%ao_shell == [1, 2]), "AO shell map changed in TREXIO round trip")
   if (allocated(error)) return

   call check(error, data%nmo, 4)
   if (allocated(error)) return
   call check(error, all(data%mo_spin == [0, 0, 1, 1]), "MO spins changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, all(abs(data%mo_energy - [-0.5_wp, 0.2_wp, -0.4_wp, 0.3_wp]) <= epsilon(1.0_wp)), &
      & "MO energies changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, all(abs(data%mo_occupation - [1.0_wp, 0.5_wp, 1.0_wp, 0.0_wp]) <= epsilon(1.0_wp)), &
      & "MO occupations changed in TREXIO round trip")
   if (allocated(error)) return
   call check(error, all(abs(data%mo_coefficient - reshape([ &
      & 1.0_wp, 0.0_wp, &
      & 0.6_wp, 0.8_wp, &
      & 0.0_wp, 1.0_wp, &
      & 0.8_wp, -0.6_wp], [4, 2])) <= epsilon(1.0_wp)), &
      & "MO coefficients changed in TREXIO round trip")
   if (allocated(error)) return

   call remove_trexio_output(filename)
end subroutine check_roundtrip

subroutine make_roundtrip_data(mol, bas, wfn)
   type(structure_type), intent(out) :: mol
   type(basis_type), intent(out) :: bas
   type(wavefunction_type), intent(out) :: wfn

   call new(mol, [1, 8], reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 1.4_wp], [3, 2]), charge=0.0_wp, uhf=1)

   bas%nsh = 2
   bas%nao = 2
   allocate(bas%nao_sh(2), bas%iao_sh(2), bas%sh2at(2))
   bas%nao_sh = [1, 1]
   bas%iao_sh = [0, 1]
   bas%sh2at = [1, 2]

   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, 2, 300.0_wp)
   wfn%nocc = 3.0_wp
   wfn%nuhf = 1.0_wp
   wfn%nel = [2.0_wp, 1.0_wp]
   wfn%emo(:, 1) = [-0.5_wp, 0.2_wp]
   wfn%emo(:, 2) = [-0.4_wp, 0.3_wp]
   wfn%focc(:, 1) = [1.0_wp, 0.5_wp]
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
