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

#ifndef TBLITE_HAS_HDF5
#define TBLITE_HAS_HDF5 0
#endif
#ifndef TBLITE_HAS_TREXIO
#define TBLITE_HAS_TREXIO 0
#endif

!> @file tblite/io/trexio.F90
!> Provides optional TREXIO output support.

!> TREXIO writer for tblite singlepoint results.
module tblite_io_trexio
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_wavefunction_type, only : wavefunction_type, get_alpha_beta_occupation
#if TBLITE_HAS_TREXIO
   use, intrinsic :: iso_c_binding, only : c_double
   use trexio, only : trexio_t, trexio_back_end_t, trexio_exit_code, &
      & TREXIO_SUCCESS, TREXIO_TEXT, TREXIO_HDF5, trexio_has_backend, &
      & trexio_open, trexio_close, trexio_string_of_error, &
      & trexio_write_nucleus_num, trexio_write_nucleus_charge, &
      & trexio_write_nucleus_coord, trexio_write_nucleus_label, &
      & trexio_write_electron_num, trexio_write_electron_up_num, &
      & trexio_write_electron_dn_num, trexio_write_state_num, &
      & trexio_write_state_id, trexio_write_state_energy, &
      & trexio_write_basis_type, trexio_write_basis_shell_num, &
      & trexio_write_basis_nucleus_index, trexio_write_basis_shell_ang_mom, &
      & trexio_write_ao_cartesian, trexio_write_ao_num, trexio_write_ao_shell, &
      & trexio_write_mo_type, trexio_write_mo_num, trexio_write_mo_coefficient, &
      & trexio_write_mo_occupation, trexio_write_mo_energy, trexio_write_mo_spin
#endif
   implicit none
   private

   public :: save_trexio

contains

!> Write tblite singlepoint data to a TREXIO file.
subroutine save_trexio(filename, mol, bas, wfn, energy, error)
   !> Output file or directory name
   character(len=*), intent(in) :: filename
   !> Molecular structure
   type(structure_type), intent(in) :: mol
   !> Basis metadata used by the calculation
   type(basis_type), intent(in) :: bas
   !> Converged wavefunction
   type(wavefunction_type), intent(in) :: wfn
   !> Total energy in atomic units
   real(wp), intent(in) :: energy
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

#if TBLITE_HAS_TREXIO
   integer(trexio_t) :: trex_file
   integer(trexio_back_end_t) :: backend
   integer(trexio_exit_code) :: rc

   call get_backend(filename, backend, error)
   if (allocated(error)) return

   if (.not.trexio_has_backend(backend)) then
      call fatal_error(error, "Requested TREXIO backend is not available in the linked TREXIO library")
      return
   end if

   trex_file = trexio_open(filename, 'w', backend, rc)
   if (rc /= TREXIO_SUCCESS) then
      call fatal_trexio(error, rc, "Failed to open TREXIO output '"//filename//"'")
      return
   end if

   call write_nucleus(trex_file, mol, error)
   if (.not.allocated(error)) call write_electrons(trex_file, wfn, error)
   if (.not.allocated(error)) call write_state(trex_file, energy, error)
   if (.not.allocated(error)) call write_basis(trex_file, bas, error)
   if (.not.allocated(error)) call write_mos(trex_file, wfn, error)

   rc = trexio_close(trex_file)
   if (rc /= TREXIO_SUCCESS .and. .not.allocated(error)) then
      call fatal_trexio(error, rc, "Failed to close TREXIO output '"//filename//"'")
   end if
#else
   call fatal_error(error, "TREXIO support is not available in this build of tblite.")
#endif
end subroutine save_trexio

#if TBLITE_HAS_TREXIO
subroutine get_backend(filename, backend, error)
   !> Output file or directory name
   character(len=*), intent(in) :: filename
   !> TREXIO backend selected from the output extension
   integer(trexio_back_end_t), intent(out) :: backend
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   if (is_hdf5_file(filename)) then
#if TBLITE_HAS_HDF5
      backend = TREXIO_HDF5
#else
      call fatal_error(error, "TREXIO HDF5 output requires HDF5 support in tblite.")
      backend = TREXIO_TEXT
#endif
   else
      backend = TREXIO_TEXT
   end if
end subroutine get_backend

subroutine write_nucleus(trex_file, mol, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Molecular structure
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(c_double), allocatable :: charge(:), coord(:, :)
   character(len=8), allocatable :: label(:)
   integer(trexio_exit_code) :: rc
   integer :: iat

   allocate(charge(mol%nat), coord(3, mol%nat), label(mol%nat))
   do iat = 1, mol%nat
      charge(iat) = real(mol%num(mol%id(iat)), c_double)
      label(iat) = mol%sym(mol%id(iat))
   end do
   coord = real(mol%xyz, c_double)

   rc = trexio_write_nucleus_num(trex_file, mol%nat)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO nucleus count")
   if (allocated(error)) return

   rc = trexio_write_nucleus_charge(trex_file, charge)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO nuclear charges")
   if (allocated(error)) return

   rc = trexio_write_nucleus_coord(trex_file, coord)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO nuclear coordinates")
   if (allocated(error)) return

   rc = trexio_write_nucleus_label(trex_file, label, len(label))
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO nuclear labels")
end subroutine write_nucleus

subroutine write_electrons(trex_file, wfn, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Converged wavefunction
   type(wavefunction_type), intent(in) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer(trexio_exit_code) :: rc
   integer :: electron_num, up_num, dn_num
   real(wp) :: nalp, nbet

   electron_num = nint(wfn%nocc)
   if (wfn%nspin == 2 .and. size(wfn%nel) >= 2) then
      up_num = nint(wfn%nel(1))
      dn_num = nint(wfn%nel(2))
   else
      call get_alpha_beta_occupation(wfn%nocc, wfn%nuhf, nalp, nbet)
      up_num = nint(nalp)
      dn_num = nint(nbet)
   end if

   rc = trexio_write_electron_num(trex_file, electron_num)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO electron count")
   if (allocated(error)) return

   rc = trexio_write_electron_up_num(trex_file, up_num)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO alpha electron count")
   if (allocated(error)) return

   rc = trexio_write_electron_dn_num(trex_file, dn_num)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO beta electron count")
end subroutine write_electrons

subroutine write_state(trex_file, energy, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Total energy in atomic units
   real(wp), intent(in) :: energy
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer(trexio_exit_code) :: rc

   rc = trexio_write_state_num(trex_file, 1)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO state count")
   if (allocated(error)) return

   rc = trexio_write_state_id(trex_file, 1)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO state id")
   if (allocated(error)) return

   rc = trexio_write_state_energy(trex_file, real(energy, c_double))
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO state energy")
end subroutine write_state

subroutine write_basis(trex_file, bas, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Basis metadata used by the calculation
   type(basis_type), intent(in) :: bas
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, allocatable :: shell_ang_mom(:), ao_shell(:)
   integer(trexio_exit_code) :: rc
   integer :: iao, ish, offset

   allocate(shell_ang_mom(bas%nsh), ao_shell(bas%nao))
   do ish = 1, bas%nsh
      shell_ang_mom(ish) = (bas%nao_sh(ish) - 1) / 2
      offset = bas%iao_sh(ish)
      do iao = 1, bas%nao_sh(ish)
         ao_shell(offset + iao) = ish
      end do
   end do

   rc = trexio_write_basis_type(trex_file, "Gaussian", 16)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO basis type")
   if (allocated(error)) return

   rc = trexio_write_basis_shell_num(trex_file, bas%nsh)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO shell count")
   if (allocated(error)) return

   rc = trexio_write_basis_nucleus_index(trex_file, bas%sh2at)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO shell atom map")
   if (allocated(error)) return

   rc = trexio_write_basis_shell_ang_mom(trex_file, shell_ang_mom)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO shell angular momenta")
   if (allocated(error)) return

   rc = trexio_write_ao_cartesian(trex_file, 0)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO AO cartesian flag")
   if (allocated(error)) return

   rc = trexio_write_ao_num(trex_file, bas%nao)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO AO count")
   if (allocated(error)) return

   rc = trexio_write_ao_shell(trex_file, ao_shell)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO AO shell map")
end subroutine write_basis

subroutine write_mos(trex_file, wfn, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Converged wavefunction
   type(wavefunction_type), intent(in) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(c_double), allocatable :: mo_coeff(:, :), mo_energy(:), mo_occ(:)
   integer, allocatable :: mo_spin(:)
   integer(trexio_exit_code) :: rc
   integer :: iao, imo, imoflat, ispin, nao, nspin

   nao = size(wfn%coeff, 1)
   nspin = size(wfn%coeff, 3)
   allocate(mo_coeff(nao*nspin, nao), mo_energy(nao*nspin), mo_occ(nao*nspin), mo_spin(nao*nspin))

   imoflat = 0
   do ispin = 1, nspin
      do imo = 1, nao
         imoflat = imoflat + 1
         mo_energy(imoflat) = real(wfn%emo(imo, ispin), c_double)
         mo_occ(imoflat) = real(wfn%focc(imo, min(ispin, size(wfn%focc, 2))), c_double)
         mo_spin(imoflat) = ispin - 1
         do iao = 1, nao
            mo_coeff(imoflat, iao) = real(wfn%coeff(iao, imo, ispin), c_double)
         end do
      end do
   end do

   rc = trexio_write_mo_type(trex_file, "Tight-binding", 32)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO MO type")
   if (allocated(error)) return

   rc = trexio_write_mo_num(trex_file, nao*nspin)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO MO count")
   if (allocated(error)) return

   rc = trexio_write_mo_coefficient(trex_file, mo_coeff)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO MO coefficients")
   if (allocated(error)) return

   rc = trexio_write_mo_occupation(trex_file, mo_occ)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO MO occupations")
   if (allocated(error)) return

   rc = trexio_write_mo_energy(trex_file, mo_energy)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO MO energies")
   if (allocated(error)) return

   rc = trexio_write_mo_spin(trex_file, mo_spin)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO MO spins")
end subroutine write_mos

subroutine fatal_trexio(error, rc, context)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> TREXIO return code
   integer(trexio_exit_code), intent(in) :: rc
   !> Operation that failed
   character(len=*), intent(in) :: context

   character(len=256) :: message

   call trexio_string_of_error(rc, message)
   call fatal_error(error, context//": "//trim(message))
end subroutine fatal_trexio
#endif

pure function is_hdf5_file(filename) result(is_hdf5)
   !> File name to inspect
   character(len=*), intent(in) :: filename
   !> Whether the file name uses an HDF5 extension
   logical :: is_hdf5

   is_hdf5 = len(filename) >= 3 .and. filename(len(filename)-2:) == ".h5"
   is_hdf5 = is_hdf5 .or. &
      & (len(filename) >= 5 .and. filename(len(filename)-4:) == ".hdf5")
end function is_hdf5_file

end module tblite_io_trexio