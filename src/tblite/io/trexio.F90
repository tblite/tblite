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

#define trexio tblite_ext_trexio
#if TBLITE_HAS_TREXIO
#include <trexio_f.f90>
#endif
#undef trexio

!> TREXIO writer for tblite singlepoint results.
module tblite_io_trexio
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type, new
   use tblite_basis_type, only : basis_type, cgto_type, new_basis, new_cgto
   use tblite_integral_trafo, only : transform0
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction, &
      & get_alpha_beta_occupation, get_density_matrix
#if TBLITE_HAS_TREXIO
   use, intrinsic :: iso_c_binding, only : c_double
   use tblite_ext_trexio, only : trexio_t, trexio_back_end_t, trexio_exit_code, &
      & TREXIO_SUCCESS, TREXIO_TEXT, TREXIO_HDF5, trexio_has_backend, &
      & trexio_open, trexio_close, trexio_string_of_error, &
      & trexio_read_nucleus_num, trexio_read_nucleus_charge, &
      & trexio_read_nucleus_coord, trexio_read_nucleus_label, &
      & trexio_read_electron_num, trexio_read_electron_up_num, &
      & trexio_read_electron_dn_num, trexio_read_state_energy, &
      & trexio_read_basis_shell_num, trexio_read_basis_nucleus_index, &
      & trexio_read_basis_shell_ang_mom, trexio_read_ecp_z_core, &
      & trexio_read_ao_num, trexio_read_ao_cartesian, trexio_read_ao_shell, &
      & trexio_read_mo_num, trexio_read_mo_coefficient, &
      & trexio_read_mo_occupation, trexio_read_mo_energy, trexio_read_mo_spin, &
      & trexio_has_mo_spin, trexio_read_basis_prim_num, trexio_read_basis_exponent, &
      & trexio_read_basis_coefficient, trexio_read_basis_shell_factor, &
      & trexio_read_basis_shell_index, trexio_read_basis_prim_factor, &
      & trexio_write_nucleus_num, trexio_write_nucleus_charge, &
      & trexio_write_nucleus_coord, trexio_write_nucleus_label, &
      & trexio_write_electron_num, trexio_write_electron_up_num, &
      & trexio_write_electron_dn_num, trexio_write_state_num, &
      & trexio_write_state_id, trexio_write_state_energy, &
      & trexio_write_basis_type, trexio_write_basis_shell_num, &
      & trexio_write_basis_nucleus_index, trexio_write_basis_shell_ang_mom, &
      & trexio_write_ecp_z_core, trexio_write_ao_cartesian, trexio_write_ao_num, &
      & trexio_write_ao_shell, trexio_write_mo_type, trexio_write_mo_num, &
      & trexio_write_mo_coefficient, trexio_write_mo_occupation, &
      & trexio_write_mo_energy, trexio_write_mo_spin, &
      & trexio_write_basis_prim_num, trexio_write_basis_exponent, &
      & trexio_write_basis_coefficient, trexio_write_basis_shell_factor, &
      & trexio_write_basis_shell_index, trexio_write_basis_prim_factor, &
      & trexio_write_cell_a, trexio_write_cell_b, trexio_write_cell_c, &
      & trexio_read_cell_a, trexio_read_cell_b, trexio_read_cell_c
#endif
   implicit none
   private

   public :: load_trexio, save_trexio

   !> Maximum contraction length allowed in TREXIO basis sets
   integer, parameter :: maxg = 12

contains

!> Read tblite TREXIO output into tblite objects
subroutine load_trexio(filename, mol, bas, wfn, energy, error)
   !> Input file or directory name
   character(len=*), intent(in) :: filename
   !> Molecular structure data loaded from TREXIO
   type(structure_type), intent(out) :: mol
   !> Basis set information loaded from TREXIO
   type(basis_type), intent(out) :: bas
   !> Wavefunction data loaded from TREXIO
   type(wavefunction_type), intent(out) :: wfn
   !> Total energy
   real(wp), intent(out) :: energy
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

   trex_file = trexio_open(filename, "r", backend, rc)
   if (rc /= TREXIO_SUCCESS) then
      call fatal_trexio(error, rc, "Failed to open TREXIO input '"//filename//"'")
      return
   end if

   call load_structure(trex_file, mol, error)
   if (.not.allocated(error)) call read_state(trex_file, energy, error)
   if (.not.allocated(error)) call load_basis(trex_file, mol, bas, error)
   if (.not.allocated(error)) call load_wavefunction(trex_file, mol, bas, wfn, error)

   rc = trexio_close(trex_file)
   if (rc /= TREXIO_SUCCESS .and. .not.allocated(error)) then
      call fatal_trexio(error, rc, "Failed to close TREXIO input '"//filename//"'")
   end if
#else
   call fatal_error(error, "TREXIO support is not available in this build of tblite.")
#endif
end subroutine load_trexio

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

   trex_file = trexio_open(filename, "w", backend, rc)
   if (rc /= TREXIO_SUCCESS) then
      call fatal_trexio(error, rc, "Failed to open TREXIO output '"//filename//"'")
      return
   end if

   call write_nucleus(trex_file, mol, error)
   if (.not.allocated(error)) call write_cell(trex_file, mol, error)
   if (.not.allocated(error)) call write_electron(trex_file, mol, wfn, error)
   if (.not.allocated(error)) call write_state(trex_file, energy, error)
   if (.not.allocated(error)) call write_basis(trex_file, mol, bas, error)
   if (.not.allocated(error)) call write_ecp(trex_file, mol, wfn, error)
   if (.not.allocated(error)) call write_ao(trex_file, bas, error)
   if (.not.allocated(error)) call write_mo(trex_file, mol, bas, wfn, error)

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

!> Load molecular structure object from TREXIO file
subroutine load_structure(trex_file, mol, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Molecular structure data loaded from TREXIO
   type(structure_type), intent(out) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: nat, electron_num, electron_up_num, electron_dn_num
   real(wp), allocatable :: nuclear_charge(:), nuclear_coord(:, :), lattice(:, :)
   character(len=8), allocatable :: nuclear_label(:)
   logical, allocatable :: periodic(:)

   integer :: uhf
   real(wp) :: charge

   ! Read structure data from TREXIO file
   call read_nucleus(trex_file, nat, nuclear_charge, nuclear_coord, nuclear_label, &
      & error)
   if (.not.allocated(error)) call read_cell(trex_file, periodic, lattice, error)
   if (.not.allocated(error)) call read_electron(trex_file, electron_num, &
      & electron_up_num, electron_dn_num, error)
   if (allocated(error)) return

   ! Total charge and multiplicity
   charge = sum(nuclear_charge) - real(electron_num, wp)
   uhf = electron_up_num - electron_dn_num

   ! Set up new molecular structure object
   call new(mol, nint(nuclear_charge), nuclear_label, nuclear_coord, charge=charge, &
      & uhf=uhf, lattice=lattice, periodic=periodic)
end subroutine load_structure

!> Read total nuclear charge, atomic symbols and coordinates
subroutine read_nucleus(trex_file, nat, nuclear_charge, nuclear_coord, &
   & nuclear_label, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Number of atoms
   integer, intent(out) :: nat
   !> Nuclear charges
   real(wp), allocatable, intent(out) :: nuclear_charge(:)
   !> Nuclear coordinates in Bohr
   real(wp), allocatable, intent(out) :: nuclear_coord(:, :)
   !> Nuclear labels
   character(len=8), allocatable, intent(out) :: nuclear_label(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(c_double), allocatable :: charge(:), coord(:, :)
   integer(trexio_exit_code) :: rc

   rc = trexio_read_nucleus_num(trex_file, nat)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO nucleus count")
   if (allocated(error)) return

   allocate(charge(nat), coord(3, nat), nuclear_label(nat))
   rc = trexio_read_nucleus_charge(trex_file, charge)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO nuclear charges")
   if (allocated(error)) return

   rc = trexio_read_nucleus_coord(trex_file, coord)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO nuclear coordinates")
   if (allocated(error)) return

   rc = trexio_read_nucleus_label(trex_file, nuclear_label, len(nuclear_label))
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO nuclear labels")
   if (allocated(error)) return

   nuclear_charge = real(charge, wp)
   nuclear_coord = real(coord, wp)
end subroutine read_nucleus

!> Read lattice vectors and periodicity information
subroutine read_cell(trex_file, periodic, lattice, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Periodic directions
   logical, allocatable, intent(out) :: periodic(:)
   !> Lattice vectors in Bohr
   real(wp), allocatable, intent(out) :: lattice(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(c_double) :: a(3), b(3), c(3)
   logical :: tmp_periodic(3)
   integer(trexio_exit_code) :: rc

   tmp_periodic(:) = .true.
   rc = trexio_read_cell_a(trex_file, a)
   if (rc /= TREXIO_SUCCESS) tmp_periodic(1) = .false.

   rc = trexio_read_cell_b(trex_file, b)
   if (rc /= TREXIO_SUCCESS) tmp_periodic(2) = .false.

   rc = trexio_read_cell_c(trex_file, c)
   if (rc /= TREXIO_SUCCESS) tmp_periodic(3) = .false.

   if (any(tmp_periodic)) then
      periodic = tmp_periodic

      allocate(lattice(3, 3), source=0.0_wp)
      if (periodic(1)) lattice(:, 1) = real(a, wp)
      if (periodic(2)) lattice(:, 2) = real(b, wp)
      if (periodic(3)) lattice(:, 3) = real(c, wp)
   end if
end subroutine read_cell

!> Read total number of alpha/beta electrons including the core electrons
subroutine read_electron(trex_file, electron_num, electron_up_num, electron_dn_num, &
   & error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Number of electrons
   integer, intent(out) :: electron_num
   !> Number of alpha electrons
   integer, intent(out) :: electron_up_num
   !> Number of beta electrons
   integer, intent(out) :: electron_dn_num
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer(trexio_exit_code) :: rc

   rc = trexio_read_electron_num(trex_file, electron_num)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO electron count")
   if (allocated(error)) return

   rc = trexio_read_electron_up_num(trex_file, electron_up_num)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO alpha electron count")
   if (allocated(error)) return

   rc = trexio_read_electron_dn_num(trex_file, electron_dn_num)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO beta electron count")
end subroutine read_electron

!> Read total energy of the state in atomic units
subroutine read_state(trex_file, energy, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Total energy in atomic units
   real(wp), intent(out) :: energy
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(c_double) :: state_energy
   integer(trexio_exit_code) :: rc

   rc = trexio_read_state_energy(trex_file, state_energy)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO state energy")
   if (allocated(error)) return

   energy = real(state_energy, wp)
end subroutine read_state


!> Load basis set information object from TREXIO file
subroutine load_basis(trex_file, mol, bas, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Molecular structure data loaded from TREXIO
   type(structure_type), intent(in) :: mol
   !> Basis set information loaded from TREXIO
   type(basis_type), intent(out) :: bas
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: nsh, nprim
   integer, allocatable :: basis_nucleus_index(:), basis_shell_index(:)
   integer, allocatable :: basis_shell_ang_mom(:)
   real(wp), allocatable :: basis_shell_factor(:), basis_primitive_factor(:)
   real(wp), allocatable :: basis_primitive_alpha(:), basis_primitive_coeff(:)

   type(cgto_type), allocatable :: cgto(:, :)
   type(cgto_type) :: tmp_cgto
   logical, allocatable :: seen_cgto(:, :)
   integer, allocatable :: nsh_id(:), nsh_at(:), sh2at(:), nsh_work(:)
   real(wp), allocatable :: alpha(:), coeff(:)
   integer :: iat, iprim, ish, isp, jsh, l, ng

   call read_basis(trex_file, nsh, nprim, basis_nucleus_index, basis_shell_ang_mom, &
      & basis_shell_factor, basis_shell_index, basis_primitive_alpha, &
      & basis_primitive_coeff, basis_primitive_factor, error)
   if (allocated(error)) return

   if (any(basis_shell_index < 1) .or. any(basis_shell_index > nsh)) then
      call fatal_error(error, "TREXIO primitive references an unknown shell")
      return
   end if

   ! Shell count per atom and shell to atom mapping
   allocate(nsh_at(mol%nat), sh2at(nsh), source=0)
   do ish = 1, nsh
      iat = basis_nucleus_index(ish)
      if (iat < 1 .or. iat > mol%nat) then
         call fatal_error(error, "TREXIO shell references an unknown atom")
         return
      end if
      sh2at(ish) = iat
      nsh_at(iat) = nsh_at(iat) + 1
   end do

   ! Shell count per species requiring identical shell counts for all atoms of a species
   allocate(nsh_id(mol%nid))
   nsh_id(:) = -1
   do iat = 1, mol%nat
      isp = mol%id(iat)
      if (nsh_id(isp) < 0) then
         nsh_id(isp) = nsh_at(iat)
      else if (nsh_id(isp) /= nsh_at(iat)) then
         call fatal_error(error, "TREXIO basis has inconsistent shell counts for one species")
         return
      end if
   end do
   if (any(nsh_id <= 0)) then
      call fatal_error(error, "TREXIO missing shells for at least one species")
      return
   end if

   ! Build CGTOs for each species
   allocate(cgto(maxval(nsh_id), mol%nid))
   allocate(seen_cgto(maxval(nsh_id), mol%nid), source=.false.)
   allocate(alpha(maxg), coeff(maxg))
   allocate(nsh_work(mol%nat), source=0)
   do ish = 1, nsh
      iat = sh2at(ish)
      isp = mol%id(iat)
      ! Count shells at the current atom
      nsh_work(iat) = nsh_work(iat) + 1
      jsh = nsh_work(iat)
      ! Shell angular momentum
      l = basis_shell_ang_mom(ish)

      alpha(:) = 0.0_wp
      coeff(:) = 0.0_wp
      ng = 0
      do iprim = 1, nprim
         if (basis_shell_index(iprim) /= ish) cycle
         ng = ng + 1
         alpha(ng) = basis_primitive_alpha(iprim)
         coeff(ng) = basis_shell_factor(ish) * basis_primitive_factor(iprim) &
            & * basis_primitive_coeff(iprim)
      end do
      if (ng <= 0 .or. ng > maxg) then
         call fatal_error(error, "TREXIO shell has no primitives or more primitives than supported")
         return
      end if
      ! Apply no additional normalization, already part of the coefficents and primitive factor
      call new_cgto(tmp_cgto, ng, l, alpha, coeff, .false.)

      ! Check if current CGTO matches earlier CGTOs for the same species
      if (.not. seen_cgto(jsh, isp)) then
         cgto(jsh, isp) = tmp_cgto
         seen_cgto(jsh, isp) = .true.
      else if (.not.cgto(jsh, isp)%compare(tmp_cgto)) then
         call fatal_error(error, "TREXIO basis differs between atoms of the same species")
         return
      end if
   end do

   call new_basis(bas, mol, nsh_id, cgto, 1.0_wp)
end subroutine load_basis

!> Read shell mappings, angular momenta, contraction coefficients and exponents
subroutine read_basis(trex_file, nsh, nprim, basis_nucleus_index, basis_shell_ang_mom, &
   & basis_shell_factor, basis_shell_index, basis_primitive_alpha, &
   & basis_primitive_coeff, basis_primitive_factor, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Number of shells
   integer, intent(out) :: nsh
   !> Number of primitives
   integer, intent(out) :: nprim
   !> Mapping from shell index to nucleus index
   integer, allocatable, intent(out) :: basis_nucleus_index(:)
   !> Angular momentum of each shell
   integer, allocatable, intent(out) :: basis_shell_ang_mom(:)
   !> Shell normalization factors
   real(wp), allocatable, intent(out) :: basis_shell_factor(:)
   !> Mapping from primitive index to shell index
   integer, allocatable, intent(out) :: basis_shell_index(:)
   !> Exponent of Gaussian primitive functions
   real(wp), allocatable, intent(out) :: basis_primitive_alpha(:)
   !> Contraction coefficient of each primitive
   real(wp), allocatable, intent(out) :: basis_primitive_coeff(:)
   !> Primitive normalization factors
   real(wp), allocatable, intent(out) :: basis_primitive_factor(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer(trexio_exit_code) :: rc
   real(c_double), allocatable :: d_shell_factor(:), d_prim_factor(:)
   real(c_double), allocatable :: d_alpha(:), d_coeff(:)

   rc = trexio_read_basis_shell_num(trex_file, nsh)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO shell count")
   if (allocated(error)) return

   allocate(basis_nucleus_index(nsh))
   rc = trexio_read_basis_nucleus_index(trex_file, basis_nucleus_index)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO shell atom map")
   if (allocated(error)) return

   allocate(basis_shell_ang_mom(nsh))
   rc = trexio_read_basis_shell_ang_mom(trex_file, basis_shell_ang_mom)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO shell angular momenta")
   if (allocated(error)) return

   allocate(d_shell_factor(nsh), basis_shell_factor(nsh))
   rc = trexio_read_basis_shell_factor(trex_file, d_shell_factor)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO shell factor")
   if (allocated(error)) return
   basis_shell_factor = real(d_shell_factor, wp)

   rc = trexio_read_basis_prim_num(trex_file, nprim)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO primitive count")
   if (allocated(error)) return

   allocate(basis_shell_index(nprim))
   rc = trexio_read_basis_shell_index(trex_file, basis_shell_index)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO primitive shell map")
   if (allocated(error)) return

   allocate(d_alpha(nprim), basis_primitive_alpha(nprim))
   rc = trexio_read_basis_exponent(trex_file, d_alpha)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO primitive exponent")
   if (allocated(error)) return
   basis_primitive_alpha = real(d_alpha, wp)

   allocate(d_coeff(nprim), basis_primitive_coeff(nprim))
   rc = trexio_read_basis_coefficient(trex_file, d_coeff)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO primitive coefficient")
   if (allocated(error)) return
   basis_primitive_coeff = real(d_coeff, wp)

   allocate(d_prim_factor(nprim), basis_primitive_factor(nprim))
   rc = trexio_read_basis_prim_factor(trex_file, d_prim_factor)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO primitive factor")
   if (allocated(error)) return
   basis_primitive_factor = real(d_prim_factor, wp)
end subroutine read_basis


subroutine load_wavefunction(trex_file, mol, bas, wfn, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Molecular structure data loaded from TREXIO
   type(structure_type), intent(in) :: mol
   !> Basis set information loaded from TREXIO
   type(basis_type), intent(in) :: bas
   !> Wavefunction data loaded from TREXIO
   type(wavefunction_type), intent(out) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: nao, nmo
   logical :: cartesian
   real(wp), allocatable :: mo_occupation(:), mo_energy(:), mo_coefficient(:, :)
   integer, allocatable :: ecp_z_core(:), ao_shell(:), mo_spin(:)

   integer :: nspin, spin, imo, ish, iao, jsh, jat, jsp, js, jj, lj, imoflat
   integer :: ncart, icart, isphr, nsphr
   integer, allocatable :: ao_pos(:, :), ao_count(:), ao_expect(:), spin_count(:)
   integer :: perm_cart(15), perm_sphr(9)
   real(wp) :: ccart(15, 1), csphr(9, 1)
   real(wp), allocatable :: tmp_focc(:)

   call read_ecp(trex_file, mol%nat, ecp_z_core, error)
   if (.not.allocated(error)) call read_ao(trex_file, cartesian, nao, ao_shell, error)
   if (.not.allocated(error)) call read_mo(trex_file, nao, nmo, mo_coefficient, &
      & mo_occupation, mo_energy, mo_spin, error)
   if (allocated(error)) return

   if (cartesian) then
      if (nao /= bas%nao_cart) then
         call fatal_error(error, "TREXIO cartesian coefficient count does not match the basis")
         return
      end if
   else
      if (nao /= bas%nao) then
         call fatal_error(error, "TREXIO spherical coefficient count does not match the basis")
         return
      end if
   end if
   if (nmo == bas%nao) then
      ! Restricted TREXIO representation with one spatial MO list
      nspin = 1
   else if (nmo == 2 * bas%nao) then
      ! UHF TREXIO representation with alpha and beta spin-orbitals
      nspin = 2
      if (.not. any(mo_spin == 0) .or. .not. any(mo_spin == 1)) then
         call fatal_error(error, "TREXIO UHF MO data is missing alpha or beta spin labels")
         return
      end if
   else
      call fatal_error(error, "TREXIO MO count does not match the spherical AO basis")
      return
   end if

   ! Construct and check the AO to shell mapping
   allocate(ao_expect(bas%nsh))
   if (cartesian) then
      ao_expect(:) = bas%nao_cart_sh
   else
      ao_expect(:) = bas%nao_sh
   end if
   allocate(ao_pos(maxval(ao_expect), bas%nsh), ao_count(bas%nsh))
   ao_pos(:, :) = 0
   ao_count(:) = 0
   do iao = 1, nao
      ish = ao_shell(iao)
      if (ish < 1 .or. ish > bas%nsh) then
         call fatal_error(error, "TREXIO shell map references an unknown shell")
         return
      end if

      ao_count(ish) = ao_count(ish) + 1
      ao_pos(ao_count(ish), ish) = iao
   end do
   if (any(ao_count /= ao_expect)) then
      call fatal_error(error, "TREXIO shell map does not match the basis set")
      return
   end if

   ! Set up new wavefunction object (by default 0K electronic temperature)
   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, nspin, 0.0_wp)

   ! Set MO occupation and unpaired electrons
   wfn%nocc = sum(mo_occupation)
   wfn%nuhf = real(mol%uhf, wp)
   if (nspin == 1) then
      call get_alpha_beta_occupation(wfn%nocc, wfn%nuhf, wfn%nel(1), wfn%nel(2))
   else
      wfn%nel(1) = sum(mo_occupation, mask=mo_spin == 0)
      wfn%nel(2) = sum(mo_occupation, mask=mo_spin == 1)
   end if
   if (abs(wfn%nuhf - abs(wfn%nel(1) - wfn%nel(2))) > epsilon(1.0_wp)) then
      call fatal_error(error, "TREXIO occupation data has an inconsistent number of unpaired electrons")
      return
   end if

   ! Set atomic reference occupation from ECP core orbitals
   wfn%n0at(:) = real(mol%num(mol%id) - ecp_z_core, wp)

   ! MO coefficients, energy, and occupation numbers
   allocate(spin_count(nspin), source=0)
   wfn%coeff(:, :, :) = 0.0_wp
   do imoflat = 1, nmo
      spin = min(max(mo_spin(imoflat) + 1, 1), nspin)
      spin_count(spin) = spin_count(spin) + 1
      imo = spin_count(spin)
      if (imo > bas%nao) then
         call fatal_error(error, "TREXIO contains too many orbitals for one spin channel")
         return
      end if

      wfn%emo(imo, spin) = mo_energy(imoflat)
      if (nspin == 1) then
         ! Distribute restricted occupation equally between alpha and beta
         wfn%focc(imo, 1) = 0.5_wp * mo_occupation(imoflat)
         wfn%focc(imo, 2) = 0.5_wp * mo_occupation(imoflat)
      else
         wfn%focc(imo, spin) = mo_occupation(imoflat)
      end if

      do jsh = 1, bas%nsh
         jat = bas%sh2at(jsh)
         jsp = mol%id(jat)
         js = bas%ish_at(jat)
         lj = bas%cgto(jsh-js, jsp)%ang
         jj = bas%iao_sh(jsh)

         nsphr = bas%nao_sh(jsh)
         if (cartesian) then
            ! Reorder MO coefficients to match Cartesian AO ordering and transform
            call get_trexio_to_tblite_cart_perm(lj, perm_cart, error)
            if (allocated(error)) return

            ncart = bas%nao_cart_sh(jsh)
            ccart(:, :) = 0.0_wp
            csphr(:, :) = 0.0_wp
            do icart = 1, ncart
               iao = ao_pos(icart, jsh)
               ccart(perm_cart(icart), 1) = mo_coefficient(iao, imoflat)
            end do

            call transform0(lj, 0, ccart(:ncart, :1), csphr(:nsphr, :1), &
               & .true., .false.)
            wfn%coeff(jj+1:jj+nsphr, imo, spin) = csphr(:nsphr, 1)
         else
            ! Reorder MO coefficients to match spherical AO ordering
            call get_trexio_to_tblite_sphr_perm(lj, perm_sphr, error)
            if (allocated(error)) return

            do isphr = 1, nsphr
               iao = ao_pos(isphr, jsh)
               wfn%coeff(jj + perm_sphr(isphr), imo, spin) = &
                  & mo_coefficient(iao, imoflat)
            end do
         end if
      end do
   end do
   if (any(spin_count /= bas%nao)) then
      call fatal_error(error, "TREXIO data has an inconsistent number of spin orbitals")
      return
   end if

   ! Calculate the density matrix from MO coefficients and occupation numbers
   if (nspin == 1) then
      allocate(tmp_focc(size(wfn%focc, 1)))
      tmp_focc(:) = wfn%focc(:, 1) + wfn%focc(:, 2)
      call get_density_matrix(tmp_focc, wfn%coeff(:, :, 1), wfn%density(:, :, 1))
   else
      do spin = 1, nspin
         call get_density_matrix(wfn%focc(:, spin), wfn%coeff(:, :, spin), &
            & wfn%density(:, :, spin))
      end do
   end if
end subroutine load_wavefunction

subroutine read_ecp(trex_file, nat, ecp_z_core, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Number of atoms
   integer, intent(in) :: nat
   !> Number of core electrons removed per atom
   integer, allocatable, intent(out) :: ecp_z_core(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer(trexio_exit_code) :: rc

   allocate(ecp_z_core(nat))
   rc = trexio_read_ecp_z_core(trex_file, ecp_z_core)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO Z ECP core")
end subroutine read_ecp

subroutine read_ao(trex_file, cartesian, nao, ao_shell, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Whether the AOs are in Cartesian (true) or spherical (false) form
   logical, intent(out) :: cartesian
   !> Number of atomic orbitals
   integer, intent(out) :: nao
   !> Shell index for each AO
   integer, allocatable, intent(out) :: ao_shell(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer(trexio_exit_code) :: rc
   integer :: cart

   rc = trexio_read_ao_cartesian(trex_file, cart)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO AO cartesian flag")
   if (allocated(error)) return
   cartesian = cart /= 0

   rc = trexio_read_ao_num(trex_file, nao)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO AO count")
   if (allocated(error)) return

   allocate(ao_shell(nao))
   rc = trexio_read_ao_shell(trex_file, ao_shell)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO AO shell map")
end subroutine read_ao

subroutine read_mo(trex_file, nao, nmo, mo_coefficient, mo_occupation, mo_energy, &
   & mo_spin, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Number of atomic orbitals
   integer, intent(in) :: nao
   !> Number of molecular orbitals
   integer, intent(out) :: nmo
   !> Molecular orbital coefficients
   real(wp), allocatable, intent(out) :: mo_coefficient(:, :)
   !> Molecular orbital occupations
   real(wp), allocatable, intent(out) :: mo_occupation(:)
   !> Molecular orbital energies
   real(wp), allocatable, intent(out) :: mo_energy(:)
   !> Molecular orbital spins
   integer, allocatable, intent(out) :: mo_spin(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(c_double), allocatable :: mo_coeff(:, :), mo_ene(:), mo_occ(:)
   integer(trexio_exit_code) :: rc
   integer :: imo

   rc = trexio_read_mo_num(trex_file, nmo)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO MO count")
   if (allocated(error)) return

   allocate(mo_coeff(nao, nmo), mo_coefficient(nao, nmo))
   rc = trexio_read_mo_coefficient(trex_file, mo_coeff)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO MO coefficients")
   if (allocated(error)) return

   allocate(mo_occ(nmo), mo_occupation(nmo))
   rc = trexio_read_mo_occupation(trex_file, mo_occ)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO MO occupations")
   if (allocated(error)) return

   allocate(mo_ene(nmo), mo_energy(nmo))
   rc = trexio_read_mo_energy(trex_file, mo_ene)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO MO energies")
   if (allocated(error)) return

   allocate(mo_spin(nmo), source=0)
   rc = trexio_has_mo_spin(trex_file)
   if (rc == TREXIO_SUCCESS) then
      rc = trexio_read_mo_spin(trex_file, mo_spin)
      if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to read TREXIO MO spins")
      if (allocated(error)) return
   end if

   do imo = 1, nmo
      mo_coefficient(:, imo) = real(mo_coeff(:, imo), wp)
   end do
   mo_occupation = real(mo_occ, wp)
   mo_energy = real(mo_ene, wp)
end subroutine read_mo

!> Spherical AO order permutation from TREXIO (0, ... l, -l) to tblite (-l, ..., 0, ..., l)
subroutine get_trexio_to_tblite_sphr_perm(l, perm, error)
   !> Angular momentum quantum number
   integer, intent(in)  :: l
   !> Permutation array
   integer, intent(out) :: perm(9)
   !> Error handling
   type(error_type), intent(out), allocatable :: error

   perm = 0
   select case(l)
   case(0)
      perm(1) = 1
   case(1)
      perm(1:3) = [2, 3, 1]
   case(2)
      perm(1:5) = [3, 4, 2, 5, 1]
   case(3)
      perm(1:7) = [4, 5, 3, 6, 2, 7, 1]
   case(4)
      perm(1:9) = [5, 6, 4, 7, 3, 8, 2, 9, 1]
   case default
      call fatal_error(error, "TREXIO reader only supports angular momenta up to g")
      return
   end select
end subroutine get_trexio_to_tblite_sphr_perm

!> Cartesian AO order permutation from TREXIO (alphabetical) to tblite
subroutine get_trexio_to_tblite_cart_perm(l, perm, error)
   !> Angular momentum quantum number
   integer, intent(in) :: l
   !> Permutation array
   integer, intent(out) :: perm(15)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   perm = 0
   select case(l)
   case(0)
      perm(1) = 1
   case(1)
      perm(1:3) = [1, 2, 3]
   case(2)
      perm(1:6) = [1, 2, 3, 4, 5, 6]
   case(3)
      perm(1:10) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
   case(4)
      perm(1:15) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
   case default
      call fatal_error(error, "TREXIO reader only supports angular momenta up to g")
      return
   end select
end subroutine get_trexio_to_tblite_cart_perm


subroutine write_nucleus(trex_file, mol, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(c_double), allocatable :: charge(:), coord(:, :)
   character(len=8), allocatable :: label(:)
   integer(trexio_exit_code) :: rc
   integer :: iat, isp

   allocate(charge(mol%nat), coord(3, mol%nat), label(mol%nat))
   do iat = 1, mol%nat
      isp = mol%id(iat)
      charge(iat) = real(mol%num(isp), c_double)
      label(iat) = mol%sym(isp)
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

subroutine write_cell(trex_file, mol, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(c_double) :: a(3), b(3), c(3)
   integer(trexio_exit_code) :: rc

   if (.not. allocated(mol%periodic) .or. .not. allocated(mol%lattice)) return
   if (.not. any(mol%periodic)) return

   if (mol%periodic(1)) then
      a = real(mol%lattice(:, 1), c_double)
      rc = trexio_write_cell_a(trex_file, a)
      if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO cell_a")
      if (allocated(error)) return
   end if

   if (mol%periodic(2)) then
      b = real(mol%lattice(:, 2), c_double)
      rc = trexio_write_cell_b(trex_file, b)
      if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO cell_b")
      if (allocated(error)) return
   end if

   if (mol%periodic(3)) then
      c = real(mol%lattice(:, 3), c_double)
      rc = trexio_write_cell_c(trex_file, c)
      if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO cell_c")
      if (allocated(error)) return
   end if
end subroutine write_cell

subroutine write_electron(trex_file, mol, wfn, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Converged wavefunction
   type(wavefunction_type), intent(in) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer(trexio_exit_code) :: rc
   integer :: core_num, core_up_num, core_dn_num, electron_num, up_num, dn_num
   real(wp) :: nalp, nbet

   ! Add back the neglected core electrons
   core_num = nint(sum(mol%num(mol%id) - wfn%n0at))
   core_up_num = core_num / 2
   core_dn_num = core_num / 2

   electron_num = nint(wfn%nocc) + core_num
   if (wfn%nspin == 2 .and. size(wfn%nel) >= 2) then
      up_num = nint(wfn%nel(1)) + core_up_num
      dn_num = nint(wfn%nel(2)) + core_dn_num
   else
      call get_alpha_beta_occupation(wfn%nocc, wfn%nuhf, nalp, nbet)
      up_num = nint(nalp) + core_up_num
      dn_num = nint(nbet) + core_dn_num
   end if

   rc = trexio_write_electron_num(trex_file, electron_num)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO electron count")
   if (allocated(error)) return

   rc = trexio_write_electron_up_num(trex_file, up_num)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO alpha electron count")
   if (allocated(error)) return

   rc = trexio_write_electron_dn_num(trex_file, dn_num)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO beta electron count")
end subroutine write_electron

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

subroutine write_basis(trex_file, mol, bas, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis metadata used by the calculation
   type(basis_type), intent(in) :: bas
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer(trexio_exit_code) :: rc
   integer, allocatable :: shell_ang_mom(:), basis_shell_index(:)
   real(c_double), allocatable :: d_alpha(:), d_coeff(:), d_s_factor(:), d_p_factor(:)
   integer :: iat, isp, is, ish, iprim, nprim, jprim

   allocate(shell_ang_mom(bas%nsh), d_s_factor(bas%nsh))
   nprim = 0
   do ish = 1, bas%nsh
      iat = bas%sh2at(ish)
      isp = mol%id(iat)
      is = bas%ish_at(iat)
      shell_ang_mom(ish) = bas%cgto(ish - is, isp)%ang
      nprim = nprim + bas%cgto(ish - is, isp)%nprim
      ! Normalization of the contracted basis function is included in the coefficients
      d_s_factor(ish) = 1.0_wp
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

   rc = trexio_write_basis_shell_factor(trex_file, d_s_factor)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO shell factor")
   if (allocated(error)) return

   rc = trexio_write_basis_prim_num(trex_file, nprim)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO primitive count")
   if (allocated(error)) return

   allocate(basis_shell_index(nprim), d_alpha(nprim), d_coeff(nprim), d_p_factor(nprim))
   jprim = 1
   do ish = 1, bas%nsh
      iat = bas%sh2at(ish)
      isp = mol%id(iat)
      is = bas%ish_at(iat)
      associate(p_cgto => bas%cgto(ish - is, isp))
         do iprim = 1, p_cgto%nprim
            basis_shell_index(jprim) = ish
            d_alpha(jprim) = real(p_cgto%alpha(iprim), c_double)
            d_coeff(jprim) = real(p_cgto%coeff(iprim), c_double)
            ! Coefficients already include normalization factors
            d_p_factor(jprim) = 1.0_c_double
            jprim = jprim + 1
         end do
      end associate
   end do

   rc = trexio_write_basis_shell_index(trex_file, basis_shell_index)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO primitive shell map")
   if (allocated(error)) return

   rc = trexio_write_basis_exponent(trex_file, d_alpha)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO primitive exponent")
   if (allocated(error)) return

   rc = trexio_write_basis_coefficient(trex_file, d_coeff)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO primitive coefficient")
   if (allocated(error)) return

   rc = trexio_write_basis_prim_factor(trex_file, d_p_factor)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO primitive factor")
end subroutine write_basis

subroutine write_ecp(trex_file, mol, wfn, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Converged wavefunction
   type(wavefunction_type), intent(in) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer(trexio_exit_code) :: rc
   integer, allocatable :: ecp_z_core(:)

   allocate(ecp_z_core(mol%nat))
   ecp_z_core = nint(mol%num(mol%id) - wfn%n0at)
   rc = trexio_write_ecp_z_core(trex_file, ecp_z_core)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO Z ECP core")
end subroutine write_ecp

subroutine write_ao(trex_file, bas, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Basis metadata used by the calculation
   type(basis_type), intent(in) :: bas
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer(trexio_exit_code) :: rc

   rc = trexio_write_ao_cartesian(trex_file, 0)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO AO cartesian flag")
   if (allocated(error)) return

   rc = trexio_write_ao_num(trex_file, bas%nao)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO AO count")
   if (allocated(error)) return

   rc = trexio_write_ao_shell(trex_file, bas%ao2sh)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO AO shell map")
end subroutine write_ao

subroutine write_mo(trex_file, mol, bas, wfn, error)
   !> Open TREXIO file handle
   integer(trexio_t), intent(in) :: trex_file
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Converged wavefunction
   type(wavefunction_type), intent(in) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(c_double), allocatable :: mo_coeff(:, :), mo_energy(:), mo_occ(:)
   integer, allocatable :: mo_spin(:)
   integer(trexio_exit_code) :: rc
   integer :: nmo, imo, imoflat, jsh, jat, jsp, js, jao_sphr, nsphr
   integer :: nspin_mo, spin, cspin, idx_trexio, idx_tblite, perm(15)
   logical :: spin_resolved_mo

   ! Represent only restricted closed-shell MOs as one spin-summed MO list
   ! and both restricted open-shell and unrestricted MOs with two spin channels
   spin_resolved_mo = wfn%nspin > 1 .or. abs(wfn%nuhf) > epsilon(1.0_wp)
   if (spin_resolved_mo) then
      nspin_mo = 2
      nmo = bas%nao * nspin_mo
      allocate(mo_coeff(bas%nao, nmo), mo_energy(nmo), mo_occ(nmo), mo_spin(nmo))
   else
      nspin_mo = 1
      nmo = bas%nao
      allocate(mo_coeff(bas%nao, nmo), mo_energy(nmo), mo_occ(nmo))
   end if

   mo_coeff(:, :) = 0.0_c_double
   mo_energy(:) = 0.0_c_double
   mo_occ(:) = 0.0_c_double
   if (allocated(mo_spin)) mo_spin(:) = 0

   ! Reorder MO coefficients to TREXIO ordering
   imoflat = 0
   do spin = 1, nspin_mo
      cspin = min(spin, wfn%nspin)
      do imo = 1, bas%nao
         imoflat = imoflat + 1
         mo_energy(imoflat) = real(wfn%emo(imo, cspin), c_double)
         if (spin_resolved_mo) then
            mo_occ(imoflat) = real(wfn%focc(imo, spin), c_double)
            mo_spin(imoflat) = spin - 1
         else
            mo_occ(imoflat) = real(wfn%focc(imo, 1) + wfn%focc(imo, 2), c_double)
         end if
         idx_trexio = 0
         do jsh = 1, bas%nsh
            nsphr = bas%nao_sh(jsh)
            jat = bas%sh2at(jsh)
            jsp = mol%id(jat)
            js = bas%ish_at(jat)
            call get_trexio_to_tblite_sphr_perm(bas%cgto(jsh - js, jsp)%ang, perm, error)
            if (allocated(error)) return
            do jao_sphr = 1, nsphr
               idx_trexio = idx_trexio + 1
               idx_tblite = bas%iao_sh(jsh) + perm(jao_sphr)
               mo_coeff(idx_trexio, imoflat) = real(wfn%coeff(idx_tblite, imo, cspin), &
                  & c_double)
            end do
         end do
      end do
   end do

   rc = trexio_write_mo_type(trex_file, "Tight-binding", 32)
   if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO MO type")
   if (allocated(error)) return

   rc = trexio_write_mo_num(trex_file, nmo)
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

   if (allocated(mo_spin)) then
      rc = trexio_write_mo_spin(trex_file, mo_spin)
      if (rc /= TREXIO_SUCCESS) call fatal_trexio(error, rc, "Failed to write TREXIO MO spins")
   end if
end subroutine write_mo

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
