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

!> @file tblite/scf/iterator.f90
!> Provides the implementation of the actual self-consistent field iteractions

!> Iterator for evaluating the Hamiltonian self-consistently
module tblite_scf_iterator
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_container, only : container_cache, container_list
   use tblite_disp, only : dispersion_type
   use tblite_integral_type, only : integral_type
   use tblite_mpi_utils, only : mpi_allreduce_sum, mpi_sync_error, new_mpi_work_partition
   use tblite_partition, only : work_partition, owns_index
   use tblite_scf_info, only : scf_info
   use tblite_scf_mixer, only : mixer_type
   use tblite_scf_potential, only : potential_type, mpi_add_pot_to_h1, reduce_potential
   use tblite_scf_solver, only : solver_type
   use tblite_wavefunction_fermi, only : get_fermi_filling
   use tblite_wavefunction_mulliken, only : get_mulliken_shell_charges, &
      & get_mulliken_atomic_multipoles
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_coulomb, only : tb_coulomb
   implicit none
   private

   public :: next_scf, get_mixer_dimension, get_electronic_energy, reduce
   public :: next_density, get_qat_from_qsh

contains

!> Evaluate self-consistent iteration for the density-dependent Hamiltonian
subroutine next_scf(iscf, mol, bas, wfn, solver, mixer, info, coulomb, dispersion, &
      & interactions, ints, pot, ccache, dcache, icache, &
      & energies, error, comm)
   !> Current iteration count
   integer, intent(inout) :: iscf
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Solver for the general eigenvalue problem
   class(solver_type), intent(inout) :: solver
   !> Convergence accelerator
   class(mixer_type), intent(inout) :: mixer
   !> Information on wavefunction data used to construct Hamiltonian
   type(scf_info), intent(in) :: info
   !> Container for coulombic interactions
   type(tb_coulomb), intent(in), optional :: coulomb
   !> Container for dispersion interactions
   class(dispersion_type), intent(in), optional :: dispersion
   !> Container for general interactions
   type(container_list), intent(in), optional :: interactions

   !> Integral container
   type(integral_type), intent(in) :: ints
   !> Density dependent potential shifts
   type(potential_type), intent(inout) :: pot
   !> Restart data for coulombic interactions
   type(container_cache), intent(inout), optional :: ccache
   !> Restart data for dispersion interactions
   type(container_cache), intent(inout), optional :: dcache
   !> Restart data for interaction containers
   type(container_cache), intent(inout), optional :: icache

   !> Self-consistent energy
   real(wp), intent(inout) :: energies(:)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Communicator to reduce the partitioned contributions over
   integer, intent(in), optional :: comm

   real(wp), allocatable :: eao(:)
   integer, allocatable :: columns(:)
   real(wp) :: ts
   type(work_partition) :: partition

   if (present(comm)) then
      call new_mpi_work_partition(error, partition, comm)
      if (allocated(error)) return
   end if

   if (iscf > 0) then
      call mixer%next(error)
      call mpi_sync_error(error, comm)
      if (allocated(error)) return

      call get_mixer(mixer, bas, wfn, info)
   end if

   iscf = iscf + 1
   call pot%reset
   if (present(coulomb) .and. present(ccache)) then
      call coulomb%get_potential(mol, ccache, wfn, pot)
   end if
   if (present(dispersion) .and. present(dcache)) then
      call dispersion%get_potential(mol, dcache, wfn, pot)
   end if
   if (present(interactions) .and. present(icache)) then
      call interactions%get_potential(mol, icache, wfn, pot)
   end if

   call reduce_potential(error, pot, comm)
   if (allocated(error)) return

   call mpi_add_pot_to_h1(error, bas, ints, pot, wfn%coeff, comm)
   if (allocated(error)) return

   call set_mixer(mixer, wfn, info)

   call next_density(wfn, solver, ints, ts, error)
   call mpi_sync_error(error, comm)
   if (allocated(error)) return

   if (ints%local) columns = ints%columns
   call get_mulliken_shell_charges(bas, ints%overlap, wfn%density, wfn%n0sh, wfn%qsh, partition, columns)
   call get_mulliken_atomic_multipoles(bas, ints%dipole, wfn%density, wfn%dpat, partition, columns)
   call get_mulliken_atomic_multipoles(bas, ints%quadrupole, wfn%density, wfn%qpat, partition, columns)
   if (partition%get_nparts() > 1) then
      call reduce_mulliken(wfn, error, comm)
      if (allocated(error)) return
   end if
   call get_qat_from_qsh(bas, wfn%qsh, wfn%qat)

   call diff_mixer(mixer, wfn, info)

   ! Assemble partitioned container and one-electron energies together.
   energies(:) = 0.0_wp
   if (present(coulomb) .and. present(ccache)) then
      call coulomb%get_energy(mol, ccache, wfn, energies)
   end if
   if (present(dispersion) .and. present(dcache)) then
      call dispersion%get_energy(mol, dcache, wfn, energies)
   end if
   if (present(interactions) .and. present(icache)) then
      call interactions%get_energy(mol, icache, wfn, energies)
   end if
   ! Preserve serial summation order; the entropy is replicated, not additive.
   if (partition%get_nparts() == 1) energies(:) = energies + ts / size(energies)
   allocate(eao(bas%nao), source=0.0_wp)
   call get_electronic_energy(ints%hamiltonian, wfn%density, eao, partition, columns)
   call reduce(energies, eao, bas%ao2at)
   call mpi_allreduce_sum(error, energies, comm)
   if (allocated(error)) return
   if (partition%get_nparts() > 1) energies(:) = energies + ts / size(energies)
end subroutine next_scf


!> Assemble all population components in one small collective, for either
!> column-distributed matrices or owner-computed contractions of full matrices.
subroutine reduce_mulliken(wfn, error, comm)
   type(wavefunction_type), intent(inout) :: wfn
   type(error_type), allocatable, intent(inout) :: error
   integer, intent(in) :: comm
   real(wp) :: buffer(size(wfn%qsh) + size(wfn%dpat) + size(wfn%qpat))
   integer :: nq, nd, nt

   nq = size(wfn%qsh)
   nd = size(wfn%dpat)
   nt = size(wfn%qpat)
   buffer(:nq) = reshape(wfn%qsh, [nq])
   buffer(nq+1:nq+nd) = reshape(wfn%dpat, [nd])
   buffer(nq+nd+1:) = reshape(wfn%qpat, [nt])
   call mpi_allreduce_sum(error, buffer, comm)
   if (allocated(error)) return
   wfn%qsh(:, :) = reshape(buffer(:nq), shape(wfn%qsh))
   wfn%dpat(:, :, :) = reshape(buffer(nq+1:nq+nd), shape(wfn%dpat))
   wfn%qpat(:, :, :) = reshape(buffer(nq+nd+1:), shape(wfn%qpat))
end subroutine reduce_mulliken


subroutine get_electronic_energy(h0, density, energies, partition, columns)
   real(wp), intent(in) :: h0(:, :)
   real(wp), intent(in) :: density(:, :, :)
   real(wp), intent(inout) :: energies(:)
   type(work_partition), intent(in), optional :: partition
   integer, intent(in), optional :: columns(2)

   integer :: iao, jao, spin, offset
   real(wp) :: eiao

   offset = 0
   if (present(columns)) offset = columns(1) - 1

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(h0, density, energies, partition, columns, offset) private(spin, iao, jao, eiao)
   do iao = 1, size(density, 2)
      if (.not.present(columns)) then
         if (.not.owns_index(partition, iao)) cycle
      end if
      eiao = 0.0_wp
      do spin = 1, size(density, 3)
         do jao = 1, size(density, 1)
            eiao = eiao + h0(jao, iao) * density(jao, iao, spin)
         end do
      end do
      energies(iao+offset) = energies(iao+offset) + eiao
   end do
end subroutine get_electronic_energy


subroutine reduce(reduced, full, map)
   real(wp), intent(inout) :: reduced(:)
   real(wp), intent(in) :: full(:)
   integer, intent(in) :: map(:)

   integer :: ix

   do ix = 1, size(map)
      reduced(map(ix)) = reduced(map(ix)) + full(ix)
   end do
end subroutine reduce


subroutine get_qat_from_qsh(bas, qsh, qat)
   type(basis_type), intent(in) :: bas
   real(wp), intent(in) :: qsh(:, :)
   real(wp), intent(out) :: qat(:, :)

   integer :: ish, ispin

   qat(:, :) = 0.0_wp
   !$omp parallel do schedule(runtime) collapse(2) default(none) &
   !$omp reduction(+:qat) shared(bas, qsh) private(ish)
   do ispin = 1, size(qsh, 2)
      do ish = 1, size(qsh, 1)
         qat(bas%sh2at(ish), ispin) = qat(bas%sh2at(ish), ispin) + qsh(ish, ispin)
      end do
   end do
end subroutine get_qat_from_qsh


function get_mixer_dimension(mol, bas, info) result(ndim)
   use tblite_scf_info, only : atom_resolved, shell_resolved
   type(structure_type), intent(in) :: mol
   type(basis_type), intent(in) :: bas
   type(scf_info), intent(in) :: info
   integer :: ndim

   ndim = 0

   select case(info%charge)
   case(atom_resolved)
      ndim = ndim + mol%nat
   case(shell_resolved)
      ndim = ndim + bas%nsh
   case default
      continue
   end select

   select case(info%dipole)
   case(atom_resolved)
      ndim = ndim + 3*mol%nat
   case default
      continue
   end select

   select case(info%quadrupole)
   case(atom_resolved)
      ndim = ndim + 6*mol%nat
   case default
      continue
   end select
end function get_mixer_dimension

subroutine set_mixer(mixer, wfn, info)
   use tblite_scf_info, only : atom_resolved, shell_resolved
   class(mixer_type), intent(inout) :: mixer
   type(wavefunction_type), intent(inout) :: wfn
   type(scf_info), intent(in) :: info

   select case(info%charge)
   case(atom_resolved)
      call mixer%set(wfn%qat)
   case(shell_resolved)
      call mixer%set(wfn%qsh)
   case default
      continue
   end select

   select case(info%dipole)
   case(atom_resolved)
      call mixer%set(wfn%dpat)
   case default
      continue
   end select

   select case(info%quadrupole)
   case(atom_resolved)
      ! Scale off-diagonal elements of packed quadrupole tensor in place, since
      ! they appear twice in the full tensor and enter the residual norm as such
      wfn%qpat([2, 4, 5], :, :) = wfn%qpat([2, 4, 5], :, :) * sqrt(2.0_wp)
      call mixer%set(wfn%qpat)
      wfn%qpat([2, 4, 5], :, :) = wfn%qpat([2, 4, 5], :, :) / sqrt(2.0_wp)
   case default
      continue
   end select
end subroutine set_mixer

subroutine diff_mixer(mixer, wfn, info)
   use tblite_scf_info, only : atom_resolved, shell_resolved
   class(mixer_type), intent(inout) :: mixer
   type(wavefunction_type), intent(inout) :: wfn
   type(scf_info), intent(in) :: info

   select case(info%charge)
   case(atom_resolved)
      call mixer%diff(wfn%qat)
   case(shell_resolved)
      call mixer%diff(wfn%qsh)
   case default
      continue
   end select

   select case(info%dipole)
   case(atom_resolved)
      call mixer%diff(wfn%dpat)
   case default
      continue
   end select

   select case(info%quadrupole)
   case(atom_resolved)
      ! Scale off-diagonal elements of packed quadrupole tensor in place, since
      ! they appear twice in the full tensor and enter the residual norm as such
      wfn%qpat([2, 4, 5], :, :) = wfn%qpat([2, 4, 5], :, :) * sqrt(2.0_wp)
      call mixer%diff(wfn%qpat)
      wfn%qpat([2, 4, 5], :, :) = wfn%qpat([2, 4, 5], :, :) / sqrt(2.0_wp)
   case default
      continue
   end select
end subroutine diff_mixer

subroutine get_mixer(mixer, bas, wfn, info)
   use tblite_scf_info, only : atom_resolved, shell_resolved
   class(mixer_type), intent(inout) :: mixer
   type(basis_type), intent(in) :: bas
   type(wavefunction_type), intent(inout) :: wfn
   type(scf_info), intent(in) :: info

   select case(info%charge)
   case(atom_resolved)
      call mixer%get(wfn%qat)
   case(shell_resolved)
      call mixer%get(wfn%qsh)
      call get_qat_from_qsh(bas, wfn%qsh, wfn%qat)
   case default
      continue
   end select

   select case(info%dipole)
   case(atom_resolved)
      call mixer%get(wfn%dpat)
   case default
      continue
   end select

   select case(info%quadrupole)
   case(atom_resolved)
      call mixer%get(wfn%qpat)
      ! Remove scaling of off-diagonal elements of packed quadrupole tensor
      wfn%qpat([2, 4, 5], :, :) = wfn%qpat([2, 4, 5], :, :) &
         & / sqrt(2.0_wp)
   case default
      continue
   end select
end subroutine get_mixer


subroutine next_density(wfn, solver, ints, ts, error)
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Solver for the general eigenvalue problem
   class(solver_type), intent(inout) :: solver
   !> Integral container
   type(integral_type), intent(in) :: ints
   !> Electronic entropy
   real(wp), intent(out) :: ts
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp) :: e_fermi, stmp(2)
   real(wp), allocatable :: focc(:)
   integer :: spin

   call solver%get_density(wfn%coeff, ints%overlap, wfn%emo, wfn%focc, wfn%density, error)
   do spin = 1, 2
      call get_electronic_entropy(wfn%focc(:, spin), wfn%kt, stmp(spin))
   end do
   ts = sum(stmp)
end subroutine next_density

subroutine get_electronic_entropy(occ, kt, s)
   real(wp), intent(in) :: occ(:)
   real(wp), intent(in) :: kt
   real(wp), intent(out) :: s

   s = sum(log(occ ** occ * (1 - occ) ** (1 - occ))) * kt
end subroutine get_electronic_entropy

end module tblite_scf_iterator
