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
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_container, only : container_cache, container_list
   use tblite_disp, only : dispersion_type
   use tblite_integral_type, only : integral_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_wavefunction_fermi, only : get_fermi_filling
   use tblite_wavefunction_mulliken, only : get_mulliken_shell_charges, &
      & get_mulliken_atomic_multipoles
   use tblite_xtb_coulomb, only : tb_coulomb
   use tblite_scf_mixer, only : mixer_type
   use tblite_scf_mixer_input, only : trial_version
   use tblite_scf_info, only : scf_info
   use tblite_scf_potential, only : potential_type, add_pot_to_h1
   use tblite_scf_solver, only : solver_type
   implicit none
   private

   public :: next_outer_scf, next_inner_scf, get_mixer_dimension, get_electronic_energy, reduce
   public :: next_density, get_qat_from_qsh, get_mixer_vector, set_mixer_vector

contains

!> Evaluate an outer self-consistent iteration for the density-dependent Hamiltonian
subroutine next_outer_scf(iscf, mol, bas, wfn, solver, mixer, info, coulomb, dispersion, &
      & interactions, ints, pot, ccache, dcache, icache, energies, error)
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

   if (mixer%use_inner_scf() .and. iscf >= mixer%trial_start) then
      call next_inner_scf(iscf, mol, bas, wfn, solver, mixer, info, coulomb, dispersion, &
         & interactions, ints, pot, ccache, dcache, icache, energies, error)
   else
      call next_standard_scf(iscf, mol, bas, wfn, solver, mixer, info, coulomb, dispersion, &
         & interactions, ints, pot, ccache, dcache, icache, energies, error)
   end if
end subroutine next_outer_scf


!> Evaluate one standard self-consistent iteration
subroutine next_standard_scf(iscf, mol, bas, wfn, solver, mixer, info, coulomb, dispersion, &
      & interactions, ints, pot, ccache, dcache, icache, energies, error)
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

   if (iscf > 0) then
      call mixer%next(error)
      if (allocated(error)) return

      call get_mixer(mixer, bas, wfn, info)
   end if

   iscf = iscf + 1
   call set_mixer(mixer, wfn, info)
   mixer%trial_evaluations = mixer%trial_evaluations + 1
   call evaluate_scf_map(mol, bas, wfn, solver, info, coulomb, dispersion, &
      & interactions, ints, pot, ccache, dcache, icache, energies, error)
   if (allocated(error)) return
   call diff_mixer(mixer, wfn, info)
end subroutine next_standard_scf


!> Evaluate a self-consistent iteration with trial damping of the mixer step
subroutine next_inner_scf(iscf, mol, bas, wfn, solver, mixer, info, coulomb, dispersion, &
      & interactions, ints, pot, ccache, dcache, icache, energies, error)
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

   integer :: ilambda, best
   real(wp) :: lambda, residual, score, best_score, energy, best_energy
   real(wp) :: start_residual, previous_energy
   real(wp), allocatable :: qstart(:), qmix(:), qtrial(:), qout(:), qbest(:), &
      & trial_energies(:), best_energies(:)
   type(wavefunction_type) :: wfn_start, wfn_trial, wfn_best, wfn_input
   logical :: early_accept
   real(wp), parameter :: oda_lambdas(*) = [1.0_wp, 0.5_wp, 0.25_wp, 0.1_wp]
   real(wp), parameter :: mesa_lambdas(*) = [1.0_wp, 0.5_wp, 0.25_wp, 0.1_wp]
   real(wp), parameter :: early_accept_factor = 0.8_wp
   real(wp), parameter :: early_accept_residual = 1.0e-2_wp

   if (iscf == 0) then
      call next_standard_scf(iscf, mol, bas, wfn, solver, mixer, info, coulomb, dispersion, &
         & interactions, ints, pot, ccache, dcache, icache, energies, error)
      return
   end if

   wfn_start = wfn
   allocate(qstart(wfn%nspin*get_mixer_dimension(mol, bas, info)))
   allocate(qmix(size(qstart)), qtrial(size(qstart)), qout(size(qstart)), qbest(size(qstart)))
   allocate(trial_energies(size(energies)), best_energies(size(energies)))
   call get_mixer_vector(qstart, bas, wfn_start, info)
   start_residual = max(mixer%get_error(), epsilon(1.0_wp))
   previous_energy = sum(energies)

   call mixer%next(error)
   if (allocated(error)) return
   call get_mixer(mixer, bas, wfn, info)
   call get_mixer_vector(qmix, bas, wfn, info)

   best = 0
   best_score = huge(best_score)
   best_energy = huge(best_energy)
   early_accept = .false.

   select case(mixer%trial)
   case(trial_version%oda)
      do ilambda = 1, size(oda_lambdas)
         lambda = oda_lambdas(ilambda)
         call evaluate_trial(lambda)
         if (allocated(error)) exit
         if (early_accept) exit
      end do
   case(trial_version%mesa)
      do ilambda = 1, size(mesa_lambdas)
         lambda = mesa_lambdas(ilambda)
         call evaluate_trial(lambda)
         if (allocated(error)) exit
         if (early_accept) exit
      end do
   case(trial_version%default)
      call evaluate_trial(1.0_wp)
   case default
      call fatal_error(error, "Unknown trial version selected for inner SCF")
   end select
   if (allocated(error)) return

   if (best == 0) return

   iscf = iscf + 1
   wfn = wfn_best
   energies(:) = best_energies

   wfn_input = wfn_best
   call set_mixer_vector(qbest, bas, wfn_input, info)
   call set_mixer(mixer, wfn_input, info)
   call diff_mixer(mixer, wfn, info)

contains

subroutine evaluate_trial(lambda)
   real(wp), intent(in) :: lambda

   wfn_trial = wfn_start
   qtrial(:) = qstart + lambda*(qmix - qstart)
   call set_mixer_vector(qtrial, bas, wfn_trial, info)
   trial_energies(:) = 0.0_wp
   mixer%trial_evaluations = mixer%trial_evaluations + 1
   call evaluate_scf_map(mol, bas, wfn_trial, solver, info, coulomb, dispersion, &
      & interactions, ints, pot, ccache, dcache, icache, trial_energies, error)
   if (allocated(error)) return

   call get_mixer_vector(qout, bas, wfn_trial, info)
   residual = sqrt(sum((qout - qtrial)**2)/real(size(qtrial), wp))
   energy = sum(trial_energies)
   score = residual
   if (mixer%trial == trial_version%mesa) then
      score = residual + max(0.0_wp, energy - sum(energies)) / &
         & max(1.0_wp, abs(sum(energies)))
   end if

   if (score < best_score .or. &
      & (abs(score - best_score) <= epsilon(1.0_wp) .and. energy < best_energy)) then
      best = ilambda
      best_score = score
      best_energy = energy
      qbest(:) = qtrial
      best_energies(:) = trial_energies
      wfn_best = wfn_trial
      early_accept = residual <= early_accept_residual .and. &
         & score <= early_accept_factor*start_residual
      if (mixer%trial == trial_version%mesa) then
         early_accept = early_accept .and. &
            & energy <= previous_energy + max(1.0e-8_wp, 1.0e-10_wp*abs(previous_energy))
      end if
   end if
end subroutine evaluate_trial

end subroutine next_inner_scf


!> Evaluate the SCF map for the density-dependent Hamiltonian
subroutine evaluate_scf_map(mol, bas, wfn, solver, info, coulomb, dispersion, &
      & interactions, ints, pot, ccache, dcache, icache, energies, error)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Solver for the general eigenvalue problem
   class(solver_type), intent(inout) :: solver
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

   real(wp), allocatable :: eao(:)
   real(wp) :: ts

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

   call add_pot_to_h1(bas, ints, pot, wfn%coeff)

   call next_density(wfn, solver, ints, ts, error)
   if (allocated(error)) return

   call get_mulliken_shell_charges(bas, ints%overlap, wfn%density, wfn%n0sh, &
      & wfn%qsh)
   call get_qat_from_qsh(bas, wfn%qsh, wfn%qat)

   call get_mulliken_atomic_multipoles(bas, ints%dipole, wfn%density, &
      & wfn%dpat)
   call get_mulliken_atomic_multipoles(bas, ints%quadrupole, wfn%density, &
      & wfn%qpat)

   allocate(eao(bas%nao), source=0.0_wp)
   call get_electronic_energy(ints%hamiltonian, wfn%density, eao)

   energies(:) = ts / size(energies)
   call reduce(energies, eao, bas%ao2at)

   if (present(coulomb) .and. present(ccache)) then
      call coulomb%get_energy(mol, ccache, wfn, energies)
   end if
   if (present(dispersion) .and. present(dcache)) then
      call dispersion%get_energy(mol, dcache, wfn, energies)
   end if
   if (present(interactions) .and. present(icache)) then
      call interactions%get_energy(mol, icache, wfn, energies)
   end if
end subroutine evaluate_scf_map


subroutine get_electronic_energy(h0, density, energies)
   real(wp), intent(in) :: h0(:, :)
   real(wp), intent(in) :: density(:, :, :)
   real(wp), intent(inout) :: energies(:)

   integer :: iao, jao, spin
   real(wp) :: eiao

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(h0, density, energies) private(spin, iao, jao, eiao)
   do iao = 1, size(density, 2)
      eiao = 0.0_wp
      do spin = 1, size(density, 3)
         do jao = 1, size(density, 1)
            eiao = eiao + h0(jao, iao) * density(jao, iao, spin)
         end do
      end do
      energies(iao) = energies(iao) + eiao
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
   end select

   select case(info%dipole)
   case(atom_resolved)
      ndim = ndim + 3*mol%nat
   end select

   select case(info%quadrupole)
   case(atom_resolved)
      ndim = ndim + 6*mol%nat
   end select
end function get_mixer_dimension

subroutine set_mixer(mixer, wfn, info)
   use tblite_scf_info, only : atom_resolved, shell_resolved
   class(mixer_type), intent(inout) :: mixer
   type(wavefunction_type), intent(in) :: wfn
   type(scf_info), intent(in) :: info

   select case(info%charge)
   case(atom_resolved)
      call mixer%set(wfn%qat)
   case(shell_resolved)
      call mixer%set(wfn%qsh)
   end select

   select case(info%dipole)
   case(atom_resolved)
      call mixer%set(wfn%dpat)
   end select

   select case(info%quadrupole)
   case(atom_resolved)
      call mixer%set(wfn%qpat)
   end select
end subroutine set_mixer

subroutine diff_mixer(mixer, wfn, info)
   use tblite_scf_info, only : atom_resolved, shell_resolved
   class(mixer_type), intent(inout) :: mixer
   type(wavefunction_type), intent(in) :: wfn
   type(scf_info), intent(in) :: info

   select case(info%charge)
   case(atom_resolved)
      call mixer%diff(wfn%qat)
   case(shell_resolved)
      call mixer%diff(wfn%qsh)
   end select

   select case(info%dipole)
   case(atom_resolved)
      call mixer%diff(wfn%dpat)
   end select

   select case(info%quadrupole)
   case(atom_resolved)
      call mixer%diff(wfn%qpat)
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
   end select

   select case(info%dipole)
   case(atom_resolved)
      call mixer%get(wfn%dpat)
   end select

   select case(info%quadrupole)
   case(atom_resolved)
      call mixer%get(wfn%qpat)
   end select
end subroutine get_mixer


subroutine get_mixer_vector(qvec, bas, wfn, info)
   use tblite_scf_info, only : atom_resolved, shell_resolved
   real(wp), intent(out) :: qvec(:)
   type(basis_type), intent(in) :: bas
   type(wavefunction_type), intent(in) :: wfn
   type(scf_info), intent(in) :: info

   integer :: offset

   offset = 0
   select case(info%charge)
   case(atom_resolved)
      call pack_2d(qvec, offset, wfn%qat)
   case(shell_resolved)
      call pack_2d(qvec, offset, wfn%qsh)
   end select

   select case(info%dipole)
   case(atom_resolved)
      call pack_3d(qvec, offset, wfn%dpat)
   end select

   select case(info%quadrupole)
   case(atom_resolved)
      call pack_3d(qvec, offset, wfn%qpat)
   end select

contains

subroutine pack_2d(qvec, offset, array)
   real(wp), intent(inout) :: qvec(:)
   integer, intent(inout) :: offset
   real(wp), intent(in) :: array(:, :)

   qvec(offset+1:offset+size(array)) = reshape(array, [size(array)])
   offset = offset + size(array)
end subroutine pack_2d

subroutine pack_3d(qvec, offset, array)
   real(wp), intent(inout) :: qvec(:)
   integer, intent(inout) :: offset
   real(wp), intent(in) :: array(:, :, :)

   qvec(offset+1:offset+size(array)) = reshape(array, [size(array)])
   offset = offset + size(array)
end subroutine pack_3d

end subroutine get_mixer_vector


subroutine set_mixer_vector(qvec, bas, wfn, info)
   use tblite_scf_info, only : atom_resolved, shell_resolved
   real(wp), intent(in) :: qvec(:)
   type(basis_type), intent(in) :: bas
   type(wavefunction_type), intent(inout) :: wfn
   type(scf_info), intent(in) :: info

   integer :: offset

   offset = 0
   select case(info%charge)
   case(atom_resolved)
      call unpack_2d(qvec, offset, wfn%qat)
   case(shell_resolved)
      call unpack_2d(qvec, offset, wfn%qsh)
      call get_qat_from_qsh(bas, wfn%qsh, wfn%qat)
   end select

   select case(info%dipole)
   case(atom_resolved)
      call unpack_3d(qvec, offset, wfn%dpat)
   end select

   select case(info%quadrupole)
   case(atom_resolved)
      call unpack_3d(qvec, offset, wfn%qpat)
   end select

contains

subroutine unpack_2d(qvec, offset, array)
   real(wp), intent(in) :: qvec(:)
   integer, intent(inout) :: offset
   real(wp), intent(out) :: array(:, :)

   array(:, :) = reshape(qvec(offset+1:offset+size(array)), shape(array))
   offset = offset + size(array)
end subroutine unpack_2d

subroutine unpack_3d(qvec, offset, array)
   real(wp), intent(in) :: qvec(:)
   integer, intent(inout) :: offset
   real(wp), intent(out) :: array(:, :, :)

   array(:, :, :) = reshape(qvec(offset+1:offset+size(array)), shape(array))
   offset = offset + size(array)
end subroutine unpack_3d

end subroutine set_mixer_vector


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
