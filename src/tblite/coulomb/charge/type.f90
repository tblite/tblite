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

!> Isotropic second-order electrostatics using the DFTB Coulomb functional.
module tblite_coulomb_charge_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_blas, only : dot, gemv, symv
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_type, only : coulomb_type
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: coulomb_charge_type


   !> General second-order electrostatics
   type, extends(coulomb_type), abstract :: coulomb_charge_type
      !> Number of shells for each atom
      integer, allocatable :: nshell(:)
      !> Index offset for each shell
      integer, allocatable :: offset(:)
   contains
      !> Update container cache
      procedure :: update
      !> Get information about density dependent quantities used in the energy
      procedure :: variable_info
      !> Evaluate selfconsistent energy of the interaction
      procedure :: get_energy
      !> Evaluate charge dependent potential shift from the interaction
      procedure :: get_potential
      !> Evaluate gradient contributions from the selfconsistent interaction
      procedure :: get_gradient
      !> Evaluate Coulomb matrix
      procedure(get_coulomb_matrix), deferred :: get_coulomb_matrix
      !> Evaluate uncontracted derivatives of Coulomb matrix
      procedure(get_coulomb_derivs), deferred :: get_coulomb_derivs
   end type coulomb_charge_type

   abstract interface
      !> Evaluate coulomb matrix
      subroutine get_coulomb_matrix(self, mol, cache, amat)
         import :: coulomb_charge_type, structure_type, coulomb_cache, wp
         !> Instance of the electrostatic container
         class(coulomb_charge_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Reusable data container
         type(coulomb_cache), intent(inout) :: cache
         !> Coulomb matrix
         real(wp), contiguous, intent(out) :: amat(:, :)
      end subroutine get_coulomb_matrix

      !> Evaluate uncontracted derivatives of Coulomb matrix
      subroutine get_coulomb_derivs(self, mol, cache, qat, qsh, dadr, dadL, atrace)
         import :: coulomb_charge_type, structure_type, coulomb_cache, wp
         !> Instance of the electrostatic container
         class(coulomb_charge_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Reusable data container
         type(coulomb_cache), intent(inout) :: cache
         !> Atomic partial charges
         real(wp), intent(in) :: qat(:)
         !> Shell-resolved partial charges
         real(wp), intent(in) :: qsh(:)
         !> Derivative of interactions with respect to cartesian displacements
         real(wp), contiguous, intent(out) :: dadr(:, :, :)
         !> Derivative of interactions with respect to strain deformations
         real(wp), contiguous, intent(out) :: dadL(:, :, :)
         !> On-site derivatives with respect to cartesian displacements
         real(wp), contiguous, intent(out) :: atrace(:, :)
      end subroutine get_coulomb_derivs
   end interface

contains


!> Update container cache
subroutine update(self, mol, cache)
   !> Instance of the electrostatic container
   class(coulomb_charge_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache

   call cache%update(mol)

   if (.not.allocated(cache%amat)) then
      allocate(cache%amat(sum(self%nshell), sum(self%nshell)))
   end if
   call self%get_coulomb_matrix(mol, cache, cache%amat)

   if (.not.allocated(cache%vvec)) then
      allocate(cache%vvec(sum(self%nshell)))
   end if

end subroutine update


!> Evaluate selfconsistent energy of the interaction
subroutine get_energy(self, mol, cache, wfn, energy)
   !> Instance of the electrostatic container
   class(coulomb_charge_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Electrostatic energy
   real(wp), intent(inout) :: energy

   call symv(cache%amat, wfn%qsh(:, 1), cache%vvec, alpha=0.5_wp)
   energy = energy + dot(cache%vvec, wfn%qsh(:, 1))

end subroutine get_energy


!> Evaluate charge dependent potential shift from the interaction
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the electrostatic container
   class(coulomb_charge_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   call symv(cache%amat, wfn%qsh(:, 1), pot%vsh(:, 1), beta=1.0_wp)

end subroutine get_potential


!> Evaluate gradient contributions from the selfconsistent interaction
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the electrostatic container
   class(coulomb_charge_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular gradient of the repulsion energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the repulsion energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   integer :: ndim
   real(wp), allocatable :: dadr(:, :, :), dadL(:, :, :), datr(:, :)

   ndim = sum(self%nshell)
   allocate(dadr(3, mol%nat, ndim), dadL(3, 3, ndim), datr(3, ndim))

   call self%get_coulomb_derivs(mol, cache, wfn%qat(:, 1), wfn%qsh(:, 1), dadr, dadL, datr)

   call gemv(dadr, wfn%qsh(:, 1), gradient, beta=1.0_wp)
   call gemv(dadL, wfn%qsh(:, 1), sigma, beta=1.0_wp, alpha=0.5_wp)

end subroutine get_gradient


!> Get information about density dependent quantities used in the energy
pure function variable_info(self) result(info)
   use tblite_scf_info, only : scf_info, shell_resolved
   !> Instance of the electrostatic container
   class(coulomb_charge_type), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=shell_resolved)
end function variable_info


end module tblite_coulomb_charge_type
