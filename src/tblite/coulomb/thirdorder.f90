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

!> Isotropic third-order onsite correction
module tblite_coulomb_thirdorder
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_container_cache, only : container_cache
   use tblite_coulomb_charge, only : effective_coulomb
   use tblite_coulomb_multipole, only : damped_multipole
   use tblite_coulomb_type, only : coulomb_type
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: onsite_thirdorder, new_onsite_thirdorder

   !> Onsite correction for third-order charge expansion
   type, extends(coulomb_type) :: onsite_thirdorder
      !> Whether the third order contribution is shell-dependent
      logical :: shell_resolved
      !> Number of shell for each atom
      integer, allocatable :: nsh_at(:)
      !> Shell offset for each atom
      integer, allocatable :: ish_at(:)
      !> Hubbard derivatives for each species
      real(wp), allocatable :: hubbard_derivs(:, :)
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
   end type onsite_thirdorder

   character(len=*), parameter :: label = "onsite third-order electrostatics"

contains


!> Create new onsite third-order contribution
subroutine new_onsite_thirdorder(self, mol, hubbard_derivs, nshell)
   !> Instance of the electrostatic container
   type(onsite_thirdorder), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Hubbard derivatives
   real(wp), intent(in) :: hubbard_derivs(:, :)
   !> Number of shells for each species
   integer, intent(in), optional :: nshell(:)

   integer :: ind, iat

   self%label = label
   self%hubbard_derivs = hubbard_derivs

   self%shell_resolved = present(nshell)
   if (present(nshell)) then
      self%nsh_at = nshell(mol%id)

      allocate(self%ish_at(mol%nat))
      ind = 0
      do iat = 1, mol%nat
         self%ish_at(iat) = ind
         ind = ind + self%nsh_at(iat)
      end do
   end if

end subroutine new_onsite_thirdorder


!> Update container cache
subroutine update(self, mol, cache)
   !> Instance of the electrostatic container
   class(onsite_thirdorder), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

end subroutine update


!> Evaluate selfconsistent energy of the interaction
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the electrostatic container
   class(onsite_thirdorder), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Electrostatic energy
   real(wp), intent(inout) :: energies(:)

   integer :: iat, izp, ii, ish

   if (self%shell_resolved) then
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ii = self%ish_at(iat)
         do ish = 1, self%nsh_at(iat)
            energies(iat) = energies(iat) &
               & + wfn%qsh(ii+ish, 1)**3 * self%hubbard_derivs(ish, izp) / 3.0_wp
         end do
      end do
   else
      do iat = 1, mol%nat
         izp = mol%id(iat)
         energies(iat) = energies(iat) &
            & + wfn%qat(iat, 1)**3 * self%hubbard_derivs(1, izp) / 3.0_wp
      end do
   end if
end subroutine get_energy


!> Evaluate charge dependent potential shift from the interaction
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the electrostatic container
   class(onsite_thirdorder), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   integer :: iat, izp, ii, ish

   if (self%shell_resolved) then
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ii = self%ish_at(iat)
         do ish = 1, self%nsh_at(iat)
            pot%vsh(ii+ish, 1) = pot%vsh(ii+ish, 1) &
               & + wfn%qsh(ii+ish, 1)**2 * self%hubbard_derivs(ish, izp)
         end do
      end do
   else
      do iat = 1, mol%nat
         izp = mol%id(iat)
         pot%vat(iat, 1) = pot%vat(iat, 1) + wfn%qat(iat, 1)**2 * self%hubbard_derivs(1, izp)
      end do
   end if
end subroutine get_potential


!> Evaluate gradient contributions from the selfconsistent interaction
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the electrostatic container
   class(onsite_thirdorder), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular gradient of the repulsion energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the repulsion energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

end subroutine get_gradient


!> Get information about density dependent quantities used in the energy
pure function variable_info(self) result(info)
   use tblite_scf_info, only : scf_info, atom_resolved, shell_resolved
   !> Instance of the electrostatic container
   class(onsite_thirdorder), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=merge(shell_resolved, atom_resolved, self%shell_resolved))
end function variable_info

end module tblite_coulomb_thirdorder
