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

!> Implementation of the density dependent potential and its contribution
!> to the effective Hamiltonian
module tblite_scf_potential
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_integral_type, only : integral_type
   use tblite_wavefunction_spin, only : magnet_to_updown
   implicit none
   private

   public :: potential_type, new_potential, add_pot_to_h1


   !> Container for density dependent potential-shifts
   type :: potential_type
      !> Atom-resolved charge-dependent potential shift
      real(wp), allocatable :: vat(:, :)
      !> Shell-resolved charge-dependent potential shift
      real(wp), allocatable :: vsh(:, :)
      !> Orbital-resolved charge-dependent potential shift
      real(wp), allocatable :: vao(:, :)

      !> Atom-resolved dipolar potential
      real(wp), allocatable :: vdp(:, :, :)
      !> Atom-resolved quadrupolar potential
      real(wp), allocatable :: vqp(:, :, :)
   contains
      !> Reset the density dependent potential
      procedure :: reset
   end type potential_type


contains


!> Create a new potential object
subroutine new_potential(self, mol, bas, nspin)
   !> Instance of the density dependent potential
   type(potential_type), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Description of the basis set
   type(basis_type), intent(in) :: bas
   !> Number of spin channels
   integer, intent(in) :: nspin

   allocate(self%vat(mol%nat, nspin))
   allocate(self%vsh(bas%nsh, nspin))
   allocate(self%vao(bas%nao, nspin))

   allocate(self%vdp(3, mol%nat, nspin))
   allocate(self%vqp(6, mol%nat, nspin))
end subroutine new_potential

!> Reset the density dependent potential
subroutine reset(self)
   !> Instance of the density dependent potential
   class(potential_type), intent(inout) :: self

   self%vat(:, :) = 0.0_wp
   self%vsh(:, :) = 0.0_wp
   self%vao(:, :) = 0.0_wp
   self%vdp(:, :, :) = 0.0_wp
   self%vqp(:, :, :) = 0.0_wp
end subroutine reset

!> Add the collected potential shifts to the effective Hamiltonian
subroutine add_pot_to_h1(bas, ints, pot, h1)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Integral container
   type(integral_type), intent(in) :: ints
   !> Density dependent potential-shifts
   type(potential_type), intent(inout) :: pot
   !> Effective Hamiltonian
   real(wp), intent(out) :: h1(:, :, :)

   h1(:, :, 1) = ints%hamiltonian
   if (size(h1, 3) > 1) h1(:, :, 2:) = 0.0_wp

   call add_vat_to_vsh(bas, pot%vat, pot%vsh)
   call add_vsh_to_vao(bas, pot%vsh, pot%vao)
   call add_vao_to_h1(bas, ints%overlap, pot%vao, h1)
   call add_vmp_to_h1(bas, ints%dipole, pot%vdp, h1)
   call add_vmp_to_h1(bas, ints%quadrupole, pot%vqp, h1)

   call magnet_to_updown(h1)
end subroutine add_pot_to_h1

!> Expand an atom-resolved potential shift to a shell-resolved potential shift
subroutine add_vat_to_vsh(bas, vat, vsh)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Atom-resolved charge-dependent potential shift
   real(wp), intent(in) :: vat(:, :)
   !> Shell-resolved charge-dependent potential shift
   real(wp), intent(inout) :: vsh(:, :)

   integer :: iat, ish, ii, spin

   !$omp parallel do schedule(runtime) collapse(2) default(none) &
   !$omp reduction(+:vsh) shared(bas, vat) private(spin, ii, ish, iat)
   do spin = 1, size(vat, 2)
      do iat = 1, size(vat, 1)
         ii = bas%ish_at(iat)
         do ish = 1, bas%nsh_at(iat)
            vsh(ii+ish, spin) = vsh(ii+ish, spin) + vat(iat, spin)
         end do
      end do
   end do
end subroutine add_vat_to_vsh

!> Expand a shell-resolved potential shift to an orbital-resolved potential shift
subroutine add_vsh_to_vao(bas, vsh, vao)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Shell-resolved charge-dependent potential shift
   real(wp), intent(in) :: vsh(:, :)
   !> Orbital-resolved charge-dependent potential shift
   real(wp), intent(inout) :: vao(:, :)

   integer :: ish, iao, ii, spin

   !$omp parallel do schedule(runtime) collapse(2) default(none) &
   !$omp reduction(+:vao) shared(bas, vsh) private(ii, iao, ish)
   do spin = 1, size(vsh, 2)
      do ish = 1, size(vsh, 1)
         ii = bas%iao_sh(ish)
         do iao = 1, bas%nao_sh(ish)
            vao(ii+iao, spin) = vao(ii+iao, spin) + vsh(ish, spin)
         end do
      end do
   end do
end subroutine add_vsh_to_vao


!> Add a charge-dependent potential to the Hamiltonian
subroutine add_vao_to_h1(bas, sint, vao, h1)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Overlap integrals
   real(wp), intent(in) :: sint(:, :)
   !> Orbital-resolved charge-dependent potential shift
   real(wp), intent(in) :: vao(:, :)
   !> Effective Hamiltonian
   real(wp), intent(inout) :: h1(:, :, :)

   integer :: iao, jao, spin

   !$omp parallel do collapse(3) schedule(runtime) default(none) &
   !$omp shared(h1, bas, sint, vao) private(spin, iao, jao)
   do spin = 1, size(h1, 3)
      do iao = 1, bas%nao
         do jao = 1, bas%nao
            h1(jao, iao, spin) = h1(jao, iao, spin) &
               & - sint(jao, iao) * 0.5_wp * (vao(jao, spin) + vao(iao, spin))
         end do
      end do
   end do
end subroutine add_vao_to_h1

!> Add a multipolar potential to the Hamiltonian
subroutine add_vmp_to_h1(bas, mpint, vmp, h1)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Multipole integrals, multipole operator is always centered on last index
   real(wp), intent(in) :: mpint(:, :, :)
   !> Multipole potential
   real(wp), intent(in) :: vmp(:, :, :)
   !> Effective Hamiltonian
   real(wp), intent(inout) :: h1(:, :, :)

   integer :: iao, jao, iat, jat, nmp, spin

   nmp = min(size(mpint, 1), size(vmp, 1))

   !$omp parallel do collapse(3) schedule(runtime) default(none) &
   !$omp shared(h1, bas, mpint, vmp, nmp) private(spin, iao, jao, iat, jat)
   do spin = 1, size(h1, 3)
      do iao = 1, bas%nao
         do jao = 1, bas%nao
            iat = bas%ao2at(iao)
            jat = bas%ao2at(jao)
            h1(jao, iao, spin) = h1(jao, iao, spin) &
               & - 0.5_wp * dot_product(mpint(:nmp, jao, iao), vmp(:nmp, iat, spin)) &
               & - 0.5_wp * dot_product(mpint(:nmp, iao, jao), vmp(:nmp, jat, spin))
         end do
      end do
   end do
end subroutine add_vmp_to_h1

end module tblite_scf_potential
