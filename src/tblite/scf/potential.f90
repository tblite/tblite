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
   implicit none
   private

   public :: potential_type, new_potential, add_pot_to_h1


   !> Container for density dependent potential-shifts
   type :: potential_type
      !> Atom-resolved charge-dependent potential shift
      real(wp), allocatable :: vat(:)
      !> Shell-resolved charge-dependent potential shift
      real(wp), allocatable :: vsh(:)
      !> Orbital-resolved charge-dependent potential shift
      real(wp), allocatable :: vao(:)

      !> Atom-resolved dipolar potential
      real(wp), allocatable :: vdp(:, :)
      !> Atom-resolved quadrupolar potential
      real(wp), allocatable :: vqp(:, :)
   contains
      !> Reset the density dependent potential
      procedure :: reset
   end type potential_type


contains


!> Create a new potential object
subroutine new_potential(self, mol, bas)
   !> Instance of the density dependent potential
   type(potential_type), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Description of the basis set
   type(basis_type), intent(in) :: bas

   allocate(self%vat(mol%nat))
   allocate(self%vsh(bas%nsh))
   allocate(self%vao(bas%nao))

   allocate(self%vdp(3, mol%nat))
   allocate(self%vqp(6, mol%nat))
end subroutine new_potential

!> Reset the density dependent potential
subroutine reset(self)
   !> Instance of the density dependent potential
   class(potential_type), intent(inout) :: self

   self%vat(:) = 0.0_wp
   self%vsh(:) = 0.0_wp
   self%vao(:) = 0.0_wp
   self%vdp(:, :) = 0.0_wp
   self%vqp(:, :) = 0.0_wp
end subroutine reset

!> Add the collected potential shifts to the effective Hamiltonian
subroutine add_pot_to_h1(bas, h0, sint, dpint, qpint, pot, h1)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Effective one-electron Hamiltonian
   real(wp), intent(in) :: h0(:, :)
   !> Overlap integrals
   real(wp), intent(in) :: sint(:, :)
   !> Dipole moment integrals, moment operator is centered on last index
   real(wp), intent(in) :: dpint(:, :, :)
   !> Quadrupole moment integrals, moment operator is centered on last index
   real(wp), intent(in) :: qpint(:, :, :)
   !> Density dependent potential-shifts
   type(potential_type), intent(inout) :: pot
   !> Effective Hamiltonian
   real(wp), intent(out) :: h1(:, :)

   call add_vat_to_vsh(bas, pot%vat, pot%vsh)
   call add_vsh_to_vao(bas, pot%vsh, pot%vao)
   call add_vao_to_h1(bas, h0, sint, pot%vao, h1)
   call add_vmp_to_h1(bas, dpint, pot%vdp, h1)
   call add_vmp_to_h1(bas, qpint, pot%vqp, h1)
end subroutine add_pot_to_h1

!> Expand an atom-resolved potential shift to a shell-resolved potential shift
subroutine add_vat_to_vsh(bas, vat, vsh)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Atom-resolved charge-dependent potential shift
   real(wp), intent(in) :: vat(:)
   !> Shell-resolved charge-dependent potential shift
   real(wp), intent(inout) :: vsh(:)

   integer :: iat, ish, ii

   !$omp parallel do schedule(runtime) default(none) &
   !$omp reduction(+:vsh) shared(bas, vat) private(ii, ish, iat)
   do iat = 1, size(vat)
      ii = bas%ish_at(iat)
      do ish = 1, bas%nsh_at(iat)
         vsh(ii+ish) = vsh(ii+ish) + vat(iat)
      end do
   end do
end subroutine add_vat_to_vsh

!> Expand a shell-resolved potential shift to an orbital-resolved potential shift
subroutine add_vsh_to_vao(bas, vsh, vao)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Shell-resolved charge-dependent potential shift
   real(wp), intent(in) :: vsh(:)
   !> Orbital-resolved charge-dependent potential shift
   real(wp), intent(inout) :: vao(:)

   integer :: ish, iao, ii

   !$omp parallel do schedule(runtime) default(none) &
   !$omp reduction(+:vao) shared(bas, vsh) private(ii, iao, ish)
   do ish = 1, size(vsh)
      ii = bas%iao_sh(ish)
      do iao = 1, bas%nao_sh(ish)
         vao(ii+iao) = vao(ii+iao) + vsh(ish)
      end do
   end do
end subroutine add_vsh_to_vao


!> Add a charge-dependent potential to the Hamiltonian
subroutine add_vao_to_h1(bas, h0, sint, vao, h1)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Effective one-electron Hamiltonian
   real(wp), intent(in) :: h0(:, :)
   !> Overlap integrals
   real(wp), intent(in) :: sint(:, :)
   !> Orbital-resolved charge-dependent potential shift
   real(wp), intent(in) :: vao(:)
   !> Effective Hamiltonian
   real(wp), intent(out) :: h1(:, :)

   integer :: iao, jao

   !$omp parallel do collapse(2) schedule(runtime) default(none) &
   !$omp shared(h1, h0, bas, sint, vao) private(iao, jao)
   do iao = 1, bas%nao
      do jao = 1, bas%nao
         h1(jao, iao) = h0(jao, iao) - sint(jao, iao) * 0.5_wp * (vao(jao) + vao(iao))
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
   real(wp), intent(in) :: vmp(:, :)
   !> Effective Hamiltonian
   real(wp), intent(inout) :: h1(:, :)

   integer :: iao, jao

   !$omp parallel do collapse(2) schedule(runtime) default(none) &
   !$omp shared(h1, bas, mpint, vmp) private(iao, jao)
   do iao = 1, bas%nao
      do jao = 1, bas%nao
         h1(jao, iao) = h1(jao, iao) &
            & - 0.5_wp * dot_product(mpint(:, jao, iao), vmp(:, bas%ao2at(iao))) &
            & - 0.5_wp * dot_product(mpint(:, iao, jao), vmp(:, bas%ao2at(jao)))
      end do
   end do
end subroutine add_vmp_to_h1

end module tblite_scf_potential
