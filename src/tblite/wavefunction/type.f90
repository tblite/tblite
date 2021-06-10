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

!> Base model for the extended tight binding hamiltonian
module tblite_wavefunction_type
   use mctc_env, only : wp
   use tblite_blas, only : gemm
   implicit none
   private

   public :: wavefunction_type, new_wavefunction, get_density_matrix

   type :: wavefunction_type
      !> Electronic temperature
      real(wp) :: kt = 0.0_wp
      !> Number of electrons in this wavefunction
      real(wp) :: nocc = 0.0_wp
      !> Number of unpaired electrons in this wavefunction
      real(wp) :: nuhf = 0.0_wp
      !> Index of the highest occupied molecular orbitals (alpha space)
      integer :: homoa = 0
      !> Index of the highest occupied molecular orbitals (beta space)
      integer :: homob = 0
      !> Reference occupation number for each atom, shape: [nat]
      real(wp), allocatable :: n0at(:)
      !> Reference occupation number for each shell, shape: [nsh]
      real(wp), allocatable :: n0sh(:)

      !> Density matrix, shape: [nao, nao]
      real(wp), allocatable :: density(:, :)
      !> Orbital coefficients, shape: [nao, nao]
      real(wp), allocatable :: coeff(:, :)
      !> Orbital energies, eigenvalues, shape: [nao]
      real(wp), allocatable :: emo(:)
      !> Occupation numbers, shape: [nao]
      real(wp), allocatable :: focc(:)

      !> Number of electrons for each atom, shape: [nat]
      real(wp), allocatable :: qat(:)
      !> Number of electrons for each shell, shape: [nsh]
      real(wp), allocatable :: qsh(:)

      !> Atomic dipole moments for each atom, shape: [3, nat]
      real(wp), allocatable :: dpat(:, :)
      !> Atomic quadrupole moments for each atom, shape: [5, nat]
      real(wp), allocatable :: qpat(:, :)
   end type wavefunction_type

contains


subroutine new_wavefunction(self, nat, nsh, nao, kt)
   type(wavefunction_type), intent(out) :: self
   integer, intent(in) :: nat
   integer, intent(in) :: nsh
   integer, intent(in) :: nao
   real(wp), intent(in) :: kt

   self%kt = kt

   allocate(self%n0at(nat))
   allocate(self%n0sh(nsh))

   allocate(self%density(nao, nao))
   allocate(self%coeff(nao, nao))
   allocate(self%emo(nao))
   allocate(self%focc(nao))

   allocate(self%qat(nat))
   allocate(self%qsh(nsh))

   allocate(self%dpat(3, nat))
   allocate(self%qpat(6, nat))

   self%qat(:) = 0.0_wp
   self%qsh(:) = 0.0_wp
   self%dpat(:, :) = 0.0_wp
   self%qpat(:, :) = 0.0_wp
end subroutine new_wavefunction


subroutine get_density_matrix(focc, coeff, pmat)
   real(wp), intent(in) :: focc(:)
   real(wp), contiguous, intent(in) :: coeff(:, :)
   real(wp), contiguous, intent(out) :: pmat(:, :)

   real(wp), allocatable :: scratch(:, :)
   integer :: iao, jao

   allocate(scratch(size(pmat, 1), size(pmat, 2)))
   do iao = 1, size(pmat, 1)
      do jao = 1, size(pmat, 2)
         scratch(jao, iao) = coeff(jao, iao) * focc(iao)
      end do
   end do
   call gemm(scratch, coeff, pmat, transb='t')
end subroutine get_density_matrix


end module tblite_wavefunction_type
