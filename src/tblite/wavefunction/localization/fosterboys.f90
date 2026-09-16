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

!> @file tblite/wavefunction/localization/fosterboys.f90
!> Implements the Foster-Boys orbital localization method used by xTB
!>
!> The localization criterion follows:
!>
!> S. F. Boys, "Construction of Some Molecular Orbitals to Be Approximately
!> Invariant for Changes from One Molecule to Another",
!> Rev. Mod. Phys. 32, 296 (1960). DOI: 10.1103/RevModPhys.32.296
!>
!> J. M. Foster, S. F. Boys, "Canonical Configurational Interaction Procedure",
!> Rev. Mod. Phys. 32, 300 (1960). DOI: 10.1103/RevModPhys.32.300
!>
!> Rotates occupied orbitals of each spin channel to maximize charge centers separation.
module tblite_wavefunction_localization_fosterboys
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_blas, only : gemm, swap
   use tblite_wavefunction_localization_jacobi, only : new_jacobi_optimizer
   use tblite_wavefunction_localization_type, only : localization_type, &
      & get_common_origin_position
   implicit none
   private

   public :: fosterboys_localization_type, new_fosterboys_localization

   !> Foster-Boys localization
   type, extends(localization_type) :: fosterboys_localization_type
   contains
      !> Prepare common-origin position operator projected into the occupied space
      procedure :: prepare
      !> Reorder the optimized transformation by the localized orbital energy
      procedure :: postprocess
   end type fosterboys_localization_type

contains


!> Create a new Foster-Boys orbital localization method
subroutine new_fosterboys_localization(self)
   !> Instance of the Foster-Boys localization method
   type(fosterboys_localization_type), intent(out) :: self

   call new_jacobi_optimizer(self%optimizer)

end subroutine new_fosterboys_localization


!> Prepare the Foster-Boys localization objective for one spin channel
subroutine prepare(self, mol, bas, overlap, dipole, coeff_occ, opmat, error)
   !> Instance of the Foster-Boys localization method
   class(fosterboys_localization_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Overlap integrals
   real(wp), intent(in) :: overlap(:, :)
   !> Dipole integrals with moment operator centered on last index
   real(wp), intent(in) :: dipole(:, :, :)
   !> Canonical occupied orbital coefficients
   real(wp), intent(in) :: coeff_occ(:, :)
   !> Common-origin position operators projected into the occupied space
   real(wp), allocatable, intent(out) :: opmat(:, :, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), allocatable :: position(:, :, :), tmp(:, :)
   integer :: nao, nocc, k

   if (allocated(mol%periodic)) then
      if (any(mol%periodic)) then
         call fatal_error(error, "Foster-Boys orbital localization is not "// &
            & "implemented for periodic systems")
         return
      end if
   end if

   ! Build the common-origin position operator
   call get_common_origin_position(mol, bas, overlap, dipole, position, error)
   if (allocated(error)) return

   nao = size(coeff_occ, 1)
   nocc = size(coeff_occ, 2)

   allocate(opmat(nocc, nocc, 3), tmp(nao, nocc))

   ! Project the common-origin position operators into the occupied MO space
   do k = 1, 3
      call gemm(position(:, :, k), coeff_occ, tmp)
      call gemm(coeff_occ, tmp, opmat(:, :, k), transa="T")
   end do

end subroutine prepare


!> Sort the localized orbitals according to their localized (rotated) orbital energy
subroutine postprocess(self, emo_occ, trafo)
   !> Instance of the Foster-Boys localization method
   class(fosterboys_localization_type), intent(in) :: self
   !> Canonical occupied orbital energies
   real(wp), intent(in) :: emo_occ(:)
   !> Occupied-space orthogonal transformation
   real(wp), intent(inout) :: trafo(:, :)

   real(wp), allocatable :: eps_local(:)
   integer :: nocc, i, j, k
   real(wp) :: eps_k

   nocc = size(emo_occ)
   allocate(eps_local(nocc), source=0.0_wp)
   do k = 1, nocc
      eps_local(:) = eps_local(:) + emo_occ(k) * trafo(k, :)**2
   end do

   ! Reorder the localized orbitals by their localized orbital energy
   nocc = size(eps_local)
   do i = 1, nocc - 1
      k = i
      eps_k = eps_local(i)
      do j = i + 1, nocc
         if (eps_local(j) < eps_k) then
            k = j
            eps_k = eps_local(j)
         end if
      end do
      if (k == i) cycle

      eps_local(k) = eps_local(i)
      eps_local(i) = eps_k
      call swap(trafo(:, i), trafo(:, k))
   end do

end subroutine postprocess

end module tblite_wavefunction_localization_fosterboys
