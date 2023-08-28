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

!> @file tblite/ncoord/ceh_std.f90
!> Provides a (standard) coordination number implementation for the CEH method

!> Coordination number implementation for the CEH method in the usual form.
module tblite_ncoord_ceh_std
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_data_covrad_ceh, only : get_covalent_cehrad
   use tblite_ncoord_type, only : ncoord_type
   implicit none
   private

   public :: new_ceh_std_ncoord
   public :: get_coordination_number

   !> Coordination number evaluator
   type, public, extends(ncoord_type) :: ceh_std_ncoord_type
      real(wp) :: cutoff
      real(wp), allocatable :: rcov(:)
   contains
      procedure :: get_cn
   end type ceh_std_ncoord_type

   real(wp), parameter :: default_cutoff = 25.0_wp

   !> Steepness of counting function
   real(wp), parameter :: kcn = 3.09_wp

contains

   subroutine new_ceh_std_ncoord(self, mol, cutoff, rcov)
      !> Coordination number container
      type(ceh_std_ncoord_type), intent(out) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Real space cutoff
      real(wp), intent(in), optional :: cutoff
      !> Covalent radii
      real(wp), intent(in), optional :: rcov(:)

      if (present(cutoff)) then
         self%cutoff = cutoff
      else
         self%cutoff = default_cutoff
      end if

      allocate(self%rcov(mol%nid))
      if (present(rcov)) then
         self%rcov(:) = rcov
      else
         self%rcov(:) = get_covalent_cehrad(mol%num)
      end if
   end subroutine new_ceh_std_ncoord


   subroutine get_cn(self, mol, cn, dcndr, dcndL)
      !> Coordination number container
      class(ceh_std_ncoord_type), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Error function coordination number.
      real(wp), intent(out) :: cn(:)
      !> Derivative of the CN with respect to the Cartesian coordinates.
      real(wp), intent(out), optional :: dcndr(:, :, :)
      !> Derivative of the CN with respect to strain deformations.
      real(wp), intent(out), optional :: dcndL(:, :, :)

      real(wp), allocatable :: lattr(:, :)

      call get_lattice_points(mol%periodic, mol%lattice, self%cutoff, lattr)
      call get_coordination_number(mol, lattr, self%cutoff, self%rcov, cn, &
      & dcndr, dcndL)
   end subroutine get_cn


!> Geometric fractional coordination number, supports exponential counting functions.
   subroutine get_coordination_number(mol, trans, cutoff, rcov, cn, &
      & dcndr, dcndL)

      !> Molecular structure data
      type(structure_type), intent(in) :: mol

      !> Lattice points
      real(wp), intent(in) :: trans(:, :)

      !> Real space cutoff
      real(wp), intent(in) :: cutoff

      !> Covalent radius
      real(wp), intent(in) :: rcov(:)

      !> Error function coordination number.
      real(wp), intent(out) :: cn(:)

      !> Derivative of the CN with respect to the Cartesian coordinates.
      real(wp), intent(out), optional :: dcndr(:, :, :)

      !> Derivative of the CN with respect to strain deformations.
      real(wp), intent(out), optional :: dcndL(:, :, :)

      call ncoord_erf(mol, trans, cutoff, rcov, cn)

   end subroutine get_coordination_number


   subroutine ncoord_erf(mol, trans, cutoff, rcov, cn)

      !> Molecular structure data
      type(structure_type), intent(in) :: mol

      !> Lattice points
      real(wp), intent(in) :: trans(:, :)

      !> Real space cutoff
      real(wp), intent(in) :: cutoff

      !> Covalent radius
      real(wp), intent(in) :: rcov(:)

      !> Error function coordination number.
      real(wp), intent(out) :: cn(:)

      integer :: iat, jat, izp, jzp, itr
      real(wp) :: r2, r1, rc, rij(3), countf, cutoff2

      cn(:) = 0.0_wp
      cutoff2 = cutoff**2

      !$omp parallel do default(none) reduction(+:cn) &
      !$omp shared(mol, trans, cutoff2, rcov) &
      !$omp private(jat, itr, izp, jzp, r2, rij, r1, rc, countf)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)
            do itr = 1, size(trans, dim=2)
               rij = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, itr))
               r2 = sum(rij**2)
               if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
               r1 = sqrt(r2)

               rc = rcov(izp) + rcov(jzp)

               countf = erf_count(kcn, r1, rc)

               cn(iat) = cn(iat) + countf
               if (iat /= jat) then
                  cn(jat) = cn(jat) + countf
               end if
            end do
         end do
      end do

   end subroutine ncoord_erf

!> Error function counting function for coordination number contributions.
   pure function erf_count(k, r, r0) result(count)

      !> Steepness of the counting function.
      real(wp), intent(in) :: k

      !> Current distance.
      real(wp), intent(in) :: r

      !> Sum of covalent radii.
      real(wp), intent(in) :: r0

      real(wp) :: count

      count = 0.5_wp * (1.0_wp + erf(-k*(r-r0)/r0))

   end function erf_count

end module tblite_ncoord_ceh_std
