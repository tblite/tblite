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

!> @file tblite/ncoord/type.f90
!> Provides a coordination number evalulator base class

!> Declaration of base class for coordination number evalulations
module tblite_ncoord_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_cutoff, only : get_lattice_points

   implicit none
   private

   public :: get_coordination_number

   !> Abstract base class for coordination number evaluator
   type, public, abstract :: ncoord_type
      real(wp)  :: cutoff
      !> Steepness of counting function
      real(wp)  :: kcn 
      !> Factor determining whether the CN is evaluated with direction
      !> if +1 the CN contribution is added equally to both partners
      !> if -1 (i.e. with the EN-dep.) it is added to one and subtracted from the other
      real(wp)  :: directed_factor
   contains
      !> Obtains lattice information and calls get_coordination number
      procedure :: get_cn
      !> Decides whether the energy or gradient is calculated
      procedure :: get_coordination_number
      !> Evaluates the CN from the specific counting function
      procedure :: ncoord
      !> Evaluates derivative of the CN from the specific counting function
      procedure :: ncoord_d
      !> Evaluates the counting function (exp, erf, ...)
      procedure(ncoord_count),  deferred :: ncoord_count
      !> Evaluates the derivative of the counting function (exp, erf, ...)
      procedure(ncoord_dcount), deferred :: ncoord_dcount
   end type ncoord_type

   abstract interface

      !> Abstract counting function
      elemental function ncoord_count(self, mol, izp, jzp, r) result(count)
         import :: ncoord_type, structure_type, wp
         !> Instance of coordination number container
         class(ncoord_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Atom i index
         integer, intent(in)  :: izp
         !> Atom j index
         integer, intent(in)  :: jzp
         !> Current distance.
         real(wp), intent(in) :: r

         real(wp) :: count
      end function ncoord_count

      !> Abstract derivative of the counting function w.r.t. the distance.
      elemental function ncoord_dcount(self, mol, izp, jzp, r) result(count)
         import :: ncoord_type, structure_type, wp
         !> Instance of coordination number container
         class(ncoord_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Atom i index
         integer, intent(in)  :: izp
         !> Atom j index
         integer, intent(in)  :: jzp
         !> Current distance.
         real(wp), intent(in) :: r

         real(wp) :: count
      end function ncoord_dcount

   end interface

contains


   subroutine get_cn(self, mol, cn, dcndr, dcndL)
      !> Coordination number container
      class(ncoord_type), intent(in) :: self
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
      call get_coordination_number(self, mol, lattr, self%cutoff, cn, dcndr, dcndL)
   end subroutine get_cn

   !> Geometric fractional coordination number, supports exponential and error counting functions.
   subroutine get_coordination_number(self, mol, trans, cutoff, cn, dcndr, dcndL)

      !> Coordination number container
      class(ncoord_type), intent(in) :: self

      !> Molecular structure data
      type(structure_type), intent(in) :: mol

      !> Lattice points
      real(wp), intent(in) :: trans(:, :)

      !> Real space cutoff
      real(wp), intent(in) :: cutoff

      !> Error function coordination number.
      real(wp), intent(out) :: cn(:)

      !> Derivative of the CN with respect to the Cartesian coordinates.
      real(wp), intent(out), optional :: dcndr(:, :, :)

      !> Derivative of the CN with respect to strain deformations.
      real(wp), intent(out), optional :: dcndL(:, :, :)

      if (present(dcndr) .and. present(dcndL)) then
         call ncoord_d(self, mol, trans, cutoff, cn, dcndr, dcndL)
      else
         call ncoord(self, mol, trans, cutoff, cn)
      end if

   end subroutine get_coordination_number

   subroutine ncoord(self, mol, trans, cutoff, cn)
      !> Coordination number container
      class(ncoord_type), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Lattice points
      real(wp), intent(in) :: trans(:, :)
      !> Real space cutoff
      real(wp), intent(in) :: cutoff
      !> Error function coordination number.
      real(wp), intent(out) :: cn(:)

      integer :: iat, jat, izp, jzp, itr
      real(wp) :: r2, r1, rij(3), countf, cutoff2

      cn(:) = 0.0_wp
      cutoff2 = cutoff**2

      !$omp parallel do schedule(runtime) default(none) reduction(+:cn) &
      !$omp shared(self, mol, trans, cutoff2) &
      !$omp private(jat, itr, izp, jzp, r2, rij, r1, countf)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)

            do itr = 1, size(trans, dim=2)
               rij = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, itr))
               r2 = sum(rij**2)
               if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
               r1 = sqrt(r2)

               countf = self%ncoord_count(mol, izp, jzp, r1)

               cn(iat) = cn(iat) + countf
               if (iat /= jat) then
                  cn(jat) = cn(jat) + countf * self%directed_factor
               end if

            end do
         end do
      end do

   end subroutine ncoord

   subroutine ncoord_d(self, mol, trans, cutoff, cn, dcndr, dcndL)
      !> Coordination number container
      class(ncoord_type), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Lattice points
      real(wp), intent(in) :: trans(:, :)
      !> Real space cutoff
      real(wp), intent(in) :: cutoff
      !> Error function coordination number.
      real(wp), intent(out) :: cn(:)
      !> Derivative of the CN with respect to the Cartesian coordinates.
      real(wp), intent(out) :: dcndr(:, :, :)
      !> Derivative of the CN with respect to strain deformations.
      real(wp), intent(out) :: dcndL(:, :, :)

      integer :: iat, jat, izp, jzp, itr
      real(wp) :: r2, r1, rij(3), countf, countd(3), sigma(3, 3), cutoff2

      cn(:) = 0.0_wp
      dcndr(:, :, :) = 0.0_wp
      dcndL(:, :, :) = 0.0_wp
      cutoff2 = cutoff**2

      !$omp parallel do schedule(runtime) default(none) &
      !$omp reduction(+:cn, dcndr, dcndL) shared(mol, trans, cutoff2) &
      !$omp shared(self) &
      !$omp private(jat, itr, izp, jzp, r2, rij, r1, countf, countd, sigma)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)

            do itr = 1, size(trans, dim=2)
               rij = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, itr))
               r2 = sum(rij**2)
               if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
               r1 = sqrt(r2)

               countf = self%ncoord_count(mol, izp, jzp, r1)
               countd = self%ncoord_dcount(mol, izp, jzp, r1) * rij/r1

               cn(iat) = cn(iat) + countf
               if (iat /= jat) then
                  cn(jat) = cn(jat) + countf * self%directed_factor
               end if

               dcndr(:, iat, iat) = dcndr(:, iat, iat) + countd
               dcndr(:, jat, jat) = dcndr(:, jat, jat) - countd * self%directed_factor
               dcndr(:, iat, jat) = dcndr(:, iat, jat) + countd * self%directed_factor
               dcndr(:, jat, iat) = dcndr(:, jat, iat) - countd

               sigma = spread(countd, 1, 3) * spread(rij, 2, 3)

               dcndL(:, :, iat) = dcndL(:, :, iat) + sigma
               if (iat /= jat) then
                  dcndL(:, :, jat) = dcndL(:, :, jat) + sigma * self%directed_factor
               end if

            end do
         end do
      end do

   end subroutine ncoord_d


end module tblite_ncoord_type
