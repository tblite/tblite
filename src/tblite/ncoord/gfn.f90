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

!> @file tblite/ncoord/gfn.f90
!> Provides a coordination number implementation with double exponential counting function

!> Coordination number implementation using a double exponential counting function
module tblite_ncoord_gfn
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_data_covrad, only : get_covalent_rad
   use tblite_ncoord_type, only : ncoord_type
   implicit none
   private

   public :: new_gfn_ncoord
   public :: get_coordination_number

   !> Coordination number evaluator
   type, public, extends(ncoord_type) :: gfn_ncoord_type
      real(wp) :: cutoff
      real(wp), allocatable :: rcov(:)
   contains
      procedure :: get_cn
   end type gfn_ncoord_type

   !> Steepness of first counting function
   real(wp),parameter :: ka = 10.0_wp
   !> Steepness of second counting function
   real(wp),parameter :: kb = 20.0_wp
   !> Offset of the second counting function
   real(wp),parameter :: r_shift = 2.0_wp

   real(wp), parameter :: default_cutoff = 25.0_wp


contains


subroutine new_gfn_ncoord(self, mol, cutoff, rcov)
   !> Coordination number container
   type(gfn_ncoord_type), intent(out) :: self
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
      self%rcov(:) = get_covalent_rad(mol%num)
   end if
end subroutine new_gfn_ncoord


subroutine get_cn(self, mol, cn, dcndr, dcndL)
   !> Coordination number container
   class(gfn_ncoord_type), intent(in) :: self
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
   call get_coordination_number(mol, lattr, self%cutoff, self%rcov, cn, dcndr, dcndL)
end subroutine get_cn


!> Geometric fractional coordination number, supports exponential counting functions.
subroutine get_coordination_number(mol, trans, cutoff, rcov, cn, dcndr, dcndL)

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

   if (present(dcndr) .and. present(dcndL)) then
      call ncoord_dexp(mol, trans, cutoff, rcov, cn, dcndr, dcndL)
   else
      call ncoord_exp(mol, trans, cutoff, rcov, cn)
   end if

end subroutine get_coordination_number


subroutine ncoord_exp(mol, trans, cutoff, rcov, cn)

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

   !$omp parallel do schedule(runtime) default(none) reduction(+:cn) &
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

            countf = exp_count(ka, r1, rc) * exp_count(kb, r1, rc + r_shift)

            cn(iat) = cn(iat) + countf
            if (iat /= jat) then
               cn(jat) = cn(jat) + countf
            end if

         end do
      end do
   end do

end subroutine ncoord_exp


subroutine ncoord_dexp(mol, trans, cutoff, rcov, cn, dcndr, dcndL)

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
   real(wp), intent(out) :: dcndr(:, :, :)

   !> Derivative of the CN with respect to strain deformations.
   real(wp), intent(out) :: dcndL(:, :, :)

   integer :: iat, jat, izp, jzp, itr
   real(wp) :: r2, r1, rc, rij(3), countf, countd(3), sigma(3, 3), cutoff2

   cn(:) = 0.0_wp
   dcndr(:, :, :) = 0.0_wp
   dcndL(:, :, :) = 0.0_wp
   cutoff2 = cutoff**2

   !$omp parallel do schedule(runtime) default(none) &
   !$omp reduction(+:cn, dcndr, dcndL) shared(mol, trans, cutoff2, rcov) &
   !$omp private(jat, itr, izp, jzp, r2, rij, r1, rc, countf, countd, sigma)
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

            countf = exp_count(ka, r1, rc) * exp_count(kb, r1, rc + r_shift)
            countd = (dexp_count(ka, r1, rc) * exp_count(kb, r1, rc + r_shift) &
               & + exp_count(ka, r1, rc) * dexp_count(kb, r1, rc + r_shift)) * rij/r1

            cn(iat) = cn(iat) + countf
            if (iat /= jat) then
               cn(jat) = cn(jat) + countf
            end if

            dcndr(:, iat, iat) = dcndr(:, iat, iat) + countd
            dcndr(:, jat, jat) = dcndr(:, jat, jat) - countd
            dcndr(:, iat, jat) = dcndr(:, iat, jat) + countd
            dcndr(:, jat, iat) = dcndr(:, jat, iat) - countd

            sigma = spread(countd, 1, 3) * spread(rij, 2, 3)

            dcndL(:, :, iat) = dcndL(:, :, iat) + sigma
            if (iat /= jat) then
               dcndL(:, :, jat) = dcndL(:, :, jat) + sigma
            end if

         end do
      end do
   end do

end subroutine ncoord_dexp


!> Exponential counting function for coordination number contributions.
elemental function exp_count(k, r, r0) result(count)

   !> Steepness of the counting function.
   real(wp), intent(in) :: k

   !> Current distance.
   real(wp), intent(in) :: r

   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp) :: count

   count =1.0_wp/(1.0_wp+exp(-k*(r0/r-1.0_wp)))

end function exp_count


!> Derivative of the counting function w.r.t. the distance.
elemental function dexp_count(k, r, r0) result(count)

   !> Steepness of the counting function.
   real(wp), intent(in) :: k

   !> Current distance.
   real(wp), intent(in) :: r

   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp) :: count
   real(wp) :: expterm

   expterm = exp(-k*(r0/r-1._wp))

   count = (-k*r0*expterm)/(r**2*((expterm+1._wp)**2))

end function dexp_count


end module tblite_ncoord_gfn
