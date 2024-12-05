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

!> @file tblite/ncoord/exp.f90
!> Provides a coordination number implementation with exponential counting function

!> Coordination number implementation using an exponential counting function, 
!> same as in dftd3.
module tblite_ncoord_exp
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_data_covrad, only : get_covalent_rad
   use tblite_ncoord_type, only : ncoord_type
   implicit none
   private

   public :: new_exp_ncoord

   !> Coordination number evaluator
   type, public, extends(ncoord_type) :: exp_ncoord_type
      real(wp), allocatable :: rcov(:)
   contains
      !> Evaluates the exponential counting function
      procedure :: ncoord_count
      !> Evaluates the derivative of the exponential counting function
      procedure :: ncoord_dcount
   end type exp_ncoord_type

   !> Steepness of counting function
   real(wp), parameter :: default_kcn = 16.0_wp

   real(wp), parameter :: default_cutoff = 25.0_wp

contains


   subroutine new_exp_ncoord(self, mol, kcn, cutoff, rcov)
      !> Coordination number container
      type(exp_ncoord_type), intent(out) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Steepness of counting function
      real(wp), optional :: kcn
      !> Real space cutoff
      real(wp), intent(in), optional :: cutoff
      !> Covalent radii
      real(wp), intent(in), optional :: rcov(:)

      if(present(kcn)) then
         self%kcn = kcn
      else
         self%kcn = default_kcn
      end if

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

      self%directed_factor = 1.0_wp

   end subroutine new_exp_ncoord
   

   !> Exponential counting function for coordination number contributions.
   elemental function ncoord_count(self, mol, izp, jzp, r) result(count)
      !> Coordination number container
      class(exp_ncoord_type), intent(in) :: self
      !> Molecular structure data (not used in exp)
      type(structure_type), intent(in) :: mol
      !> Atom i index
      integer, intent(in)  :: izp
      !> Atom j index
      integer, intent(in)  :: jzp
      !> Current distance.
      real(wp), intent(in) :: r

      real(wp) :: rc, count

      rc = self%rcov(izp) + self%rcov(jzp)
      ! exponential function based counting function
      count =1.0_wp/(1.0_wp+exp(-self%kcn*(rc/r-1.0_wp)))

   end function ncoord_count


   !> Derivative of the exponential counting function w.r.t. the distance.
   elemental function ncoord_dcount(self, mol, izp, jzp, r) result(count)
      !> Coordination number container
      class(exp_ncoord_type), intent(in) :: self
      !> Molecular structure data (not used in gfn)
      type(structure_type), intent(in) :: mol
      !> Atom i index
      integer, intent(in)  :: izp
      !> Atom j index
      integer, intent(in)  :: jzp
      !> Current distance.
      real(wp), intent(in) :: r

      real(wp) :: rc, expterm, count

      rc = self%rcov(izp) + self%rcov(jzp)
      ! exponential function based counting function derivative
      expterm = exp(-self%kcn*(rc/r-1._wp))
      count = (-self%kcn*rc*expterm)/(r**2*((expterm+1._wp)**2))

   end function ncoord_dcount


end module tblite_ncoord_exp
