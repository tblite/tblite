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
   use tblite_data_covrad, only : get_covalent_rad
   use tblite_ncoord_type, only : ncoord_type
   implicit none
   private

   public :: new_gfn_ncoord
   

   !> Coordination number evaluator
   type, public, extends(ncoord_type) :: gfn_ncoord_type
      real(wp), allocatable :: rcov(:)
   contains
      !> Evaluates the gfn counting function 
      procedure :: ncoord_count
      !> Evaluates the derivative of the gfn counting function
      procedure :: ncoord_dcount
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
   
   self%directed_factor = 1.0_wp

end subroutine new_gfn_ncoord

!> Double-exponential counting function for coordination number contributions.
elemental function ncoord_count(self, mol, izp, jzp, r) result(count)
   !> Coordination number container
   class(gfn_ncoord_type), intent(in) :: self
   !> Molecular structure data (not used in gfn)
   type(structure_type), intent(in) :: mol
   !> Atom i index
   integer, intent(in)  :: izp
   !> Atom j index
   integer, intent(in)  :: jzp
   !> Current distance.
   real(wp), intent(in) :: r

   real(wp) :: rc, count

   rc = self%rcov(izp) + self%rcov(jzp)
   ! double exponential function based counting function
   count = gfn_count(ka, r, rc) * gfn_count(kb, r, rc + r_shift)

end function ncoord_count

!> Derivative of the double-exponential counting function w.r.t. the distance.
elemental function ncoord_dcount(self, mol, izp, jzp, r) result(count)
   !> Coordination number container
   class(gfn_ncoord_type), intent(in) :: self
   !> Molecular structure data (not used in gfn)
   type(structure_type), intent(in) :: mol
   !> Atom i index
   integer, intent(in)  :: izp
   !> Atom j index
   integer, intent(in)  :: jzp
   !> Current distance.
   real(wp), intent(in) :: r

   real(wp) :: rc, count
   
   rc = self%rcov(izp) + self%rcov(jzp)
   ! double exponential function based counting function derivative
   count = (gfn_dcount(ka, r, rc) * gfn_count(kb, r, rc + r_shift) &
               & + gfn_count(ka, r, rc) * gfn_dcount(kb, r, rc + r_shift))

end function ncoord_dcount

!> Mono-exponential counting function for coordination number contributions.
elemental function gfn_count(k, r, r0) result(count)
   !> Steepness of the counting function.
   real(wp), intent(in) :: k
   !> Current distance.
   real(wp), intent(in) :: r
   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp) :: count

   count = 1.0_wp/(1.0_wp+exp(-k*(r0/r-1.0_wp)))

end function gfn_count


!> Derivative of the mono-exponential counting function w.r.t. the distance.
elemental function gfn_dcount(k, r, r0) result(count)
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

end function gfn_dcount

end module tblite_ncoord_gfn
