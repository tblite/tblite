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

!> @file tblite/ncoord_ceh/type_ceh.f90
!> Provides a coordination number evalulator base class

!> Declaration of base class for coordination number evalulations
module tblite_ncoord_type_ceh
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none
   private

   !> Abstract base class for coordination number evaluator
   type, public, abstract :: ncoord_type_ceh
   contains
      procedure(get_cn), deferred :: get_cn
   end type ncoord_type_ceh

   abstract interface
      subroutine get_cn(self, mol, cn, cn_en, dcndr, dcndL, dcnendr, dcnendL)
         import :: ncoord_type_ceh, structure_type, wp
         !> Coordination number container
         class(ncoord_type_ceh), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Error function coordination number.
         real(wp), intent(out) :: cn(:)
         !> Error function coordination number.
         real(wp), intent(out) :: cn_en(:)
         !> Derivative of the CN with respect to the Cartesian coordinates.
         real(wp), intent(out), optional :: dcndr(:, :, :)
         !> Derivative of the CN with respect to strain deformations.
         real(wp), intent(out), optional :: dcndL(:, :, :)
         !> Derivative of the CN_EN with respect to the Cartesian coordinates.
         real(wp), intent(out), optional :: dcnendr(:, :, :)
         !> Derivative of the CN_EN with respect to strain deformations.
         real(wp), intent(out), optional :: dcnendL(:, :, :)
      end subroutine get_cn
   end interface

end module tblite_ncoord_type_ceh
