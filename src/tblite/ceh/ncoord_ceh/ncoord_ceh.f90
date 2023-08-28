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

!> @dir tblite/ceh/ncoord
!> Contains the implementation of the CEH coordination number.

!> @file tblite/ceh/ncoord_ceh/ncoord_ceh.f90
!> Reexports the coordination number evaluation modules.

!> Proxy module to expose coordination number containers
module tblite_ncoord_ceh
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_ncoord_type_ceh, only : ncoord_type_ceh
   use tblite_ncoord_ceh_standard, only : ceh_ncoord_type, new_ceh_ncoord
   implicit none
   private

   public :: ncoord_type_ceh, new_ncoord


contains

!> Create a new generic coordination number container
subroutine new_ncoord(self, mol, cn_type)
   !> Instance of the coordination number container
   class(ncoord_type_ceh), allocatable, intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Coordination number type
   character(len=*), intent(in) :: cn_type

   select case(cn_type)
   case("ceh")
      block
         type(ceh_ncoord_type), allocatable :: tmp
         allocate(tmp)
         call new_ceh_ncoord(tmp, mol)
         call move_alloc(tmp, self)
      end block
   end select
end subroutine new_ncoord

end module tblite_ncoord_ceh
