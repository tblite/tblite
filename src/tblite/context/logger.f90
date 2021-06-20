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

!> Logger to display strings
module tblite_context_logger
   implicit none
   private

   public :: context_logger


   type, abstract :: context_logger
   contains
      procedure(message), deferred :: message
   end type context_logger


   abstract interface
      subroutine message(self, msg)
         import :: context_logger
         class(context_logger), intent(inout) :: self
         character(len=*), intent(in) :: msg
      end subroutine message
   end interface


end module tblite_context_logger
