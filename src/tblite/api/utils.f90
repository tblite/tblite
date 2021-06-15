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

module tblite_api_utils
   use, intrinsic :: iso_c_binding
   implicit none
   private

   public :: f_c_character


contains


subroutine f_c_character(rhs, lhs, len)
   character(kind=c_char), intent(out) :: lhs(*)
   character(len=*), intent(in) :: rhs
   integer, intent(in) :: len
   integer :: length
   length = min(len-1, len_trim(rhs))

   lhs(1:length) = transfer(rhs(1:length), lhs(1:length))
   lhs(length+1:length+1) = c_null_char

end subroutine f_c_character

end module tblite_api_utils
