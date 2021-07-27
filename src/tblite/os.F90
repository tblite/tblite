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

module tblite_os
   use, intrinsic :: iso_c_binding, only : c_char, c_null_char, c_int
   implicit none
   private

   public :: setenv
   public :: file_exists, delete_file


   interface
#ifndef _WIN32
      function sys_setenv(name, value, overwrite) result(stat) &
            & bind(c, name="setenv")
         import :: c_char, c_int
         character(len=c_char), intent(in) :: name(*)
         character(len=c_char), intent(in) :: value(*)
         integer(c_int), value :: overwrite
         integer(c_int) :: stat
      end function
#else
      function win_setenv(name, value) result(stat) &
            & bind(c, name="SetEnvironmentVariable")
         import :: c_char, c_int
         character(len=c_char), intent(in) :: name(*)
         character(len=c_char), intent(in) :: value(*)
         integer(c_int) :: stat
      end function
#endif
   end interface


contains


#ifdef _WIN32
function sys_setenv(name, value, overwrite) result(stat)
   character(len=c_char), intent(in) :: name(*)
   character(len=c_char), intent(in) :: value(*)
   integer(c_int), value :: overwrite
   integer(c_int) :: stat

   stat = win_setenv(name, value)
end function
#endif


subroutine setenv(name, value, stat)
   character(len=*), intent(in) :: name
   character(len=*), intent(in) :: value
   integer, intent(out) :: stat
   integer(c_int), parameter :: overwrite = 1_c_int

   stat = sys_setenv(as_c_char(name), as_c_char(value), overwrite)
end subroutine setenv


pure function as_c_char(str) result(res)
   character(len=*), intent(in) :: str
   character(kind=c_char) :: res(len(str)+1)
   res = transfer(str // c_null_char, res)
end function as_c_char


function file_exists(file) result(exist)
   character(len=*), intent(in) :: file
   logical :: exist
   inquire(file=file, exist=exist)
end function file_exists


subroutine delete_file(file)
   character(len=*), intent(in) :: file
   integer :: unit
   if (file_exists(file)) then
      open(file=file, newunit=unit)
      close(unit, status="delete")
   end if
end subroutine delete_file


end module tblite_os
