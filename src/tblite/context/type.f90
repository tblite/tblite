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

!> Calculation context for storing and communicating with the environment
module tblite_context_type
   use iso_fortran_env, only : output_unit
   use mctc_env, only : error_type
   implicit none
   private

   public :: context_type

   type :: context_type
      integer :: unit = output_unit
      integer :: verbosity = 1
      type(error_type), allocatable :: error_log(:)
   contains
      procedure :: message
      procedure :: set_error
      procedure :: get_error
      procedure :: failed
   end type context_type

contains

subroutine set_error(self, error)
   class(context_type), intent(inout) :: self
   type(error_type), intent(in), optional :: error

   if (present(error)) then
      if (.not.allocated(self%error_log)) allocate(self%error_log(0))
      self%error_log = [self%error_log, error]
   end if
end subroutine set_error

subroutine get_error(self, error)
   class(context_type), intent(inout) :: self
   type(error_type), allocatable, intent(out) :: error

   if (.not.allocated(self%error_log)) allocate(self%error_log(0))

   if (size(self%error_log) > 0) then
      error = self%error_log(size(self%error_log))
      self%error_log = self%error_log(:size(self%error_log)-1)
   end if
end subroutine get_error

subroutine message(self, msg)
   class(context_type), intent(inout) :: self
   character(len=*), intent(in) :: msg

   write(self%unit, '(a)') msg
end subroutine message

pure function failed(self)
   class(context_type), intent(in) :: self
   logical :: failed

   failed = .false.
   if (allocated(self%error_log)) then
      failed = size(self%error_log) > 0
   end if
end function failed


end module tblite_context_type
