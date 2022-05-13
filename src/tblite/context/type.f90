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
   use tblite_context_logger, only : context_logger
   use tblite_context_terminal, only : context_terminal
   implicit none
   private

   public :: context_type


   !> Calculation context type for error handling and output messages
   type :: context_type
      !> Default output unit for this context
      integer :: unit = output_unit
      !> Default verbosity for procedures using this context
      integer :: verbosity = 1
      !> Stack containing the error messages of this context
      type(error_type), allocatable :: error_log(:)
      !> Optional logger to be used for writing messages
      class(context_logger), allocatable :: io
      !> Color support for output
      type(context_terminal) :: terminal = context_terminal()
   contains
      !> Write a message to the output
      procedure :: message
      !> Push an error message to the context
      procedure :: set_error
      !> Pop an error message from the context
      procedure :: get_error
      !> Query the context for errors
      procedure :: failed
   end type context_type


contains


!> Add an error message to the context
subroutine set_error(self, error)
   !> Instance of the calculation context
   class(context_type), intent(inout) :: self
   !> Error handling
   type(error_type), intent(in), optional :: error

   if (present(error)) then
      if (.not.allocated(self%error_log)) allocate(self%error_log(0))
      self%error_log = [self%error_log, error]
   end if
end subroutine set_error


!> Pop an error message from the context
subroutine get_error(self, error)
   !> Instance of the calculation context
   class(context_type), intent(inout) :: self
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   if (.not.allocated(self%error_log)) allocate(self%error_log(0))

   if (size(self%error_log) > 0) then
      error = self%error_log(size(self%error_log))
      self%error_log = self%error_log(:size(self%error_log)-1)
   end if
end subroutine get_error


!> Write a message to the output
subroutine message(self, msg)
   !> Instance of the calculation context
   class(context_type), intent(inout) :: self
   !> Message to write
   character(len=*), intent(in) :: msg

   if (allocated(self%io)) then
      call self%io%message(msg)
   else
      write(self%unit, '(a)') msg
   end if
end subroutine message


!> Query the context for errors
pure function failed(self)
   !> Instance of the calculation context
   class(context_type), intent(in) :: self
   !> Error status of the context
   logical :: failed

   failed = .false.
   if (allocated(self%error_log)) then
      failed = size(self%error_log) > 0
   end if
end function failed


end module tblite_context_type
