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

!> @dir tblite/api/post_processing.f90
!> Allows API access to post_processing container, to create a new instance from a string.

!> @file tblite/api/post_processing.f90
!> Implements post processing container API access
module tblite_api_post_processing
   use mctc_env, only : error_type, fatal_error
   use tblite_post_processing_list, only : post_processing_list, post_processing_type
   use tblite_api_version, only : namespace
   use tblite_post_processing_list, only : add_post_processing
   use tblite_api_utils, only : c_f_character, f_c_character
   use iso_c_binding
   implicit none
   private

   public :: vp_post_processing, new_post_processing_api, delete_post_processing_api
   public :: push_back_post_processing_api

!> Void pointer to a container instance
   type :: vp_post_processing
      !> Actual container
      type(post_processing_list) :: ptr
   end type vp_post_processing

   logical, parameter :: debug = .true.

contains

function new_post_processing_api() result(vpost_proc)&
      & bind(C, name=namespace//"new_post_processing")
   type(vp_post_processing), pointer :: post_proc
   type(c_ptr) :: vpost_proc
   
   if (debug) print '("[Info]", 1x, a)', "new_post_processing"

   allocate(post_proc)
   
   vpost_proc = c_loc(post_proc)

end function

subroutine push_back_post_processing_api(vpost_proc, charptr) &
      & bind(C, name=namespace//"push_back_post_processing")
   character(kind=c_char), intent(in) :: charptr(*)
   type(c_ptr) :: vpost_proc
   character(len=:), allocatable :: config_str
   class(post_processing_type), allocatable :: pproc
   type(vp_post_processing), pointer :: post_proc
   type(error_type), allocatable :: error

   if (debug) print '("[Info]", 1x, a)', "push_back_post_processing"
   
   call c_f_character(charptr, config_str)

   if (.not.c_associated(vpost_proc)) then 
      call fatal_error(error, "Post Processor container is missing")
      return
   end if
   call c_f_pointer(vpost_proc, post_proc)

   call add_post_processing(post_proc%ptr, config_str, error)

end subroutine

subroutine delete_post_processing_api(vpost_proc) &
      & bind(C, name=namespace//"delete_post_processing")
   type(c_ptr), intent(inout) :: vpost_proc
   type(vp_post_processing), pointer :: post_proc

   if (debug) print '("[Info]", 1x, a)', "delete_post_processing"

   if (c_associated(vpost_proc)) then
      call c_f_pointer(vpost_proc, post_proc)

      deallocate(post_proc)
      vpost_proc = c_null_ptr
   end if
end subroutine delete_post_processing_api

end module
