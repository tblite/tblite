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

!> @file tblite_lapack/solver.f90
!> Provides a wrapper for the eigenvalue solvers provided by LAPACK

!> LAPACK based eigenvalue solvers
module tblite_lapack_solver
   use mctc_env, only : wp
   use iso_c_binding, only : c_ptr, c_null_ptr, c_associated, c_size_t
   use tblite_context_solver, only : context_solver
   use tblite_lapack_sygvd, only : sygvd_solver, new_sygvd
   use tblite_lapack_sygvr, only : sygvr_solver, new_sygvr
   use tblite_cusolver_sygvd, only : sygvd_cusolver, new_sygvd_gpu
   use tblite_scf_solver, only : solver_type
   implicit none
   private

   public :: solver_type, lapack_algorithm


   !> Possible solvers provided by LAPACK
   type :: enum_lapack
      !> Divide-and-conquer solver
      integer :: gvd = 1
      !> Relatively robust solver
      integer :: gvr = 2
      !> Divide-and-conquer solver cuSolver implementation
      integer :: gvd_cusolver = 3 
      
   end type enum_lapack

   !> Actual enumerator of possible solvers
   type(enum_lapack), parameter :: lapack_algorithm = enum_lapack()


   !> Generator for LAPACK based electronic solvers
   type, public, extends(context_solver) :: lapack_solver
      !> Selected electronic solver algorithm
      integer :: algorithm = lapack_algorithm%gvd
      !> Pointer to store the C++ solver instance
      type(c_ptr) :: ptr = c_null_ptr
      contains
      !> Create new instance of electronic solver
      procedure :: new
      !> Delete an electronic solver instance
      procedure :: delete
   end type lapack_solver

   interface lapack_solver
      procedure :: new_lapack_solver
   end interface lapack_solver

contains

type(lapack_solver) function new_lapack_solver(algorithm)
   
   !> Selected electronic solver algorithm
   integer, intent(in), optional :: algorithm
   write(*,*) "Creating new LAPACK solver with algorithm: "
   if (present(algorithm)) new_lapack_solver%algorithm = algorithm
   
end function new_lapack_solver


!> Create new electronic solver
subroutine new(self, overlap, nel, kt)
   !> Instance of the solver factory
   class(lapack_solver), intent(inout) :: self
   !> Overlap matrix
   real(wp), intent(in) :: overlap(:, :)
   !> Number of electrons per spin channel
   real(wp), intent(in) :: nel(:)
   !> Electronic temperature
   real(wp), intent(in) :: kt


   if (self%ndim /= size(overlap, 1) .or. .not.(self%reuse)) then
      self%ndim = size(overlap, 1)
      if (allocated(self%solver)) call self%delete()
      select case(self%algorithm)
      case(lapack_algorithm%gvd)
         block
            type(sygvd_solver), allocatable :: tmp
            allocate(tmp)
            call new_sygvd(tmp, overlap, nel, kt)
            call move_alloc(tmp, self%solver)
         end block
      case(lapack_algorithm%gvr)
         block
            type(sygvr_solver), allocatable :: tmp
            allocate(tmp)
            call new_sygvr(tmp, overlap, nel, kt)
            call move_alloc(tmp, self%solver)
         end block
      case(lapack_algorithm%gvd_cusolver)
         block
           
            type(sygvd_cusolver), allocatable :: tmp
            allocate(tmp)
            call new_sygvd_gpu(tmp, overlap, nel, kt, self%ptr)
            call move_alloc(tmp, self%solver)
         end block
      end select
   end if
end subroutine new


!> Delete electronic solver instance
subroutine delete(self)
   !> Instance of the solver factory
   class(lapack_solver), intent(inout) :: self
   

   if (allocated(self%solver)) then
      if (self%reuse .and. c_associated(self%ptr)) then
         call self%solver%delete()
      else
         call self%solver%delete(self%ptr) 
      end if
      deallocate(self%solver)
   end if
end subroutine delete


end module tblite_lapack_solver
