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
   use tblite_context_solver, only : context_solver
   use tblite_lapack_sygvd, only : sygvd_solver
   use tblite_lapack_sygvr, only : sygvr_solver
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
   end type enum_lapack

   !> Actual enumerator of possible solvers
   type(enum_lapack), parameter :: lapack_algorithm = enum_lapack()


   !> Generator for LAPACK based electronic solvers
   type, public, extends(context_solver) :: lapack_solver
      !> Selected electronic solver algorithm
      integer :: algorithm = lapack_algorithm%gvd
   contains
      !> Create new instance of electronic solver
      procedure :: new_solver
      !> Delete an electronic solver instance
      procedure :: delete_solver
   end type lapack_solver


contains


!> Create new electronic solver
subroutine new_solver(self, solver, ndim)
   !> Instance of the solver factory
   class(lapack_solver), intent(inout) :: self
   !> New electronic solver
   class(solver_type), allocatable, intent(out) :: solver
   !> Dimension of the eigenvalue problem
   integer, intent(in) :: ndim

   select case(self%algorithm)
   case(lapack_algorithm%gvd)
      solver = sygvd_solver()
   case(lapack_algorithm%gvr)
      solver = sygvr_solver()
   end select
end subroutine new_solver


!> Delete electronic solver instance
subroutine delete_solver(self, solver)
   !> Instance of the solver factory
   class(lapack_solver), intent(inout) :: self
   !> Electronic solver instance
   class(solver_type), allocatable, intent(inout) :: solver

   if (allocated(solver)) deallocate(solver)
end subroutine delete_solver


end module tblite_lapack_solver
