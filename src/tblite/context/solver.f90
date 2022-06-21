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

!> @file tblite/context/solver.f90
!> Provides an abstract base class for an electronic solver

!> Abstract base class for an electronic solver
module tblite_context_solver
   use tblite_scf_solver, only : solver_type
   implicit none
   private


   !> Abstract base class for creating electronic solver instances
   type, public, abstract :: context_solver
   contains
      !> Create new instance of electronic solver
      procedure(new), deferred :: new
      !> Delete an electronic solver instance
      procedure(delete), deferred :: delete
   end type context_solver


   abstract interface
      !> Create new electronic solver
      subroutine new(self, solver, ndim)
         import :: context_solver, solver_type
         !> Instance of the solver factory
         class(context_solver), intent(inout) :: self
         !> New electronic solver
         class(solver_type), allocatable, intent(out) :: solver
         !> Dimension of the eigenvalue problem
         integer, intent(in) :: ndim
      end subroutine new

      !> Delete electronic solver instance
      subroutine delete(self, solver)
         import :: context_solver, solver_type
         !> Instance of the solver factory
         class(context_solver), intent(inout) :: self
         !> Electronic solver instance
         class(solver_type), allocatable, intent(inout) :: solver
      end subroutine delete
   end interface

end module tblite_context_solver
