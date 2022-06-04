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

!> @file tblite/scf/solver.f90
!> Provides a base class for defining electronic solvers

!> Declaration of the abstract base class for electronic solvers
module tblite_scf_solver
   use mctc_env, only : sp, dp, error_type
   implicit none
   private

   !> Abstract base class for electronic solvers
   type, public, abstract :: solver_type
   contains
      generic :: solve => solve_sp, solve_dp
      procedure(solve_sp), deferred :: solve_sp
      procedure(solve_dp), deferred :: solve_dp
   end type solver_type

   abstract interface
      subroutine solve_sp(self, hmat, smat, eval, error)
         import :: solver_type, error_type, sp
         class(solver_type), intent(inout) :: self
         real(sp), contiguous, intent(inout) :: hmat(:, :)
         real(sp), contiguous, intent(in) :: smat(:, :)
         real(sp), contiguous, intent(inout) :: eval(:)
         type(error_type), allocatable, intent(out) :: error
      end subroutine solve_sp
      subroutine solve_dp(self, hmat, smat, eval, error)
         import :: solver_type, error_type, dp
         class(solver_type), intent(inout) :: self
         real(dp), contiguous, intent(inout) :: hmat(:, :)
         real(dp), contiguous, intent(in) :: smat(:, :)
         real(dp), contiguous, intent(inout) :: eval(:)
         type(error_type), allocatable, intent(out) :: error
      end subroutine solve_dp
   end interface


end module tblite_scf_solver
