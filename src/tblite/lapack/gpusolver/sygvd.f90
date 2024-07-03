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

!> @file tblite/lapack/sygvd.f90
!> Provides an inverface to symmetric divide-and-conquer solver

!> Wrapper to symmetric divide-and-conquer solver for general eigenvalue problems
module tblite_cusolver_sygvd
   use mctc_env, only : sp, dp, error_type, fatal_error, wp
   use tblite_output_format, only : format_string
   use tblite_scf_solver, only : solver_type
   use iso_c_binding
   implicit none
   private

   public :: new_sygvd_gpu


   interface
      type(c_ptr) function cusolver_setup_dp(ndim, coeff_return) bind(C, name="cuSolverSetupDP")
         use iso_c_binding
         integer(c_size_t), value , intent(in) :: ndim
         integer(c_size_t), value :: coeff_return
      end function
      type(c_ptr) function cusolver_setup_sp(ndim, coeff_return) bind(C, name="cuSolverSetupSP")
         use iso_c_binding
         integer(c_size_t), value , intent(in) :: ndim
         integer(c_size_t), value :: coeff_return
      end function
      subroutine cusolve_dp(ptr, Fock, Overlap, orbitalenergies, info) bind(C, name="GPUSolverSolveDP")
         use iso_c_binding
         type(c_ptr), value :: ptr
         real(c_double) :: Fock(*)
         real(c_double) :: Overlap(*)
         real(c_double) :: orbitalenergies(*)
         integer(c_int) :: info
      end subroutine
      subroutine cusolve_sp(ptr, Fock, Overlap, orbitalenergies, info) bind(C, name="GPUSolverSolveSP")
         use iso_c_binding
         type(c_ptr), value :: ptr
         real(c_float) :: Fock(*)
         real(c_float) :: Overlap(*)
         real(c_float) :: orbitalenergies(*)
         integer(c_int) :: info
      end subroutine
      subroutine cusolve_get_density(ptr, coeff, focc, Density) bind(C, name="GPUSolverGetDensityDP")
         use iso_c_binding
         type(c_ptr), value :: ptr
         real(c_double) :: coeff(*)
         real(c_double) :: focc(*)
         real(c_double) :: Density(*)
      end subroutine
      subroutine delete_ptr(ptr) bind(C, name="GPUSolverDelete")
         use iso_c_binding
         type(c_ptr), value :: ptr
      end subroutine
   end interface 


   !> Wrapper class for solving symmetric general eigenvalue problems
   type, public, extends(solver_type) :: sygvd_cusolver
      private
      integer(c_size_t) :: n = 0
      type(c_ptr) :: ptr = c_null_ptr
      integer(c_size_t) :: return_coeff = 0
   contains
      procedure :: solve_sp
      procedure :: solve_dp
      procedure :: delete
      procedure :: get_density_matrix
   end type sygvd_cusolver

contains

subroutine new_sygvd_gpu(self, ndim, return_coeff)
   type(sygvd_cusolver), intent(out) :: self
   integer, intent(in) :: ndim
   logical, intent(in), optional :: return_coeff
   self%n = ndim
   
   if (present(return_coeff)) then
      if (return_coeff) then
         self%return_coeff = 1
      else
         self%return_coeff = 0
      end if
   end if
   
  
end subroutine new_sygvd_gpu

subroutine solve_sp(self, hmat, smat, eval, error)
   class(sygvd_cusolver), intent(inout) :: self
   real(sp), contiguous, intent(inout) :: hmat(:, :)
   real(sp), contiguous, intent(in) :: smat(:, :)
   real(sp), contiguous, intent(inout) :: eval(:)
   type(error_type), allocatable, intent(out) :: error
   integer :: info = 0
   if (.not.c_associated(self%ptr))  self%ptr = cusolver_setup_sp(self%n, self%return_coeff)

   call cusolve_sp(self%ptr, hmat, smat, eval, info)

   call handle_info(error, info)

end subroutine solve_sp

subroutine solve_dp(self, hmat, smat, eval, error)
   class(sygvd_cusolver), intent(inout) :: self
   real(c_double), contiguous, intent(inout) :: hmat(:, :)
   real(c_double), contiguous, intent(in) :: smat(:, :)
   real(c_double), contiguous, intent(inout) :: eval(:)
   type(error_type), allocatable, intent(out) :: error
   integer :: info = 0

   if (.not.c_associated(self%ptr))  self%ptr = cusolver_setup_dp(self%n, self%return_coeff)

   call cusolve_dp(self%ptr, hmat, smat, eval, info)
   
   call handle_info(error, info)

end subroutine solve_dp

subroutine delete(self)
   class(sygvd_cusolver) :: self
   call delete_ptr(self%ptr)

end subroutine delete

subroutine handle_info(error, info)
   type(error_type), allocatable, intent(out) :: error
   integer, intent(in) :: info

   if (info /= 0) then
      call fatal_error(error, "(sygvd) failed to solve eigenvalue problem.&
         & info="//format_string(info, '(i0)'))
   end if
end subroutine handle_info

subroutine get_density_matrix(self, focc, coeff, pmat)
   class(sygvd_cusolver) :: self
   real(wp), intent(in) :: focc(:)
   real(wp), contiguous, intent(in) :: coeff(:, :)
   real(wp), contiguous, intent(out) :: pmat(:, :)

   call cusolve_get_density(self%ptr, coeff, focc, pmat);
   
end subroutine get_density_matrix

end module tblite_cusolver_sygvd
