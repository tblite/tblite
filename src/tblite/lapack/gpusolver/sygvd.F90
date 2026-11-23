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

!> @file tblite/lapack/gpusolver/sygvd.F90
!> Provides a generalized symmetric eigensolver backed by cuSOLVER.

!> Wrapper to symmetric divide-and-conquer solver for general eigenvalue problems cuSOLVER implementation
module tblite_cusolver_sygvd
   use, intrinsic :: iso_c_binding, only : c_associated, c_double, c_int, c_null_ptr, c_ptr
   use mctc_env, only : sp, dp, error_type, fatal_error, wp
   use tblite_output_format, only : format_string
   use tblite_scf_diag, only : diag_solver_type
   implicit none
   private

   public :: new_sygvd_gpu

#ifdef USE_CUDA
   interface
      integer(c_int) function cuda_sygvd(ptr, n, single, hmat, smat, eval) bind(C, name="tblite_cuda_sygvd")
         import :: c_int, c_ptr
         type(c_ptr), intent(inout) :: ptr
         integer(c_int), value :: n, single
         type(*), intent(inout) :: hmat(*), eval(*)
         type(*), intent(in) :: smat(*)
      end function cuda_sygvd
      integer(c_int) function cuda_density(ptr, coeff, focc, density) bind(C, name="tblite_cuda_density")
         import :: c_double, c_int, c_ptr
         type(c_ptr), value :: ptr
         real(c_double), intent(in) :: coeff(*), focc(*)
         real(c_double), intent(out) :: density(*)
      end function cuda_density
      subroutine delete_cusolver_ptr(ptr) bind(C, name="tblite_cuda_delete")
         import :: c_ptr
         type(c_ptr), value :: ptr
      end subroutine delete_cusolver_ptr
   end interface


#endif

   type, public, extends(diag_solver_type) :: sygvd_cusolver
      private
      type(c_ptr) :: double_solver = c_null_ptr
      type(c_ptr) :: single_solver = c_null_ptr
   contains
      procedure :: solve_sp
      procedure :: solve_dp
#ifdef USE_CUDA
      procedure :: get_density_matrix
#endif
      procedure :: delete
      final :: finalize
   end type sygvd_cusolver

contains

subroutine new_sygvd_gpu(self, overlap, nel, kt)
   type(sygvd_cusolver), intent(out) :: self
   real(wp), intent(in) :: overlap(:, :), nel(:), kt

   self%nel = nel
   self%kt = kt
end subroutine new_sygvd_gpu

subroutine solve_sp(self, hmat, smat, eval, error)
   class(sygvd_cusolver), intent(inout) :: self
   real(sp), contiguous, intent(inout) :: hmat(:, :), eval(:)
   real(sp), contiguous, intent(in) :: smat(:, :)
   type(error_type), allocatable, intent(out) :: error
#ifdef USE_CUDA
   integer(c_int) :: info

   info = cuda_sygvd(self%single_solver, size(hmat, 1, kind=c_int), 1_c_int, hmat, smat, eval)
   if (info /= 0) call fatal_error(error, "cuSOLVER generalized diagonalization failed, info="//format_string(info, "(i0)"))
#else
   call fatal_error(error, "This tblite build does not support GPU acceleration.")
#endif
end subroutine solve_sp

subroutine solve_dp(self, hmat, smat, eval, error)
   class(sygvd_cusolver), intent(inout) :: self
   real(dp), contiguous, intent(inout) :: hmat(:, :), eval(:)
   real(dp), contiguous, intent(in) :: smat(:, :)
   type(error_type), allocatable, intent(out) :: error
#ifdef USE_CUDA
   integer(c_int) :: info

   info = cuda_sygvd(self%double_solver, size(hmat, 1, kind=c_int), 0_c_int, hmat, smat, eval)
   if (info /= 0) call fatal_error(error, "cuSOLVER generalized diagonalization failed, info="//format_string(info, "(i0)"))
#else
   call fatal_error(error, "This tblite build does not support GPU acceleration.")
#endif
end subroutine solve_dp

#ifdef USE_CUDA
subroutine get_density_matrix(self, focc, coeff, pmat, error)
   class(sygvd_cusolver), intent(inout) :: self
   real(wp), intent(in) :: focc(:)
   real(wp), contiguous, intent(in) :: coeff(:, :)
   real(wp), contiguous, intent(inout) :: pmat(:, :)
   type(error_type), allocatable, intent(inout) :: error
   integer(c_int) :: info

   info = cuda_density(self%double_solver, coeff, focc, pmat)
   if (info /= 0) call fatal_error(error, "GPU density construction failed, info="//format_string(info, "(i0)"))
end subroutine get_density_matrix
#endif

subroutine delete(self)
   class(sygvd_cusolver), intent(inout) :: self
#ifdef USE_CUDA
   if (c_associated(self%double_solver)) call delete_cusolver_ptr(self%double_solver)
   if (c_associated(self%single_solver)) call delete_cusolver_ptr(self%single_solver)
#endif
   self%double_solver = c_null_ptr
   self%single_solver = c_null_ptr
end subroutine delete

subroutine finalize(self)
   type(sygvd_cusolver), intent(inout) :: self
   call self%delete()
end subroutine finalize

end module tblite_cusolver_sygvd
