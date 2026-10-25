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
module tblite_lapack_sygvd
   use mctc_env, only : sp, dp, wp, error_type, fatal_error
   use tblite_blas_level3, only : wrap_trsm
   use tblite_lapack_potrf, only : update_potrf
   use tblite_lapack_sygst, only : wrap_sygst
   use tblite_output_format, only : format_string
   use tblite_scf_diag, only : diag_solver_type
   implicit none
   private

   public :: new_sygvd


   interface lapack_syevd
      pure subroutine ssyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)
         import :: sp
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(out) :: w(*)
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
         integer, intent(out) :: info
      end subroutine ssyevd
      pure subroutine dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)
         import :: dp
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(out) :: w(*)
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
         integer, intent(out) :: info
      end subroutine dsyevd
   end interface lapack_syevd


   !> Wrapper class for solving symmetric general eigenvalue problems
   type, public, extends(diag_solver_type) :: sygvd_solver
      private
      integer :: n = 0
      integer, allocatable :: iwork(:)
      real(sp), allocatable :: swork(:)
      real(sp), allocatable :: soverlap(:)
      real(sp), allocatable :: schol(:, :)
      real(dp), allocatable :: dwork(:)
      real(dp), allocatable :: doverlap(:)
      real(dp), allocatable :: dchol(:, :)
   contains
      procedure :: solve_sp
      procedure :: solve_dp
   end type sygvd_solver

contains

subroutine new_sygvd(self, overlap, nel, kt)
   type(sygvd_solver), intent(out) :: self
   real(wp), intent(in) :: overlap(:, :)
   real(wp), intent(in) :: nel(:)
   real(wp), intent(in) :: kt
   integer :: info
   self%n = size(overlap, 1)
   self%nel = nel
   self%kt = kt
   call update_potrf(overlap, self%doverlap, self%dchol, info)
end subroutine new_sygvd

subroutine solve_sp(self, hmat, smat, eval, error)
   class(sygvd_solver), intent(inout) :: self
   real(sp), contiguous, intent(inout) :: hmat(:, :)
   real(sp), contiguous, intent(in) :: smat(:, :)
   real(sp), contiguous, intent(inout) :: eval(:)
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: itype = 1
   character(len=1), parameter :: uplo = "U"
   integer :: info, lswork, liwork

   if (self%n /= size(hmat, 1)) then
      self%n = size(hmat, 1)
      if (allocated(self%swork)) deallocate(self%swork)
      if (allocated(self%dwork)) deallocate(self%dwork)
      if (allocated(self%iwork)) deallocate(self%iwork)
   end if
   if (.not.allocated(self%swork)) then
      allocate(self%swork(1 + 6*self%n + 2*self%n**2))
   end if
   if (.not.allocated(self%iwork)) then
      allocate(self%iwork(3 + 5*self%n))
   end if
   lswork = size(self%swork)
   liwork = size(self%iwork)
   info = 0

   call update_potrf(smat, self%soverlap, self%schol, info)
   call handle_info(error, info)
   if (allocated(error)) return

   call wrap_sygst(hmat, self%schol, info, itype=itype, uplo=uplo)
   call handle_info(error, info)
   if (allocated(error)) return

   call lapack_syevd("V", uplo, self%n, hmat, max(1, self%n), eval, self%swork, lswork, &
      & self%iwork, liwork, info)
   call handle_info(error, info)
   if (allocated(error)) return

   if (self%n > 0) then
      call wrap_trsm(self%schol, hmat, side="L", uplo=uplo, transa="N", diag="N")
   end if

end subroutine solve_sp

subroutine solve_dp(self, hmat, smat, eval, error)
   class(sygvd_solver), intent(inout) :: self
   real(dp), contiguous, intent(inout) :: hmat(:, :)
   real(dp), contiguous, intent(in) :: smat(:, :)
   real(dp), contiguous, intent(inout) :: eval(:)
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: itype = 1
   character(len=1), parameter :: uplo = "U"
   integer :: info, ldwork, liwork

   if (self%n /= size(hmat, 1)) then
      self%n = size(hmat, 1)
      if (allocated(self%swork)) deallocate(self%swork)
      if (allocated(self%dwork)) deallocate(self%dwork)
      if (allocated(self%iwork)) deallocate(self%iwork)
   end if
   if (.not.allocated(self%dwork)) then
      allocate(self%dwork(1 + 6*self%n + 2*self%n**2))
   end if
   if (.not.allocated(self%iwork)) then
      allocate(self%iwork(3 + 5*self%n))
   end if
   ldwork = size(self%dwork)
   liwork = size(self%iwork)
   info = 0

   call update_potrf(smat, self%doverlap, self%dchol, info)
   call handle_info(error, info)
   if (allocated(error)) return

   call wrap_sygst(hmat, self%dchol, info, itype=itype, uplo=uplo)
   call handle_info(error, info)
   if (allocated(error)) return

   call lapack_syevd("V", uplo, self%n, hmat, max(1, self%n), eval, self%dwork, ldwork, &
      & self%iwork, liwork, info)
   call handle_info(error, info)
   if (allocated(error)) return

   if (self%n > 0) then
      call wrap_trsm(self%dchol, hmat, side="L", uplo=uplo, transa="N", diag="N")
   end if

end subroutine solve_dp

subroutine handle_info(error, info)
   type(error_type), allocatable, intent(out) :: error
   integer, intent(in) :: info

   if (info /= 0) then
      call fatal_error(error, "(sygvd) failed to solve eigenvalue problem.&
         & info="//format_string(info, "(i0)"))
   end if
end subroutine handle_info

end module tblite_lapack_sygvd
