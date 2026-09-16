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

!> @file tblite/lapack/geqp3.f90
!> Provides wrappers for column-pivoted QR decompositions

!> Wrapper routines for computing a rank-revealing, column-pivoted economy
!> QR decomposition of a general rectangular matrix
module tblite_lapack_geqp3
   use mctc_env, only : sp, dp
   implicit none
   private

   public :: wrap_geqp3

   !> Computes the column-pivoted economy QR decomposition of a real
   !> M-by-N matrix A
   !>
   !>    A * P = Q * R
   !>
   !> where P is the column permutation described by jpvt, Q has
   !> orthonormal columns, and R is upper trapezoidal with |R(1,1)| >=
   !> |R(2,2)| >= ... along columns selected by decreasing pivoted norm.
   !> For an M-by-N matrix, Q has dimensions M-by-min(M,N), while
   !> R has dimensions min(M,N)-by-N.
   interface wrap_geqp3
      module procedure :: wrap_sgeqp3
      module procedure :: wrap_dgeqp3
   end interface wrap_geqp3

   !> Computes a QR factorization with column pivoting of a real
   !> M-by-N matrix A:
   !>
   !>    A * P = Q * ( R ),
   !>                ( 0 )
   !>
   !> where Q is a M-by-M orthogonal matrix, R is upper-triangular, and P
   !> is represented by jpvt: column jpvt(j) of A is moved to position j.
   interface lapack_geqp3
      pure subroutine sgeqp3(m, n, a, lda, jpvt, tau, work, lwork, info)
         import :: sp
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: lwork
         real(sp), intent(inout) :: a(lda, *)
         integer, intent(inout) :: jpvt(*)
         real(sp), intent(out) :: tau(*)
         real(sp), intent(inout) :: work(*)
         integer, intent(out) :: info
      end subroutine sgeqp3

      pure subroutine dgeqp3(m, n, a, lda, jpvt, tau, work, lwork, info)
         import :: dp
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: lwork
         real(dp), intent(inout) :: a(lda, *)
         integer, intent(inout) :: jpvt(*)
         real(dp), intent(out) :: tau(*)
         real(dp), intent(inout) :: work(*)
         integer, intent(out) :: info
      end subroutine dgeqp3
   end interface lapack_geqp3

   !> Generates an M-by-N real matrix Q with orthonormal columns,
   !> which is defined as the first N columns of a product of K elementary
   !> reflectors of order M
   !>
   !>       Q  =  H(1) H(2) . . . H(k)
   !>
   !> as returned by DGEQP3.
   interface lapack_orgqr
      pure subroutine sorgqr(m, n, k, a, lda, tau, work, lwork, info)
         import :: sp
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: lwork
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(in) :: tau(*)
         real(sp), intent(inout) :: work(*)
         integer, intent(out) :: info
      end subroutine sorgqr

      pure subroutine dorgqr(m, n, k, a, lda, tau, work, lwork, info)
         import :: dp
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: lwork
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(in) :: tau(*)
         real(dp), intent(inout) :: work(*)
         integer, intent(out) :: info
      end subroutine dorgqr
   end interface lapack_orgqr

contains

!> Compute a single-precision column-pivoted economy QR decomposition
subroutine wrap_sgeqp3(amat, jpvt, qmat, rmat, info)
   !> Input matrix on entry and QR factorization storage on exit
   real(sp), contiguous, intent(inout) :: amat(:, :)
   !> Permutation so column jpvt(j) of the input ends up in position j
   integer, contiguous, intent(out) :: jpvt(:)
   !> Economy orthogonal matrix Q
   real(sp), contiguous, intent(out) :: qmat(:, :)
   !> Upper trapezoidal matrix R in pivoted column order
   real(sp), contiguous, intent(out) :: rmat(:, :)
   !> Error handling
   integer, intent(out) :: info

   integer :: m, n, k, lda, ldq, lwork, irow, icol
   real(sp), allocatable :: tau(:), work(:)
   real(sp) :: query_geqp3(1), query_orgqr(1)

   m = size(amat, 1)
   n = size(amat, 2)
   k = min(m, n)

   lda = max(1, m)
   ldq = max(1, size(qmat, 1))

   allocate(tau(max(1, k)))

   ! Every column is free to be pivoted
   jpvt(:) = 0

   ! Query the optimal workspace for the pivoted QR factorization
   call lapack_geqp3(m, n, amat, lda, jpvt, tau, query_geqp3, -1, info)
   if (info /= 0) return
   ! ORGQR does not inspect the reflector data during workspace query
   call lapack_orgqr(m, k, k, qmat, ldq, tau, query_orgqr, -1, info)
   if (info /= 0) return

   lwork = max(1, ceiling(query_geqp3(1)), ceiling(query_orgqr(1)))
   allocate(work(lwork))

   ! Compute the column-pivoted QR factorization
   jpvt(:) = 0
   call lapack_geqp3(m, n, amat, lda, jpvt, tau, work, lwork, info)
   if (info /= 0) return

   ! Extract the economy upper trapezoidal matrix R
   rmat(:, :) = 0.0_sp
   do icol = 1, n
      do irow = 1, min(icol, k)
         rmat(irow, icol) = amat(irow, icol)
      end do
   end do

   ! Copy the Householder vectors required to generate the economy Q
   qmat(:, :) = amat(:, 1:k)

   ! Construct the economy orthogonal matrix Q explicitly
   call lapack_orgqr(m, k, k, qmat, ldq, tau, work, lwork, info)

end subroutine wrap_sgeqp3


!> Compute a double-precision column-pivoted economy QR decomposition
subroutine wrap_dgeqp3(amat, jpvt, qmat, rmat, info)
   !> Input matrix on entry and QR factorization storage on exit
   real(dp), contiguous, intent(inout) :: amat(:, :)
   !> Permutation so column jpvt(j) of the input ends up in position j
   integer, contiguous, intent(out) :: jpvt(:)
   !> Economy orthogonal matrix Q
   real(dp), contiguous, intent(out) :: qmat(:, :)
   !> Upper trapezoidal matrix R in pivoted column order
   real(dp), contiguous, intent(out) :: rmat(:, :)
   !> Error handling
   integer, intent(out) :: info

   integer :: m, n, k, lda, ldq, lwork, irow, icol
   real(dp), allocatable :: tau(:), work(:)
   real(dp) :: query_geqp3(1), query_orgqr(1)

   m = size(amat, 1)
   n = size(amat, 2)
   k = min(m, n)

   lda = max(1, m)
   ldq = max(1, size(qmat, 1))

   allocate(tau(max(1, k)))

   ! Every column is free to be pivoted
   jpvt(:) = 0

   ! Query the optimal workspace for the pivoted QR factorization
   call lapack_geqp3(m, n, amat, lda, jpvt, tau, query_geqp3, -1, info)
   if (info /= 0) return
   ! ORGQR does not inspect the reflector data during workspace query
   call lapack_orgqr(m, k, k, qmat, ldq, tau, query_orgqr, -1, info)
   if (info /= 0) return

   lwork = max(1, ceiling(query_geqp3(1)), ceiling(query_orgqr(1)))
   allocate(work(lwork))

   ! Compute the column-pivoted QR factorization
   jpvt(:) = 0
   call lapack_geqp3(m, n, amat, lda, jpvt, tau, work, lwork, info)
   if (info /= 0) return

   ! Extract the economy upper trapezoidal matrix R
   rmat(:, :) = 0.0_dp
   do icol = 1, n
      do irow = 1, min(icol, k)
         rmat(irow, icol) = amat(irow, icol)
      end do
   end do

   ! Copy the Householder vectors required to generate the economy Q
   qmat(:, :) = amat(:, 1:k)

   ! Construct the economy orthogonal matrix Q explicitly
   call lapack_orgqr(m, k, k, qmat, ldq, tau, work, lwork, info)

end subroutine wrap_dgeqp3

end module tblite_lapack_geqp3