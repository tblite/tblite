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

!> @file tblite/lapack/getri.f90
!> Provides wrappers for computing a matrix inverse

!> Wrappers to obtain the inverse of a matrix

#ifndef IK
#define IK i4
#endif

module tblite_lapack_getri
   use mctc_env, only : sp, dp, ik => IK
   implicit none
   private

   public :: wrap_getri


   !> Computes the inverse of a matrix using the LU factorization
   !> computed by ?GETRF.
   !>
   !> This method inverts U and then computes inv(A) by solving the system
   !> inv(A)*L = inv(U) for inv(A).
   interface wrap_getri
      module procedure :: wrap_sgetri
      module procedure :: wrap_dgetri
   end interface wrap_getri


   !> Computes the inverse of a matrix using the LU factorization
   !> computed by ?GETRF.
   !>
   !> This method inverts U and then computes inv(A) by solving the system
   !> inv(A)*L = inv(U) for inv(A).
   interface lapack_getri
      pure subroutine sgetri(n, a, lda, ipiv, work, lwork, info)
         import :: sp, ik
         integer(ik), intent(in) :: lda
         real(sp), intent(inout) :: a(lda, *)
         integer(ik), intent(in) :: ipiv(*)
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         real(sp), intent(inout) :: work(*)
         integer(ik), intent(in) :: lwork
      end subroutine sgetri
      pure subroutine dgetri(n, a, lda, ipiv, work, lwork, info)
         import :: dp, ik
         integer(ik), intent(in) :: lda
         real(dp), intent(inout) :: a(lda, *)
         integer(ik), intent(in) :: ipiv(*)
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         real(dp), intent(inout) :: work(*)
         integer(ik), intent(in) :: lwork
      end subroutine dgetri
   end interface lapack_getri

contains

subroutine wrap_sgetri(amat, ipiv, info)
   real(sp), intent(inout) :: amat(:, :)
   integer(ik), intent(in) :: ipiv(:)
   integer(ik), intent(out) :: info
   integer(ik) :: n, lda, lwork
   real(sp), allocatable :: work(:)
   real(sp) :: test(1)
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   lwork = -1
   call lapack_getri(n, amat, lda, ipiv, test, lwork, info)
   if (info == 0) then
      lwork = nint(test(1))
      allocate(work(lwork))
      call lapack_getri(n, amat, lda, ipiv, work, lwork, info)
   end if
end subroutine wrap_sgetri


subroutine wrap_dgetri(amat, ipiv, info)
   real(dp), intent(inout) :: amat(:, :)
   integer(ik), intent(in) :: ipiv(:)
   integer(ik), intent(out) :: info
   integer(ik) :: n, lda, lwork
   real(dp), allocatable :: work(:)
   real(dp) :: test(1)
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   lwork = -1
   call lapack_getri(n, amat, lda, ipiv, test, lwork, info)
   if (info == 0) then
      lwork = nint(test(1))
      allocate(work(lwork))
      call lapack_getri(n, amat, lda, ipiv, work, lwork, info)
   end if
end subroutine wrap_dgetri

end module tblite_lapack_getri
