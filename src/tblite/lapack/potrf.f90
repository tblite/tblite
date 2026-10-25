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

!> Computes the Cholesky factorization of a real symmetric positive definite matrix A.
module tblite_lapack_potrf
   use mctc_env, only : sp, dp
   implicit none
   private

   public :: wrap_potrf, update_potrf


   !> Computes the Cholesky factorization of a real symmetric
   !> positive definite matrix A.
   !>
   !> The factorization has the form
   !>    A = U**T * U,  if UPLO = 'U', or
   !>    A = L  * L**T,  if UPLO = 'L',
   !> where U is an upper triangular matrix and L is lower triangular.
   !>
   !> This is the block version of the algorithm, calling Level 3 BLAS.
   interface wrap_potrf
      module procedure :: wrap_spotrf
      module procedure :: wrap_dpotrf
   end interface wrap_potrf

   !> Forms or reuses an upper-triangular Cholesky factorization when the
   !> packed upper triangle of A is unchanged.
   interface update_potrf
      module procedure :: update_spotrf
      module procedure :: update_dpotrf
   end interface update_potrf

   !> Computes the Cholesky factorization of a real symmetric
   !> positive definite matrix A.
   !>
   !> The factorization has the form
   !>    A = U**T * U,  if UPLO = 'U', or
   !>    A = L  * L**T,  if UPLO = 'L',
   !> where U is an upper triangular matrix and L is lower triangular.
   !>
   !> This is the block version of the algorithm, calling Level 3 BLAS.
   interface lapack_potrf
      pure subroutine spotrf(uplo, n, a, lda, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine spotrf
      pure subroutine dpotrf(uplo, n, a, lda, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dpotrf
   end interface lapack_potrf

contains

subroutine wrap_spotrf(amat, info, uplo)
   real(sp), intent(inout) :: amat(:, :)
   integer, intent(out) :: info
   character(len=1), intent(in), optional :: uplo
   integer :: n, lda
   character(len=1) :: ula

   ula = "u"
   if (present(uplo)) ula = uplo
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   call lapack_potrf(ula, n, amat, lda, info)
end subroutine wrap_spotrf

subroutine wrap_dpotrf(amat, info, uplo)
   real(dp), intent(inout) :: amat(:, :)
   integer, intent(out) :: info
   character(len=1), intent(in), optional :: uplo
   integer :: n, lda
   character(len=1) :: ula

   ula = "u"
   if (present(uplo)) ula = uplo
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   call lapack_potrf(ula, n, amat, lda, info)
end subroutine wrap_dpotrf

subroutine update_spotrf(amat, packed, chol, info)
   !> Symmetric positive-definite matrix (upper triangle read)
   real(sp), intent(in) :: amat(:, :)
   !> Packed upper-triangle snapshot of the valid cached factorization
   real(sp), allocatable, intent(inout) :: packed(:)
   !> Cached upper-triangular Cholesky factor
   real(sp), allocatable, intent(inout) :: chol(:, :)
   !> Status flag from spotrf (0 when cached or successfully factored)
   integer, intent(out) :: info

   integer :: n, npk, ii, jj, k0
   logical :: same

   info = 0
   n = size(amat, 1)
   npk = n * (n + 1) / 2
   same = allocated(packed) .and. allocated(chol)
   if (same) same = (size(packed) == npk .and. size(chol, 1) == n)
   if (same) then
      outer: do jj = 1, n
         k0 = jj * (jj - 1) / 2
         do ii = 1, jj
            if (packed(k0 + ii) /= amat(ii, jj)) then
               same = .false.
               exit outer
            end if
         end do
      end do outer
   end if
   if (.not. same) then
      if (allocated(packed)) then
         if (size(packed) /= npk) deallocate(packed)
      end if
      if (allocated(chol)) then
         if (size(chol, 1) /= n) deallocate(chol)
      end if
      if (.not. allocated(packed)) allocate(packed(npk))
      if (.not. allocated(chol)) allocate(chol(n, n))
      do jj = 1, n
         k0 = jj * (jj - 1) / 2
         packed(k0 + 1:k0 + jj) = amat(1:jj, jj)
         chol(1:jj, jj) = amat(1:jj, jj)
      end do
      call wrap_spotrf(chol, info, uplo="U")
      if (info /= 0) deallocate(packed)
   end if
end subroutine update_spotrf

subroutine update_dpotrf(amat, packed, chol, info)
   !> Symmetric positive-definite matrix (upper triangle read)
   real(dp), intent(in) :: amat(:, :)
   !> Packed upper-triangle snapshot of the valid cached factorization
   real(dp), allocatable, intent(inout) :: packed(:)
   !> Cached upper-triangular Cholesky factor
   real(dp), allocatable, intent(inout) :: chol(:, :)
   !> Status flag from dpotrf (0 when cached or successfully factored)
   integer, intent(out) :: info

   integer :: n, npk, ii, jj, k0
   logical :: same

   info = 0
   n = size(amat, 1)
   npk = n * (n + 1) / 2
   same = allocated(packed) .and. allocated(chol)
   if (same) same = (size(packed) == npk .and. size(chol, 1) == n)
   if (same) then
      outer: do jj = 1, n
         k0 = jj * (jj - 1) / 2
         do ii = 1, jj
            if (packed(k0 + ii) /= amat(ii, jj)) then
               same = .false.
               exit outer
            end if
         end do
      end do outer
   end if
   if (.not. same) then
      if (allocated(packed)) then
         if (size(packed) /= npk) deallocate(packed)
      end if
      if (allocated(chol)) then
         if (size(chol, 1) /= n) deallocate(chol)
      end if
      if (.not. allocated(packed)) allocate(packed(npk))
      if (.not. allocated(chol)) allocate(chol(n, n))
      do jj = 1, n
         k0 = jj * (jj - 1) / 2
         packed(k0 + 1:k0 + jj) = amat(1:jj, jj)
         chol(1:jj, jj) = amat(1:jj, jj)
      end do
      call wrap_dpotrf(chol, info, uplo="U")
      if (info /= 0) deallocate(packed)
   end if
end subroutine update_dpotrf

end module tblite_lapack_potrf
