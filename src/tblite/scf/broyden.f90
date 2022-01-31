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

module tblite_scf_broyden
   use mctc_env, only : wp
   use tblite_lapack, only : getrf, getrs
   implicit none
   private

   public :: broyden_mixer, new_broyden, broyden

   type :: broyden_mixer
      integer :: ndim
      integer :: memory
      integer :: iter
      integer :: iset
      integer :: idif
      integer :: iget
      real(wp) :: damp
      real(wp), allocatable :: df(:, :)
      real(wp), allocatable :: u(:, :)
      real(wp), allocatable :: a(:, :)
      real(wp), allocatable :: dq(:)
      real(wp), allocatable :: dqlast(:)
      real(wp), allocatable :: qlast_in(:)
      real(wp), allocatable :: omega(:)
      real(wp), allocatable :: q_in(:)
   contains
      procedure :: next
      generic :: set => set_1d, set_2d, set_3d
      procedure :: set_1d
      procedure :: set_2d
      procedure :: set_3d
      generic :: diff => diff_1d, diff_2d, diff_3d
      procedure :: diff_1d
      procedure :: diff_2d
      procedure :: diff_3d
      generic :: get => get_1d, get_2d, get_3d
      procedure :: get_1d
      procedure :: get_2d
      procedure :: get_3d
      procedure :: get_error
   end type broyden_mixer

contains

subroutine new_broyden(self, memory, ndim, damp)
   type(broyden_mixer), intent(out) :: self
   integer, intent(in) :: memory
   integer, intent(in) :: ndim
   real(wp), intent(in) :: damp

   self%ndim = ndim
   self%memory = memory
   self%iter = 0
   self%iset = 0
   self%idif = 0
   self%iget = 0
   self%damp = damp
   allocate(self%df(ndim, memory))
   allocate(self%u(ndim, memory))
   allocate(self%a(memory, memory))
   allocate(self%dq(ndim))
   allocate(self%dqlast(ndim))
   allocate(self%qlast_in(ndim))
   allocate(self%omega(memory))
   allocate(self%q_in(ndim))
end subroutine new_broyden

subroutine set_2d(self, qvec)
   class(broyden_mixer), intent(inout) :: self
   real(wp), contiguous, intent(in), target :: qvec(:, :)
   real(wp), pointer :: qptr(:)
   qptr(1:size(qvec)) => qvec
   call self%set(qptr)
end subroutine set_2d

subroutine set_3d(self, qvec)
   class(broyden_mixer), intent(inout) :: self
   real(wp), contiguous, intent(in), target :: qvec(:, :, :)
   real(wp), pointer :: qptr(:)
   qptr(1:size(qvec)) => qvec
   call self%set(qptr)
end subroutine set_3d

subroutine set_1d(self, qvec)
   class(broyden_mixer), intent(inout) :: self
   real(wp), intent(in) :: qvec(:)
   self%q_in(self%iset+1:self%iset+size(qvec)) = qvec
   self%iset = self%iset + size(qvec)
end subroutine set_1d

subroutine diff_2d(self, qvec)
   class(broyden_mixer), intent(inout) :: self
   real(wp), contiguous, intent(in), target :: qvec(:, :)
   real(wp), pointer :: qptr(:)
   qptr(1:size(qvec)) => qvec
   call self%diff(qptr)
end subroutine diff_2d

subroutine diff_3d(self, qvec)
   class(broyden_mixer), intent(inout) :: self
   real(wp), contiguous, intent(in), target :: qvec(:, :, :)
   real(wp), pointer :: qptr(:)
   qptr(1:size(qvec)) => qvec
   call self%diff(qptr)
end subroutine diff_3d

subroutine diff_1d(self, qvec)
   class(broyden_mixer), intent(inout) :: self
   real(wp), intent(in) :: qvec(:)
   self%dq(self%idif+1:self%idif+size(qvec)) = qvec &
      & - self%q_in(self%idif+1:self%idif+size(qvec))
   self%idif = self%idif + size(qvec)
end subroutine diff_1d

subroutine next(self)
   class(broyden_mixer), intent(inout) :: self
   self%iset = 0
   self%idif = 0
   self%iget = 0
   self%iter = self%iter + 1
   call broyden(self%ndim, self%q_in, self%qlast_in, self%dq, self%dqlast, &
      & self%iter, self%memory, self%damp, self%omega, self%df, self%u, self%a)
end subroutine next

subroutine get_2d(self, qvec)
   class(broyden_mixer), intent(inout) :: self
   real(wp), contiguous, intent(out), target :: qvec(:, :)
   real(wp), pointer :: qptr(:)
   qptr(1:size(qvec)) => qvec
   call self%get(qptr)
end subroutine get_2d

subroutine get_3d(self, qvec)
   class(broyden_mixer), intent(inout) :: self
   real(wp), contiguous, intent(out), target :: qvec(:, :, :)
   real(wp), pointer :: qptr(:)
   qptr(1:size(qvec)) => qvec
   call self%get(qptr)
end subroutine get_3d

subroutine get_1d(self, qvec)
   class(broyden_mixer), intent(inout) :: self
   real(wp), intent(out) :: qvec(:)
   qvec(:) = self%q_in(self%iget+1:self%iget+size(qvec))
   self%iget = self%iget + size(qvec)
end subroutine get_1d

subroutine broyden(n, q, qlast, dq, dqlast, iter, memory, alpha, omega, df, u, a)
   integer, intent(in) :: n
   integer, intent(in) :: iter
   integer, intent(in) :: memory
   real(wp), intent(inout) :: q(n)
   real(wp), intent(inout) :: qlast(n)
   real(wp), intent(in) :: dq(n)
   real(wp), intent(inout) :: dqlast(n)
   real(wp), intent(inout) :: df(n, memory)
   real(wp), intent(inout) :: u(n, memory)
   real(wp), intent(inout) :: a(memory, memory)
   real(wp), intent(inout) :: omega(memory)
   real(wp), intent(in) :: alpha

   real(wp), allocatable :: beta(:,:), c(:, :)
   integer :: i, it1
   real(wp) :: inv, omega0, minw, maxw, wfac

   it1 = iter - 1

   ! set parameters
   ! alpha = 0.25d0
   omega0 = 0.01d0
   minw = 1.0d0
   maxw = 100000.0d0
   wfac = 0.01d0
   ! wfac = 0.05d0

   ! if case for first iteration: simple damping
   if (iter == 1) then
      dqlast(:) = dq
      qlast(:) = q
      q(:) = q + alpha * dq
      return
   end if

   allocate(beta(it1,it1),c(it1, 1))

   ! create omega (weight) for the current iteration
   omega(it1) = sqrt(dot_product(dq, dq))
   if (omega(it1) > (wfac / maxw)) then
      omega(it1) = wfac / omega(it1)
   else
      omega(it1) = maxw
   end if
   if (omega(it1) < minw) then
      omega(it1) = minw
   end if

   ! Build dF(iter-1)
   df(:, it1) = dq - dqlast
   inv = max(sqrt(dot_product(df(:, it1), df(:, it1))), epsilon(1.0_wp))
   inv = 1.0_wp / inv
   df(:, it1) = inv*df(:, it1)

   ! Next: build a, beta, c, gamma
   do i = 1, it1
      a(i, it1) = dot_product(df(:, i), df(:, it1))
      a(it1, i) = a(i, it1)
      c(i, 1) = omega(i) * dot_product(df(:, i), dq)
   end do

   ! Build beta from a and omega
   do i = 1, it1
      beta(:it1, i) = omega(:it1) * omega(i) * a(:it1, i)
      beta(i, i) = beta(i, i) + omega0*omega0
   end do

   ! build beta^-1
   call lineq(beta, c)

   ! Build |u>
   u(:, it1) = alpha * df(:, it1) + inv * (q-qlast) !!!

   ! save charges and deltas
   dqlast(:) = dq
   qlast(:) = q

   ! calculate new charges
   q(:) = q + alpha * dq

   do i = 1, it1
      q(:) = q - omega(i) * c(i, 1) * u(:, i)
   end do

end subroutine broyden

subroutine lineq(a, c)
   real(wp), intent(inout) :: a(:, :)
   real(wp), intent(inout) :: c(:, :)

   integer info
   integer, allocatable :: ipiv(:)

   allocate(ipiv(size(a, 1)))
   ! LU decomoposition of a general matrix
   call getrf(a, ipiv, info)
   if (info == 0) then
      ! generate inverse of a matrix given its LU decomposition
      call getrs(a, c, ipiv, info, trans="t")
   endif
   if (info /= 0)then
      error stop "Error in Broyden matrix inversion!"
   endif

end subroutine lineq

pure function get_error(self) result(error)
   class(broyden_mixer), intent(in) :: self
   real(wp) :: error
   integer :: i
   error = 0.0_wp
   do i = 1, size(self%dq)
      error = error + self%dq(i)**2 / size(self%dq)
   end do
   error = sqrt(error)
end function get_error

end module tblite_scf_broyden
