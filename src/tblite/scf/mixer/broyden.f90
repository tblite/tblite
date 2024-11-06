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

!> @file tblite/scf/mixer/broyden.f90
!> Provides an electronic mixer implementation

!> Implementation of a modified Broyden mixing
module tblite_scf_mixer_broyden
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_lapack, only : getrf, getrs
   use tblite_scf_mixer_type, only : mixer_type
   implicit none
   private

   public :: new_broyden


   !> Configuration for the Broyden mixer
   type, public :: broyden_input
      !> Number of steps to keep in memory
      integer :: memory
      !> Damping parameter
      real(wp) :: damp
   end type broyden_input

   !> Electronic mixer using modified Broyden scheme
   type, public, extends(mixer_type) :: broyden_mixer
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
      !> Apply mixing to the density
      procedure :: next
      !> Set new density from 1D array
      procedure :: set_1d
      !> Set difference between new and old density from 1D array
      procedure :: diff_1d
      !> Get density as 1D array
      procedure :: get_1d
      !> Get error metric from mixing
      procedure :: get_error
   end type broyden_mixer

contains

!> Create new instance of electronic mixer
subroutine new_broyden(self, ndim, input)
   !> Instance of the mixer
   type(broyden_mixer), intent(out) :: self
   !> Number of variables to consider
   integer, intent(in) :: ndim
   !> Configuration of the Broyden mixer
   type(broyden_input), intent(in) :: input

   self%ndim = ndim
   self%memory = input%memory
   self%iter = 0
   self%iset = 0
   self%idif = 0
   self%iget = 0
   self%damp = input%damp
   allocate(self%df(ndim, input%memory))
   allocate(self%u(ndim, input%memory))
   allocate(self%a(input%memory, input%memory))
   allocate(self%dq(ndim))
   allocate(self%dqlast(ndim))
   allocate(self%qlast_in(ndim))
   allocate(self%omega(input%memory))
   allocate(self%q_in(ndim))
end subroutine new_broyden

!> Set new density from 1D array
subroutine set_1d(self, qvec)
   !> Instance of the mixer
   class(broyden_mixer), intent(inout) :: self
   !> Density vector
   real(wp), intent(in) :: qvec(:)
   self%q_in(self%iset+1:self%iset+size(qvec)) = qvec
   self%iset = self%iset + size(qvec)
end subroutine set_1d

!> Set difference between new and old density from 1D array
subroutine diff_1d(self, qvec)
   !> Instance of the mixer
   class(broyden_mixer), intent(inout) :: self
   !> Density vector
   real(wp), intent(in) :: qvec(:)
   self%dq(self%idif+1:self%idif+size(qvec)) = qvec &
      & - self%q_in(self%idif+1:self%idif+size(qvec))
   self%idif = self%idif + size(qvec)
end subroutine diff_1d

!> Apply mixing to the density
subroutine next(self, error)
   !> Instance of the mixer
   class(broyden_mixer), intent(inout) :: self
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: info

   self%iset = 0
   self%idif = 0
   self%iget = 0
   self%iter = self%iter + 1
   call broyden(self%ndim, self%q_in, self%qlast_in, self%dq, self%dqlast, &
      & self%iter, self%memory, self%damp, self%omega, self%df, self%u, self%a, info)
   if (info /= 0) then
      call fatal_error(error, "Broyden mixing failed to obtain next iteration")
   end if
end subroutine next

!> Get density as 1D array
subroutine get_1d(self, qvec)
   !> Instance of the mixer
   class(broyden_mixer), intent(inout) :: self
   !> Density vector
   real(wp), intent(out) :: qvec(:)
   qvec(:) = self%q_in(self%iget+1:self%iget+size(qvec))
   self%iget = self%iget + size(qvec)
end subroutine get_1d

subroutine broyden(n, q, qlast, dq, dqlast, iter, memory, alpha, omega, df, u, a, info)
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
   integer, intent(out) :: info

   real(wp), allocatable :: beta(:,:), c(:, :)
   integer :: i, j, it1, itn
   real(wp) :: inv, omega0, minw, maxw, wfac

   info = 0
   itn = iter - 1
   it1 = mod(itn - 1, memory) + 1

   ! set parameters
   ! alpha = 0.25e0_wp
   omega0 = 0.01e0_wp
   minw = 1.0e0_wp
   maxw = 100000.0e0_wp
   wfac = 0.01e0_wp
   ! wfac = 0.05e0_wp

   ! if case for first iteration: simple damping
   if (iter == 1) then
      dqlast(:) = dq
      qlast(:) = q
      q(:) = q + alpha * dq
      return
   end if

   allocate(beta(min(memory, itn), min(memory, itn)), c(min(memory, itn), 1))

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
   do j = max(1, itn - memory + 1), itn
      i = mod(j - 1, memory) + 1
      a(i, it1) = dot_product(df(:, i), df(:, it1))
      a(it1, i) = a(i, it1)
      c(i, 1) = omega(i) * dot_product(df(:, i), dq)
   end do

   ! Build beta from a and omega
   do j = max(1, itn - memory + 1), itn
      i = mod(j - 1, memory) + 1
      beta(:it1, i) = omega(:it1) * omega(i) * a(:it1, i)
      beta(i, i) = beta(i, i) + omega0*omega0
   end do

   ! build beta^-1
   call lineq(beta, c, info)
   if (info /= 0) return

   ! Build |u>
   u(:, it1) = alpha * df(:, it1) + inv * (q-qlast) !!!

   ! save charges and deltas
   dqlast(:) = dq
   qlast(:) = q

   ! calculate new charges
   q(:) = q + alpha * dq

   do j = max(1, itn - memory + 1), itn
      i = mod(j - 1, memory) + 1
      q(:) = q - omega(i) * c(i, 1) * u(:, i)
   end do

end subroutine broyden

subroutine lineq(a, c, info)
   real(wp), intent(inout) :: a(:, :)
   real(wp), intent(inout) :: c(:, :)
   integer, intent(out) :: info

   integer, allocatable :: ipiv(:)

   allocate(ipiv(size(a, 1)))
   ! LU decomoposition of a general matrix
   call getrf(a, ipiv, info)
   if (info == 0) then
      ! generate inverse of a matrix given its LU decomposition
      call getrs(a, c, ipiv, info, trans="t")
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

end module tblite_scf_mixer_broyden
