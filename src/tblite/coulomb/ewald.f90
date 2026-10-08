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

!> @file tblite/coulomb/ewald.f90
!> Provides an utilities for implementing Ewald summations

!> Helper tools for dealing with Ewald summation related calculations
module tblite_coulomb_ewald
   use mctc_env, only : wp
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   use tblite_cutoff, only : get_lattice_points
   implicit none
   private

   public :: get_alpha, get_dir_cutoff, get_rec_cutoff, ewald_cache

   real(wp), parameter :: twopi = 2 * pi
   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

   !> Cell-dependent reciprocal coefficients, shared by energies and derivatives.
   !> These do not depend on coordinates or on the work partition.
   type :: ewald_cache
      real(wp), private :: lattice(3, 3), alpha, tolerance
      logical, private :: ko = .false.
      real(wp), allocatable :: vec(:, :), weight(:), strain(:, :, :)
      real(wp), allocatable :: weight3(:), strain3(:, :, :)
   contains
      procedure :: update => update_ewald_cache
   end type ewald_cache

   real(wp), parameter :: euler_gamma = 0.57721566490153286061_wp

   !> Evaluator for interaction term
   type, abstract :: term_type
   contains
      !> Interaction at specified distance
      procedure(get_value), deferred :: get_value
   end type term_type

   abstract interface
      !> Interaction at specified distance
      pure function get_value(self, dist, alpha, vol) result(val)
         import :: term_type, wp
         !> Instance of interaction
         class(term_type), intent(in) :: self
         !> Distance between the two atoms
         real(wp), intent(in) :: dist
         !> Parameter of the Ewald summation
         real(wp), intent(in) :: alpha
         !> Volume of the real space unit cell
         real(wp), intent(in) :: vol
         !> Value of the interaction
         real(wp) :: val
      end function get_value
   end interface

   !> Real space Coulombic interaction 1/R
   type, extends(term_type) :: dir_term
   contains
      procedure :: get_value => dir_value
   end type dir_term

   !> Reciprocal space Coulombic interaction 1/R in 3D
   type, extends(term_type) :: rec_3d_term
   contains
      procedure :: get_value => rec_3d_value
   end type rec_3d_term

   !> Real space Coulombic interaction 1/R^2
   type, extends(term_type) :: dir_mp_term
   contains
      procedure :: get_value => dir_mp_value
   end type dir_mp_term

   !> Reciprocal space Coulombic interaction 1/R^2 in 3D
   type, extends(term_type) :: rec_3d_mp_term
   contains
      procedure :: get_value => rec_3d_mp_value
   end type rec_3d_mp_term

contains


!> Rebuild only when the cell, splitting parameter, tolerance or kernel changes.
subroutine update_ewald_cache(self, lattice, alpha, tolerance, ko)
   class(ewald_cache), intent(inout) :: self
   real(wp), intent(in) :: lattice(3, 3), alpha, tolerance
   logical, intent(in), optional :: ko

   logical :: with_ko
   integer :: itr, a, b, nvec
   real(wp) :: volume, rec_lat(3, 3), g2, x, fac, weight, scale
   real(wp), allocatable :: trans(:, :)

   with_ko = .false.
   if (present(ko)) with_ko = ko
   if (allocated(self%vec)) then
      if (all(self%lattice == lattice) .and. self%alpha == alpha &
         & .and. self%tolerance == tolerance .and. (self%ko .eqv. with_ko)) return
      deallocate(self%vec, self%weight, self%strain)
      if (allocated(self%weight3)) deallocate(self%weight3, self%strain3)
   end if
   self%lattice = lattice
   self%alpha = alpha
   self%tolerance = tolerance
   self%ko = with_ko
   volume = abs(matdet_3x3(lattice))
   rec_lat = twopi*transpose(matinv_3x3(lattice))
   call get_lattice_points([.true.], rec_lat, get_rec_cutoff(alpha, volume, tolerance), trans)
   nvec = size(trans, 2) - 1
   allocate(self%vec(3, nvec), self%weight(nvec), self%strain(3, 3, nvec))
   self%vec(:, :) = trans(:, 2:)
   if (with_ko) allocate(self%weight3(nvec), self%strain3(3, 3, nvec))
   fac = 4*pi/volume
   do itr = 1, nvec
      g2 = dot_product(self%vec(:, itr), self%vec(:, itr))
      x = g2/(4*alpha*alpha)
      weight = fac*exp(-x)/g2
      scale = 2/g2 + 0.5_wp/(alpha*alpha)
      self%weight(itr) = weight
      do b = 1, 3
         do a = 1, 3
            self%strain(a, b, itr) = weight*scale*self%vec(a, itr)*self%vec(b, itr)
         end do
         self%strain(b, b, itr) = self%strain(b, b, itr) - weight
      end do
      ! Preserve the small-G exclusion of the Coulomb 1/r kernels.
      if (g2 < eps) then
         self%weight(itr) = 0.0_wp
         self%strain(:, :, itr) = 0.0_wp
      end if
      if (.not.with_ko) cycle
      self%weight3(itr) = 0.5_wp*fac*expint_e1(x)
      do b = 1, 3
         do a = 1, 3
            self%strain3(a, b, itr) = weight*self%vec(a, itr)*self%vec(b, itr)
         end do
         self%strain3(b, b, itr) = self%strain3(b, b, itr) - self%weight3(itr)
      end do
   end do
end subroutine update_ewald_cache

!> Exponential integral E1(x) for positive x
pure function expint_e1(x) result(e1)
   real(wp), intent(in) :: x
   real(wp) :: e1

   integer :: iter
   real(wp) :: term, sum, a, b, c, d, delta, h
   real(wp), parameter :: fpmin = 10.0_wp*tiny(1.0_wp)

   if (x <= 1.0_wp) then
      term = 1.0_wp
      sum = 0.0_wp
      do iter = 1, 100
         term = -term*x/real(iter, wp)
         sum = sum + term/real(iter, wp)
         if (abs(term) < epsilon(1.0_wp)*abs(sum)) exit
      end do
      e1 = -euler_gamma - log(x) - sum
   else
      b = x + 1.0_wp
      c = 1.0_wp/fpmin
      d = 1.0_wp/b
      h = d
      do iter = 1, 100
         a = -real(iter*iter, wp)
         b = b + 2.0_wp
         d = 1.0_wp/max(abs(a*d + b), fpmin)*sign(1.0_wp, a*d + b)
         c = b + a/c
         if (abs(c) < fpmin) c = fpmin
         delta = c*d
         h = h*delta
         if (abs(delta - 1.0_wp) < epsilon(1.0_wp)) exit
      end do
      e1 = exp(-x)*h
   end if
end function expint_e1


!> Convenience interface to determine Ewald splitting parameter
subroutine get_alpha(lattice, alpha, multipole)
   !> Lattice vectors
   real(wp), intent(in) :: lattice(:, :)
   !> Estimated Ewald splitting parameter
   real(wp), intent(out) :: alpha
   !> Multipole expansion is used
   logical, intent(in) :: multipole

   real(wp) :: vol, rec_lat(3, 3)
   class(term_type), allocatable :: dirv, recv

   vol = abs(matdet_3x3(lattice))
   rec_lat = twopi*transpose(matinv_3x3(lattice))
   if (multipole) then
      dirv = dir_mp_term()
      recv = rec_3d_mp_term()
   else
      dirv = dir_term()
      recv = rec_3d_term()
   end if

   call search_alpha(lattice, rec_lat, vol, eps, dirv, recv, alpha)
end subroutine get_alpha


!> Get optimal alpha-parameter for the Ewald summation by finding alpha, where
!> decline of real and reciprocal part of Ewald are equal.
subroutine search_alpha(lattice, rec_lat, volume, tolerance, dirv, recv, alpha)
   !> Lattice vectors
   real(wp), intent(in) :: lattice(:,:)
   !> Reciprocal vectors
   real(wp), intent(in) :: rec_lat(:,:)
   !> Volume of the unit cell
   real(wp), intent(in) :: volume
   !> Tolerance for difference in real and rec. part
   real(wp), intent(in) :: tolerance
   !> Real-space interaction term
   class(term_type), intent(in) :: dirv
   !> Reciprocal space interaction term
   class(term_type), intent(in) :: recv
   !> Optimal alpha
   real(wp), intent(out) :: alpha

   real(wp) :: alpl, alpr, rlen, dlen, diff
   real(wp), parameter :: alpha0 = sqrt(epsilon(0.0_wp))
   integer, parameter :: niter = 30
   integer :: ibs, stat

   rlen = sqrt(minval(sum(rec_lat(:,:)**2, dim=1)))
   dlen = sqrt(minval(sum(lattice(:,:)**2, dim=1)))

   stat = 0
   alpha = alpha0
   diff = rec_dir_diff(alpha, dirv, recv, rlen, dlen, volume)
   do while (diff < -tolerance .and. alpha <= huge(1.0_wp))
      alpha = 2.0_wp * alpha
      diff = rec_dir_diff(alpha, dirv, recv, rlen, dlen, volume)
   end do
   if (alpha > huge(1.0_wp)) then
      stat = 1
   else if (alpha == alpha0) then
      stat = 2
   end if

   if (stat == 0) then
      alpl = 0.5_wp * alpha
      do while (diff < tolerance .and. alpha <= huge(1.0_wp))
         alpha = 2.0_wp * alpha
         diff = rec_dir_diff(alpha, dirv, recv, rlen, dlen, volume)
      end do
      if (alpha > huge(1.0_wp)) then
         stat = 3
      end if
   end if

   if (stat == 0) then
      alpr = alpha
      alpha = (alpl + alpr) * 0.5_wp
      ibs = 0
      diff = rec_dir_diff(alpha, dirv, recv, rlen, dlen, volume)
      do while (abs(diff) > tolerance .and. ibs <= niter)
         if (diff < 0) then
            alpl = alpha
         else
            alpr = alpha
         end if
         alpha = (alpl + alpr) * 0.5_wp
         diff = rec_dir_diff(alpha, dirv, recv, rlen, dlen, volume)
         ibs = ibs + 1
      end do
      if (ibs > niter) then
         stat = 4
      end if
   end if

   if (stat /= 0) then
      alpha = 0.25_wp
   end if

end subroutine search_alpha


!> Return cutoff for reciprocal contributions
function get_rec_cutoff(alpha, volume, conv) result(x)
   !> Parameter of Ewald summation
   real(wp), intent(in) :: alpha
   !> Volume of the unit cell
   real(wp), intent(in) :: volume
   !> Tolerance value
   real(wp), intent(in) :: conv
   !> Magnitude of reciprocal vector
   real(wp) :: x

   class(term_type), allocatable :: term

   term = rec_3d_term()
   x = search_cutoff(term, alpha, volume, conv)

end function get_rec_cutoff


!> Return cutoff for real-space contributions
function get_dir_cutoff(alpha, conv) result(x)
   !> Parameter of Ewald summation
   real(wp), intent(in) :: alpha
   !> Tolerance for Ewald summation
   real(wp), intent(in) :: conv
   !> Magnitude of real-space vector
   real(wp) :: x

   real(wp), parameter :: volume = 0.0_wp
   class(term_type), allocatable :: term

   term = dir_term()
   x = search_cutoff(term, alpha, volume, conv)

end function get_dir_cutoff


!> Search for cutoff value of interaction term
function search_cutoff(term, alpha, volume, conv) result(x)
   !> Interaction term
   class(term_type), intent(in) :: term
   !> Parameter of Ewald summation
   real(wp), intent(in) :: alpha
   !> Volume of the unit cell
   real(wp), intent(in) :: volume
   !> Tolerance value
   real(wp), intent(in) :: conv
   !> Magnitude of reciprocal vector
   real(wp) :: x

   real(wp), parameter :: init = sqrt(epsilon(0.0_wp))
   integer, parameter :: miter = 30
   real(wp) :: xl, xr, yl, yr, y
   integer :: iter

   x = init
   y = term%get_value(x, alpha, volume)
   do while (y > conv .and. x <= huge(1.0_wp))
      x = 2.0_wp * x
      y = term%get_value(x, alpha, volume)
   end do

   xl = 0.5_wp * x
   yl = term%get_value(xl, alpha, volume)
   xr = x
   yr = y

   do iter = 1, miter
      if (yl - yr <= conv) exit
      x = 0.5_wp * (xl + xr)
      y = term%get_value(x, alpha, volume)
      if (y >= conv) then
         xl = x
         yl = y
      else
         xr = x
         yr = y
      end if
   end do

end function search_cutoff


!> Returns the difference in the decrease of the real and reciprocal parts of the
!> Ewald sum. In order to make the real space part shorter than the reciprocal
!> space part, the values are taken at different distances for the real and the
!> reciprocal space parts.
pure function rec_dir_diff(alpha, dirv, recv, rlen, dlen, volume) result(diff)
   !> Parameter for the Ewald summation
   real(wp), intent(in) :: alpha
   !> Procedure pointer to real-space routine
   class(term_type), intent(in) :: dirv
   !> Procedure pointer to reciprocal routine
   class(term_type), intent(in) :: recv
   !> Length of the shortest reciprocal space vector in the sum
   real(wp), intent(in) :: rlen
   !> Length of the shortest real space vector in the sum
   real(wp), intent(in) :: dlen
   !> Volume of the real space unit cell
   real(wp), intent(in) :: volume
   !> Difference between changes in the two terms
   real(wp) :: diff

   diff = ((recv%get_value(4*rlen, alpha, volume) - recv%get_value(5*rlen, alpha, volume))) &
      & - (dirv%get_value(2*dlen, alpha, volume) - dirv%get_value(3*dlen, alpha, volume))

end function rec_dir_diff


!> Returns the max. value of a term in the reciprocal space part of the Ewald
!> summation for a given vector length.
pure function rec_3d_value(self, dist, alpha, vol) result(rval)
   !> Instance of interaction
   class(rec_3d_term), intent(in) :: self
   !> Length of the reciprocal space vector
   real(wp), intent(in) :: dist
   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha
   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol
   !> Reciprocal term
   real(wp) :: rval

   rval = 4.0_wp*pi*(exp(-0.25_wp*dist*dist/(alpha**2))/(vol*dist*dist))

end function rec_3d_value


!> Returns the max. value of a term in the reciprocal space part of the Ewald
!> summation for a given vector length.
pure function rec_3d_mp_value(self, dist, alpha, vol) result(rval)
   !> Instance of interaction
   class(rec_3d_mp_term), intent(in) :: self
   !> Length of the reciprocal space vector
   real(wp), intent(in) :: dist
   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha
   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol
   !> Reciprocal term
   real(wp) :: rval

   real(wp) :: g2

   g2 = dist*dist

   rval = 4.0_wp*pi*(exp(-0.25_wp*g2/(alpha**2))/vol)

end function rec_3d_mp_value


!> Direct space interaction at specified distance
pure function dir_value(self, dist, alpha, vol) result(val)
   !> Instance of interaction
   class(dir_term), intent(in) :: self
   !> Distance between the two atoms
   real(wp), intent(in) :: dist
   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha
   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol
   !> Value of the interaction
   real(wp) :: val

   val = erfc(alpha*dist)/dist
end function dir_value

!> Direct space interaction at specified distance
pure function dir_mp_value(self, dist, alpha, vol) result(val)
   !> Instance of interaction
   class(dir_mp_term), intent(in) :: self
   !> Distance between the two atoms
   real(wp), intent(in) :: dist
   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha
   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol
   !> Value of the interaction
   real(wp) :: val

   real(wp) :: arg

   arg = alpha * dist

   val = (2/sqrtpi*arg*exp(-arg*arg) + erfc(arg))/dist**3
end function dir_mp_value

end module tblite_coulomb_ewald
