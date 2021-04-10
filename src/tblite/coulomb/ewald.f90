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

!> Helper tools for dealing with Ewald summation related calculations
module tblite_coulomb_ewald
   use mctc_env, only : wp
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   implicit none
   private

   public :: get_alpha

   real(wp), parameter :: twopi = 2 * pi
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

   abstract interface
      !> Returns the max. value of a term in the reciprocal space part of the Ewald
      !> summation for a given vector length.
      pure function get_rec_term_gen(gg, alpha, vol) result(gTerm)
         import :: wp

         !> Length of the reciprocal space vector
         real(wp), intent(in) :: gg

         !> Parameter of the Ewald summation
         real(wp), intent(in) :: alpha

         !> Volume of the real space unit cell
         real(wp), intent(in) :: vol

         !> Reciprocal term
         real(wp) :: gTerm

      end function get_rec_term_gen
   end interface

contains


subroutine get_alpha(lattice, alpha)
   real(wp), intent(in) :: lattice(:, :)
   real(wp), intent(out) :: alpha
   real(wp) :: vol, rec_lat(3, 3)

   vol = abs(matdet_3x3(lattice))
   rec_lat = twopi*transpose(matinv_3x3(lattice))

   call search_alpha(lattice, rec_lat, vol, eps, alpha)

end subroutine get_alpha


!> Get optimal alpha-parameter for the Ewald summation by finding alpha, where
!> decline of real and reciprocal part of Ewald are equal.
subroutine search_alpha(lattice, rec_lat, volume, tolerance, alpha)
   !> Lattice vectors
   real(wp), intent(in) :: lattice(:,:)
   !> Reciprocal vectors
   real(wp), intent(in) :: rec_lat(:,:)
   !> Volume of the unit cell
   real(wp), intent(in) :: volume
   !> Tolerance for difference in real and rec. part
   real(wp), intent(in) :: tolerance
   !> Optimal alpha
   real(wp), intent(out) :: alpha

   real(wp) :: alpl, alpr, rlen, dlen, diff
   real(wp), parameter :: alpha0 = 1.0e-8_wp
   integer, parameter :: niter = 30
   integer :: ibs, stat

   rlen = sqrt(minval(sum(rec_lat(:,:)**2, dim=1)))
   dlen = sqrt(minval(sum(lattice(:,:)**2, dim=1)))

   stat = 0
   alpha = alpha0
   diff = rec_dir_diff(alpha, get_rec_term_3d, rlen, dlen, volume)
   do while (diff < -tolerance .and. alpha <= huge(1.0_wp))
      alpha = 2.0_wp * alpha
      diff = rec_dir_diff(alpha, get_rec_term_3d, rlen, dlen, volume)
   end do
   if (alpha > huge(1.0_wp)) then
      stat = 1
   elseif (alpha == alpha0) then
      stat = 2
   end if

   if (stat == 0) then
      alpl = 0.5_wp * alpha
      do while (diff < tolerance .and. alpha <= huge(1.0_wp))
         alpha = 2.0_wp * alpha
         diff = rec_dir_diff(alpha, get_rec_term_3d, rlen, dlen, volume)
      end do
      if (alpha > huge(1.0_wp)) then
         stat = 3
      end if
   end if

   if (stat == 0) then
      alpr = alpha
      alpha = (alpl + alpr) * 0.5_wp
      ibs = 0
      diff = rec_dir_diff(alpha, get_rec_term_3d, rlen, dlen, volume)
      do while (abs(diff) > tolerance .and. ibs <= niter)
         if (diff < 0) then
            alpl = alpha
         else
            alpr = alpha
         end if
         alpha = (alpl + alpr) * 0.5_wp
         diff = rec_dir_diff(alpha, get_rec_term_3d, rlen, dlen, volume)
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

!> Returns the difference in the decrease of the real and reciprocal parts of the
!> Ewald sum. In order to make the real space part shorter than the reciprocal
!> space part, the values are taken at different distances for the real and the
!> reciprocal space parts.
pure function rec_dir_diff(alpha, get_rec_term, rlen, dlen, volume) result(diff)

   !> Parameter for the Ewald summation
   real(wp), intent(in) :: alpha

   !> Procedure pointer to reciprocal routine
   procedure(get_rec_term_gen) :: get_rec_term

   !> Length of the shortest reciprocal space vector in the sum
   real(wp), intent(in) :: rlen

   !> Length of the shortest real space vector in the sum
   real(wp), intent(in) :: dlen

   !> Volume of the real space unit cell
   real(wp), intent(in) :: volume

   !> Difference between changes in the two terms
   real(wp) :: diff

   diff = ((get_rec_term(4*rlen, alpha, volume) &
      & - get_rec_term(5*rlen, alpha, volume))) &
      & - (get_dir_term(2*dlen, alpha) - get_dir_term(3*dlen, alpha))

end function rec_dir_diff

!> Returns the max. value of a term in the real space part of the Ewald summation
!> for a given vector length.
pure function get_dir_term(rr, alpha) result(dval)

   !> Length of the real space vector
   real(wp), intent(in) :: rr

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Real space term
   real(wp) :: dval

   dval = erfc(alpha*rr)/rr

end function get_dir_term


!> Returns the max. value of a term in the reciprocal space part of the Ewald
!> summation for a given vector length.
pure function get_rec_term_3d(gg, alpha, vol) result(rval)

   !> Length of the reciprocal space vector
   real(wp), intent(in) :: gg

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol

   !> Reciprocal term
   real(wp) :: rval

   rval = 4.0_wp*pi*(exp(-0.25_wp*gg*gg/(alpha**2))/(vol*gg*gg))

end function get_rec_term_3d

end module tblite_coulomb_ewald
