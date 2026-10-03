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

!> @file tblite/integral/diat_trafo.f90
!> Evaluation of the diatomic scaled overlap based on projection
!> of the sigma, pi, and delta contributions onto the bond axis.
module tblite_integral_diat_trafo
   use mctc_env, only : wp
   implicit none
   private

   public :: setup_diat_trafo, diat_trafo_cache, diat_trafo

   integer, parameter :: max_diat_l = 2
   integer, parameter :: max_diat_dim = 9
   integer, parameter :: sdim(0:max_diat_l) = [1, 4, 9]

   !> Thread-private reusable cache for one atom-pair scaling
   type :: diat_trafo_cache
      !> Maximum angular momentum on atom j
      integer :: maxlj = 0
      !> Maximum angular momentum on atom i
      integer :: maxli = 0
      !> Dimension of the block on atom j
      integer :: dimj = 1
      !> Dimension of the block on atom i
      integer :: dimi = 1
      !> Angular momentum of the smaller side used for projection
      integer :: active_l = 0
      !> Apply the operator from the left (atom j) or from the right (atom i)
      logical :: left = .true.
      !> Bond axis exists, otherwise no scaling is applied
      logical :: active = .false.

      !> Unit vector along the bond axis
      real(wp) :: n(3) = 0.0_wp
      !> Derivative of the unit vector, dn(:, ic) = d n / d vec(ic)
      real(wp) :: dn(3, 3) = 0.0_wp
      !> Scaling operator
      real(wp) :: op(max_diat_dim, max_diat_dim) = 0.0_wp
      !> Derivative of the scaling operator for one Cartesian component
      real(wp) :: dop(max_diat_dim, max_diat_dim) = 0.0_wp
      !> General temporary multiplication intermediate to avoid repeated allocation
      real(wp) :: tmp(max_diat_dim, max_diat_dim) = 0.0_wp
   end type diat_trafo_cache

   real(wp), parameter :: sqrtthotw = sqrt(1.5_wp)
   real(wp), parameter :: sqrttwo = sqrt(2.0_wp)
   real(wp), parameter :: sqrtsix = sqrt(6.0_wp)

contains

!> Prepare the bond axis and its derivative once per atom pair reusing memory
pure subroutine setup_diat_trafo(cache, vec, maxlj, maxli, grad)
   !> Cache for intermediate scaling data
   type(diat_trafo_cache), intent(inout) :: cache
   !> Interatomic vector defining the bond axis
   real(wp), intent(in) :: vec(3)
   !> Maximum angular momentum on atom j
   integer, intent(in) :: maxlj
   !> Maximum angular momentum on atom i
   integer, intent(in) :: maxli
   !> Optional flag to prepare the derivative of the bond axis
   logical, intent(in), optional :: grad

   integer :: ic, jc
   real(wp) :: r2, r
   logical :: gradient

   cache%maxlj = min(max(maxlj, 0), max_diat_l)
   cache%maxli = min(max(maxli, 0), max_diat_l)
   cache%dimj = sdim(cache%maxlj)
   cache%dimi = sdim(cache%maxli)
   cache%active_l = min(cache%maxlj, cache%maxli)
   cache%left = cache%dimj <= cache%dimi

   gradient = .false.
   if (present(grad)) gradient = grad

   ! Check for vanishing bond vectors
   r2 = sum(vec**2)
   cache%active = r2 /= 0.0_wp
   cache%n = 0.0_wp
   cache%dn = 0.0_wp
   if (.not. cache%active) return

   ! Compute the unit vector along the bond axis and its derivative
   r = sqrt(r2)
   cache%n(:) = vec / r
   if (gradient) then
      do ic = 1, 3
         do jc = 1, 3
            cache%dn(jc, ic) = -cache%n(jc) * cache%n(ic) / r
         end do
         cache%dn(ic, ic) = cache%dn(ic, ic) + 1.0_wp / r
      end do
   end if
end subroutine setup_diat_trafo


!> Apply the diatomic sigma, pi and delta scaling to a shell-pair overlap block
!> and optionally calculate the Cartesian derivatives
pure subroutine diat_trafo(cache, ksig, kpi, kdel, block_overlap, block_doverlap)
   !> Cache for intermediate diatomic scaling data
   type(diat_trafo_cache), intent(inout) :: cache
   !> Scaling parameters for sigma bonding contributions
   real(wp), intent(in) :: ksig
   !> Scaling parameters for pi bonding contributions
   real(wp), intent(in) :: kpi
   !> Scaling parameters for delta bonding contributions
   real(wp), intent(in) :: kdel
   !> Diatomic block of CGTO overlap to be scaled
   real(wp), intent(inout), contiguous :: block_overlap(:, :)
   !> Derivative of diatomic block of CGTO overlap to be scaled
   real(wp), intent(inout), optional, contiguous :: block_doverlap(:, :, :)

   integer :: ic, dimj, dimi

   if (.not. cache%active) return

   dimj = cache%dimj
   dimi = cache%dimi

   ! For spherically symmetric s functions and equal scaling for both sides,
   ! the operator reduces to a constant factor.
   if ((cache%active_l < 1 .or. kpi == ksig) &
      & .and. (cache%active_l < 2 .or. kdel == ksig)) then
      if (ksig == 1.0_wp) return
      block_overlap(:dimj, :dimi) = ksig * block_overlap(:dimj, :dimi)
      if (present(block_doverlap)) then
         block_doverlap(:dimj, :dimi, :) = ksig * block_doverlap(:dimj, :dimi, :)
      end if
      return
   end if

   ! Compute block diagonal scaling operator for the current shell pair
   call get_scale_operator(cache%active_l, cache%n, ksig, kpi, kdel, cache%op)

   ! Evaluate derivatives before the overlap is overwritten.
   if (present(block_doverlap)) then
      do ic = 1, 3
         ! Derivative of block diagaonal scaling operator
         call get_scale_operator_deriv(cache%active_l, cache%n, cache%dn(:, ic), &
            & ksig, kpi, kdel, cache%dop)

         ! Scaling derivative dS' = dK S + K dS or dS' = dS K + S dK
         if (cache%left) then
            ! s block
            block_doverlap(1, :dimi, ic) = ksig * block_doverlap(1, :dimi, ic)
            ! p block
            if (cache%active_l >= 1) then
               cache%tmp(2:4, :dimi) = &
                  & matmul(cache%dop(2:4, 2:4), block_overlap(2:4, :dimi)) &
                  & + matmul(cache%op(2:4, 2:4), block_doverlap(2:4, :dimi, ic))
               block_doverlap(2:4, :dimi, ic) = cache%tmp(2:4, :dimi)
            end if
            ! d block
            if (cache%active_l >= 2) then
               cache%tmp(5:9, :dimi) = &
                  & matmul(cache%dop(5:9, 5:9), block_overlap(5:9, :dimi)) &
                  & + matmul(cache%op(5:9, 5:9), block_doverlap(5:9, :dimi, ic))
               block_doverlap(5:9, :dimi, ic) = cache%tmp(5:9, :dimi)
            end if
         else
            ! s block
            block_doverlap(:dimj, 1, ic) = ksig * block_doverlap(:dimj, 1, ic)
            ! p block
            if (cache%active_l >= 1) then
               cache%tmp(:dimj, 2:4) = &
                  & matmul(block_overlap(:dimj, 2:4), cache%dop(2:4, 2:4)) &
                  & + matmul(block_doverlap(:dimj, 2:4, ic), cache%op(2:4, 2:4))
               block_doverlap(:dimj, 2:4, ic) = cache%tmp(:dimj, 2:4)
            end if
            ! d block
            if (cache%active_l >= 2) then
               cache%tmp(:dimj, 5:9) = &
                  & matmul(block_overlap(:dimj, 5:9), cache%dop(5:9, 5:9)) &
                  & + matmul(block_doverlap(:dimj, 5:9, ic), cache%op(5:9, 5:9))
               block_doverlap(:dimj, 5:9, ic) = cache%tmp(:dimj, 5:9)
            end if
         end if
      end do
   end if

   ! Scale the overlap itself as S' = K S or S' = S K
   if (cache%left) then
      ! s block
      block_overlap(1, :dimi) = ksig * block_overlap(1, :dimi)
      ! p block
      if (cache%active_l >= 1) then
         cache%tmp(2:4, :dimi) = &
            matmul(cache%op(2:4, 2:4), block_overlap(2:4, :dimi))
         block_overlap(2:4, :dimi) = cache%tmp(2:4, :dimi)
      end if
      ! d block
      if (cache%active_l >= 2) then
         cache%tmp(5:9, :dimi) = &
            matmul(cache%op(5:9, 5:9), block_overlap(5:9, :dimi))
         block_overlap(5:9, :dimi) = cache%tmp(5:9, :dimi)
      end if
   else
      ! s block
      block_overlap(:dimj, 1) = ksig * block_overlap(:dimj, 1)
      ! p block
      if (cache%active_l >= 1) then
         cache%tmp(:dimj, 2:4) = &
            matmul(block_overlap(:dimj, 2:4), cache%op(2:4, 2:4))
         block_overlap(:dimj, 2:4) = cache%tmp(:dimj, 2:4)
      end if
      ! d block
      if (cache%active_l >= 2) then
         cache%tmp(:dimj, 5:9) = &
            matmul(block_overlap(:dimj, 5:9), cache%op(5:9, 5:9))
         block_overlap(:dimj, 5:9) = cache%tmp(:dimj, 5:9)
      end if
   end if

end subroutine diat_trafo


!> Scaling operator of the sigma, pi, and delta channels for all shells up to maxl
pure subroutine get_scale_operator(maxl, n, ksig, kpi, kdel, op)
   !> Maximum angular momentum
   integer, intent(in) :: maxl
   !> Unit vector along the bond axis
   real(wp), intent(in) :: n(3)
   !> Scaling parameters for sigma contributions
   real(wp), intent(in) :: ksig
   !> Scaling parameters for pi contributions
   real(wp), intent(in) :: kpi
   !> Scaling parameters for delta contributions
   real(wp), intent(in) :: kdel
   !> Block diagonal scaling operator in the order s, p, d
   real(wp), intent(out) :: op(:, :)

   integer :: ic
   real(wp) :: q(5), b(3, 5), btb(5, 5), wsig, wpi

   op = 0.0_wp

   ! s functions: only sigma character
   op(1, 1) = ksig

   if (maxl < 1) return

   ! p functions: q is the orbital components pointing along the bond axis.
   ! Projectors: P_sig = q q^T and P_pi = 1 - q q^T
   q(:3) = n([2, 3, 1])
   do ic = 1, 3
      op(2:4, 1+ic) = (ksig - kpi) * q(:3) * q(ic)
      op(1+ic, 1+ic) = op(1+ic, 1+ic) + kpi
   end do

   if (maxl < 2) return

   ! d functions: B are vectors coupling via the d orbitals to the bond axis direction.
   ! Projectors: P_sig = q q^T, P_pi = 2 B^T B - 4/3 P_sig, and P_del = 1 - P_sig - P_pi
   call get_d_tensor_vector(n, b)
   btb = matmul(transpose(b), b)
   q = sqrtthotw * matmul(n, b)
   wsig = (ksig - kdel) - (4.0_wp / 3.0_wp) * (kpi - kdel)
   wpi = 2.0_wp * (kpi - kdel)
   do ic = 1, 5
      op(5:9, 4+ic) = wsig * q * q(ic) + wpi * btb(:, ic)
      op(4+ic, 4+ic) = op(4+ic, 4+ic) + kdel
   end do

end subroutine get_scale_operator

!> Derivative of the scaling operator along a change of the unit vector
pure subroutine get_scale_operator_deriv(maxl, n, dn, ksig, kpi, kdel, dop)
   !> Maximum angular momentum
   integer, intent(in) :: maxl
   !> Unit vector along the bond axis
   real(wp), intent(in) :: n(3)
   !> Change of the unit vector
   real(wp), intent(in) :: dn(3)
   !> Scaling parameters for sigma contributions
   real(wp), intent(in) :: ksig
   !> Scaling parameters for pi contributions
   real(wp), intent(in) :: kpi
   !> Scaling parameters for delta contributions
   real(wp), intent(in) :: kdel
   !> Derivative of the block diagonal scaling operator
   real(wp), intent(out) :: dop(:, :)

   integer :: ic
   real(wp) :: q(5), dq(5), b(3, 5), db(3, 5), dbtb(5, 5), wsig, wpi

   dop = 0.0_wp

   if (maxl < 1) return

   q(:3) = n([2, 3, 1])
   dq(:3) = dn([2, 3, 1])
   do ic = 1, 3
      dop(2:4, 1+ic) = (ksig - kpi) * (dq(:3) * q(ic) + q(:3) * dq(ic))
   end do

   if (maxl < 2) return

   ! B is linear in n and q quadratic in n
   call get_d_tensor_vector(n, b)
   call get_d_tensor_vector(dn, db)
   dbtb = matmul(transpose(db), b) + matmul(transpose(b), db)
   q = sqrtthotw * matmul(n, b)
   dq = sqrtsix * matmul(dn, b)
   wsig = (ksig - kdel) - (4.0_wp / 3.0_wp) * (kpi - kdel)
   wpi = 2.0_wp * (kpi - kdel)
   do ic = 1, 5
      dop(5:9, 4+ic) = wsig * (dq * q(ic) + q * dq(ic)) + wpi * dbtb(:, ic)
   end do

end subroutine get_scale_operator_deriv

!> Apply the normalized traceless tensors of the d functions to a vector
pure subroutine get_d_tensor_vector(n, b)
   !> Vector the tensors are applied to
   real(wp), intent(in) :: n(3)
   !> Tensors for (xy, yz, z2, xz, x2-y2) applied to the vector
   real(wp), intent(out) :: b(3, 5)

   b(:, 1) = [n(2), n(1), 0.0_wp] / sqrttwo
   b(:, 2) = [0.0_wp, n(3), n(2)] / sqrttwo
   b(:, 3) = [-n(1), -n(2), 2.0_wp*n(3)] / sqrtsix
   b(:, 4) = [n(3), 0.0_wp, n(1)] / sqrttwo
   b(:, 5) = [n(1), -n(2), 0.0_wp] / sqrttwo

end subroutine get_d_tensor_vector

end module tblite_integral_diat_trafo
