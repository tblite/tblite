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
!> Evaluation of the diatomic scaled overlap
module tblite_integral_diat_trafo
   use mctc_env, only : wp
   use tblite_blas, only : gemm
   implicit none
   private

   public :: setup_diat_trafo, diat_trafo_cache, diat_trafo

   integer, parameter :: max_diat_l = 2
   integer, parameter :: max_diat_dim = 9
   integer, parameter :: sdim(0:max_diat_l) = [1, 4, 9]

   !> Thread-private reusable cache for one atom-pair transformation
   type :: diat_trafo_cache
      !> Maximum angular momentum on atom j
      integer :: maxlj = 0
      !> Maximum angular momentum on atom i
      integer :: maxli = 0
      !> Dimension of the block on atom j
      integer :: dimj = 1
      !> Dimension of the block on atom i
      integer :: dimi = 1

      !> Ordinary rotation matrix used away from the z-axis singularity
      !> and for the x and z component at the z-axis singularity.
      real(wp) :: rot(max_diat_dim, max_diat_dim) = 0.0_wp
      !> Alternative y-direction limiting rotation at vec || z.
      real(wp) :: rot_ylim(max_diat_dim, max_diat_dim) = 0.0_wp
      !> Derivatives of the rotation matrix.
      real(wp) :: drot(max_diat_dim, max_diat_dim, 3) = 0.0_wp

      !> General temporary multiplication intermediate
      real(wp) :: tmpprod(max_diat_dim, max_diat_dim) = 0.0_wp
      !> Overlap matrix tranformed to the diatomic frame with optional scaling
      real(wp) :: diat_overlap(max_diat_dim, max_diat_dim) = 0.0_wp
      !> Temporary storage for derivatives in the diatomic frame
      real(wp) :: diat_deriv(max_diat_dim, max_diat_dim) = 0.0_wp
   end type diat_trafo_cache

contains

!> Prepare the diatomic frame rotation matrices once per atom pair reusing memory.
pure subroutine setup_diat_trafo(cache, vec, maxlj, maxli, grad)
   !> Cache for intermediate diatomic frame transformation data
   type(diat_trafo_cache), intent(inout) :: cache
   !> Interatomic vector defining the diatomic frame z-axis
   real(wp), intent(in) :: vec(3)
   !> Maximum angular momentum on atom j
   integer, intent(in) :: maxlj
   !> Maximum angular momentum on atom i
   integer, intent(in) :: maxli
   !> Optional flag to prepare diatomic frame transformation gradient matrices
   logical, intent(in), optional :: grad

   integer :: maxl
   logical :: gradient

   cache%maxlj = min(max(maxlj, 0), max_diat_l)
   cache%maxli = min(max(maxli, 0), max_diat_l)
   cache%dimj = sdim(cache%maxlj)
   cache%dimi = sdim(cache%maxli)
   maxl = max(cache%maxlj, cache%maxli)

   if (present(grad)) then
      gradient = grad
   else
      gradient = .false.
   end if

   ! Setup transformation matrix with optional gradient
   if (gradient) then
      call d_harmtr(maxl, vec, cache%rot, cache%rot_ylim, cache%drot)
   else
      call harmtr(maxl, vec, cache%rot)
      cache%rot_ylim = cache%rot
   end if
end subroutine setup_diat_trafo

!> Apply a prepared diatomic frame transformation to a shell-pair overlap block
!> and optionally calaculate the Cartesian derivatives
pure subroutine diat_trafo(cache, ksig, kpi, kdel, block_overlap, block_doverlap)
   !> Cache for intermediate diatomic frame transformation data
   type(diat_trafo_cache), intent(inout) :: cache
   !> Scaling parameters for sigma bonding contributions
   real(wp), intent(in) :: ksig
   !> Scaling parameters for pi bonding contributions
   real(wp), intent(in) :: kpi
   !> Scaling parameters for delta bonding contributions
   real(wp), intent(in) :: kdel
   !> Diatomic block of CGTO overlap to be transformed and scaled
   real(wp), intent(inout), contiguous :: block_overlap(:, :)
   !> Derivative of diatomic block of CGTO overlap to be transformed and scaled
   real(wp), intent(inout), optional, contiguous :: block_doverlap(:, :, :)

   ! Transform only for a non-trivial scaling
   if (ksig == 1.0_wp .and. kpi == 1.0_wp .and. kdel == 1.0_wp) return

   ! Scale s-s block without transformation
   if (cache%dimj == 1 .and. cache%dimi == 1) then
      block_overlap(1, 1) = ksig * block_overlap(1, 1)
      if (present(block_doverlap)) then
         block_doverlap(1, 1, :) = ksig * block_doverlap(1, 1, :)
      end if
      return
   end if

   ! Transform the overlap submatrix to the diatomic frame: S' = O_j^T * S * O_i
   call gemm(amat=cache%rot, bmat=block_overlap, cmat=cache%tmpprod, &
      & transa="T", m=cache%dimj, n=cache%dimi, k=cache%dimj)
   call gemm(amat=cache%tmpprod, bmat=cache%rot, cmat=cache%diat_overlap, &
      & m=cache%dimj, n=cache%dimi, k=cache%dimi)

   ! Scale diatomic frame overlap: Ssc' = K(S')
   call scale_diatomic_frame(cache%diat_overlap, ksig, kpi, kdel, &
      & cache%maxlj, cache%maxli)

   if (present(block_doverlap)) then
      ! Gradient x-component of diatomic frame scaled overlap
      call diat_trafo_grad_component(block_overlap, block_doverlap(:, :, 1), &
         & cache%diat_overlap, cache%rot, cache%drot(:, :, 1), cache%tmpprod, &
         & cache%diat_deriv, cache%dimj, cache%dimi, cache%maxlj, cache%maxli, &
         & ksig, kpi, kdel)
      ! Gradient y-component of diatomic frame scaled overlap
      call diat_trafo_grad_component(block_overlap, block_doverlap(:, :, 2), &
         & cache%diat_overlap, cache%rot_ylim, cache%drot(:, :, 2), cache%tmpprod, &
         & cache%diat_deriv, cache%dimj, cache%dimi, cache%maxlj, cache%maxli, &
         & ksig, kpi, kdel)
      ! Gradient z-component of diatomic frame scaled overlap
      call diat_trafo_grad_component(block_overlap, block_doverlap(:, :, 3), &
         & cache%diat_overlap, cache%rot, cache%drot(:, :, 3), cache%tmpprod, &
         & cache%diat_deriv, cache%dimj, cache%dimi, cache%maxlj, cache%maxli, &
         & ksig, kpi, kdel)
   end if

   ! Transform the diatomic-frame-scaled overlap back to the original frame:
   ! Ssc = O_j Ssc' O_i^T
   call gemm(amat=cache%rot, bmat=cache%diat_overlap, cmat=cache%tmpprod, &
      & m=cache%dimj, n=cache%dimi, k=cache%dimj)
   call gemm(amat=cache%tmpprod, bmat=cache%rot, cmat=block_overlap, &
      & transb="T", m=cache%dimj, n=cache%dimi, k=cache%dimi)

end subroutine diat_trafo

! Calculate diatomic frame scaled overlap gradient for one cartesian component
pure subroutine diat_trafo_grad_component(block_overlap, block_doverlap, &
   & scaled_diat_overlap, rotation, drotation, tmpprod, diat_deriv, &
   & dimj, dimi, maxlj, maxli, ksig, kpi, kdel)
   !> Diatomic block of CGTO overlap to be transformed and scaled
   real(wp), intent(in), contiguous :: block_overlap(:, :)
   !> One cartesian component of the derivative of the diatomic block of CGTO overlap
   real(wp), intent(inout), contiguous :: block_doverlap(:, :)
   !> Scaled overlap in the diatomic frame
   real(wp), intent(in), contiguous :: scaled_diat_overlap(:, :)
   !> Rotation matrix for the current cartesian component
   real(wp), intent(in), contiguous :: rotation(:, :)
   !> Derivative of the rotation matrix for the current cartesian component
   real(wp), intent(in), contiguous :: drotation(:, :)
   !> Temporary storage for intermediate products, reused across components
   real(wp), intent(inout), contiguous :: tmpprod(:, :)
   !> Temporary storage for derivatives in the diatomic frame
   real(wp), intent(inout), contiguous :: diat_deriv(:, :)
   !> Dimension of the block on atom j
   integer, intent(in) :: dimj
   !> Dimension of the block on atom i
   integer, intent(in) :: dimi
   !> Maximum angular momentum on atom j
   integer, intent(in) :: maxlj
   !> Maximum angular momentum on atom i
   integer, intent(in) :: maxli
   !> Scaling parameters for sigma bonding contributions
   real(wp), intent(in) :: ksig
   !> Scaling parameters for pi bonding contributions
   real(wp), intent(in) :: kpi
   !> Scaling parameters for delta bonding contributions
   real(wp), intent(in) :: kdel

   ! diat_deriv = dO_j^T S O_i
   call gemm(amat=block_overlap, bmat=rotation, cmat=tmpprod, &
      & m=dimj, n=dimi, k=dimi)
   call gemm(amat=drotation, bmat=tmpprod, cmat=diat_deriv, &
      & transa="T", m=dimj, n=dimi, k=dimj)

   ! diat_deriv += O_j^T (dS O_i + S dO_i)
   call gemm(amat=block_doverlap, bmat=rotation, cmat=tmpprod, &
      & m=dimj, n=dimi, k=dimi)
   call gemm(amat=block_overlap, bmat=drotation, cmat=tmpprod, &
      & beta=1.0_wp, m=dimj, n=dimi, k=dimi)
   call gemm(amat=rotation, bmat=tmpprod, cmat=diat_deriv, &
      & transa="T", beta=1.0_wp, m=dimj, n=dimi, k=dimj)

   ! Scale derivative in the diatomic frame
   call scale_diatomic_frame(diat_deriv, ksig, kpi, kdel, maxlj, maxli)

   ! dSsc = [dO_j (Ssc') + O_j (dSsc')] O_i^T
   call gemm(amat=drotation, bmat=scaled_diat_overlap, cmat=tmpprod, &
      & m=dimj, n=dimi, k=dimj)
   call gemm(amat=rotation, bmat=diat_deriv, cmat=tmpprod, &
      & beta=1.0_wp, m=dimj, n=dimi, k=dimj)
   call gemm(amat=tmpprod, bmat=rotation, cmat=block_doverlap, &
      & transb="T", m=dimj, n=dimi, k=dimi)

   ! dSsc += O_j (Ssc') dO_i^T
   call gemm(amat=rotation, bmat=scaled_diat_overlap, cmat=tmpprod, &
      & m=dimj, n=dimi, k=dimj)
   call gemm(amat=tmpprod, bmat=drotation, cmat=block_doverlap, &
      & transb="T", beta=1.0_wp, m=dimj, n=dimi, k=dimi)

end subroutine diat_trafo_grad_component

!> Setup diatomic frame transformation matrix for up to d-functions.
pure subroutine harmtr(maxl, vec, trafomat)
   !> Maximum angular momentum
   integer, intent(in)  :: maxl
   !> Interatomic vector defining the diatomic frame z-axis
   real(wp), intent(in) :: vec(3)
   !> Transformation matrix
   real(wp), intent(out) :: trafomat(:, :)

   real(wp) :: cos2p, cos2t, cosp, cost, sin2p, sin2t, sinp, sint, sqrt3, len
   real(wp) :: norm_vec(3)

   trafomat = 0.0_wp

   ! -----------------------------
   !  s functions (trafomat(1x1))
   ! -----------------------------

   trafomat(1,1) = 1.0_wp

   if ( maxl == 0 ) return

   ! Normalize the vector
   len = sqrt(sum(vec**2))
   norm_vec = vec / len

   ! Prepare spherical coordinats
   cost = norm_vec(3)
   if ( abs(cost) == 1.0_wp ) then
      ! Here, phi is arbitrary as the vector is parallel to the z-axis.
      ! We choose the x-axis as the arbitrary direction.
      sint = 0.0_wp
      cosp = 1.0_wp
      sinp = 0.0_wp
   else if ( abs(cost) == 0.0_wp ) then
      sint = 1.0_wp
      cosp = norm_vec(1)
      sinp = norm_vec(2)
   else
      sint = sqrt(norm_vec(1)**2 + norm_vec(2)**2)
      cosp = norm_vec(1)/sint
      sinp = norm_vec(2)/sint
   end if

   ! -----------------------------
   !  p functions (trafomat(4x4))
   ! -----------------------------

   ! Adapted to tblite ordering from MSINDO
   ! 1st index:
   ! (2,:)_MSINDO -> (px,:) -> (4,:)_tblite
   ! (3,:)_MSINDO -> (py,:) -> (2,:)_tblite
   ! (4,:)_MSINDO -> (pz,:) -> (3,:)_tblite
   trafomat(2, 2) = cosp
   trafomat(3, 2) = 0.0_wp
   trafomat(4, 2) = -sinp
   trafomat(2, 3) = sint*sinp
   trafomat(3, 3) = cost
   trafomat(4, 3) = sint*cosp
   trafomat(2, 4) = cost*sinp
   trafomat(3, 4) = -sint
   trafomat(4, 4) = cost*cosp

   if ( maxl <= 1 ) return

   ! -----------------------------
   !  d functions (trafomat(9x9))
   ! -----------------------------

   cos2t = cost**2 - sint**2
   sin2t = 2.0_wp * sint*cost
   cos2p = cosp**2 - sinp**2
   sin2p = 2.0_wp * sinp*cosp
   sqrt3 = sqrt(3.0_wp)

   ! Changed from MSINDO ordering (0,-1,1,-2,2) to tblite ordering (-2,-1,0,1,2) of d-functions
   ! (5,:)_MSINDO -> (dz2,:) -> (7,:)_tblite
   ! (6,:)_MSINDO -> (dxz,:) -> (8,:)_tblite
   ! (7,:)_MSINDO -> (dyz,:) -> (6,:)_tblite
   ! (8,:)_MSINDO -> (dx2-y2,:) -> (5,:)_tblite
   ! (9,:)_MSINDO -> (dxy,:) -> (9,:)_tblite

   trafomat(5, 5) = cost*cos2p
   trafomat(6, 5) = -sint*cosp
   trafomat(7, 5) = 0.0_wp
   trafomat(8, 5) = sint*sinp
   trafomat(9, 5) = -cost*sin2p
   trafomat(5, 6) = sint*cos2p
   trafomat(6, 6) = cost*cosp
   trafomat(7, 6) = 0.0_wp
   trafomat(8, 6) = -cost*sinp
   trafomat(9, 6) = -sint*sin2p
   trafomat(5, 7) = sqrt3*sint**2*sin2p*0.5_wp
   trafomat(6, 7) = sqrt3*sin2t*sinp*0.5_wp
   trafomat(7, 7) = (3.0_wp * cost**2 - 1.0_wp) * 0.5_wp
   trafomat(8, 7) = sqrt3*sin2t*cosp*0.5_wp
   trafomat(9, 7) = sqrt3*sint**2*cos2p*0.5_wp
   trafomat(5, 8) = sin2t*sin2p*0.5_wp
   trafomat(6, 8) = cos2t*sinp
   trafomat(7, 8) = -sqrt3*sin2t*0.5_wp
   trafomat(8, 8) = cos2t*cosp
   trafomat(9, 8) = sin2t*cos2p*0.5_wp
   trafomat(5, 9) = (1.0_wp + cost**2) * sin2p * 0.5_wp
   trafomat(6, 9) = -sin2t*sinp*0.5_wp
   trafomat(7, 9) = sqrt3*sint**2 * 0.5_wp
   trafomat(8, 9) = -sin2t*cosp*0.5_wp
   trafomat(9, 9) = (1.0_wp + cost**2) * cos2p * 0.5_wp

end subroutine harmtr

!> Setup diatomic frame transformation matrix and gradient for up to d-functions.
pure subroutine d_harmtr(maxl, vec, trafomat, trafomat_y, dtrafomat)
   !> Maximum angular momentum
   integer, intent(in) :: maxl
   !> Interatomic vector defining the diatomic frame z-axis
   real(wp), intent(in) :: vec(3)
   !> Transformation matrix for all components away from the z-axis singularity
   !> and for x and z component at the z-axis singularity.
   real(wp), intent(out) :: trafomat(:, :)
   !> Transformation matrix for the y component at the z-axis singularity.
   real(wp), intent(out) :: trafomat_y(:, :)
   !> Derivative of transformation matrix
   real(wp), intent(out) :: dtrafomat(:, :, :)

   ! Intermediate variables for the trigonometric functions
   ! Separte version for y for the case: vec || z-axis
   real(wp) :: cos2p, cos2t, cosp, cost, sin2p, sin2t, sinp
   real(wp) :: cos2py, cospy, sin2py, sinpy
   real(wp) :: dcos2t, dsin2t, dcos2p, dsin2p
   real(wp) :: dcos2py, dsin2py
   real(wp) :: dpdx, dpdy, dpdz, dtdx, dtdy, dtdz
   real(wp) :: norm_vec(3), sint, sqrt3, len

   trafomat = 0.0_wp
   trafomat_y = 0.0_wp
   dtrafomat = 0.0_wp

   ! -----------------------------
   !  s functions (trafomat(1x1))
   ! -----------------------------

   trafomat(1, 1) = 1.0_wp
   trafomat_y(1, 1) = 1.0_wp
   dtrafomat(1, 1, :) = 0.0_wp

   if ( maxl == 0 ) return

   ! Normalize the vector
   len = sqrt(sum(vec**2))
   norm_vec = vec / len

   ! Prepare spherical coordinats
   cost = norm_vec(3)
   if ( abs(cost) == 1.0_wp ) then
      sint = 0.0_wp
      ! Here, phi is arbitrary as the vector is parallel to the z-axis.
      ! In turn, the derivative is ill defined and has to be evaluated
      ! assuming either x or y orientation for the infinitesimal change.
      ! This is only require for the p- and d-functions, which are transformed.
      cosp = 1.0_wp
      sinp = 0.0_wp
      cospy = 0.0_wp
      sinpy = -1.0_wp
   else if ( abs(cost) == 0.0_wp ) then
      sint = 1.0_wp
      cosp = norm_vec(1)
      sinp = norm_vec(2)
      cospy = cosp
      sinpy = sinp
   else
      sint = sqrt(norm_vec(1)**2 + norm_vec(2)**2)
      cosp = norm_vec(1)/sint
      sinp = norm_vec(2)/sint
      cospy = cosp
      sinpy = sinp
   end if

   ! Prepare sperical coordinate derivative
   ! In the case of (exactly) vec || z-axis, the phi derivative vanishes.
   if( norm_vec(1)**2 + norm_vec(2)**2 == 0.0_wp ) then
      dpdx = 0.0_wp
      dpdy = 0.0_wp
   else
      dpdx = -sinp / sqrt(vec(1)**2 + vec(2)**2)
      dpdy = cospy / sqrt(vec(1)**2 + vec(2)**2)
   end if
   dpdz = 0.0_wp
   dtdx = cost * cosp / len
   dtdy = cost * sinpy / len
   dtdz = -sint / len

   ! -----------------------------
   !  p functions (trafomat(4x4))
   ! -----------------------------

   ! Adapted to tblite ordering from MSINDO
   ! 1st index:
   ! (2,:)_MSINDO -> (px,:) -> (4,:)_tblite
   ! (3,:)_MSINDO -> (py,:) -> (2,:)_tblite
   ! (4,:)_MSINDO -> (pz,:) -> (3,:)_tblite

   ! x- and z-direction
   trafomat(2, 2) = cosp
   trafomat(3, 2) = 0.0_wp
   trafomat(4, 2) = -sinp
   trafomat(2, 3) = sint*sinp
   trafomat(3, 3) = cost
   trafomat(4, 3) = sint*cosp
   trafomat(2, 4) = cost*sinp
   trafomat(3, 4) = -sint
   trafomat(4, 4) = cost*cosp
   ! y-direction
   trafomat_y(2, 2) = cospy
   trafomat_y(3, 2) = 0.0_wp
   trafomat_y(4, 2) = -sinpy
   trafomat_y(2, 3) = sint*sinpy
   trafomat_y(3, 3) = cost
   trafomat_y(4, 3) = sint*cospy
   trafomat_y(2, 4) = cost*sinpy
   trafomat_y(3, 4) = -sint
   trafomat_y(4, 4) = cost*cospy

   ! x-direction derivative w.r.t. phi and w.r.t. theta
   dtrafomat(2, 2, 1) = -sinp * dpdx
   dtrafomat(3, 2, 1) = 0.0_wp
   dtrafomat(4, 2, 1) = -cosp * dpdx
   dtrafomat(2, 3, 1) = sint*cosp * dpdx + cost*sinp * dtdx
   dtrafomat(3, 3, 1) = -sint * dtdx
   dtrafomat(4, 3, 1) = -sint*sinp * dpdx + cost*cosp * dtdx
   dtrafomat(2, 4, 1) = cost*cosp * dpdx - sint*sinp * dtdx
   dtrafomat(3, 4, 1) = -cost * dtdx
   dtrafomat(4, 4, 1) = -cost*sinp * dpdx - sint*cosp * dtdx
   ! y-direction derivative w.r.t. phi and w.r.t. theta
   dtrafomat(2, 2, 2) = -sinpy * dpdy
   dtrafomat(3, 2, 2) = 0.0_wp
   dtrafomat(4, 2, 2) = -cospy * dpdy
   dtrafomat(2, 3, 2) = sint*cospy * dpdy + cost*sinpy * dtdy
   dtrafomat(3, 3, 2) = -sint * dtdy
   dtrafomat(4, 3, 2) = -sint*sinpy * dpdy + cost*cospy * dtdy
   dtrafomat(2, 4, 2) = cost*cospy * dpdy - sint*sinpy * dtdy
   dtrafomat(3, 4, 2) = -cost * dtdy
   dtrafomat(4, 4, 2) = -cost*sinpy * dpdy - sint*cospy * dtdy
   ! z-direaction derivative w.r.t. theta
   dtrafomat(2, 2, 3) = 0.0_wp
   dtrafomat(3, 2, 3) = 0.0_wp
   dtrafomat(4, 2, 3) = 0.0_wp
   dtrafomat(2, 3, 3) = cost*sinp * dtdz
   dtrafomat(3, 3, 3) = -sint * dtdz
   dtrafomat(4, 3, 3) = cost*cosp * dtdz
   dtrafomat(2, 4, 3) = -sint*sinp * dtdz
   dtrafomat(3, 4, 3) = -cost * dtdz
   dtrafomat(4, 4, 3) = -sint*cosp * dtdz

   if ( maxl <= 1 ) return

   ! -----------------------------
   !  d functions (trafomat(9x9))
   ! -----------------------------

   sqrt3 = sqrt(3.0_wp)
   cos2t = cost**2 - sint**2
   sin2t = 2.0_wp * sint*cost
   cos2p = cosp**2 - sinp**2
   sin2p = 2.0_wp * sinp*cosp
   cos2py = cospy**2 - sinpy**2
   sin2py = 2.0_wp * sinpy*cospy

   dcos2t = -2.0_wp * sin2t
   dsin2t =  2.0_wp * cos2t
   dcos2p = -2.0_wp * sin2p
   dsin2p =  2.0_wp * cos2p
   dcos2py = -2.0_wp * sin2py
   dsin2py =  2.0_wp * cos2py

   ! Changed from MSINDO ordering (0,-1,1,-2,2) to tblite ordering (-2,-1,0,1,2) of d-functions
   ! (5,:)_MSINDO -> (dz2,:) -> (7,:)_tblite
   ! (6,:)_MSINDO -> (dxz,:) -> (8,:)_tblite
   ! (7,:)_MSINDO -> (dyz,:) -> (6,:)_tblite
   ! (8,:)_MSINDO -> (dx2-y2,:) -> (5,:)_tblite
   ! (9,:)_MSINDO -> (dxy,:) -> (9,:)_tblite
   ! x- and z-direction
   trafomat(5, 5) = cost*cos2p
   trafomat(6, 5) = -sint*cosp
   trafomat(7, 5) = 0.0_wp
   trafomat(8, 5) = sint*sinp
   trafomat(9, 5) = -cost*sin2p
   trafomat(5, 6) = sint*cos2p
   trafomat(6, 6) = cost*cosp
   trafomat(7, 6) = 0.0_wp
   trafomat(8, 6) = -cost*sinp
   trafomat(9, 6) = -sint*sin2p
   trafomat(5, 7) = sqrt3*sint**2*sin2p*0.5_wp
   trafomat(6, 7) = sqrt3*sin2t*sinp*0.5_wp
   trafomat(7, 7) = (3.0_wp * cost**2 - 1.0_wp) * 0.5_wp
   trafomat(8, 7) = sqrt3*sin2t*cosp*0.5_wp
   trafomat(9, 7) = sqrt3*sint**2*cos2p*0.5_wp
   trafomat(5, 8) = sin2t*sin2p*0.5_wp
   trafomat(6, 8) = cos2t*sinp
   trafomat(7, 8) = -sqrt3*sin2t*0.5_wp
   trafomat(8, 8) = cos2t*cosp
   trafomat(9, 8) = sin2t*cos2p*0.5_wp
   trafomat(5, 9) = (1.0_wp + cost**2) * sin2p * 0.5_wp
   trafomat(6, 9) = -sin2t*sinp*0.5_wp
   trafomat(7, 9) = sqrt3*sint**2 * 0.5_wp
   trafomat(8, 9) = -sin2t*cosp*0.5_wp
   trafomat(9, 9) = (1.0_wp + cost**2) * cos2p * 0.5_wp
   ! y-direction
   trafomat_y(5, 5) = cost*cos2py
   trafomat_y(6, 5) = -sint*cospy
   trafomat_y(7, 5) = 0.0_wp
   trafomat_y(8, 5) = sint*sinpy
   trafomat_y(9, 5) = -cost*sin2py
   trafomat_y(5, 6) = sint*cos2py
   trafomat_y(6, 6) = cost*cospy
   trafomat_y(7, 6) = 0.0_wp
   trafomat_y(8, 6) = -cost*sinpy
   trafomat_y(9, 6) = -sint*sin2py
   trafomat_y(5, 7) = sqrt3*sint**2*sin2py*0.5_wp
   trafomat_y(6, 7) = sqrt3*sin2t*sinpy*0.5_wp
   trafomat_y(7, 7) = (3.0_wp * cost**2 - 1.0_wp) * 0.5_wp
   trafomat_y(8, 7) = sqrt3*sin2t*cospy*0.5_wp
   trafomat_y(9, 7) = sqrt3*sint**2*cos2py*0.5_wp
   trafomat_y(5, 8) = sin2t*sin2py*0.5_wp
   trafomat_y(6, 8) = cos2t*sinpy
   trafomat_y(7, 8) = -sqrt3*sin2t*0.5_wp
   trafomat_y(8, 8) = cos2t*cospy
   trafomat_y(9, 8) = sin2t*cos2py*0.5_wp
   trafomat_y(5, 9) = (1.0_wp + cost**2) * sin2py * 0.5_wp
   trafomat_y(6, 9) = -sin2t*sinpy*0.5_wp
   trafomat_y(7, 9) = sqrt3*sint**2 * 0.5_wp
   trafomat_y(8, 9) = -sin2t*cospy*0.5_wp
   trafomat_y(9, 9) = (1.0_wp + cost**2) * cos2py * 0.5_wp

   ! x-direction derivative w.r.t. phi and w.r.t. theta
   dtrafomat(5, 5, 1) = cost*dcos2p * dpdx - sint*cos2p * dtdx
   dtrafomat(6, 5, 1) = sint*sinp * dpdx - cost*cosp * dtdx
   dtrafomat(7, 5, 1) = 0.0_wp
   dtrafomat(8, 5, 1) = sint*cosp * dpdx + cost*sinp * dtdx
   dtrafomat(9, 5, 1) = -cost*dsin2p * dpdx + sint*sin2p * dtdx
   dtrafomat(5, 6, 1) = sint*dcos2p * dpdx + cost*cos2p * dtdx
   dtrafomat(6, 6, 1) = -cost*sinp * dpdx - sint*cosp * dtdx
   dtrafomat(7, 6, 1) = 0.0_wp
   dtrafomat(8, 6, 1) = -cost*cosp * dpdx + sint*sinp * dtdx
   dtrafomat(9, 6, 1) = -sint*dsin2p * dpdx - cost*sin2p * dtdx
   dtrafomat(5, 7, 1) = 0.5_wp*sqrt3 * (sint**2*dsin2p * dpdx + sin2t*sin2p * dtdx)
   dtrafomat(6, 7, 1) = 0.5_wp*sqrt3 * (sin2t*cosp * dpdx + dsin2t*sinp * dtdx)
   dtrafomat(7, 7, 1) = -3.0_wp*sin2t*0.5_wp * dtdx
   dtrafomat(8, 7, 1) = 0.5_wp*sqrt3 * (-sin2t*sinp * dpdx + dsin2t*cosp * dtdx)
   dtrafomat(9, 7, 1) = 0.5_wp*sqrt3 * (sint**2*dcos2p * dpdx + sin2t*cos2p * dtdx)
   dtrafomat(5, 8, 1) = 0.5_wp * (sin2t*dsin2p * dpdx + dsin2t*sin2p * dtdx)
   dtrafomat(6, 8, 1) = cos2t*cosp * dpdx + dcos2t*sinp * dtdx
   dtrafomat(7, 8, 1) = -sqrt3*dsin2t*0.5_wp * dtdx
   dtrafomat(8, 8, 1) = -cos2t*sinp * dpdx + dcos2t*cosp * dtdx
   dtrafomat(9, 8, 1) = 0.5_wp * (sin2t*dcos2p * dpdx + dsin2t*cos2p * dtdx)
   dtrafomat(5, 9, 1) = 0.5_wp * ((1.0_wp + cost**2)*dsin2p * dpdx - sin2t*sin2p * dtdx)
   dtrafomat(6, 9, 1) = 0.5_wp * (-sin2t*cosp * dpdx - dsin2t*sinp * dtdx)
   dtrafomat(7, 9, 1) = sqrt3*sin2t*0.5_wp * dtdx
   dtrafomat(8, 9, 1) = 0.5_wp * (sin2t*sinp * dpdx - dsin2t*cosp * dtdx)
   dtrafomat(9, 9, 1) = 0.5_wp * ((1.0_wp + cost**2)*dcos2p * dpdx - sin2t*cos2p * dtdx)
   ! y-direction derivative w.r.t. phi and w.r.t. theta
   dtrafomat(5, 5, 2) = cost*dcos2py * dpdy - sint*cos2py * dtdy
   dtrafomat(6, 5, 2) = sint*sinpy * dpdy - cost*cospy * dtdy
   dtrafomat(7, 5, 2) = 0.0_wp
   dtrafomat(8, 5, 2) = sint*cospy * dpdy + cost*sinpy * dtdy
   dtrafomat(9, 5, 2) = -cost*dsin2py * dpdy + sint*sin2py * dtdy
   dtrafomat(5, 6, 2) = sint*dcos2py * dpdy + cost*cos2py * dtdy
   dtrafomat(6, 6, 2) = -cost*sinpy * dpdy - sint*cospy * dtdy
   dtrafomat(7, 6, 2) = 0.0_wp
   dtrafomat(8, 6, 2) = -cost*cospy * dpdy + sint*sinpy * dtdy
   dtrafomat(9, 6, 2) = -sint*dsin2py * dpdy - cost*sin2py * dtdy
   dtrafomat(5, 7, 2) = 0.5_wp*sqrt3 * (sint**2*dsin2py * dpdy + sin2t*sin2py * dtdy)
   dtrafomat(6, 7, 2) = 0.5_wp*sqrt3 * (sin2t*cospy * dpdy + dsin2t*sinpy * dtdy)
   dtrafomat(7, 7, 2) = -3.0_wp*sin2t*0.5_wp * dtdy
   dtrafomat(8, 7, 2) = 0.5_wp*sqrt3 * (-sin2t*sinpy * dpdy + dsin2t*cospy * dtdy)
   dtrafomat(9, 7, 2) = 0.5_wp*sqrt3 * (sint**2*dcos2py * dpdy + sin2t*cos2py * dtdy)
   dtrafomat(5, 8, 2) = 0.5_wp * (sin2t*dsin2py * dpdy + dsin2t*sin2py * dtdy)
   dtrafomat(6, 8, 2) = cos2t*cospy * dpdy + dcos2t*sinpy * dtdy
   dtrafomat(7, 8, 2) = -sqrt3*dsin2t*0.5_wp * dtdy
   dtrafomat(8, 8, 2) = -cos2t*sinpy * dpdy + dcos2t*cospy * dtdy
   dtrafomat(9, 8, 2) = 0.5_wp * (sin2t*dcos2py * dpdy + dsin2t*cos2py * dtdy)
   dtrafomat(5, 9, 2) = 0.5_wp * ((1.0_wp + cost**2)*dsin2py * dpdy - sin2t*sin2py * dtdy)
   dtrafomat(6, 9, 2) = 0.5_wp * (-sin2t*cospy * dpdy - dsin2t*sinpy * dtdy)
   dtrafomat(7, 9, 2) = sqrt3*sin2t*0.5_wp * dtdy
   dtrafomat(8, 9, 2) = 0.5_wp * (sin2t*sinpy * dpdy - dsin2t*cospy * dtdy)
   dtrafomat(9, 9, 2) = 0.5_wp * ((1.0_wp + cost**2)*dcos2py * dpdy - sin2t*cos2py * dtdy)
   ! z-direaction derivative w.r.t. theta
   dtrafomat(5, 5, 3) = -sint*cos2p * dtdz
   dtrafomat(6, 5, 3) = -cost*cosp * dtdz
   dtrafomat(7, 5, 3) = 0.0_wp
   dtrafomat(8, 5, 3) = cost*sinp * dtdz
   dtrafomat(9, 5, 3) = sint*sin2p * dtdz
   dtrafomat(5, 6, 3) = cost*cos2p * dtdz
   dtrafomat(6, 6, 3) = -sint*cosp * dtdz
   dtrafomat(7, 6, 3) = 0.0_wp
   dtrafomat(8, 6, 3) = sint*sinp * dtdz
   dtrafomat(9, 6, 3) = -cost*sin2p * dtdz
   dtrafomat(5, 7, 3) = sqrt3*sin2t*sin2p*0.5_wp * dtdz
   dtrafomat(6, 7, 3) = sqrt3*dsin2t*sinp*0.5_wp * dtdz
   dtrafomat(7, 7, 3) = -3.0_wp*sin2t*0.5_wp * dtdz
   dtrafomat(8, 7, 3) = sqrt3*dsin2t*cosp*0.5_wp * dtdz
   dtrafomat(9, 7, 3) = sqrt3*sin2t*cos2p*0.5_wp * dtdz
   dtrafomat(5, 8, 3) = dsin2t*sin2p*0.5_wp * dtdz
   dtrafomat(6, 8, 3) = dcos2t*sinp * dtdz
   dtrafomat(7, 8, 3) = -sqrt3*dsin2t*0.5_wp * dtdz
   dtrafomat(8, 8, 3) = dcos2t*cosp * dtdz
   dtrafomat(9, 8, 3) = dsin2t*cos2p*0.5_wp * dtdz
   dtrafomat(5, 9, 3) = -sin2t*sin2p*0.5_wp * dtdz
   dtrafomat(6, 9, 3) = -dsin2t*sinp*0.5_wp * dtdz
   dtrafomat(7, 9, 3) = sqrt3*sin2t*0.5_wp * dtdz
   dtrafomat(8, 9, 3) = -dsin2t*cosp*0.5_wp * dtdz
   dtrafomat(9, 9, 3) = -sin2t*cos2p*0.5_wp * dtdz

end subroutine d_harmtr

pure subroutine scale_diatomic_frame(diat_mat, ksig, kpi, kdel, maxlj, maxli)
   !> Block matrix in the diatomic frame to be scaled
   real(wp),intent(inout)    :: diat_mat(:,:)
   !> Scaling parameters for different bonding contributions
   real(wp),intent(in)       :: ksig, kpi, kdel
   !> Highest angular momentum of atom j (first index)
   integer,intent(in)        :: maxlj
   !> Highest angular momentum of atom i (second index)
   integer,intent(in)        :: maxli

   diat_mat(1,1) = diat_mat(1,1)*ksig ! Sigma bond s   <-> s
   if(maxlj > 0) then
      diat_mat(3,1) = diat_mat(3,1)*ksig ! Sigma bond pz  <-> s
   end if
   if(maxli > 0)  then
      diat_mat(1,3) = diat_mat(1,3)*ksig ! Sigma bond s   <-> pz
   end if
   if(maxlj > 0 .and. maxli > 0) then
      diat_mat(3,3) = diat_mat(3,3)*ksig ! Sigma bond pz  <-> pz
      diat_mat(4,4) = diat_mat(4,4)*kpi  ! Pi    bond px  <-> px
      diat_mat(2,2) = diat_mat(2,2)*kpi  ! Pi    bond py  <-> py
      if(maxlj > 1) then
         diat_mat(7,3) = diat_mat(7,3)*ksig ! Sigma bond dz2 <-> pz
         diat_mat(8,4) = diat_mat(8,4)*kpi  ! Pi    bond dxz <-> px
         diat_mat(6,2) = diat_mat(6,2)*kpi  ! Pi    bond dyz <-> py
      end if
      if(maxli > 1) then
         diat_mat(3,7) = diat_mat(3,7)*ksig ! Sigma bond pz  <-> dz2
         diat_mat(4,8) = diat_mat(4,8)*kpi  ! Pi    bond px  <-> dxz
         diat_mat(2,6) = diat_mat(2,6)*kpi  ! Pi    bond py  <-> dyz
      end if
   end if
   if (maxlj > 1) then
      diat_mat(7,1) = diat_mat(7,1)*ksig ! Sigma bond dz2 <-> s
   end if
   if (maxli > 1) then
      diat_mat(1,7) = diat_mat(1,7)*ksig ! Sigma bond s   <-> dz2
   end if
   if (maxlj > 1 .and. maxli > 1) then
      diat_mat(7,7) = diat_mat(7,7)*ksig ! Sigma bond dz2 <-> dz2
      diat_mat(8,8) = diat_mat(8,8)*kpi  ! Pi    bond dxz <-> dxz
      diat_mat(6,6) = diat_mat(6,6)*kpi  ! Pi    bond dyz <-> dyz
      diat_mat(9,9) = diat_mat(9,9)*kdel ! Delta bond dxy <-> dxy
      diat_mat(5,5) = diat_mat(5,5)*kdel ! Delta bond dx2-y2 <-> dx2-y2
   end if
   ! f- and g-functions remain unscaled

end subroutine scale_diatomic_frame

end module tblite_integral_diat_trafo
