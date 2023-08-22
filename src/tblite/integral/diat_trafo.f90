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
!> Provides help tools for evaluation of the diatomic overlap
module tblite_integral_diat_trafo
   use mctc_env, only : wp
   use tblite_blas, only: gemm

   implicit none
   private

   public :: relvec, diat_trafo

contains

   subroutine diat_trafo(block_overlap, vec_diat_trafo, ksig, kpi, kdel)
      !> Transformation vector for the diatomic frame
      real(wp),intent(inout)    :: block_overlap(9,9)
      real(wp),intent(in)       :: vec_diat_trafo(3)
      real(wp),intent(in)       :: ksig, kpi, kdel
      real(wp) :: trafomat(9,9), trans_block_s(9,9),tmp(9,9)

      !> 1. Calculate the transformation matrix
      call harmtr(2, vec_diat_trafo, trafomat)
      ! Doing it only for the maxl_ish_jsh elements (for s,p -> 4x4) would in
      ! theory be more efficient, but is not implemented yet:
      ! call harmtr(maxl_ish_jsh, vec_diat_trafo, trafomat)

      !> 2. Transform the submatrix
      ! trans_block_s = matmul(matmul(transpose(trafomat), block_overlap),trafomat)
      ! trans_block_s = O^T * S * O
      call gemm(amat=trafomat,bmat=block_overlap,cmat=tmp,transa='T',transb='N')
      call gemm(amat=tmp,bmat=trafomat,cmat=trans_block_s,transa='N',transb='N')

      !> 3.1. Scale elements in diatomic frame
      !> 3.2. Scale elements with equivalent bonding situation in the
      !>      diatomic frame.
      trans_block_s(1,1) = trans_block_s(1,1)*ksig ! Sigma bond s   <-> s
      trans_block_s(1,3) = trans_block_s(1,3)*ksig ! Sigma bond s   <-> pz
      trans_block_s(3,1) = trans_block_s(3,1)*ksig ! Sigma bond pz  <-> s
      trans_block_s(1,5) = trans_block_s(1,5)*ksig ! Sigma bond s   <-> dz2
      trans_block_s(5,1) = trans_block_s(5,1)*ksig ! Sigma bond dz2 <-> s
      trans_block_s(3,3) = trans_block_s(3,3)*ksig ! Sigma bond pz  <-> pz
      trans_block_s(3,5) = trans_block_s(3,5)*ksig ! Sigma bond pz  <-> dz2
      trans_block_s(5,3) = trans_block_s(5,3)*ksig ! Sigma bond dz2 <-> pz
      trans_block_s(5,5) = trans_block_s(5,5)*ksig ! Sigma bond dz2 <-> dz2
      trans_block_s(4,4) = trans_block_s(4,4)*kpi  ! Pi    bond px  <-> px
      trans_block_s(4,6) = trans_block_s(4,6)*kpi  ! Pi    bond px  <-> dxz
      trans_block_s(6,4) = trans_block_s(6,4)*kpi  ! Pi    bond dxz <-> px
      trans_block_s(6,6) = trans_block_s(6,6)*kpi  ! Pi    bond dxz <-> dxz
      trans_block_s(2,2) = trans_block_s(2,2)*kpi  ! Pi    bond py  <-> py
      trans_block_s(2,7) = trans_block_s(2,7)*kpi  ! Pi    bond py  <-> dyz
      trans_block_s(7,2) = trans_block_s(7,2)*kpi  ! Pi    bond dyz <-> py
      trans_block_s(7,7) = trans_block_s(7,7)*kpi  ! Pi    bond dyz <-> dyz
      trans_block_s(8,8) = trans_block_s(8,8)*kdel ! Delta bond dx2-y2 <-> dx2-y2
      trans_block_s(9,9) = trans_block_s(9,9)*kdel ! Delta bond dxy <-> dxy

      !> 4. Transform back to original frame
      ! block_overlap = matmul(matmul(trafomat, trans_block_s),transpose(trafomat))
      ! block_overlap = O * S * O^T
      call gemm(amat=trafomat,bmat=trans_block_s,cmat=tmp,transa='N',transb='N')
      call gemm(amat=tmp,bmat=trafomat,cmat=block_overlap,transa='N',transb='T')
   end subroutine diat_trafo

   subroutine relvec(vec, rkl, veckl)

      !> Original vector between atoms A and B
      real(wp), intent(in)             :: vec(3)
      !> Distance between the two atoms
      real(wp), intent(in)             :: rkl
      !> Normalized vector from atom k to atom l
      real(wp), intent(out)            :: veckl(3)

      real(wp), parameter              :: eps = 4.0e-08_wp

      real(wp)                         :: sq

      veckl(1:3) = vec(1:3) / rkl
      if ( abs(1.0_wp-abs(veckl(1))) .lt. eps ) then
         veckl(1) = sign(1.0_wp,veckl(1))
         veckl(2) = 0.0_wp
         veckl(3) = 0.0_wp
      else if ( abs(1.0_wp-abs(veckl(2))) .lt. eps ) then
         veckl(1) = 0.0_wp
         veckl(2) = sign(1.0_wp,veckl(2))
         veckl(3) = 0.0_wp
      else if ( abs(1.0_wp-abs(veckl(3))) .lt. eps ) then
         veckl(1) = 0.0_wp
         veckl(2) = 0.0_wp
         veckl(3) = sign(1.0_wp,veckl(3))
      else if ( (abs(veckl(1)) .lt. eps) .and. .not. eff_equality(veckl(1),0.0_wp) ) then
         veckl(1) = 0.0_wp
         sq = sqrt( veckl(2)**2 + veckl(3)**2 )
         veckl(2) = veckl(2)/sq
         veckl(3) = veckl(3)/sq
      else if ( (abs(veckl(2)) .lt. eps) .and. .not. eff_equality(veckl(2),0.0_wp) ) then
         veckl(2) = 0.0_wp
         sq = sqrt( veckl(1)**2 + veckl(3)**2 )
         veckl(1) = veckl(1)/sq
         veckl(3) = veckl(3)/sq
      else if ( (abs(veckl(3)) .lt. eps) .and. .not. eff_equality(veckl(3),0.0_wp) ) then
         veckl(3) = 0.0_wp
         sq = sqrt(veckl(1)**2 + veckl(2)**2)
         veckl(1) = veckl(1)/sq
         veckl(2) = veckl(2)/sq
      endif

   end subroutine relvec

   logical pure function eff_equality(num1, num2)
      !> Numbers to compare
      real(wp), intent(in) :: num1, num2
      !> Logical deciding if numbers are (almost) equal or not
      eff_equality = (abs( num1 - num2 ) .le. 1.0e-12_wp)

   end function eff_equality

   pure subroutine harmtr(maxkl,veckl,trafomat)
      !> Maximum angular momentum
      integer, intent(in)  :: maxkl
      !> Normalized vector from atom k to atom l
      real(wp), intent(in) :: veckl(3)
      !> Transformation matrix
      real(wp), intent(out) :: trafomat(9,9)

      real(wp) :: cos2p, cos2t, cosp, cost, sin2p, sin2t, sinp, sint, sqrt3

      !     ------------------------------------------------------------------
      if (maxkl > 2) then
         error stop "ERROR: f function or higher ang. mom. not implemented in harmtr"
      endif

      trafomat = 0.0_wp
      !     -----------------------------
      !     *** s functions (trafomat(1x1)) ***
      !     -----------------------------
      trafomat(1,1) = 1.0

      if ( maxkl == 0 ) return
      !     -----------------------------
      !     *** p functions (trafomat(4x4)) ***
      !     -----------------------------

      cost = veckl(3)
      if ( abs(cost) .eq. 1.0_wp ) then
         sint = 0.0_wp
         cosp = 1.0_wp
         sinp = 0.0_wp
      else if ( abs(cost) .eq. 0.0_wp ) then
         sint = 1.0_wp
         cosp = veckl(1)
         sinp = veckl(2)
      else
         sint = SQRT(1.0_wp-COST**2)
         cosp = veckl(1)/SINT
         sinp = veckl(2)/SINT
      endif

      !> tblite ordering with adapted column ordering
      ! 1st index:
      ! MSINDO defintion of p function ordering is converted to
      ! tblite definition of p function ordering. E.g. for first entry:
      ! trafomat(2,:)_MSINDO -> trafomat(px,:) -> trafomat(4:)_tblite
      ! 2nd index:
      ! Final ordering of p functions (see below) corresponds to the
      ! tblite ordering of p functions. For the second index, the ordering
      ! 3, 4, 2 holds, corresponding (in MSINDO convention)
      ! to the tblite ordering of p functions, i.e.
      ! y, z, x
      trafomat(4,3) = SINT*COSP
      trafomat(2,3) = SINT*SINP
      trafomat(3,3) = COST
      trafomat(4,4) = COST*COSP
      trafomat(2,4) = COST*SINP
      trafomat(3,4) = -SINT
      trafomat(4,2) = -SINP
      trafomat(2,2) = COSP
      trafomat(3,2) = 0.0_wp

      if ( maxkl <= 1 ) return

!     -----------------------------
!     *** d functions (trafomat(9x9)) ***
!     -----------------------------

      COS2T = COST**2 - SINT**2
      SIN2T = 2.0_wp * SINT*COST
      COS2P = COSP**2 - SINP**2
      SIN2P = 2.0_wp * SINP*COSP
      SQRT3 = SQRT(3.0_wp)

      !> Original MSINDO ordering
      !> The MSINDO d SAO ordering corresponds to the
      !> tblite ordering of d SAOs
      trafomat(5,5) = (3.0_wp * COST**2 - 1.0_wp) * 0.5_wp
      trafomat(6,5) = SQRT3*SIN2T*COSP*0.5_wp
      trafomat(7,5) = SQRT3*SIN2T*SINP*0.5_wp
      trafomat(8,5) = SQRT3*SINT**2*COS2P*0.5_wp
      trafomat(9,5) = SQRT3*SINT**2*SIN2P*0.5_wp
      trafomat(5,6) = -SQRT3*SIN2T*0.5_wp
      trafomat(6,6) = COS2T*COSP
      trafomat(7,6) = COS2T*SINP
      trafomat(8,6) = SIN2T*COS2P*0.5_wp
      trafomat(9,6) = SIN2T*SIN2P*0.5_wp
      trafomat(5,7) = 0.0_wp
      trafomat(6,7) = -COST*SINP
      trafomat(7,7) = COST*COSP
      trafomat(8,7) = -SINT*SIN2P
      trafomat(9,7) = SINT*COS2P
      trafomat(5,8) = SQRT3*SINT**2 * 0.5_wp
      trafomat(6,8) = -SIN2T*COSP*0.5_wp
      trafomat(7,8) = -SIN2T*SINP*0.5_wp
      trafomat(8,8) = (1.0_wp + COST**2) * COS2P * 0.5_wp
      trafomat(9,8) = (1.0_wp + COST**2) * SIN2P * 0.5_wp
      trafomat(5,9) = 0.0_wp
      trafomat(6,9) = SINT*SINP
      trafomat(7,9) = -SINT*COSP
      trafomat(8,9) = -COST*SIN2P
      trafomat(9,9) = COST*COS2P

   end subroutine harmtr
end module tblite_integral_diat_trafo
