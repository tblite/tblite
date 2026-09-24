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

!> @file tblite/integral/trafo.f90
!> Provides transformation from cartesian to spherical harmonic basis functions

!> Implementation of transformations from cartesian to spherical harmonic basis functions
!> and adjoint transformation for contravariant vectors from spherical harmonic
!> to cartesian basis functions, as well as the inverse of the latter.
!>
!> Spherical harmonics use standard ordering, *i.e.* [-l, ..., 0, ..., l].
module tblite_integral_trafo
   use mctc_env, only : wp
   implicit none
   private

   public :: transform0, transform1, transform2
   public :: adjoint_transform0, adjoint_transform1, adjoint_transform2
   public :: contravariant_transform0, contravariant_transform1, contravariant_transform2


   real(wp), parameter :: d1_3 = 1.0_wp/3.0_wp
   real(wp), parameter :: d2_3 = 2.0_wp/3.0_wp
   real(wp), parameter :: s3 = sqrt(3.0_wp)
   real(wp), parameter :: s1_3 = s3/3.0_wp
   real(wp), parameter :: s3_4 = s3 * 0.5_wp
   real(wp), parameter :: dtrafo(5, 6) = reshape([&
      ! -2      -1       0       1       2
      & 0.0_wp, 0.0_wp, -0.5_wp, 0.0_wp,   s3_4, & ! xx
      &     s3, 0.0_wp,  0.0_wp, 0.0_wp, 0.0_wp, & ! xy
      & 0.0_wp, 0.0_wp,  0.0_wp,     s3, 0.0_wp, & ! xz
      & 0.0_wp, 0.0_wp, -0.5_wp, 0.0_wp,  -s3_4, & ! yy
      & 0.0_wp,     s3,  0.0_wp, 0.0_wp, 0.0_wp, & ! yz
      & 0.0_wp, 0.0_wp,  1.0_wp, 0.0_wp, 0.0_wp],& ! zz
      & shape(dtrafo))
   real(wp), parameter :: dtrafo_pinv(5, 6) = reshape([&
      ! -2      -1      0       1       2
      & 0.0_wp, 0.0_wp,  -d1_3, 0.0_wp,   s1_3, & ! xx
      &   s1_3, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, & ! xy
      & 0.0_wp, 0.0_wp, 0.0_wp,   s1_3, 0.0_wp, & ! xz
      & 0.0_wp, 0.0_wp,  -d1_3, 0.0_wp,  -s1_3, & ! yy
      & 0.0_wp,   s1_3, 0.0_wp, 0.0_wp, 0.0_wp, & ! yz
      & 0.0_wp, 0.0_wp,   d2_3, 0.0_wp, 0.0_wp],& ! zz
      & shape(dtrafo_pinv))

   real(wp), parameter :: d3_2 = 3.0_wp/2.0_wp
   real(wp), parameter :: d3_11 = 3.0_wp/11.0_wp
   real(wp), parameter :: d2_11 = 2.0_wp/11.0_wp
   real(wp), parameter :: s3_8 = sqrt(3.0_wp/8.0_wp)
   real(wp), parameter :: s5_8 = sqrt(5.0_wp/8.0_wp)
   real(wp), parameter :: s6 = sqrt(6.0_wp)
   real(wp), parameter :: s3_242 = s6/22.0_wp
   real(wp), parameter :: s1_726 = s6/66.0_wp
   real(wp), parameter :: s50_363 = 5.0_wp*s6/33.0_wp
   real(wp), parameter :: s10 = sqrt(10.0_wp)
   real(wp), parameter :: s5_242 = s10/22.0_wp
   real(wp), parameter :: s2_605 = s10/55.0_wp
   real(wp), parameter :: s15 = sqrt(15.0_wp)
   real(wp), parameter :: s15_4 = s15/2.0_wp
   real(wp), parameter :: s1_15 = s15/15.0_wp
   real(wp), parameter :: s45_8 = 3.0_wp*s10/4.0_wp
   real(wp), parameter :: s169_1210 = 13.0_wp*s10/110.0_wp
   real(wp), parameter :: ftrafo(7, 10) = reshape([&
      ! -3       -2       -1       0        1         2         3
      &  0.0_wp,  0.0_wp,  0.0_wp, 0.0_wp,   -s3_8,   0.0_wp,     s5_8, & ! xxx
      &   s45_8,  0.0_wp,   -s3_8, 0.0_wp,  0.0_wp,   0.0_wp,   0.0_wp, & ! xxy
      &  0.0_wp,  0.0_wp,  0.0_wp,  -d3_2,  0.0_wp,    s15_4,   0.0_wp, & ! xxz
      &  0.0_wp,  0.0_wp,  0.0_wp, 0.0_wp,   -s3_8,   0.0_wp,   -s45_8, & ! xyy
      &  0.0_wp,     s15,  0.0_wp, 0.0_wp,  0.0_wp,   0.0_wp,   0.0_wp, & ! xyz
      &  0.0_wp,  0.0_wp,  0.0_wp, 0.0_wp,      s6,   0.0_wp,   0.0_wp, & ! xzz
      &   -s5_8,  0.0_wp,   -s3_8, 0.0_wp,  0.0_wp,   0.0_wp,   0.0_wp, & ! yyy
      &  0.0_wp,  0.0_wp,  0.0_wp,  -d3_2,  0.0_wp,   -s15_4,   0.0_wp, & ! yyz
      &  0.0_wp,  0.0_wp,      s6, 0.0_wp,  0.0_wp,   0.0_wp,   0.0_wp, & ! yzz
      &  0.0_wp,  0.0_wp,  0.0_wp, 1.0_wp,  0.0_wp,   0.0_wp,   0.0_wp],& ! zzz
      & shape(ftrafo))
   real(wp), parameter :: ftrafo_pinv(7, 10) = reshape([&
      ! -3         -2      -1       0       1        2       3
      &    0.0_wp, 0.0_wp,  0.0_wp, 0.0_wp, -s3_242, 0.0_wp,     s5_242, & ! xxx
      & s169_1210, 0.0_wp, -s1_726, 0.0_wp,  0.0_wp, 0.0_wp,     0.0_wp, & ! xxy
      &    0.0_wp, 0.0_wp,  0.0_wp, -d3_11,  0.0_wp,  s1_15,     0.0_wp, & ! xxz
      &    0.0_wp, 0.0_wp,  0.0_wp, 0.0_wp, -s1_726, 0.0_wp, -s169_1210, & ! xyy
      &    0.0_wp,  s1_15,  0.0_wp, 0.0_wp,  0.0_wp, 0.0_wp,     0.0_wp, & ! xyz
      &    0.0_wp, 0.0_wp,  0.0_wp, 0.0_wp, s50_363, 0.0_wp,    -s2_605, & ! xzz
      &   -s5_242, 0.0_wp, -s3_242, 0.0_wp,  0.0_wp, 0.0_wp,     0.0_wp, & ! yyy
      &    0.0_wp, 0.0_wp,  0.0_wp, -d3_11,  0.0_wp, -s1_15,     0.0_wp, & ! yyz
      &    s2_605, 0.0_wp, s50_363, 0.0_wp,  0.0_wp, 0.0_wp,     0.0_wp, & ! yzz
      &    0.0_wp, 0.0_wp,  0.0_wp,  d2_11,  0.0_wp, 0.0_wp,     0.0_wp],& ! zzz
      & shape(ftrafo_pinv))

   real(wp), parameter :: d3_8 = 3.0_wp/8.0_wp
   real(wp), parameter :: d3_4 = 3.0_wp/4.0_wp
   real(wp), parameter :: d3_370 = 3.0_wp/370.0_wp
   real(wp), parameter :: d9_370 = 9.0_wp/370.0_wp
   real(wp), parameter :: d19_370 = 19.0_wp/370.0_wp
   real(wp), parameter :: d57_370 = 57.0_wp/370.0_wp
   real(wp), parameter :: s5 = sqrt(5.0_wp)
   real(wp), parameter :: s5_4 = s5/2.0_wp
   real(wp), parameter :: s5_16 = s5/4.0_wp
   real(wp), parameter :: s4_6845 = 2.0_wp*s5/185.0_wp
   real(wp), parameter :: s144_6845 = 12.0_wp*s5/185.0_wp
   real(wp), parameter :: s1_1805 = s5/95.0_wp
   real(wp), parameter :: s36_1805 = 6.0_wp*s5/95.0_wp
   real(wp), parameter :: s9_3610 = 3.0_wp*s10/190.0_wp
   real(wp), parameter :: s81_3610 = 9.0_wp*s10/190.0_wp
   real(wp), parameter :: s10_361 = s10/19.0_wp
   real(wp), parameter :: s35 = sqrt(35.0_wp)
   real(wp), parameter :: s35_4 = s35/2.0_wp
   real(wp), parameter :: s35_8 = sqrt(35.0_wp/8.0_wp)
   real(wp), parameter :: s35_64 = s35/8.0_wp
   real(wp), parameter :: s1_35 = s35/35.0_wp
   real(wp), parameter :: s1_191660 = s35/2590.0_wp
   real(wp), parameter :: s9_191660 = 3.0_wp*s1_191660
   real(wp), parameter :: s81_191660 = 9.0_wp*s1_191660
   real(wp), parameter :: s289_191660 = 17.0_wp*s1_191660
   real(wp), parameter :: s8649_191660 = 93.0_wp*s1_191660
   real(wp), parameter :: s315_16 = 3.0_wp*s35/4.0_wp
   real(wp), parameter :: s45 = sqrt(45.0_wp)
   real(wp), parameter :: s45_4 = s45/2.0_wp
   real(wp), parameter :: s70 = sqrt(70.0_wp)
   real(wp), parameter :: s63_3610 = 3.0_wp*s70/190.0_wp
   real(wp), parameter :: s18_12635 = 3.0_wp*s70/665.0_wp
   real(wp), parameter :: s169_25270 = 13.0_wp*s70/1330.0_wp
   real(wp), parameter :: s315_8 = 3.0_wp*s70/4.0_wp
   real(wp), parameter :: gtrafo(9, 15) = reshape([&
      !  -4      -3      -2      -1      0       1       2        3         4
      &  0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,   d3_8, 0.0_wp, -s5_16,  0.0_wp,   s35_64, & ! xxxx
      &   s35_4, 0.0_wp,  -s5_4, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,  0.0_wp,   0.0_wp, & ! xxxy
      &  0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -s45_8, 0.0_wp,   s35_8,   0.0_wp, & ! xxxz
      &  0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,   d3_4, 0.0_wp, 0.0_wp,  0.0_wp, -s315_16, & ! xxyy
      &  0.0_wp, s315_8, 0.0_wp, -s45_8, 0.0_wp, 0.0_wp, 0.0_wp,  0.0_wp,   0.0_wp, & ! xxyz
      &  0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,-3.0_wp, 0.0_wp,  s45_4,  0.0_wp,   0.0_wp, & ! xxzz
      &  -s35_4, 0.0_wp,  -s5_4, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,  0.0_wp,   0.0_wp, & ! xyyy
      &  0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -s45_8, 0.0_wp, -s315_8,   0.0_wp, & ! xyyz
      &  0.0_wp, 0.0_wp,    s45, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,  0.0_wp,   0.0_wp, & ! xyzz
      &  0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,    s10, 0.0_wp,  0.0_wp,   0.0_wp, & ! xzzz
      &  0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,   d3_8, 0.0_wp,  s5_16,  0.0_wp,   s35_64, & ! yyyy
      &  0.0_wp, -s35_8, 0.0_wp, -s45_8, 0.0_wp, 0.0_wp, 0.0_wp,  0.0_wp,   0.0_wp, & ! yyyz
      &  0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,-3.0_wp, 0.0_wp, -s45_4,  0.0_wp,   0.0_wp, & ! yyzz
      &  0.0_wp, 0.0_wp, 0.0_wp,    s10, 0.0_wp, 0.0_wp, 0.0_wp,  0.0_wp,   0.0_wp, & ! yzzz
      &  0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp,  0.0_wp,   0.0_wp],& ! zzzz
      & shape(gtrafo))
   real(wp), parameter :: gtrafo_pinv(9, 15) = reshape([&
      ! -4       -3           -2        -1         0         1          2           3           4
      &  0.0_wp,      0.0_wp,   0.0_wp,    0.0_wp,   d9_370,    0.0_wp,   -s4_6845,     0.0_wp,   s289_191660, & ! xxxx
      &   s1_35,      0.0_wp, -s1_1805,    0.0_wp,   0.0_wp,    0.0_wp,     0.0_wp,     0.0_wp,        0.0_wp, & ! xxxy
      &  0.0_wp,      0.0_wp,   0.0_wp,    0.0_wp,   0.0_wp, -s81_3610,     0.0_wp, s169_25270,        0.0_wp, & ! xxxz
      &  0.0_wp,      0.0_wp,   0.0_wp,    0.0_wp,   d3_370,    0.0_wp,     0.0_wp,     0.0_wp, -s8649_191660, & ! xxyy
      &  0.0_wp,    s63_3610,   0.0_wp,  -s9_3610,   0.0_wp,    0.0_wp,     0.0_wp,     0.0_wp,        0.0_wp, & ! xxyz
      &  0.0_wp,      0.0_wp,   0.0_wp,    0.0_wp, -d57_370,    0.0_wp,  s144_6845,     0.0_wp,   -s81_191660, & ! xxzz
      &  -s1_35,      0.0_wp, -s1_1805,    0.0_wp,   0.0_wp,    0.0_wp,     0.0_wp,     0.0_wp,        0.0_wp, & ! xyyy
      &  0.0_wp,      0.0_wp,   0.0_wp,    0.0_wp,   0.0_wp,  -s9_3610,     0.0_wp,  -s63_3610,        0.0_wp, & ! xyyz
      &  0.0_wp,      0.0_wp, s36_1805,    0.0_wp,   0.0_wp,    0.0_wp,     0.0_wp,     0.0_wp,        0.0_wp, & ! xyzz
      &  0.0_wp,      0.0_wp,   0.0_wp,    0.0_wp,   0.0_wp,   s10_361,     0.0_wp, -s18_12635,        0.0_wp, & ! xzzz
      &  0.0_wp,      0.0_wp,   0.0_wp,    0.0_wp,   d9_370,    0.0_wp,    s4_6845,     0.0_wp,   s289_191660, & ! yyyy
      &  0.0_wp, -s169_25270,   0.0_wp, -s81_3610,   0.0_wp,    0.0_wp,     0.0_wp,     0.0_wp,        0.0_wp, & ! yyyz
      &  0.0_wp,      0.0_wp,   0.0_wp,    0.0_wp, -d57_370,    0.0_wp, -s144_6845,     0.0_wp,   -s81_191660, & ! yyzz
      &  0.0_wp,   s18_12635,   0.0_wp,   s10_361,   0.0_wp,    0.0_wp,     0.0_wp,     0.0_wp,        0.0_wp, & ! yzzz
      &  0.0_wp,      0.0_wp,   0.0_wp,    0.0_wp,  d19_370,    0.0_wp,     0.0_wp,     0.0_wp,     s9_191660],& ! zzzz
      & shape(gtrafo_pinv))

contains



!> Transformation from the cartesian to the spherical harmonic basis
!> for a shell pair block.
pure subroutine transform0(lj, li, cart, sphr, bra, ket)
   !> Angular momentum of ket shell i
   integer, intent(in) :: li
   !> Angular momentum of bra shell j
   integer, intent(in) :: lj
   !> Cartesian representation of the integral [bra j, ket i]
   real(wp), intent(in) :: cart(:, :)
   !> Spherical harmonic representation of the integral [bra j, ket i]
   real(wp), intent(out) :: sphr(:, :)
   !> Flag for transformation of the bra dimension
   logical, intent(in) :: bra
   !> Flag for transformation of the ket dimension
   logical, intent(in) :: ket
   if (.not. bra .and. .not. ket) then
      sphr = cart
      if (bra .and. lj == 1) sphr = sphr([2, 3, 1], :)
      if (ket .and. li == 1) sphr = sphr(:, [2, 3, 1])
      return
   end if

   if (.not. bra .and. ket) then
      select case(li)
      case(0, 1)
         sphr = cart
      case(2)
         sphr = matmul(cart, transpose(dtrafo))
      case(3)
         sphr = matmul(cart, transpose(ftrafo))
      case(4)
         sphr = matmul(cart, transpose(gtrafo))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select
      if (bra .and. lj == 1) sphr = sphr([2, 3, 1], :)
      if (ket .and. li == 1) sphr = sphr(:, [2, 3, 1])
      return
   end if

   if (bra .and. .not. ket) then
      select case(lj)
      case(0, 1)
         sphr = cart
      case(2)
         sphr = matmul(dtrafo, cart)
      case(3)
         sphr = matmul(ftrafo, cart)
      case(4)
         sphr = matmul(gtrafo, cart)
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select
      if (bra .and. lj == 1) sphr = sphr([2, 3, 1], :)
      if (ket .and. li == 1) sphr = sphr(:, [2, 3, 1])
      return
   end if

   ! Transform both dimensions
   select case(li)
   case(0, 1)
      select case(lj)
      case(0, 1)
         sphr = cart
      case(2)
         !sphr = matmul(dtrafo, cart)
         sphr(3, :) = cart(6, :) - 0.5_wp * (cart(1, :) + cart(4, :))
         sphr(4, :) = s3 * cart(3, :)
         sphr(2, :) = s3 * cart(5, :)
         sphr(5, :) = s3_4 * (cart(1, :) - cart(4, :))
         sphr(1, :) = s3 * cart(2, :)
      case(3)
         sphr = matmul(ftrafo, cart)
      case(4)
         sphr = matmul(gtrafo, cart)
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case(2)
      select case(lj)
      case(0, 1)
         !sphr = matmul(cart, transpose(dtrafo))
         sphr(:, 3) = cart(:, 6) - 0.5_wp * (cart(:, 1) + cart(:, 4))
         sphr(:, 4) = s3 * cart(:, 3)
         sphr(:, 2) = s3 * cart(:, 5)
         sphr(:, 5) = s3_4 * (cart(:, 1) - cart(:, 4))
         sphr(:, 1) = s3 * cart(:, 2)
      case(2)
         !sphr = matmul(dtrafo, matmul(cart, transpose(dtrafo)))
         sphr(3, 3) = cart(6, 6) &
            & - 0.5_wp * (cart(6, 1) + cart(6, 4) + cart(1, 6) + cart(4, 6)) &
            & + 0.25_wp * (cart(1, 1) + cart(1, 4) + cart(4, 1) + cart(4, 4))
         sphr([4, 2, 1], 3) = s3 * cart([3, 5, 2], 6) &
            & - s3_4 * (cart([3, 5, 2], 1) + cart([3, 5, 2], 4))
         sphr(5, 3) = s3_4 * (cart(1, 6) - cart(4, 6)) &
            & - s3 * 0.25_wp * (cart(1, 1) - cart(4, 1) + cart(1, 4) - cart(4, 4))
         sphr(3, 4) = s3 * cart(6, 3) - s3_4 * (cart(1, 3) + cart(4, 3))
         sphr([4, 2, 1], 4) = 3 * cart([3, 5, 2], 3)
         sphr(5, 4) = 1.5_wp * (cart(1, 3) - cart(4, 3))
         sphr(3, 2) = s3 * cart(6, 5) - s3_4 * (cart(1, 5) + cart(4, 5))
         sphr([4, 2, 1], 2) = 3 * cart([3, 5, 2], 5)
         sphr(5, 2) = 1.5_wp * (cart(1, 5) - cart(4, 5))
         sphr(3, 5) = s3_4 * (cart(6, 1) - cart(6, 4)) &
            & - s3 * 0.25_wp * (cart(1, 1) - cart(1, 4) + cart(4, 1) - cart(4, 4))
         sphr([4, 2, 1], 5) = 1.5_wp * (cart([3, 5, 2], 1) - cart([3, 5, 2], 4))
         sphr(5, 5) = 0.75_wp * (cart(1, 1) - cart(4, 1) - cart(1, 4) + cart(4, 4))
         sphr(3, 1) = s3 * cart(6, 2) - s3_4 * (cart(1, 2) + cart(4, 2))
         sphr([4, 2, 1], 1) = 3 * cart([3, 5, 2], 2)
         sphr(5, 1) = 1.5_wp * (cart(1, 2) - cart(4, 2))
      case(3)
         sphr = matmul(ftrafo, matmul(cart, transpose(dtrafo)))
      case(4)
         sphr = matmul(gtrafo, matmul(cart, transpose(dtrafo)))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case(3)
      select case(lj)
      case(0, 1)
         sphr = matmul(cart, transpose(ftrafo))
      case(2)
         sphr = matmul(dtrafo, matmul(cart, transpose(ftrafo)))
      case(3)
         sphr = matmul(ftrafo, matmul(cart, transpose(ftrafo)))
      case(4)
         sphr = matmul(gtrafo, matmul(cart, transpose(ftrafo)))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case(4)
      select case(lj)
      case(0, 1)
         sphr = matmul(cart, transpose(gtrafo))
      case(2)
         sphr = matmul(dtrafo, matmul(cart, transpose(gtrafo)))
      case(3)
         sphr = matmul(ftrafo, matmul(cart, transpose(gtrafo)))
      case(4)
         sphr = matmul(gtrafo, matmul(cart, transpose(gtrafo)))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case default
      error stop "[Fatal] Moments higher than g are not supported"
   end select

   if (bra .and. lj == 1) sphr = sphr([2, 3, 1], :)
   if (ket .and. li == 1) sphr = sphr(:, [2, 3, 1])

end subroutine transform0

!> Transformation from the cartesian to the spherical harmonic basis
!> for a vector of shell pair block.
pure subroutine transform1(lj, li, cart, sphr, bra, ket)
   !> Angular momentum of ket shell i
   integer, intent(in) :: li
   !> Angular momentum of bra shell j
   integer, intent(in) :: lj
   !> Cartesian representation of the integral [:, bra j, ket i]
   real(wp), intent(in) :: cart(:, :, :)
   !> Spherical harmonic representation of the integral [:, bra j, ket i]
   real(wp), intent(out) :: sphr(:, :, :)
   !> Flag for transformation of the bra dimension
   logical, intent(in) :: bra
   !> Flag for transformation of the ket dimension
   logical, intent(in) :: ket

   integer :: k

   do k = 1, size(cart, 1)
      call transform0(lj, li, cart(k, :, :), sphr(k, :, :), bra, ket)
   end do
end subroutine transform1

!> Transformation from the cartesian to the spherical harmonic basis
!> for a matrix of shell pair block.
pure subroutine transform2(lj, li, cart, sphr, bra, ket)
   !> Angular momentum of ket shell i
   integer, intent(in) :: li
   !> Angular momentum of bra shell j
   integer, intent(in) :: lj
   !> Cartesian representation of the integral [:, :, bra j, ket i]
   real(wp), intent(in) :: cart(:, :, :, :)
   !> Spherical harmonic representation of the integral [:, :, bra j, ket i]
   real(wp), intent(out) :: sphr(:, :, :, :)
   !> Flag for transformation of the bra dimension
   logical, intent(in) :: bra
   !> Flag for transformation of the ket dimension
   logical, intent(in) :: ket

   integer :: k, l

   do l = 1, size(cart, 2)
      do k = 1, size(cart, 1)
         call transform0(lj, li, cart(k, l, :, :), sphr(k, l, :, :), bra, ket)
      end do
   end do
end subroutine transform2


!> Adjoint transformation from the spherical harmonic to the cartesian basis
!> for a shell pair block. Applies quantities which behave contravariant w.r.t.
!> the basis functions (i.e. MO expansion coefficients or the density matrix)
pure subroutine adjoint_transform0(lj, li, sphr, cart, bra, ket)
   !> Angular momentum of ket shell i
   integer, intent(in) :: li
   !> Angular momentum of bra shell j
   integer, intent(in) :: lj
   !> Spherical harmonic representation of the integral [bra j, ket i]
   real(wp), intent(in) :: sphr(:, :)
   !> Cartesian representation of the integral [bra j, ket i]
   real(wp), intent(out) :: cart(:, :)
   !> Flag for transformation of the bra dimension
   logical, intent(in) :: bra
   !> Flag for transformation of the ket dimension
   logical, intent(in) :: ket
   if (.not. bra .and. .not. ket) then
      cart = sphr
      if (bra .and. lj == 1) cart = cart([3, 1, 2], :)
      if (ket .and. li == 1) cart = cart(:, [3, 1, 2])
      return
   end if

   ! Transform only ket dimension
   if (.not. bra .and. ket) then
      select case(li)
      case(0, 1)
         cart = sphr
      case(2)
         cart = matmul(sphr, dtrafo)
      case(3)
         cart = matmul(sphr, ftrafo)
      case(4)
         cart = matmul(sphr, gtrafo)
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select
      if (bra .and. lj == 1) cart = cart([3, 1, 2], :)
      if (ket .and. li == 1) cart = cart(:, [3, 1, 2])
      return
   end if

   ! Transform only bra dimension
   if (bra .and. .not. ket) then
      select case(lj)
      case(0, 1)
         cart = sphr
      case(2)
         cart = matmul(transpose(dtrafo), sphr)
      case(3)
         cart = matmul(transpose(ftrafo), sphr)
      case(4)
         cart = matmul(transpose(gtrafo), sphr)
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select
      if (bra .and. lj == 1) cart = cart([3, 1, 2], :)
      if (ket .and. li == 1) cart = cart(:, [3, 1, 2])
      return
   end if

   ! Transform both dimensions
   select case(lj)
   case(0, 1)
      select case(li)
      case(0, 1)
         cart = sphr
      case(2)
         cart = matmul(sphr, dtrafo)
      case(3)
         cart = matmul(sphr, ftrafo)
      case(4)
         cart = matmul(sphr, gtrafo)
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case(2)
      select case(li)
      case(0, 1)
         cart = matmul(transpose(dtrafo), sphr)
      case(2)
         cart = matmul(transpose(dtrafo), matmul(sphr, dtrafo))
      case(3)
         cart = matmul(transpose(dtrafo), matmul(sphr, ftrafo))
      case(4)
         cart = matmul(transpose(dtrafo), matmul(sphr, gtrafo))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case(3)
      select case(li)
      case(0, 1)
         cart = matmul(transpose(ftrafo), sphr)
      case(2)
         cart = matmul(transpose(ftrafo), matmul(sphr, dtrafo))
      case(3)
         cart = matmul(transpose(ftrafo), matmul(sphr, ftrafo))
      case(4)
         cart = matmul(transpose(ftrafo), matmul(sphr, gtrafo))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case(4)
      select case(li)
      case(0, 1)
         cart = matmul(transpose(gtrafo), sphr)
      case(2)
         cart = matmul(transpose(gtrafo), matmul(sphr, dtrafo))
      case(3)
         cart = matmul(transpose(gtrafo), matmul(sphr, ftrafo))
      case(4)
         cart = matmul(transpose(gtrafo), matmul(sphr, gtrafo))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case default
      error stop "[Fatal] Moments higher than g are not supported"
   end select
   if (bra .and. lj == 1) cart = cart([3, 1, 2], :)
   if (ket .and. li == 1) cart = cart(:, [3, 1, 2])

end subroutine adjoint_transform0

!> Adjoint transformation from the spherical harmonic to the cartesian basis
!> for a vector of shell pair block. Applies quantities which behave contravariant
!> w.r.t. the basis functions (i.e. MO expansion coefficients or the density matrix)
pure subroutine adjoint_transform1(lj, li, sphr, cart, bra, ket)
   !> Angular momentum of ket shell i
   integer, intent(in) :: li
   !> Angular momentum of bra shell j
   integer, intent(in) :: lj
   !> Spherical harmonic representation of the integral [:, bra j, ket i]
   real(wp), intent(in) :: sphr(:, :, :)
   !> Cartesian representation of the integral [:, bra j, ket i]
   real(wp), intent(out) :: cart(:, :, :)
   !> Flag for transformation of the bra dimension
   logical, intent(in) :: bra
   !> Flag for transformation of the ket dimension
   logical, intent(in) :: ket

   integer :: k

   do k = 1, size(sphr, 1)
      call adjoint_transform0(lj, li, sphr(k, :, :), cart(k, :, :), bra, ket)
   end do
end subroutine adjoint_transform1

!> Adjoint transformation from the spherical harmonic to the cartesian basis
!> for a matrix of shell pair block. Applies quantities which behave contravariant
!> w.r.t. the basis functions (i.e. MO expansion coefficients or the density matrix)
pure subroutine adjoint_transform2(lj, li, sphr, cart, bra, ket)
   !> Angular momentum of ket shell i
   integer, intent(in) :: li
   !> Angular momentum of bra shell j
   integer, intent(in) :: lj
   !> Spherical harmonic representation of the integral [:, :, bra j, ket i]
   real(wp), intent(in) :: sphr(:, :, :, :)
   !> Cartesian representation of the integral [:, :, bra j, ket i]
   real(wp), intent(out) :: cart(:, :, :, :)
   !> Flag for transformation of the bra dimension
   logical, intent(in) :: bra
   !> Flag for transformation of the ket dimension
   logical, intent(in) :: ket

   integer :: k, l

   do l = 1, size(sphr, 2)
      do k = 1, size(sphr, 1)
         call adjoint_transform0(lj, li, sphr(k, l, :, :), cart(k, l, :, :), bra, ket)
      end do
   end do
end subroutine adjoint_transform2


!> Inverse of the adjoint transformation, from the cartesian back to the spherical
!> harmonic basis for a shell pair block. Applies to quantities which behave contravariant
!> w.r.t. the basis functions (i.e. MO expansion coefficients or the density matrix)
pure subroutine contravariant_transform0(lj, li, cart, sphr, bra, ket)
   !> Angular momentum of ket shell i
   integer, intent(in) :: li
   !> Angular momentum of bra shell j
   integer, intent(in) :: lj
   !> Cartesian representation of the integral [bra j, ket i]
   real(wp), intent(in) :: cart(:, :)
   !> Spherical harmonic representation of the integral [bra j, ket i]
   real(wp), intent(out) :: sphr(:, :)
   !> Flag for transformation of the bra dimension
   logical, intent(in) :: bra
   !> Flag for transformation of the ket dimension
   logical, intent(in) :: ket
   if (.not. bra .and. .not. ket) then
      sphr = cart
      if (bra .and. lj == 1) sphr = sphr([2, 3, 1], :)
      if (ket .and. li == 1) sphr = sphr(:, [2, 3, 1])
      return
   end if

   ! Transform only ket dimension
   if (.not. bra .and. ket) then
      select case(li)
      case(0, 1)
         sphr = cart
      case(2)
         sphr = matmul(cart, transpose(dtrafo_pinv))
      case(3)
         sphr = matmul(cart, transpose(ftrafo_pinv))
      case(4)
         sphr = matmul(cart, transpose(gtrafo_pinv))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select
      if (bra .and. lj == 1) sphr = sphr([2, 3, 1], :)
      if (ket .and. li == 1) sphr = sphr(:, [2, 3, 1])
      return
   end if

   ! Transform only bra dimension
   if (bra .and. .not. ket) then
      select case(lj)
      case(0, 1)
         sphr = cart
      case(2)
         sphr = matmul(dtrafo_pinv, cart)
      case(3)
         sphr = matmul(ftrafo_pinv, cart)
      case(4)
         sphr = matmul(gtrafo_pinv, cart)
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select
      if (bra .and. lj == 1) sphr = sphr([2, 3, 1], :)
      if (ket .and. li == 1) sphr = sphr(:, [2, 3, 1])
      return
   end if

   ! Transform both dimensions
   select case(lj)
   case(0, 1)
      select case(li)
      case(0, 1)
         sphr = cart
      case(2)
         sphr = matmul(cart, transpose(dtrafo_pinv))
      case(3)
         sphr = matmul(cart, transpose(ftrafo_pinv))
      case(4)
         sphr = matmul(cart, transpose(gtrafo_pinv))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case(2)
      select case(li)
      case(0, 1)
         sphr = matmul(dtrafo_pinv, cart)
      case(2)
         sphr = matmul(dtrafo_pinv, matmul(cart, transpose(dtrafo_pinv)))
      case(3)
         sphr = matmul(dtrafo_pinv, matmul(cart, transpose(ftrafo_pinv)))
      case(4)
         sphr = matmul(dtrafo_pinv, matmul(cart, transpose(gtrafo_pinv)))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case(3)
      select case(li)
      case(0, 1)
         sphr = matmul(ftrafo_pinv, cart)
      case(2)
         sphr = matmul(ftrafo_pinv, matmul(cart, transpose(dtrafo_pinv)))
      case(3)
         sphr = matmul(ftrafo_pinv, matmul(cart, transpose(ftrafo_pinv)))
      case(4)
         sphr = matmul(ftrafo_pinv, matmul(cart, transpose(gtrafo_pinv)))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case(4)
      select case(li)
      case(0, 1)
         sphr = matmul(gtrafo_pinv, cart)
      case(2)
         sphr = matmul(gtrafo_pinv, matmul(cart, transpose(dtrafo_pinv)))
      case(3)
         sphr = matmul(gtrafo_pinv, matmul(cart, transpose(ftrafo_pinv)))
      case(4)
         sphr = matmul(gtrafo_pinv, matmul(cart, transpose(gtrafo_pinv)))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case default
      error stop "[Fatal] Moments higher than g are not supported"
   end select
   if (bra .and. lj == 1) sphr = sphr([2, 3, 1], :)
   if (ket .and. li == 1) sphr = sphr(:, [2, 3, 1])

end subroutine contravariant_transform0

!> Inverse of the adjoint transformation, from the cartesian back to the spherical
!> harmonic basis for a vector of shell pair block. Applies to quantities which behave
!> contravariant w.r.t. the basis functions (i.e. MO expansion coefficients or the
!> density matrix)
pure subroutine contravariant_transform1(lj, li, cart, sphr, bra, ket)
   !> Angular momentum of ket shell i
   integer, intent(in) :: li
   !> Angular momentum of bra shell j
   integer, intent(in) :: lj
   !> Cartesian representation of the integral [:, bra j, ket i]
   real(wp), intent(in) :: cart(:, :, :)
   !> Spherical harmonic representation of the integral [:, bra j, ket i]
   real(wp), intent(out) :: sphr(:, :, :)
   !> Flag for transformation of the bra dimension
   logical, intent(in) :: bra
   !> Flag for transformation of the ket dimension
   logical, intent(in) :: ket

   integer :: k

   do k = 1, size(cart, 1)
      call contravariant_transform0(lj, li, cart(k, :, :), sphr(k, :, :), bra, ket)
   end do
end subroutine contravariant_transform1

!> Inverse of the adjoint transformation, from the cartesian back to the spherical
!> harmonic basis for a matrix of shell pair block. Applies to quantities which behave
!> contravariant w.r.t. the basis functions (i.e. MO expansion coefficients or the
!> density matrix)
pure subroutine contravariant_transform2(lj, li, cart, sphr, bra, ket)
   !> Angular momentum of ket shell i
   integer, intent(in) :: li
   !> Angular momentum of bra shell j
   integer, intent(in) :: lj
   !> Cartesian representation of the integral [:, :, bra j, ket i]
   real(wp), intent(in) :: cart(:, :, :, :)
   !> Spherical harmonic representation of the integral [:, :, bra j, ket i]
   real(wp), intent(out) :: sphr(:, :, :, :)
   !> Flag for transformation of the bra dimension
   logical, intent(in) :: bra
   !> Flag for transformation of the ket dimension
   logical, intent(in) :: ket

   integer :: k, l

   do l = 1, size(cart, 2)
      do k = 1, size(cart, 1)
         call contravariant_transform0(lj, li, cart(k, l, :, :), sphr(k, l, :, :), bra, ket)
      end do
   end do
end subroutine contravariant_transform2

end module tblite_integral_trafo
