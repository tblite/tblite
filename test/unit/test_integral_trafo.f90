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

module test_integral_trafo
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_basis_type
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_cutoff, only : get_lattice_points
   use tblite_integral_overlap
   use tblite_integral_trafo, only : transform0, transform1, transform2, &
      & adjoint_transform0, adjoint_transform1, adjoint_transform2

   implicit none
   private

   public :: collect_integral_trafo

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains


!> Collect all exported unit tests
subroutine collect_integral_trafo(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("trafo-ss", test_trafo_ss), &
      new_unittest("trafo-sp", test_trafo_sp), &
      new_unittest("trafo-sd", test_trafo_sd), &
      new_unittest("trafo-sf", test_trafo_sf), &
      new_unittest("trafo-sg", test_trafo_sg), &
      new_unittest("trafo-pp", test_trafo_pp), &
      new_unittest("trafo-pd", test_trafo_pd), &
      new_unittest("trafo-pf", test_trafo_pf), &
      new_unittest("trafo-pg", test_trafo_pg), &
      new_unittest("trafo-dd", test_trafo_dd), &
      new_unittest("trafo-df", test_trafo_df), &
      new_unittest("trafo-dg", test_trafo_dg), &
      new_unittest("trafo-ff", test_trafo_ff), &
      new_unittest("trafo-fg", test_trafo_fg), &
      new_unittest("trafo-gg", test_trafo_gg) &
      ]

end subroutine collect_integral_trafo


subroutine test_trafo_pair(error, lj, li)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Angular momentum of bra
   integer, intent(in) :: lj
   !> Angular momentum of ket
   integer, intent(in) :: li

   ! Transform neither bra nor ket
   call test_trafo_adjoint(error, lj, li, .false., .false.)
   if (allocated(error)) return

   ! Transform bra only
   call test_trafo_adjoint(error, lj, li, .true., .false.)
   if (allocated(error)) return

   ! Transform ket only
   call test_trafo_adjoint(error, lj, li, .false., .true.)
   if (allocated(error)) return

   ! Transform both bra and ket
   call test_trafo_adjoint(error, lj, li, .true., .true.)
   if (allocated(error)) return

end subroutine test_trafo_pair

subroutine test_trafo_adjoint(error, lj, li, bra, ket)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Angular momentum of bra
   integer, intent(in) :: lj
   !> Angular momentum of ket
   integer, intent(in) :: li
   !> Whether to transform the bra
   logical, intent(in) :: bra
   !> Whether to transform the ket
   logical, intent(in) :: ket

   real(wp), allocatable :: cart1(:, :), cart2(:, :)
   real(wp), allocatable :: sphr1(:, :), sphr2(:, :)
   real(wp) :: prod_sphr, prod_cart
   integer :: ncj, nci, nsj, nsi, nrow_cart, ncol_cart

   ! Number of Cartesian functions for bra and ket
   ncj = (lj + 1) * (lj + 2) / 2
   nci = (li + 1) * (li + 2) / 2

   ! Number of spherical functions for bra and ket
   nsj = 2*lj + 1
   nsi = 2*li + 1

   if (bra) then
      nrow_cart = ncj
   else
      nrow_cart = nsj
   end if

   if (ket) then
      ncol_cart = nci
   else
      ncol_cart = nsi
   end if

   allocate(cart1(nrow_cart, ncol_cart), cart2(nrow_cart, ncol_cart))
   allocate(sphr1(nsj, nsi), sphr2(nsj, nsi))

   ! Random coefficients around around 0.0
   call random_number(cart1)
   call random_number(sphr1)
   cart1 = cart1 - 0.5_wp
   sphr1 = sphr1 - 0.5_wp

   ! Transformation from cartesian to spherical
   if (bra .and. ket) then
      call transform0(lj, li, cart1, sphr2)

   else if (bra .and. .not. ket) then
      call transform0(lj, 0, cart1, sphr2)

   else if (.not. bra .and. ket) then
      call transform0(0, li, cart1, sphr2)

   else
      sphr2 = cart1
   end if

   ! Adjoint transformation from spherical to cartesian
   call adjoint_transform0(lj, li, sphr1, cart2, bra, ket)

   ! Hadamard product in spherical and cartesian basis
   prod_sphr = sum(sphr1 * sphr2)
   prod_cart = sum(cart1 * cart2)

   call check(error, prod_sphr, prod_cart, thr=thr)

end subroutine test_trafo_adjoint


subroutine test_trafo_ss(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 0, 0)

end subroutine test_trafo_ss

subroutine test_trafo_sp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 0, 1)
   call test_trafo_pair(error, 1, 0)

end subroutine test_trafo_sp

subroutine test_trafo_sd(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 0, 2)
   call test_trafo_pair(error, 2, 0)

end subroutine test_trafo_sd

subroutine test_trafo_sf(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 0, 3)
   call test_trafo_pair(error, 3, 0)

end subroutine test_trafo_sf

subroutine test_trafo_sg(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 0, 4)
   call test_trafo_pair(error, 4, 0)

end subroutine test_trafo_sg

subroutine test_trafo_pp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 1, 1)

end subroutine test_trafo_pp

subroutine test_trafo_pd(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 1, 2)
   call test_trafo_pair(error, 2, 1)

end subroutine test_trafo_pd

subroutine test_trafo_pf(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 1, 3)
   call test_trafo_pair(error, 3, 1)

end subroutine test_trafo_pf

subroutine test_trafo_pg(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 1, 4)
   call test_trafo_pair(error, 4, 1)

end subroutine test_trafo_pg

subroutine test_trafo_dd(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 2, 2)

end subroutine test_trafo_dd

subroutine test_trafo_df(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 2, 3)
   call test_trafo_pair(error, 3, 2)

end subroutine test_trafo_df

subroutine test_trafo_dg(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 2, 4)
   call test_trafo_pair(error, 4, 2)

end subroutine test_trafo_dg

subroutine test_trafo_ff(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 3, 3)

end subroutine test_trafo_ff

subroutine test_trafo_fg(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 3, 4)
   call test_trafo_pair(error, 4, 3)

end subroutine test_trafo_fg

subroutine test_trafo_gg(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 4, 4)

end subroutine test_trafo_gg

end module test_integral_trafo
