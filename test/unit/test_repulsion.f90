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

module test_repulsion
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_container_cache, only : container_cache
   use tblite_repulsion, only : tb_repulsion, new_repulsion
   implicit none
   private

   public :: collect_repulsion

   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   abstract interface
      subroutine repulsion_maker(rep, mol)
         import :: tb_repulsion, structure_type
         type(tb_repulsion), intent(out) :: rep
         type(structure_type), intent(in) :: mol
      end subroutine repulsion_maker
   end interface

contains


!> Collect all exported unit tests
subroutine collect_repulsion(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("energy-1", test_e_effective_m01), &
      new_unittest("energy-2", test_e_effective_m02), &
      new_unittest("energy-pbc", test_e_effective_uracil), &
      new_unittest("gradient-1", test_g_effective_m03), &
      new_unittest("gradient-2", test_g_effective_m04), &
      new_unittest("gradient-pbc", test_g_effective_urea), &
      new_unittest("sigma-1", test_s_effective_m05), &
      new_unittest("sigma-2", test_s_effective_m06), &
      new_unittest("sigma-pbc", test_s_effective_succinic) &
      ]

end subroutine collect_repulsion


!> Factory to create repulsion objects based on GFN1-xTB
subroutine make_repulsion1(rep, mol)
   type(tb_repulsion), intent(out) :: rep
   type(structure_type), intent(in) :: mol

   real(wp), parameter :: alpha_gfn1(1:20) = [&
      & 2.209700_wp, 1.382907_wp, 0.671797_wp, 0.865377_wp, 1.093544_wp, &
      & 1.281954_wp, 1.727773_wp, 2.004253_wp, 2.507078_wp, 3.038727_wp, &
      & 0.704472_wp, 0.862629_wp, 0.929219_wp, 0.948165_wp, 1.067197_wp, &
      & 1.200803_wp, 1.404155_wp, 1.323756_wp, 0.581529_wp, 0.665588_wp]
   real(wp), parameter :: zeff_gfn1(1:20) = [&
      &  1.116244_wp,  0.440231_wp,  2.747587_wp,  4.076830_wp,  4.458376_wp, &
      &  4.428763_wp,  5.498808_wp,  5.171786_wp,  6.931741_wp,  9.102523_wp, &
      & 10.591259_wp, 15.238107_wp, 16.283595_wp, 16.898359_wp, 15.249559_wp, &
      & 15.100323_wp, 17.000000_wp, 17.153132_wp, 20.831436_wp, 19.840212_wp]
   real(wp), allocatable :: alpha(:), zeff(:), kexp

   alpha = alpha_gfn1(mol%num)
   zeff = zeff_gfn1(mol%num)
   kexp = 1.5_wp
   call new_repulsion(rep, mol, alpha, zeff, 1.5_wp, kexp, 1.0_wp)

end subroutine make_repulsion1

!> Factory to create repulsion objects based on GFN2-xTB
subroutine make_repulsion2(rep, mol)
   type(tb_repulsion), intent(out) :: rep
   type(structure_type), intent(in) :: mol

   real(wp), parameter :: alpha_gfn2(1:20) = [&
      & 2.213717_wp, 3.604670_wp, 0.475307_wp, 0.939696_wp, 1.373856_wp, &
      & 1.247655_wp, 1.682689_wp, 2.165712_wp, 2.421394_wp, 3.318479_wp, &
      & 0.572728_wp, 0.917975_wp, 0.876623_wp, 1.187323_wp, 1.143343_wp, &
      & 1.214553_wp, 1.577144_wp, 0.896198_wp, 0.482206_wp, 0.683051_wp]
   real(wp), parameter :: zeff_gfn2(1:20) = [&
      &  1.105388_wp,  1.094283_wp,  1.289367_wp,  4.221216_wp,  7.192431_wp, &
      &  4.231078_wp,  5.242592_wp,  5.784415_wp,  7.021486_wp, 11.041068_wp, &
      &  5.244917_wp, 18.083164_wp, 17.867328_wp, 40.001111_wp, 19.683502_wp, &
      & 14.995090_wp, 17.353134_wp,  7.266606_wp, 10.439482_wp, 14.786701_wp]
   real(wp), allocatable :: alpha(:), zeff(:), kexp

   alpha = alpha_gfn2(mol%num)
   zeff = zeff_gfn2(mol%num)
   kexp = 1.0_wp
   call new_repulsion(rep, mol, alpha, zeff, 1.5_wp, kexp, 1.0_wp)

end subroutine make_repulsion2


subroutine test_generic(error, mol, make_repulsion, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Factory to produce repulsion objects
   procedure(repulsion_maker) :: make_repulsion

   !> Reference value to compare against
   real(wp), intent(in) :: ref

   type(tb_repulsion) :: rep
   type(container_cache) :: cache
   real(wp) :: energy(mol%nat), sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), numgrad(:, :)
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   call make_repulsion(rep, mol)
   call rep%get_engrad(mol, cache, energy)

   call check(error, sum(energy), ref, thr=thr)
   if (allocated(error)) then
      print*,sum(energy)
   end if

end subroutine test_generic


subroutine test_numgrad(error, mol, make_repulsion)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Factory to create a repulsion object
   procedure(repulsion_maker) :: make_repulsion

   integer :: iat, ic
   type(tb_repulsion) :: rep
   type(container_cache) :: cache
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), numgrad(:, :)
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   call make_repulsion(rep, mol)

   do iat = 1, mol%nat
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call rep%get_engrad(mol, cache, er)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call rep%get_engrad(mol, cache, el)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   call rep%get_engrad(mol, cache, energy, gradient, sigma)

   if (any(abs(gradient - numgrad) > thr2)) then
      call test_failed(error, "Gradient of dispersion energy does not match")
      print'(3es21.14)', gradient-numgrad
   end if

end subroutine test_numgrad


subroutine test_numsigma(error, mol, make_repulsion)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Factory to create a repulsion object
   procedure(repulsion_maker) :: make_repulsion

   integer :: ic, jc
   type(tb_repulsion) :: rep
   type(container_cache) :: cache
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3), eps(3, 3), numsigma(3, 3)
   real(wp), allocatable :: gradient(:, :), xyz(:, :), lat(:, :)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(gradient(3, mol%nat), xyz(3, mol%nat), lat(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   call make_repulsion(rep, mol)

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   if (any(mol%periodic)) &
   lat(:, :) = mol%lattice
   do ic = 1, 3
      do jc = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (any(mol%periodic)) &
         mol%lattice(:, :) = matmul(eps, lat)
         call rep%get_engrad(mol, cache, er)
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (any(mol%periodic)) &
         mol%lattice(:, :) = matmul(eps, lat)
         call rep%get_engrad(mol, cache, el)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         if (any(mol%periodic)) &
         mol%lattice(:, :) = lat
         numsigma(jc, ic) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   call rep%get_engrad(mol, cache, energy, gradient, sigma)

   if (any(abs(sigma - numsigma) > thr2)) then
      call test_failed(error, "Strain derivatives do not match")
      print'(3es21.14)', sigma-numsigma
   end if

end subroutine test_numsigma


subroutine test_e_effective_m01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call test_generic(error, mol, make_repulsion1, 0.16777923624986593_wp)

end subroutine test_e_effective_m01


subroutine test_e_effective_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "02")
   call test_generic(error, mol, make_repulsion2, 0.10745931926703985_wp)

end subroutine test_e_effective_m02


subroutine test_e_effective_uracil(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "uracil")
   call test_generic(error, mol, make_repulsion2, 1.0401472262740301_wp)

end subroutine test_e_effective_uracil


subroutine test_g_effective_m03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "03")
   call test_numgrad(error, mol, make_repulsion1)

end subroutine test_g_effective_m03


subroutine test_g_effective_m04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "04")
   call test_numgrad(error, mol, make_repulsion2)

end subroutine test_g_effective_m04


subroutine test_g_effective_urea(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "urea")
   call test_numgrad(error, mol, make_repulsion2)

end subroutine test_g_effective_urea


subroutine test_s_effective_m05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_numsigma(error, mol, make_repulsion1)

end subroutine test_s_effective_m05


subroutine test_s_effective_m06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "06")
   call test_numsigma(error, mol, make_repulsion2)

end subroutine test_s_effective_m06


subroutine test_s_effective_succinic(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "succinic")
   call test_numsigma(error, mol, make_repulsion2)

end subroutine test_s_effective_succinic


end module test_repulsion
