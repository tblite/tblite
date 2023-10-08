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

module test_ceh
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
   & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_context_type, only : context_type
   use tblite_lapack_solver, only : lapack_solver, lapack_algorithm
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_ceh_calculator, only : ceh_calculator
   use tblite_ceh_ceh, only: ceh_guess, new_ceh_calculator

   use tblite_blas, only: gemv

   use tblite_container, only : container_type, container_cache
   use tblite_external_field, only : electric_field
   implicit none
   private

   public :: collect_ceh

   real(wp), parameter :: kt = 300.0_wp * 3.166808578545117e-06_wp

contains

   !> Collect all exported unit tests
   subroutine collect_ceh(testsuite)

      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("cn-mol", test_cn_mb12), &
         new_unittest("q-mol", test_q_mb01), &
         new_unittest("q-chrgd-efield-mol", test_q_ef_chrg_mb01), &
         new_unittest("d-mol", test_d_mb01), &
         new_unittest("d-field-mol", test_d_mb04), &
         new_unittest("d-field-change-mol", test_d_hcn) &
         ]

   end subroutine collect_ceh

   subroutine test_cn_mb12(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(structure_type) :: mol
      type(ceh_calculator) :: calc
      type(wavefunction_type) :: wfn
      real(wp), allocatable :: cn(:), cn_en(:)
      real(wp), parameter :: ref1(16) = reshape([ & ! calculated with GP3 standalone
      1.31825913_wp, &
      1.47303090_wp, &
      0.46520171_wp, &
      1.01541807_wp, &
      0.48417989_wp, &
      1.00793531_wp, &
      0.42618861_wp, &
      0.23257833_wp, &
      0.73395771_wp, &
      1.02341545_wp, &
      1.22594362_wp, &
      0.26119136_wp, &
      1.55186287_wp, &
      1.40171890_wp, &
      0.44331842_wp, &
      0.50585021_wp], shape(ref1))
      real(wp), parameter :: ref2(16) = reshape([ & ! calculated with GP3 standalone
       0.30785024_wp, &
      -0.02792822_wp, &
      -0.11358521_wp, &
       0.03271103_wp, &
       0.10217515_wp, &
       0.25825112_wp, &
      -0.20769464_wp, &
      -0.03447769_wp, &
      -0.25817610_wp, &
       0.24460001_wp, &
       0.26974150_wp, &
      -0.01050015_wp, &
      -0.28007123_wp, &
      -0.42862573_wp, &
       0.03898681_wp, &
       0.10674311_wp], shape(ref2))
      integer :: i

      call get_structure(mol, "MB16-43", "12")
      call new_ceh_calculator(calc, mol)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      ctx%verbosity = 0
      call ceh_guess(ctx, calc, mol, error, wfn)
      allocate(cn(mol%nat), cn_en(mol%nat))
      call calc%ncoordstd%get_cn(mol, cn)
      call calc%ncoorden%get_cn(mol, cn_en)
      do i = 1, mol%nat
         call check(error, cn(i), ref1(i), thr=1e-7_wp)
         if (allocated(error)) return
         call check(error, cn_en(i), ref2(i), thr=1e-7_wp)
         if (allocated(error)) return
      enddo

   end subroutine test_cn_mb12

   subroutine test_q_mb01(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(structure_type) :: mol
      type(ceh_calculator) :: calc
      type(wavefunction_type) :: wfn
      real(wp), parameter :: ref(16) = reshape([ & ! calculated with GP3 standalone 
       0.5041712306_wp, & 
      -0.0768741000_wp, & 
      -0.4935157669_wp, & 
      -0.0831876027_wp, & 
      -0.2122917586_wp, & 
       0.1274119295_wp, & 
      -0.0434563264_wp, & 
      -0.3788163344_wp, & 
      -0.3016588466_wp, & 
       0.1576514268_wp, & 
       0.1353213766_wp, & 
       0.0150687156_wp, & 
       0.0511522155_wp, & 
       0.1399127014_wp, & 
      -0.0749090701_wp, & 
       0.5340202097_wp], shape(ref))
      integer :: i

      call get_structure(mol, "MB16-43", "01")
      call new_ceh_calculator(calc, mol)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      ctx%verbosity = 0
      call ceh_guess(ctx, calc, mol, error, wfn)
      do i = 1, mol%nat
         call check(error, wfn%qat(i,1), ref(i), thr=1e-6_wp)
         if (allocated(error)) return
      enddo

   end subroutine test_q_mb01

   subroutine test_q_ef_chrg_mb01(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(structure_type) :: mol
      type(ceh_calculator) :: calc
      type(wavefunction_type) :: wfn
      class(container_type), allocatable :: cont
      real(wp), parameter :: ref(16) = reshape([ & ! calculated with GP3 standalone
      -6.1090763982_wp, &
      -0.9999265530_wp, &
       3.7865813535_wp, &
      -0.9753348768_wp, &
       6.9999970090_wp, &
       0.6329394986_wp, &
      -0.6638462346_wp, &
       5.2331540651_wp, &
       1.3850165624_wp, &
       0.7793708797_wp, &
       0.9967264400_wp, &
      -10.7642634986_wp, &
      -4.9181285308_wp, &
       2.4379171576_wp, &
       4.2321691647_wp, &
      -0.0532960387_wp &
       & ], shape(ref))
      real(wp) :: efield(3)
      integer :: i

      efield = 0.0_wp
      efield(3) = 0.2_wp

      call get_structure(mol, "MB16-43", "01")
      mol%charge = 2.0_wp
      call new_ceh_calculator(calc, mol)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      cont = electric_field(efield)
      call calc%push_back(cont)
      ctx%verbosity = 0
      call ceh_guess(ctx, calc, mol, error, wfn)
      do i = 1, mol%nat
         call check(error, wfn%qat(i,1), ref(i), thr=5e-6_wp, message="Calculated charge& 
         & does not match reference")
         if (allocated(error)) return
      enddo

   end subroutine test_q_ef_chrg_mb01

   subroutine test_d_mb01(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(structure_type) :: mol
      type(ceh_calculator) :: calc
      type(wavefunction_type) :: wfn
      real(wp) :: dipole(3), tmp(3)
      real(wp), parameter :: ref(3) = reshape([ & ! calculated with GP3 standalone
         0.584361099036660_wp, &
         -1.47304239280996_wp, &
         -2.25861915370679_wp], shape(ref))

      call get_structure(mol, "MB16-43", "01")

      call new_ceh_calculator(calc, mol)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      ctx%verbosity = 0
      call ceh_guess(ctx, calc, mol, error, wfn)
      tmp = 0.0_wp
      dipole = 0.0_wp
      call gemv(mol%xyz, wfn%qat(:, 1), tmp)
      dipole(:) = tmp + sum(wfn%dpat(:, :, 1), 2)
      if (any(abs(dipole - ref) > 1e-5_wp)) then
         call test_failed(error, "Numerical dipole moment does not match")
         print '(3es21.14)', dipole
         print '("---")'
         print '(3es21.14)', ref 
         print '("---")'
         print '(3es21.14)', dipole - ref 
      end if

   end subroutine test_d_mb01

   subroutine test_d_mb04(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(structure_type) :: mol
      type(ceh_calculator) :: calc
      type(wavefunction_type) :: wfn
      class(container_type), allocatable :: cont
      real(wp) :: energy, efield(3), dipole(3), tmp(3)
      real(wp), parameter :: ref(3) = reshape([ & ! calculated with GP3 standalone
         -16.4396031161495_wp, &
         90.2215123832578_wp, &
         -8.00262461340548_wp], shape(ref))

      call get_structure(mol, "MB16-43", "04")
      energy = 0.0_wp
      efield = 0.0_wp
      efield(2) = 0.2_wp

      call new_ceh_calculator(calc, mol)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)

      cont = electric_field(efield)
      call calc%push_back(cont)

      ctx%verbosity = 0
      call ceh_guess(ctx, calc, mol, error, wfn)
      tmp = 0.0_wp
      dipole = 0.0_wp
      call gemv(mol%xyz, wfn%qat(:, 1), tmp)
      dipole(:) = tmp + sum(wfn%dpat(:, :, 1), 2)

      if (any(abs(dipole - ref) > 1e-5_wp)) then
         call test_failed(error, "Numerical dipole moment does not match")
         print '(3es21.14)', dipole
         print '("---")'
         print '(3es21.14)', ref 
         print '("---")'
         print '(3es21.14)', dipole - ref
      end if

   end subroutine test_d_mb04

   subroutine test_d_hcn(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(structure_type) :: mol1,mol2
      type(ceh_calculator) :: calc1,calc2
      type(wavefunction_type) :: wfn1,wfn2
      class(container_type), allocatable :: cont1,cont2
      real(wp) :: efield(3), dip1(3), dip2(3), tmp(3)
      integer, parameter :: num(3) = reshape([ &
         7, &
         6, &
         1], shape(num))
      integer, parameter :: nat = 3
      real(wp) :: xyz(3, nat) = reshape([ &
      & -0.09604091224796_wp,  0.0_wp, 0.0_wp, &
      &  2.09604091224796_wp,  0.0_wp, 0.0_wp, &
      &  4.10859879422050_wp,  0.0_wp, 0.0_wp], &
      & shape(xyz))

      ctx%verbosity = 0
      call new(mol1, num, xyz) 
      efield = 0.0_wp
      efield(1) = -0.1_wp
      call new_ceh_calculator(calc1, mol1)
      call new_wavefunction(wfn1, mol1%nat, calc1%bas%nsh, calc1%bas%nao, 1, kt)
      cont1 = electric_field(efield)
      call calc1%push_back(cont1)
      call ceh_guess(ctx, calc1, mol1, error, wfn1)
      tmp = 0.0_wp
      dip1 = 0.0_wp
      call gemv(mol1%xyz, wfn1%qat(:, 1), tmp)
      dip1(:) = tmp + sum(wfn1%dpat(:, :, 1), 2)

      xyz(1, :) = xyz(1, :) - 1.0_wp
      call new(mol2, num, xyz) 
      call new_ceh_calculator(calc2, mol2)
      call new_wavefunction(wfn2, mol2%nat, calc2%bas%nsh, calc2%bas%nao, 1, kt)
      cont2 = electric_field(efield)
      call calc2%push_back(cont2)
      call ceh_guess(ctx, calc2, mol2, error, wfn2)
      tmp = 0.0_wp
      dip2 = 0.0_wp
      call gemv(mol2%xyz, wfn2%qat(:, 1), tmp)
      dip2(:) = tmp + sum(wfn2%dpat(:, :, 1), 2)

      if (any(abs(dip1 - dip2) > 1e-7_wp)) then
         call test_failed(error, "Numerical dipole moment does not match")
         print '(3es21.14)', dip1
         print '("---")'
         print '(3es21.14)', dip2
         print '("---")'
         print '(3es21.14)', dip1 - dip2
      end if

   end subroutine test_d_hcn
end module test_ceh
