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

module test_ncoord_gfn
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, &
      & test_failed
   use mctc_io_structure, only : structure_type
   use mstore, only : get_structure
   use tblite_cutoff, only : get_lattice_points
   use tblite_data_covrad, only : get_covalent_rad
   use tblite_ncoord_gfn
   implicit none
   private

   public :: collect_ncoord_gfn

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))


contains


!> Collect all exported unit tests
subroutine collect_ncoord_gfn(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("cn-mb01", test_cn_mb01), &
      & new_unittest("cn-mb02", test_cn_mb02), &
      & new_unittest("cn-mb03", test_cn_mb03), &
      & new_unittest("cn-acetic", test_cn_acetic), &
      & new_unittest("dcndr-mb04", test_dcndr_mb04), &
      & new_unittest("dcndr-mb05", test_dcndr_mb05), &
      & new_unittest("dcndr-ammonia", test_dcndr_ammonia), &
      & new_unittest("dcndL-mb06", test_dcndL_mb06), &
      & new_unittest("dcndL-mb07", test_dcndL_mb07), &
      & new_unittest("dcndL-antracene", test_dcndL_anthracene) &
      & ]

end subroutine collect_ncoord_gfn


subroutine test_cn_gen(error, mol, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol

   !> Reference CNs
   real(wp), intent(in) :: ref(:)

   real(wp), allocatable :: cn(:), rcov(:)
   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), allocatable :: lattr(:, :)

   allocate(rcov(mol%nid), cn(mol%nat))
   rcov(:) = get_covalent_rad(mol%num)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call get_coordination_number(mol, lattr, cutoff, rcov, cn)

   if (any(abs(cn - ref) > thr)) then
      call test_failed(error, "Coordination numbers do not match")
      print'(3es21.14)', cn
   end if

end subroutine test_cn_gen


subroutine test_numgrad(error, mol)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   integer :: iat, ic
   real(wp), allocatable :: cn(:), rcov(:), cnr(:), cnl(:)
   real(wp), allocatable :: dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: numdr(:, :, :)
   real(wp), allocatable :: lattr(:, :)
   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(rcov(mol%nid), cn(mol%nat), cnr(mol%nat), cnl(mol%nat), &
      & dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & numdr(3, mol%nat, mol%nat))
   rcov(:) = get_covalent_rad(mol%num)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

   do iat = 1, mol%nat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_coordination_number(mol, lattr, cutoff, rcov, cnr)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call get_coordination_number(mol, lattr, cutoff, rcov, cnl)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numdr(ic, iat, :) = 0.5_wp*(cnr - cnl)/step
      end do
   end do

   call get_coordination_number(mol, lattr, cutoff, rcov, cn, dcndr, dcndL)

   if (any(abs(dcndr - numdr) > thr2)) then
      call test_failed(error, "Derivative of coordination number does not match")
   end if

end subroutine test_numgrad


subroutine test_numsigma(error, mol)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   integer :: ic, jc
   real(wp) :: eps(3, 3)
   real(wp), allocatable :: cn(:), rcov(:), cnr(:), cnl(:), xyz(:, :)
   real(wp), allocatable :: dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: numdL(:, :, :)
   real(wp), allocatable :: lattr(:, :), trans(:, :)
   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(rcov(mol%nid), cn(mol%nat), cnr(mol%nat), cnl(mol%nat), &
      & dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), xyz(3, mol%nat), &
      & numdL(3, 3, mol%nat))
   rcov(:) = get_covalent_rad(mol%num)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   trans = lattr
   do ic = 1, 3
      do jc = 1, 3
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         lattr(:, :) = matmul(eps, trans)
         call get_coordination_number(mol, lattr, cutoff, rcov, cnr)
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         lattr(:, :) = matmul(eps, trans)
         call get_coordination_number(mol, lattr, cutoff, rcov, cnl)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         lattr(:, :) = trans
         numdL(jc, ic, :) = 0.5_wp*(cnr - cnl)/step
      end do
   end do

   call get_coordination_number(mol, lattr, cutoff, rcov, cn, dcndr, dcndL)

   if (any(abs(dcndL - numdL) > thr2)) then
      call test_failed(error, "Derivative of coordination number does not match")
   end if

end subroutine test_numsigma


subroutine test_cn_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = [&
      & 4.11453659059991E+0_wp, 9.32058998762811E-1_wp, 2.03554597140311E+0_wp, &
      & 1.42227835389358E+0_wp, 1.12812426574031E+0_wp, 1.05491602558828E+0_wp, &
      & 1.52709064704269E+0_wp, 1.95070367247232E+0_wp, 3.83759889196540E+0_wp, &
      & 1.09388314007182E+0_wp, 1.07090773695340E+0_wp, 2.00285254082830E+0_wp, &
      & 4.36400837813955E+0_wp, 3.83469860546080E+0_wp, 3.91542517673963E+0_wp, &
      & 5.58571682419960E+0_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_cn_gen(error, mol, ref)

end subroutine test_cn_mb01


subroutine test_cn_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = [&
      & 9.16720394154780E-1_wp, 3.81678481096046E+0_wp, 3.58669504034282E+0_wp, &
      & 2.90017371823910E+0_wp, 5.09896950232009E+0_wp, 1.13810142820063E+0_wp, &
      & 9.35735928561309E-1_wp, 9.43324181507483E-1_wp, 4.86704300506033E+0_wp, &
      & 1.22024405676881E+0_wp, 3.86661472534511E+0_wp, 4.02897853934353E+0_wp, &
      & 1.96869999324996E+0_wp, 9.04495366214972E-1_wp, 1.49977675358986E+0_wp, &
      & 2.04937665719480E+0_wp]

   call get_structure(mol, "MB16-43", "02")
   call test_cn_gen(error, mol, ref)

end subroutine test_cn_mb02


subroutine test_cn_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = [&
      & 4.03003984525826E+0_wp, 2.88404829383272E+0_wp, 1.12030447706437E+0_wp, &
      & 4.89326653038565E+0_wp, 6.97177319120287E+0_wp, 4.43722632690247E+0_wp, &
      & 4.74866591812568E+0_wp, 1.30338583081493E+0_wp, 1.06130190313902E+0_wp, &
      & 1.04048373066340E+0_wp, 1.92270060977307E+0_wp, 2.98737904327655E+0_wp, &
      & 4.69081207436918E+0_wp, 1.26260125612717E+0_wp, 3.36368006785498E+0_wp, &
      & 1.35743886389476E+0_wp]

   call get_structure(mol, "MB16-43", "03")
   call test_cn_gen(error, mol, ref)

end subroutine test_cn_mb03


subroutine test_cn_acetic(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(32) = [&
      & 1.11869327906258E+0_wp, 1.11873798419384E+0_wp, 1.11867539605108E+0_wp, &
      & 1.11871888406156E+0_wp, 1.00298061967398E+0_wp, 1.00296246055984E+0_wp, &
      & 1.00292045987951E+0_wp, 9.96672431444661E-1_wp, 9.96665461112044E-1_wp, &
      & 9.96644836940640E-1_wp, 9.96639762485971E-1_wp, 9.96009549496349E-1_wp, &
      & 9.96081138973917E-1_wp, 9.96044567123842E-1_wp, 9.96030395875155E-1_wp, &
      & 1.00297417207048E+0_wp, 3.11025027428883E+0_wp, 3.11027217276691E+0_wp, &
      & 3.11021682807939E+0_wp, 3.11025492840344E+0_wp, 4.00458350822475E+0_wp, &
      & 4.00464380514876E+0_wp, 4.00453849863709E+0_wp, 4.00457231615400E+0_wp, &
      & 2.09575556191806E+0_wp, 2.09584291489233E+0_wp, 2.09574612300023E+0_wp, &
      & 2.09582178519930E+0_wp, 1.31402642902757E+0_wp, 1.31407686428049E+0_wp, &
      & 1.31404047639461E+0_wp, 1.31407349269578E+0_wp]

   call get_structure(mol, "X23", "acetic")
   call test_cn_gen(error, mol, ref)

end subroutine test_cn_acetic


subroutine test_dcndr_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "04")
   call test_numgrad(error, mol)

end subroutine test_dcndr_mb04


subroutine test_dcndr_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_numgrad(error, mol)

end subroutine test_dcndr_mb05


subroutine test_dcndr_ammonia(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "ammonia")
   call test_numgrad(error, mol)

end subroutine test_dcndr_ammonia


subroutine test_dcndL_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "06")
   call test_numsigma(error, mol)

end subroutine test_dcndL_mb06


subroutine test_dcndL_mb07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "07")
   call test_numsigma(error, mol)

end subroutine test_dcndL_mb07


subroutine test_dcndL_anthracene(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "anthracene")
   call test_numsigma(error, mol)

end subroutine test_dcndL_anthracene


end module test_ncoord_gfn
