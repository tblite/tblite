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
   real(wp), parameter :: ref(16) = &
      &[4.11453332407988E+0_wp, 9.32058698646034E-1_wp, 2.03554558097438E+0_wp, &
      & 1.42227746710936E+0_wp, 1.12812389434909E+0_wp, 1.05491579723435E+0_wp, &
      & 1.52708946810529E+0_wp, 1.95070344717242E+0_wp, 3.83759768609002E+0_wp, &
      & 1.09388281571287E+0_wp, 1.07090747617095E+0_wp, 2.00285147985041E+0_wp, &
      & 4.36400613276303E+0_wp, 3.83469756749583E+0_wp, 3.91542430689253E+0_wp, &
      & 5.58571407311513E+0_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_cn_gen(error, mol, ref)

end subroutine test_cn_mb01


subroutine test_cn_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = &
      &[9.16720040324353E-1_wp, 3.81678274047125E+0_wp, 3.58669371363926E+0_wp, &
      & 2.90017307806127E+0_wp, 5.09896725493552E+0_wp, 1.13810099508394E+0_wp, &
      & 9.35735601902201E-1_wp, 9.43323828566168E-1_wp, 4.86704045877793E+0_wp, &
      & 1.22024345749220E+0_wp, 3.86661271986982E+0_wp, 4.02897680486562E+0_wp, &
      & 1.96869979721683E+0_wp, 9.04495111104734E-1_wp, 1.49977570995823E+0_wp, &
      & 2.04937594337046E+0_wp]

   call get_structure(mol, "MB16-43", "02")
   call test_cn_gen(error, mol, ref)

end subroutine test_cn_mb02


subroutine test_cn_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = &
      &[4.03003841929577E+0_wp, 2.88404791395922E+0_wp, 1.12030403562302E+0_wp, &
      & 4.89326350512164E+0_wp, 6.97176967437239E+0_wp, 4.43722359168413E+0_wp, &
      & 4.74866431687226E+0_wp, 1.30338515919807E+0_wp, 1.06130157537971E+0_wp, &
      & 1.04048346760957E+0_wp, 1.92270034504840E+0_wp, 2.98737782950765E+0_wp, &
      & 4.69081042627144E+0_wp, 1.26260066175568E+0_wp, 3.36367832701891E+0_wp, &
      & 1.35743816568803E+0_wp]

   call get_structure(mol, "MB16-43", "03")
   call test_cn_gen(error, mol, ref)

end subroutine test_cn_mb03


subroutine test_cn_acetic(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(32) = &
      &[1.11869274777858E+0_wp, 1.11873745311598E+0_wp, 1.11867486488261E+0_wp, &
      & 1.11871835301562E+0_wp, 1.00298041850431E+0_wp, 1.00296225934508E+0_wp, &
      & 1.00292025857609E+0_wp, 9.96672235171986E-1_wp, 9.96665264804772E-1_wp, &
      & 9.96644640634157E-1_wp, 9.96639566141031E-1_wp, 9.96009354820402E-1_wp, &
      & 9.96080944431250E-1_wp, 9.96044372561812E-1_wp, 9.96030201248580E-1_wp, &
      & 1.00297397087696E+0_wp, 3.11024982090751E+0_wp, 3.11027171938537E+0_wp, &
      & 3.11021637465367E+0_wp, 3.11025447500635E+0_wp, 4.00458274208821E+0_wp, &
      & 4.00464303907056E+0_wp, 4.00453773239264E+0_wp, 4.00457154999029E+0_wp, &
      & 2.09575496538156E+0_wp, 2.09584231844821E+0_wp, 2.09574552648371E+0_wp, &
      & 2.09582118878727E+0_wp, 1.31402587717708E+0_wp, 1.31407631240435E+0_wp, &
      & 1.31403992456289E+0_wp, 1.31407294078440E+0_wp]

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
