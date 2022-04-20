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

module test_integral_overlap
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_basis_type
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_cutoff, only : get_lattice_points
   use tblite_integral_overlap
   implicit none
   private

   public :: collect_integral_overlap

   real(wp), parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains


!> Collect all exported unit tests
subroutine collect_integral_overlap(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("overlap-1", test_overlap_alh3), &
      new_unittest("overlap-2", test_overlap_bh3), &
      new_unittest("overlap-3", test_overlap_beh2), &
      new_unittest("overlap-4", test_overlap_ch4), &
      new_unittest("overlap-5", test_overlap_cl2), &
      new_unittest("overlap-6", test_overlap_f2), &
      new_unittest("overlap-7", test_overlap_h2), &
      new_unittest("overlap-8", test_overlap_lih), &
      new_unittest("overlap-grad-ss", test_overlap_grad_ss), &
      new_unittest("overlap-grad-pp", test_overlap_grad_pp), &
      new_unittest("overlap-grad-dd", test_overlap_grad_dd) &
      ]

end subroutine collect_integral_overlap


subroutine make_basis(bas, mol, ng)
   type(basis_type), intent(out) :: bas
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: ng

   integer, parameter :: nsh(20) = [&
      & 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3]
   integer, parameter :: lsh(3, 20) = reshape([&
      & 0, 0, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0, &
      & 0, 1, 0,  0, 1, 0,  0, 1, 2,  0, 1, 0,  0, 1, 2,  0, 1, 2,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 0,  0, 1, 2], &
      & shape(lsh))
   integer, parameter :: pqn(3, 20) = reshape([&
      & 1, 0, 0,  1, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0, &
      & 2, 2, 0,  2, 2, 0,  2, 2, 3,  3, 3, 0,  3, 3, 3,  3, 3, 3,  3, 3, 3, &
      & 3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  4, 4, 0,  4, 4, 3], &
      & shape(pqn))
   real(wp), parameter :: zeta(3, 20) = reshape([&
      & 1.230000_wp, 0.000000_wp, 0.000000_wp, 1.669667_wp, 1.500000_wp, 0.000000_wp, &
      & 0.750060_wp, 0.557848_wp, 0.000000_wp, 1.034720_wp, 0.949332_wp, 0.000000_wp, &
      & 1.479444_wp, 1.479805_wp, 0.000000_wp, 2.096432_wp, 1.800000_wp, 0.000000_wp, &
      & 2.339881_wp, 2.014332_wp, 0.000000_wp, 2.439742_wp, 2.137023_wp, 0.000000_wp, &
      & 2.416361_wp, 2.308399_wp, 0.000000_wp, 3.084104_wp, 2.312051_wp, 2.815609_wp, &
      & 0.763787_wp, 0.573553_wp, 0.000000_wp, 1.184203_wp, 0.717769_wp, 1.300000_wp, &
      & 1.352531_wp, 1.391201_wp, 1.000000_wp, 1.773917_wp, 1.718996_wp, 1.250000_wp, &
      & 1.816945_wp, 1.903247_wp, 1.167533_wp, 1.981333_wp, 2.025643_wp, 1.702555_wp, &
      & 2.485265_wp, 2.199650_wp, 2.476089_wp, 2.329679_wp, 2.149419_wp, 1.950531_wp, &
      & 0.875961_wp, 0.631694_wp, 0.000000_wp, 1.267130_wp, 0.786247_wp, 1.380000_wp],&
      & shape(zeta))

   integer :: isp, izp, ish, stat
   integer, allocatable :: nshell(:)
   type(cgto_type), allocatable :: cgto(:, :)

   nshell = nsh(mol%num)
   allocate(cgto(maxval(nshell), mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, nshell(isp)
         call slater_to_gauss(ng, pqn(ish, izp), lsh(ish, izp), zeta(ish, izp), &
            & cgto(ish, isp), .true., stat)
      end do
   end do

   call new_basis(bas, mol, nshell, cgto, 1.0_wp)

end subroutine make_basis


subroutine test_overlap_mol(error, mol, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: ref(:, :)

   type(basis_type) :: bas
   real(wp), allocatable :: lattr(:, :), overlap(:, :)
   real(wp) :: cutoff
   integer :: ii, jj

   call make_basis(bas, mol, 6)
   call check(error, bas%nao, size(ref, 1))
   if (allocated(error)) return

   cutoff = get_cutoff(bas)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

   allocate(overlap(bas%nao, bas%nao))
   call get_overlap(mol, lattr, cutoff, bas, overlap)

   !where(abs(overlap) < thr) overlap = 0.0_wp
   !print '(*(6x,"&", 3(es20.14e1, "_wp":, ","), "&", /))', overlap

   do ii = 1, size(overlap, 2)
      do jj = 1, size(overlap, 1)
         call check(error, overlap(jj, ii), ref(jj, ii), thr=thr)
         if (allocated(error)) return
      end do
   end do


end subroutine test_overlap_mol

subroutine test_overlap_alh3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 12
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 9.99999999869333E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 4.02664545809937E-1_wp, 4.02664545809937E-1_wp, 4.02664545809939E-1_wp,&
      & 0.00000000000000E+0_wp, 9.99999999998060E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 4.32775379754938E-1_wp,-4.32775379754938E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999998060E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999998060E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.49862982000156E-1_wp,-2.49862982000156E-1_wp, 4.99725964000314E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830206E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.26789082148469E-1_wp,-2.26789082148469E-1_wp,-2.26789082148469E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830206E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999830206E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830206E-1_wp, 0.00000000000000E+0_wp,&
      &-1.96405106441530E-1_wp,-1.96405106441530E-1_wp, 3.92810212883060E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830206E-1_wp,&
      &-3.40183623222704E-1_wp, 3.40183623222704E-1_wp, 0.00000000000000E+0_wp,&
      & 4.02664545809937E-1_wp, 4.32775379754938E-1_wp, 0.00000000000000E+0_wp,&
      &-2.49862982000156E-1_wp,-2.26789082148469E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.96405106441530E-1_wp,-3.40183623222704E-1_wp,&
      & 9.99999999881495E-1_wp, 3.54600353330803E-2_wp, 3.54600353330805E-2_wp,&
      & 4.02664545809937E-1_wp,-4.32775379754938E-1_wp, 0.00000000000000E+0_wp,&
      &-2.49862982000156E-1_wp,-2.26789082148469E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.96405106441530E-1_wp, 3.40183623222704E-1_wp,&
      & 3.54600353330803E-2_wp, 9.99999999881495E-1_wp, 3.54600353330805E-2_wp,&
      & 4.02664545809939E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 4.99725964000314E-1_wp,-2.26789082148469E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 3.92810212883060E-1_wp, 0.00000000000000E+0_wp,&
      & 3.54600353330805E-2_wp, 3.54600353330805E-2_wp, 9.99999999881495E-1_wp],&
      & shape(overlap))
   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "AlH3")
   call test_overlap_mol(error, mol, overlap)

end subroutine test_overlap_alh3

subroutine test_overlap_bh3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 7
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 1.00000000000600E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 4.80899928496809E-1_wp, 4.80899928496809E-1_wp,&
      & 4.80899928496810E-1_wp, 0.00000000000000E+0_wp, 9.99999999925689E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 4.26420603078617E-1_wp,&
      &-4.26420603078617E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999925689E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp,-2.46194049975443E-1_wp,-2.46194049975443E-1_wp,&
      & 4.92388099950885E-1_wp, 4.80899928496809E-1_wp, 4.26420603078617E-1_wp,&
      & 0.00000000000000E+0_wp,-2.46194049975443E-1_wp, 9.99999999881495E-1_wp,&
      & 1.10583333710332E-1_wp, 1.10583333710331E-1_wp, 4.80899928496809E-1_wp,&
      &-4.26420603078617E-1_wp, 0.00000000000000E+0_wp,-2.46194049975443E-1_wp,&
      & 1.10583333710332E-1_wp, 9.99999999881495E-1_wp, 1.10583333710331E-1_wp,&
      & 4.80899928496810E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 4.92388099950885E-1_wp, 1.10583333710331E-1_wp, 1.10583333710331E-1_wp,&
      & 9.99999999881495E-1_wp],shape(overlap))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "BH3")
   call test_overlap_mol(error, mol, overlap)

end subroutine test_overlap_bh3

subroutine test_overlap_beh2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 6
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 1.00000000000600E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 4.75332491800271E-1_wp, 4.75332491800271E-1_wp,&
      & 0.00000000000000E+0_wp, 9.99999999925689E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999925689E-1_wp,&
      & 0.00000000000000E+0_wp,-5.55623013601629E-1_wp, 5.55623013601629E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 4.75332491800271E-1_wp, 0.00000000000000E+0_wp,-5.55623013601629E-1_wp,&
      & 0.00000000000000E+0_wp, 9.99999999881495E-1_wp, 4.00645037406879E-2_wp,&
      & 4.75332491800271E-1_wp, 0.00000000000000E+0_wp, 5.55623013601629E-1_wp,&
      & 0.00000000000000E+0_wp, 4.00645037406879E-2_wp, 9.99999999881495E-1_wp],&
      & shape(overlap))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "BeH2")
   call test_overlap_mol(error, mol, overlap)

end subroutine test_overlap_beh2

subroutine test_overlap_ch4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 8
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 1.00000000000600E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 4.22144465267144E-1_wp, 4.22144465267144E-1_wp,&
      & 4.22144465267144E-1_wp, 4.22144465267144E-1_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 2.57682345348632E-1_wp,-2.57682345348632E-1_wp,-2.57682345348632E-1_wp,&
      & 2.57682345348632E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 2.57682345348632E-1_wp,&
      & 2.57682345348632E-1_wp,-2.57682345348632E-1_wp,-2.57682345348632E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp,-2.57682345348632E-1_wp, 2.57682345348632E-1_wp,&
      &-2.57682345348632E-1_wp, 2.57682345348632E-1_wp, 4.22144465267144E-1_wp,&
      & 2.57682345348632E-1_wp, 2.57682345348632E-1_wp,-2.57682345348632E-1_wp,&
      & 9.99999999881495E-1_wp, 1.71953274424808E-1_wp, 1.71953274424808E-1_wp,&
      & 1.71953274424808E-1_wp, 4.22144465267144E-1_wp,-2.57682345348632E-1_wp,&
      & 2.57682345348632E-1_wp, 2.57682345348632E-1_wp, 1.71953274424808E-1_wp,&
      & 9.99999999881495E-1_wp, 1.71953274424808E-1_wp, 1.71953274424808E-1_wp,&
      & 4.22144465267144E-1_wp,-2.57682345348632E-1_wp,-2.57682345348632E-1_wp,&
      &-2.57682345348632E-1_wp, 1.71953274424808E-1_wp, 1.71953274424808E-1_wp,&
      & 9.99999999881495E-1_wp, 1.71953274424808E-1_wp, 4.22144465267144E-1_wp,&
      & 2.57682345348632E-1_wp,-2.57682345348632E-1_wp, 2.57682345348632E-1_wp,&
      & 1.71953274424808E-1_wp, 1.71953274424808E-1_wp, 1.71953274424808E-1_wp,&
      & 9.99999999881495E-1_wp],shape(overlap))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "CH4")
   call test_overlap_mol(error, mol, overlap)

end subroutine test_overlap_ch4

subroutine test_overlap_cl2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 18
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 9.99999999869332E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.17210205959034E-2_wp, 0.00000000000000E+0_wp,-1.66880110378166E-1_wp,&
      & 0.00000000000000E+0_wp, 1.08719248478375E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999998060E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 6.66413853324784E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-7.28795473739356E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999998060E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.66880110378166E-1_wp, 0.00000000000000E+0_wp,-2.43857348208300E-1_wp,&
      & 0.00000000000000E+0_wp, 1.54878516243462E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999998060E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 6.66413853324784E-2_wp, 0.00000000000000E+0_wp,-7.28795473739356E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830206E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.08719248478375E-1_wp, 0.00000000000000E+0_wp,-1.54878516243462E-1_wp,&
      & 0.00000000000000E+0_wp, 1.20587356016693E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830205E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 7.28795473739356E-2_wp, 0.00000000000000E+0_wp,-8.83217466681701E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999830205E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 7.28795473739356E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-8.83217466681701E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830206E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 1.65407586749691E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.0000000000000E+0_wp, 9.99999999830205E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.65407586749691E-2_wp,&
      & 9.17210205959034E-2_wp, 0.00000000000000E+0_wp, 1.66880110378166E-1_wp,&
      & 0.00000000000000E+0_wp, 1.08719248478375E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999869332E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 6.66413853324784E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 7.28795473739356E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999998060E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.66880110378166E-1_wp, 0.00000000000000E+0_wp,-2.43857348208300E-1_wp,&
      & 0.00000000000000E+0_wp,-1.54878516243462E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999998060E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 6.66413853324784E-2_wp, 0.00000000000000E+0_wp, 7.28795473739356E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999998060E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.08719248478375E-1_wp, 0.00000000000000E+0_wp, 1.54878516243462E-1_wp,&
      & 0.00000000000000E+0_wp, 1.20587356016693E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830206E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-7.28795473739356E-2_wp, 0.00000000000000E+0_wp,-8.83217466681701E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830205E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-7.28795473739356E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-8.83217466681701E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999830205E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 1.65407586749691E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830206E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.65407586749691E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830205E-1_wp],&
      & shape(overlap))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "Cl2")
   call test_overlap_mol(error, mol, overlap)

end subroutine test_overlap_cl2

subroutine test_overlap_f2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 8
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 1.00000000000600E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 1.47958427103274E-1_wp, 0.00000000000000E+0_wp,&
      &-2.00620208590218E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 7.94087596423658E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 2.00620208590218E-1_wp,&
      & 0.00000000000000E+0_wp,-2.36369347922990E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 7.94087596423658E-2_wp, 1.47958427103274E-1_wp,&
      & 0.00000000000000E+0_wp, 2.00620208590218E-1_wp, 0.00000000000000E+0_wp,&
      & 1.00000000000600E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 7.94087596423658E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.00620208590218E-1_wp, 0.00000000000000E+0_wp,-2.36369347922990E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 7.94087596423658E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp],shape(overlap))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "F2")
   call test_overlap_mol(error, mol, overlap)

end subroutine test_overlap_f2

subroutine test_overlap_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 2
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 9.99999999881495E-1_wp, 6.61346655776026E-1_wp, 6.61346655776026E-1_wp,&
      & 9.99999999881495E-1_wp],shape(overlap))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "H2")
   call test_overlap_mol(error, mol, overlap)

end subroutine test_overlap_h2

subroutine test_overlap_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 5
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 1.00000000000600E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 3.99089038384911E-1_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 4.65790780903622E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 3.99089038384911E-1_wp,&
      & 0.00000000000000E+0_wp, 4.65790780903622E-1_wp, 0.00000000000000E+0_wp,&
      & 9.99999999881495E-1_wp],shape(overlap))


   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_overlap_mol(error, mol, overlap)

end subroutine test_overlap_lih


subroutine test_overlap_grad_ss(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6
   integer :: stat, i
   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3), r2
   real(wp) :: overlap(1, 1), doverlapi(3, 1, 1), doverlapj(3, 1, 1), sr(1, 1), sl(1, 1)
   real(wp), parameter :: step = 1.0e-6_wp

   call slater_to_gauss(ng, 2, 0, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 1, 0, 1.0_wp, cgtoj, .true., stat)

   call random_number(vec)
   vec = vec - 0.5_wp
   r2 = sum(vec**2)

   call overlap_grad_cgto(cgtoi, cgtoj, r2, vec, 100.0_wp, overlap, doverlapj)

   vec(:) = -vec

   call overlap_grad_cgto(cgtoj, cgtoi, r2, vec, 100.0_wp, overlap, doverlapi)

   do i = 1, 3
      call check(error, doverlapi(i, 1, 1), -doverlapj(i, 1, 1), thr=thr)
      if (allocated(error)) exit
   end do
   if (allocated(error)) return

   do i = 1, 3
      vec(i) = vec(i) + step
      r2 = sum(vec**2)
      call overlap_cgto(cgtoj, cgtoi, r2, vec, 100.0_wp, sr)
      vec(i) = vec(i) - 2*step
      r2 = sum(vec**2)
      call overlap_cgto(cgtoj, cgtoi, r2, vec, 100.0_wp, sl)
      vec(i) = vec(i) + step
      doverlapj(i, :, :) = 0.5_wp * (sr - sl) / step
   end do

   do i = 1, 3
      call check(error, doverlapi(i, 1, 1), doverlapj(i, 1, 1), thr=thr)
      if (allocated(error)) exit
   end do
   if (allocated(error)) return

end subroutine test_overlap_grad_ss


subroutine test_overlap_grad_pp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6
   integer :: stat, i, j, k
   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3), r2
   real(wp) :: overlap(3, 3), doverlapi(3, 3, 3), doverlapj(3, 3, 3), sr(3, 3), sl(3, 3)
   real(wp), parameter :: step = 1.0e-6_wp

   call slater_to_gauss(ng, 3, 1, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 2, 1, 1.0_wp, cgtoj, .true., stat)

   call random_number(vec)
   vec = vec - 0.5_wp
   r2 = sum(vec**2)

   call overlap_grad_cgto(cgtoi, cgtoj, r2, vec, 100.0_wp, overlap, doverlapj)

   vec(:) = -vec

   call overlap_grad_cgto(cgtoj, cgtoi, r2, vec, 100.0_wp, overlap, doverlapi)

   lp: do i = 1, 3
      do j = 1, 3
         do k = 1, 3
            call check(error, doverlapi(i, j, k), -doverlapj(i, k, j), thr=thr)
            if (allocated(error)) exit
         end do
      end do
   end do lp
   if (allocated(error)) return

   do i = 1, 3
      vec(i) = vec(i) + step
      r2 = sum(vec**2)
      call overlap_cgto(cgtoj, cgtoi, r2, vec, 100.0_wp, sr)
      vec(i) = vec(i) - 2*step
      r2 = sum(vec**2)
      call overlap_cgto(cgtoj, cgtoi, r2, vec, 100.0_wp, sl)
      vec(i) = vec(i) + step
      doverlapj(i, :, :) = 0.5_wp * (sr - sl) / step
   end do

   num: do i = 1, 3
      do j = 1, 3
         do k = 1, 3
            call check(error, doverlapi(i, j, k), doverlapj(i, j, k), thr=thr)
            if (allocated(error)) exit num
         end do
      end do
   end do num
   if (allocated(error)) return

end subroutine test_overlap_grad_pp


subroutine test_overlap_grad_dd(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6
   integer :: stat, i, j, k
   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3), r2
   real(wp) :: overlap(5, 5), doverlapi(3, 5, 5), doverlapj(3, 5, 5), sl(5, 5), sr(5, 5)
   real(wp), parameter :: step = 1.0e-6_wp

   call slater_to_gauss(ng, 4, 2, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 3, 2, 1.0_wp, cgtoj, .true., stat)

   call random_number(vec)
   vec = vec - 0.5_wp
   r2 = sum(vec**2)

   call overlap_grad_cgto(cgtoi, cgtoj, r2, vec, 100.0_wp, overlap, doverlapj)

   vec(:) = -vec

   call overlap_grad_cgto(cgtoj, cgtoi, r2, vec, 100.0_wp, overlap, doverlapi)

   lp: do i = 1, 3
      do j = 1, 5
         do k = 1, 5
            call check(error, doverlapi(i, j, k), -doverlapj(i, k, j), thr=thr)
            if (allocated(error)) exit lp
         end do
      end do
   end do lp
   if (allocated(error)) return

   do i = 1, 3
      vec(i) = vec(i) + step
      r2 = sum(vec**2)
      call overlap_cgto(cgtoj, cgtoi, r2, vec, 100.0_wp, sr)
      vec(i) = vec(i) - 2*step
      r2 = sum(vec**2)
      call overlap_cgto(cgtoj, cgtoi, r2, vec, 100.0_wp, sl)
      vec(i) = vec(i) + step
      doverlapj(i, :, :) = 0.5_wp * (sr - sl) / step
   end do

   num: do i = 1, 3
      do j = 1, 5
         do k = 1, 5
            call check(error, doverlapi(i, j, k), doverlapj(i, j, k), thr=thr)
            if (allocated(error)) exit num
         end do
      end do
   end do num
   if (allocated(error)) return

end subroutine test_overlap_grad_dd


end module test_integral_overlap
