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

module test_ncoord
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, &
   & test_failed
   use mctc_io_structure, only : structure_type
   use mstore, only : get_structure
   use tblite_cutoff, only : get_lattice_points
   use tblite_data_covrad, only : get_covalent_rad
   use tblite_data_paulingen, only : get_pauling_en
   use tblite_ncoord_gfn
   use tblite_ncoord_exp
   use tblite_ncoord_erf
   use tblite_ncoord_erf_en
   use tblite_ncoord_type !, only : get_coordination_number
   implicit none
   private

   public :: collect_ncoord

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))


contains


   !> Collect all exported unit tests
   subroutine collect_ncoord(testsuite)

      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
      & new_unittest("cn-mb01_gfn", test_cn_mb01_gfn), &
      & new_unittest("cn-mb02_gfn", test_cn_mb02_gfn), &
      & new_unittest("cn-mb03_gfn", test_cn_mb03_gfn), &
      & new_unittest("cn-acetic_gfn", test_cn_acetic_gfn), &
      & new_unittest("dcndr-mb04_gfn", test_dcndr_mb04_gfn), &
      & new_unittest("dcndr-mb05_gfn", test_dcndr_mb05_gfn), &
      & new_unittest("dcndr-ammonia_gfn", test_dcndr_ammonia_gfn), &
      & new_unittest("dcndL-mb06_gfn", test_dcndL_mb06_gfn), &
      & new_unittest("dcndL-mb07_gfn", test_dcndL_mb07_gfn), &
      & new_unittest("dcndL-antracene_gfn", test_dcndL_anthracene_gfn), &
      & new_unittest("cn-mb01_exp", test_cn_mb01_exp), &
      & new_unittest("cn-mb02_exp", test_cn_mb02_exp), &
      & new_unittest("cn-mb03_exp", test_cn_mb03_exp), &
      & new_unittest("cn-acetic_exp", test_cn_acetic_exp), &
      & new_unittest("dcndr-mb04_exp", test_dcndr_mb04_exp), &
      & new_unittest("dcndr-mb05_exp", test_dcndr_mb05_exp), &
      & new_unittest("dcndr-ammonia_exp", test_dcndr_ammonia_exp), &
      & new_unittest("dcndL-mb06_exp", test_dcndL_mb06_exp), &
      & new_unittest("dcndL-mb07_exp", test_dcndL_mb07_exp), &
      & new_unittest("dcndL-antracene_exp", test_dcndL_anthracene_exp), &
      & new_unittest("cn-mb01_erf", test_cn_mb01_erf), &
      & new_unittest("cn-mb02_erf", test_cn_mb02_erf), &
      & new_unittest("cn-mb03_erf", test_cn_mb03_erf), &
      & new_unittest("cn-acetic_erf", test_cn_acetic_erf), &
      & new_unittest("dcndr-mb04_erf", test_dcndr_mb04_erf), &
      & new_unittest("dcndr-mb05_erf", test_dcndr_mb05_erf), &
      & new_unittest("dcndr-ammonia_erf", test_dcndr_ammonia_erf), &
      & new_unittest("dcndL-mb06_erf", test_dcndL_mb06_erf), &
      & new_unittest("dcndL-mb07_erf", test_dcndL_mb07_erf), &
      & new_unittest("dcndL-antracene_erf", test_dcndL_anthracene_erf), &
      & new_unittest("cn-mb01_erf_en", test_cn_mb01_erf_en), &
      & new_unittest("cn-mb02_erf_en", test_cn_mb02_erf_en), &
      & new_unittest("cn-mb03_erf_en", test_cn_mb03_erf_en), &
      & new_unittest("cn-acetic_erf_en", test_cn_acetic_erf_en), &
      & new_unittest("dcndr-mb04_erf_en", test_dcndr_mb04_erf_en), &
      & new_unittest("dcndr-mb05_erf_en", test_dcndr_mb05_erf_en), &
      & new_unittest("dcndr-ammonia_erf_en", test_dcndr_ammonia_erf_en), &
      & new_unittest("dcndL-mb06_erf_en", test_dcndL_mb06_erf_en), &
      & new_unittest("dcndL-mb07_erf_en", test_dcndL_mb07_erf_en), &
      & new_unittest("dcndL-antracene_erf_en", test_dcndL_anthracene_erf_en) &
      & ]

   end subroutine collect_ncoord


   subroutine test_cn_gen(error, mol, ncoord, cutoff, ref)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type) :: mol

      !> Coordination number type
      class(ncoord_type)   :: ncoord

      real(wp), intent(in) :: cutoff

      !> Reference CNs
      real(wp), intent(in) :: ref(:)
      real(wp), allocatable :: cn(:)
      real(wp), allocatable :: lattr(:, :)

      allocate(cn(mol%nat))

      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      call get_coordination_number(ncoord, mol, lattr, cutoff, cn)

      if (any(abs(cn - ref) > thr)) then
         call test_failed(error, "Coordination numbers do not match")
         print'(3es21.14)', cn
      end if

   end subroutine test_cn_gen


   subroutine test_numgrad(error, mol, ncoord, cutoff)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type), intent(inout) :: mol

      !> Coordination number type
      class(ncoord_type)   :: ncoord

      real(wp), intent(in) :: cutoff

      integer :: iat, ic
      real(wp), allocatable :: cn(:), cnr(:), cnl(:)
      real(wp), allocatable :: dcndr(:, :, :), dcndL(:, :, :)
      real(wp), allocatable :: numdr(:, :, :)
      real(wp), allocatable :: lattr(:, :)
      real(wp), parameter :: step = 1.0e-6_wp

      allocate(cn(mol%nat), cnr(mol%nat), cnl(mol%nat), &
      & dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & numdr(3, mol%nat, mol%nat))

      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

      do iat = 1, mol%nat
         do ic = 1, 3
            mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
            call get_coordination_number(ncoord, mol, lattr, cutoff, cnr)
            mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
            call get_coordination_number(ncoord, mol, lattr, cutoff, cnl)
            mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
            numdr(ic, iat, :) = 0.5_wp*(cnr - cnl)/step
         end do
      end do

      call get_coordination_number(ncoord, mol, lattr, cutoff, cn, dcndr, dcndL)

      if (any(abs(dcndr - numdr) > thr2)) then
         call test_failed(error, "Derivative of coordination number does not match")
      end if

   end subroutine test_numgrad


   subroutine test_numsigma(error, mol, ncoord, cutoff)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type), intent(inout) :: mol

      !> Coordination number type
      class(ncoord_type)   :: ncoord

      real(wp), intent(in) :: cutoff

      integer :: ic, jc
      real(wp) :: eps(3, 3)
      real(wp), allocatable :: cn(:), cnr(:), cnl(:), xyz(:, :)
      real(wp), allocatable :: dcndr(:, :, :), dcndL(:, :, :)
      real(wp), allocatable :: numdL(:, :, :)
      real(wp), allocatable :: lattr(:, :), trans(:, :)
      real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
      real(wp), parameter :: step = 1.0e-6_wp

      allocate(cn(mol%nat), cnr(mol%nat), cnl(mol%nat), &
      & dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), xyz(3, mol%nat), &
      & numdL(3, 3, mol%nat))
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

      eps(:, :) = unity
      xyz(:, :) = mol%xyz
      trans = lattr
      do ic = 1, 3
         do jc = 1, 3
            eps(jc, ic) = eps(jc, ic) + step
            mol%xyz(:, :) = matmul(eps, xyz)
            lattr(:, :) = matmul(eps, trans)
            call get_coordination_number(ncoord, mol, lattr, cutoff, cnr)
            eps(jc, ic) = eps(jc, ic) - 2*step
            mol%xyz(:, :) = matmul(eps, xyz)
            lattr(:, :) = matmul(eps, trans)
            call get_coordination_number(ncoord, mol, lattr, cutoff, cnl)
            eps(jc, ic) = eps(jc, ic) + step
            mol%xyz(:, :) = xyz
            lattr(:, :) = trans
            numdL(jc, ic, :) = 0.5_wp*(cnr - cnl)/step
         end do
      end do

      call get_coordination_number(ncoord, mol, lattr, cutoff, cn, dcndr, dcndL)

      if (any(abs(dcndL - numdL) > thr2)) then
         call test_failed(error, "Derivative of coordination number does not match")
      end if

   end subroutine test_numsigma


   !> ----------------------------------------------------
   !> Tests for double-exponential (GFN) coordination number
   !> ----------------------------------------------------
   subroutine test_cn_mb01_gfn(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(gfn_ncoord_type) :: gfn_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp
      real(wp), parameter :: ref(16) = [&
      & 4.11453659059991E+0_wp, 9.32058998762811E-1_wp, 2.03554597140311E+0_wp, &
      & 1.42227835389358E+0_wp, 1.12812426574031E+0_wp, 1.05491602558828E+0_wp, &
      & 1.52709064704269E+0_wp, 1.95070367247232E+0_wp, 3.83759889196540E+0_wp, &
      & 1.09388314007182E+0_wp, 1.07090773695340E+0_wp, 2.00285254082830E+0_wp, &
      & 4.36400837813955E+0_wp, 3.83469860546080E+0_wp, 3.91542517673963E+0_wp, &
      & 5.58571682419960E+0_wp]

      call get_structure(mol, "MB16-43", "01")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_gfn_ncoord(gfn_ncoord, mol, cutoff, rcov)
      call test_cn_gen(error, mol, gfn_ncoord, cutoff, ref)

   end subroutine test_cn_mb01_gfn


   subroutine test_cn_mb02_gfn(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(gfn_ncoord_type) :: gfn_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp
      real(wp), parameter :: ref(16) = [&
      & 9.16720394154780E-1_wp, 3.81678481096046E+0_wp, 3.58669504034282E+0_wp, &
      & 2.90017371823910E+0_wp, 5.09896950232009E+0_wp, 1.13810142820063E+0_wp, &
      & 9.35735928561309E-1_wp, 9.43324181507483E-1_wp, 4.86704300506033E+0_wp, &
      & 1.22024405676881E+0_wp, 3.86661472534511E+0_wp, 4.02897853934353E+0_wp, &
      & 1.96869999324996E+0_wp, 9.04495366214972E-1_wp, 1.49977675358986E+0_wp, &
      & 2.04937665719480E+0_wp]

      call get_structure(mol, "MB16-43", "02")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_gfn_ncoord(gfn_ncoord, mol, cutoff, rcov)
      call test_cn_gen(error, mol,  gfn_ncoord, cutoff, ref)

   end subroutine test_cn_mb02_gfn


   subroutine test_cn_mb03_gfn(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(gfn_ncoord_type) :: gfn_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp
      real(wp), parameter :: ref(16) = [&
      & 4.03003984525826E+0_wp, 2.88404829383272E+0_wp, 1.12030447706437E+0_wp, &
      & 4.89326653038565E+0_wp, 6.97177319120287E+0_wp, 4.43722632690247E+0_wp, &
      & 4.74866591812568E+0_wp, 1.30338583081493E+0_wp, 1.06130190313902E+0_wp, &
      & 1.04048373066340E+0_wp, 1.92270060977307E+0_wp, 2.98737904327655E+0_wp, &
      & 4.69081207436918E+0_wp, 1.26260125612717E+0_wp, 3.36368006785498E+0_wp, &
      & 1.35743886389476E+0_wp]

      call get_structure(mol, "MB16-43", "03")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_gfn_ncoord(gfn_ncoord, mol, cutoff, rcov)
      call test_cn_gen(error, mol, gfn_ncoord, cutoff, ref)

   end subroutine test_cn_mb03_gfn


   subroutine test_cn_acetic_gfn(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(gfn_ncoord_type) :: gfn_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp
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

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_gfn_ncoord(gfn_ncoord, mol, cutoff, rcov)
      call test_cn_gen(error, mol, gfn_ncoord, cutoff, ref)

   end subroutine test_cn_acetic_gfn


   subroutine test_dcndr_mb04_gfn(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(gfn_ncoord_type) :: gfn_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "MB16-43", "04")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_gfn_ncoord(gfn_ncoord, mol, cutoff, rcov)
      call test_numgrad(error, mol, gfn_ncoord, cutoff)

   end subroutine test_dcndr_mb04_gfn


   subroutine test_dcndr_mb05_gfn(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(gfn_ncoord_type) :: gfn_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "MB16-43", "05")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_gfn_ncoord(gfn_ncoord, mol, cutoff, rcov)
      call test_numgrad(error, mol, gfn_ncoord, cutoff)

   end subroutine test_dcndr_mb05_gfn


   subroutine test_dcndr_ammonia_gfn(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(gfn_ncoord_type) :: gfn_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "X23", "ammonia")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_gfn_ncoord(gfn_ncoord, mol, cutoff, rcov)
      call test_numgrad(error, mol, gfn_ncoord, cutoff)

   end subroutine test_dcndr_ammonia_gfn


   subroutine test_dcndL_mb06_gfn(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(gfn_ncoord_type) :: gfn_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "MB16-43", "06")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_gfn_ncoord(gfn_ncoord, mol, cutoff, rcov)
      call test_numsigma(error, mol, gfn_ncoord, cutoff)

   end subroutine test_dcndL_mb06_gfn


   subroutine test_dcndL_mb07_gfn(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(gfn_ncoord_type) :: gfn_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "MB16-43", "07")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_gfn_ncoord(gfn_ncoord, mol, cutoff, rcov)
      call test_numsigma(error, mol, gfn_ncoord, cutoff)

   end subroutine test_dcndL_mb07_gfn


   subroutine test_dcndL_anthracene_gfn(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(gfn_ncoord_type) :: gfn_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "X23", "anthracene")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_gfn_ncoord(gfn_ncoord, mol, cutoff, rcov)
      call test_numsigma(error, mol, gfn_ncoord, cutoff)

   end subroutine test_dcndL_anthracene_gfn


   !> ----------------------------------------------------
   !> Tests for mono-exponential coordination number
   !> ----------------------------------------------------
   subroutine test_cn_mb01_exp(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(exp_ncoord_type) :: exp_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp
      real(wp), parameter :: ref(16) = [&
      & 4.15066368951397E+0_wp, 9.78868026389781E-1_wp, 2.01080985633859E+0_wp, &
      & 1.47865697827818E+0_wp, 1.03577822442117E+0_wp, 1.01206994314781E+0_wp, &
      & 1.50329777127401E+0_wp, 1.99858468272609E+0_wp, 3.89181927539324E+0_wp, &
      & 1.04323373360740E+0_wp, 1.01526584450636E+0_wp, 1.99315213227354E+0_wp, &
      & 4.63526560889683E+0_wp, 3.87312260639335E+0_wp, 3.99316800677884E+0_wp, &
      & 5.45068226903888E+0_wp]

      call get_structure(mol, "MB16-43", "01")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_exp_ncoord(exp_ncoord, mol, cutoff, rcov)
      call test_cn_gen(error, mol, exp_ncoord, cutoff, ref)

   end subroutine test_cn_mb01_exp


   subroutine test_cn_mb02_exp(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(exp_ncoord_type) :: exp_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp
      real(wp), parameter :: ref(16) = [&
      & 9.68434565947196E-1_wp, 3.94488220883276E+0_wp, 3.82701677880409E+0_wp, &
      & 2.99161201243234E+0_wp, 5.46904892971914E+0_wp, 1.06438740698832E+0_wp, &
      & 9.77627896999762E-1_wp, 9.81715643929621E-1_wp, 5.06702924169705E+0_wp, &
      & 1.08324093335500E+0_wp, 4.00161384251868E+0_wp, 4.03067393311321E+0_wp, &
      & 2.00249301823990E+0_wp, 9.73399178780401E-1_wp, 1.66868575900646E+0_wp, &
      & 2.03273936244268E+0_wp]

      call get_structure(mol, "MB16-43", "02")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_exp_ncoord(exp_ncoord, mol, cutoff, rcov)
      call test_cn_gen(error, mol,  exp_ncoord, cutoff, ref)

   end subroutine test_cn_mb02_exp


   subroutine test_cn_mb03_exp(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(exp_ncoord_type) :: exp_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp
      real(wp), parameter :: ref(16) = [&
      & 4.04992403668766E+0_wp, 2.98487291056010E+0_wp, 1.03236609075831E+0_wp, &
      & 4.86782876076800E+0_wp, 7.48154833549122E+0_wp, 4.76588477343835E+0_wp, &
      & 4.92432613481557E+0_wp, 1.19761963808520E+0_wp, 1.01809037129368E+0_wp, &
      & 1.01042713594272E+0_wp, 1.99186789992505E+0_wp, 3.03157225205500E+0_wp, &
      & 4.89702969622395E+0_wp, 1.11663432026701E+0_wp, 3.52903544775683E+0_wp, &
      & 1.33563958333433E+0_wp]

      call get_structure(mol, "MB16-43", "03")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_exp_ncoord(exp_ncoord, mol, cutoff, rcov)
      call test_cn_gen(error, mol, exp_ncoord, cutoff, ref)

   end subroutine test_cn_mb03_exp


   subroutine test_cn_acetic_exp(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(exp_ncoord_type) :: exp_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp
      real(wp), parameter :: ref(32) = [&
      & 1.04073967324021E+0_wp, 1.04075472444211E+0_wp, 1.04072486016808E+0_wp, &
      & 1.04074563131443E+0_wp, 1.00106790331471E+0_wp, 1.00106709316340E+0_wp, &
      & 1.00105733071182E+0_wp, 1.00011861283032E+0_wp, 1.00011458903338E+0_wp, &
      & 1.00010904222439E+0_wp, 1.00010907333542E+0_wp, 9.99919377862809E-1_wp, &
      & 9.99937754940446E-1_wp, 9.99930065513371E-1_wp, 9.99925395685880E-1_wp, &
      & 1.00107017611203E+0_wp, 3.03354090037860E+0_wp, 3.03354424850001E+0_wp, &
      & 3.03353363463649E+0_wp, 3.03354145642558E+0_wp, 4.03132899427226E+0_wp, &
      & 4.03134766043645E+0_wp, 4.03131729363474E+0_wp, 4.03132452876117E+0_wp, &
      & 2.03473414573198E+0_wp, 2.03477399887022E+0_wp, 2.03473692087023E+0_wp, &
      & 2.03476613472564E+0_wp, 1.09702037696180E+0_wp, 1.09703703566506E+0_wp, &
      & 1.09702583457622E+0_wp, 1.09704757516935E+0_wp]

      call get_structure(mol, "X23", "acetic")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_exp_ncoord(exp_ncoord, mol, cutoff, rcov)
      call test_cn_gen(error, mol, exp_ncoord, cutoff, ref)

   end subroutine test_cn_acetic_exp


   subroutine test_dcndr_mb04_exp(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(exp_ncoord_type) :: exp_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "MB16-43", "04")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_exp_ncoord(exp_ncoord, mol, cutoff, rcov)
      call test_numgrad(error, mol, exp_ncoord, cutoff)

   end subroutine test_dcndr_mb04_exp


   subroutine test_dcndr_mb05_exp(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(exp_ncoord_type) :: exp_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "MB16-43", "05")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_exp_ncoord(exp_ncoord, mol, cutoff, rcov)
      call test_numgrad(error, mol, exp_ncoord, cutoff)

   end subroutine test_dcndr_mb05_exp


   subroutine test_dcndr_ammonia_exp(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(exp_ncoord_type) :: exp_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "X23", "ammonia")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_exp_ncoord(exp_ncoord, mol, cutoff, rcov)
      call test_numgrad(error, mol, exp_ncoord, cutoff)

   end subroutine test_dcndr_ammonia_exp


   subroutine test_dcndL_mb06_exp(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(exp_ncoord_type) :: exp_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "MB16-43", "06")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_exp_ncoord(exp_ncoord, mol, cutoff, rcov)
      call test_numsigma(error, mol, exp_ncoord, cutoff)

   end subroutine test_dcndL_mb06_exp


   subroutine test_dcndL_mb07_exp(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(exp_ncoord_type) :: exp_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "MB16-43", "07")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_exp_ncoord(exp_ncoord, mol, cutoff, rcov)
      call test_numsigma(error, mol, exp_ncoord, cutoff)

   end subroutine test_dcndL_mb07_exp


   subroutine test_dcndL_anthracene_exp(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(exp_ncoord_type) :: exp_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "X23", "anthracene")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_exp_ncoord(exp_ncoord, mol, cutoff, rcov)
      call test_numsigma(error, mol, exp_ncoord, cutoff)

   end subroutine test_dcndL_anthracene_exp


   !> ----------------------------------------------------
   !> Tests for error-function based CEH/GP3 coordination number 
   !> using the Pyykko covalent radii and Pauling EN
   !> ----------------------------------------------------
   subroutine test_cn_mb01_erf(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(erf_ncoord_type) :: erf_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp
      real(wp), parameter :: ref(16) = [&
      & +4.0293255193985615E+00_wp, +7.8510334142664440E-01_wp, +1.8101410224986185E+00_wp, & 
      & +1.2757066882291497E+00_wp, +1.0781936213942493E+00_wp, +9.2557169711499776E-01_wp, & 
      & +1.4615446456579280E+00_wp, +1.6847684089029940E+00_wp, +3.3648064939330400E+00_wp, & 
      & +1.0124657146010481E+00_wp, +9.4167385098988499E-01_wp, +1.8661189314109503E+00_wp, & 
      & +3.8220710498421062E+00_wp, +3.4377504655667011E+00_wp, +3.4349836685142909E+00_wp, & 
      & +5.2383444295555552E+00_wp]

      call get_structure(mol, "MB16-43", "01")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_erf_ncoord(erf_ncoord, mol, cutoff, rcov)
      call test_cn_gen(error, mol, erf_ncoord, cutoff, ref)

   end subroutine test_cn_mb01_erf


   subroutine test_cn_mb02_erf(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(erf_ncoord_type) :: erf_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp
      real(wp), parameter :: ref(16) = [&
      & +7.6933845548840352E-01_wp, +3.4390628220068522E+00_wp, +3.0968216100725039E+00_wp, & 
      & +2.5188744074433327E+00_wp, +4.4847318088928416E+00_wp, +1.0442798803389632E+00_wp, & 
      & +7.9199922696084635E-01_wp, +8.1060504520365595E-01_wp, +4.3849962975469783E+00_wp, & 
      & +1.1815409616134993E+00_wp, +3.5066198235979593E+00_wp, +3.7980320477763003E+00_wp, & 
      & +1.7195072406902749E+00_wp, +7.5226274752252320E-01_wp, +1.3098148757645638E+00_wp, & 
      & +1.9302354819530141E+00_wp]

      call get_structure(mol, "MB16-43", "02")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_erf_ncoord(erf_ncoord, mol, cutoff, rcov)
      call test_cn_gen(error, mol, erf_ncoord, cutoff, ref)

   end subroutine test_cn_mb02_erf


   subroutine test_cn_mb03_erf(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(erf_ncoord_type) :: erf_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp
      real(wp), parameter :: ref(16) = [&
      & +3.6560485845339090E+00_wp, +2.4744404243444635E+00_wp, +1.0472994524445995E+00_wp, & 
      & +4.6902866431890162E+00_wp, +6.1333128089596878E+00_wp, +3.9670698548854917E+00_wp, & 
      & +4.1909395790329649E+00_wp, +1.2475870668105131E+00_wp, +9.5657288286854358E-01_wp, & 
      & +9.1546525520671240E-01_wp, +1.6538274001105115E+00_wp, +2.7161929906912068E+00_wp, & 
      & +4.1336603661167155E+00_wp, +1.2590885278639488E+00_wp, +3.1241492930194554E+00_wp, & 
      & +1.2455016734877240E+00_wp]

      call get_structure(mol, "MB16-43", "03")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_erf_ncoord(erf_ncoord, mol, cutoff, rcov)
      call test_cn_gen(error, mol, erf_ncoord, cutoff, ref)

   end subroutine test_cn_mb03_erf


   subroutine test_cn_acetic_erf(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(erf_ncoord_type) :: erf_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp
      real(wp), parameter :: ref(32) = [&
      & +1.0461853752186683E+00_wp, +1.0462144049836850E+00_wp, +1.0461526901127018E+00_wp, & 
      & +1.0461927068096812E+00_wp, +8.3953739457313425E-01_wp, +8.3951267982957412E-01_wp, & 
      & +8.3945549246201678E-01_wp, +8.3657043332907832E-01_wp, +8.3655989299146305E-01_wp, & 
      & +8.3653217035350158E-01_wp, +8.3652859765160237E-01_wp, +8.3501289305917448E-01_wp, & 
      & +8.3510827350536454E-01_wp, +8.3505828052912601E-01_wp, +8.3504133253387924E-01_wp, & 
      & +8.3952949590967618E-01_wp, +2.7521054643316014E+00_wp, +2.7521582090229151E+00_wp, & 
      & +2.7520567053245850E+00_wp, +2.7521059694337313E+00_wp, +3.5004593976447507E+00_wp, & 
      & +3.5005482072646514E+00_wp, +3.5003932049077329E+00_wp, +3.5004460582287180E+00_wp, & 
      & +1.8964535656891905E+00_wp, +1.8965371449325730E+00_wp, +1.8964092180523289E+00_wp, & 
      & +1.8965025248107723E+00_wp, +1.3405486119794887E+00_wp, +1.3406581897362071E+00_wp, & 
      & +1.3405893806587921E+00_wp, +1.3406201938022535E+00_wp]

      call get_structure(mol, "X23", "acetic")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_erf_ncoord(erf_ncoord, mol, cutoff, rcov)
      call test_cn_gen(error, mol, erf_ncoord, cutoff, ref)

   end subroutine test_cn_acetic_erf


   subroutine test_dcndr_mb04_erf(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(erf_ncoord_type) :: erf_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "MB16-43", "04")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_erf_ncoord(erf_ncoord, mol, cutoff, rcov)
      call test_numgrad(error, mol, erf_ncoord, cutoff)

   end subroutine test_dcndr_mb04_erf


   subroutine test_dcndr_mb05_erf(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(erf_ncoord_type) :: erf_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "MB16-43", "05")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_erf_ncoord(erf_ncoord, mol, cutoff, rcov)
      call test_numgrad(error, mol, erf_ncoord, cutoff)

   end subroutine test_dcndr_mb05_erf


   subroutine test_dcndr_ammonia_erf(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(erf_ncoord_type) :: erf_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "X23", "ammonia")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_erf_ncoord(erf_ncoord, mol, cutoff, rcov)
      call test_numgrad(error, mol, erf_ncoord, cutoff)

   end subroutine test_dcndr_ammonia_erf


   subroutine test_dcndL_mb06_erf(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(erf_ncoord_type) :: erf_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "MB16-43", "06")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_erf_ncoord(erf_ncoord, mol, cutoff, rcov)
      call test_numsigma(error, mol, erf_ncoord, cutoff)

   end subroutine test_dcndL_mb06_erf


   subroutine test_dcndL_mb07_erf(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(erf_ncoord_type) :: erf_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "MB16-43", "07")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_erf_ncoord(erf_ncoord, mol, cutoff, rcov)
      call test_numsigma(error, mol, erf_ncoord, cutoff)

   end subroutine test_dcndL_mb07_erf


   subroutine test_dcndL_anthracene_erf(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(erf_ncoord_type) :: erf_ncoord
      real(wp), allocatable :: rcov(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "X23", "anthracene")

      allocate(rcov(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)

      call new_erf_ncoord(erf_ncoord, mol, cutoff, rcov)
      call test_numsigma(error, mol, erf_ncoord, cutoff)

   end subroutine test_dcndL_anthracene_erf


   !> ----------------------------------------------------
   !> Tests for electronegativity-weighted
   !> error-function-based CEH/GP3 coordination number
   !> using the Pyykko covalent radii and Pauling EN
   !> ----------------------------------------------------
   subroutine test_cn_mb01_erf_en(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(erf_en_ncoord_type) :: erf_en_ncoord
      real(wp), allocatable :: rcov(:)
      real(wp), allocatable :: en(:)

      real(wp), parameter :: cutoff = 30.0_wp
      real(wp), parameter :: ref(16) = [&
      & +6.4983577121808374E+00_wp, -9.2478042638094324E-02_wp, -3.1940402380604640E+00_wp, & 
      & -7.2646340531332720E-01_wp, -1.9217155297919475E+00_wp, +6.7170334215331562E-01_wp, & 
      & -6.4786183247720208E-01_wp, -2.5626344794490645E+00_wp, -3.2866057193592120E+00_wp, & 
      & +8.8640430680121118E-01_wp, +6.6508475327474137E-01_wp, -2.5117202580702962E+00_wp, & 
      & +4.2433691465974610E-01_wp, +2.9891464457177497E+00_wp, -2.7299215313690528E+00_wp, & 
      & +5.5384075617410602E+00_wp]

      call get_structure(mol, "MB16-43", "01")

      allocate(rcov(mol%nid), en(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)

      call new_erf_en_ncoord(erf_en_ncoord, mol, cutoff, rcov, en)
      call test_cn_gen(error, mol, erf_en_ncoord, cutoff, ref)

   end subroutine test_cn_mb01_erf_en


   subroutine test_cn_mb02_erf_en(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(erf_en_ncoord_type) :: erf_en_ncoord
      real(wp), allocatable :: rcov(:)
      real(wp), allocatable :: en(:)

      real(wp), parameter :: cutoff = 30.0_wp
      real(wp), parameter :: ref(16) = [&
      & -1.1223458428612135E-01_wp, -2.9400422417672871E+00_wp, -5.5856244559830737E-02_wp, & 
      & -4.1405542269891660E+00_wp, +4.9411871126140534E+00_wp, +6.9932942264825770E-01_wp, & 
      & -1.0382623569172383E-01_wp, -1.3791391815223347E-01_wp, +9.8051892904817461E-01_wp, & 
      & +6.6915562378679438E-01_wp, +5.4858837079167211E-01_wp, +6.7088630384026269E+00_wp, & 
      & -4.8169615059643691E+00_wp, -2.2159428443733206E-01_wp, -2.7932235946222728E-01_wp, & 
      & -1.7393368959812894E+00_wp]

      call get_structure(mol, "MB16-43", "02")

      allocate(rcov(mol%nid), en(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)

      call new_erf_en_ncoord(erf_en_ncoord, mol, cutoff, rcov, en)
      call test_cn_gen(error, mol,  erf_en_ncoord, cutoff, ref)

   end subroutine test_cn_mb02_erf_en


   subroutine test_cn_mb03_erf_en(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(erf_en_ncoord_type) :: erf_en_ncoord
      real(wp), allocatable :: rcov(:)
      real(wp), allocatable :: en(:)

      real(wp), parameter :: cutoff = 30.0_wp
      real(wp), parameter :: ref(16) = [&
      & -1.5712931068031073E+00_wp, -5.1253864860013216E+00_wp, +3.4440401067550655E-02_wp, & 
      & +5.5388805704334949E+00_wp, +6.0840392187133130E+00_wp, +1.7048440119191477E+00_wp, & 
      & -3.0579764303975088E+00_wp, -1.7733925517625523E-01_wp, +1.7710135398417912E-01_wp, & 
      & +3.1161198244520993E-01_wp, -4.7200177640622911E+00_wp, -1.7554540348635912E+00_wp, & 
      & -2.8757203214343594E+00_wp, +4.4516409543356017E-01_wp, +5.1226022769395847E+00_wp, & 
      & -1.3549651219760761E-01_wp]

      call get_structure(mol, "MB16-43", "03")

      allocate(rcov(mol%nid), en(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)

      call new_erf_en_ncoord(erf_en_ncoord, mol, cutoff, rcov, en)
      call test_cn_gen(error, mol, erf_en_ncoord, cutoff, ref)

   end subroutine test_cn_mb03_erf_en


   subroutine test_cn_acetic_erf_en(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(erf_en_ncoord_type) :: erf_en_ncoord
      real(wp), allocatable :: rcov(:)
      real(wp), allocatable :: en(:)

      real(wp), parameter :: cutoff = 30.0_wp
      real(wp), parameter :: ref(32) = [&
      & +1.2145797661349247E+00_wp, +1.2146193518354280E+00_wp, +1.2145352947774128E+00_wp, & 
      & +1.2145940886306013E+00_wp, +2.9497726378417549E-01_wp, +2.9496846852927938E-01_wp, & 
      & +2.9494914420211849E-01_wp, +2.9279771664404630E-01_wp, +2.9279399798719413E-01_wp, & 
      & +2.9278432351907768E-01_wp, +2.9278307996425129E-01_wp, +2.9233223687373783E-01_wp, & 
      & +2.9236558667241352E-01_wp, +2.9234812736512095E-01_wp, +2.9234223550974719E-01_wp, & 
      & +2.9497445400262740E-01_wp, +1.4789496078983553E+00_wp, +1.4789637312132844E+00_wp, & 
      & +1.4789530579350458E+00_wp, +1.4789353061979014E+00_wp, -6.0467609654193188E-01_wp, & 
      & -6.0465825267909690E-01_wp, -6.0466025694059378E-01_wp, -6.0466687608819336E-01_wp, & 
      & -1.8310642730391724E+00_wp, -1.8311341950273163E+00_wp, -1.8310560154433486E+00_wp, & 
      & -1.8311340900309643E+00_wp, -1.1378198331905931E+00_wp, -1.1378945508115248E+00_wp, & 
      & -1.1378778934968552E+00_wp, -1.1379045063871540E+00_wp]

      call get_structure(mol, "X23", "acetic")

      allocate(rcov(mol%nid), en(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)

      call new_erf_en_ncoord(erf_en_ncoord, mol, cutoff, rcov, en)
      call test_cn_gen(error, mol, erf_en_ncoord, cutoff, ref)

   end subroutine test_cn_acetic_erf_en


   subroutine test_dcndr_mb04_erf_en(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(erf_en_ncoord_type) :: erf_en_ncoord
      real(wp), allocatable :: rcov(:)
      real(wp), allocatable :: en(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "MB16-43", "04")

      allocate(rcov(mol%nid), en(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)

      call new_erf_en_ncoord(erf_en_ncoord, mol, cutoff, rcov, en)
      call test_numgrad(error, mol, erf_en_ncoord, cutoff)

   end subroutine test_dcndr_mb04_erf_en


   subroutine test_dcndr_mb05_erf_en(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(erf_en_ncoord_type) :: erf_en_ncoord
      real(wp), allocatable :: rcov(:)
      real(wp), allocatable :: en(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "MB16-43", "05")

      allocate(rcov(mol%nid), en(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)

      call new_erf_en_ncoord(erf_en_ncoord, mol, cutoff, rcov, en)
      call test_numgrad(error, mol, erf_en_ncoord, cutoff)

   end subroutine test_dcndr_mb05_erf_en


   subroutine test_dcndr_ammonia_erf_en(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(erf_en_ncoord_type) :: erf_en_ncoord
      real(wp), allocatable :: rcov(:)
      real(wp), allocatable :: en(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "X23", "ammonia")

      allocate(rcov(mol%nid), en(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)

      call new_erf_en_ncoord(erf_en_ncoord, mol, cutoff, rcov, en)
      call test_numgrad(error, mol, erf_en_ncoord, cutoff)

   end subroutine test_dcndr_ammonia_erf_en


   subroutine test_dcndL_mb06_erf_en(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(erf_en_ncoord_type) :: erf_en_ncoord
      real(wp), allocatable :: rcov(:)
      real(wp), allocatable :: en(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "MB16-43", "06")

      allocate(rcov(mol%nid), en(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)

      call new_erf_en_ncoord(erf_en_ncoord, mol, cutoff, rcov, en)
      call test_numsigma(error, mol, erf_en_ncoord, cutoff)

   end subroutine test_dcndL_mb06_erf_en


   subroutine test_dcndL_mb07_erf_en(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(erf_en_ncoord_type) :: erf_en_ncoord
      real(wp), allocatable :: rcov(:)
      real(wp), allocatable :: en(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "MB16-43", "07")

      allocate(rcov(mol%nid), en(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)

      call new_erf_en_ncoord(erf_en_ncoord, mol, cutoff, rcov, en)
      call test_numsigma(error, mol, erf_en_ncoord, cutoff)

   end subroutine test_dcndL_mb07_erf_en


   subroutine test_dcndL_anthracene_erf_en(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(erf_en_ncoord_type) :: erf_en_ncoord
      real(wp), allocatable :: rcov(:)
      real(wp), allocatable :: en(:)

      real(wp), parameter :: cutoff = 30.0_wp

      call get_structure(mol, "X23", "anthracene")

      allocate(rcov(mol%nid), en(mol%nid))
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)

      call new_erf_en_ncoord(erf_en_ncoord, mol, cutoff, rcov, en)
      call test_numsigma(error, mol, erf_en_ncoord, cutoff)

   end subroutine test_dcndL_anthracene_erf_en


end module test_ncoord
