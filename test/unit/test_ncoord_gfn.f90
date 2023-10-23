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
   use tblite_data_covrad_ceh, only : get_covalent_cehrad
   use tblite_data_paulingen_ceh, only : get_pauling_en_ceh
   use tblite_ncoord_gfn
   use tblite_ncoord_exp
   use tblite_ncoord_ceh_std
   use tblite_ncoord_ceh_en
   use tblite_ncoord_type !, only : get_coordination_number
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
      & new_unittest("cn-mb01_ceh_std", test_cn_mb01_ceh_std), &
      & new_unittest("cn-mb02_ceh_std", test_cn_mb02_ceh_std), &
      & new_unittest("cn-mb03_ceh_std", test_cn_mb03_ceh_std), &
      & new_unittest("cn-acetic_ceh_std", test_cn_acetic_ceh_std), &
      & new_unittest("dcndr-mb04_ceh_std", test_dcndr_mb04_ceh_std), &
      & new_unittest("dcndr-mb05_ceh_std", test_dcndr_mb05_ceh_std), &
      & new_unittest("dcndr-ammonia_ceh_std", test_dcndr_ammonia_ceh_std), &
      & new_unittest("dcndL-mb06_ceh_std", test_dcndL_mb06_ceh_std), &
      & new_unittest("dcndL-mb07_ceh_std", test_dcndL_mb07_ceh_std), &
      & new_unittest("dcndL-antracene_ceh_std", test_dcndL_anthracene_ceh_std), &
      & new_unittest("cn-mb01_ceh_en", test_cn_mb01_ceh_en), &
      & new_unittest("cn-mb02_ceh_en", test_cn_mb02_ceh_en), &
      & new_unittest("cn-mb03_ceh_en", test_cn_mb03_ceh_en), &
      & new_unittest("cn-acetic_ceh_en", test_cn_acetic_ceh_en), &
      & new_unittest("dcndr-mb04_ceh_en", test_dcndr_mb04_ceh_en), &
      & new_unittest("dcndr-mb05_ceh_en", test_dcndr_mb05_ceh_en), &
      & new_unittest("dcndr-ammonia_ceh_en", test_dcndr_ammonia_ceh_en), &
      & new_unittest("dcndL-mb06_ceh_en", test_dcndL_mb06_ceh_en), &
      & new_unittest("dcndL-mb07_ceh_en", test_dcndL_mb07_ceh_en), &
      & new_unittest("dcndL-antracene_ceh_en", test_dcndL_anthracene_ceh_en) &
      & ]

end subroutine collect_ncoord_gfn


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
!> ----------------------------------------------------
subroutine test_cn_mb01_ceh_std(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(ceh_std_ncoord_type) :: ceh_std_ncoord
   real(wp), allocatable :: rcov(:)

   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), parameter :: ref(16) = [&
      & 2.99602033361644E-1_wp, 2.94632924190307E-1_wp, 8.39338029333359E-1_wp, &
      & 2.72853670043490E-1_wp, 4.65294250591486E-1_wp, 4.96615809565662E-1_wp, &
      & 2.57299033286322E-1_wp, 1.02366732559721E+0_wp, 1.37742318797612E+0_wp, &
      & 5.07478856979855E-1_wp, 4.92545892335530E-1_wp, 5.12888339327036E-1_wp, &
      & 1.04817884216878E+0_wp, 1.62461738249920E+0_wp, 1.69546464698267E+0_wp, &
      & 1.67615854937172E+0_wp]
  
   call get_structure(mol, "MB16-43", "01")
 
   allocate(rcov(mol%nid))
   rcov(:) = get_covalent_cehrad(mol%num)

   call new_ceh_std_ncoord(ceh_std_ncoord, mol, cutoff, rcov)
   call test_cn_gen(error, mol, ceh_std_ncoord, cutoff, ref)

end subroutine test_cn_mb01_ceh_std


subroutine test_cn_mb02_ceh_std(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(ceh_std_ncoord_type) :: ceh_std_ncoord
   real(wp), allocatable :: rcov(:)

   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), parameter :: ref(16) = [&
      & 2.64224330027233E-1_wp, 9.45142882305708E-1_wp, 1.17350403979181E+0_wp, &
      & 1.06367969241646E+0_wp, 9.31411718959526E-1_wp, 4.78516391083669E-1_wp, &
      & 2.86744968440476E-1_wp, 2.92618519639304E-1_wp, 1.83617445072684E+0_wp, &
      & 4.64315766512035E-1_wp, 1.09500586721714E+0_wp, 3.69669889756518E-1_wp, &
      & 4.91050486360982E-1_wp, 4.25925275638909E-1_wp, 3.89724987256085E-1_wp, &
      & 8.17338263052412E-1_wp]

   call get_structure(mol, "MB16-43", "02")

   allocate(rcov(mol%nid))
   rcov(:) = get_covalent_cehrad(mol%num)

   call new_ceh_std_ncoord(ceh_std_ncoord, mol, cutoff, rcov)
   call test_cn_gen(error, mol, ceh_std_ncoord, cutoff, ref)

end subroutine test_cn_mb02_ceh_std


subroutine test_cn_mb03_ceh_std(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(ceh_std_ncoord_type) :: ceh_std_ncoord
   real(wp), allocatable :: rcov(:)

   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), parameter :: ref(16) = [&
      & 1.70130387316316E+0_wp, 9.39803385893297E-1_wp, 4.35943196783110E-1_wp, &
      & 3.63546168569012E-1_wp, 1.08989121337887E+0_wp, 6.36196297091134E-1_wp, &
      & 1.50314678323341E+0_wp, 4.34744185826058E-1_wp, 4.38356709382285E-1_wp, &
      & 4.51355458073888E-1_wp, 4.77403137241262E-1_wp, 9.60062266192319E-1_wp, &
      & 1.56766012525407E+0_wp, 5.04574836212309E-1_wp, 3.53605176448097E-1_wp, &
      & 4.12701662714700E-1_wp]

   call get_structure(mol, "MB16-43", "03")

   allocate(rcov(mol%nid))
   rcov(:) = get_covalent_cehrad(mol%num)

   call new_ceh_std_ncoord(ceh_std_ncoord, mol, cutoff, rcov)
   call test_cn_gen(error, mol, ceh_std_ncoord, cutoff, ref)

end subroutine test_cn_mb03_ceh_std


subroutine test_cn_acetic_ceh_std(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(ceh_std_ncoord_type) :: ceh_std_ncoord
   real(wp), allocatable :: rcov(:)

   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), parameter :: ref(32) = [&
      & 4.04183870924080E-1_wp, 4.04324894002753E-1_wp, 4.04226264330234E-1_wp, &
      & 4.04323764557260E-1_wp, 4.47560166869845E-1_wp, 4.47501204622850E-1_wp, &
      & 4.47379658202106E-1_wp, 4.39737194647952E-1_wp, 4.39700299538792E-1_wp, &
      & 4.39679475831586E-1_wp, 4.39642582687677E-1_wp, 4.41022406717090E-1_wp, &
      & 4.41210571027497E-1_wp, 4.41155394171328E-1_wp, 4.41086851760964E-1_wp, &
      & 4.47532374597696E-1_wp, 1.85712514293316E+0_wp, 1.85721806394664E+0_wp, &
      & 1.85706937710972E+0_wp, 1.85711765923651E+0_wp, 1.85060931071132E+0_wp, &
      & 1.85073101217519E+0_wp, 1.85045570461204E+0_wp, 1.85057980248076E+0_wp, &
      & 1.02678810459768E+0_wp, 1.02686905781610E+0_wp, 1.02672544845891E+0_wp, &
      & 1.02686783987009E+0_wp, 7.20763207476270E-1_wp, 7.20893527835198E-1_wp, &
      & 7.20864010565967E-1_wp, 7.20794506792329E-1_wp]

   call get_structure(mol, "X23", "acetic")

   allocate(rcov(mol%nid))
   rcov(:) = get_covalent_cehrad(mol%num)

   call new_ceh_std_ncoord(ceh_std_ncoord, mol, cutoff, rcov)
   call test_cn_gen(error, mol, ceh_std_ncoord, cutoff, ref)

end subroutine test_cn_acetic_ceh_std


subroutine test_dcndr_mb04_ceh_std(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(ceh_std_ncoord_type) :: ceh_std_ncoord
   real(wp), allocatable :: rcov(:)

   real(wp), parameter :: cutoff = 30.0_wp

   call get_structure(mol, "MB16-43", "04")

   allocate(rcov(mol%nid))
   rcov(:) = get_covalent_cehrad(mol%num)

   call new_ceh_std_ncoord(ceh_std_ncoord, mol, cutoff, rcov)
   call test_numgrad(error, mol, ceh_std_ncoord, cutoff)

end subroutine test_dcndr_mb04_ceh_std


subroutine test_dcndr_mb05_ceh_std(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(ceh_std_ncoord_type) :: ceh_std_ncoord
   real(wp), allocatable :: rcov(:)

   real(wp), parameter :: cutoff = 30.0_wp

   call get_structure(mol, "MB16-43", "05")

   allocate(rcov(mol%nid))
   rcov(:) = get_covalent_cehrad(mol%num)

   call new_ceh_std_ncoord(ceh_std_ncoord, mol, cutoff, rcov)
   call test_numgrad(error, mol, ceh_std_ncoord, cutoff)

end subroutine test_dcndr_mb05_ceh_std


subroutine test_dcndr_ammonia_ceh_std(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(ceh_std_ncoord_type) :: ceh_std_ncoord
   real(wp), allocatable :: rcov(:)

   real(wp), parameter :: cutoff = 30.0_wp

   call get_structure(mol, "X23", "ammonia")

   allocate(rcov(mol%nid))
   rcov(:) = get_covalent_cehrad(mol%num)

   call new_ceh_std_ncoord(ceh_std_ncoord, mol, cutoff, rcov)
   call test_numgrad(error, mol, ceh_std_ncoord, cutoff)

end subroutine test_dcndr_ammonia_ceh_std


subroutine test_dcndL_mb06_ceh_std(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(ceh_std_ncoord_type) :: ceh_std_ncoord
   real(wp), allocatable :: rcov(:)

   real(wp), parameter :: cutoff = 30.0_wp

   call get_structure(mol, "MB16-43", "06")
      
   allocate(rcov(mol%nid))
   rcov(:) = get_covalent_cehrad(mol%num)

   call new_ceh_std_ncoord(ceh_std_ncoord, mol, cutoff, rcov)
   call test_numsigma(error, mol, ceh_std_ncoord, cutoff)

end subroutine test_dcndL_mb06_ceh_std


subroutine test_dcndL_mb07_ceh_std(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(ceh_std_ncoord_type) :: ceh_std_ncoord
   real(wp), allocatable :: rcov(:)

   real(wp), parameter :: cutoff = 30.0_wp

   call get_structure(mol, "MB16-43", "07")
         
   allocate(rcov(mol%nid))
   rcov(:) = get_covalent_cehrad(mol%num)

   call new_ceh_std_ncoord(ceh_std_ncoord, mol, cutoff, rcov)
   call test_numsigma(error, mol, ceh_std_ncoord, cutoff)

end subroutine test_dcndL_mb07_ceh_std


subroutine test_dcndL_anthracene_ceh_std(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(ceh_std_ncoord_type) :: ceh_std_ncoord
   real(wp), allocatable :: rcov(:)

   real(wp), parameter :: cutoff = 30.0_wp

   call get_structure(mol, "X23", "anthracene")
         
   allocate(rcov(mol%nid))
   rcov(:) = get_covalent_cehrad(mol%num)

   call new_ceh_std_ncoord(ceh_std_ncoord, mol, cutoff, rcov)
   call test_numsigma(error, mol, ceh_std_ncoord, cutoff)

end subroutine test_dcndL_anthracene_ceh_std


!> ----------------------------------------------------
!> Tests for electronegativity-weighted  
!> error-function-based CEH/GP3 coordination number 
!> ----------------------------------------------------
subroutine test_cn_mb01_ceh_en(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(ceh_en_ncoord_type) :: ceh_en_ncoord
   real(wp), allocatable :: rcov(:)
   real(wp), allocatable :: en(:)

   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), parameter :: ref(16) = [&
      & 1.64837690033030E-1_wp, -1.18442609679655E-2_wp, -3.49242315441712E-1_wp, &
      &-1.20043854727766E-2_wp, -2.26741511076709E-1_wp,  1.04810774566835E-1_wp, &
      &-1.37200255399630E-2_wp, -3.95528202678488E-1_wp, -3.16687108761406E-1_wp, &
      & 1.57870443276754E-1_wp,  1.03948705629303E-1_wp, -1.75836699870889E-1_wp, &
      & 8.91343409428138E-2_wp,  5.62967900802320E-1_wp, -3.26991579947931E-1_wp, &
      & 6.45026234506784E-1_wp]
  
   call get_structure(mol, "MB16-43", "01")
 
   allocate(rcov(mol%nid), en(mol%nid))
   rcov(:) = get_covalent_cehrad(mol%num)
   en(:) = get_pauling_en_ceh(mol%num)

   call new_ceh_en_ncoord(ceh_en_ncoord, mol, cutoff, rcov, en)
   call test_cn_gen(error, mol, ceh_en_ncoord, cutoff, ref)

end subroutine test_cn_mb01_ceh_en


subroutine test_cn_mb02_ceh_en(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(ceh_en_ncoord_type) :: ceh_en_ncoord
   real(wp), allocatable :: rcov(:)
   real(wp), allocatable :: en(:)

   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), parameter :: ref(16) = [&
      &-1.06213013717222E-2_wp, -1.59361567146995E-1_wp,  3.41448491614630E-2_wp, &
      &-3.68480554591165E-1_wp,  3.67376861584839E-1_wp,  1.49084995118337E-1_wp, &
      &-1.15254384363943E-2_wp, -1.17632512371385E-2_wp,  1.51678937490452E-1_wp, &
      & 1.44627771030999E-1_wp,  7.54651838437863E-2_wp,  2.20643899491151E-1_wp, &
      &-3.42942435301755E-1_wp, -3.21049158305882E-2_wp, -2.48842854873811E-2_wp, &
      &-1.81338748317887E-1_wp]

   call get_structure(mol, "MB16-43", "02")

   allocate(rcov(mol%nid), en(mol%nid))
   rcov(:) = get_covalent_cehrad(mol%num)
   en(:) = get_pauling_en_ceh(mol%num)
   
   call new_ceh_en_ncoord(ceh_en_ncoord, mol, cutoff, rcov, en)
   call test_cn_gen(error, mol,  ceh_en_ncoord, cutoff, ref)

end subroutine test_cn_mb02_ceh_en


subroutine test_cn_mb03_ceh_en(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(ceh_en_ncoord_type) :: ceh_en_ncoord
   real(wp), allocatable :: rcov(:)
   real(wp), allocatable :: en(:)

   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), parameter :: ref(16) = [&
      &-5.99035988018810E-2_wp, -4.28671245878299E-1_wp,  3.83357862396699E-2_wp, &
      & 1.85747489975467E-1_wp,  4.27210508080517E-1_wp,  1.04606875688658E-1_wp, &
      &-1.56423676461829E-1_wp,  3.82214556265731E-2_wp,  3.85369480640030E-2_wp, &
      & 3.96928112487448E-2_wp, -3.34027743669695E-1_wp, -1.28376323193933E-1_wp, &
      &-1.91517181561922E-1_wp,  1.57195970331334E-1_wp,  2.34735906485441E-1_wp, &
      & 3.46360178271516E-2_wp]

   call get_structure(mol, "MB16-43", "03")

   allocate(rcov(mol%nid), en(mol%nid))
   rcov(:) = get_covalent_cehrad(mol%num)
   en(:) = get_pauling_en_ceh(mol%num)

   call new_ceh_en_ncoord(ceh_en_ncoord, mol, cutoff, rcov, en)
   call test_cn_gen(error, mol, ceh_en_ncoord, cutoff, ref)

end subroutine test_cn_mb03_ceh_en


subroutine test_cn_acetic_ceh_en(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(ceh_en_ncoord_type) :: ceh_en_ncoord
   real(wp), allocatable :: rcov(:)
   real(wp), allocatable :: en(:)

   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), parameter :: ref(32) = [&

      & 1.25896629244905E-1_wp,  1.25940569890756E-1_wp,  1.25909831705899E-1_wp, &
      & 1.25940218003629E-1_wp,  3.93583062789118E-2_wp,  3.93531211566207E-2_wp, &
      & 3.93424324012873E-2_wp,  3.86703563132633E-2_wp,  3.86671117684880E-2_wp, &
      & 3.86652805379547E-2_wp,  3.86620361660030E-2_wp,  3.87833774751703E-2_wp, &
      & 3.87999245878945E-2_wp,  3.87950723517993E-2_wp,  3.87890447528980E-2_wp, &
      & 3.93558622348648E-2_wp,  2.99092397546000E-1_wp,  2.99106816436381E-1_wp, &
      & 2.99090746382539E-1_wp,  2.99084366202306E-1_wp, -1.16162802460569E-1_wp, &
      &-1.16170706669792E-1_wp, -1.16153544831317E-1_wp, -1.16157705381245E-1_wp, &
      &-2.64504473005920E-1_wp, -2.64534916941950E-1_wp, -2.64494669249799E-1_wp, &
      &-2.64534916941947E-1_wp, -1.61132693878179E-1_wp, -1.61161568341170E-1_wp, &
      &-1.61155501185592E-1_wp, -1.61140002550092E-1_wp]

   call get_structure(mol, "X23", "acetic")

   allocate(rcov(mol%nid), en(mol%nid))
   rcov(:) = get_covalent_cehrad(mol%num)
   en(:) = get_pauling_en_ceh(mol%num)

   call new_ceh_en_ncoord(ceh_en_ncoord, mol, cutoff, rcov, en)
   call test_cn_gen(error, mol, ceh_en_ncoord, cutoff, ref)

end subroutine test_cn_acetic_ceh_en


subroutine test_dcndr_mb04_ceh_en(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(ceh_en_ncoord_type) :: ceh_en_ncoord
   real(wp), allocatable :: rcov(:)
   real(wp), allocatable :: en(:)

   real(wp), parameter :: cutoff = 30.0_wp

   call get_structure(mol, "MB16-43", "04")

   allocate(rcov(mol%nid), en(mol%nid))
   rcov(:) = get_covalent_cehrad(mol%num)
   en(:) = get_pauling_en_ceh(mol%num)

   call new_ceh_en_ncoord(ceh_en_ncoord, mol, cutoff, rcov, en)
   call test_numgrad(error, mol, ceh_en_ncoord, cutoff)

end subroutine test_dcndr_mb04_ceh_en


subroutine test_dcndr_mb05_ceh_en(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(ceh_en_ncoord_type) :: ceh_en_ncoord
   real(wp), allocatable :: rcov(:)
   real(wp), allocatable :: en(:)

   real(wp), parameter :: cutoff = 30.0_wp

   call get_structure(mol, "MB16-43", "05")

   allocate(rcov(mol%nid), en(mol%nid))
   rcov(:) = get_covalent_cehrad(mol%num)
   en(:) = get_pauling_en_ceh(mol%num)

   call new_ceh_en_ncoord(ceh_en_ncoord, mol, cutoff, rcov, en)
   call test_numgrad(error, mol, ceh_en_ncoord, cutoff)

end subroutine test_dcndr_mb05_ceh_en


subroutine test_dcndr_ammonia_ceh_en(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(ceh_en_ncoord_type) :: ceh_en_ncoord
   real(wp), allocatable :: rcov(:)
   real(wp), allocatable :: en(:)

   real(wp), parameter :: cutoff = 30.0_wp

   call get_structure(mol, "X23", "ammonia")

   allocate(rcov(mol%nid), en(mol%nid))
   rcov(:) = get_covalent_cehrad(mol%num)
   en(:) = get_pauling_en_ceh(mol%num)

   call new_ceh_en_ncoord(ceh_en_ncoord, mol, cutoff, rcov, en)
   call test_numgrad(error, mol, ceh_en_ncoord, cutoff)

end subroutine test_dcndr_ammonia_ceh_en


subroutine test_dcndL_mb06_ceh_en(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(ceh_en_ncoord_type) :: ceh_en_ncoord
   real(wp), allocatable :: rcov(:)
   real(wp), allocatable :: en(:)

   real(wp), parameter :: cutoff = 30.0_wp

   call get_structure(mol, "MB16-43", "06")
      
   allocate(rcov(mol%nid), en(mol%nid))
   rcov(:) = get_covalent_cehrad(mol%num)
   en(:) = get_pauling_en_ceh(mol%num)

   call new_ceh_en_ncoord(ceh_en_ncoord, mol, cutoff, rcov, en)
   call test_numsigma(error, mol, ceh_en_ncoord, cutoff)

end subroutine test_dcndL_mb06_ceh_en


subroutine test_dcndL_mb07_ceh_en(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(ceh_en_ncoord_type) :: ceh_en_ncoord
   real(wp), allocatable :: rcov(:)
   real(wp), allocatable :: en(:)

   real(wp), parameter :: cutoff = 30.0_wp

   call get_structure(mol, "MB16-43", "07")
         
   allocate(rcov(mol%nid), en(mol%nid))
   rcov(:) = get_covalent_cehrad(mol%num)
   en(:) = get_pauling_en_ceh(mol%num)

   call new_ceh_en_ncoord(ceh_en_ncoord, mol, cutoff, rcov, en)
   call test_numsigma(error, mol, ceh_en_ncoord, cutoff)

end subroutine test_dcndL_mb07_ceh_en


subroutine test_dcndL_anthracene_ceh_en(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(ceh_en_ncoord_type) :: ceh_en_ncoord
   real(wp), allocatable :: rcov(:)
   real(wp), allocatable :: en(:)

   real(wp), parameter :: cutoff = 30.0_wp

   call get_structure(mol, "X23", "anthracene")
         
   allocate(rcov(mol%nid), en(mol%nid))
   rcov(:) = get_covalent_rad(mol%num)
   en(:) = get_pauling_en_ceh(mol%num)

   call new_ceh_en_ncoord(ceh_en_ncoord, mol, cutoff, rcov, en)
   call test_numsigma(error, mol, ceh_en_ncoord, cutoff)

end subroutine test_dcndL_anthracene_ceh_en


end module test_ncoord_gfn
