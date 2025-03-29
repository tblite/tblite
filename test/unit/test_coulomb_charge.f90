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

module test_coulomb_charge
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_basis_type
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_container_cache, only : container_cache
   use tblite_coulomb_charge, only : coulomb_charge_type, effective_coulomb, &
      & new_effective_coulomb, harmonic_average, arithmetic_average, &
      & gamma_coulomb, new_gamma_coulomb
   use tblite_scf, only: new_potential, potential_type
   use tblite_ceh_ceh, only : get_effective_qat
   use tblite_ncoord_erf_en
   use tblite_ncoord_type, only : get_coordination_number
   use tblite_cutoff, only : get_lattice_points
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   implicit none
   private

   public :: collect_coulomb_charge

   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   abstract interface
      subroutine coulomb_maker(coulomb, mol, shell)
         import :: coulomb_charge_type, structure_type
         class(coulomb_charge_type), allocatable, intent(out) :: coulomb
         type(structure_type), intent(in) :: mol
         logical, intent(in) :: shell
      end subroutine coulomb_maker
   end interface

   abstract interface
      subroutine charge_maker(wfn, mol, nshell)
         import :: wavefunction_type, structure_type
         class(wavefunction_type), intent(inout) :: wfn
         type(structure_type), intent(in) :: mol
         integer, optional, intent(in) :: nshell(:)
      end subroutine charge_maker
   end interface

contains


!> Collect all exported unit tests
subroutine collect_coulomb_charge(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("energy-atom-e1", test_e_effective_m01), &
      new_unittest("energy-atom-e2", test_e_effective_m02), &
      new_unittest("energy-shell-e1", test_e_effective_m07), &
      new_unittest("energy-atom-pbc-e2", test_e_effective_oxacb), &
      new_unittest("energy-atom-sc-e2", test_e_effective_oxacb_sc), &
      new_unittest("energy-atom-g1", test_e_effective_m10), &
      new_unittest("energy-shell-g2", test_e_effective_m13), &
      new_unittest("gradient-atom-e1", test_g_effective_m03), &
      new_unittest("gradient-atom-e2", test_g_effective_m04), &
      new_unittest("gradient-shell-e1", test_g_effective_m08), &
      new_unittest("gradient-atom-pbc-e2", test_g_effective_co2), &
      new_unittest("gradient-atom-g1", test_g_effective_m11), &
      new_unittest("gradient-shell-g2", test_g_effective_m14), &
      new_unittest("gradient-atom-pbc-g1", test_g_effective_urea), &
      new_unittest("sigma-atom-e1", test_s_effective_m05), &
      new_unittest("sigma-atom-e2", test_s_effective_m06), &
      new_unittest("sigma-shell-e1", test_s_effective_m09), &
      new_unittest("sigma-atom-pbc-e2", test_s_effective_ammonia), &
      new_unittest("sigma-atom-g1", test_s_effective_m12), &
      new_unittest("sigma-shell-g2", test_s_effective_m15), &
      new_unittest("sigma-atom-pbc-g2", test_s_effective_pyrazine), &
      new_unittest("potential-gradient-lih-effceh", test_ceh_potgrad_lih), &
      new_unittest("potential-gradient-mb15-effceh", test_ceh_potgrad_mb15), &
      new_unittest("potential-gradient-mb16-shell-effceh", test_ceh_potgrad_mb16), &
      new_unittest("potential-sigma-lih-effceh", test_ceh_potsigma_lih), &
      new_unittest("potential-sigma-co2-effceh", test_ceh_potsigma_co2), &
      new_unittest("potential-sigma-mb05-effceh", test_ceh_potsigma_mb05), &
      new_unittest("potential-sigma-mb17-shell-effceh", test_ceh_potsigma_mb17) &
      ]

end subroutine collect_coulomb_charge


!> Factory to setup the CEH basis set for testing of the potential (gradient)
subroutine make_basis(bas, mol, ng)
   type(basis_type), intent(out) :: bas
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: ng

   integer, parameter :: nsh(20) = [&
   & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3]
   integer, parameter :: lsh(3, 20) = reshape([&
   & 0, 0, 0,  0, 0, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0, &
   & 0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 2, &
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, &
   & 0, 1, 0,  0, 1, 2], shape(lsh))

   integer, parameter :: pqn(3, 20) = reshape([&
   & 1, 0, 0,  1, 0, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0, &
   & 2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  3, 3, 0,  3, 3, 3, &
   & 3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3, &
   & 4, 4, 0,  4, 4, 3], shape(pqn))

   real(wp), parameter :: zeta(3, 20) = reshape([&
   & 1.23363166_wp, 0.00000000_wp, 0.00000000_wp, 2.27004605_wp, 0.00000000_wp, 0.00000000_wp, &
   & 0.86185456_wp, 1.42017184_wp, 0.00000000_wp, 1.76817995_wp, 1.44095844_wp, 0.00000000_wp, &
   & 2.06339837_wp, 1.52051807_wp, 0.00000000_wp, 2.56058582_wp, 1.86484737_wp, 0.00000000_wp, &
   & 2.71233631_wp, 2.19848968_wp, 0.00000000_wp, 3.21585650_wp, 2.41309737_wp, 0.00000000_wp, &
   & 3.82146807_wp, 2.63063636_wp, 0.00000000_wp, 4.62721228_wp, 2.53599954_wp, 0.00000000_wp, &
   & 0.93221172_wp, 1.55333839_wp, 0.00000000_wp, 1.77220557_wp, 1.59942632_wp, 2.98596647_wp, &
   & 2.26040231_wp, 1.78718151_wp, 2.00990188_wp, 1.85259089_wp, 1.81733349_wp, 1.65269988_wp, &
   & 2.65701241_wp, 2.03189759_wp, 2.03883661_wp, 2.60609998_wp, 2.16530440_wp, 2.41888232_wp, &
   & 2.78818934_wp, 2.24732894_wp, 1.99081182_wp, 2.55424399_wp, 2.20946190_wp, 1.93619550_wp, &
   & 1.73713827_wp, 1.33788617_wp, 0.00000000_wp, 2.47982574_wp, 1.07250770_wp, 2.11920764_wp],&
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


!> Factory to create electrostatic objects based on GFN1-xTB values
subroutine make_coulomb_e1(coulomb, mol, shell)

   !> New electrostatic object
   class(coulomb_charge_type), allocatable, intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return a shell resolved object
   logical, intent(in) :: shell

   real(wp), parameter :: hubbard_parameter(20) = [&
      & 0.470099_wp, 1.441379_wp, 0.205342_wp, 0.274022_wp, 0.340530_wp, &
      & 0.479988_wp, 0.476106_wp, 0.583349_wp, 0.788194_wp, 0.612878_wp, &
      & 0.165908_wp, 0.354151_wp, 0.221658_wp, 0.438331_wp, 0.798319_wp, &
      & 0.643959_wp, 0.519712_wp, 0.529906_wp, 0.114358_wp, 0.134187_wp]
   integer, parameter :: shell_count(20) = [&
      & 2, 1, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 3, 3, 3, 3, 3, 3, 2, 3]
   real(wp), parameter :: shell_scale(3, 20) = reshape([&
      & 0.0_wp, 0.0000000_wp, 0.0000000_wp,  0.0_wp, 0.0000000_wp, 0.0000000_wp, &
      & 0.0_wp,-0.0772012_wp, 0.0000000_wp,  0.0_wp, 0.1113005_wp, 0.0000000_wp, &
      & 0.0_wp, 0.0165643_wp, 0.0000000_wp,  0.0_wp,-0.0471181_wp, 0.0000000_wp, &
      & 0.0_wp, 0.0315090_wp, 0.0000000_wp,  0.0_wp, 0.0374608_wp, 0.0000000_wp, &
      & 0.0_wp,-0.0827352_wp, 0.0000000_wp,  0.0_wp,-0.3892542_wp, 0.0000000_wp, &
      & 0.0_wp,-0.3004391_wp, 0.0000000_wp,  0.0_wp, 0.0674819_wp, 0.0000000_wp, &
      & 0.0_wp, 0.0503564_wp, 0.0000000_wp,  0.0_wp,-0.5925834_wp, 0.0000000_wp, &
      & 0.0_wp,-0.2530875_wp, 0.0000000_wp,  0.0_wp,-0.1678147_wp, 0.0000000_wp, &
      & 0.0_wp,-0.4481841_wp, 0.0000000_wp,  0.0_wp,-0.1450000_wp, 0.0000000_wp, &
      & 0.0_wp,-0.5332978_wp, 0.0000000_wp,  0.0_wp, 1.1522018_wp, 0.0000000_wp],&
      & shape(shell_scale)) + 1.0_wp
   real(wp), parameter :: gexp = 2.0_wp
   real(wp), allocatable :: hubbard(:, :)
   integer :: isp, izp, ish
   type(effective_coulomb), allocatable :: tmp

   allocate(tmp)
   if (shell) then
      allocate(hubbard(3, mol%nid))
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, shell_count(izp)
            hubbard(ish, isp) = hubbard_parameter(izp) * shell_scale(ish, izp)
         end do
      end do
      call new_effective_coulomb(tmp, mol, gexp, hubbard, harmonic_average, &
         & shell_count(mol%num))
   else
      hubbard = reshape(hubbard_parameter(mol%num), [1, mol%nid])
      call new_effective_coulomb(tmp, mol, gexp, hubbard, harmonic_average)
   end if
   call move_alloc(tmp, coulomb)

end subroutine make_coulomb_e1

!> Factory to create electrostatic objects based on GFN2-xTB values
subroutine make_coulomb_e2(coulomb, mol, shell)

   !> New electrostatic object
   class(coulomb_charge_type), allocatable, intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return a shell resolved object
   logical, intent(in) :: shell

   real(wp), parameter :: hubbard_parameter(20) = [&
      & 0.405771_wp, 0.642029_wp, 0.245006_wp, 0.684789_wp, 0.513556_wp, &
      & 0.538015_wp, 0.461493_wp, 0.451896_wp, 0.531518_wp, 0.850000_wp, &
      & 0.271056_wp, 0.344822_wp, 0.364801_wp, 0.720000_wp, 0.297739_wp, &
      & 0.339971_wp, 0.248514_wp, 0.502376_wp, 0.247602_wp, 0.320378_wp]
   integer, parameter :: shell_count(20) = [&
      & 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3]
   real(wp), parameter :: shell_scale(3, 20) = reshape([&
      & 0.0_wp, 0.0000000_wp, 0.0000000_wp, 0.0_wp, 0.0000000_wp, 0.0000000_wp, &
      & 0.0_wp, 0.1972612_wp, 0.0000000_wp, 0.0_wp, 0.9658467_wp, 0.0000000_wp, &
      & 0.0_wp, 0.3994080_wp, 0.0000000_wp, 0.0_wp, 0.1056358_wp, 0.0000000_wp, &
      & 0.0_wp, 0.1164892_wp, 0.0000000_wp, 0.0_wp, 0.1497020_wp, 0.0000000_wp, &
      & 0.0_wp, 0.1677376_wp, 0.0000000_wp, 0.0_wp, 0.1190576_wp,-0.3200000_wp, &
      & 0.0_wp, 0.1018894_wp, 0.0000000_wp, 0.0_wp, 1.4000000_wp,-0.0500000_wp, &
      & 0.0_wp,-0.0603699_wp, 0.2000000_wp, 0.0_wp,-0.5580042_wp,-0.2300000_wp, &
      & 0.0_wp,-0.1558060_wp,-0.3500000_wp, 0.0_wp,-0.1085866_wp,-0.2500000_wp, &
      & 0.0_wp, 0.4989400_wp, 0.5000000_wp, 0.0_wp,-0.0461133_wp,-0.0100000_wp, &
      & 0.0_wp, 0.3483655_wp, 0.0000000_wp, 0.0_wp, 1.5000000_wp,-0.2500000_wp],&
      & shape(shell_scale)) + 1.0_wp
   real(wp), parameter :: gexp = 2.0_wp
   real(wp), allocatable :: hubbard(:, :)
   integer :: isp, izp, ish
   type(effective_coulomb), allocatable :: tmp

   allocate(tmp)
   if (shell) then
      allocate(hubbard(3, mol%nid))
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, shell_count(izp)
            hubbard(ish, isp) = hubbard_parameter(izp) * shell_scale(ish, izp)
         end do
      end do
      call new_effective_coulomb(tmp, mol, gexp, hubbard, arithmetic_average, &
         & shell_count(mol%num))
   else
      hubbard = reshape(hubbard_parameter(mol%num), [1, mol%nid])
      call new_effective_coulomb(tmp, mol, gexp, hubbard, arithmetic_average)
   end if
   call move_alloc(tmp, coulomb)

end subroutine make_coulomb_e2

!> Factory to create electrostatic objects based on CEH values
subroutine make_coulomb_eceh(coulomb, mol, shell)

   !> New electrostatic object
   class(coulomb_charge_type), allocatable, intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return a shell resolved object
   logical, intent(in) :: shell

   real(wp), parameter :: hubbard_parameter(20) = [&
      & 0.47259288_wp, 0.92203391_wp, 0.17452888_wp, 0.25700733_wp, &
      & 0.33949086_wp, 0.42195412_wp, 0.50438193_wp, 0.58691863_wp, &
      & 0.66931351_wp, 0.75191607_wp, 0.17964105_wp, 0.22157276_wp, &
      & 0.26348578_wp, 0.30539645_wp, 0.34734014_wp, 0.38924725_wp, &
      & 0.43115670_wp, 0.47308269_wp, 0.17105469_wp, 0.20276244_wp]
   integer, parameter :: shell_count(20) = [&
      & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3]

   real(wp), parameter :: gexp = 1.0_wp
   real(wp), allocatable :: hubbard(:, :)
   integer :: isp, izp, ish
   type(effective_coulomb), allocatable :: tmp

   allocate(tmp)
   if (shell) then
      ! If shell-resolved, we use the atomic parameter for each shell
      allocate(hubbard(3, mol%nid))
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, shell_count(izp)
            hubbard(ish, isp) = hubbard_parameter(izp)
         end do
      end do
      call new_effective_coulomb(tmp, mol, gexp, hubbard, arithmetic_average, &
         & shell_count(mol%num))
   else
      hubbard = reshape(hubbard_parameter(mol%num), [1, mol%nid])
      call new_effective_coulomb(tmp, mol, gexp, hubbard, arithmetic_average)
   end if
   call move_alloc(tmp, coulomb)

end subroutine make_coulomb_eceh

!> Factory to create electrostatic objects based on GFN1-xTB values
subroutine make_coulomb_g1(coulomb, mol, shell)

   !> New electrostatic object
   class(coulomb_charge_type), allocatable, intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return a shell resolved object
   logical, intent(in) :: shell

   real(wp), parameter :: hubbard_parameter(20) = [&
      & 0.470099_wp, 1.441379_wp, 0.205342_wp, 0.274022_wp, 0.340530_wp, &
      & 0.479988_wp, 0.476106_wp, 0.583349_wp, 0.788194_wp, 0.612878_wp, &
      & 0.165908_wp, 0.354151_wp, 0.221658_wp, 0.438331_wp, 0.798319_wp, &
      & 0.643959_wp, 0.519712_wp, 0.529906_wp, 0.114358_wp, 0.134187_wp]
   integer, parameter :: shell_count(20) = [&
      & 2, 1, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 3, 3, 3, 3, 3, 3, 2, 3]
   real(wp), parameter :: shell_scale(3, 20) = reshape([&
      & 0.0_wp, 0.0000000_wp, 0.0000000_wp,  0.0_wp, 0.0000000_wp, 0.0000000_wp, &
      & 0.0_wp,-0.0772012_wp, 0.0000000_wp,  0.0_wp, 0.1113005_wp, 0.0000000_wp, &
      & 0.0_wp, 0.0165643_wp, 0.0000000_wp,  0.0_wp,-0.0471181_wp, 0.0000000_wp, &
      & 0.0_wp, 0.0315090_wp, 0.0000000_wp,  0.0_wp, 0.0374608_wp, 0.0000000_wp, &
      & 0.0_wp,-0.0827352_wp, 0.0000000_wp,  0.0_wp,-0.3892542_wp, 0.0000000_wp, &
      & 0.0_wp,-0.3004391_wp, 0.0000000_wp,  0.0_wp, 0.0674819_wp, 0.0000000_wp, &
      & 0.0_wp, 0.0503564_wp, 0.0000000_wp,  0.0_wp,-0.5925834_wp, 0.0000000_wp, &
      & 0.0_wp,-0.2530875_wp, 0.0000000_wp,  0.0_wp,-0.1678147_wp, 0.0000000_wp, &
      & 0.0_wp,-0.4481841_wp, 0.0000000_wp,  0.0_wp,-0.1450000_wp, 0.0000000_wp, &
      & 0.0_wp,-0.5332978_wp, 0.0000000_wp,  0.0_wp, 1.1522018_wp, 0.0000000_wp],&
      & shape(shell_scale)) + 1.0_wp
   real(wp), parameter :: gexp = 2.0_wp
   real(wp), allocatable :: hubbard(:, :)
   integer :: isp, izp, ish
   type(gamma_coulomb), allocatable :: tmp

   allocate(tmp)
   if (shell) then
      allocate(hubbard(3, mol%nid))
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, shell_count(izp)
            hubbard(ish, isp) = hubbard_parameter(izp) * shell_scale(ish, izp)
         end do
      end do
      call new_gamma_coulomb(tmp, mol, hubbard, shell_count(mol%num))
   else
      hubbard = reshape(hubbard_parameter(mol%num), [1, mol%nid])
      call new_gamma_coulomb(tmp, mol, hubbard)
   end if
   call move_alloc(tmp, coulomb)

end subroutine make_coulomb_g1


!> Factory to create electrostatic objects based on GFN1-xTB values
subroutine make_coulomb_g2(coulomb, mol, shell)

   !> New electrostatic object
   class(coulomb_charge_type), allocatable, intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return a shell resolved object
   logical, intent(in) :: shell

   real(wp), parameter :: hubbard_parameter(20) = [&
      & 0.405771_wp, 0.642029_wp, 0.245006_wp, 0.684789_wp, 0.513556_wp, &
      & 0.538015_wp, 0.461493_wp, 0.451896_wp, 0.531518_wp, 0.850000_wp, &
      & 0.271056_wp, 0.344822_wp, 0.364801_wp, 0.720000_wp, 0.297739_wp, &
      & 0.339971_wp, 0.248514_wp, 0.502376_wp, 0.247602_wp, 0.320378_wp]
   integer, parameter :: shell_count(20) = [&
      & 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3]
   real(wp), parameter :: shell_scale(3, 20) = reshape([&
      & 0.0_wp, 0.0000000_wp, 0.0000000_wp, 0.0_wp, 0.0000000_wp, 0.0000000_wp, &
      & 0.0_wp, 0.1972612_wp, 0.0000000_wp, 0.0_wp, 0.9658467_wp, 0.0000000_wp, &
      & 0.0_wp, 0.3994080_wp, 0.0000000_wp, 0.0_wp, 0.1056358_wp, 0.0000000_wp, &
      & 0.0_wp, 0.1164892_wp, 0.0000000_wp, 0.0_wp, 0.1497020_wp, 0.0000000_wp, &
      & 0.0_wp, 0.1677376_wp, 0.0000000_wp, 0.0_wp, 0.1190576_wp,-0.3200000_wp, &
      & 0.0_wp, 0.1018894_wp, 0.0000000_wp, 0.0_wp, 1.4000000_wp,-0.0500000_wp, &
      & 0.0_wp,-0.0603699_wp, 0.2000000_wp, 0.0_wp,-0.5580042_wp,-0.2300000_wp, &
      & 0.0_wp,-0.1558060_wp,-0.3500000_wp, 0.0_wp,-0.1085866_wp,-0.2500000_wp, &
      & 0.0_wp, 0.4989400_wp, 0.5000000_wp, 0.0_wp,-0.0461133_wp,-0.0100000_wp, &
      & 0.0_wp, 0.3483655_wp, 0.0000000_wp, 0.0_wp, 1.5000000_wp,-0.2500000_wp],&
      & shape(shell_scale)) + 1.0_wp
   real(wp), parameter :: gexp = 2.0_wp
   real(wp), allocatable :: hubbard(:, :)
   integer :: isp, izp, ish
   type(gamma_coulomb), allocatable :: tmp

   allocate(tmp)
   if (shell) then
      allocate(hubbard(3, mol%nid))
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, shell_count(izp)
            hubbard(ish, isp) = hubbard_parameter(izp) * shell_scale(ish, izp)
         end do
      end do
      call new_gamma_coulomb(tmp, mol, hubbard, shell_count(mol%num))
   else
      hubbard = reshape(hubbard_parameter(mol%num), [1, mol%nid])
      call new_gamma_coulomb(tmp, mol, hubbard)
   end if
   call move_alloc(tmp, coulomb)

end subroutine make_coulomb_g2


!> Procedure to create CN based effective charges and gradients from CEH
subroutine get_charges_effceh(wfn, mol, nshell)

   !> New wavefunction object
   class(wavefunction_type), intent(inout) :: wfn

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return shell-resolved charges
   integer, optional, intent(in) :: nshell(:)

   real(wp), parameter :: ceh_cov_radii(20) = 0.5 * [&
   &  2.4040551903_wp,  1.8947380542_wp,  3.4227634078_wp,  3.5225408137_wp, &
   &  3.6150631704_wp,  2.8649682108_wp,  2.4695867541_wp,  2.3533691180_wp, &
   &  2.4992147462_wp,  3.3442607441_wp,  4.4665909451_wp,  4.3877250907_wp, &
   &  4.6647077385_wp,  4.2086223530_wp,  4.4750280107_wp,  4.2847281423_wp, &
   &  3.8560304959_wp,  3.9017061017_wp,  5.2392192639_wp,  5.1872031383_wp]

   real(wp), parameter :: pauling_en_ceh(20) = (1.0_wp/3.98_wp) * [ &
   &  1.9435211923_wp,  3.6116085622_wp,  2.4630915335_wp,  2.0658837656_wp, &
   &  2.3619778807_wp,  2.9484294262_wp,  3.8753937411_wp,  4.6235054741_wp, &
   &  3.9800000000_wp,  3.5124865506_wp,  2.3578254072_wp,  2.4225832022_wp, &
   &  2.1120078826_wp,  2.4607564741_wp,  2.7410779326_wp,  3.3517034720_wp, &
   &  4.1093492601_wp,  3.7979559518_wp,  2.4147937668_wp,  2.1974781961_wp]
   
   real(wp), allocatable :: lattr(:, :), cn_en(:), dcn_endr(:, :, :), dcn_endL(:, :, :)
   type(erf_en_ncoord_type) :: ncoord_en
   integer :: iat, ii, ish

   allocate(cn_en(mol%nat), dcn_endr(3, mol%nat, mol%nat), dcn_endL(3, 3, mol%nat))

   ! Get electronegativity-weighted coordination number
   call new_erf_en_ncoord(ncoord_en, mol, &
   & rcov=ceh_cov_radii(mol%num), en=pauling_en_ceh(mol%num))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call get_coordination_number(ncoord_en, mol, lattr, cutoff, &
      & cn_en, dcn_endr, dcn_endL)

   ! Get effective charges and their gradients
   call get_effective_qat(mol, cn_en, wfn%qat, &
      & dcn_endr, dcn_endL, wfn%dqatdr, wfn%dqatdL)

   ! If shell-resolved charges are requested, partition atomic charges on shells
   if(present(nshell)) then
      ii = 0
      do iat = 1, mol%nat
         do ish = 1, nshell(iat)
            wfn%qsh(ii+ish, :) = wfn%qat(iat, :) / real(nshell(iat), wp)
            wfn%dqshdr(:, :, ii+ish, :) = wfn%dqatdr(:, :, iat, :) / real(nshell(iat), wp)
            wfn%dqshdL(:, :, ii+ish, :) = wfn%dqatdL(:, :, iat, :) / real(nshell(iat), wp)
         end do 
         ii = ii + nshell(iat)
      end do 
   end if

end subroutine get_charges_effceh


subroutine test_generic(error, mol, qat, qsh, make_coulomb, ref, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges for this structure
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges for this structure
   real(wp), intent(in), optional :: qsh(:)

   !> Factory to create new electrostatic objects
   procedure(coulomb_maker) :: make_coulomb

   !> Reference value to check against
   real(wp), intent(in) :: ref

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   class(coulomb_charge_type), allocatable :: coulomb
   type(container_cache) :: cache
   real(wp) :: energy(mol%nat), sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), numgrad(:, :)
   real(wp) :: thr_
   type(wavefunction_type) :: wfn

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   if (present(qsh)) then
      wfn%qsh = reshape(qsh, [size(qsh), 1])
   else
      wfn%qsh = reshape(qat, [size(qat), 1])
   end if
   wfn%qat = reshape(qat, [size(qat), 1])
   call make_coulomb(coulomb, mol, present(qsh))
   call coulomb%update(mol, cache)
   call coulomb%get_energy(mol, cache, wfn, energy)

   call check(error, sum(energy), ref, thr=thr_)
   if (allocated(error)) then
      print*,ref, sum(energy)
   end if

end subroutine test_generic


subroutine test_numgrad(error, mol, qat, qsh, make_coulomb)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges for this structure
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges for this structure
   real(wp), intent(in), optional :: qsh(:)

   !> Factory to create new electrostatic objects
   procedure(coulomb_maker) :: make_coulomb

   integer :: iat, ic
   class(coulomb_charge_type), allocatable :: coulomb
   type(container_cache) :: cache
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), numgrad(:, :)
   real(wp), parameter :: step = 5.0e-5_wp
   type(wavefunction_type) :: wfn

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   if (present(qsh)) then
      wfn%qsh = reshape(qsh, [size(qsh), 1])
   else
      wfn%qsh = reshape(qat, [size(qat), 1])
   end if
   wfn%qat = reshape(qat, [size(qat), 1])
   call make_coulomb(coulomb, mol, present(qsh))

   do iat = 1, mol%nat
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call coulomb%update(mol, cache)
         call coulomb%get_energy(mol, cache, wfn, er)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call coulomb%update(mol, cache)
         call coulomb%get_energy(mol, cache, wfn, el)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   call coulomb%update(mol, cache)
   call coulomb%get_gradient(mol, cache, wfn, gradient, sigma)

   if (any(abs(gradient - numgrad) > thr2)) then
      call test_failed(error, "Gradient of energy does not match")
      print'(3es21.14)', gradient-numgrad
   end if

end subroutine test_numgrad


subroutine test_numpotgrad(error, mol, get_charges, make_coulomb, shell)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Procedure to provide atomic charges (and their gradients)
   procedure(charge_maker) :: get_charges

   !> Factory to create new electrostatic objects
   procedure(coulomb_maker) :: make_coulomb

   !> Test shell-resolved
   logical, intent(in) :: shell

   integer :: iat, ic
   type(basis_type) :: bas
   type(potential_type) :: potl, potr
   class(coulomb_charge_type), allocatable :: coulomb
   type(container_cache) :: cache
   real(wp), allocatable :: numpotgrad(:, :, :)
   real(wp), parameter :: step = 5.0e-5_wp
   type(wavefunction_type) :: wfn

   ! Setup potentials and wavefunction with dummy basis set 
   call make_basis(bas, mol, 6)
   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, 1, 0.0_wp, .true.)
   call new_potential(potr, mol, bas, 1, .true.)
   call new_potential(potl, mol, bas, 1, .true.)
   call make_coulomb(coulomb, mol, shell)

   if(shell) then
      allocate(numpotgrad(3, mol%nat, bas%nsh), source=0.0_wp)
   else
      allocate(numpotgrad(3, mol%nat, mol%nat), source=0.0_wp)
   end if

   do iat = 1, mol%nat
      do ic = 1, 3
         call potr%reset
         call potl%reset
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_charges(wfn, mol, coulomb%nshell)
         call coulomb%update(mol, cache)
         call coulomb%get_potential(mol, cache, wfn, potr)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call get_charges(wfn, mol, coulomb%nshell)
         call coulomb%update(mol, cache)
         call coulomb%get_potential(mol, cache, wfn, potl)
         
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         if(shell) then
            numpotgrad(ic, iat, :) = 0.5_wp*(potr%vsh(:,1) - potl%vsh(:,1))/step
         else
            numpotgrad(ic, iat, :) = 0.5_wp*(potr%vat(:,1) - potl%vat(:,1))/step
         end if
      end do
   end do

   call get_charges(wfn, mol, coulomb%nshell)
   call coulomb%update(mol, cache)
   call coulomb%get_potential(mol, cache, wfn, potl)
   call coulomb%get_potential_gradient(mol, cache, wfn, potl)
   
   if(shell) then
      if (any(abs(potl%dvshdr(:,:,:,1) - numpotgrad) > thr2)) then
         call test_failed(error, "Gradient of shell-resolved potential does not match")
         print'(3es21.14)', potl%dvshdr(:,:,:,1) - numpotgrad
      end if
   else
      if (any(abs(potl%dvatdr(:,:,:,1) - numpotgrad) > thr2)) then
         call test_failed(error, "Gradient of atom-resolved potential does not match")
         print'(3es21.14)', potl%dvatdr(:,:,:,1) - numpotgrad
      end if
   end if

end subroutine test_numpotgrad


subroutine test_numsigma(error, mol, qat, qsh, make_coulomb)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges for this structure
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges for this structure
   real(wp), intent(in), optional :: qsh(:)

   !> Factory to create new electrostatic objects
   procedure(coulomb_maker) :: make_coulomb

   integer :: ic, jc
   class(coulomb_charge_type), allocatable :: coulomb
   type(container_cache) :: cache
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3), eps(3, 3), numsigma(3, 3)
   real(wp), allocatable :: gradient(:, :), xyz(:, :), lattice(:, :)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 1.0e-5_wp
   type(wavefunction_type) :: wfn

   allocate(gradient(3, mol%nat), xyz(3, mol%nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   if (present(qsh)) then
      wfn%qsh = reshape(qsh, [size(qsh), 1])
   else
      wfn%qsh = reshape(qat, [size(qat), 1])
   end if
   wfn%qat = reshape(qat, [size(qat), 1])
   call make_coulomb(coulomb, mol, present(qsh))

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   if (any(mol%periodic)) lattice = mol%lattice
   do ic = 1, 3
      do jc = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call coulomb%update(mol, cache)
         call coulomb%get_energy(mol, cache, wfn, er)
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call coulomb%update(mol, cache)
         call coulomb%get_energy(mol, cache, wfn, el)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         if (allocated(lattice)) mol%lattice = lattice
         numsigma(jc, ic) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   call coulomb%update(mol, cache)
   call coulomb%get_gradient(mol, cache, wfn, gradient, sigma)

   if (any(abs(sigma - numsigma) > thr2)) then
      call test_failed(error, "Strain derivatives do not match")
      print'(3es21.14)', sigma-numsigma
   end if

end subroutine test_numsigma


subroutine test_numpotsigma(error, mol, get_charges, make_coulomb, shell)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Procedure to provide atomic charges (and their gradients)
   procedure(charge_maker) :: get_charges

   !> Factory to create new electrostatic objects
   procedure(coulomb_maker) :: make_coulomb

   !> Test shell-resolved
   logical, intent(in) :: shell

   integer :: ic, jc
   type(basis_type) :: bas
   type(potential_type) :: potl, potr
   class(coulomb_charge_type), allocatable :: coulomb
   type(container_cache) :: cache
   real(wp) :: eps(3, 3)
   real(wp), allocatable :: numpotsigma(:, :, :), xyz(:, :), lattice(:, :)
   real(wp), parameter :: unity(3, 3) = reshape(&
   & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 5.0e-5_wp
   type(wavefunction_type) :: wfn

   allocate(xyz(3, mol%nat), source=0.0_wp)

   ! Setup potentials and wavefunction with dummy basis set 
   call make_basis(bas, mol, 6)
   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, 1, 0.0_wp, .true.)
   call new_potential(potr, mol, bas, 1, .true.)
   call new_potential(potl, mol, bas, 1, .true.)
   call make_coulomb(coulomb, mol, shell)

   if(shell) then
      allocate(numpotsigma(3, 3, bas%nsh), source=0.0_wp)
   else
      allocate(numpotsigma(3, 3, mol%nat), source=0.0_wp)
   end if

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   if (any(mol%periodic)) lattice = mol%lattice
   do ic = 1, 3
      do jc = 1, 3
         call potr%reset
         call potl%reset
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call get_charges(wfn, mol, coulomb%nshell)
         call coulomb%update(mol, cache)
         call coulomb%get_potential(mol, cache, wfn, potr)

         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call get_charges(wfn, mol, coulomb%nshell)
         call coulomb%update(mol, cache)
         call coulomb%get_potential(mol, cache, wfn, potl)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         if (allocated(lattice)) mol%lattice = lattice
         if(shell) then
            numpotsigma(jc, ic, :) = 0.5_wp*(potr%vsh(:,1) - potl%vsh(:,1))/step
         else
            numpotsigma(jc, ic, :) = 0.5_wp*(potr%vat(:,1) - potl%vat(:,1))/step
         end if
      end do
   end do

   call get_charges(wfn, mol, coulomb%nshell)
   call coulomb%update(mol, cache)
   call coulomb%get_potential(mol, cache, wfn, potl)
   call coulomb%get_potential_gradient(mol, cache, wfn, potl)

   if(shell) then
      if (any(abs(potl%dvshdL(:,:,:,1) - numpotsigma) > thr2)) then
         call test_failed(error, "Sigma of shell-resolved potential does not match")
         print'(3es21.14)', potl%dvshdL(:,:,:,1) - numpotsigma
      end if
   else
      if (any(abs(potl%dvatdL(:,:,:,1) - numpotsigma) > thr2)) then
         call test_failed(error, "Sigma of atom-resolved potential does not match")
         print'(3es21.14)', potl%dvatdL(:,:,:,1) - numpotsigma
      end if
   end if

end subroutine test_numpotsigma


subroutine test_e_effective_m01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 7.73347900345264E-1_wp, 1.07626888948184E-1_wp,-3.66999593831010E-1_wp,&
      & 4.92833325937897E-2_wp,-1.83332156197733E-1_wp, 2.33302086605469E-1_wp,&
      & 6.61837152062315E-2_wp,-5.43944165050002E-1_wp,-2.70264356583716E-1_wp,&
      & 2.66618968841682E-1_wp, 2.62725033202480E-1_wp,-7.15315510172571E-2_wp,&
      &-3.73300777019193E-1_wp, 3.84585237785621E-2_wp,-5.05851088366940E-1_wp,&
      & 5.17677238544189E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "01")
   call test_generic(error, mol, qat, qsh, make_coulomb_e1, 0.10952019883948200_wp)

end subroutine test_e_effective_m01


subroutine test_e_effective_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 7.38394711236234E-2_wp,-1.68354976558608E-1_wp,-3.47642833746823E-1_wp,&
      &-7.05489267186003E-1_wp, 7.73548301641266E-1_wp, 2.30207581365386E-1_wp,&
      & 1.02748501676354E-1_wp, 9.47818107467040E-2_wp, 2.44260351729187E-2_wp,&
      & 2.34984927037408E-1_wp,-3.17839896393030E-1_wp, 6.67112994818879E-1_wp,&
      &-4.78119977010488E-1_wp, 6.57536027459275E-2_wp, 1.08259054549882E-1_wp,&
      &-3.58215329983396E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "02")
   call test_generic(error, mol, qat, qsh, make_coulomb_e2, 0.10635843572138280_wp)

end subroutine test_e_effective_m02


subroutine test_e_effective_m10(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-1.68014634286730E-2_wp, 4.42632083957210E-1_wp, 5.27534110370006E-3_wp, &
      &-3.49386920985384E-1_wp,-2.13178401684440E-1_wp, 4.07942205949245E-1_wp, &
      &-4.08514260972236E-1_wp, 1.34978625814380E-1_wp,-2.48254330281858E-1_wp, &
      &-3.65112235756872E-1_wp, 3.19858617682441E-1_wp, 2.96731604233838E-2_wp, &
      & 2.84061022228221E-1_wp,-3.25028474749853E-2_wp,-5.33914616408729E-3_wp, &
      & 1.46685494419848E-2_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "10")
   call test_generic(error, mol, qat, qsh, make_coulomb_g1, 7.964758847421499E-2_wp, thr2)

end subroutine test_e_effective_m10


subroutine test_e_effective_oxacb(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 3.41731844312030E-1_wp, 3.41716020106239E-1_wp, 3.41730526585671E-1_wp,&
      & 3.41714427217954E-1_wp, 3.80996046757999E-1_wp, 3.80989821246195E-1_wp,&
      & 3.81000747720282E-1_wp, 3.80990494183703E-1_wp,-3.70406587264474E-1_wp,&
      &-3.70407565207006E-1_wp,-3.70417590212352E-1_wp,-3.70399716470705E-1_wp,&
      &-3.52322260586075E-1_wp,-3.52304269439196E-1_wp,-3.52313440903261E-1_wp,&
      &-3.52298498047004E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "oxacb")
   call test_generic(error, mol, qat, qsh, make_coulomb_e2, 0.10130450083781417_wp)

end subroutine test_e_effective_oxacb


subroutine test_e_effective_oxacb_sc(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat1(*) = [&
      & 3.41731844312030E-1_wp, 3.41716020106239E-1_wp, 3.41730526585671E-1_wp,&
      & 3.41714427217954E-1_wp, 3.80996046757999E-1_wp, 3.80989821246195E-1_wp,&
      & 3.81000747720282E-1_wp, 3.80990494183703E-1_wp,-3.70406587264474E-1_wp,&
      &-3.70407565207006E-1_wp,-3.70417590212352E-1_wp,-3.70399716470705E-1_wp,&
      &-3.52322260586075E-1_wp,-3.52304269439196E-1_wp,-3.52313440903261E-1_wp,&
      &-3.52298498047004E-1_wp]
   integer, parameter :: supercell(*) = [2, 2, 2]
   real(wp), parameter :: qat(*) = [spread(qat1, 2, product(supercell))]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "oxacb")
   call make_supercell(mol, supercell)
   call test_generic(error, mol, qat, qsh, make_coulomb_e2, &
      & 0.10130450083781417_wp*product(supercell), 1.0e-7_wp)

end subroutine test_e_effective_oxacb_sc


subroutine make_supercell(mol, rep)
   type(structure_type), intent(inout) :: mol
   integer, intent(in) :: rep(3)

   real(wp), allocatable :: xyz(:, :), lattice(:, :)
   integer, allocatable :: num(:)
   integer :: i, j, k, c

   num = reshape(spread([mol%num(mol%id)], 2, product(rep)), [product(rep)*mol%nat])
   lattice = reshape(&
      [rep(1)*mol%lattice(:, 1), rep(2)*mol%lattice(:, 2), rep(3)*mol%lattice(:, 3)], &
      shape(mol%lattice))
   allocate(xyz(3, product(rep)*mol%nat))
   c = 0
   do i = 0, rep(1)-1
      do j = 0, rep(2)-1
         do k = 0, rep(3)-1
            xyz(:, c+1:c+mol%nat) = mol%xyz &
               & + spread(matmul(mol%lattice, [real(wp):: i, j, k]), 2, mol%nat)
            c = c + mol%nat
         end do
      end do
   end do

   call new(mol, num, xyz, lattice=lattice)
end subroutine make_supercell


subroutine test_g_effective_m03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-1.77788256288236E-1_wp,-8.22943267808161E-1_wp, 4.04578389873281E-2_wp,&
      & 5.79710531992282E-1_wp, 6.99601887637659E-1_wp, 6.84309612639107E-2_wp,&
      &-3.42971414989811E-1_wp, 4.64954031865410E-2_wp, 6.77012204116428E-2_wp,&
      & 8.49931225363225E-2_wp,-5.22285304699699E-1_wp,-2.92515001764712E-1_wp,&
      &-3.98375452377043E-1_wp, 2.09769668297792E-1_wp, 7.23140464830357E-1_wp,&
      & 3.65775987838250E-2_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "03")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_e1)

end subroutine test_g_effective_m03


subroutine test_g_effective_m04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 9.33596160193497E-2_wp,-3.41088061922851E-1_wp, 7.32474961830646E-2_wp,&
      &-2.21649975471802E-1_wp, 6.24413528413759E-3_wp, 1.07366683260668E-1_wp,&
      & 1.25982547197317E-1_wp, 9.65935501843890E-2_wp, 1.02704543049803E-1_wp,&
      & 1.45380937882263E-1_wp,-1.55978251071729E-1_wp, 3.42948437914661E-1_wp,&
      & 5.65504846503244E-2_wp,-3.37789986050220E-1_wp, 1.13510089629769E-1_wp,&
      &-2.07382246739143E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "04")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_e2)

end subroutine test_g_effective_m04


subroutine test_g_effective_m11(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-2.90839692052395E-1_wp,-2.66254004019676E-1_wp,-1.14585349859674E-1_wp, &
      &-1.94168606876184E-1_wp, 1.69243400195097E-1_wp, 5.94099012700995E-2_wp, &
      & 3.58048105537982E-1_wp, 3.65662316444953E-2_wp, 3.38204991437465E-1_wp, &
      &-4.07570397699211E-1_wp, 5.27525279437458E-1_wp,-2.17656283937311E-1_wp, &
      &-3.01540791618602E-1_wp, 2.89587569744254E-1_wp,-1.87215698022911E-2_wp, &
      & 3.27512165984881E-2_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "11")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_g1)

end subroutine test_g_effective_m11


subroutine test_g_effective_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 4.56275672862067E-1_wp, 4.56284770386671E-1_wp, 4.56284770386671E-1_wp,&
      & 4.56284770386671E-1_wp,-2.28127680925611E-1_wp,-2.28138283131909E-1_wp,&
      &-2.28145770512561E-1_wp,-2.28145770512561E-1_wp,-2.28150142163058E-1_wp,&
      &-2.28145770512561E-1_wp,-2.28138283131909E-1_wp,-2.28138283131909E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "CO2")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_e2)

end subroutine test_g_effective_co2


subroutine test_g_effective_urea(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 5.55723890858218E-1_wp, 5.55765354442035E-1_wp, 2.50200231242017E-1_wp,&
      & 2.50282053284422E-1_wp, 2.39786980460652E-1_wp, 2.39895142481200E-1_wp,&
      & 2.50103678240412E-1_wp, 2.50425041601730E-1_wp, 2.39464477136495E-1_wp,&
      & 2.40360053062669E-1_wp,-4.38369096728919E-1_wp,-4.38451412936599E-1_wp,&
      &-4.38310020776279E-1_wp,-4.38617373848238E-1_wp,-6.59141030224988E-1_wp,&
      &-6.59117968294813E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "urea")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_g1)

end subroutine test_g_effective_urea


subroutine test_s_effective_m05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-2.01138111283277E-1_wp, 1.30358706339300E-1_wp, 9.38825924720944E-2_wp,&
      & 8.92795900801844E-2_wp, 5.13625440660610E-2_wp,-2.65500121876709E-2_wp,&
      & 9.26496972837658E-2_wp,-9.61095258223972E-2_wp,-4.92845009674246E-1_wp,&
      & 2.66730531684206E-1_wp, 3.37256104303071E-2_wp, 1.63170419985976E-1_wp,&
      & 6.91343155032824E-2_wp, 1.04287482572171E-1_wp, 6.09307909835941E-2_wp,&
      &-3.38869622433350E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "05")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_e1)

end subroutine test_s_effective_m05


subroutine test_s_effective_m06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-2.13983049532933E-1_wp,-5.10521279217923E-1_wp, 7.70190120699491E-2_wp,&
      &-3.68835155548212E-1_wp,-4.08747874260092E-1_wp,-4.09471309598929E-2_wp,&
      & 2.94164204769172E-1_wp, 9.76819709672870E-2_wp,-7.84337476935767E-3_wp,&
      & 7.07702520795024E-1_wp, 2.38774840136381E-1_wp, 1.08934666297455E-1_wp,&
      & 1.10156911889136E-1_wp, 9.25098455002779E-2_wp,-1.96776817442259E-1_wp,&
      & 2.07107093059868E-2_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "06")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_e2)

end subroutine test_s_effective_m06


subroutine test_s_effective_m12(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 2.57523069224332E-1_wp,-7.13834599743419E-2_wp,-9.43945788149514E-2_wp, &
      & 1.02376579062895E-2_wp, 1.97960723912756E-1_wp, 3.16253846282241E-1_wp, &
      &-3.91548233613895E-1_wp,-1.68829398890385E-1_wp,-3.99798824173873E-1_wp, &
      & 4.22333212859043E-1_wp, 3.16980455307683E-1_wp,-1.09800808615267E-1_wp, &
      &-1.07582049146789E-1_wp,-4.48115660454027E-1_wp, 7.19397904194672E-2_wp, &
      & 1.98224257774052E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "12")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_g1)

end subroutine test_s_effective_m12


subroutine test_s_effective_ammonia(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 2.95376975876519E-1_wp, 2.95376975876519E-1_wp, 2.95376975876519E-1_wp,&
      & 2.95329109335847E-1_wp, 2.95332441468412E-1_wp, 2.95347202855778E-1_wp,&
      & 2.95347202855779E-1_wp, 2.95329109335848E-1_wp, 2.95332441468411E-1_wp,&
      & 2.95347202855777E-1_wp, 2.95329109335847E-1_wp, 2.95332441468412E-1_wp,&
      &-8.86118742099358E-1_wp,-8.86012815503436E-1_wp,-8.86012815503437E-1_wp,&
      &-8.86012815503434E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "ammonia")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_e2)

end subroutine test_s_effective_ammonia


subroutine test_s_effective_pyrazine(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 1.23705026512437E-1_wp, 1.22537765989959E-1_wp, 1.23932626831231E-1_wp,&
      & 1.22418959326822E-1_wp, 1.23788033684569E-1_wp, 1.23643389058068E-1_wp,&
      & 1.22389811551880E-1_wp, 1.22402837718155E-1_wp, 3.17174594622166E-2_wp,&
      & 3.15585817789390E-2_wp, 3.14290005061543E-2_wp, 3.10358506314526E-2_wp,&
      & 3.11084225356749E-2_wp, 3.11187325528499E-2_wp, 3.10707965296438E-2_wp,&
      & 3.13508824430484E-2_wp,-3.09107070033859E-1_wp,-3.09144123845085E-1_wp,&
      &-3.08473365847409E-1_wp,-3.08483617386745E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "pyrazine")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_g2)

end subroutine test_s_effective_pyrazine

subroutine write_charges(mol)
   use dftd4_charge, only : get_charges
   type(structure_type) :: mol
   real(wp), allocatable :: qat(:)
   allocate(qat(mol%nat))
   call get_charges(mol, qat)
   write(*, '(3x,a)') "real(wp), parameter :: qat(*) = [&"
   write(*, '(*(6x,"&",3(es20.14e1, "_wp":, ","),"&", /))', advance='no') qat
   write(*, '(a)') "]"
end subroutine

subroutine test_e_effective_m07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-1.49712318034775E-1_wp, 2.12665850975202E-1_wp, 3.35977061494489E-1_wp, &
      & 3.16737890491354E-2_wp, 4.12434432866631E-2_wp,-3.21014009885608E-1_wp, &
      &-3.06535419089495E-1_wp,-5.36251066565321E-1_wp, 4.48758364798896E-1_wp, &
      & 6.00309584480896E-2_wp,-2.75470557482709E-1_wp, 3.60263594022551E-1_wp, &
      & 3.77425314022796E-1_wp,-6.30561365518420E-1_wp,-2.50675164323255E-1_wp, &
      & 6.02181524801775E-1_wp]
   real(wp), parameter :: qsh(*) = [&
      & 8.85960229060055E-1_wp,-1.03567241653662E+0_wp, 2.34499192077770E-1_wp, &
      &-2.18333480864186E-2_wp, 1.09026104661485E+0_wp,-7.54283954798938E-1_wp, &
      & 4.12740327203921E-2_wp,-9.60021563849638E-3_wp, 5.17672944681095E-2_wp, &
      &-1.05238375989861E-2_wp, 5.94332546515988E-2_wp,-3.94897989828280E-1_wp, &
      & 1.44506731071946E-2_wp, 1.57870128213110E-1_wp,-4.64405557396352E-1_wp, &
      & 4.78122334280047E-1_wp,-1.01437364107707E+0_wp, 9.10337331767967E-1_wp, &
      &-4.61579000227231E-1_wp, 9.07619848805192E-2_wp,-3.07310018122722E-2_wp, &
      & 1.13955875471381E-1_wp,-3.99913576087036E-1_wp, 1.04872002787662E-2_wp, &
      & 4.12951024314537E-1_wp,-5.26874026571100E-2_wp, 4.04435991881125E-1_wp, &
      &-2.70107073714884E-2_wp, 3.13675308978710E-1_wp,-9.44236655190031E-1_wp, &
      & 1.75329569882602E-1_wp,-4.26004749886597E-1_wp, 1.24860566181157E+0_wp, &
      &-6.46424080267374E-1_wp]
   real(wp), allocatable :: qsh0(:)

   call get_structure(mol, "MB16-43", "07")
   call test_generic(error, mol, qat, qsh0, make_coulomb_e1, 0.13650692645610521_wp)
   call test_generic(error, mol, qat, qsh, make_coulomb_e1, 0.12017418620257683_wp)

end subroutine test_e_effective_m07

subroutine test_g_effective_m08(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-2.11048312695985E-1_wp,-5.02011645803230E-1_wp, 4.15238062649689E-1_wp, &
      &-3.25959600753673E-1_wp, 2.51473641195433E-2_wp, 2.93748490123740E-1_wp, &
      & 2.56736194030896E-2_wp, 2.38762690307426E-2_wp,-6.03118603733083E-1_wp, &
      & 3.91990240426822E-1_wp, 8.97114734113785E-1_wp, 1.93532936362436E-1_wp, &
      &-1.03136268223866E-1_wp,-1.04447608767710E-1_wp,-2.64818891498402E-2_wp, &
      &-3.90117787102468E-1_wp]
   real(wp), parameter :: qsh(*) = [&
      & 9.06259904944829E-1_wp,-1.11730821567902E+0_wp, 2.78017329305492E-1_wp, &
      &-7.80028989546297E-1_wp, 1.11352815063389E+0_wp,-6.98290073981154E-1_wp, &
      & 2.03943236255318E-1_wp,-5.29902840441233E-1_wp, 4.38219939650397E-2_wp, &
      &-1.86746328945826E-2_wp, 4.65996457236599E-2_wp, 4.97590807484258E-1_wp, &
      &-2.50441962186972E-1_wp, 4.83295451755440E-2_wp,-2.26559244782012E-2_wp, &
      & 4.50331992744248E-2_wp,-2.11569328297532E-2_wp, 3.12470620007346E-1_wp, &
      &-9.15589243372491E-1_wp, 1.06394261835743E+0_wp,-6.71952361588756E-1_wp, &
      & 1.82322476598938E+0_wp,-9.26110009158329E-1_wp, 9.78357111140355E-1_wp, &
      &-7.84824170464332E-1_wp,-9.43549308434806E-2_wp,-8.78133979988158E-3_wp, &
      &-7.07783143624696E-2_wp,-3.36692984665466E-2_wp, 6.75375129657761E-1_wp, &
      &-7.01857024422455E-1_wp, 2.11598132242645E-1_wp,-6.01715925641418E-1_wp]

   call get_structure(mol, "MB16-43", "08")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_e1)

end subroutine test_g_effective_m08

subroutine test_s_effective_m09(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-1.13260038539900E-2_wp, 1.10070523471231E-2_wp,-9.97165474630829E-2_wp, &
      &-8.78527301724521E-2_wp, 2.89049242695863E-1_wp, 3.57284006856323E-2_wp, &
      &-1.73226219187217E-1_wp, 1.61174372420268E-1_wp,-8.89089419183055E-2_wp, &
      & 3.23950178196666E-2_wp, 1.88420688366637E-1_wp, 4.14882523279327E-2_wp, &
      &-2.23498403532295E-1_wp,-3.55334728213004E-1_wp,-7.15753987897201E-2_wp, &
      & 3.52175946466941E-1_wp]
   real(wp), parameter :: qsh(*) = [&
      & 1.40887956776581E-3_wp,-1.27348820058716E-2_wp, 2.63961183739554E-2_wp, &
      &-1.53890176131402E-2_wp,-8.73648546608390E-2_wp,-1.23517435478204E-2_wp, &
      &-8.71021014527735E-2_wp,-7.50559382736492E-4_wp, 7.82044211296174E-1_wp, &
      &-4.92995083533018E-1_wp, 4.84143136555792E-2_wp,-1.26858387490357E-2_wp, &
      & 9.72488073646510E-1_wp,-1.14571448042039E+0_wp, 1.07574874045191E+0_wp, &
      &-9.14574293473561E-1_wp,-7.63358458276189E-2_wp,-1.25730981035572E-2_wp, &
      & 4.44349073468088E-2_wp,-1.20397879426510E-2_wp, 5.20245311277456E-1_wp, &
      & 1.92282483805197E-1_wp,-5.24107355799204E-1_wp, 5.39382871928999E-2_wp, &
      &-1.24499232808976E-2_wp, 7.97368410133983E-2_wp,-3.13209082796440E-1_wp, &
      & 9.97387287057362E-3_wp, 1.94446888375020E-1_wp,-5.49781435696375E-1_wp, &
      &-6.89789344411558E-2_wp,-2.59643153694089E-3_wp, 1.09519797601190E+0_wp, &
      &-7.43022154621128E-1_wp]

   call get_structure(mol, "MB16-43", "09")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_e1)

end subroutine test_s_effective_m09

subroutine test_e_effective_m13(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-2.62608233282119E-1_wp, 3.73633121487967E-1_wp, 1.51424532948944E-1_wp, &
      & 8.11274419840145E-2_wp, 4.55582555217907E-1_wp, 1.89469664895825E-1_wp, &
      &-3.59350817183894E-1_wp,-1.38911850317377E-1_wp,-1.83689392824396E-1_wp, &
      &-1.88906495161279E-1_wp, 5.33440028285669E-2_wp, 1.94112134916556E-1_wp, &
      & 2.02080948377078E-1_wp, 1.74595453525400E-1_wp,-4.46124496927388E-1_wp, &
      &-2.95778570663624E-1_wp]
   real(wp), parameter :: qsh(*) = [&
      & 5.72134559421376E-2_wp,-2.68193499948548E-1_wp,-5.17064903935052E-2_wp, &
      & 3.73632853173886E-1_wp, 1.51424477665324E-1_wp, 8.11205008366953E-2_wp, &
      & 1.05453876337982E+0_wp,-4.64589774617786E-1_wp,-1.34371775944868E-1_wp, &
      & 3.21958772020979E-1_wp,-1.32435004307411E-1_wp, 2.88638899705825E-1_wp, &
      &-6.47972813769995E-1_wp, 6.82118177705109E-2_wp,-1.10631364729565E-1_wp, &
      &-9.64955180671905E-2_wp, 1.27185941911165E-1_wp,-3.10873201558534E-1_wp, &
      & 9.97036415531523E-2_wp,-2.88615729133477E-1_wp,-1.09656595674679E-1_wp, &
      & 1.63000176490660E-1_wp, 1.94112048312228E-1_wp,-4.48012133376332E-2_wp, &
      & 2.46906120671464E-1_wp, 1.74594629792853E-1_wp, 2.06932598673206E-1_wp, &
      &-6.52990922441455E-1_wp,-6.98603488893812E-3_wp,-2.88854759086314E-1_wp]
   real(wp), allocatable :: qsh0(:)

   call get_structure(mol, "MB16-43", "13")
   call test_generic(error, mol, qat, qsh0, make_coulomb_g2, 4.4263535114062461E-2_wp, thr2)
   if (allocated(error)) return
   call test_generic(error, mol, qat, qsh, make_coulomb_g2, 6.1195162961497636E-2_wp, thr2)

end subroutine test_e_effective_m13

subroutine test_g_effective_m14(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-5.35371225694038E-1_wp, 1.27905155882876E-1_wp, 2.06910619292535E-1_wp, &
      &-1.93061647443670E-1_wp, 5.46833043573218E-1_wp, 2.98577669101319E-1_wp, &
      &-3.62405585534705E-1_wp, 2.07231134137244E-1_wp, 2.85826164709174E-1_wp, &
      &-1.76518940177473E-1_wp, 9.44972704818130E-2_wp,-1.17451405142691E-1_wp, &
      &-1.41198286268662E-1_wp, 1.05227974737201E-2_wp,-1.31666840327078E-1_wp, &
      &-1.20629924063582E-1_wp]
   real(wp), parameter :: qsh(*) = [&
      & 2.92177048596496E-1_wp,-8.27551283559270E-1_wp, 3.81811514612779E-1_wp, &
      & 4.17784916263666E-1_wp,-6.71683789364610E-1_wp,-5.11576334445887E-2_wp, &
      & 2.68488368058409E-1_wp,-1.04391705441501E-2_wp,-1.93062932974169E-1_wp, &
      & 7.61952673849415E-1_wp,-2.15114295745990E-1_wp, 2.98579359296096E-1_wp, &
      & 2.43594536377056E-2_wp,-3.79828052486015E-1_wp,-6.93861559657517E-3_wp, &
      & 9.40499498626918E-1_wp, 3.55733525506819E-1_wp,-1.08899859046528E+0_wp, &
      & 2.85827068858603E-1_wp,-1.76521898491791E-1_wp, 1.16110895118352E+0_wp, &
      &-1.06660369382576E+0_wp,-1.17453146019547E-1_wp,-1.41199843368407E-1_wp, &
      & 1.11757931464608E+0_wp,-1.10704872438469E+0_wp,-1.31668946428877E-1_wp, &
      &-1.20631076436794E-1_wp]

   call get_structure(mol, "MB16-43", "14")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_g2)

end subroutine test_g_effective_m14

subroutine test_s_effective_m15(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-5.41402496268596E-2_wp,-1.33777153976276E-1_wp, 4.14313829600631E-1_wp, &
      &-1.16641170075389E-1_wp, 4.56021377424607E-1_wp,-5.20378766868989E-1_wp, &
      &-1.63965423099635E-1_wp,-7.65345311273482E-2_wp, 2.08304494730413E-1_wp, &
      &-1.71827679329874E-1_wp,-3.30458156481514E-1_wp, 5.58638267294323E-1_wp, &
      &-3.10094747569162E-1_wp, 1.56794592474036E-1_wp, 2.13459748796815E-1_wp, &
      &-1.29714432165776E-1_wp]
   real(wp), parameter :: qsh(*) = [&
      &-5.41503330947710E-2_wp,-1.33779040525986E-1_wp, 6.99211871952128E-2_wp, &
      & 5.19410210372243E-1_wp,-1.75021497979320E-1_wp,-1.16646288193380E-1_wp, &
      & 5.88579806566877E-1_wp,-1.32544517759254E-1_wp, 7.82133136453700E-4_wp, &
      &-5.23133226533954E-1_wp, 1.95314227063571E-3_wp,-1.63971802434084E-1_wp, &
      &-7.65416499768423E-2_wp, 6.26092659280207E-1_wp, 3.00350609071998E-1_wp, &
      &-7.18137148569625E-1_wp,-1.71830433005059E-1_wp,-1.31761941444373E-1_wp, &
      &-1.98713319565700E-1_wp, 8.86264348375974E-2_wp, 7.46929250011706E-1_wp, &
      &-2.76884353276538E-1_wp, 4.27025703238206E-2_wp,-3.52590124894769E-1_wp, &
      &-2.13102162478342E-4_wp, 1.13142747674328E+0_wp,-9.74609683930292E-1_wp, &
      & 1.02802427493832E+0_wp,-8.14556120567196E-1_wp,-1.29715170834719E-1_wp]

   call get_structure(mol, "MB16-43", "15")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_g2)

end subroutine test_s_effective_m15

subroutine test_ceh_potgrad_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_eceh, .false.)

end subroutine test_ceh_potgrad_lih

subroutine test_ceh_potgrad_mb15(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "15")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_eceh, .false.)

end subroutine test_ceh_potgrad_mb15

subroutine test_ceh_potgrad_mb16(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "16")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_eceh, .true.)

end subroutine test_ceh_potgrad_mb16

subroutine test_ceh_potsigma_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_eceh, .false.)

end subroutine test_ceh_potsigma_lih

subroutine test_ceh_potsigma_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "CO2")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_eceh, .false.)

end subroutine test_ceh_potsigma_co2

subroutine test_ceh_potsigma_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_eceh, .false.)

end subroutine test_ceh_potsigma_mb05

subroutine test_ceh_potsigma_mb17(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "17")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_eceh, .true.)

end subroutine test_ceh_potsigma_mb17

end module test_coulomb_charge
