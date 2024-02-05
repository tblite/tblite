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
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_cutoff, only : get_lattice_points
   use tblite_basis_type
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_xtb_h0, only : tb_hamiltonian, new_hamiltonian
   use tblite_data_covrad, only : get_covalent_rad
   use tblite_data_paulingen, only : get_pauling_en
   use tblite_ncoord_erf
   use tblite_ncoord_erf_en
   use tblite_ncoord_type, only : get_coordination_number

   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_ceh_singlepoint, only : ceh_guess
   use tblite_ceh_ceh, only : ceh_h0spec, new_ceh_calculator
   use tblite_ceh_h0, only : get_scaled_selfenergy, get_hamiltonian

   use tblite_blas, only: gemv

   use tblite_container, only : container_type, container_cache
   use tblite_external_field, only : electric_field
   implicit none
   private

   public :: collect_ceh

   real(wp), parameter :: kt = 5000.0_wp * 3.166808578545117e-06_wp
   real(wp), parameter :: thr2 = 1.0e2_wp*sqrt(epsilon(1.0_wp))

contains

   !> Collect all exported unit tests
   subroutine collect_ceh(testsuite)

      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("scaled-selfenergy-H2", test_scaled_selfenergy_h2), &
         new_unittest("scaled-selfenergy-LiH", test_scaled_selfenergy_lih), &
         new_unittest("scaled-selfenergy-S2", test_scaled_selfenergy_s2), &
         new_unittest("scaled-selfenergy-SiH4", test_scaled_selfenergy_sih4), &
         new_unittest("hamiltonian-H2", test_hamiltonian_h2), &
         new_unittest("hamiltonian-LiH", test_hamiltonian_lih), &
         new_unittest("hamiltonian-S2", test_hamiltonian_s2), &
         new_unittest("hamiltonian-SiH4", test_hamiltonian_sih4), &
         new_unittest("overlap_diat-H2", test_overlap_diat_h2), &
         new_unittest("overlap_diat-LiH", test_overlap_diat_lih), &
         new_unittest("overlap_diat-S2", test_overlap_diat_s2), &
         new_unittest("overlap_diat-SiH4", test_overlap_diat_sih4), &
         new_unittest("q-mol-h2", test_q_h2), &
         new_unittest("q-mol-lih", test_q_lih), &
         new_unittest("q-mol-1", test_q_mb01), &
         new_unittest("q-mol-2", test_q_mb02), &
         new_unittest("q-mol-3", test_q_mb03), &
         new_unittest("q-mol-4", test_q_mb04), &
         new_unittest("q-chrgd-efield-mol", test_q_ef_chrg_mb01), &
         new_unittest("d-mol", test_d_mb01), &
         new_unittest("d-field-mol", test_d_field_mb04), &
         new_unittest("d-field-change-mol", test_d_hcn) &
         ]

   end subroutine collect_ceh


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

   subroutine test_scaled_selfenergy_mol(error, mol, ref)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: ref(:)

      type(basis_type) :: bas
      type(tb_hamiltonian) :: h0
      type(erf_ncoord_type) :: ncoord
      type(erf_en_ncoord_type) :: ncoord_en
      type(adjacency_list) :: list
      real(wp), parameter :: cn_cutoff = 30.0_wp
      real(wp), allocatable :: lattr(:, :), cn(:), cn_en(:), rcov(:), en(:)
      real(wp), allocatable :: selfenergy(:)
      real(wp) :: cutoff
      integer :: ii, jj

      call make_basis(bas, mol, 6)

      call check(error, bas%nsh, size(ref, 1))
      if (allocated(error)) return

      call new_hamiltonian(h0, mol, bas, ceh_h0spec(mol))

      allocate(cn(mol%nat), cn_en(mol%nat), rcov(mol%nid), en(mol%nid))
      ! test with the standard Pyykko radii and Pauling EN (not as in CEH parametrization)
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)
      call new_erf_ncoord(ncoord, mol, cn_cutoff, rcov)
      call new_erf_en_ncoord(ncoord_en, mol, cn_cutoff, rcov)
      call get_lattice_points(mol%periodic, mol%lattice, cn_cutoff, lattr)
      call get_coordination_number(ncoord , mol, lattr, cn_cutoff, cn)
      call get_coordination_number(ncoord_en, mol, lattr, cn_cutoff, cn_en)

      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      call new_adjacency_list(list, mol, lattr, cutoff)

      allocate(selfenergy(bas%nsh))
      call get_scaled_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, cn=cn, cn_en=cn_en, &
      & selfenergy=selfenergy)

      do ii = 1, size(selfenergy, 1)
        call check(error, selfenergy(ii), ref(ii), thr=thr2)
        if (allocated(error)) then
           print '(2es20.13)', selfenergy(ii), ref(ii)
           return
        end if
      end do

   end subroutine test_scaled_selfenergy_mol


   subroutine test_hamiltonian_mol(error, mol, ref)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: ref(:, :)

      type(basis_type) :: bas
      type(tb_hamiltonian) :: h0
      type(erf_ncoord_type) :: ncoord
      type(erf_en_ncoord_type) :: ncoord_en
      type(adjacency_list) :: list
      real(wp), parameter :: cn_cutoff = 30.0_wp
      real(wp), allocatable :: lattr(:, :), cn(:), cn_en(:), rcov(:), en(:)
      real(wp), allocatable :: overlap(:, :), overlap_diat(:, :), dpint(:, :, :)
      real(wp), allocatable :: hamiltonian(:, :), selfenergy(:)
      real(wp) :: cutoff
      integer :: ii, jj

      call make_basis(bas, mol, 6)

      call check(error, bas%nao, size(ref, 1))
      if (allocated(error)) return

      call new_hamiltonian(h0, mol, bas, ceh_h0spec(mol))

      allocate(cn(mol%nat), cn_en(mol%nat), rcov(mol%nid), en(mol%nid))
      ! test with the standard Pyykko radii and Pauling EN (not as in CEH parametrization)
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)
      call new_erf_ncoord(ncoord, mol, cn_cutoff, rcov)
      call new_erf_en_ncoord(ncoord_en, mol, cn_cutoff, rcov)
      call get_lattice_points(mol%periodic, mol%lattice, cn_cutoff, lattr)
      call get_coordination_number(ncoord , mol, lattr, cn_cutoff, cn)
      call get_coordination_number(ncoord_en, mol, lattr, cn_cutoff, cn_en)

      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      call new_adjacency_list(list, mol, lattr, cutoff)

      allocate(selfenergy(bas%nsh))
      call get_scaled_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, cn=cn, cn_en=cn_en, &
      & selfenergy=selfenergy)

      allocate(overlap(bas%nao, bas%nao), overlap_diat(bas%nao, bas%nao), &
      & dpint(3, bas%nao, bas%nao), hamiltonian(bas%nao, bas%nao))

      call get_hamiltonian(mol, lattr, list, bas, h0, selfenergy, &
      & overlap, overlap_diat, dpint, hamiltonian)

      do ii = 1, size(hamiltonian, 2)
         do jj = 1, size(hamiltonian, 1)
            call check(error, hamiltonian(jj, ii), ref(jj, ii), thr=thr2)
            if (allocated(error)) then
               print '(3es21.13)', hamiltonian(jj, ii), ref(jj, ii), &
               & hamiltonian(jj, ii) - ref(jj, ii)
               return
            end if
         end do
      end do

   end subroutine test_hamiltonian_mol

   subroutine test_overlap_diat_mol(error, mol, ref)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: ref(:, :)

      type(basis_type) :: bas
      type(tb_hamiltonian) :: h0
      type(erf_ncoord_type) :: ncoord
      type(erf_en_ncoord_type) :: ncoord_en
      type(adjacency_list) :: list
      real(wp), parameter :: cn_cutoff = 30.0_wp
      real(wp), allocatable :: lattr(:, :), cn(:), cn_en(:), rcov(:), en(:)
      real(wp), allocatable :: overlap(:, :), overlap_diat(:, :), dpint(:, :, :)
      real(wp), allocatable :: hamiltonian(:, :), selfenergy(:)
      real(wp) :: cutoff
      integer :: ii, jj

      call make_basis(bas, mol, 6)

      call check(error, bas%nao, size(ref, 1))
      if (allocated(error)) return

      call new_hamiltonian(h0, mol, bas, ceh_h0spec(mol))

      allocate(cn(mol%nat), cn_en(mol%nat), rcov(mol%nid), en(mol%nid))
      ! test with the standard Pyykko radii and Pauling EN (not as in CEH parametrization)
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)
      call new_erf_ncoord(ncoord, mol, cn_cutoff, rcov)
      call new_erf_en_ncoord(ncoord_en, mol, cn_cutoff, rcov)
      call get_lattice_points(mol%periodic, mol%lattice, cn_cutoff, lattr)
      call get_coordination_number(ncoord , mol, lattr, cn_cutoff, cn)
      call get_coordination_number(ncoord_en, mol, lattr, cn_cutoff, cn_en)

      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      call new_adjacency_list(list, mol, lattr, cutoff)

      allocate(selfenergy(bas%nsh))
      call get_scaled_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, cn=cn, cn_en=cn_en, &
      & selfenergy=selfenergy)

      allocate(overlap(bas%nao, bas%nao), overlap_diat(bas%nao, bas%nao), &
      & dpint(3, bas%nao, bas%nao), hamiltonian(bas%nao, bas%nao))

      call get_hamiltonian(mol, lattr, list, bas, h0, selfenergy, &
      & overlap, overlap_diat, dpint, hamiltonian)

      do ii = 1, size(overlap_diat, 2)
         do jj = 1, size(overlap_diat, 1)
            call check(error, overlap_diat(jj, ii), ref(jj, ii), thr=thr2)
            if (allocated(error)) then
               print '(3es21.13)', overlap_diat(jj, ii), ref(jj, ii), &
               & overlap_diat(jj, ii) - ref(jj, ii)
               return
            end if
         end do
      end do

   end subroutine test_overlap_diat_mol

   subroutine test_q_gen(error, mol, ref)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type), intent(in) :: mol

      !> Reference CEH charges
      real(wp), intent(in) :: ref(:)

      type(context_type) :: ctx
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
      real(wp), allocatable :: cn(:), cn_en(:)
      real(wp), parameter :: accuracy = 1e-8_wp
      real(wp), allocatable :: lattr(:, :)
      integer :: i
      allocate(cn(mol%nat))

      call new_ceh_calculator(calc, mol)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      ctx%verbosity = 0
      call ceh_guess(ctx, calc, mol, error, wfn, accuracy)

      do i = 1, mol%nat
         call check(error, wfn%qat(i,1), ref(i), thr=1e-6_wp)
         if (allocated(error)) then
            print '(3es21.13)',  wfn%qat(i,1), ref(i), &
            & wfn%qat(i,1) - ref(i)
            return
         end if
      enddo

   end subroutine test_q_gen


   subroutine test_scaled_selfenergy_h2(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nsh = 2
      real(wp), parameter :: scaled_selfenergy(nsh) = reshape([&
      & -5.1041627058615E-01_wp, -5.1041627058615E-01_wp &
      &],shape(scaled_selfenergy))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "H2")
      call test_scaled_selfenergy_mol(error, mol, scaled_selfenergy)

   end subroutine test_scaled_selfenergy_h2

   subroutine test_scaled_selfenergy_lih(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nsh = 3
      real(wp), parameter :: scaled_selfenergy(nsh) = reshape([&
      & -3.7307764740843E-01_wp, -3.9446056938638E-01_wp, -3.3732801377653E-01_wp &
      &],shape(scaled_selfenergy))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "LiH")
      call test_scaled_selfenergy_mol(error, mol, scaled_selfenergy)

   end subroutine test_scaled_selfenergy_lih

   subroutine test_scaled_selfenergy_s2(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nsh = 6
      real(wp), parameter :: scaled_selfenergy(nsh) = reshape([&
      & -5.7642093144110E-01_wp, -5.3180811904793E-01_wp, -2.9175444046022E-01_wp, & 
      & -5.7642093144110E-01_wp, -5.3180811904793E-01_wp, -2.9175444046022E-01_wp &
      &], shape(scaled_selfenergy))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "S2")
      call test_scaled_selfenergy_mol(error, mol, scaled_selfenergy)

   end subroutine test_scaled_selfenergy_s2

   subroutine test_scaled_selfenergy_sih4(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nsh = 7
      real(wp), parameter :: scaled_selfenergy(nsh) = reshape([&
      & -5.2140420559246E-01_wp, -4.8639401015524E-01_wp, -2.4597945091348E-01_wp, & 
      & -4.7036092200061E-01_wp, -4.7036092200061E-01_wp, -4.7036092200061E-01_wp, & 
      & -4.7036092200061E-01_wp], shape(scaled_selfenergy))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "SiH4")
      call test_scaled_selfenergy_mol(error, mol, scaled_selfenergy)

   end subroutine test_scaled_selfenergy_sih4


   subroutine test_hamiltonian_h2(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 2
      real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -5.1041627058615E-01_wp, -3.9128731178966E-01_wp, -3.9128731178966E-01_wp, & 
      & -5.1041627058615E-01_wp],shape(hamiltonian))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "H2")
      call test_hamiltonian_mol(error, mol, hamiltonian)

   end subroutine test_hamiltonian_h2

   subroutine test_hamiltonian_lih(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 5
      real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -3.7307764740843E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, -1.7777739269090E-01_wp, +0.0000000000000E+00_wp, & 
      & -3.9446056938638E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -3.9446056938638E-01_wp, +0.0000000000000E+00_wp, -2.3777030579175E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -3.9446056938638E-01_wp, +0.0000000000000E+00_wp, -1.7777739269090E-01_wp, & 
      & +0.0000000000000E+00_wp, -2.3777030579175E-01_wp, +0.0000000000000E+00_wp, & 
      & -3.3732801377653E-01_wp],shape(hamiltonian))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "LiH")
      call test_hamiltonian_mol(error, mol, hamiltonian)

   end subroutine test_hamiltonian_lih

   subroutine test_hamiltonian_s2(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 18
      real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -5.7642093144110E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -1.8327174566928E-01_wp, +0.0000000000000E+00_wp, +2.2661674580606E-01_wp, & 
      & +0.0000000000000E+00_wp, -2.1711090569016E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, -5.3180811904793E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, -1.1516226432909E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +2.1604179999700E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -5.3180811904793E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -2.2661674580606E-01_wp, +0.0000000000000E+00_wp, +2.6953901229836E-01_wp, & 
      & +0.0000000000000E+00_wp, -1.8797977314300E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -5.3180811904793E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -1.1516226432909E-01_wp, +0.0000000000000E+00_wp, +2.1604179999700E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, -2.9175444046022E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -2.1711090569016E-01_wp, +0.0000000000000E+00_wp, +1.8797977314300E-01_wp, & 
      & +0.0000000000000E+00_wp, -7.5823435266772E-02_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -2.9175444046022E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -2.1604179999700E-01_wp, +0.0000000000000E+00_wp, +2.3389877153471E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -2.9175444046022E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, -2.1604179999700E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +2.3389877153471E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, -2.9175444046022E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, -1.6808397034405E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -2.9175444046022E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -1.6808397034405E-01_wp, & 
      & -1.8327174566928E-01_wp, +0.0000000000000E+00_wp, -2.2661674580606E-01_wp, & 
      & +0.0000000000000E+00_wp, -2.1711090569016E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -5.7642093144110E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, -1.1516226432909E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -2.1604179999700E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, -5.3180811904793E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +2.2661674580606E-01_wp, +0.0000000000000E+00_wp, +2.6953901229836E-01_wp, & 
      & +0.0000000000000E+00_wp, +1.8797977314300E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -5.3180811904793E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -1.1516226432909E-01_wp, +0.0000000000000E+00_wp, -2.1604179999700E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -5.3180811904793E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -2.1711090569016E-01_wp, +0.0000000000000E+00_wp, -1.8797977314300E-01_wp, & 
      & +0.0000000000000E+00_wp, -7.5823435266772E-02_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, -2.9175444046022E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +2.1604179999700E-01_wp, +0.0000000000000E+00_wp, +2.3389877153471E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -2.9175444046022E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +2.1604179999700E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +2.3389877153471E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -2.9175444046022E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, -1.6808397034405E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, -2.9175444046022E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -1.6808397034405E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -2.9175444046022E-01_wp],&
      & shape(hamiltonian))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "S2")
      call test_hamiltonian_mol(error, mol, hamiltonian)

   end subroutine test_hamiltonian_s2

   subroutine test_hamiltonian_sih4(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 13
      real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -5.2140420559246E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -2.2761846976247E-01_wp, -2.2761846976247E-01_wp, -2.2761846976247E-01_wp, & 
      & -2.2761846976247E-01_wp, +0.0000000000000E+00_wp, -4.8639401015524E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, -1.6954715848881E-01_wp, +1.6954715848881E-01_wp, & 
      & +1.6954715848881E-01_wp, -1.6954715848881E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, -4.8639401015524E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +1.6954715848881E-01_wp, & 
      & +1.6954715848881E-01_wp, -1.6954715848881E-01_wp, -1.6954715848881E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -4.8639401015524E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -1.6954715848881E-01_wp, +1.6954715848881E-01_wp, -1.6954715848881E-01_wp, & 
      & +1.6954715848881E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -2.4597945091348E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, -1.0301163014236E-17_wp, -1.0301163014236E-17_wp, & 
      & -1.0301163014236E-17_wp, -1.0301163014236E-17_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, -2.4597945091348E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +1.5146408363682E-01_wp, & 
      & -1.5146408363682E-01_wp, -1.5146408363682E-01_wp, +1.5146408363682E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -2.4597945091348E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +1.5146408363682E-01_wp, -1.5146408363682E-01_wp, +1.5146408363682E-01_wp, & 
      & -1.5146408363682E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -2.4597945091348E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, -2.4597945091348E-01_wp, -1.5146408363682E-01_wp, & 
      & -1.5146408363682E-01_wp, +1.5146408363682E-01_wp, +1.5146408363682E-01_wp, & 
      & -2.2761846976247E-01_wp, -1.6954715848881E-01_wp, +1.6954715848881E-01_wp, & 
      & -1.6954715848881E-01_wp, -1.0301163014236E-17_wp, +1.5146408363682E-01_wp, & 
      & +1.5146408363682E-01_wp, +0.0000000000000E+00_wp, -1.5146408363682E-01_wp, & 
      & -4.7036092200061E-01_wp, -3.3693071890466E-02_wp, -3.3693071890466E-02_wp, & 
      & -3.3693071890466E-02_wp, -2.2761846976247E-01_wp, +1.6954715848881E-01_wp, & 
      & +1.6954715848881E-01_wp, +1.6954715848881E-01_wp, -1.0301163014236E-17_wp, & 
      & -1.5146408363682E-01_wp, -1.5146408363682E-01_wp, +0.0000000000000E+00_wp, & 
      & -1.5146408363682E-01_wp, -3.3693071890466E-02_wp, -4.7036092200061E-01_wp, & 
      & -3.3693071890466E-02_wp, -3.3693071890466E-02_wp, -2.2761846976247E-01_wp, & 
      & +1.6954715848881E-01_wp, -1.6954715848881E-01_wp, -1.6954715848881E-01_wp, & 
      & -1.0301163014236E-17_wp, -1.5146408363682E-01_wp, +1.5146408363682E-01_wp, & 
      & +0.0000000000000E+00_wp, +1.5146408363682E-01_wp, -3.3693071890466E-02_wp, & 
      & -3.3693071890466E-02_wp, -4.7036092200061E-01_wp, -3.3693071890466E-02_wp, & 
      & -2.2761846976247E-01_wp, -1.6954715848881E-01_wp, -1.6954715848881E-01_wp, & 
      & +1.6954715848881E-01_wp, -1.0301163014236E-17_wp, +1.5146408363682E-01_wp, & 
      & -1.5146408363682E-01_wp, +0.0000000000000E+00_wp, +1.5146408363682E-01_wp, & 
      & -3.3693071890466E-02_wp, -3.3693071890466E-02_wp, -3.3693071890466E-02_wp, & 
      & -4.7036092200061E-01_wp],shape(hamiltonian))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "SiH4")
      call test_hamiltonian_mol(error, mol, hamiltonian)

   end subroutine test_hamiltonian_sih4


   subroutine test_overlap_diat_h2(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 2
      real(wp), parameter :: overlap_diat(nao, nao) = reshape([&
      & +9.9999999988150E-01_wp, +1.1625785680539E+00_wp, +1.1625785680539E+00_wp, & 
      & +9.9999999988150E-01_wp],shape(overlap_diat))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "H2")
      call test_overlap_diat_mol(error, mol, overlap_diat)

   end subroutine test_overlap_diat_h2

   subroutine test_overlap_diat_lih(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 5
      real(wp), parameter :: overlap_diat(nao, nao) = reshape([&
      & +1.0000000000060E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +7.5901642084662E-01_wp, +0.0000000000000E+00_wp, & 
      & +9.9999999992569E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +9.9999999992569E-01_wp, +0.0000000000000E+00_wp, +8.8587462290517E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +9.9999999992569E-01_wp, +0.0000000000000E+00_wp, +7.5901642084662E-01_wp, & 
      & +0.0000000000000E+00_wp, +8.8587462290517E-01_wp, +0.0000000000000E+00_wp, & 
      & +9.9999999988150E-01_wp],shape(overlap_diat))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "LiH")
      call test_overlap_diat_mol(error, mol, overlap_diat)

   end subroutine test_overlap_diat_lih

   subroutine test_overlap_diat_s2(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 18
      real(wp), parameter :: overlap_diat(nao, nao) = reshape([&
      & +9.9999999986933E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +4.8217737784357E-01_wp, +0.0000000000000E+00_wp, -5.5752294134525E-01_wp, & 
      & +0.0000000000000E+00_wp, +5.5764799996744E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +9.9999999999806E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +2.6810514555495E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -5.4029331669609E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +9.9999999999806E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +5.5752294134525E-01_wp, +0.0000000000000E+00_wp, -6.2750412685952E-01_wp, & 
      & +0.0000000000000E+00_wp, +4.7011372384704E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +9.9999999999806E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +2.6810514555495E-01_wp, +0.0000000000000E+00_wp, -5.4029331669609E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +5.5764799996744E-01_wp, +0.0000000000000E+00_wp, -4.7011372384704E-01_wp, & 
      & +0.0000000000000E+00_wp, +2.2909719098522E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +5.4029331669609E-01_wp, +0.0000000000000E+00_wp, -7.0671490080821E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +9.9999999983021E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +5.4029331669609E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -7.0671490080821E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +5.0785835962169E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +5.0785835962169E-01_wp, & 
      & +4.8217737784357E-01_wp, +0.0000000000000E+00_wp, +5.5752294134525E-01_wp, & 
      & +0.0000000000000E+00_wp, +5.5764799996744E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +9.9999999986933E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +2.6810514555495E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +5.4029331669609E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +9.9999999999806E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -5.5752294134525E-01_wp, +0.0000000000000E+00_wp, -6.2750412685952E-01_wp, & 
      & +0.0000000000000E+00_wp, -4.7011372384704E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +9.9999999999806E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +2.6810514555495E-01_wp, +0.0000000000000E+00_wp, +5.4029331669609E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +9.9999999999806E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +5.5764799996744E-01_wp, +0.0000000000000E+00_wp, +4.7011372384704E-01_wp, & 
      & +0.0000000000000E+00_wp, +2.2909719098522E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -5.4029331669609E-01_wp, +0.0000000000000E+00_wp, -7.0671490080821E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, -5.4029331669609E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -7.0671490080821E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +9.9999999983021E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +5.0785835962169E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +5.0785835962169E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +9.9999999983021E-01_wp],&
      & shape(overlap_diat))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "S2")
      call test_overlap_diat_mol(error, mol, overlap_diat)

   end subroutine test_overlap_diat_s2

   subroutine test_overlap_diat_sih4(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 13
      real(wp), parameter :: overlap_diat(nao, nao) = reshape([&
      & +9.9999999986933E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +6.9611295875884E-01_wp, +6.9611295875884E-01_wp, +6.9611295875884E-01_wp, & 
      & +6.9611295875884E-01_wp, +0.0000000000000E+00_wp, +9.9999999999806E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +4.8315900715143E-01_wp, -4.8315900715143E-01_wp, & 
      & -4.8315900715143E-01_wp, +4.8315900715143E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +9.9999999999806E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -4.8315900715143E-01_wp, & 
      & -4.8315900715143E-01_wp, +4.8315900715143E-01_wp, +4.8315900715143E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +9.9999999999806E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +4.8315900715143E-01_wp, -4.8315900715143E-01_wp, +4.8315900715143E-01_wp, & 
      & -4.8315900715143E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +3.2066591762028E-17_wp, +3.2066591762028E-17_wp, & 
      & +3.2066591762028E-17_wp, +3.2066591762028E-17_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -4.7149403711788E-01_wp, & 
      & +4.7149403711788E-01_wp, +4.7149403711788E-01_wp, -4.7149403711788E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +9.9999999983021E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & -4.7149403711788E-01_wp, +4.7149403711788E-01_wp, -4.7149403711788E-01_wp, & 
      & +4.7149403711788E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, & 
      & +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, +4.7149403711788E-01_wp, & 
      & +4.7149403711788E-01_wp, -4.7149403711788E-01_wp, -4.7149403711788E-01_wp, & 
      & +6.9611295875884E-01_wp, +4.8315900715143E-01_wp, -4.8315900715143E-01_wp, & 
      & +4.8315900715143E-01_wp, +3.2066591762028E-17_wp, -4.7149403711788E-01_wp, & 
      & -4.7149403711788E-01_wp, +0.0000000000000E+00_wp, +4.7149403711788E-01_wp, & 
      & +9.9999999988150E-01_wp, +1.0863266473352E-01_wp, +1.0863266473352E-01_wp, & 
      & +1.0863266473352E-01_wp, +6.9611295875884E-01_wp, -4.8315900715143E-01_wp, & 
      & -4.8315900715143E-01_wp, -4.8315900715143E-01_wp, +3.2066591762028E-17_wp, & 
      & +4.7149403711788E-01_wp, +4.7149403711788E-01_wp, +0.0000000000000E+00_wp, & 
      & +4.7149403711788E-01_wp, +1.0863266473352E-01_wp, +9.9999999988150E-01_wp, & 
      & +1.0863266473352E-01_wp, +1.0863266473352E-01_wp, +6.9611295875884E-01_wp, & 
      & -4.8315900715143E-01_wp, +4.8315900715143E-01_wp, +4.8315900715143E-01_wp, & 
      & +3.2066591762028E-17_wp, +4.7149403711788E-01_wp, -4.7149403711788E-01_wp, & 
      & +0.0000000000000E+00_wp, -4.7149403711788E-01_wp, +1.0863266473352E-01_wp, & 
      & +1.0863266473352E-01_wp, +9.9999999988150E-01_wp, +1.0863266473352E-01_wp, & 
      & +6.9611295875884E-01_wp, +4.8315900715143E-01_wp, +4.8315900715143E-01_wp, & 
      & -4.8315900715143E-01_wp, +3.2066591762028E-17_wp, -4.7149403711788E-01_wp, & 
      & +4.7149403711788E-01_wp, +0.0000000000000E+00_wp, -4.7149403711788E-01_wp, & 
      & +1.0863266473352E-01_wp, +1.0863266473352E-01_wp, +1.0863266473352E-01_wp, & 
      & +9.9999999988150E-01_wp],shape(overlap_diat))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "SiH4")
      call test_overlap_diat_mol(error, mol, overlap_diat)

   end subroutine test_overlap_diat_sih4


   subroutine test_q_h2(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(2) = reshape([ & 
      & 0.0000000000000_wp, 0.00000000000000_wp &
      &], shape(charges))

      call get_structure(mol, "MB16-43", "H2")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_h2


   subroutine test_q_lih(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(2) = reshape([ &
      & 0.383963042442916_wp, -0.383963042442917_wp &
      &], shape(charges))

      call get_structure(mol, "MB16-43", "LiH")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_lih

   subroutine test_q_mb01(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(16) = reshape([ &
      & 0.504694175287870_wp, -0.049320603215405_wp, -0.442887082747601_wp, &     
      &-0.040897715947291_wp, -0.225355610004319_wp,  0.079338661074347_wp, &
      &-0.012184763492912_wp, -0.354104643499813_wp, -0.225243243808309_wp, &     
      & 0.087276255050057_wp,  0.085471100209451_wp, -0.038779744241030_wp, &
      & 0.038616630415300_wp,  0.127389947162501_wp, -0.039323802393347_wp, &
      & 0.505310440150523_wp], shape(charges))

      call get_structure(mol, "MB16-43", "01")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_mb01

   subroutine test_q_mb02(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(16) = reshape([ &
      &-0.059404506571487_wp, -0.118702191502332_wp, -0.195235491894287_wp, &
      &-0.210021399475394_wp,  0.496594666771036_wp,  0.158190196939209_wp, &
      &-0.075318639484230_wp, -0.045052163348822_wp,  0.286946319499636_wp, &
      & 0.157881990181518_wp, -0.029236175252485_wp,  0.281560060657950_wp, &
      &-0.312312030061058_wp, -0.064265894380938_wp,  0.066269279024472_wp, &
      &-0.337894021272329_wp], shape(charges))

      call get_structure(mol, "MB16-43", "02")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_mb02

   subroutine test_q_mb03(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(16) = reshape([ &
      & 0.047060647546981_wp, -0.452721646051320_wp,  0.000750024797478_wp, &
      & 0.303998197523052_wp,  0.479998103945730_wp,  0.099926859160474_wp, &
      &-0.229477809903280_wp,  0.016858971022812_wp,  0.008564268752175_wp, &
      & 0.004606322535074_wp, -0.354502824466586_wp, -0.187631548264506_wp, &     
      &-0.300801554641486_wp,  0.041952785802619_wp,  0.482401552816662_wp, &     
      & 0.039017649441556_wp], shape(charges))

      call get_structure(mol, "MB16-43", "03")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_mb03

   subroutine test_q_mb04(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(16) = reshape([ &
      & 0.017867676479831_wp, -0.103259042275922_wp, -0.035078199415400_wp, &
      &-0.183433026064262_wp,  0.263017445728753_wp,  0.012000420842077_wp, &
      & 0.040327964255115_wp, -0.119339642424892_wp, -0.010898440812181_wp, &
      &-0.056108839369092_wp, -0.137112266394061_wp,  0.334034245022961_wp, &
      & 0.147772637227760_wp, -0.242482916232017_wp,  0.033049020618869_wp, &
      & 0.039642962813086_wp], shape(charges))

      call get_structure(mol, "MB16-43", "04")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_mb04


   subroutine test_q_ef_chrg_mb01(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
      real(wp), parameter :: accuracy = 1e-8_wp
      class(container_type), allocatable :: cont      
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: ref(16) = reshape([ &
      &-6.90111445726635_wp, -0.99979611590118_wp,   4.2837802019751_wp, &
      &-0.99980021634224_wp,  7.00023420054907_wp,   0.9454454365410_wp, &
      &-0.98804624196187_wp,  5.84908677934913_wp,   2.2839742642425_wp, & 
      & 0.95737100853797_wp,  0.99908251026128_wp, -10.9475535476315_wp, &
      &-4.99687086866152_wp,  2.83244134467160_wp,   4.8453679895474_wp, &
      &-2.16360228791069_wp], shape(ref))

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
      call ceh_guess(ctx, calc, mol, error, wfn, accuracy)

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
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
      real(wp) :: dipole(3), tmp(3)
      real(wp), parameter :: accuracy = 1e-8_wp
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: ref(3) = reshape([ &
      & 0.438025381623586_wp, -0.735884148841272_wp, -2.76541717331434_wp &
      &], shape(ref))

      call get_structure(mol, "MB16-43", "01")

      call new_ceh_calculator(calc, mol)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      ctx%verbosity = 0
      call ceh_guess(ctx, calc, mol, error, wfn, accuracy)
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

   subroutine test_d_field_mb04(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
      class(container_type), allocatable :: cont
      real(wp) :: energy, efield(3), dipole(3), tmp(3)
      real(wp), parameter :: accuracy = 1e-8_wp
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: ref(3) = reshape([ & 
         & -20.835631606789_wp,  97.0349021135889_wp, -7.70521258527074_wp &
         ], shape(ref))

      call get_structure(mol, "MB16-43", "04")
      energy = 0.0_wp
      efield = 0.0_wp
      efield(2) = 0.2_wp

      call new_ceh_calculator(calc, mol)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)

      cont = electric_field(efield)
      call calc%push_back(cont)

      ctx%verbosity = 0
      call ceh_guess(ctx, calc, mol, error, wfn, accuracy)
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

   end subroutine test_d_field_mb04

   subroutine test_d_hcn(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(structure_type) :: mol1,mol2
      type(xtb_calculator) :: calc1,calc2
      type(wavefunction_type) :: wfn1,wfn2
      class(container_type), allocatable :: cont1,cont2
      real(wp), parameter :: accuracy = 1e-8_wp
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
      call ceh_guess(ctx, calc1, mol1, error, wfn1, accuracy)
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
      call ceh_guess(ctx, calc2, mol2, error, wfn2, accuracy)
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
