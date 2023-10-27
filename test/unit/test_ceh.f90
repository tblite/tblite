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
   use tblite_data_covrad_ceh, only : get_covalent_cehrad
   use tblite_data_paulingen_ceh, only : get_pauling_en_ceh
   use tblite_ncoord_ceh_std
   use tblite_ncoord_ceh_en
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

   real(wp), parameter :: kt = 300.0_wp * 3.166808578545117e-06_wp
   real(wp), parameter :: thr2 = 1.0e2_wp*sqrt(epsilon(1.0_wp))

contains

   !> Collect all exported unit tests
   subroutine collect_ceh(testsuite)

      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("scaled-selfenergy-1", test_scaled_selfenergy_h2), &
         new_unittest("scaled-selfenergy-2", test_scaled_selfenergy_lih), &
         new_unittest("scaled-selfenergy-3", test_scaled_selfenergy_s2), &
         new_unittest("scaled-selfenergy-4", test_scaled_selfenergy_sih4), &
         new_unittest("hamiltonian-1", test_hamiltonian_h2), &
         new_unittest("hamiltonian-2", test_hamiltonian_lih), &
         new_unittest("hamiltonian-3", test_hamiltonian_s2), &
         new_unittest("hamiltonian-4", test_hamiltonian_sih4), &
         new_unittest("overlap_diat-1", test_overlap_diat_h2), &
         new_unittest("overlap_diat-2", test_overlap_diat_lih), &
         new_unittest("overlap_diat-3", test_overlap_diat_s2), &
         new_unittest("overlap_diat-4", test_overlap_diat_sih4), &
         new_unittest("q-mol-1", test_q_mb01), &
         new_unittest("q-mol-2", test_q_mb02), &
         new_unittest("q-mol-3", test_q_mb03), &
         new_unittest("q-mol-4", test_q_mb04), &
         new_unittest("q-chrgd-efield-mol", test_q_ef_chrg_mb01), &
         new_unittest("d-mol", test_d_mb01), &
         new_unittest("d-field-mol", test_d_mb04), &
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
   type(ceh_std_ncoord_type) :: ncoord
   type(ceh_en_ncoord_type) :: ncoord_en
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
   rcov(:) = get_covalent_cehrad(mol%num)
   en(:) = get_pauling_en_ceh(mol%num)
   call new_ceh_std_ncoord(ncoord, mol, cn_cutoff, rcov)
   call new_ceh_en_ncoord(ncoord_en, mol, cn_cutoff, rcov)
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
   type(ceh_std_ncoord_type) :: ncoord
   type(ceh_en_ncoord_type) :: ncoord_en
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
   rcov(:) = get_covalent_cehrad(mol%num)
   en(:) = get_pauling_en_ceh(mol%num)
   call new_ceh_std_ncoord(ncoord, mol, cn_cutoff, rcov)
   call new_ceh_en_ncoord(ncoord_en, mol, cn_cutoff, rcov)
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
   type(ceh_std_ncoord_type) :: ncoord
   type(ceh_en_ncoord_type) :: ncoord_en
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
   rcov(:) = get_covalent_cehrad(mol%num)
   en(:) = get_pauling_en_ceh(mol%num)
   call new_ceh_std_ncoord(ncoord, mol, cn_cutoff, rcov)
   call new_ceh_en_ncoord(ncoord_en, mol, cn_cutoff, rcov)
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
      &-5.02531612491150E-1_wp,-5.02531612491150E-1_wp &
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
      &-3.95120002130243E-1_wp,-1.68766212568570E-1_wp,-5.00535900276072E-1_wp &
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
      &-5.49096214933834E-1_wp,-4.93891369275332E-1_wp,-1.52761590706418E-1_wp,&
      &-5.49096214933834E-1_wp,-4.93891369275332E-1_wp,-1.52761590706418E-1_wp &
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
      &-5.87966583817092E-1_wp,-4.27252582137204E-1_wp,-1.16479369796453E-1_wp,&
      &-5.06480584862866E-1_wp,-5.06480584862866E-1_wp,-5.06480584862866E-1_wp,&
      &-5.06480584862866E-1_wp], shape(scaled_selfenergy))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "SiH4")
   call test_scaled_selfenergy_mol(error, mol, scaled_selfenergy)

end subroutine test_scaled_selfenergy_sih4


subroutine test_hamiltonian_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 2
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      &-5.0253161249115E-1_wp, -4.8328023191842E-1_wp, -4.8328023191842E-1_wp, &
      &-5.0253161249115E-1_wp],shape(hamiltonian))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "H2")
   call test_hamiltonian_mol(error, mol, hamiltonian)

end subroutine test_hamiltonian_h2

subroutine test_hamiltonian_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 5
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      &-3.9512000213024E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp,-2.9461923997699E-1_wp, 0.0000000000000E+0_wp, &
      &-1.6876621256857E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      &-1.6876621256857E-1_wp, 0.0000000000000E+0_wp,-3.2190465369231E-1_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      &-1.6876621256857E-1_wp, 0.0000000000000E+0_wp,-2.9461923997699E-1_wp, &
      & 0.0000000000000E+0_wp,-3.2190465369231E-1_wp, 0.0000000000000E+0_wp, &
      &-5.0053590027607E-1_wp],shape(hamiltonian))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_hamiltonian_mol(error, mol, hamiltonian)

end subroutine test_hamiltonian_lih

subroutine test_hamiltonian_s2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 18
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      &-5.4909621493383E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      &-1.5557648632663E-1_wp, 0.0000000000000E+0_wp, 2.1402505681592E-1_wp, & 
      & 0.0000000000000E+0_wp,-1.6876705809971E-1_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp,-4.9389136927533E-1_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp,-1.1575948410446E-1_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 1.7451690255052E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp,-4.9389136927533E-1_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      &-2.1402505681592E-1_wp, 0.0000000000000E+0_wp, 2.7416806818271E-1_wp, & 
      & 0.0000000000000E+0_wp,-1.5365959440037E-1_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      &-4.9389136927533E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      &-1.1575948410446E-1_wp, 0.0000000000000E+0_wp, 1.7451690255052E-1_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp,-1.5276159070642E-1_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      &-1.6876705809971E-1_wp, 0.0000000000000E+0_wp, 1.5365959440037E-1_wp, & 
      & 0.0000000000000E+0_wp,-3.9798498518230E-2_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp,-1.5276159070642E-1_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      &-1.7451690255052E-1_wp, 0.0000000000000E+0_wp, 1.2132277062713E-1_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      &-1.5276159070642E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp,-1.7451690255052E-1_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 1.2132277062713E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp,-1.5276159070642E-1_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp,-1.4337015160658E-1_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp,-1.5276159070642E-1_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp,-1.4337015160658E-1_wp, &
      &-1.5557648632663E-1_wp, 0.0000000000000E+0_wp,-2.1402505681592E-1_wp, & 
      & 0.0000000000000E+0_wp,-1.6876705809971E-1_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      &-5.4909621493383E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp,-1.1575948410446E-1_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      &-1.7451690255052E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp,-4.9389136927533E-1_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 2.1402505681592E-1_wp, 0.0000000000000E+0_wp, 2.7416806818271E-1_wp, & 
      & 0.0000000000000E+0_wp, 1.5365959440037E-1_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp,-4.9389136927533E-1_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      &-1.1575948410446E-1_wp, 0.0000000000000E+0_wp,-1.7451690255052E-1_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      &-4.9389136927533E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      &-1.6876705809971E-1_wp, 0.0000000000000E+0_wp,-1.5365959440037E-1_wp, & 
      & 0.0000000000000E+0_wp,-3.9798498518230E-2_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp,-1.5276159070642E-1_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 1.7451690255052E-1_wp, 0.0000000000000E+0_wp, 1.2132277062713E-1_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp,-1.5276159070642E-1_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 1.7451690255052E-1_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 1.2132277062713E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      &-1.5276159070642E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp,-1.4337015160658E-1_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp,-1.5276159070642E-1_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp,-1.4337015160658E-1_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, & 
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp,-1.5276159070642E-1_wp],&
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
      &-5.8796658381709E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, & 
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, & 
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      &-2.8957077919196E-01_wp,-2.8957077919196E-01_wp,-2.8957077919196E-01_wp, &
      &-2.8957077919196E-01_wp, 0.0000000000000E+00_wp,-4.2725258213720E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, & 
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, & 
      & 0.0000000000000E+00_wp,-2.1481130685659E-01_wp, 2.1481130685659E-01_wp, & 
      & 2.1481130685659E-01_wp,-2.1481130685659E-01_wp, 0.0000000000000E+00_wp, & 
      & 0.0000000000000E+00_wp,-4.2725258213720E-01_wp, 0.0000000000000E+00_wp, & 
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, & 
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 2.1481130685659E-01_wp, & 
      & 2.1481130685659E-01_wp,-2.1481130685659E-01_wp,-2.1481130685659E-01_wp, & 
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      &-4.2725258213720E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, & 
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      &-2.1481130685659E-01_wp, 2.1481130685659E-01_wp,-2.1481130685659E-01_wp, & 
      & 2.1481130685659E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, & 
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp,-1.1647936979645E-01_wp, & 
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, & 
      & 0.0000000000000E+00_wp,-1.4453496267419E-17_wp,-1.4453496267419E-17_wp, &
      &-1.4453496267419E-17_wp,-1.4453496267419E-17_wp, 0.0000000000000E+00_wp, & 
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, & 
      & 0.0000000000000E+00_wp,-1.1647936979645E-01_wp, 0.0000000000000E+00_wp, & 
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.6384615846636E-01_wp, &
      &-1.6384615846636E-01_wp,-1.6384615846636E-01_wp, 1.6384615846636E-01_wp, & 
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, & 
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      &-1.1647936979645E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, & 
      & 1.6384615846636E-01_wp,-1.6384615846636E-01_wp, 1.6384615846636E-01_wp, &
      &-1.6384615846636E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, & 
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, & 
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp,-1.1647936979645E-01_wp, & 
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, & 
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, & 
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, & 
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, & 
      & 0.0000000000000E+00_wp,-1.1647936979645E-01_wp,-1.6384615846636E-01_wp, &
      &-1.6384615846636E-01_wp, 1.6384615846636E-01_wp, 1.6384615846636E-01_wp, &
      &-2.8957077919196E-01_wp,-2.1481130685659E-01_wp, 2.1481130685659E-01_wp, &
      &-2.1481130685659E-01_wp,-1.4453496267419E-17_wp, 1.6384615846636E-01_wp, & 
      & 1.6384615846636E-01_wp, 0.0000000000000E+00_wp,-1.6384615846636E-01_wp, &
      &-5.0648058486287E-01_wp,-4.5513115639586E-02_wp,-4.5513115639586E-02_wp, &
      &-4.5513115639586E-02_wp,-2.8957077919196E-01_wp, 2.1481130685659E-01_wp, & 
      & 2.1481130685659E-01_wp, 2.1481130685659E-01_wp,-1.4453496267419E-17_wp, &
      &-1.6384615846636E-01_wp,-1.6384615846636E-01_wp, 0.0000000000000E+00_wp, &
      &-1.6384615846636E-01_wp,-4.5513115639586E-02_wp,-5.0648058486287E-01_wp, &
      &-4.5513115639586E-02_wp,-4.5513115639586E-02_wp,-2.8957077919196E-01_wp, & 
      & 2.1481130685659E-01_wp,-2.1481130685659E-01_wp,-2.1481130685659E-01_wp, &
      &-1.4453496267419E-17_wp,-1.6384615846636E-01_wp, 1.6384615846636E-01_wp, & 
      & 0.0000000000000E+00_wp, 1.6384615846636E-01_wp,-4.5513115639586E-02_wp, &
      &-4.5513115639586E-02_wp,-5.0648058486287E-01_wp,-4.5513115639586E-02_wp, &
      &-2.8957077919196E-01_wp,-2.1481130685659E-01_wp,-2.1481130685659E-01_wp, & 
      & 2.1481130685659E-01_wp,-1.4453496267419E-17_wp, 1.6384615846636E-01_wp, &
      &-1.6384615846636E-01_wp, 0.0000000000000E+00_wp, 1.6384615846636E-01_wp, &
      &-4.5513115639586E-02_wp,-4.5513115639586E-02_wp,-4.5513115639586E-02_wp, &
      &-5.0648058486287E-01_wp],shape(hamiltonian))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "SiH4")
   call test_hamiltonian_mol(error, mol, hamiltonian)

end subroutine test_hamiltonian_sih4


subroutine test_overlap_diat_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 2
   real(wp), parameter :: overlap_diat(nao, nao) = reshape([&
      & 9.9999999988150E-1_wp, 1.5106679310445E+0_wp, 1.5106679310445E+0_wp, &
      & 9.9999999988150E-1_wp],shape(overlap_diat))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "H2")
   call test_overlap_diat_mol(error, mol, overlap_diat)

end subroutine test_overlap_diat_h2

subroutine test_overlap_diat_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 5
   real(wp), parameter :: overlap_diat(nao, nao) = reshape([&
      & 1.0000000000060E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 1.0334352569168E+0_wp, 0.0000000000000E+0_wp, &
      & 9.9999999992569E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 9.9999999992569E-1_wp, 0.0000000000000E+0_wp, 1.2061584484521E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 9.9999999992569E-1_wp, 0.0000000000000E+0_wp, 1.0334352569168E+0_wp, &
      & 0.0000000000000E+0_wp, 1.2061584484521E+0_wp, 0.0000000000000E+0_wp, &
      & 9.9999999988150E-1_wp,-3.9512000213024E-1_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp],shape(overlap_diat))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_overlap_diat_mol(error, mol, overlap_diat)

end subroutine test_overlap_diat_lih

subroutine test_overlap_diat_s2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 18
   real(wp), parameter :: overlap_diat(nao, nao) = reshape([&
      & 9.9999999986933E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 4.4507056660994E-1_wp, 0.0000000000000E+0_wp,-5.1461777927516E-1_wp, &
      & 0.0000000000000E+0_wp, 5.1473321379033E-1_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 9.9999999999806E-1_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 2.4455601027678E-1_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      &-4.9283641176265E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 9.9999999999806E-1_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 5.1461777927516E-1_wp, 0.0000000000000E+0_wp,-5.7921343912997E-1_wp, &
      & 0.0000000000000E+0_wp, 4.3393529240103E-1_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 9.9999999999806E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 2.4455601027678E-1_wp, 0.0000000000000E+0_wp,-4.9283641176265E-1_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 9.9999999983021E-1_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 5.1473321379033E-1_wp, 0.0000000000000E+0_wp,-4.3393529240103E-1_wp, &
      & 0.0000000000000E+0_wp, 2.1146661225906E-1_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 9.9999999983021E-1_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 4.9283641176265E-1_wp, 0.0000000000000E+0_wp,-6.4464028165914E-1_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 9.9999999983021E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 4.9283641176265E-1_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      &-6.4464028165914E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 9.9999999983021E-1_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 7.6178753943254E-1_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 9.9999999983021E-1_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 7.6178753943254E-1_wp, &
      & 4.4507056660994E-1_wp, 0.0000000000000E+0_wp, 5.1461777927516E-1_wp, &
      & 0.0000000000000E+0_wp, 5.1473321379033E-1_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 9.9999999986933E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 2.4455601027678E-1_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 4.9283641176265E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 9.9999999999806E-1_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      &-5.1461777927516E-1_wp, 0.0000000000000E+0_wp,-5.7921343912997E-1_wp, &
      & 0.0000000000000E+0_wp,-4.3393529240103E-1_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 9.9999999999806E-1_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 2.4455601027678E-1_wp, 0.0000000000000E+0_wp, 4.9283641176265E-1_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 9.9999999999806E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 5.1473321379033E-1_wp, 0.0000000000000E+0_wp, 4.3393529240103E-1_wp, &
      & 0.0000000000000E+0_wp, 2.1146661225906E-1_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 9.9999999983021E-1_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      &-4.9283641176265E-1_wp, 0.0000000000000E+0_wp,-6.4464028165914E-1_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 9.9999999983021E-1_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp,-4.9283641176265E-1_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      &-6.4464028165914E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 9.9999999983021E-1_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 7.6178753943254E-1_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 9.9999999983021E-1_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 7.6178753943254E-1_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, &
      & 0.0000000000000E+0_wp, 0.0000000000000E+0_wp, 9.9999999983021E-1_wp],&
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
      & 9.9999999986933E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 8.3123399075353E-01_wp, 8.3123399075353E-01_wp, 8.3123399075353E-01_wp, &
      & 8.3123399075353E-01_wp, 0.0000000000000E+00_wp, 9.9999999999806E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 5.7694399253689E-01_wp,-5.7694399253689E-01_wp, &
      &-5.7694399253689E-01_wp, 5.7694399253689E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 9.9999999999806E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp,-5.7694399253689E-01_wp, &
      &-5.7694399253689E-01_wp, 5.7694399253689E-01_wp, 5.7694399253689E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 9.9999999999806E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 5.7694399253689E-01_wp,-5.7694399253689E-01_wp, 5.7694399253689E-01_wp, &
      &-5.7694399253689E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 9.9999999983021E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 4.9665685187358E-17_wp, 4.9665685187358E-17_wp, &
      & 4.9665685187358E-17_wp, 4.9665685187358E-17_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 9.9999999983021E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp,-5.6301475954244E-01_wp, &
      & 5.6301475954244E-01_wp, 5.6301475954244E-01_wp,-5.6301475954244E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 9.9999999983021E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      &-5.6301475954244E-01_wp, 5.6301475954244E-01_wp,-5.6301475954244E-01_wp, &
      & 5.6301475954244E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 9.9999999983021E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 9.9999999983021E-01_wp, 5.6301475954244E-01_wp, &
      & 5.6301475954244E-01_wp,-5.6301475954244E-01_wp,-5.6301475954244E-01_wp, &
      & 8.3123399075353E-01_wp, 5.7694399253689E-01_wp,-5.7694399253689E-01_wp, &
      & 5.7694399253689E-01_wp, 4.9665685187358E-17_wp,-5.6301475954244E-01_wp, &
      &-5.6301475954244E-01_wp, 0.0000000000000E+00_wp, 5.6301475954244E-01_wp, &
      & 9.9999999988150E-01_wp, 1.4115853103291E-01_wp, 1.4115853103291E-01_wp, &
      & 1.4115853103291E-01_wp, 8.3123399075353E-01_wp,-5.7694399253689E-01_wp, &
      &-5.7694399253689E-01_wp,-5.7694399253689E-01_wp, 4.9665685187358E-17_wp, &
      & 5.6301475954244E-01_wp, 5.6301475954244E-01_wp, 0.0000000000000E+00_wp, &
      & 5.6301475954244E-01_wp, 1.4115853103291E-01_wp, 9.9999999988150E-01_wp, &
      & 1.4115853103291E-01_wp, 1.4115853103291E-01_wp, 8.3123399075353E-01_wp, &
      &-5.7694399253689E-01_wp, 5.7694399253689E-01_wp, 5.7694399253689E-01_wp, &
      & 4.9665685187358E-17_wp, 5.6301475954244E-01_wp,-5.6301475954244E-01_wp, &
      & 0.0000000000000E+00_wp,-5.6301475954244E-01_wp, 1.4115853103291E-01_wp, &
      & 1.4115853103291E-01_wp, 9.9999999988150E-01_wp, 1.4115853103291E-01_wp, &
      & 8.3123399075353E-01_wp, 5.7694399253689E-01_wp, 5.7694399253689E-01_wp, &
      &-5.7694399253689E-01_wp, 4.9665685187358E-17_wp,-5.6301475954244E-01_wp, &
      & 5.6301475954244E-01_wp, 0.0000000000000E+00_wp,-5.6301475954244E-01_wp, &
      & 1.4115853103291E-01_wp, 1.4115853103291E-01_wp, 1.4115853103291E-01_wp, &
      & 9.9999999988150E-01_wp],shape(overlap_diat))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "SiH4")
   call test_overlap_diat_mol(error, mol, overlap_diat)

end subroutine test_overlap_diat_sih4


  subroutine test_q_mb01(error)
     !> Error handling
     type(error_type), allocatable, intent(out) :: error

     type(structure_type) :: mol
     real(wp), parameter :: charges(16) = reshape([ & ! calculated with GP3 standalone 
        & 0.5041712306_wp, -0.0768741000_wp, -0.4935157669_wp, &
        &-0.0831876027_wp, -0.2122917586_wp,  0.1274119295_wp, &
        &-0.0434563264_wp, -0.3788163344_wp, -0.3016588466_wp, &
        & 0.1576514268_wp,  0.1353213766_wp,  0.0150687156_wp, &
        & 0.0511522155_wp,  0.1399127014_wp, -0.0749090701_wp, & 
        & 0.5340202097_wp], shape(charges))

     call get_structure(mol, "MB16-43", "01")
     call test_q_gen(error, mol, charges)
     
  end subroutine test_q_mb01

  subroutine test_q_mb02(error)
     !> Error handling
     type(error_type), allocatable, intent(out) :: error

     type(structure_type) :: mol
     real(wp), parameter :: charges(16) = reshape([ &  
        &-8.2903701616729E-2_wp, -1.4466316351265E-1_wp, -1.1895317453167E-1_wp, &
        &-2.1055463325104E-1_wp,  4.6672727684789E-1_wp,  1.8631662619422E-1_wp, &
        &-9.6524530569012E-2_wp, -5.5370442753631E-2_wp,  2.7239683511750E-1_wp, & 
        & 1.8069642306914E-1_wp,  7.6749322633498E-3_wp,  3.6631643797735E-1_wp, &
        &-3.8041142618418E-1_wp, -9.0237343862766E-2_wp,  3.2741009140768E-2_wp, &
        &-3.3325112432852E-1_wp], shape(charges))

     call get_structure(mol, "MB16-43", "02")
     call test_q_gen(error, mol, charges)
     
  end subroutine test_q_mb02

  subroutine test_q_mb03(error)
     !> Error handling
     type(error_type), allocatable, intent(out) :: error

     type(structure_type) :: mol
     real(wp), parameter :: charges(16) = reshape([ & 
        & 1.0164987232184E-1_wp, -5.2353191449427E-1_wp,  7.8468915936174E-3_wp, &
        & 3.3621716693231E-1_wp,  5.3343641670527E-1_wp,  3.1220707646058E-2_wp, &
        &-2.5184973957319E-1_wp, -8.8831770719766E-4_wp,  1.3902266486746E-2_wp, &
        & 2.8751747800250E-3_wp, -3.9246331080380E-1_wp, -2.2054378319140E-1_wp, &
        &-3.1296934118473E-1_wp,  1.0756899130705E-1_wp,  5.3832038987584E-1_wp, &
        & 2.9208529305834E-2_wp], shape(charges))

     call get_structure(mol, "MB16-43", "03")
     call test_q_gen(error, mol, charges)
     
  end subroutine test_q_mb03

  subroutine test_q_mb04(error)
     !> Error handling
     type(error_type), allocatable, intent(out) :: error

     type(structure_type) :: mol
     real(wp), parameter :: charges(16) = reshape([ & 
        & 7.6319227329846E-4_wp, -6.3602035431309E-2_wp, -6.5671587004899E-2_wp, &
        &-1.8864318926079E-1_wp,  2.8955146220578E-1_wp,  2.6724125028655E-3_wp, & 
        & 2.9922719357434E-2_wp, -1.5947752989537E-1_wp, -5.3993174966949E-2_wp, &
        &-5.1908771187738E-2_wp, -1.5497804319660E-1_wp,  4.8902197420826E-1_wp, & 
        & 1.0421863678148E-1_wp, -3.2971896505580E-1_wp,  3.8648437944251E-2_wp, & 
        & 1.1319446072610E-1_wp], shape(charges))

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
      real(wp), parameter :: ref(3) = reshape([ & ! calculated with GP3 standalone
         0.584361099036660_wp, &
         -1.47304239280996_wp, &
         -2.25861915370679_wp], shape(ref))

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

   subroutine test_d_mb04(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
      class(container_type), allocatable :: cont
      real(wp) :: energy, efield(3), dipole(3), tmp(3)
      real(wp), parameter :: accuracy = 1e-8_wp
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

   end subroutine test_d_mb04

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
