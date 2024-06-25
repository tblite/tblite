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

module test_hamiltonian
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_basis_type
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_cutoff, only : get_lattice_points
   use tblite_data_covrad, only : get_covalent_rad
   use tblite_lapack_sygvd, only : sygvd_solver
   use tblite_integral_overlap
   use tblite_ncoord_gfn, only : gfn_ncoord_type, new_gfn_ncoord
   use tblite_ncoord_type, only : get_coordination_number
   use tblite_xtb_gfn2
   use tblite_xtb_h0
   implicit none
   private

   public :: collect_hamiltonian

   real(wp), parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = 1.0e2_wp*sqrt(epsilon(1.0_wp))

contains


!> Collect all exported unit tests
subroutine collect_hamiltonian(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("hamiltonian-h2", test_hamiltonian_h2), &
      new_unittest("hamiltonian-lih", test_hamiltonian_lih), &
      new_unittest("hamiltonian-s2", test_hamiltonian_s2), &
      new_unittest("hamiltonian-sih4", test_hamiltonian_sih4) &
      ]

end subroutine collect_hamiltonian


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


subroutine test_hamiltonian_mol(error, mol, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: ref(:, :)

   type(basis_type) :: bas
   type(tb_hamiltonian) :: h0
   type(gfn_ncoord_type) :: gfn_ncoord
   type(adjacency_list) :: list
   real(wp), parameter :: cn_cutoff = 30.0_wp
   real(wp), allocatable :: lattr(:, :), cn(:), rcov(:)
   real(wp), allocatable :: overlap(:, :), hamiltonian(:, :), selfenergy(:)
   real(wp), allocatable :: dpint(:, :, :), qpint(:, :, :)
   real(wp) :: cutoff
   integer :: ii, jj

   call make_basis(bas, mol, 6)
   !print*, bas%nao
   call check(error, bas%nao, size(ref, 1))
   if (allocated(error)) return

   call new_hamiltonian(h0, mol, bas, gfn2_h0spec(mol))

   allocate(cn(mol%nat), rcov(mol%nid))
   rcov(:) = get_covalent_rad(mol%num)
   call new_gfn_ncoord(gfn_ncoord, mol, cn_cutoff, rcov)
   call get_lattice_points(mol%periodic, mol%lattice, cn_cutoff, lattr)
   call get_coordination_number(gfn_ncoord, mol, lattr, cn_cutoff, cn)

   cutoff = get_cutoff(bas)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call new_adjacency_list(list, mol, lattr, cutoff)

   allocate(selfenergy(bas%nsh))
   call get_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, cn=cn, selfenergy=selfenergy)

   allocate(overlap(bas%nao, bas%nao), dpint(3, bas%nao, bas%nao), &
      & qpint(6, bas%nao, bas%nao), hamiltonian(bas%nao, bas%nao))
   call get_hamiltonian(mol, lattr, list, bas, h0, selfenergy, overlap, dpint, qpint, &
      & hamiltonian)

   !where(abs(hamiltonian) < thr) hamiltonian = 0.0_wp
   !print '(*(6x,"&", 3(es20.14e1, "_wp":, ","), "&", /))', hamiltonian

   do ii = 1, size(hamiltonian, 2)
      do jj = 1, size(hamiltonian, 1)
         call check(error, hamiltonian(jj, ii), ref(jj, ii), thr=thr2)
         if (allocated(error)) then
            print '(2es20.13)', hamiltonian(jj, ii), ref(jj, ii), &
               & hamiltonian(jj, ii) - ref(jj, ii)
            return
         end if
      end do
   end do

   !allocate(eigval(bas%nao))
   !call sygvd%solve(hamiltonian, overlap, eigval, error)
   !if (allocated(error)) return

   !print '(*("&", 3(es20.14e1, "_wp":, ","), "&", /))', eigval

end subroutine test_hamiltonian_mol

subroutine test_hamiltonian_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 2
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      &-3.91986875628795E-1_wp,-4.69784163992013E-1_wp,-4.69784163992013E-1_wp,&
      &-3.91986875628795E-1_wp],shape(hamiltonian))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "H2")
   call test_hamiltonian_mol(error, mol, hamiltonian)

end subroutine test_hamiltonian_h2

subroutine test_hamiltonian_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 5
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      &-1.85652592301171E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-2.04060202125475E-1_wp, 0.00000000000000E+0_wp,&
      &-7.93540995798542E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-7.93540995798542E-2_wp, 0.00000000000000E+0_wp,-2.64332069820779E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-7.93540995798542E-2_wp, 0.00000000000000E+0_wp,-2.04060202125475E-1_wp,&
      & 0.00000000000000E+0_wp,-2.64332069820779E-1_wp, 0.00000000000000E+0_wp,&
      &-3.91761150560104E-1_wp],shape(hamiltonian))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_hamiltonian_mol(error, mol, hamiltonian)

end subroutine test_hamiltonian_lih

subroutine test_hamiltonian_s2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 18
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      &-7.35145168755863E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.92782989368414E-1_wp, 0.00000000000000E+0_wp, 2.36427133955232E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.05951879709695E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-4.17765769244149E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-9.33756600277693E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.20176383632229E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-4.17765769244149E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.36427133955232E-1_wp, 0.00000000000000E+0_wp, 2.58607490728104E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.23733826451155E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-4.17765769244149E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-9.33756600277693E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 1.20176383632229E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-2.27200179621343E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.00344447885664E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.27200179621343E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.20176383632229E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.45142819682488E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.27200179621343E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.05951879709695E-1_wp, 0.00000000000000E+0_wp, 1.23733826451155E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-9.40352474081307E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-2.27200179621343E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.20176383632229E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 2.45142819682488E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.27200179621343E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.00344447885664E-2_wp,&
      &-1.92782989368414E-1_wp, 0.00000000000000E+0_wp,-2.36427133955232E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.05951879709695E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-7.35145168755863E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-9.33756600277693E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.20176383632229E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-4.17765769244149E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 2.36427133955232E-1_wp, 0.00000000000000E+0_wp, 2.58607490728104E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.23733826451155E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-4.17765769244149E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-9.33756600277693E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.20176383632229E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-4.17765769244149E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.00344447885664E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-2.27200179621343E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 1.20176383632229E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.45142819682488E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.27200179621343E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.05951879709695E-1_wp, 0.00000000000000E+0_wp,-1.23733826451155E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-9.40352474081307E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.27200179621343E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.20176383632229E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 2.45142819682488E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-2.27200179621343E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.00344447885664E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.27200179621343E-2_wp],&
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
      &-5.52421014706411E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-3.36004320989251E-1_wp,-3.36004320989251E-1_wp,-3.36004320989251E-1_wp,&
      &-3.36004320989251E-1_wp, 0.00000000000000E+0_wp,-2.35769691508991E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.53874698572085E-1_wp, 1.53874698572085E-1_wp,&
      & 1.53874698572085E-1_wp,-1.53874698572085E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-2.35769691508991E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.53874698572085E-1_wp,&
      & 1.53874698572085E-1_wp,-1.53874698572085E-1_wp,-1.53874698572085E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.35769691508991E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.53874698572085E-1_wp, 1.53874698572085E-1_wp,-1.53874698572085E-1_wp,&
      & 1.53874698572085E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-4.13801903077918E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.23912378597848E-1_wp,-1.23912378597848E-1_wp,&
      & 1.23912378597848E-1_wp, 1.23912378597848E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-4.13801903077918E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.23912378597848E-1_wp,&
      &-1.23912378597848E-1_wp, 1.23912378597848E-1_wp,-1.23912378597848E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-4.13801903077918E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-4.13801903077918E-2_wp,&
      & 0.00000000000000E+0_wp, 1.23912378597848E-1_wp,-1.23912378597848E-1_wp,&
      &-1.23912378597848E-1_wp, 1.23912378597848E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-4.13801903077918E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-3.36004320989251E-1_wp,-1.53874698572085E-1_wp, 1.53874698572085E-1_wp,&
      &-1.53874698572085E-1_wp,-1.23912378597848E-1_wp, 1.23912378597848E-1_wp,&
      & 0.00000000000000E+0_wp, 1.23912378597848E-1_wp, 0.00000000000000E+0_wp,&
      &-3.91823589869151E-1_wp,-4.31486730191658E-2_wp,-4.31486730191658E-2_wp,&
      &-4.31486730191658E-2_wp,-3.36004320989251E-1_wp, 1.53874698572085E-1_wp,&
      & 1.53874698572085E-1_wp, 1.53874698572085E-1_wp,-1.23912378597848E-1_wp,&
      &-1.23912378597848E-1_wp, 0.00000000000000E+0_wp,-1.23912378597848E-1_wp,&
      & 0.00000000000000E+0_wp,-4.31486730191658E-2_wp,-3.91823589869151E-1_wp,&
      &-4.31486730191658E-2_wp,-4.31486730191658E-2_wp,-3.36004320989251E-1_wp,&
      & 1.53874698572085E-1_wp,-1.53874698572085E-1_wp,-1.53874698572085E-1_wp,&
      & 1.23912378597848E-1_wp, 1.23912378597848E-1_wp, 0.00000000000000E+0_wp,&
      &-1.23912378597848E-1_wp, 0.00000000000000E+0_wp,-4.31486730191658E-2_wp,&
      &-4.31486730191658E-2_wp,-3.91823589869151E-1_wp,-4.31486730191658E-2_wp,&
      &-3.36004320989251E-1_wp,-1.53874698572085E-1_wp,-1.53874698572085E-1_wp,&
      & 1.53874698572085E-1_wp, 1.23912378597848E-1_wp,-1.23912378597848E-1_wp,&
      & 0.00000000000000E+0_wp, 1.23912378597848E-1_wp, 0.00000000000000E+0_wp,&
      &-4.31486730191658E-2_wp,-4.31486730191658E-2_wp,-4.31486730191658E-2_wp,&
      &-3.91823589869151E-1_wp],shape(hamiltonian))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "SiH4")
   call test_hamiltonian_mol(error, mol, hamiltonian)

end subroutine test_hamiltonian_sih4


end module test_hamiltonian
