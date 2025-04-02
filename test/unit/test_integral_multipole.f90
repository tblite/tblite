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

module test_integral_multipole
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_basis_type
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_cutoff, only : get_lattice_points
   use tblite_integral_dipole
   use tblite_integral_multipole
   implicit none
   private

   public :: collect_integral_multipole

   real(wp), parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains


!> Collect all exported unit tests
subroutine collect_integral_multipole(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("overlap-dipole-diat-alh3", test_overlap_dipole_diat_alh3), &
      new_unittest("overlap-multipole-diat-alh3", test_overlap_multipole_diat_alh3), &
      new_unittest("dipole-trans-ss", test_dipole_ss), &
      new_unittest("dipole-trans-pp", test_dipole_pp), &
      new_unittest("dipole-trans-dd", test_dipole_dd), &
      new_unittest("dipole-grad-ss", test_dipole_grad_ss) &
      ]

end subroutine collect_integral_multipole


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


subroutine test_dipole_ss(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6
   integer :: stat, i
   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3), r2
   real(wp) :: overlap(1, 1), dipolei(3, 1, 1), dipolej(3, 1, 1)

   call slater_to_gauss(ng, 2, 0, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 1, 0, 1.0_wp, cgtoj, .true., stat)

   call random_number(vec)
   vec = vec - 0.5_wp
   r2 = sum(vec**2)

   call dipole_cgto(cgtoi, cgtoj, r2, vec, 100.0_wp, overlap, dipolej)

   vec(:) = -vec

   call dipole_cgto(cgtoj, cgtoi, r2, vec, 100.0_wp, overlap, dipolei)

   do i = 1, 3
      call check(error, dipolei(i, 1, 1) + vec(i) * overlap(1, 1), dipolej(i, 1, 1), thr=thr)
      if (allocated(error)) return
   end do

end subroutine test_dipole_ss


subroutine test_dipole_pp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6, n = 2

   integer :: stat, i, j, k
   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3), r2
   real(wp) :: overlap(3, 3), dipolei(3, 3, 3), dipolej(3, 3, 3)

   call slater_to_gauss(ng, 2, 1, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 3, 1, 1.0_wp, cgtoj, .true., stat)

   call random_number(vec)
   vec = vec - 0.5_wp
   r2 = sum(vec**2)

   call dipole_cgto(cgtoi, cgtoj, r2, vec, 100.0_wp, overlap, dipolej)

   vec(:) = -vec

   call dipole_cgto(cgtoj, cgtoi, r2, vec, 100.0_wp, overlap, dipolei)

   do i = 1, 3
      do j = 1, 3
         do k = 1, 3
            call check(error, dipolei(k, j, i) + vec(k) * overlap(j, i), dipolej(k, i, j), thr=thr)
            if (allocated(error)) return
         end do
      end do
   end do

end subroutine test_dipole_pp


subroutine test_dipole_dd(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6, n = 2

   integer :: stat, i, j, k
   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3), r2
   real(wp) :: overlap(5, 5), dipolei(3, 5, 5), dipolej(3, 5, 5)

   call slater_to_gauss(ng, 3, 2, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 4, 2, 1.0_wp, cgtoj, .true., stat)

   call random_number(vec)
   vec = vec - 0.5_wp
   r2 = sum(vec**2)

   call dipole_cgto(cgtoi, cgtoj, r2, vec, 100.0_wp, overlap, dipolej)

   vec(:) = -vec

   call dipole_cgto(cgtoj, cgtoi, r2, vec, 100.0_wp, overlap, dipolei)

   do i = 1, 5
      do j = 1, 5
         do k = 1, 3
            call check(error, dipolei(k, j, i) + vec(k) * overlap(j, i), dipolej(k, i, j), thr=thr)
            if (allocated(error)) return
         end do
      end do
   end do

end subroutine test_dipole_dd


subroutine test_dipole_grad_ss(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6
   integer :: stat, i, j
   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3), r2, zero(3)
   real(wp) :: overlap(1, 1), doverlapi(3, 1, 1)
   real(wp) :: dipole(3, 1, 1), ddipolei(3, 3, 1, 1), ddipolej(3, 3, 1, 1)
   real(wp) :: quadrupole(6, 1, 1), dquadrupolei(3, 6, 1, 1), dquadrupolej(3, 6, 1, 1)
   real(wp) :: sr(1, 1), sl(1, 1), dr(3, 1, 1), dl(3, 1, 1), qr(6, 1, 1), ql(6, 1, 1)
   real(wp), parameter :: step = 1.0e-6_wp

   call slater_to_gauss(ng, 2, 0, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 1, 0, 1.0_wp, cgtoj, .true., stat)

   call random_number(vec)
   vec = vec - 0.5_wp
   zero = 0
   r2 = sum(vec**2)

   call multipole_grad_cgto(cgtoj, cgtoi, r2, vec, 100.0_wp, overlap, dipole, quadrupole, &
      & doverlapi, ddipolej, dquadrupolej, ddipolei, dquadrupolei)

   do i = 1, 3
      vec(i) = vec(i) + step
      r2 = sum(vec**2)
      call multipole_cgto(cgtoj, cgtoi, r2, vec, 100.0_wp, sr, dr, qr)
      vec(i) = vec(i) - 2*step
      r2 = sum(vec**2)
      call multipole_cgto(cgtoj, cgtoi, r2, vec, 100.0_wp, sl, dl, ql)
      vec(i) = vec(i) + step
      ddipolej(i, :, :, :) = 0.5_wp * (dr - dl) / step
      dquadrupolej(i, :, :, :) = 0.5_wp * (qr - ql) / step
   end do

   num: do i = 1, 3
      do j = 1, 3
         call check(error, ddipolei(i, j, 1, 1), ddipolej(i, j, 1, 1), thr=thr)
         if (allocated(error)) exit num
      end do
      do j = 1, 6
         call check(error, dquadrupolei(i, j, 1, 1), dquadrupolej(i, j, 1, 1), thr=thr)
         if (allocated(error)) exit num
      end do
   end do num
   if (allocated(error)) return

end subroutine test_dipole_grad_ss

subroutine test_overlap_dipole_diat_mol(error, mol, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: ref(:, :)

   type(basis_type) :: bas
   real(wp), allocatable :: lattr(:, :), overlap(:, :), overlap_diat(:, :)
   real(wp), allocatable :: dipole(:, :, :)
   real(wp) :: cutoff
   integer :: ii, jj
   real(wp) :: scalfac(3,86)

   scalfac = 1.2_wp

   call make_basis(bas, mol, 6)
   call check(error, bas%nao, size(ref, 1))
   if (allocated(error)) return

   cutoff = get_cutoff(bas)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

   allocate(overlap(bas%nao, bas%nao), overlap_diat(bas%nao, bas%nao))
   allocate(dipole(3, bas%nao, bas%nao))
   call get_dipole_integrals(mol, lattr, cutoff, bas, scalfac, overlap, overlap_diat, dipole)

   do ii = 1, size(overlap_diat, 2)
      do jj = 1, size(overlap_diat, 1)
         call check(error, overlap_diat(jj, ii), ref(jj, ii), thr=thr)
         if (allocated(error)) return
      end do
   end do


end subroutine test_overlap_dipole_diat_mol

subroutine test_overlap_dipole_diat_alh3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 12
   real(wp), parameter :: overlap_diat(nao, nao) = reshape([&
      & 9.99999999869333E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 4.83197454971925E-1_wp, 4.83197454971925E-1_wp, 4.83197454971926E-1_wp,&
      & 0.00000000000000E+0_wp, 9.99999999998060E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 5.19330455705926E-1_wp,-5.19330455705926E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999998060E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999998060E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.99835578400188E-1_wp,-2.99835578400188E-1_wp, 5.99671156800377E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830206E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-4.08220347867244E-1_wp, 4.08220347867244E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830206E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999830206E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.72146898578163E-1_wp,-2.72146898578163E-1_wp,-2.72146898578163E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830206E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830206E-1_wp,&
      &-2.35686127729836E-1_wp,-2.35686127729836E-1_wp, 4.71372255459672E-1_wp,&
      & 4.83197454971925E-1_wp, 5.19330455705926E-1_wp, 0.00000000000000E+0_wp,&
      &-2.99835578400188E-1_wp,-4.08220347867244E-1_wp, 0.00000000000000E+0_wp,&
      &-2.72146898578163E-1_wp, 0.00000000000000E+0_wp,-2.35686127729836E-1_wp,&
      & 9.99999999881495E-1_wp, 4.25520423996964E-2_wp, 4.25520423996966E-2_wp,&
      & 4.83197454971925E-1_wp,-5.19330455705926E-1_wp, 0.00000000000000E+0_wp,&
      &-2.99835578400188E-1_wp, 4.08220347867244E-1_wp, 0.00000000000000E+0_wp,&
      &-2.72146898578163E-1_wp, 0.00000000000000E+0_wp,-2.35686127729836E-1_wp,&
      & 4.25520423996964E-2_wp, 9.99999999881495E-1_wp, 4.25520423996966E-2_wp,&
      & 4.83197454971927E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 5.99671156800377E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.72146898578163E-1_wp, 0.00000000000000E+0_wp, 4.71372255459672E-1_wp,&
      & 4.25520423996966E-2_wp, 4.25520423996966E-2_wp, 9.99999999881495E-1_wp],&
      & shape(overlap_diat))
   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "AlH3")
   call test_overlap_dipole_diat_mol(error, mol, overlap_diat)

end subroutine test_overlap_dipole_diat_alh3

subroutine test_overlap_multipole_diat_mol(error, mol, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: ref(:, :)

   type(basis_type) :: bas
   real(wp), allocatable :: lattr(:, :), overlap(:, :), overlap_diat(:, :)
   real(wp), allocatable :: dipole(:, :, :), quadrupole(:, :, :)
   real(wp) :: cutoff
   integer :: ii, jj
   real(wp) :: scalfac(3,86)

   scalfac = 1.2_wp

   call make_basis(bas, mol, 6)
   call check(error, bas%nao, size(ref, 1))
   if (allocated(error)) return

   cutoff = get_cutoff(bas)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

   allocate(overlap(bas%nao, bas%nao), overlap_diat(bas%nao, bas%nao))
   allocate(dipole(3, bas%nao, bas%nao), quadrupole(6, bas%nao, bas%nao))
   call get_multipole_integrals(mol, lattr, cutoff, bas, scalfac, overlap, overlap_diat, &
      & dipole, quadrupole)

   do ii = 1, size(overlap_diat, 2)
      do jj = 1, size(overlap_diat, 1)
         call check(error, overlap_diat(jj, ii), ref(jj, ii), thr=thr)
         if (allocated(error)) return
      end do
   end do


end subroutine test_overlap_multipole_diat_mol

subroutine test_overlap_multipole_diat_alh3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 12
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 9.99999999869333E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 4.83197454971925E-1_wp, 4.83197454971925E-1_wp, 4.83197454971926E-1_wp,&
      & 0.00000000000000E+0_wp, 9.99999999998060E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 5.19330455705926E-1_wp,-5.19330455705926E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999998060E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999998060E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.99835578400188E-1_wp,-2.99835578400188E-1_wp, 5.99671156800377E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830206E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-4.08220347867244E-1_wp, 4.08220347867244E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830206E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999830206E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.72146898578163E-1_wp,-2.72146898578163E-1_wp,-2.72146898578163E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830206E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830206E-1_wp,&
      &-2.35686127729836E-1_wp,-2.35686127729836E-1_wp, 4.71372255459672E-1_wp,&
      & 4.83197454971925E-1_wp, 5.19330455705926E-1_wp, 0.00000000000000E+0_wp,&
      &-2.99835578400188E-1_wp,-4.08220347867244E-1_wp, 0.00000000000000E+0_wp,&
      &-2.72146898578163E-1_wp, 0.00000000000000E+0_wp,-2.35686127729836E-1_wp,&
      & 9.99999999881495E-1_wp, 4.25520423996964E-2_wp, 4.25520423996966E-2_wp,&
      & 4.83197454971925E-1_wp,-5.19330455705926E-1_wp, 0.00000000000000E+0_wp,&
      &-2.99835578400188E-1_wp, 4.08220347867244E-1_wp, 0.00000000000000E+0_wp,&
      &-2.72146898578163E-1_wp, 0.00000000000000E+0_wp,-2.35686127729836E-1_wp,&
      & 4.25520423996964E-2_wp, 9.99999999881495E-1_wp, 4.25520423996966E-2_wp,&
      & 4.83197454971927E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 5.99671156800377E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.72146898578163E-1_wp, 0.00000000000000E+0_wp, 4.71372255459672E-1_wp,&
      & 4.25520423996966E-2_wp, 4.25520423996966E-2_wp, 9.99999999881495E-1_wp],&
      & shape(overlap))
   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "AlH3")
   call test_overlap_multipole_diat_mol(error, mol, overlap)

end subroutine test_overlap_multipole_diat_alh3

end module test_integral_multipole
