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

module test_solvation_surface
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mctc_io_convert, only : aatoau, kcaltoau
   use mstore, only : get_structure
   use tblite_container, only : container_cache
   use tblite_mesh_lebedev, only : grid_size
   use tblite_scf_potential, only : potential_type
   use tblite_solvation_data, only : get_vdw_rad_bondi, get_vdw_rad_cosmo, get_vdw_rad_d3
   use tblite_solvation_surface, only : new_surface_integrator, surface_integrator
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: collect_solvation_surface

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains


!> Collect all exported unit tests
subroutine collect_solvation_surface(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("surface-1", test_mb01), &
      new_unittest("surface-2", test_mb02), &
      new_unittest("surface-3", test_mb03), &
      new_unittest("surface-4", test_mb04), &
      new_unittest("sigma-surface-1", test_s_mb01), &
      new_unittest("sigma-surface-2", test_s_mb02), &
      new_unittest("sigma-surface-3", test_s_mb03), &
      new_unittest("sigma-surface-asymmetry", test_s_asymmetry) &
      ]

end subroutine collect_solvation_surface


subroutine test_numg(error, sasa, mol)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Surface integrator
   type(surface_integrator), intent(inout) :: sasa
   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   integer :: iat, ic
   real(wp), allocatable :: surface(:), sr(:), sl(:), dsdr(:, :, :), dsdL(:, :, :)
   real(wp), allocatable :: numg(:, :, :)
   real(wp), parameter :: step = 1.0e-5_wp

   allocate(surface(mol%nat), sr(mol%nat), sl(mol%nat))
   allocate(dsdr(3, mol%nat, mol%nat), numg(3, mol%nat, mol%nat))
   allocate(dsdL(3, 3, mol%nat))

   call sasa%get_surface(mol, surface, dsdr, dsdL)

   do iat = 1, mol%nat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call sasa%get_surface(mol, sr)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call sasa%get_surface(mol, sl)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step

         numg(ic, iat, :) = 0.5_wp * (sr - sl) / step
      end do
   end do

   if (any(abs(numg - dsdr) > thr2)) then
      call test_failed(error, "Surface derivative does not much finite difference solution")
   end if
end subroutine test_numg


subroutine test_nums(error, sasa, mol)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Surface integrator
   type(surface_integrator), intent(inout) :: sasa
   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   integer :: ic, jc
   real(wp), allocatable :: surface(:), sr(:), sl(:), dsdr(:, :, :), dsdL(:, :, :)
   real(wp), allocatable :: nums(:, :, :), xyz(:, :)
   real(wp) :: eps(3, 3)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(surface(mol%nat), sr(mol%nat), sl(mol%nat))
   allocate(dsdr(3, mol%nat, mol%nat), dsdL(3, 3, mol%nat))
   allocate(nums(3, 3, mol%nat), xyz(3, mol%nat))

   call sasa%get_surface(mol, surface, dsdr, dsdL)

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   do ic = 1, 3
      do jc = 1, 3
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         call sasa%get_surface(mol, sr)

         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         call sasa%get_surface(mol, sl)

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz

         ! the angular grid is not rotationally invariant, therefore the strain
         ! derivative is not symmetric
         nums(jc, ic, :) = 0.5_wp * (sr - sl) / step
      end do
   end do

   ! the smoothed surface is only once continuously differentiable and the
   ! atomic surfaces are large, which limits the finite difference accuracy
   if (any(abs(nums - dsdL) > 10*thr2)) then
      call test_failed(error, "Surface strain derivative does not match")
   end if
end subroutine test_nums


subroutine test_mb01(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(surface_integrator) :: sasa
   real(wp), allocatable :: rad(:), surface(:), dsdr(:, :, :), dsdL(:, :, :)
   real(wp), parameter :: probe = 1.4_wp * aatoau
   integer, parameter :: nang = 110
   real(wp), parameter :: ref(16) = [&
      & 1.98249964603498E+2_wp, 9.34967918541344E+1_wp, 7.26746425976157E+1_wp, &
      & 3.72308705072405E+1_wp, 1.00057039498616E+2_wp, 8.72703799995796E+1_wp, &
      & 1.75563553107864E+1_wp, 5.79324044295481E+1_wp, 9.81701754804677E-3_wp, &
      & 1.05256238904348E+2_wp, 6.62363240313345E+1_wp, 1.44944528018566E+2_wp, &
      & 3.33346853562456E+1_wp, 5.79746582175529E+1_wp, 6.69252984752073E+0_wp, &
      & 4.86484694486336E+1_wp]

   call get_structure(mol, "MB16-43", "01")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat), dsdL(3, 3, mol%nat))
   rad = get_vdw_rad_d3(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)
   call sasa%get_surface(mol, surface, dsdr, dsdL)

   if (any(abs(surface - ref) > thr2)) then
      call test_failed(error, "Surface area values do not match")
      print "(es20.14e1)", surface
      return
   end if

   call test_numg(error, sasa, mol)

end subroutine test_mb01


subroutine test_mb02(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(surface_integrator) :: sasa
   real(wp), allocatable :: rad(:), surface(:), dsdr(:, :, :), dsdL(:, :, :)
   real(wp), parameter :: probe = 1.2_wp * aatoau
   integer, parameter :: nang = 230
   real(wp), parameter :: ref(16) = [&
      & 2.86084867868854E+1_wp, 7.50937555534059E+1_wp, 8.05879869880977E+1_wp, &
      & 8.24020440962820E+1_wp, 6.48136730299052E+1_wp, 1.97586791688521E+1_wp, &
      & 4.90632288004349E+1_wp, 5.29220735596789E+1_wp, 9.14599031786151E+1_wp, &
      & 1.38294851260743E+1_wp, 9.02032751808618E+1_wp, 1.13713659875286E+2_wp, &
      & 9.83820274680035E+1_wp, 5.95926059359978E+1_wp, 2.96614646358023E+0_wp, &
      & 1.44874751490690E+2_wp]

   call get_structure(mol, "MB16-43", "02")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat), dsdL(3, 3, mol%nat))
   rad = get_vdw_rad_bondi(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)
   call sasa%get_surface(mol, surface, dsdr, dsdL)

   if (any(abs(surface - ref) > thr2)) then
      call test_failed(error, "Surface area values do not match")
      print "(es20.14e1)", surface
      return
   end if

   call test_numg(error, sasa, mol)

end subroutine test_mb02


subroutine test_mb03(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(surface_integrator) :: sasa
   real(wp), allocatable :: rad(:), surface(:), dsdr(:, :, :), dsdL(:, :, :)
   real(wp), parameter :: probe = 0.2_wp * aatoau
   integer, parameter :: nang = 111
   real(wp), parameter :: ref(16) = [&
      & 4.93447390726497E+1_wp, 5.42387849176901E+1_wp, 2.58043997374119E+1_wp, &
      & 3.26892803192176E+1_wp, 1.27988010759842E+1_wp, 9.45810634518707E+1_wp, &
      & 3.43532470377123E+1_wp, 2.76341416140764E+1_wp, 2.74903764017798E+1_wp, &
      & 2.85813017859723E+1_wp, 7.99313005786035E+1_wp, 1.26258175473983E+2_wp, &
      & 5.38016574162998E+1_wp, 4.16287245622076E+1_wp, 9.95930646536509E+1_wp, &
      & 2.36024718294637E+1_wp]

   call get_structure(mol, "MB16-43", "03")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat), dsdL(3, 3, mol%nat))
   rad = get_vdw_rad_cosmo(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)
   call sasa%get_surface(mol, surface, dsdr, dsdL)

   if (any(abs(surface - ref) > thr2)) then
      call test_failed(error, "Surface area values do not match")
      print "(es20.14e1)", surface
      return
   end if

   call test_numg(error, sasa, mol)

end subroutine test_mb03


subroutine test_mb04(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(surface_integrator) :: sasa
   integer :: iang
   real(wp) :: this_thr
   real(wp), allocatable :: rad(:), surface(:), dsdr(:, :, :), dsdL(:, :, :)
   real(wp), parameter :: probe = 1.4_wp * aatoau
   real(wp), parameter :: ref(16) = [&
      & 3.17888247611432E+1_wp, 1.10843192696983E+2_wp, 6.88322052202638E+1_wp, &
      & 1.14544539499935E+2_wp, 1.70720777217273E+2_wp, 3.13678106939536E+1_wp, &
      & 4.58475696363698E+1_wp, 1.93179973492600E+2_wp, 6.00038960472752E+1_wp, &
      & 6.11241830969292E+1_wp, 4.51433358912678E+1_wp, 9.79240755827738E+0_wp, &
      & 1.11790314288253E+2_wp, 3.26024198194955E+1_wp, 7.04914426556603E+1_wp, &
      & 7.70033482758105E+1_wp]

   call get_structure(mol, "MB16-43", "04")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat), dsdL(3, 3, mol%nat))
   rad = get_vdw_rad_bondi(mol%num)

   do iang = 1, size(grid_size)
      call new_surface_integrator(sasa, mol%id, rad, probe, grid_size(iang))
      call sasa%get_surface(mol, surface, dsdr, dsdL)

      if (grid_size(iang) > 1000) then
         this_thr = 0.1_wp
      else if (grid_size(iang) > 250) then
         this_thr = 1.0_wp
      else if (grid_size(iang) > 75) then
         this_thr = 10.0_wp
      else
         cycle
      end if
      call check(error, 0.0_wp, norm2(abs(surface - ref)), thr=this_thr)
      if (allocated(error)) return
   end do

   if (any(abs(surface - ref) > thr2)) then
      call test_failed(error, "Surface area values do not match")
      print "(es20.14e1)", surface
      return
   end if

end subroutine test_mb04


subroutine test_s_mb01(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(surface_integrator) :: sasa
   real(wp), allocatable :: rad(:)
   real(wp), parameter :: probe = 1.4_wp * aatoau
   integer, parameter :: nang = 110

   call get_structure(mol, "MB16-43", "01")
   rad = get_vdw_rad_d3(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)

   call test_nums(error, sasa, mol)

end subroutine test_s_mb01


subroutine test_s_mb02(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(surface_integrator) :: sasa
   real(wp), allocatable :: rad(:)
   real(wp), parameter :: probe = 1.2_wp * aatoau
   integer, parameter :: nang = 230

   call get_structure(mol, "MB16-43", "02")
   rad = get_vdw_rad_bondi(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)

   call test_nums(error, sasa, mol)

end subroutine test_s_mb02


subroutine test_s_mb03(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(surface_integrator) :: sasa
   real(wp), allocatable :: rad(:)
   real(wp), parameter :: probe = 0.2_wp * aatoau
   integer, parameter :: nang = 111

   call get_structure(mol, "MB16-43", "03")
   rad = get_vdw_rad_cosmo(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)

   call test_nums(error, sasa, mol)

end subroutine test_s_mb03


subroutine test_s_asymmetry(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(surface_integrator) :: sasa
   real(wp), allocatable :: rad(:), surface(:), dsdr(:, :, :), dsdL(:, :, :)
   real(wp) :: sigma(3, 3), asym_small, asym_large
   real(wp), parameter :: probe = 1.4_wp * aatoau

   call get_structure(mol, "MB16-43", "04")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat), dsdL(3, 3, mol%nat))
   rad = get_vdw_rad_bondi(mol%num)

   ! the antisymmetric strain is the torque of the fixed angular grid
   call new_surface_integrator(sasa, mol%id, rad, probe, 110)
   call sasa%get_surface(mol, surface, dsdr, dsdL)
   sigma(:, :) = sum(dsdL, 3)
   asym_small = norm2(0.5_wp * (sigma - transpose(sigma))) / norm2(sigma)

   call new_surface_integrator(sasa, mol%id, rad, probe, 5810)
   call sasa%get_surface(mol, surface, dsdr, dsdL)
   sigma(:, :) = sum(dsdL, 3)
   asym_large = norm2(0.5_wp * (sigma - transpose(sigma))) / norm2(sigma)

   if (asym_small < 1.0e-2_wp) then
      call test_failed(error, &
         & "Surface strain derivative is not asymmetric for small grid")
      print "(es20.13)", asym_small
      return
   end if

   if (asym_large > 1.0e-3_wp) then
      call test_failed(error, "Surface strain asymmetry does not decrease for large grid")
      print "(es20.13)", asym_large
      return
   end if

end subroutine test_s_asymmetry

end module test_solvation_surface
