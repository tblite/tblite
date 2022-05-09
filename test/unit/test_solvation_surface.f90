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
   use mctc_io_constants, only : pi
   use mctc_io_convert, only : aatoau, kcaltoau
   use mstore, only : get_structure
   use tblite_container, only : container_cache
   use tblite_mesh_lebedev, only : grid_size
   use tblite_scf_potential, only : potential_type
   use tblite_solvation_cds
   use tblite_solvation_data
   use tblite_solvation_surface
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: collect_solvation_surface

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   real(wp), parameter :: tension_water(20) = 1.0e-5_wp * [&
      &-0.08499967_wp, 0.46780225_wp,-2.87013596_wp,-3.95935069_wp,-0.29783987_wp, &
      &-0.48323273_wp, 0.00133622_wp, 0.20448945_wp, 0.20150600_wp, 0.36379863_wp, &
      &-3.47082133_wp,-0.93451053_wp,-1.46342018_wp,-0.32774697_wp,-0.38015204_wp, &
      &-0.35311116_wp,-0.19972593_wp,-0.12891363_wp,-1.19450558_wp,-1.61289300_wp]
   real(wp), parameter :: hbond_water(20) = -kcaltoau * [&
      & 6.70894947_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, &
      & 1.26459036_wp, 3.52206160_wp, 2.30440543_wp, 1.98829409_wp, 0.00000000_wp, &
      & 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, 2.68116653_wp, &
      & 0.38262428_wp, 1.02948365_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp]**2
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
      new_unittest("sasa-e", test_e_sasa), &
      new_unittest("sasa-g", test_g_sasa), &
      new_unittest("sasa-p", test_p_sasa) &
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
   real(wp), allocatable :: surface(:), sr(:), sl(:), dsdr(:, :, :), numg(:, :, :)
   real(wp), parameter :: step = 1.0e-5_wp

   allocate(surface(mol%nat), sr(mol%nat), sl(mol%nat))
   allocate(dsdr(3, mol%nat, mol%nat), numg(3, mol%nat, mol%nat))

   call sasa%get_surface(mol, surface, dsdr)

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


subroutine test_mb01(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(surface_integrator) :: sasa
   real(wp), allocatable :: rad(:), surface(:), dsdr(:, :, :)
   real(wp), parameter :: probe = 1.4_wp * aatoau
   integer, parameter :: nang = 110
   real(wp), parameter :: ref(16) = [&
      & 1.98249964115938E+2_wp, 9.34967917451750E+1_wp, 7.26746426166825E+1_wp, &
      & 3.72308704194178E+1_wp, 1.00057039382294E+2_wp, 8.72703799067120E+1_wp, &
      & 1.75563552659476E+1_wp, 5.79324044451170E+1_wp, 9.81702069377092E-3_wp, &
      & 1.05256238709304E+2_wp, 6.62363239088925E+1_wp, 1.44944527871945E+2_wp, &
      & 3.33346854609299E+1_wp, 5.79746583890200E+1_wp, 6.69252987306994E+0_wp, &
      & 4.86484695123678E+1_wp]

   call get_structure(mol, "MB16-43", "01")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat))
   rad = get_vdw_rad_d3(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)
   call sasa%get_surface(mol, surface, dsdr)

   if (any(abs(surface - ref) > thr2)) then
      call test_failed(error, "Surface area values do not match")
      print '(es20.14e1)', surface
      return
   end if

   call test_numg(error, sasa, mol)

end subroutine test_mb01


subroutine test_mb02(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(surface_integrator) :: sasa
   real(wp), allocatable :: rad(:), surface(:), dsdr(:, :, :)
   real(wp), parameter :: probe = 1.2_wp * aatoau
   integer, parameter :: nang = 230
   real(wp), parameter :: ref(16) = [&
      & 2.86084867409084E+1_wp, 7.50937554738823E+1_wp, 8.05879870633801E+1_wp, &
      & 8.24020440071839E+1_wp, 6.48136730333555E+1_wp, 1.97586790838718E+1_wp, &
      & 4.90632286612891E+1_wp, 5.29220735196209E+1_wp, 9.14599027767086E+1_wp, &
      & 1.38294851038605E+1_wp, 9.02032750404975E+1_wp, 1.13713659886359E+2_wp, &
      & 9.83820273872444E+1_wp, 5.95926059064588E+1_wp, 2.96614651477179E+0_wp, &
      & 1.44874751678823E+2_wp]

   call get_structure(mol, "MB16-43", "02")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat))
   rad = get_vdw_rad_bondi(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)
   call sasa%get_surface(mol, surface, dsdr)

   if (any(abs(surface - ref) > thr2)) then
      call test_failed(error, "Surface area values do not match")
      print '(es20.14e1)', surface
      return
   end if

   call test_numg(error, sasa, mol)

end subroutine test_mb02


subroutine test_mb03(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(surface_integrator) :: sasa
   real(wp), allocatable :: rad(:), surface(:), dsdr(:, :, :)
   real(wp), parameter :: probe = 0.2_wp * aatoau
   integer, parameter :: nang = 111
   real(wp), parameter :: ref(16) = [&
      & 4.93447390062137E+1_wp, 5.42387848008597E+1_wp, 2.58043996919982E+1_wp, &
      & 3.26892802882332E+1_wp, 1.27988011902719E+1_wp, 9.45810634328472E+1_wp, &
      & 3.43532468932011E+1_wp, 2.76341415226557E+1_wp, 2.74903763271451E+1_wp, &
      & 2.85813017380994E+1_wp, 7.99313006099443E+1_wp, 1.26258175368922E+2_wp, &
      & 5.38016574327888E+1_wp, 4.16287245551897E+1_wp, 9.95930646163274E+1_wp, &
      & 2.36024718010776E+1_wp]

   call get_structure(mol, "MB16-43", "03")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat))
   rad = get_vdw_rad_cosmo(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)
   call sasa%get_surface(mol, surface, dsdr)

   if (any(abs(surface - ref) > thr2)) then
      call test_failed(error, "Surface area values do not match")
      print '(es20.14e1)', surface
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
   real(wp), allocatable :: rad(:), surface(:), dsdr(:, :, :)
   real(wp), parameter :: probe = 1.4_wp * aatoau
   real(wp), parameter :: ref(16) = [&
      & 3.17888246831669E+1_wp, 1.10843192647448E+2_wp, 6.88322051937029E+1_wp, &
      & 1.14544539434107E+2_wp, 1.70720777045875E+2_wp, 3.13678105884819E+1_wp, &
      & 4.58475695169198E+1_wp, 1.93179973322729E+2_wp, 6.00038959489535E+1_wp, &
      & 6.11241830519643E+1_wp, 4.51433360077588E+1_wp, 9.79240756407461E+0_wp, &
      & 1.11790314128626E+2_wp, 3.26024197541769E+1_wp, 7.04914424301540E+1_wp, &
      & 7.70033481249225E+1_wp]

   call get_structure(mol, "MB16-43", "04")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat))
   rad = get_vdw_rad_bondi(mol%num)

   do iang = 1, size(grid_size)
      call new_surface_integrator(sasa, mol%id, rad, probe, grid_size(iang))
      call sasa%get_surface(mol, surface, dsdr)

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
      print '(es20.14e1)', surface
      return
   end if

end subroutine test_mb04


subroutine test_e(error, mol, input, qat, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Solvation model input
   type(cds_input), intent(in) :: input

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Reference energy
   real(wp), intent(in) :: ref

   type(cds_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp) :: energy(mol%nat)

   energy = 0.0_wp
   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1))

   solv = cds_solvation(mol, input)

   call solv%update(mol, cache)
   call solv%get_engrad(mol, cache, energy)
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)

   if (abs(sum(energy) - ref) > thr) then
      call test_failed(error, "Energy does not match reference")
      print *, sum(energy)
   end if
end subroutine test_e


subroutine test_g(error, mol, input, qat)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Solvation model input
   type(cds_input), intent(in) :: input

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   type(cds_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp), allocatable :: gradient(:, :), numg(:, :)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3)
   integer :: ii, ic

   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1))

   solv = cds_solvation(mol, input)

   allocate(numg(3, mol%nat), gradient(3, mol%nat))
   do ii = 1, mol%nat
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(ic, ii) = mol%xyz(ic, ii) + step
         call solv%update(mol, cache)
         call solv%get_engrad(mol, cache, er)
         call solv%get_potential(mol, cache, wfn, pot)
         call solv%get_energy(mol, cache, wfn, er)

         mol%xyz(ic, ii) = mol%xyz(ic, ii) - 2*step
         call solv%update(mol, cache)
         call solv%get_engrad(mol, cache, el)
         call solv%get_potential(mol, cache, wfn, pot)
         call solv%get_energy(mol, cache, wfn, el)

         mol%xyz(ic, ii) = mol%xyz(ic, ii) + step
         numg(ic, ii) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call solv%update(mol, cache)
   call solv%get_engrad(mol, cache, energy, gradient, sigma)
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)
   call solv%get_gradient(mol, cache, wfn, gradient, sigma)

   if (any(abs(gradient - numg) > thr)) then
      call test_failed(error, "Gradient does not match")
      print '(3es20.13)', gradient
      print '(a)', "---"
      print '(3es20.13)', numg
      print '(a)', "---"
      print '(3es20.13)', gradient - numg
   end if
end subroutine test_g


subroutine test_p(error, mol, input, qat)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Solvation model input
   type(cds_input), intent(in) :: input

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   type(cds_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), parameter :: thr = 1e+3_wp*sqrt(epsilon(1.0_wp))
   real(wp), allocatable :: vat(:)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat)
   integer :: ii

   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1))

   solv = cds_solvation(mol, input)

   call solv%update(mol, cache)

   allocate(vat(mol%nat))
   do ii = 1, mol%nat
      er = 0.0_wp
      el = 0.0_wp
      wfn%qat(ii, 1) = wfn%qat(ii, 1) + step
      call solv%get_potential(mol, cache, wfn, pot)
      call solv%get_energy(mol, cache, wfn, er)

      wfn%qat(ii, 1) = wfn%qat(ii, 1) - 2*step
      call solv%get_potential(mol, cache, wfn, pot)
      call solv%get_energy(mol, cache, wfn, el)

      wfn%qat(ii, 1) = wfn%qat(ii, 1) + step
      vat(ii) = 0.5_wp*(sum(er) - sum(el))/step
   end do

   energy = 0.0_wp
   pot%vat(:, :) = 0.0_wp
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)

   if (any(abs([pot%vat] - vat) > thr)) then
      call test_failed(error, "Potential does not match")
      print '(3es20.13)', pot%vat
      print '(a)', "---"
      print '(3es20.13)', vat
      print '(a)', "---"
      print '(3es20.13)', [pot%vat] - vat
   end if
end subroutine test_p


subroutine test_e_sasa(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-8.99890404486076E-2_wp, 9.42168087556583E-2_wp,-1.49387631509499E-1_wp, &
      &-2.99114121895542E-1_wp, 4.85527734875224E-1_wp,-6.83156326406137E-2_wp, &
      & 1.50011293889337E-2_wp, 2.79368544459465E-1_wp,-1.24072452878322E-1_wp, &
      &-9.36760994051244E-2_wp,-2.19062123031622E-1_wp, 2.14538817685587E-1_wp, &
      & 3.06156726072831E-1_wp,-3.86105514712244E-1_wp,-1.51265171389388E-3_wp, &
      & 3.64255069977693E-2_wp]
   real(wp), allocatable :: rad(:), tension(:), hbond(:)

   call get_structure(mol, "MB16-43", "04")
   rad = get_vdw_rad_cosmo(mol%num)
   tension = tension_water(mol%num)
   hbond = hbond_water(mol%num)

   call test_e(error, mol, cds_input(probe=0.3_wp, nang=110, rad=rad, tension=tension), &
      & qat, -2.0260925782007388E-3_wp)
   if (allocated(error)) return

   call test_e(error, mol, cds_input(probe=0.3_wp, nang=110, rad=rad, tension=tension, &
      & hbond=hbond), qat, -0.21782880741530361_wp)

end subroutine test_e_sasa


subroutine test_g_sasa(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 2.08159387594211E-1_wp,-3.78010519998818E-1_wp, 3.36498247356244E-2_wp, &
      &-4.11556158912895E-1_wp, 8.14928196660512E-2_wp,-2.00886649303053E-1_wp, &
      & 2.44756994282684E-1_wp, 2.54580499189089E-2_wp, 2.59835128092562E-1_wp, &
      & 4.21683321877209E-1_wp, 1.37097163086023E-1_wp, 4.06951664942900E-2_wp, &
      &-1.10955378625897E-1_wp,-6.44033540918074E-2_wp,-1.91525919028143E-1_wp, &
      &-9.54898757869102E-2_wp]
   real(wp), allocatable :: rad(:), tension(:), hbond(:)

   call get_structure(mol, "MB16-43", "06")
   rad = get_vdw_rad_bondi(mol%num)
   tension = tension_water(mol%num)
   hbond = hbond_water(mol%num)
   call test_g(error, mol, cds_input(probe=2.0_wp, nang=110, rad=rad, tension=tension, &
      & hbond=hbond), qat)

end subroutine test_g_sasa


subroutine test_p_sasa(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-2.05668345919710E-1_wp,-3.99553123071811E-1_wp, 3.29242774348191E-1_wp, &
      &-3.11737933844111E-1_wp, 3.58851882478133E-2_wp, 3.21886835736497E-1_wp, &
      & 4.14743455841314E-2_wp, 2.95727359478547E-2_wp,-5.06347224522431E-1_wp, &
      & 3.43067182413129E-1_wp, 6.88373767679515E-1_wp, 7.03357390141253E-2_wp, &
      &-9.62424552888750E-2_wp,-1.32209348056625E-1_wp, 9.78998441832186E-2_wp, &
      &-3.05979982450903E-1_wp]
   real(wp), allocatable :: rad(:), tension(:), hbond(:)

   call get_structure(mol, "MB16-43", "08")
   rad = get_vdw_rad_d3(mol%num)
   tension = tension_water(mol%num)
   hbond = hbond_water(mol%num)
   call test_p(error, mol, cds_input(probe=2.2_wp, nang=110, rad=rad, tension=tension, &
      & hbond=hbond), qat)

end subroutine test_p_sasa


end module test_solvation_surface
