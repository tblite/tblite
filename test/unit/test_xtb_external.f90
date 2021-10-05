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

module test_xtb_external
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_container, only : container_type, container_cache
   use tblite_context_type, only : context_type
   use tblite_external_field, only : electric_field
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: collect_xtb_external

   real(wp), parameter :: acc = 0.01_wp
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: thr2 = 1e+4_wp*sqrt(epsilon(1.0_wp))
   real(wp), parameter :: kt = 300.0_wp * 3.166808578545117e-06_wp

  real(wp), parameter :: aatoau = 1.0_wp / 0.529177249_wp, &
     & ctoau = 1.0_wp / 1.60217653e-19_wp, jtoau = 1.0_wp / 4.3597441775e-18_wp
  !> Convert V/Å = J/(C·Å) to atomic units
  real(wp), parameter :: vatoau = jtoau / (ctoau * aatoau)

  type, extends(container_type) :: empty_interaction
  end type empty_interaction

contains


!> Collect all exported unit tests
subroutine collect_xtb_external(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("gfn1-efield", test_e_mb01), &
      new_unittest("gfn2-efield", test_e_mb02), &
      new_unittest("gfn1-dipole", test_d_mb03), &
      new_unittest("gfn2-dipole", test_d_mb04), &
      new_unittest("gfn1-empty", test_g_mb05), &
      new_unittest("gfn2-empty", test_g_mb06) &
      ]

end subroutine collect_xtb_external


subroutine test_e_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   class(container_type), allocatable :: cont
   real(wp) :: energy
   real(wp), parameter :: ref1 = -33.092759804236330_wp, ref0 = -33.040345086576387_wp

   call get_structure(mol, "MB16-43", "01")
   energy = 0.0_wp

   call new_gfn1_calculator(calc, mol)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, kt)

   cont = electric_field([-2.0_wp, 0.0_wp, 0.0_wp]*vatoau)
   call calc%push_back(cont)

   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)

   call check(error, energy, ref1, thr=thr)

   call calc%pop(cont)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)

   call check(error, energy, ref0, thr=thr)

end subroutine test_e_mb01


subroutine test_e_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   class(container_type), allocatable :: cont
   real(wp) :: energy
   real(wp), allocatable :: gradient(:, :), sigma(:, :)
   real(wp), parameter :: ref1 = -24.157134025160921_wp, ref0 = -24.069929678894123_wp

   call get_structure(mol, "MB16-43", "02")
   allocate(gradient(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp

   call new_gfn2_calculator(calc, mol)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, kt)

   cont = electric_field([0.0_wp, sqrt(2.0_wp), -sqrt(2.0_wp)]*vatoau)
   call calc%push_back(cont)

   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, verbosity=0)

   call check(error, energy, ref1, thr=thr)

   call calc%pop(cont)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, verbosity=0)

   call check(error, energy, ref0, thr=thr)

end subroutine test_e_mb02


subroutine test_d_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn, wfn0
   class(container_type), allocatable :: cont
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp) :: energy, efield(3), er, el, numdip(3), dipole(3)
   integer :: i

   call get_structure(mol, "MB16-43", "03")
   energy = 0.0_wp
   efield(:) = 0.0_wp

   call new_gfn1_calculator(calc, mol)
   call new_wavefunction(wfn0, mol%nat, calc%bas%nsh, calc%bas%nao, kt)

   cont = electric_field(efield)
   call calc%push_back(cont)

   call xtb_singlepoint(ctx, mol, calc, wfn0, acc, energy, verbosity=0)
   dipole(:) = matmul(mol%xyz, wfn0%qat) + sum(wfn0%dpat, 2)

   do i = 1, 3
      wfn = wfn0
      efield(i) = efield(i) + step
      call calc%pop(cont)
      cont = electric_field(efield)
      call calc%push_back(cont)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, er, verbosity=0)

      wfn = wfn0
      efield(i) = efield(i) - 2*step
      call calc%pop(cont)
      cont = electric_field(efield)
      call calc%push_back(cont)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, el, verbosity=0)

      efield(i) = efield(i) + step
      numdip(i) = -0.5_wp * (er - el) / step
   end do

   if (any(abs(dipole - numdip) > thr2)) then
      call test_failed(error, "Numerical dipole moment does not match")
      print '(3es21.14)', dipole
      print '("---")'
      print '(3es21.14)', numdip
      print '("---")'
      print '(3es21.14)', dipole - numdip
   end if

end subroutine test_d_mb03


subroutine test_d_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn, wfn0
   class(container_type), allocatable :: cont
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp) :: energy, efield(3), er, el, numdip(3), dipole(3)
   integer :: i

   call get_structure(mol, "MB16-43", "04")
   energy = 0.0_wp
   efield(:) = 0.0_wp

   call new_gfn2_calculator(calc, mol)
   call new_wavefunction(wfn0, mol%nat, calc%bas%nsh, calc%bas%nao, kt)

   cont = electric_field(efield)
   call calc%push_back(cont)

   call xtb_singlepoint(ctx, mol, calc, wfn0, acc, energy, verbosity=0)
   dipole(:) = matmul(mol%xyz, wfn0%qat) + sum(wfn0%dpat, 2)

   do i = 1, 3
      wfn = wfn0
      efield(i) = efield(i) + step
      call calc%pop(cont)
      cont = electric_field(efield)
      call calc%push_back(cont)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, er, verbosity=0)

      wfn = wfn0
      efield(i) = efield(i) - 2*step
      call calc%pop(cont)
      cont = electric_field(efield)
      call calc%push_back(cont)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, el, verbosity=0)

      efield(i) = efield(i) + step
      numdip(i) = -0.5_wp * (er - el) / step
   end do

   if (any(abs(dipole - numdip) > thr2)) then
      call test_failed(error, "Numerical dipole moment does not match")
      print '(3es21.14)', dipole
      print '("---")'
      print '(3es21.14)', numdip
      print '("---")'
      print '(3es21.14)', dipole - numdip
   end if

end subroutine test_d_mb04


subroutine test_g_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   class(container_type), allocatable :: cont
   real(wp) :: energy
   real(wp), allocatable :: gradient(:, :), sigma(:, :)
   real(wp), parameter :: eref = -27.735987610325889_wp, gref(3, 16) = reshape([&
      & -1.8235429544448E-3_wp,  1.3181894439006E-3_wp,  8.5822519510144E-3_wp, &
      &  1.9768865248121E-2_wp, -7.2945842651456E-3_wp, -6.0438086560557E-3_wp, &
      &  3.9881019940576E-3_wp,  6.1490714375871E-3_wp,  4.8038465530667E-3_wp, &
      &  9.1166248719117E-3_wp,  8.1189091588512E-4_wp, -2.0603363004095E-3_wp, &
      &  6.9329682115644E-3_wp, -9.9964557147168E-4_wp, -5.7443369652198E-3_wp, &
      & -4.8675411574309E-3_wp, -3.0627052769839E-3_wp, -2.0109062332000E-3_wp, &
      &  4.1644098339087E-3_wp, -2.0340497310258E-3_wp, -6.0780858751439E-4_wp, &
      & -2.2038175531200E-2_wp, -1.1147439653379E-4_wp,  7.8831226805861E-3_wp, &
      & -1.1352594800217E-3_wp, -9.0312601610063E-3_wp,  5.4144106194533E-3_wp, &
      & -4.3176311990670E-4_wp,  7.2864066035523E-3_wp,  1.4638861555741E-3_wp, &
      & -9.5239504517465E-3_wp, -3.0619341605843E-4_wp, -1.1800933866746E-2_wp, &
      & -3.1721506980028E-3_wp,  1.2341293868359E-2_wp, -6.9305589804901E-3_wp, &
      &  3.8470358756443E-3_wp, -5.3096535160343E-3_wp, -4.7665012734588E-3_wp, &
      &  4.0181327665139E-3_wp, -6.6099642502266E-3_wp, -4.9460830121791E-3_wp, &
      &  4.6662894517864E-3_wp,  2.5460207907372E-3_wp,  1.5317323586343E-3_wp, &
      & -1.3510044860754E-2_wp,  4.3066575244654E-3_wp,  1.5232023556944E-2_wp],&
      & shape(gref))

   call get_structure(mol, "MB16-43", "05")
   allocate(gradient(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call new_gfn2_calculator(calc, mol)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, kt)

   cont = empty_interaction()
   call calc%push_back(cont)

   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, verbosity=0)

   call check(error, energy, eref, thr=thr)
   if (allocated(error)) return
   call check(error, all(abs(gradient - gref) < thr))

end subroutine test_g_mb05


subroutine test_g_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   class(container_type), allocatable :: cont
   real(wp) :: energy
   real(wp), allocatable :: gradient(:, :), sigma(:, :)
   real(wp), parameter :: eref = -18.559531624771843_wp, gref(3, 16) = reshape([&
      & -7.6885542678065E-3_wp,  2.1063888360685E-2_wp,  2.4121523424853E-3_wp, &
      &  2.9080591624961E-3_wp, -1.7501761976798E-3_wp,  3.8558280713210E-3_wp, &
      &  1.5421922710101E-3_wp,  1.3237640678443E-3_wp,  1.0182491425554E-3_wp, &
      & -1.9145261823387E-2_wp, -2.2338434369955E-2_wp,  5.4858599573081E-3_wp, &
      & -4.3974565005313E-3_wp, -1.1052970054546E-2_wp, -2.0343428051935E-2_wp, &
      &  9.6286425110520E-3_wp,  9.3349579578424E-3_wp,  4.2851303903926E-3_wp, &
      &  3.4030110887005E-3_wp,  1.3332626054614E-2_wp,  3.8129668684993E-4_wp, &
      & -2.3215966669115E-3_wp,  8.7910099194196E-4_wp, -2.0312084579466E-3_wp, &
      &  1.7342359254211E-2_wp,  1.5273548827631E-2_wp, -5.0199979432449E-3_wp, &
      & -3.3507318065947E-3_wp, -1.5974401847186E-2_wp,  8.3152808490089E-3_wp, &
      & -1.8257802211007E-3_wp, -3.8668297167786E-3_wp, -1.4746216779591E-3_wp, &
      &  1.3359190167041E-3_wp, -9.8876405264024E-4_wp,  6.4383472418861E-3_wp, &
      & -2.8465113729483E-3_wp,  1.6299038735390E-3_wp,  1.1044803571726E-3_wp, &
      &  4.0033955200458E-3_wp, -1.0314447785249E-3_wp,  4.0504720509594E-3_wp, &
      &  1.7256805104575E-3_wp, -6.6058663378452E-3_wp, -3.9807435012693E-3_wp, &
      & -3.1336667539659E-4_wp,  7.7109722105795E-4_wp, -4.4970974575845E-3_wp],&
      & shape(gref))

   call get_structure(mol, "MB16-43", "06")
   allocate(gradient(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call new_gfn2_calculator(calc, mol)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, kt)

   cont = empty_interaction()
   call calc%push_back(cont)

   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, verbosity=0)

   call check(error, energy, eref, thr=thr)
   if (allocated(error)) return
   call check(error, all(abs(gradient - gref) < thr))

end subroutine test_g_mb06


end module test_xtb_external
