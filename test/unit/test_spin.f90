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

module test_spin
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_basis_type, only : basis_type
   use tblite_container, only : container_type, container_cache
   use tblite_context_type, only : context_type
   use tblite_data_spin, only : get_spin_constant
   use tblite_spin, only : spin_polarization, new_spin_polarization
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: collect_spin

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
subroutine collect_spin(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("gfn1-e-spin", test_e_crcp2), &
      new_unittest("gfn2-e-spin", test_e_p10), &
      new_unittest("gfn1-g-spin", test_g_p10), &
      new_unittest("gfn2-g-spin", test_g_crcp2) &
      ]

end subroutine collect_spin


subroutine test_e_p10(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   class(container_type), allocatable :: cont
   real(wp) :: energy
   real(wp), parameter :: ref1 = -10.802158458119459_wp, ref0 = -10.789711352994029_wp

   call rse43_p10(mol)
   energy = 0.0_wp

   call new_gfn2_calculator(calc, mol)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 2, kt)

   block
      type(spin_polarization), allocatable :: spin
      real(wp), allocatable :: wll(:, :, :)
      allocate(spin)
      call get_spin_constants(wll, mol, calc%bas)
      call new_spin_polarization(spin, mol, wll, calc%bas%nsh_id)
      call move_alloc(spin, cont)
      call calc%push_back(cont)
   end block

   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)

   call check(error, energy, ref1, thr=thr)
   if (allocated(error)) return

   call calc%pop(cont)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)

   call check(error, energy, ref0, thr=thr)

end subroutine test_e_p10


subroutine test_e_crcp2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   class(container_type), allocatable :: cont
   real(wp) :: energy
   real(wp), allocatable :: gradient(:, :), sigma(:, :)
   real(wp), parameter :: ref1 = -28.373975300991656_wp, ref0 = -28.349613833733063_wp

   call crcp2(mol)
   allocate(gradient(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp

   call new_gfn1_calculator(calc, mol)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 2, kt)

   block
      type(spin_polarization), allocatable :: spin
      real(wp), allocatable :: wll(:, :, :)
      allocate(spin)
      call get_spin_constants(wll, mol, calc%bas)
      call new_spin_polarization(spin, mol, wll, calc%bas%nsh_id)
      call move_alloc(spin, cont)
      call calc%push_back(cont)
   end block

   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, verbosity=0)

   call check(error, energy, ref1, thr=thr)
   if (allocated(error)) return

   mol%uhf = 0
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 2, kt)

   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, verbosity=0)

   call check(error, energy, ref0, thr=thr)
   if (allocated(error)) return

   call calc%pop(cont)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, verbosity=0)

   call check(error, energy, ref0, thr=thr)

end subroutine test_e_crcp2


subroutine test_g_p10(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   class(container_type), allocatable :: cont
   real(wp) :: energy
   real(wp), allocatable :: gradient(:, :), sigma(:, :)
   real(wp), parameter :: eref = -10.8021584679735_wp, gref(3, 8) = reshape([&
 &  4.8823767545216577E-03_wp,  2.9532959748048769E-03_wp, -3.4762232155316486E-18_wp, &
 & -3.8646457247001344E-03_wp,  8.7557577403917213E-03_wp,  8.1965684239904135E-17_wp, &
 &  7.3902847289607299E-03_wp, -8.6908539939121350E-03_wp, -1.5788694136820158E-17_wp, &
 & -7.4140027505004906E-04_wp, -6.5698737289731001E-03_wp, -2.4394933519908687E-18_wp, &
 & -5.6853104768989585E-03_wp,  1.8547919494240630E-03_wp,  7.4811174128806384E-18_wp, &
 & -1.0993280366217745E-03_wp,  4.6356506108841702E-04_wp,  1.5684980489191667E-03_wp, &
 & -1.0993280366217060E-03_wp,  4.6356506108842049E-04_wp, -1.5684980489192484E-03_wp, &
 &  2.1735106641023914E-04_wp,  7.6975193608779685E-04_wp,  6.1798729738285702E-18_wp], &
     & shape(gref))

   call rse43_p10(mol)
   allocate(gradient(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call new_gfn2_calculator(calc, mol)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 2, kt)

   block
      type(spin_polarization), allocatable :: spin
      real(wp), allocatable :: wll(:, :, :)
      allocate(spin)
      call get_spin_constants(wll, mol, calc%bas)
      call new_spin_polarization(spin, mol, wll, calc%bas%nsh_id)
      call move_alloc(spin, cont)
      call calc%push_back(cont)
   end block

   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, verbosity=0)

   call check(error, energy, eref, thr=thr)
   if (allocated(error)) return
   call check(error, all(abs(gradient - gref) < thr))

end subroutine test_g_p10


subroutine test_g_crcp2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   class(container_type), allocatable :: cont
   real(wp) :: energy
   real(wp), allocatable :: gradient(:, :), sigma(:, :)
   real(wp), parameter :: eref = -28.4689356654401_wp, gref(3, 21) = reshape([&
 &  1.2465018166957220E-14_wp,  5.2745758065061987E-16_wp,  9.2343247960639959E-04_wp, &
 & -7.5755720995450258E-15_wp, -1.4946698886261040E-03_wp,  1.4801049719676116E-02_wp, &
 &  1.2235710223571861E-02_wp,  9.6568994280106301E-05_wp,  8.4886419950261328E-04_wp, &
 & -4.8943631566281012E-05_wp,  1.9388653721408560E-03_wp, -7.4740329469348242E-03_wp, &
 &  4.8943631558773048E-05_wp,  1.9388653721439488E-03_wp, -7.4740329469184207E-03_wp, &
 & -1.2235710223562847E-02_wp,  9.6568994287620446E-05_wp,  8.4886419949170111E-04_wp, &
 &  7.8520670246701161E-16_wp, -4.7507814854341667E-03_wp,  1.7083299906172900E-03_wp, &
 &  1.8053306517207197E-03_wp, -4.7380777534232121E-03_wp, -1.3096026161590606E-04_wp, &
 &  4.5430538932696222E-04_wp, -4.6026094828900081E-03_wp, -1.7294189660027603E-03_wp, &
 & -4.5430538932641167E-04_wp, -4.6026094828900167E-03_wp, -1.7294189660045371E-03_wp, &
 & -1.8053306517222033E-03_wp, -4.7380777534236276E-03_wp, -1.3096026161485376E-04_wp, &
 &  1.2235710223570812E-02_wp, -9.6568994280520412E-05_wp,  8.4886419950127049E-04_wp, &
 & -7.5648736850663688E-15_wp,  1.4946698886260854E-03_wp,  1.4801049719675817E-02_wp, &
 & -4.8943631566335927E-05_wp, -1.9388653721411687E-03_wp, -7.4740329469333835E-03_wp, &
 &  1.8053306517208780E-03_wp,  4.7380777534232034E-03_wp, -1.3096026161570797E-04_wp, &
 & -1.2235710223562603E-02_wp, -9.6568994287553022E-05_wp,  8.4886419949212275E-04_wp, &
 &  7.7894493033091442E-16_wp,  4.7507814854342969E-03_wp,  1.7083299906173477E-03_wp, &
 &  4.8943631559400069E-05_wp, -1.9388653721439989E-03_wp, -7.4740329469187624E-03_wp, &
 &  4.5430538932713971E-04_wp,  4.6026094828899638E-03_wp, -1.7294189660029661E-03_wp, &
 & -1.8053306517220801E-03_wp,  4.7380777534236874E-03_wp, -1.3096026161482866E-04_wp, &
 & -4.5430538932661718E-04_wp,  4.6026094828900983E-03_wp, -1.7294189660046091E-03_wp], &
   & shape(gref))

   call crcp2(mol)
   allocate(gradient(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call new_gfn2_calculator(calc, mol)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 2, kt)

   block
      type(spin_polarization), allocatable :: spin
      real(wp), allocatable :: wll(:, :, :)
      allocate(spin)
      call get_spin_constants(wll, mol, calc%bas)
      call new_spin_polarization(spin, mol, wll, calc%bas%nsh_id)
      call move_alloc(spin, cont)
      call calc%push_back(cont)
   end block

   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, verbosity=0)

   call check(error, energy, eref, thr=thr)
   if (allocated(error)) return
   call check(error, all(abs(gradient - gref) < thr))

end subroutine test_g_crcp2


subroutine get_spin_constants(wll, mol, bas)
   real(wp), allocatable, intent(out) :: wll(:, :, :)
   type(structure_type), intent(in) :: mol
   type(basis_type), intent(in) :: bas

   integer :: izp, ish, jsh, il, jl

   allocate(wll(bas%nsh, bas%nsh, mol%nid), source=0.0_wp)

   do izp = 1, mol%nid
      do ish = 1, bas%nsh_id(izp)
         il = bas%cgto(ish, izp)%ang
         do jsh = 1, bas%nsh_id(izp)
            jl = bas%cgto(jsh, izp)%ang
            wll(jsh, ish, izp) = get_spin_constant(jl, il, mol%num(izp))
         end do
      end do
   end do
end subroutine get_spin_constants


subroutine rse43_p10(self)
   type(structure_type), intent(out) :: self
   integer, parameter :: nat = 8
   character(len=*), parameter :: sym(nat) = [character(len=4)::&
      & "C", "C", "O", "H", "H", "H", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & -1.97051959765227_wp, -0.86572333787476_wp,  0.00000000000000_wp, &
      &  0.35098462279192_wp,  0.68629061984403_wp,  0.00000000000000_wp, &
      &  2.50609985217434_wp, -0.93449614912242_wp,  0.00000000000000_wp, &
      & -1.83649606109455_wp, -2.90299181092582_wp,  0.00000000000000_wp, &
      & -3.80466245712260_wp,  0.03498324286025_wp,  0.00000000000000_wp, &
      &  0.37355558151150_wp,  1.94431040908593_wp, -1.66596178649581_wp, &
      &  0.37355558151150_wp,  1.94431040908593_wp,  1.66596178649581_wp, &
      &  4.00748247788017_wp,  0.09331661704686_wp,  0.00000000000000_wp],&
      & shape(xyz))
   integer, parameter :: uhf = 1
   call new(self, sym, xyz, uhf=uhf)
end subroutine rse43_p10


subroutine crcp2(self)
   type(structure_type), intent(out) :: self
   integer, parameter :: nat = 21
   character(len=*), parameter :: sym(nat) = [character(len=4)::&
      & "Cr", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "C", "C", "C", &
      & "H", "C", "H", "C", "H", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      &  0.00000000000000_wp,  0.00000000000000_wp, -0.06044684528305_wp, &
      &  0.00000000000000_wp,  3.19613712523832_wp,  2.30877824528581_wp, &
      &  2.18828801115897_wp,  3.32943780995849_wp,  0.70249948585734_wp, &
      &  1.33235791539260_wp,  3.55640652898451_wp, -1.83908673090077_wp, &
      & -1.33235791539260_wp,  3.55640652898451_wp, -1.83908673090077_wp, &
      & -2.18828801115897_wp,  3.32943780995849_wp,  0.70249948585734_wp, &
      &  0.00000000000000_wp,  3.10509505378016_wp,  4.34935395653655_wp, &
      &  4.13810718850644_wp,  3.28428734944129_wp,  1.31235006648465_wp, &
      &  2.52190264478214_wp,  3.60569548880830_wp, -3.50208900904435_wp, &
      & -2.52190264478214_wp,  3.60569548880830_wp, -3.50208900904435_wp, &
      & -4.13810718850644_wp,  3.28428734944129_wp,  1.31235006648465_wp, &
      &  2.18828801115897_wp, -3.32943780995849_wp,  0.70249948585734_wp, &
      &  0.00000000000000_wp, -3.19613712523832_wp,  2.30877824528581_wp, &
      &  1.33235791539260_wp, -3.55640652898451_wp, -1.83908673090077_wp, &
      &  4.13810718850644_wp, -3.28428734944129_wp,  1.31235006648465_wp, &
      & -2.18828801115897_wp, -3.32943780995849_wp,  0.70249948585734_wp, &
      &  0.00000000000000_wp, -3.10509505378016_wp,  4.34935395653655_wp, &
      & -1.33235791539260_wp, -3.55640652898451_wp, -1.83908673090077_wp, &
      &  2.52190264478214_wp, -3.60569548880830_wp, -3.50208900904435_wp, &
      & -4.13810718850644_wp, -3.28428734944129_wp,  1.31235006648465_wp, &
      & -2.52190264478214_wp, -3.60569548880830_wp, -3.50208900904435_wp],&
      & shape(xyz))
   integer, parameter :: uhf = 2
   call new(self, sym, xyz, uhf=uhf)
end subroutine crcp2


end module test_spin
