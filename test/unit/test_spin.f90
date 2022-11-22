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
   real(wp), parameter :: ref1 = -10.802158467973536_wp, ref0 = -10.789711352994029_wp

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
   real(wp), parameter :: eref = -11.540515632359202_wp, gref(3, 8) = reshape([&
   &  4.2839758814104774E-003_wp,   2.8415043845744620E-003_wp,   2.0735366548785272E-017_wp, &
   &  -5.6230454193085854E-003_wp,   7.3528998110541062E-003_wp,   1.7374339814080209E-017_wp, &
   &  8.4826457719626765E-003_wp,  -7.8808435812119330E-003_wp,  -3.2670907808545618E-017_wp, &
   &  1.4568794973941307E-004_wp,  -2.4493593451781650E-003_wp,   9.9156919353589606E-018_wp, &
   & -2.2504987498990882E-003_wp,   1.0183045601019182E-003_wp,  -1.7988696781190693E-018_wp, &
   & -1.2283924647833188E-003_wp,  -6.4217353074644073E-004_wp,   2.2301318032430474E-003_wp, &
   & -1.2283924647833711E-003_wp,  -6.4217353074638446E-004_wp,  -2.2301318032430760E-003_wp, &
   & -2.5819805043382172E-003_wp,   4.0184123215238727E-004_wp,   6.2793464055717471E-018_wp], &
     & shape(gref))


   call rse43_p10(mol)
   allocate(gradient(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

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
   real(wp), parameter :: eref = -28.468935665440110_wp, gref(3, 21) = reshape([&
  & -1.4111501138620863E-014_wp,   7.8152816249815737E-015_wp,   9.2358786075009362E-004_wp, &
  &  9.3150964390106904E-015_wp,  -1.4950226428480252E-003_wp,   1.4800995373278392E-002_wp, &
  &  1.2235687416282660E-002_wp,   9.6237573962382580E-005_wp,   8.4893856278399880E-004_wp, &
  & -4.8808333658941972E-005_wp,   1.9385955285310765E-003_wp,  -7.4741281680472801E-003_wp, &
  &  4.8808333668234858E-005_wp,   1.9385955285268364E-003_wp,  -7.4741281680671410E-003_wp, &
  & -1.2235687416293498E-002_wp,   9.6237573952830962E-005_wp,   8.4893856279680171E-004_wp, &
  & -1.0522051242511776E-015_wp,  -4.7507380366085435E-003_wp,   1.7083241206784932E-003_wp, &
  &  1.8053024455472434E-003_wp,  -4.7380488806736155E-003_wp,  -1.3098313338247091E-004_wp, &
  &  4.5427087875652610E-004_wp,  -4.6025889342161666E-003_wp,  -1.7293839735182659E-003_wp, &
  & -4.5427087875726422E-004_wp,  -4.6025889342162195E-003_wp,  -1.7293839735160789E-003_wp, &
  & -1.8053024455454850E-003_wp,  -4.7380488806730639E-003_wp,  -1.3098313338377962E-004_wp, &
  &  1.2235687416277187E-002_wp,  -9.6237573963227891E-005_wp,   8.4893856278902224E-004_wp, &
  &  6.8875613581106781E-015_wp,   1.4950226428423906E-003_wp,   1.4800995373269408E-002_wp, &
  & -4.8808333659240913E-005_wp,  -1.9385955285299141E-003_wp,  -7.4741281680470381E-003_wp, &
  &  1.8053024455486461E-003_wp,   4.7380488806740371E-003_wp,  -1.3098313338274719E-004_wp, &
  & -1.2235687416284242E-002_wp,  -9.6237573956721201E-005_wp,   8.4893856279851887E-004_wp, &
  & -6.4472354345589481E-016_wp,   4.7507380366092816E-003_wp,   1.7083241206810073E-003_wp, &
  &  4.8808333666354974E-005_wp,  -1.9385955285270873E-003_wp,  -7.4741281680611648E-003_wp, &
  &  4.5427087875731393E-004_wp,   4.6025889342163426E-003_wp,  -1.7293839735191541E-003_wp, &
  & -1.8053024455479106E-003_wp,   4.7380488806734438E-003_wp,  -1.3098313338341354E-004_wp, &
  & -4.5427087875786975E-004_wp,   4.6025889342161432E-003_wp,  -1.7293839735178563E-003_wp], &
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
   & -1.97051959765227E+00_wp,   -8.65723337874754E-01_wp,    0.00000000000000E+00_wp, &     
   &  3.50984622791913E-01_wp,    6.86290619844032E-01_wp,    0.00000000000000E+00_wp, &      
   &  2.50609985217434E+00_wp,   -9.34496149122418E-01_wp,    0.00000000000000E+00_wp, &      
   & -1.83649606109455E+00_wp,   -2.90299181092583E+00_wp,    0.00000000000000E+00_wp, &      
   & -3.80466245712260E+00_wp,    3.49832428602470E-02_wp,    0.00000000000000E+00_wp, &      
   &  3.73555581511497E-01_wp,    1.94431040908594E+00_wp,   -1.66596178649581E+00_wp, &      
   &  3.73555581511497E-01_wp,    1.94431040908594E+00_wp,    1.66596178649581E+00_wp, &      
   &  4.00748247788016E+00_wp,    9.33166170468600E-02_wp,    0.00000000000000E+00_wp], &      
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
  &  0.00000000000000E+00_wp,    0.00000000000000E+00_wp,   -6.04468452830504E-02_wp, &      
  &  0.00000000000000E+00_wp,    3.19613712523833E+00_wp,    2.30877824528580E+00_wp, &      
  &  2.18828801115897E+00_wp,    3.32943780995850E+00_wp,    7.02499485857345E-01_wp, &      
  &  1.33235791539260E+00_wp,    3.55640652898451E+00_wp,   -1.83908673090077E+00_wp, &      
  & -1.33235791539260E+00_wp,    3.55640652898451E+00_wp,   -1.83908673090077E+00_wp, &      
  & -2.18828801115897E+00_wp,    3.32943780995850E+00_wp,    7.02499485857345E-01_wp, &      
  &  0.00000000000000E+00_wp,    3.10509505378016E+00_wp,    4.34935395653655E+00_wp, &      
  &  4.13810718850644E+00_wp,    3.28428734944129E+00_wp,    1.31235006648465E+00_wp, &      
  &  2.52190264478215E+00_wp,    3.60569548880831E+00_wp,   -3.50208900904436E+00_wp, &      
  & -2.52190264478215E+00_wp,    3.60569548880831E+00_wp,   -3.50208900904436E+00_wp, &      
  & -4.13810718850644E+00_wp,    3.28428734944129E+00_wp,    1.31235006648465E+00_wp, &      
  &  2.18828801115897E+00_wp,   -3.32943780995850E+00_wp,    7.02499485857345E-01_wp, &      
  &  0.00000000000000E+00_wp,   -3.19613712523833E+00_wp,    2.30877824528580E+00_wp, &      
  &  1.33235791539260E+00_wp,   -3.55640652898451E+00_wp,   -1.83908673090077E+00_wp, &      
  &  4.13810718850644E+00_wp,   -3.28428734944129E+00_wp,    1.31235006648465E+00_wp, &      
  & -2.18828801115897E+00_wp,   -3.32943780995850E+00_wp,    7.02499485857345E-01_wp, &      
  &  0.00000000000000E+00_wp,   -3.10509505378016E+00_wp,    4.34935395653655E+00_wp, &      
  & -1.33235791539260E+00_wp,   -3.55640652898451E+00_wp,   -1.83908673090077E+00_wp, &      
  &  2.52190264478215E+00_wp,   -3.60569548880831E+00_wp,   -3.50208900904436E+00_wp, &      
  & -4.13810718850644E+00_wp,   -3.28428734944129E+00_wp,    1.31235006648465E+00_wp, &      
  & -2.52190264478215E+00_wp,   -3.60569548880831E+00_wp,   -3.50208900904436E+00_wp], &      
 & shape(xyz))
   integer, parameter :: uhf = 2
   call new(self, sym, xyz, uhf=uhf)
end subroutine crcp2


end module test_spin
