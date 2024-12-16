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
   real(wp), parameter :: ref1 = -10.8011953472361_wp, ref0 = -10.7897113668574_wp

   call rse43_p10(mol)
   energy = 0.0_wp

   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return
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
   real(wp), parameter :: ref1 = -28.3763344264682_wp, ref0 = -28.3496138337329_wp

   call crcp2(mol)
   allocate(gradient(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp

   call new_gfn1_calculator(calc, mol, error)
   if (allocated(error)) return
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
   real(wp), parameter :: eref = -11.5396425107754_wp, gref(3, 8) = reshape([&
  &  4.617922713076135E-003_wp,   3.055183799431468E-003_wp,  -8.709334070555372E-019_wp, &
  & -5.803572662303740E-003_wp,   7.033951465230674E-003_wp,  -4.168577567334383E-017_wp, &
  &  8.474107442195294E-003_wp,  -7.855009709204370E-003_wp,  -1.908380750365926E-017_wp, &
  &  1.673932557644197E-004_wp,  -2.701027978285800E-003_wp,  -8.348359302333630E-018_wp, &
  & -2.474761698194612E-003_wp,   1.134460559816175E-003_wp,   1.000175684686357E-017_wp, &
  & -1.201151483881359E-003_wp,  -5.344335356455324E-004_wp,   2.150192604901390E-003_wp, &
  & -1.201151483881335E-003_wp,  -5.344335356455487E-004_wp,  -2.150192604901349E-003_wp, &
  & -2.578786082774796E-003_wp,   4.013089343029468E-004_wp,   6.873111356709724E-019_wp], &
     & shape(gref))


   call rse43_p10(mol)
   allocate(gradient(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call new_gfn1_calculator(calc, mol, error)
   if (allocated(error)) return
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
   real(wp), parameter :: eref = -28.4712820267740_wp, gref(3, 21) = reshape([&
  & -3.615380752790689E-014_wp,  -6.128593609364818E-015_wp,   7.446373803562436E-004_wp, &
  &  2.255179235373212E-014_wp,  -1.387678543360229E-003_wp,   1.455234962155995E-002_wp, &
  &  1.203721801045242E-002_wp,   8.661445659815034E-005_wp,   9.573029611071851E-004_wp, &
  &  3.200216258616877E-005_wp,   1.839889792025707E-003_wp,  -7.422807303885855E-003_wp, &
  & -3.200216256614849E-005_wp,   1.839889792017230E-003_wp,  -7.422807303928467E-003_wp, &
  & -1.203721801047703E-002_wp,   8.661445657771957E-005_wp,   9.573029611334458E-004_wp, &
  & -2.268978408809382E-015_wp,  -4.763291809188703E-003_wp,   1.720252959073481E-003_wp, &
  &  1.815694181767953E-003_wp,  -4.732567673037307E-003_wp,  -1.240992083303616E-004_wp, &
  &  4.599245077263882E-004_wp,  -4.584244909899227E-003_wp,  -1.732857084288763E-003_wp, &
  & -4.599245077278532E-004_wp,  -4.584244909898971E-003_wp,  -1.732857084284529E-003_wp, &
  & -1.815694181764320E-003_wp,  -4.732567673036048E-003_wp,  -1.240992083328980E-004_wp, &
  &  1.203721801045781E-002_wp,  -8.661445659636432E-005_wp,   9.573029611068801E-004_wp, &
  &  2.332366001737109E-014_wp,   1.387678543363916E-003_wp,   1.455234962156562E-002_wp, &
  &  3.200216258490514E-005_wp,  -1.839889792026087E-003_wp,  -7.422807303889528E-003_wp, &
  &  1.815694181766321E-003_wp,   4.732567673036918E-003_wp,  -1.240992083306916E-004_wp, &
  & -1.203721801048352E-002_wp,  -8.661445657521662E-005_wp,   9.573029611335994E-004_wp, &
  & -2.412213268822006E-015_wp,   4.763291809188074E-003_wp,   1.720252959071919E-003_wp, &
  & -3.200216256423079E-005_wp,  -1.839889792017341E-003_wp,  -7.422807303933058E-003_wp, &
  &  4.599245077256235E-004_wp,   4.584244909899196E-003_wp,  -1.732857084287647E-003_wp, &
  & -1.815694181762447E-003_wp,   4.732567673035733E-003_wp,  -1.240992083334584E-004_wp, &
  & -4.599245077272394E-004_wp,   4.584244909898960E-003_wp,  -1.732857084283062E-003_wp], &
   & shape(gref))

   call crcp2(mol)
   allocate(gradient(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return
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
