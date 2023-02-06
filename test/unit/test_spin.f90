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
   real(wp), parameter :: ref1 = -10.801225962675073_wp, ref0 = -10.789711366857366_wp

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
   real(wp), parameter :: ref1 = -28.341206524033051_wp, ref0 = -28.349613833732931_wp

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
   real(wp), parameter :: eref = -11.539672597844298_wp, gref(3, 8) = reshape([&
  &  4.6249752761279780E-003_wp,   3.0612696305948269E-003_wp,  -5.6747819334249106E-017_wp, &
  & -5.8129085837906090E-003_wp,   7.0212297277206437E-003_wp,   7.1150767569361228E-017_wp, &
  &  8.4737593035400759E-003_wp,  -7.8529544361057250E-003_wp,  -1.1264691165534940E-017_wp, &
  &  1.6730068671934933E-004_wp,  -2.6983454320247020E-003_wp,   2.1951168193644288E-017_wp, &
  & -2.4723837081867704E-003_wp,   1.1334200185385455E-003_wp,   1.7120845352279658E-017_wp, &
  & -1.2010888914798831E-003_wp,  -5.3290987794842807E-004_wp,   2.1500542349721535E-003_wp, &
  & -1.2010888914799251E-003_wp,  -5.3290987794850787E-004_wp,  -2.1500542349721808E-003_wp, &
  & -2.5785651914501462E-003_wp,   4.0120024717337009E-004_wp,   1.0110039674639520E-017_wp], &
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
   real(wp), parameter :: eref = -28.436616901389360_wp, gref(3, 21) = reshape([&
  &  2.7997013886914279E-015_wp,  -4.3832260954515534E-016_wp,   5.4124614084470968E-003_wp, &
  & -1.5580186280115722E-015_wp,  -3.6934877716637000E-003_wp,   1.7078904720189136E-002_wp, &
  &  1.4193184776494644E-002_wp,  -2.9571145241722070E-004_wp,  -7.6126628794408570E-004_wp, &
  & -1.2838177112578729E-003_wp,   3.3463466951857211E-003_wp,  -7.9685459655836465E-003_wp, &
  &  1.2838177112565168E-003_wp,   3.3463466951861040E-003_wp,  -7.9685459655811624E-003_wp, &
  & -1.4193184776493045E-002_wp,  -2.9571145241611541E-004_wp,  -7.6126628794560846E-004_wp, &
  &  1.4572362198890261E-016_wp,  -4.5203718825806219E-003_wp,   1.6529187663672056E-003_wp, &
  &  1.7765817706421033E-003_wp,  -4.7435151346884637E-003_wp,  -2.1933840669344169E-004_wp, &
  &  4.1897814701754319E-004_wp,  -4.8466083415670764E-003_wp,  -1.7698764351693076E-003_wp, &
  & -4.1897814701725771E-004_wp,  -4.8466083415671267E-003_wp,  -1.7698764351694899E-003_wp, &
  & -1.7765817706424534E-003_wp,  -4.7435151346885435E-003_wp,  -2.1933840669329996E-004_wp, &
  &  1.4193184776495792E-002_wp,   2.9571145241773694E-004_wp,  -7.6126628794314028E-004_wp, &
  & -1.8043868470733449E-015_wp,   3.6934877716644481E-003_wp,   1.7078904720190281E-002_wp, &
  & -1.2838177112590803E-003_wp,  -3.3463466951858260E-003_wp,  -7.9685459655852997E-003_wp, &
  &  1.7765817706418313E-003_wp,   4.7435151346885245E-003_wp,  -2.1933840669361216E-004_wp, &
  & -1.4193184776493348E-002_wp,   2.9571145241588122E-004_wp,  -7.6126628794698334E-004_wp, &
  &  3.7213848042667519E-016_wp,   4.5203718825805005E-003_wp,   1.6529187663667968E-003_wp, &
  &  1.2838177112567434E-003_wp,  -3.3463466951864826E-003_wp,  -7.9685459655803714E-003_wp, &
  &  4.1897814701761041E-004_wp,   4.8466083415671379E-003_wp,  -1.7698764351690493E-003_wp, &
  & -1.7765817706424641E-003_wp,   4.7435151346884811E-003_wp,  -2.1933840669311868E-004_wp, &
  & -4.1897814701732786E-004_wp,   4.8466083415670764E-003_wp,  -1.7698764351696782E-003_wp], &
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
