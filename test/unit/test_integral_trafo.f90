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

module test_integral_trafo
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_blas, only : gemm
   use tblite_context_type, only : context_type
   use tblite_integral_trafo, only : transform0, adjoint_transform0, &
      & adjoint_transform1, adjoint_transform2
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction, eeq_guess
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint

   implicit none
   private

   public :: collect_integral_trafo

   real(wp), parameter :: kt = 300.0_wp * 3.166808578545117e-06_wp
   real(wp), parameter :: acc = 0.01_wp

   real(wp), parameter :: thr = 1e4*epsilon(1.0_wp)

contains


!> Collect all exported unit tests
subroutine collect_integral_trafo(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("trafo-ss", test_trafo_ss), &
      new_unittest("trafo-sp", test_trafo_sp), &
      new_unittest("trafo-sd", test_trafo_sd), &
      new_unittest("trafo-sf", test_trafo_sf), &
      new_unittest("trafo-sg", test_trafo_sg), &
      new_unittest("trafo-pp", test_trafo_pp), &
      new_unittest("trafo-pd", test_trafo_pd), &
      new_unittest("trafo-pf", test_trafo_pf), &
      new_unittest("trafo-pg", test_trafo_pg), &
      new_unittest("trafo-dd", test_trafo_dd), &
      new_unittest("trafo-df", test_trafo_df), &
      new_unittest("trafo-dg", test_trafo_dg), &
      new_unittest("trafo-ff", test_trafo_ff), &
      new_unittest("trafo-fg", test_trafo_fg), &
      new_unittest("trafo-gg", test_trafo_gg), &
      new_unittest("trafo-pcl", test_trafo_pcl) &
      ]

end subroutine collect_integral_trafo


subroutine test_trafo_pair(error, lj, li)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Angular momentum of bra
   integer, intent(in) :: lj
   !> Angular momentum of ket
   integer, intent(in) :: li

   ! Transform neither bra nor ket
   call test_trafo_adjoint(error, lj, li, .false., .false.)
   if (allocated(error)) return

   ! Transform bra only
   call test_trafo_adjoint(error, lj, li, .true., .false.)
   if (allocated(error)) return

   ! Transform ket only
   call test_trafo_adjoint(error, lj, li, .false., .true.)
   if (allocated(error)) return

   ! Transform both bra and ket
   call test_trafo_adjoint(error, lj, li, .true., .true.)
   if (allocated(error)) return

end subroutine test_trafo_pair

subroutine test_trafo_adjoint(error, lj, li, bra, ket)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Angular momentum of bra
   integer, intent(in) :: lj
   !> Angular momentum of ket
   integer, intent(in) :: li
   !> Whether to transform the bra
   logical, intent(in) :: bra
   !> Whether to transform the ket
   logical, intent(in) :: ket

   real(wp), allocatable :: cart1(:, :), cart2(:, :)
   real(wp), allocatable :: sphr1(:, :), sphr2(:, :)
   real(wp), allocatable :: cart3(:, :, :), sphr3(:, :, :)
   real(wp), allocatable :: cart4(:, :, :, :), sphr4(:, :, :, :)
   real(wp) :: prod_sphr, prod_cart
   integer :: ncj, nci, nsj, nsi, nrow_cart, ncol_cart

   ! Number of Cartesian functions for bra and ket
   ncj = (lj + 1) * (lj + 2) / 2
   nci = (li + 1) * (li + 2) / 2

   ! Number of spherical functions for bra and ket
   nsj = 2*lj + 1
   nsi = 2*li + 1

   if (bra) then
      nrow_cart = ncj
   else
      nrow_cart = nsj
   end if

   if (ket) then
      ncol_cart = nci
   else
      ncol_cart = nsi
   end if

   allocate(cart1(nrow_cart, ncol_cart), cart2(nrow_cart, ncol_cart))
   allocate(sphr1(nsj, nsi), sphr2(nsj, nsi))

   ! Random coefficients around 0.0
   call random_number(cart1)
   call random_number(sphr1)
   cart1 = cart1 - 0.5_wp
   sphr1 = sphr1 - 0.5_wp

   ! Transformation from cartesian to spherical
   call transform0(lj, li, cart1, sphr2, bra, ket)

   ! Adjoint transformation from spherical to cartesian
   call adjoint_transform0(lj, li, sphr1, cart2, bra, ket)

   ! Hadamard product in spherical and cartesian basis
   prod_sphr = sum(sphr1 * sphr2)
   prod_cart = sum(cart1 * cart2)

   call check(error, prod_sphr, prod_cart, thr=thr)
   if (allocated(error)) return

   ! Check wrapper for three dimensional arrays with leading batch dimension
   allocate(sphr3(1, nsj, nsi), source=0.0_wp)
   allocate(cart3(1, nrow_cart, ncol_cart), source=0.0_wp)
   sphr3(1, :, :) = sphr1

   call adjoint_transform1(lj, li, sphr3, cart3, bra, ket)
   prod_cart = sum(cart1 * cart3(1, :, :))

   call check(error, prod_sphr, prod_cart, thr=thr)
   if (allocated(error)) return

   ! Check wrapper for four dimensional arrays with leading batch dimensions
   allocate(sphr4(1, 1, nsj, nsi), source=0.0_wp)
   allocate(cart4(1, 1, nrow_cart, ncol_cart), source=0.0_wp)
   sphr4(1, 1, :, :) = sphr1

   call adjoint_transform2(lj, li, sphr4, cart4, bra, ket)
   prod_cart = sum(cart1 * cart4(1, 1, :, :))

   call check(error, prod_sphr, prod_cart, thr=thr)
   if (allocated(error)) return

end subroutine test_trafo_adjoint


subroutine test_trafo_ss(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 0, 0)

end subroutine test_trafo_ss

subroutine test_trafo_sp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 0, 1)
   call test_trafo_pair(error, 1, 0)

end subroutine test_trafo_sp

subroutine test_trafo_sd(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 0, 2)
   call test_trafo_pair(error, 2, 0)

end subroutine test_trafo_sd

subroutine test_trafo_sf(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 0, 3)
   call test_trafo_pair(error, 3, 0)

end subroutine test_trafo_sf

subroutine test_trafo_sg(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 0, 4)
   call test_trafo_pair(error, 4, 0)

end subroutine test_trafo_sg

subroutine test_trafo_pp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 1, 1)

end subroutine test_trafo_pp

subroutine test_trafo_pd(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 1, 2)
   call test_trafo_pair(error, 2, 1)

end subroutine test_trafo_pd

subroutine test_trafo_pf(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 1, 3)
   call test_trafo_pair(error, 3, 1)

end subroutine test_trafo_pf

subroutine test_trafo_pg(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 1, 4)
   call test_trafo_pair(error, 4, 1)

end subroutine test_trafo_pg

subroutine test_trafo_dd(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 2, 2)

end subroutine test_trafo_dd

subroutine test_trafo_df(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 2, 3)
   call test_trafo_pair(error, 3, 2)

end subroutine test_trafo_df

subroutine test_trafo_dg(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 2, 4)
   call test_trafo_pair(error, 4, 2)

end subroutine test_trafo_dg

subroutine test_trafo_ff(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 3, 3)

end subroutine test_trafo_ff

subroutine test_trafo_fg(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 3, 4)
   call test_trafo_pair(error, 4, 3)

end subroutine test_trafo_fg

subroutine test_trafo_gg(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_trafo_pair(error, 4, 4)

end subroutine test_trafo_gg


subroutine test_density_trafo(mol, calc, wfn, ref, error, thr_in)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Calculator instance
   type(xtb_calculator), intent(in) :: calc
   !> Wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Reference density matrix in Cartesian basis
   real(wp), intent(in) :: ref(:, :, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: spin
   real(wp) :: energy, thr_
   real(wp), allocatable :: coeff_cart(:, :, :), focc(:), density_cart(:, :, :)
   type(context_type) :: ctx

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   ! Perform single point calculation
   call eeq_guess(mol, calc, wfn, error)
   if (allocated(error)) return
   energy = 0.0_wp
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)

   ! Transform MO coefficients from spherical to cartesian basis
   allocate(coeff_cart(calc%bas%nao_cart, calc%bas%nao, wfn%nspin), source=0.0_wp)
   do spin = 1, wfn%nspin
      call calc%bas%spherical_to_cartesian_trafo(mol, wfn%coeff(:, :, spin), &
         & coeff_cart(:, :, spin))
   end do

   ! Calculate density matrix in the cartesian basis
   ! since the coefficients are sensitive to the phase
   allocate(focc(calc%bas%nao))
   allocate(density_cart(calc%bas%nao_cart, calc%bas%nao_cart, wfn%nspin))
   do spin = 1, wfn%nspin
      if (wfn%nspin == 1) then
         focc = wfn%focc(:, 1) + wfn%focc(:, 2)
      else
         focc = wfn%focc(:, spin)
      end if

      call gemm(coeff_cart(:, :, spin) * spread(focc, 1, calc%bas%nao_cart), &
         & transpose(coeff_cart(:, :, spin)), density_cart(:, :, spin))
   end do

   if (any(abs(density_cart - ref) > thr_)) then
      call test_failed(error, "Cartesian density matrix does not match")
      print"(3es21.14)", density_cart
      print'("---")'
      print"(3es21.14)", ref
      print'("---")'
      print"(3es21.14)", density_cart-ref
   end if

end subroutine test_density_trafo

subroutine test_trafo_pcl(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: density_cart(20, 20, 1) = reshape([&
      &  1.96022219468921E+00_wp, -4.31628878279072E-16_wp, -3.69459754888647E-01_wp, &
      & -1.11098031904690E-16_wp,  4.71236654550718E-02_wp,  4.71236654550742E-02_wp, &
      & -9.42473309101460E-02_wp,  5.62965590903630E-16_wp, -4.06828497062294E-16_wp, &
      & -8.52830081768616E-16_wp, -1.89854934222313E-01_wp,  8.09624498998321E-16_wp, &
      &  7.25134761302635E-02_wp,  1.59916712486622E-16_wp, -2.34268009256941E-02_wp, &
      & -2.34268009256942E-02_wp,  4.68536018513883E-02_wp,  3.48317988599535E-16_wp, &
      &  1.88088687294335E-16_wp, -1.76463194084002E-16_wp, -4.31628878279072E-16_wp, &
      &  1.07791095451279E+00_wp, -8.31505858137497E-16_wp,  7.22200077518664E-14_wp, &
      & -6.20379008988276E-17_wp, -3.24781092909577E-16_wp,  3.86818993808405E-16_wp, &
      &  1.07556891827876E-16_wp, -8.61810622865278E-15_wp, -5.04737692105317E-02_wp, &
      &  1.13747411177022E-16_wp,  1.91867623619473E-01_wp, -1.02920847827749E-15_wp, &
      & -2.68951527715444E-14_wp, -9.11978402679842E-17_wp,  6.66905072726581E-17_wp, &
      &  2.45073329953261E-17_wp, -2.29032624584476E-16_wp, -4.62824223390612E-15_wp, &
      & -7.61574005556389E-02_wp, -3.69459754888647E-01_wp, -8.31505858137497E-16_wp, &
      &  5.06882651413875E-01_wp,  6.61199693299378E-16_wp,  5.03840047492308E-03_wp, &
      &  5.03840047492275E-03_wp, -1.00768009498458E-02_wp, -1.14256085434253E-16_wp, &
      &  3.31421482075486E-16_wp,  1.82617289715317E-16_wp, -1.39976031367795E-02_wp, &
      &  2.04803054630828E-18_wp, -7.33326590175914E-01_wp, -4.24957452481721E-16_wp, &
      & -1.19309309178673E-03_wp, -1.19309309178663E-03_wp,  2.38618618357336E-03_wp, &
      &  3.14132994583873E-18_wp, -1.88721574223030E-16_wp, -1.26943910123494E-17_wp, &
      & -1.11098031904690E-16_wp,  7.22200077518664E-14_wp,  6.61199693299378E-16_wp, &
      &  1.07791095451263E+00_wp,  1.42412783202672E-16_wp,  1.00685168645336E-16_wp, &
      & -2.43097951848008E-16_wp, -2.48680357460767E-16_wp, -5.04737692105127E-02_wp, &
      & -8.43769498715119E-15_wp, -4.77025776053928E-17_wp, -2.65620858641569E-14_wp, &
      &  5.25687101242358E-16_wp,  1.91867623619530E-01_wp, -1.45451846493814E-16_wp, &
      &  5.70778474664071E-17_wp,  8.83739990274073E-17_wp,  2.92526454109793E-16_wp, &
      & -7.61574005556284E-02_wp, -4.41660596983695E-15_wp,  4.71236654550718E-02_wp, &
      & -6.20379008988276E-17_wp,  5.03840047492308E-03_wp,  1.42412783202672E-16_wp, &
      &  1.97110953674444E-03_wp,  1.97110953674450E-03_wp, -3.94221907348894E-03_wp, &
      &  1.32961974553750E-17_wp, -1.55896748336360E-17_wp,  6.40735691835751E-18_wp, &
      & -3.38948815549140E-02_wp,  2.16655877069712E-16_wp, -2.24138588507156E-02_wp, &
      & -5.72903643784898E-17_wp, -8.09838349982601E-04_wp, -8.09838349982605E-04_wp, &
      &  1.61967669996521E-03_wp,  9.02129887592915E-18_wp, -8.57688494343901E-18_wp, &
      & -1.97876393374441E-17_wp,  4.71236654550742E-02_wp, -3.24781092909577E-16_wp, &
      &  5.03840047492275E-03_wp,  1.00685168645336E-16_wp,  1.97110953674450E-03_wp, &
      &  1.97110953674456E-03_wp, -3.94221907348906E-03_wp,  1.32961974553757E-17_wp, &
      & -2.40519587348739E-17_wp,  1.25681864548440E-17_wp, -3.38948815549142E-02_wp, &
      &  1.24259581954180E-16_wp, -2.24138588507157E-02_wp, -1.42095173197775E-16_wp, &
      & -8.09838349982632E-04_wp, -8.09838349982636E-04_wp,  1.61967669996527E-03_wp, &
      &  9.02129887592965E-18_wp, -4.90329933150083E-18_wp, -7.96338946182486E-19_wp, &
      & -9.42473309101460E-02_wp,  3.86818993808405E-16_wp, -1.00768009498458E-02_wp, &
      & -2.43097951848008E-16_wp, -3.94221907348894E-03_wp, -3.94221907348906E-03_wp, &
      &  7.88443814697800E-03_wp, -2.65923949107507E-17_wp,  3.96416335685099E-17_wp, &
      & -1.89755433732015E-17_wp,  6.77897631098282E-02_wp, -3.40915459023892E-16_wp, &
      &  4.48277177014313E-02_wp,  1.99385537576265E-16_wp,  1.61967669996523E-03_wp, &
      &  1.61967669996524E-03_wp, -3.23935339993047E-03_wp, -1.80425977518588E-17_wp, &
      &  1.34801842749398E-17_wp,  2.05839782836266E-17_wp,  5.62965590903630E-16_wp, &
      &  1.07556891827876E-16_wp, -1.14256085434253E-16_wp, -2.48680357460767E-16_wp, &
      &  1.32961974553750E-17_wp,  1.32961974553757E-17_wp, -2.65923949107507E-17_wp, &
      &  3.03868259254195E-31_wp,  2.95679184093213E-17_wp,  4.00370232499131E-17_wp, &
      & -5.51399477298318E-17_wp,  3.53975478135562E-16_wp,  3.41688307716610E-17_wp, &
      &  8.88793920882234E-17_wp, -6.62731219503248E-18_wp, -6.62731219503245E-18_wp, &
      &  1.32546243900649E-17_wp, -1.06242698042029E-31_wp,  1.63217263204427E-17_wp, &
      & -1.07382345176700E-17_wp, -4.06828497062294E-16_wp, -8.61810622865278E-15_wp, &
      &  3.31421482075486E-16_wp, -5.04737692105127E-02_wp, -1.55896748336360E-17_wp, &
      & -2.40519587348739E-17_wp,  3.96416335685099E-17_wp,  2.95679184093213E-17_wp, &
      &  3.41889525777604E-02_wp,  1.17007098454636E-15_wp, -9.57915379548323E-17_wp, &
      &  3.47118167542959E-15_wp, -6.67167431068387E-16_wp,  2.27433048909729E-01_wp, &
      & -2.28934973120123E-18_wp,  3.10411422805550E-17_wp, -2.87517925493538E-17_wp, &
      &  1.03182578415653E-17_wp,  1.34969454111339E-03_wp,  6.11490025281825E-16_wp, &
      & -8.52830081768616E-16_wp, -5.04737692105317E-02_wp,  1.82617289715317E-16_wp, &
      & -8.43769498715119E-15_wp,  6.40735691835751E-18_wp,  1.25681864548440E-17_wp, &
      & -1.89755433732015E-17_wp,  4.00370232499131E-17_wp,  1.17007098454636E-15_wp, &
      &  3.41889525777625E-02_wp,  1.92950814569396E-16_wp,  2.27433048909735E-01_wp, &
      &  1.23395282461191E-16_wp,  4.02802791121815E-15_wp, -2.27416693403459E-17_wp, &
      &  3.25960150188036E-17_wp, -9.85434567845777E-18_wp, -7.97202226222352E-17_wp, &
      &  5.55545193181572E-16_wp,  1.34969454111463E-03_wp, -1.89854934222313E-01_wp, &
      &  1.13747411177022E-16_wp, -1.39976031367795E-02_wp, -4.77025776053928E-17_wp, &
      & -3.38948815549140E-02_wp, -3.38948815549142E-02_wp,  6.77897631098282E-02_wp, &
      & -5.51399477298318E-17_wp, -9.57915379548323E-17_wp,  1.92950814569396E-16_wp, &
      &  1.97254097945181E+00_wp,  6.32122311057884E-17_wp,  1.62414000110393E-01_wp, &
      &  2.18292545659950E-18_wp,  7.69029197968299E-03_wp,  7.69029197968316E-03_wp, &
      & -1.53805839593662E-02_wp,  6.67419266931980E-17_wp,  3.54773941002931E-18_wp, &
      &  9.94892415398064E-16_wp,  8.09624498998321E-16_wp,  1.91867623619473E-01_wp, &
      &  2.04803054630830E-18_wp, -2.65620858641569E-14_wp,  2.16655877069712E-16_wp, &
      &  1.24259581954180E-16_wp, -3.40915459023892E-16_wp,  3.53975478135562E-16_wp, &
      &  3.47118167542959E-15_wp,  2.27433048909735E-01_wp,  6.32122311057884E-17_wp, &
      &  1.79039127284740E+00_wp, -1.58698189670615E-16_wp,  9.72139035937403E-15_wp, &
      & -3.12677515642898E-16_wp,  1.81426312876657E-16_wp,  1.31251202766241E-16_wp, &
      & -7.12641717691887E-16_wp,  1.87697080100691E-15_wp, -3.00207544958856E-02_wp, &
      &  7.25134761302635E-02_wp, -1.02920847827749E-15_wp, -7.33326590175914E-01_wp, &
      &  5.25687101242358E-16_wp, -2.24138588507156E-02_wp, -2.24138588507157E-02_wp, &
      &  4.48277177014313E-02_wp,  3.41688307716610E-17_wp, -6.67167431068387E-16_wp, &
      &  1.23395282461191E-16_wp,  1.62414000110393E-01_wp, -1.58698189670615E-16_wp, &
      &  1.19108922406209E+00_wp, -2.07339616780469E-16_wp,  8.57919007600030E-03_wp, &
      &  8.57919007600018E-03_wp, -1.71583801520005E-02_wp, -9.54742728903731E-17_wp, &
      &  1.58300091699251E-16_wp,  2.85677087660719E-16_wp,  1.59916712486622E-16_wp, &
      & -2.68951527715444E-14_wp, -4.24957452481721E-16_wp,  1.91867623619530E-01_wp, &
      & -5.72903643784898E-17_wp, -1.42095173197775E-16_wp,  1.99385537576265E-16_wp, &
      &  8.88793920882234E-17_wp,  2.27433048909729E-01_wp,  4.02629318774217E-15_wp, &
      &  2.18292545659951E-18_wp,  9.72139035937403E-15_wp, -2.07339616780469E-16_wp, &
      &  1.79039127284737E+00_wp, -9.67122129575804E-17_wp,  2.57384468356189E-16_wp, &
      & -1.60672255398608E-16_wp,  2.30473471113603E-16_wp, -3.00207544958891E-02_wp, &
      &  2.29937596740726E-15_wp, -2.34268009256941E-02_wp, -9.11978402679842E-17_wp, &
      & -1.19309309178673E-03_wp, -1.45451846493814E-16_wp, -8.09838349982601E-04_wp, &
      & -8.09838349982632E-04_wp,  1.61967669996523E-03_wp, -6.62731219503248E-18_wp, &
      & -2.28934973120123E-18_wp, -2.27416693403459E-17_wp,  7.69029197968299E-03_wp, &
      & -3.12677515642898E-16_wp,  8.57919007600030E-03_wp, -9.67122129575804E-17_wp, &
      &  3.63656149333574E-04_wp,  3.63656149333575E-04_wp, -7.27312298667149E-04_wp, &
      & -4.77930967387482E-18_wp,  1.02256863982123E-17_wp,  1.53607452191457E-17_wp, &
      & -2.34268009256942E-02_wp,  6.66905072726581E-17_wp, -1.19309309178663E-03_wp, &
      &  5.70778474664071E-17_wp, -8.09838349982605E-04_wp, -8.09838349982636E-04_wp, &
      &  1.61967669996524E-03_wp, -6.62731219503245E-18_wp,  3.10411422805550E-17_wp, &
      &  3.25960150188036E-17_wp,  7.69029197968316E-03_wp,  1.81426312876657E-16_wp, &
      &  8.57919007600018E-03_wp,  2.57384468356189E-16_wp,  3.63656149333575E-04_wp, &
      &  3.63656149333576E-04_wp, -7.27312298667151E-04_wp, -4.77930967387496E-18_wp, &
      & -7.06529069837249E-18_wp, -1.63261802835476E-19_wp,  4.68536018513883E-02_wp, &
      &  2.45073329953261E-17_wp,  2.38618618357336E-03_wp,  8.83739990274073E-17_wp, &
      &  1.61967669996521E-03_wp,  1.61967669996527E-03_wp, -3.23935339993047E-03_wp, &
      &  1.32546243900649E-17_wp, -2.87517925493538E-17_wp, -9.85434567845777E-18_wp, &
      & -1.53805839593662E-02_wp,  1.31251202766241E-16_wp, -1.71583801520005E-02_wp, &
      & -1.60672255398608E-16_wp, -7.27312298667149E-04_wp, -7.27312298667151E-04_wp, &
      &  1.45462459733430E-03_wp,  9.55861934774978E-18_wp, -3.16039569983986E-18_wp, &
      & -1.51974834163102E-17_wp,  3.48317988599535E-16_wp, -2.29032624584476E-16_wp, &
      &  3.14132994583873E-18_wp,  2.92526454109793E-16_wp,  9.02129887592915E-18_wp, &
      &  9.02129887592965E-18_wp, -1.80425977518588E-17_wp, -1.06242698042029E-31_wp, &
      &  1.03182578415653E-17_wp, -7.97202226222352E-17_wp,  6.67419266931980E-17_wp, &
      & -7.12641717691887E-16_wp, -9.54742728903731E-17_wp,  2.30473471113603E-16_wp, &
      & -4.77930967387482E-18_wp, -4.77930967387496E-18_wp,  9.55861934774978E-18_wp, &
      &  4.81946035529280E-31_wp, -2.23403459632380E-17_wp,  2.24806227067167E-17_wp, &
      &  1.88088687294335E-16_wp, -4.62824223390612E-15_wp, -1.88721574223030E-16_wp, &
      & -7.61574005556284E-02_wp, -8.57688494343901E-18_wp, -4.90329933150083E-18_wp, &
      &  1.34801842749398E-17_wp,  1.63217263204427E-17_wp,  1.34969454111339E-03_wp, &
      &  5.55545193181572E-16_wp,  3.54773941002931E-18_wp,  1.87697080100691E-15_wp, &
      &  1.58300091699251E-16_wp, -3.00207544958891E-02_wp,  1.02256863982123E-17_wp, &
      & -7.06529069837249E-18_wp, -3.16039569983986E-18_wp, -2.23403459632380E-17_wp, &
      &  5.53508929793042E-03_wp,  2.74953670942324E-16_wp, -1.76463194084002E-16_wp, &
      & -7.61574005556389E-02_wp, -1.26943910123494E-17_wp, -4.42007541678890E-15_wp, &
      & -1.97876393374441E-17_wp, -7.96338946182485E-19_wp,  2.05839782836266E-17_wp, &
      & -1.07382345176700E-17_wp,  6.11490025281825E-16_wp,  1.34969454111463E-03_wp, &
      &  9.94892415398064E-16_wp, -3.00207544958856E-02_wp,  2.85677087660719E-16_wp, &
      &  2.29937596740726E-15_wp,  1.53607452191457E-17_wp, -1.63261802835476E-19_wp, &
      & -1.51974834163102E-17_wp,  2.24806227067167E-17_wp,  2.74953670942324E-16_wp, &
      &  5.53508929793114E-03_wp], shape(density_cart))

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn

   call get_structure(mol, "MB16-43", "PCl")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   call test_density_trafo(mol, calc, wfn, density_cart, error, thr_in=thr*100.0_wp)

end subroutine test_trafo_pcl



end module test_integral_trafo
