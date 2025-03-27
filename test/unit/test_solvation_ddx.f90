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

module test_solvation_ddx
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_container, only : container_cache
   use tblite_scf_potential, only : potential_type
   use tblite_solvation_ddx, only : ddx_solvation, ddx_solvation_model, ddx_input 
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction, eeq_guess
   use tblite_basis_type, only : basis_type
   use tblite_integral_type, only : integral_type
   use tblite_context_type, only : context_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint

   implicit none
   private

   public :: collect_solvation_ddx

   real(wp), parameter :: acc = 0.01_wp
   real(wp), parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = 10*sqrt(epsilon(1.0_wp))
   real(wp), parameter :: kt = 300.0_wp * 3.166808578545117e-06_wp


contains


!> Collect all exported unit tests
subroutine collect_solvation_ddx(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("energy-mol-cosmo", test_e_cosmo_m01), &
      new_unittest("energy-mol-pcm", test_e_pcm_m01), &
      ! new_unittest("energy-mol-lpb", test_e_lpb_m01), &
      new_unittest("gradient-mol-num-cosmo", test_g_num_cosmo_m02), &
      ! new_unittest("gradient-mol-cosmo", test_g_cosmo_m02), &
      new_unittest("gradient-mol-num-pcm", test_g_pcm_m02), &
      ! new_unittest("gradient-mol-lpb", test_g_lpb_m02), &
      new_unittest("potential-mol-cosmo", test_p_cosmo_m03), &
      new_unittest("potential-mol-pcm", test_p_pcm_m03) &
      ! new_unittest("potential-mol-lpb", test_p_lpb_m03) &
      ]

end subroutine collect_solvation_ddx


subroutine test_e(error, model, mol, qat, ref, kappa)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Solvation model (COSMO=11, CPCM=12, PCM=2, LPB=3)
   integer, intent(in) :: model

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Reference energy
   real(wp), intent(in) :: ref

   !> Debye-Hückel screening parameter (only used in LPB)
   real(wp), optional, intent(in) :: kappa

   type(ddx_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: feps = 80.0_wp, rscale = 1.0_wp
   integer, parameter :: nang = 110
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp) :: energy(mol%nat)

   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1), source=0.0_wp)
   energy = 0.0_wp

   if (present(kappa)) then
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale, kappa=kappa))
   else
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale))
   end if

   call solv%update(mol, cache)

   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)

   if (abs(sum(energy) - ref) > thr) then
      call test_failed(error, "Energy does not match reference")
      print *, sum(energy)
   end if
end subroutine test_e


subroutine test_g_num(error, model, mol, qat, kappa)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Solvation model (COSMO=11, CPCM=12, PCM=2, LPB=3)
   integer, intent(in) :: model

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Debye-Hückel screening parameter (only used in LPB)
   real(wp), optional, intent(in) :: kappa

   type(ddx_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: feps = 80.0_wp, rscale = 1.0_wp
   integer, parameter :: nang = 26
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp), allocatable :: gradient(:, :), numg(:, :)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3)
   integer :: ii, ic

   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1))

   if (present(kappa)) then
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale, kappa=kappa))
   else
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale))
   end if

   allocate(numg(3, mol%nat), gradient(3, mol%nat))
   do ii = 1, mol%nat
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(ic, ii) = mol%xyz(ic, ii) + step
         call solv%update(mol, cache)
         call solv%get_potential(mol, cache, wfn, pot)
         call solv%get_energy(mol, cache, wfn, er)

         mol%xyz(ic, ii) = mol%xyz(ic, ii) - 2*step
         call solv%update(mol, cache)
         call solv%get_potential(mol, cache, wfn, pot)
         call solv%get_energy(mol, cache, wfn, el)

         mol%xyz(ic, ii) = mol%xyz(ic, ii) + step
         numg(ic, ii) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp

   call solv%update(mol, cache)
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
end subroutine test_g_num

subroutine test_g(error, model, mol, qat, ref, kappa)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Solvation model (COSMO=11, CPCM=12, PCM=2, LPB=3)
   integer, intent(in) :: model

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Reference gradient
   real(wp), intent(in) :: ref(:, :)

   !> Debye-Hückel screening parameter (only used in LPB)
   real(wp), optional, intent(in) :: kappa

   type(ddx_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: feps = 80.0_wp, rscale = 1.0_wp
   integer, parameter :: nang = 26
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp), allocatable :: gradient(:, :), numg(:, :)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3)
   integer :: ii, ic

   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1))

   if (present(kappa)) then
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale, kappa=kappa))
   else
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale))
   end if

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp

   call solv%update(mol, cache)
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)
   call solv%get_gradient(mol, cache, wfn, gradient, sigma)

   if (any(abs(gradient - ref) > thr)) then
      call test_failed(error, "Gradient does not match")
      print '(3es20.13)', gradient
      print '(a)', "---"
      print '(3es20.13)', ref
      print '(a)', "---"
      print '(3es20.13)', gradient - ref
   end if
end subroutine test_g


subroutine test_p(error, model, mol, qat, kappa)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Solvation model (COSMO=11, CPCM=12, PCM=2, LPB=3)
   integer, intent(in) :: model

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Debye-Hückel screening parameter (only used in LPB)
   real(wp), optional, intent(in) :: kappa

   type(ddx_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache), allocatable :: cache
   real(wp), parameter :: feps = 80.0_wp, rscale = 1.0_wp
   integer, parameter :: nang = 26
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), parameter :: thr = 1e+3_wp*sqrt(epsilon(1.0_wp))
   real(wp), allocatable :: vat(:)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat)
   integer :: ii

   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1))

   if (present(kappa)) then
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale, kappa=kappa))
   else
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale))
   end if

   allocate(cache)
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
      print '(a)', 'analytical'
      print '(3es20.13)', pot%vat
      print '(a)', "---"
      print '(a)', 'numerical'
      print '(3es20.13)', vat
      print '(a)', "---"
      print '(a)', 'diff'
      print '(3es20.13)', [pot%vat] - vat
   end if
end subroutine test_p


subroutine test_e_cosmo_m01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 7.73347900345264E-1_wp, 1.07626888948184E-1_wp,-3.66999593831010E-1_wp,&
      & 4.92833325937897E-2_wp,-1.83332156197733E-1_wp, 2.33302086605469E-1_wp,&
      & 6.61837152062315E-2_wp,-5.43944165050002E-1_wp,-2.70264356583716E-1_wp,&
      & 2.66618968841682E-1_wp, 2.62725033202480E-1_wp,-7.15315510172571E-2_wp,&
      &-3.73300777019193E-1_wp, 3.84585237785621E-2_wp,-5.05851088366940E-1_wp,&
      & 5.17677238544189E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_e(error, ddx_solvation_model%cosmo, mol, qat, -3.4697720884118800E-2_wp)

end subroutine test_e_cosmo_m01

subroutine test_e_pcm_m01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 7.73347900345264E-1_wp, 1.07626888948184E-1_wp,-3.66999593831010E-1_wp,&
      & 4.92833325937897E-2_wp,-1.83332156197733E-1_wp, 2.33302086605469E-1_wp,&
      & 6.61837152062315E-2_wp,-5.43944165050002E-1_wp,-2.70264356583716E-1_wp,&
      & 2.66618968841682E-1_wp, 2.62725033202480E-1_wp,-7.15315510172571E-2_wp,&
      &-3.73300777019193E-1_wp, 3.84585237785621E-2_wp,-5.05851088366940E-1_wp,&
      & 5.17677238544189E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_e(error, ddx_solvation_model%pcm, mol, qat, -3.3624259293951506E-2_wp)

end subroutine test_e_pcm_m01

subroutine test_e_lpb_m01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 7.73347900345264E-1_wp, 1.07626888948184E-1_wp,-3.66999593831010E-1_wp,&
      & 4.92833325937897E-2_wp,-1.83332156197733E-1_wp, 2.33302086605469E-1_wp,&
      & 6.61837152062315E-2_wp,-5.43944165050002E-1_wp,-2.70264356583716E-1_wp,&
      & 2.66618968841682E-1_wp, 2.62725033202480E-1_wp,-7.15315510172571E-2_wp,&
      &-3.73300777019193E-1_wp, 3.84585237785621E-2_wp,-5.05851088366940E-1_wp,&
      & 5.17677238544189E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_e(error, ddx_solvation_model%lpb, mol, qat, -3.1747235888764228E-2_wp, kappa=0.5_wp)

end subroutine test_e_lpb_m01


subroutine test_g_num_cosmo_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 7.38394711236234E-2_wp,-1.68354976558608E-1_wp,-3.47642833746823E-1_wp,&
      &-7.05489267186003E-1_wp, 7.73548301641266E-1_wp, 2.30207581365386E-1_wp,&
      & 1.02748501676354E-1_wp, 9.47818107467040E-2_wp, 2.44260351729187E-2_wp,&
      & 2.34984927037408E-1_wp,-3.17839896393030E-1_wp, 6.67112994818879E-1_wp,&
      &-4.78119977010488E-1_wp, 6.57536027459275E-2_wp, 1.08259054549882E-1_wp,&
      &-3.58215329983396E-1_wp]

   call get_structure(mol, "MB16-43", "02")
   call test_g_num(error, ddx_solvation_model%cosmo, mol, qat)

end subroutine test_g_num_cosmo_m02

subroutine test_g_cosmo_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 7.38394711236234E-2_wp,-1.68354976558608E-1_wp,-3.47642833746823E-1_wp,&
      &-7.05489267186003E-1_wp, 7.73548301641266E-1_wp, 2.30207581365386E-1_wp,&
      & 1.02748501676354E-1_wp, 9.47818107467040E-2_wp, 2.44260351729187E-2_wp,&
      & 2.34984927037408E-1_wp,-3.17839896393030E-1_wp, 6.67112994818879E-1_wp,&
      &-4.78119977010488E-1_wp, 6.57536027459275E-2_wp, 1.08259054549882E-1_wp,&
      &-3.58215329983396E-1_wp]

   real(wp), parameter :: ref(3,16) = reshape([&
      & -8.72885014e-04_wp,  1.29722126e-03_wp,  4.82135743e-04_wp,  2.86976918e-03_wp, &
      &  4.91109654e-05_wp, -4.29260562e-04_wp, -5.67289804e-03_wp,  8.14714673e-04_wp, &
      & -4.12468023e-04_wp, -1.64554349e-04_wp,  1.83163634e-03_wp,  8.14645705e-04_wp, &
      & -2.92212009e-03_wp, -1.29547615e-03_wp,  1.77475155e-03_wp,  1.83559208e-03_wp, &
      & -2.82980413e-04_wp,  1.71557120e-03_wp,  2.73353837e-03_wp, -2.05688333e-03_wp, &
      &  6.38494126e-03_wp,  5.78537370e-04_wp, -2.37129787e-03_wp,  2.90753892e-04_wp, &
      &  5.08752936e-04_wp,  2.15648519e-04_wp, -1.25921100e-03_wp,  1.01160620e-03_wp, &
      & -4.55980356e-03_wp,  1.90794668e-03_wp, -5.63872868e-03_wp,  8.21396325e-04_wp, &
      &  5.02731375e-03_wp,  4.01930567e-03_wp, -1.74485028e-03_wp, -1.42778390e-02_wp, &
      & -6.01680881e-03_wp,  2.69775116e-04_wp,  1.42002404e-02_wp, -7.83354901e-04_wp, &
      & -1.47082304e-03_wp,  5.63906422e-04_wp,  6.02706612e-03_wp,  5.75784078e-04_wp, &
      & -7.23463064e-04_wp, -1.65844186e-04_wp, -5.66640254e-03_wp,  1.66094445e-04_wp &
      ], shape=[3,16])

   call get_structure(mol, "MB16-43", "02")
   call test_g(error, ddx_solvation_model%cosmo, mol, qat, ref)

end subroutine test_g_cosmo_m02

subroutine test_g_pcm_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 7.38394711236234E-2_wp,-1.68354976558608E-1_wp,-3.47642833746823E-1_wp,&
      &-7.05489267186003E-1_wp, 7.73548301641266E-1_wp, 2.30207581365386E-1_wp,&
      & 1.02748501676354E-1_wp, 9.47818107467040E-2_wp, 2.44260351729187E-2_wp,&
      & 2.34984927037408E-1_wp,-3.17839896393030E-1_wp, 6.67112994818879E-1_wp,&
      &-4.78119977010488E-1_wp, 6.57536027459275E-2_wp, 1.08259054549882E-1_wp,&
      &-3.58215329983396E-1_wp]

   call get_structure(mol, "MB16-43", "02")
   call test_g_num(error, ddx_solvation_model%pcm, mol, qat)

end subroutine test_g_pcm_m02

subroutine test_g_lpb_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 7.38394711236234E-2_wp,-1.68354976558608E-1_wp,-3.47642833746823E-1_wp,&
      &-7.05489267186003E-1_wp, 7.73548301641266E-1_wp, 2.30207581365386E-1_wp,&
      & 1.02748501676354E-1_wp, 9.47818107467040E-2_wp, 2.44260351729187E-2_wp,&
      & 2.34984927037408E-1_wp,-3.17839896393030E-1_wp, 6.67112994818879E-1_wp,&
      &-4.78119977010488E-1_wp, 6.57536027459275E-2_wp, 1.08259054549882E-1_wp,&
      &-3.58215329983396E-1_wp]

   call get_structure(mol, "MB16-43", "02")
   call test_g_num(error, ddx_solvation_model%lpb, mol, qat, kappa=0.5_wp)

end subroutine test_g_lpb_m02


subroutine test_p_cosmo_m03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-1.77788256288236E-1_wp,-8.22943267808161E-1_wp, 4.04578389873281E-2_wp,&
      & 5.79710531992282E-1_wp, 6.99601887637659E-1_wp, 6.84309612639107E-2_wp,&
      &-3.42971414989811E-1_wp, 4.64954031865410E-2_wp, 6.77012204116428E-2_wp,&
      & 8.49931225363225E-2_wp,-5.22285304699699E-1_wp,-2.92515001764712E-1_wp,&
      &-3.98375452377043E-1_wp, 2.09769668297792E-1_wp, 7.23140464830357E-1_wp,&
      & 3.65775987838250E-2_wp]

   call get_structure(mol, "MB16-43", "03")
   call test_p(error, ddx_solvation_model%cosmo, mol, qat)

end subroutine test_p_cosmo_m03

subroutine test_p_pcm_m03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-1.77788256288236E-1_wp,-8.22943267808161E-1_wp, 4.04578389873281E-2_wp,&
      & 5.79710531992282E-1_wp, 6.99601887637659E-1_wp, 6.84309612639107E-2_wp,&
      &-3.42971414989811E-1_wp, 4.64954031865410E-2_wp, 6.77012204116428E-2_wp,&
      & 8.49931225363225E-2_wp,-5.22285304699699E-1_wp,-2.92515001764712E-1_wp,&
      &-3.98375452377043E-1_wp, 2.09769668297792E-1_wp, 7.23140464830357E-1_wp,&
      & 3.65775987838250E-2_wp]

   call get_structure(mol, "MB16-43", "03")
   call test_p(error, ddx_solvation_model%pcm, mol, qat)

end subroutine test_p_pcm_m03

subroutine test_p_lpb_m03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-1.77788256288236E-1_wp,-8.22943267808161E-1_wp, 4.04578389873281E-2_wp,&
      & 5.79710531992282E-1_wp, 6.99601887637659E-1_wp, 6.84309612639107E-2_wp,&
      &-3.42971414989811E-1_wp, 4.64954031865410E-2_wp, 6.77012204116428E-2_wp,&
      & 8.49931225363225E-2_wp,-5.22285304699699E-1_wp,-2.92515001764712E-1_wp,&
      &-3.98375452377043E-1_wp, 2.09769668297792E-1_wp, 7.23140464830357E-1_wp,&
      & 3.65775987838250E-2_wp]

   call get_structure(mol, "MB16-43", "03")
   call test_p(error, ddx_solvation_model%lpb, mol, qat, kappa=0.5_wp)

end subroutine test_p_lpb_m03

end module test_solvation_ddx
