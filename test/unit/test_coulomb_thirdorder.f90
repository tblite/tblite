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

module test_coulomb_thirdorder
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_basis_type
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_container_cache, only : container_cache
   use tblite_coulomb_type, only : coulomb_type
   use tblite_coulomb_thirdorder, only : onsite_thirdorder, new_onsite_thirdorder
   use tblite_scf, only: new_potential, potential_type
   use tblite_ceh_ceh, only : get_effective_qat
   use tblite_ncoord_erf_en
   use tblite_ncoord_type, only : get_coordination_number
   use tblite_cutoff, only : get_lattice_points
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   implicit none
   private

   public :: collect_coulomb_thirdorder

   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   abstract interface
      subroutine coulomb_maker(coulomb, mol, shell)
         import :: coulomb_type, structure_type
         class(coulomb_type), allocatable, intent(out) :: coulomb
         type(structure_type), intent(in) :: mol
         logical, intent(in) :: shell
      end subroutine coulomb_maker
   end interface

   abstract interface
      subroutine charge_maker(wfn, mol, nshell)
         import :: wavefunction_type, structure_type
         class(wavefunction_type), intent(inout) :: wfn
         type(structure_type), intent(in) :: mol
         integer, optional, intent(in) :: nshell(:)
      end subroutine charge_maker
   end interface

contains


!> Collect all exported unit tests
subroutine collect_coulomb_thirdorder(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("energy-atom-gfn1", test_e_gfn1_m01), &
      new_unittest("energy-atom-gfn2", test_e_gfn2_m02), &
      new_unittest("energy-shell-gfn2", test_e_gfn2_m07), &
      new_unittest("energy-atom-pbc-gfn2", test_e_gfn2_oxacb), &
      new_unittest("energy-atom-sc-gfn2", test_e_gfn2_oxacb_sc), &
      new_unittest("gradient-atom-gfn1", test_g_gfn1_m03), &
      new_unittest("gradient-atom-gfn2", test_g_gfn2_m04), &
      new_unittest("gradient-shell-gfn2", test_g_gfn2_m08), &
      new_unittest("gradient-atom-pbc-gfn1", test_g_gfn1_urea), &
      new_unittest("sigma-shell-gfn2", test_s_gfn2_m09), &
      new_unittest("sigma-atom-pbc-gfn22", test_s_gfn2_pyrazine), &
      new_unittest("potential gradient-lih-effceh", test_ceh_potgrad_lih), &
      new_unittest("potential gradient-mb15-effceh", test_ceh_potgrad_m15), &
      new_unittest("potential gradient-mb16-shell-effceh", test_ceh_potgrad_m16), &
      new_unittest("potential-sigma-lih-effceh", test_ceh_potsigma_lih), &
      new_unittest("potential-sigma-co2-effceh", test_ceh_potsigma_co2), &
      new_unittest("potential-sigma-mb05-effceh", test_ceh_potsigma_mb05), &
      new_unittest("potential-sigma-mb17-shell-effceh", test_ceh_potsigma_mb17) &
      ]

end subroutine collect_coulomb_thirdorder


!> Factory to setup the CEH basis set for testing of the potential gradient
subroutine make_basis(bas, mol, ng)
   type(basis_type), intent(out) :: bas
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: ng

   integer, parameter :: nsh(20) = [&
   & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3]
   integer, parameter :: lsh(3, 20) = reshape([&
   & 0, 0, 0,  0, 0, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0, &
   & 0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 2, &
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, &
   & 0, 1, 0,  0, 1, 2], shape(lsh))

   integer, parameter :: pqn(3, 20) = reshape([&
   & 1, 0, 0,  1, 0, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0, &
   & 2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  3, 3, 0,  3, 3, 3, &
   & 3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3, &
   & 4, 4, 0,  4, 4, 3], shape(pqn))

   real(wp), parameter :: zeta(3, 20) = reshape([&
   & 1.23363166_wp, 0.00000000_wp, 0.00000000_wp, 2.27004605_wp, 0.00000000_wp, 0.00000000_wp, &
   & 0.86185456_wp, 1.42017184_wp, 0.00000000_wp, 1.76817995_wp, 1.44095844_wp, 0.00000000_wp, &
   & 2.06339837_wp, 1.52051807_wp, 0.00000000_wp, 2.56058582_wp, 1.86484737_wp, 0.00000000_wp, &
   & 2.71233631_wp, 2.19848968_wp, 0.00000000_wp, 3.21585650_wp, 2.41309737_wp, 0.00000000_wp, &
   & 3.82146807_wp, 2.63063636_wp, 0.00000000_wp, 4.62721228_wp, 2.53599954_wp, 0.00000000_wp, &
   & 0.93221172_wp, 1.55333839_wp, 0.00000000_wp, 1.77220557_wp, 1.59942632_wp, 2.98596647_wp, &
   & 2.26040231_wp, 1.78718151_wp, 2.00990188_wp, 1.85259089_wp, 1.81733349_wp, 1.65269988_wp, &
   & 2.65701241_wp, 2.03189759_wp, 2.03883661_wp, 2.60609998_wp, 2.16530440_wp, 2.41888232_wp, &
   & 2.78818934_wp, 2.24732894_wp, 1.99081182_wp, 2.55424399_wp, 2.20946190_wp, 1.93619550_wp, &
   & 1.73713827_wp, 1.33788617_wp, 0.00000000_wp, 2.47982574_wp, 1.07250770_wp, 2.11920764_wp],&
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


!> Factory to create onsite thirdorder objects based on GFN1-xTB values
subroutine make_coulomb_o1(coulomb, mol, shell)

   !> New coulomb object
   class(coulomb_type), allocatable, intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return a shell resolved object
   logical, intent(in) :: shell

   real(wp), parameter :: p_hubbard_derivs(20) = 0.1_wp * [&
      & 0.000000_wp, 1.500000_wp, 1.027370_wp, 0.900554_wp, 1.300000_wp, &
      & 1.053856_wp, 0.042507_wp,-0.005102_wp, 1.615037_wp, 1.600000_wp, &
      & 1.200000_wp, 1.100000_wp, 1.200000_wp, 1.500000_wp, 1.500000_wp, &
      & 1.500000_wp, 1.000000_wp, 0.829312_wp, 0.732923_wp, 1.116963_wp] 
   
   real(wp), allocatable :: hubbard_derivs(:, :)
   type(onsite_thirdorder), allocatable :: tmp

   allocate(tmp)
   hubbard_derivs = reshape(p_hubbard_derivs(mol%num), [1, mol%nid])
   call new_onsite_thirdorder(tmp, mol, hubbard_derivs)
   call move_alloc(tmp, coulomb)

end subroutine make_coulomb_o1

!> Factory to create onsite thirdorder objects based on GFN2-xTB values
subroutine make_coulomb_o2(coulomb, mol, shell)

   !> New coulomb object
   class(coulomb_type), allocatable, intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return a shell resolved object
   logical, intent(in) :: shell

   real(wp), parameter :: p_hubbard_derivs(20) = 0.1_wp * [&
      & 0.800000_wp, 2.000000_wp, 1.303821_wp, 0.574239_wp, 0.946104_wp, &
      & 1.500000_wp,-0.639780_wp,-0.517134_wp, 1.426212_wp, 0.500000_wp, &
      & 1.798727_wp, 2.349164_wp, 1.400000_wp, 1.936289_wp, 0.711291_wp, &
      &-0.501722_wp, 1.495483_wp,-0.315455_wp, 2.033085_wp, 2.006898_wp]
   real(wp), parameter :: shell_hubbard_derivs(0:4) = [1.0_wp, 0.5_wp, spread(0.25_wp, 1, 3)]
   integer, parameter :: shell_count(20) = [&
      & 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3]

   real(wp), allocatable :: hubbard_derivs(:, :)
   integer :: isp, izp, ish, il
   type(onsite_thirdorder), allocatable :: tmp

   allocate(tmp)
   if (shell) then
      allocate(hubbard_derivs(3, mol%nid))
      hubbard_derivs(:, :) = 0.0_wp
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, shell_count(izp)
            il = ish - 1
            hubbard_derivs(ish, isp) = p_hubbard_derivs(izp) * shell_hubbard_derivs(il)
         end do
      end do
      call new_onsite_thirdorder(tmp, mol, hubbard_derivs, shell_count(mol%num))
   else
      hubbard_derivs = reshape(p_hubbard_derivs(mol%num), [1, mol%nid])
      call new_onsite_thirdorder(tmp, mol, hubbard_derivs)
   end if
   call move_alloc(tmp, coulomb)

end subroutine make_coulomb_o2

!> Factory to create onsite thirdorder objects based on CEH values
subroutine make_coulomb_oceh(coulomb, mol, shell)

   !> New coulomb object
   class(coulomb_type), allocatable, intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return a shell resolved object
   logical, intent(in) :: shell

   real(wp), parameter :: p_hubbard_derivs(20) = [&
      &  0.8936213810_wp, -0.3936567743_wp, -0.7726174171_wp, -0.2849896764_wp, &
      &  0.0126634714_wp, -0.0082561791_wp,  0.0992949802_wp, -0.0267387652_wp, &
      & -0.0632999086_wp, -1.0106414497_wp, -0.3492075197_wp, -0.3191627473_wp, &
      &  0.0467483747_wp, -0.0920002125_wp, -0.0728788864_wp, -0.0213909690_wp, &
      & -0.0206065548_wp, -0.0432378706_wp, -0.0686554093_wp, -0.1671301006_wp] 
   integer, parameter :: shell_count(20) = [&
      & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3]

   real(wp), allocatable :: hubbard_derivs(:, :)
   integer :: isp, izp, ish, il
   type(onsite_thirdorder), allocatable :: tmp

   allocate(tmp)
   if (shell) then
      ! If shell-resolved, we use the atomic parameter for each shell
      allocate(hubbard_derivs(3, mol%nid))
      hubbard_derivs(:, :) = 0.0_wp
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, shell_count(izp)
            il = ish - 1
            hubbard_derivs(ish, isp) = p_hubbard_derivs(izp)
         end do
      end do
      call new_onsite_thirdorder(tmp, mol, hubbard_derivs, shell_count(mol%num))
   else
      hubbard_derivs = reshape(p_hubbard_derivs(mol%num), [1, mol%nid])
      call new_onsite_thirdorder(tmp, mol, hubbard_derivs)
   end if
   call move_alloc(tmp, coulomb)

end subroutine make_coulomb_oceh


!> Procedure to create CN based effective charges and gradients from CEH
subroutine get_charges_effceh(wfn, mol, nshell)

   !> New wavefunction object
   class(wavefunction_type), intent(inout) :: wfn

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return shell-resolved charges
   integer, optional, intent(in) :: nshell(:)

   real(wp), parameter :: ceh_cov_radii(20) = 0.5 * [&
   &  2.4040551903_wp,  1.8947380542_wp,  3.4227634078_wp,  3.5225408137_wp, &
   &  3.6150631704_wp,  2.8649682108_wp,  2.4695867541_wp,  2.3533691180_wp, &
   &  2.4992147462_wp,  3.3442607441_wp,  4.4665909451_wp,  4.3877250907_wp, &
   &  4.6647077385_wp,  4.2086223530_wp,  4.4750280107_wp,  4.2847281423_wp, &
   &  3.8560304959_wp,  3.9017061017_wp,  5.2392192639_wp,  5.1872031383_wp]

   real(wp), parameter :: pauling_en_ceh(20) = (1.0_wp/3.98_wp) * [ &
   &  1.9435211923_wp,  3.6116085622_wp,  2.4630915335_wp,  2.0658837656_wp, &
   &  2.3619778807_wp,  2.9484294262_wp,  3.8753937411_wp,  4.6235054741_wp, &
   &  3.9800000000_wp,  3.5124865506_wp,  2.3578254072_wp,  2.4225832022_wp, &
   &  2.1120078826_wp,  2.4607564741_wp,  2.7410779326_wp,  3.3517034720_wp, &
   &  4.1093492601_wp,  3.7979559518_wp,  2.4147937668_wp,  2.1974781961_wp]
   
   real(wp), allocatable :: lattr(:, :), cn_en(:), dcn_endr(:, :, :), dcn_endL(:, :, :)
   type(erf_en_ncoord_type) :: ncoord_en
   integer :: iat, ii, ish

   allocate(cn_en(mol%nat), dcn_endr(3, mol%nat, mol%nat), dcn_endL(3, 3, mol%nat))

   ! Get electronegativity-weighted coordination number
   call new_erf_en_ncoord(ncoord_en, mol, &
   & rcov=ceh_cov_radii(mol%num), en=pauling_en_ceh(mol%num))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call get_coordination_number(ncoord_en, mol, lattr, cutoff, &
      & cn_en, dcn_endr, dcn_endL)

   ! Get effective charges and their gradients
   call get_effective_qat(mol, cn_en, wfn%qat, &
      & dcn_endr, dcn_endL, wfn%dqatdr, wfn%dqatdL)

   ! If shell-resolved charges are requested, partition atomic charges on shells
   if(present(nshell)) then
      ii = 0
      do iat = 1, mol%nat
         do ish = 1, nshell(iat)
            wfn%qsh(ii+ish, :) = wfn%qat(iat, :) / real(nshell(iat), wp)
            wfn%dqshdr(:, :, ii+ish, :) = wfn%dqatdr(:, :, iat, :) / real(nshell(iat), wp)
            wfn%dqshdL(:, :, ii+ish, :) = wfn%dqatdL(:, :, iat, :) / real(nshell(iat), wp)
         end do 
         ii = ii + nshell(iat)
      end do 
   end if

end subroutine get_charges_effceh


subroutine test_generic(error, mol, qat, qsh, make_coulomb, ref, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges for this structure
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges for this structure
   real(wp), intent(in), optional :: qsh(:)

   !> Factory to create new onsite thirdorder objects
   procedure(coulomb_maker) :: make_coulomb

   !> Reference value to check against
   real(wp), intent(in) :: ref

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   class(coulomb_type), allocatable :: coulomb
   type(container_cache) :: cache
   real(wp) :: energy(mol%nat), sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), numgrad(:, :)
   real(wp) :: thr_
   type(wavefunction_type) :: wfn

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   if (present(qsh)) then
      wfn%qsh = reshape(qsh, [size(qsh), 1])
   else
      wfn%qsh = reshape(qat, [size(qat), 1])
   end if
   wfn%qat = reshape(qat, [size(qat), 1])
   call make_coulomb(coulomb, mol, present(qsh))
   call coulomb%update(mol, cache)
   call coulomb%get_energy(mol, cache, wfn, energy)

   call check(error, sum(energy), ref, thr=thr_)
   if (allocated(error)) then
      print*,ref, sum(energy)
   end if

end subroutine test_generic


subroutine test_numgrad(error, mol, qat, qsh, make_coulomb)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges for this structure
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges for this structure
   real(wp), intent(in), optional :: qsh(:)

   !> Factory to create new onsite thirdorder objects
   procedure(coulomb_maker) :: make_coulomb

   integer :: iat, ic
   class(coulomb_type), allocatable :: coulomb
   type(container_cache) :: cache
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), numgrad(:, :)
   real(wp), parameter :: step = 5.0e-5_wp
   type(wavefunction_type) :: wfn

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   if (present(qsh)) then
      wfn%qsh = reshape(qsh, [size(qsh), 1])
   else
      wfn%qsh = reshape(qat, [size(qat), 1])
   end if
   wfn%qat = reshape(qat, [size(qat), 1])
   call make_coulomb(coulomb, mol, present(qsh))

   do iat = 1, mol%nat
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call coulomb%update(mol, cache)
         call coulomb%get_energy(mol, cache, wfn, er)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call coulomb%update(mol, cache)
         call coulomb%get_energy(mol, cache, wfn, el)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   call coulomb%update(mol, cache)
   call coulomb%get_gradient(mol, cache, wfn, gradient, sigma)

   if (any(abs(gradient - numgrad) > thr2)) then
      call test_failed(error, "Gradient of energy does not match")
      print'(3es21.14)', gradient-numgrad
   end if

end subroutine test_numgrad


subroutine test_numpotgrad(error, mol, get_charges, make_coulomb, shell)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Procedure to provide atomic charges (and their gradients)
   procedure(charge_maker) :: get_charges

   !> Factory to create new onsite thirdorder objects
   procedure(coulomb_maker) :: make_coulomb

   !> Test shell-resolved
   logical, intent(in) :: shell

   integer :: iat, ic
   type(basis_type) :: bas
   type(potential_type) :: potl, potr
   class(coulomb_type), allocatable :: coulomb
   type(container_cache) :: cache
   real(wp), allocatable :: numpotgrad(:, :, :)
   integer, allocatable :: nshell(:)
   real(wp), parameter :: step = 5.0e-5_wp
   type(wavefunction_type) :: wfn

   ! Setup potentials and wavefunction with dummy basis set 
   call make_basis(bas, mol, 6)
   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, 1, 0.0_wp, .true.)
   call new_potential(potr, mol, bas, 1, .true.)
   call new_potential(potl, mol, bas, 1, .true.)
   call make_coulomb(coulomb, mol, shell)

   if(shell) then
      allocate(numpotgrad(3, mol%nat, bas%nsh), source=0.0_wp)
      allocate(nshell(mol%nat))
      nshell = bas%nsh_at
   else
      allocate(numpotgrad(3, mol%nat, mol%nat), source=0.0_wp)
   end if
   
   do iat = 1, mol%nat
      do ic = 1, 3
         call potr%reset
         call potl%reset
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_charges(wfn, mol, nshell)
         call coulomb%update(mol, cache)
         call coulomb%get_potential(mol, cache, wfn, potr)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call get_charges(wfn, mol, nshell)
         call coulomb%update(mol, cache)
         call coulomb%get_potential(mol, cache, wfn, potl)
         
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         if(shell) then
            numpotgrad(ic, iat, :) = 0.5_wp*(potr%vsh(:,1) - potl%vsh(:,1))/step
         else
            numpotgrad(ic, iat, :) = 0.5_wp*(potr%vat(:,1) - potl%vat(:,1))/step
         end if
      end do
   end do

   call get_charges(wfn, mol, nshell)
   call coulomb%update(mol, cache)
   call coulomb%get_potential(mol, cache, wfn, potl)
   call coulomb%get_potential_gradient(mol, cache, wfn, potl)

   if(shell) then
      if (any(abs(potl%dvshdr(:,:,:,1) - numpotgrad) > thr2)) then
         call test_failed(error, "Gradient of shell-resolved potential does not match")
         print'(3es21.14)', potl%dvshdr(:,:,:,1) - numpotgrad
      end if
   else
      if (any(abs(potl%dvatdr(:,:,:,1) - numpotgrad) > thr2)) then
         call test_failed(error, "Gradient of atom-resolved potential does not match")
         print'(3es21.14)', potl%dvatdr(:,:,:,1) - numpotgrad
      end if
   end if

end subroutine test_numpotgrad

subroutine test_numsigma(error, mol, qat, qsh, make_coulomb)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges for this structure
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges for this structure
   real(wp), intent(in), optional :: qsh(:)

   !> Factory to create new electrostatic objects
   procedure(coulomb_maker) :: make_coulomb

   integer :: ic, jc
   class(coulomb_type), allocatable :: coulomb
   type(container_cache) :: cache
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3), eps(3, 3), numsigma(3, 3)
   real(wp), allocatable :: gradient(:, :), xyz(:, :), lattice(:, :)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 1.0e-5_wp
   type(wavefunction_type) :: wfn

   allocate(gradient(3, mol%nat), xyz(3, mol%nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   if (present(qsh)) then
      wfn%qsh = reshape(qsh, [size(qsh), 1])
   else
      wfn%qsh = reshape(qat, [size(qat), 1])
   end if
   wfn%qat = reshape(qat, [size(qat), 1])
   call make_coulomb(coulomb, mol, present(qsh))

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   if (any(mol%periodic)) lattice = mol%lattice
   do ic = 1, 3
      do jc = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call coulomb%update(mol, cache)
         call coulomb%get_energy(mol, cache, wfn, er)
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call coulomb%update(mol, cache)
         call coulomb%get_energy(mol, cache, wfn, el)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         if (allocated(lattice)) mol%lattice = lattice
         numsigma(jc, ic) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   call coulomb%update(mol, cache)
   call coulomb%get_gradient(mol, cache, wfn, gradient, sigma)

   if (any(abs(sigma - numsigma) > thr2)) then
      call test_failed(error, "Strain derivatives do not match")
      print'(3es21.14)', sigma-numsigma
   end if

end subroutine test_numsigma


subroutine test_numpotsigma(error, mol, get_charges, make_coulomb, shell)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Procedure to provide atomic charges (and their gradients)
   procedure(charge_maker) :: get_charges

   !> Factory to create new electrostatic objects
   procedure(coulomb_maker) :: make_coulomb

   !> Test shell-resolved
   logical, intent(in) :: shell

   integer :: ic, jc
   type(basis_type) :: bas
   type(potential_type) :: potl, potr
   class(coulomb_type), allocatable :: coulomb
   type(container_cache) :: cache
   real(wp) :: eps(3, 3)
   real(wp), allocatable :: numpotsigma(:, :, :), xyz(:, :), lattice(:, :)
   real(wp), parameter :: unity(3, 3) = reshape(&
   & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 5.0e-5_wp
   type(wavefunction_type) :: wfn
   integer, allocatable :: nshell(:)

   allocate(xyz(3, mol%nat), source=0.0_wp)

   ! Setup potentials and wavefunction with dummy basis set 
   call make_basis(bas, mol, 6)
   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, 1, 0.0_wp, .true.)
   call new_potential(potr, mol, bas, 1, .true.)
   call new_potential(potl, mol, bas, 1, .true.)
   call make_coulomb(coulomb, mol, shell)

   if(shell) then
      allocate(numpotsigma(3, 3, bas%nsh), source=0.0_wp)
      allocate(nshell(mol%nat))
      nshell = bas%nsh_at
   else
      allocate(numpotsigma(3, 3, mol%nat), source=0.0_wp)
   end if

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   if (any(mol%periodic)) lattice = mol%lattice
   do ic = 1, 3
      do jc = 1, 3
         call potr%reset
         call potl%reset
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call get_charges(wfn, mol, nshell)
         call coulomb%update(mol, cache)
         call coulomb%get_potential(mol, cache, wfn, potr)

         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call get_charges(wfn, mol, nshell)
         call coulomb%update(mol, cache)
         call coulomb%get_potential(mol, cache, wfn, potl)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         if (allocated(lattice)) mol%lattice = lattice
         if(shell) then
            numpotsigma(jc, ic, :) = 0.5_wp*(potr%vsh(:,1) - potl%vsh(:,1))/step
         else
            numpotsigma(jc, ic, :) = 0.5_wp*(potr%vat(:,1) - potl%vat(:,1))/step
         end if
      end do
   end do

   call get_charges(wfn, mol, nshell)
   call coulomb%update(mol, cache)
   call coulomb%get_potential(mol, cache, wfn, potl)
   call coulomb%get_potential_gradient(mol, cache, wfn, potl)

   if(shell) then
      if (any(abs(potl%dvshdL(:,:,:,1) - numpotsigma) > thr2)) then
         call test_failed(error, "Sigma of shell-resolved potential does not match")
         print'(3es21.14)', potl%dvshdL(:,:,:,1) - numpotsigma
      end if
   else
      if (any(abs(potl%dvatdL(:,:,:,1) - numpotsigma) > thr2)) then
         call test_failed(error, "Sigma of atom-resolved potential does not match")
         print'(3es21.14)', potl%dvatdL(:,:,:,1) - numpotsigma
      end if
   end if

end subroutine test_numpotsigma


subroutine test_e_gfn1_m01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 0.80928533739041852_wp, -0.0985788702919893_wp, -1.1789498512968015_wp, & 
      &-0.07630466864269804_wp, -0.5449852206782641_wp,  0.3220779574144627_wp, & 
      &-0.02986669797918212_wp, -1.1079458022741191_wp, -0.6616262018005230_wp, & 
      & 0.56887324141215356_wp,  0.3231023443889342_wp,  0.0959376082660442_wp, &
      & 0.27886474952257534_wp,  0.8972931001216102_wp, -0.3079177622577322_wp, & 
      & 0.71074073670510574_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "01")
   call test_generic(error, mol, qat, qsh, make_coulomb_o1, 5.9181986430716059E-002_wp)

end subroutine test_e_gfn1_m01


subroutine test_e_gfn2_m02(error)

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
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "02")
   call test_generic(error, mol, qat, qsh, make_coulomb_o2, 4.9292337881703292E-002_wp)

end subroutine test_e_gfn2_m02


subroutine test_e_gfn2_oxacb(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 0.43142759318744994_wp, 0.43150967285357533_wp, 0.43146165224958366_wp, &
      & 0.43147407117646708_wp, 0.90726023139894296_wp, 0.90728336242768870_wp, &
      & 0.90737069516079072_wp, 0.90712544694131525_wp,-0.77945631396614101_wp, &
      &-0.78097308159910250_wp,-0.78002089184600454_wp,-0.77924355344945684_wp, &
      &-0.55876993965122157_wp,-0.55823164158343652_wp,-0.55927229556741787_wp, &
      &-0.55894500773303069_wp]   
   real(wp), parameter :: qsh(*)=[&
      & 0.43142759318744994_wp, 0.43150967285357533_wp, 0.43146165224958366_wp, &
      & 0.43147407117646708_wp, 0.12631625958529813_wp, 0.78094397181364483_wp, &
      & 0.12627921822540056_wp, 0.78100414420228814_wp, 0.12624063356686954_wp, &
      & 0.78113006159392118_wp, 0.12626556398329924_wp, 0.78085988295801601_wp, &
      & 0.19373236130795646_wp,-0.97318867527409747_wp, 0.19324264343503161_wp, &
      &-0.97421572503413412_wp, 0.19355130642179486_wp,-0.97357219826779939_wp, &
      & 0.19382295575165664_wp,-0.97306650920111348_wp, 0.23748355725433878_wp, &
      &-0.79625349690556035_wp, 0.23769535856264201_wp,-0.79592700014607853_wp, &
      & 0.23731240854425084_wp,-0.79658470411166871_wp, 0.23739675443019159_wp, &
      &-0.79634176216322228_wp]

   call get_structure(mol, "X23", "oxacb")
   call test_generic(error, mol, qat, qsh, make_coulomb_o2, 0.10439618995869278_wp)

end subroutine test_e_gfn2_oxacb


subroutine test_e_gfn2_oxacb_sc(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat1(*) = [&
      & 3.41731844312030E-1_wp, 3.41716020106239E-1_wp, 3.41730526585671E-1_wp,&
      & 3.41714427217954E-1_wp, 3.80996046757999E-1_wp, 3.80989821246195E-1_wp,&
      & 3.81000747720282E-1_wp, 3.80990494183703E-1_wp,-3.70406587264474E-1_wp,&
      &-3.70407565207006E-1_wp,-3.70417590212352E-1_wp,-3.70399716470705E-1_wp,&
      &-3.52322260586075E-1_wp,-3.52304269439196E-1_wp,-3.52313440903261E-1_wp,&
      &-3.52298498047004E-1_wp]
   integer, parameter :: supercell(*) = [2, 2, 2]
   real(wp), parameter :: qat(*) = [spread(qat1, 2, product(supercell))]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "oxacb")
   call make_supercell(mol, supercell)
   call test_generic(error, mol, qat, qsh, make_coulomb_o2, &
      & 0.02183660722324880_wp*product(supercell), 1.0e-7_wp)

end subroutine test_e_gfn2_oxacb_sc


subroutine make_supercell(mol, rep)
   type(structure_type), intent(inout) :: mol
   integer, intent(in) :: rep(3)

   real(wp), allocatable :: xyz(:, :), lattice(:, :)
   integer, allocatable :: num(:)
   integer :: i, j, k, c

   num = reshape(spread([mol%num(mol%id)], 2, product(rep)), [product(rep)*mol%nat])
   lattice = reshape(&
      [rep(1)*mol%lattice(:, 1), rep(2)*mol%lattice(:, 2), rep(3)*mol%lattice(:, 3)], &
      shape(mol%lattice))
   allocate(xyz(3, product(rep)*mol%nat))
   c = 0
   do i = 0, rep(1)-1
      do j = 0, rep(2)-1
         do k = 0, rep(3)-1
            xyz(:, c+1:c+mol%nat) = mol%xyz &
               & + spread(matmul(mol%lattice, [real(wp):: i, j, k]), 2, mol%nat)
            c = c + mol%nat
         end do
      end do
   end do

   call new(mol, num, xyz, lattice=lattice)
end subroutine make_supercell


subroutine test_g_gfn1_m03(error)

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
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "03")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_o1)

end subroutine test_g_gfn1_m03

subroutine test_g_gfn2_m04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 9.33596160193497E-2_wp,-3.41088061922851E-1_wp, 7.32474961830646E-2_wp,&
      &-2.21649975471802E-1_wp, 6.24413528413759E-3_wp, 1.07366683260668E-1_wp,&
      & 1.25982547197317E-1_wp, 9.65935501843890E-2_wp, 1.02704543049803E-1_wp,&
      & 1.45380937882263E-1_wp,-1.55978251071729E-1_wp, 3.42948437914661E-1_wp,&
      & 5.65504846503244E-2_wp,-3.37789986050220E-1_wp, 1.13510089629769E-1_wp,&
      &-2.07382246739143E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "04")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_o2)

end subroutine test_g_gfn2_m04

subroutine test_g_gfn1_urea(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 5.55723890858218E-1_wp, 5.55765354442035E-1_wp, 2.50200231242017E-1_wp,&
      & 2.50282053284422E-1_wp, 2.39786980460652E-1_wp, 2.39895142481200E-1_wp,&
      & 2.50103678240412E-1_wp, 2.50425041601730E-1_wp, 2.39464477136495E-1_wp,&
      & 2.40360053062669E-1_wp,-4.38369096728919E-1_wp,-4.38451412936599E-1_wp,&
      &-4.38310020776279E-1_wp,-4.38617373848238E-1_wp,-6.59141030224988E-1_wp,&
      &-6.59117968294813E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "urea")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_o1)

end subroutine test_g_gfn1_urea

subroutine test_s_gfn2_pyrazine(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 1.23705026512437E-1_wp, 1.22537765989959E-1_wp, 1.23932626831231E-1_wp,&
      & 1.22418959326822E-1_wp, 1.23788033684569E-1_wp, 1.23643389058068E-1_wp,&
      & 1.22389811551880E-1_wp, 1.22402837718155E-1_wp, 3.17174594622166E-2_wp,&
      & 3.15585817789390E-2_wp, 3.14290005061543E-2_wp, 3.10358506314526E-2_wp,&
      & 3.11084225356749E-2_wp, 3.11187325528499E-2_wp, 3.10707965296438E-2_wp,&
      & 3.13508824430484E-2_wp,-3.09107070033859E-1_wp,-3.09144123845085E-1_wp,&
      &-3.08473365847409E-1_wp,-3.08483617386745E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "pyrazine")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_o2)

end subroutine test_s_gfn2_pyrazine

subroutine test_e_gfn2_m07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-2.986820554251878_wp,  0.2133150304286247_wp,  1.2219537941124052_wp, &
      & 0.637526091539623_wp,  0.6288673640786957_wp, -0.4793698954628428_wp, &
      &-0.381215207423406_wp, -1.0863256053671664_wp,  0.8906955087366708_wp, &
      & 0.159224917897852_wp,  0.9937121523714770_wp,  0.4332827430406178_wp, &
      & 0.459237705248708_wp, -0.6369517649821776_wp, -0.3853972301505781_wp, &
      & 0.318264950183372_wp]
   real(wp), parameter :: qsh(*) = [&
      &-0.2372082175033831_wp, -2.7496123367484957_wp,  0.2133150304286247_wp, & 
      & 1.2459662589528724_wp, -0.0240124648404671_wp,  0.6375260915396232_wp, &
      & 0.6288673640786957_wp,  0.0427487776207620_wp, -0.4119239254478257_wp, &
      &-0.1101947476357790_wp,  0.1962526334754219_wp, -0.5774678408988283_wp, &
      &-0.0954706513992993_wp, -0.9908549539678670_wp,  0.1217211077680937_wp, &
      & 0.7689744009685770_wp,  0.1592249178978529_wp,  0.2104607013244239_wp, &
      & 0.7889318241389675_wp, -0.0056803730919144_wp,  0.4332827430406178_wp, &
      & 0.4592377052487081_wp,  0.2833371093194456_wp, -0.9202888743016233_wp, &
      & 0.1737487393096157_wp, -0.5591459694601939_wp,  1.0261224370878024_wp, &
      &-0.3777382313598706_wp, -0.3301192555445596_wp]
   real(wp), allocatable :: qsh0(:)

   call get_structure(mol, "MB16-43", "07")
   call test_generic(error, mol, qat, qsh0, make_coulomb_o2, -1.2136915157280399_wp)
   call test_generic(error, mol, qat, qsh, make_coulomb_o2, -0.34488498240315002_wp)

end subroutine test_e_gfn2_m07

subroutine test_g_gfn2_m08(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-2.11048312695985E-1_wp,-5.02011645803230E-1_wp, 4.15238062649689E-1_wp, &
      &-3.25959600753673E-1_wp, 2.51473641195433E-2_wp, 2.93748490123740E-1_wp, &
      & 2.56736194030896E-2_wp, 2.38762690307426E-2_wp,-6.03118603733083E-1_wp, &
      & 3.91990240426822E-1_wp, 8.97114734113785E-1_wp, 1.93532936362436E-1_wp, &
      &-1.03136268223866E-1_wp,-1.04447608767710E-1_wp,-2.64818891498402E-2_wp, &
      &-3.90117787102468E-1_wp]
   real(wp), parameter :: qsh(*) = [&
      & 9.06259904944829E-1_wp,-1.11730821567902E+0_wp, 2.78017329305492E-1_wp, &
      &-7.80028989546297E-1_wp, 1.11352815063389E+0_wp,-6.98290073981154E-1_wp, &
      & 2.03943236255318E-1_wp,-5.29902840441233E-1_wp, 4.38219939650397E-2_wp, &
      &-1.86746328945826E-2_wp, 4.65996457236599E-2_wp, 4.97590807484258E-1_wp, &
      &-2.50441962186972E-1_wp, 4.83295451755440E-2_wp,-2.26559244782012E-2_wp, &
      & 4.50331992744248E-2_wp,-2.11569328297532E-2_wp, 3.12470620007346E-1_wp, &
      &-9.15589243372491E-1_wp, 1.06394261835743E+0_wp,-6.71952361588756E-1_wp, &
      & 1.82322476598938E+0_wp,-9.26110009158329E-1_wp, 9.78357111140355E-1_wp, &
      &-7.84824170464332E-1_wp,-9.43549308434806E-2_wp,-8.78133979988158E-3_wp, &
      &-7.07783143624696E-2_wp,-3.36692984665466E-2_wp, 6.75375129657761E-1_wp, &
      &-7.01857024422455E-1_wp, 2.11598132242645E-1_wp,-6.01715925641418E-1_wp]

   call get_structure(mol, "MB16-43", "08")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_o2)

end subroutine test_g_gfn2_m08

subroutine test_s_gfn2_m09(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-1.13260038539900E-2_wp, 1.10070523471231E-2_wp,-9.97165474630829E-2_wp, &
      &-8.78527301724521E-2_wp, 2.89049242695863E-1_wp, 3.57284006856323E-2_wp, &
      &-1.73226219187217E-1_wp, 1.61174372420268E-1_wp,-8.89089419183055E-2_wp, &
      & 3.23950178196666E-2_wp, 1.88420688366637E-1_wp, 4.14882523279327E-2_wp, &
      &-2.23498403532295E-1_wp,-3.55334728213004E-1_wp,-7.15753987897201E-2_wp, &
      & 3.52175946466941E-1_wp]
   real(wp), parameter :: qsh(*) = [&
      & 1.40887956776581E-3_wp,-1.27348820058716E-2_wp, 2.63961183739554E-2_wp, &
      &-1.53890176131402E-2_wp,-8.73648546608390E-2_wp,-1.23517435478204E-2_wp, &
      &-8.71021014527735E-2_wp,-7.50559382736492E-4_wp, 7.82044211296174E-1_wp, &
      &-4.92995083533018E-1_wp, 4.84143136555792E-2_wp,-1.26858387490357E-2_wp, &
      & 9.72488073646510E-1_wp,-1.14571448042039E+0_wp, 1.07574874045191E+0_wp, &
      &-9.14574293473561E-1_wp,-7.63358458276189E-2_wp,-1.25730981035572E-2_wp, &
      & 4.44349073468088E-2_wp,-1.20397879426510E-2_wp, 5.20245311277456E-1_wp, &
      & 1.92282483805197E-1_wp,-5.24107355799204E-1_wp, 5.39382871928999E-2_wp, &
      &-1.24499232808976E-2_wp, 7.97368410133983E-2_wp,-3.13209082796440E-1_wp, &
      & 9.97387287057362E-3_wp, 1.94446888375020E-1_wp,-5.49781435696375E-1_wp, &
      &-6.89789344411558E-2_wp,-2.59643153694089E-3_wp, 1.09519797601190E+0_wp, &
      &-7.43022154621128E-1_wp]

   call get_structure(mol, "MB16-43", "09")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_o2)

end subroutine test_s_gfn2_m09

subroutine test_ceh_potgrad_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_oceh, .false.)

end subroutine test_ceh_potgrad_lih

subroutine test_ceh_potgrad_m15(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "15")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_oceh, .false.)

end subroutine test_ceh_potgrad_m15

subroutine test_ceh_potgrad_m16(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "16")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_oceh, .true.)

end subroutine test_ceh_potgrad_m16

subroutine test_ceh_potsigma_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_oceh, .false.)

end subroutine test_ceh_potsigma_lih

subroutine test_ceh_potsigma_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "CO2")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_oceh, .false.)

end subroutine test_ceh_potsigma_co2

subroutine test_ceh_potsigma_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_oceh, .false.)

end subroutine test_ceh_potsigma_mb05

subroutine test_ceh_potsigma_mb17(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "17")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_oceh, .true.)

end subroutine test_ceh_potsigma_mb17

end module test_coulomb_thirdorder
