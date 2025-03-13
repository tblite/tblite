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
   use tblite_solvation_ddx, only : ddx_solvation, new_ddx, ddx_input
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction, eeq_guess
   use tblite_basis_type, only : basis_type
   use tblite_integral_type, only : integral_type
   use tblite_context_type, only : context_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_wavefunction_mulliken, only : get_mulliken_shell_charges
   use tblite_scf_iterator, only : get_qat_from_qsh
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   use tblite_integral_type, only : new_integral
   use tblite_scf_potential, only : add_pot_to_h1, new_potential
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_xtb_h0, only : get_selfenergy, get_hamiltonian
   use tblite_cutoff, only : get_lattice_points
   use tblite_basis_type, only : get_cutoff

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
      new_unittest("energy-mol", test_e_m01), &
      new_unittest("gradient-mol", test_g_m02), &
      new_unittest("potential-mol", test_p_m03) &
      ! new_unittest("fock", test_fock_m01) &
      ]

end subroutine collect_solvation_ddx


subroutine test_e(error, mol, qat, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Reference energy
   real(wp), intent(in) :: ref

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

   solv = ddx_solvation(mol, ddx_input(feps, 1, nang=nang, rscale=rscale))

   call solv%update(mol, cache)

   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)

   if (abs(sum(energy) - ref) < thr) then
      call test_failed(error, "Energy does not match reference")
      print *, sum(energy)
   end if
end subroutine test_e


subroutine test_g(error, mol, qat)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

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

   solv = ddx_solvation(mol, ddx_input(feps, 1, nang=nang, rscale=rscale))

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
end subroutine test_g


subroutine test_p(error, mol, qat)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   type(ddx_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache), pointer :: cache
   real(wp), parameter :: feps = 80.0_wp, rscale = 1.0_wp
   integer, parameter :: nang = 26
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), parameter :: thr = 1e+3_wp*sqrt(epsilon(1.0_wp))
   real(wp), allocatable :: vat(:)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat)
   integer :: ii

   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1))

   solv = ddx_solvation(mol, ddx_input(feps, 1, nang=nang, rscale=rscale))

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


subroutine test_fock(error, mol, wfn, calc)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Basis set information
   type(xtb_calculator), intent(in) :: calc
   !> Integral container
   type(integral_type):: ints
   
   type(ddx_solvation) :: solv
   type(potential_type) :: pot
   type(container_cache), pointer :: cache
   type(adjacency_list) :: list
   
   real(wp), parameter :: feps = 80.0_wp, rscale = 1.0_wp
   integer, parameter :: nang = 26
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), parameter :: thr = 1e+3_wp*sqrt(epsilon(1.0_wp))
   real(wp) :: er, el
   real(wp), allocatable :: energies(:)
   real(wp), allocatable :: fock(:, :, :), temp(:,:,:), temp_P(:,:,:)
   real(wp), allocatable :: numfock(:,:)
   real(wp), allocatable :: lattr(:, :), selfenergy(:)
   integer :: i, ci, cj
   real(wp) :: cutoff

   solv = ddx_solvation(mol, ddx_input(feps, 1, nang=nang, rscale=rscale))

   allocate(cache)
   call solv%update(mol, cache)
   
   call new_potential(pot, mol, calc%bas, 1)
   call pot%reset
   call solv%get_potential(mol, cache, wfn, pot)

   allocate(selfenergy(calc%bas%nsh))
   call get_selfenergy(calc%h0, mol%id, calc%bas%ish_at, calc%bas%nsh_id, &
      & selfenergy=selfenergy)

   cutoff = get_cutoff(calc%bas, acc)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call new_adjacency_list(list, mol, lattr, cutoff)
   call new_integral(ints, calc%bas%nao)
   call get_hamiltonian(mol, lattr, list, calc%bas, calc%h0, selfenergy, &
      & ints%overlap, ints%dipole, ints%quadrupole, ints%hamiltonian)

   allocate(fock(size(wfn%coeff, 1), size(wfn%coeff, 2), size(wfn%coeff, 3)))
   fock = 0.0_wp
   ints%hamiltonian = 0.0_wp
   call add_pot_to_h1(calc%bas, ints, pot, fock)

   !> Get numerical Fock matrix
   allocate(energies(mol%nat))
   allocate(numfock(calc%bas%nao, calc%bas%nao), source=0.0_wp)
   allocate(temp(size(wfn%coeff, 1), size(wfn%coeff, 2), size(wfn%coeff, 3)))
   allocate(temp_P(size(wfn%density, 1), size(wfn%density, 2), size(wfn%density, 3)))
   do ci = 1, calc%bas%nao
      do cj = 1, calc%bas%nao

         temp = 0.0_wp
         temp_P = 0.0_wp

         wfn%coeff(ci, cj, 1) = wfn%coeff(ci, cj, 1) + step

         ! P = C nocc C^T
         ! Store C nocc in temp
         do i = 1, calc%bas%nao
            temp(i, i, 1) = wfn%focc(i, 1) 
         end do
         temp(:,:,1) = matmul(wfn%coeff(:, :, 1), temp(:,:,1))
         temp_P(:,:,1) = matmul(temp(:,:,1), transpose(wfn%coeff(:, :, 1)))
         wfn%density(:, :, 1) = temp_P(:, :, 1)

         wfn%qsh = 0.0_wp
         call get_mulliken_shell_charges(calc%bas, ints%overlap, wfn%density, wfn%n0sh, &
            & wfn%qsh)
         wfn%qat = 0.0_wp
         call get_qat_from_qsh(calc%bas, wfn%qsh, wfn%qat)

         energies = 0.00_wp
         call pot%reset
         ! call solv%get_potential(mol, cache, wfn, pot)
         call solv%get_energy(mol, cache, wfn, energies)
         er = 0.0_wp
         er = sum(energies)

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         temp = 0.0_wp
         temp_P = 0.0_wp

         wfn%coeff(ci, cj, 1) = wfn%coeff(ci, cj, 1) - 2*step

         ! P = CnoccC^T
         ! Store Cnocc in temp
         do i = 1, calc%bas%nao
            temp(i, i, 1) = wfn%focc(i, 1) 
         end do
         temp(:,:,1) = matmul(wfn%coeff(:, :, 1), temp(:,:,1))
         temp_P(:,:,1) = matmul(temp(:,:,1), transpose(wfn%coeff(:, :, 1)))
         wfn%density(:, :, 1) = temp_P(:, :, 1)

         wfn%qsh = 0.0_wp
         call get_mulliken_shell_charges(calc%bas, ints%overlap, wfn%density, wfn%n0sh, &
            & wfn%qsh)
         wfn%qat = 0.0_wp
         call get_qat_from_qsh(calc%bas, wfn%qsh, wfn%qat)

         energies = 0.00_wp
         call pot%reset
         call solv%get_potential(mol, cache, wfn, pot)
         call solv%get_energy(mol, cache, wfn, energies)

         el = 0.0_wp
         el = sum(energies)

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         numfock(ci, cj) = (er -  el) / (2.0_wp * step)

         wfn%coeff(ci, cj, 1) = wfn%coeff(ci, cj, 1) + step

         end do
   end do

   !> Get analytical Fock matrix
   fock(:,:,1) = matmul(fock(:,:,1), wfn%coeff(:,:,1))
   temp = 0.0_wp
   do i = 1, calc%bas%nao
      temp(i, i, 1) = wfn%focc(i, 1) 
   end do
   fock(:,:,1) = matmul(fock(:, :, 1), temp(:,:,1))
   fock(:,:,1) = fock(:,:,1) * 2.0_wp

   if (any(abs(fock(:,:,1) - numfock(:,:)) > thr)) then
      call test_failed(error, "Potential does not match")
      print '(a)', 'analytical'
      print '(3es20.13)', fock
      print '(a)', "---"
      print '(a)', 'numerical'
      print '(3es20.13)', numfock
      print '(a)', "---"
      print '(a)', 'diff'
      print '(3es20.13)', fock(:,:,1) - numfock(:,:)
   end if

end subroutine test_fock


subroutine test_e_m01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 7.73347900345264E-1_wp, 1.07626888948184E-1_wp,-3.66999593831010E-1_wp,&
      & 4.92833325937897E-2_wp,-1.83332156197733E-1_wp, 2.33302086605469E-1_wp,&
      & 6.61837152062315E-2_wp,-5.43944165050002E-1_wp,-2.70264356583716E-1_wp,&
      & 2.66618968841682E1_wp, 2.62725033202480E-1_wp,-7.15315510172571E-2_wp,&
      &-3.73300777019193E-1_wp, 3.84585237785621E-2_wp,-5.05851088366940E-1_wp,&
      & 5.17677238544189E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_e(error, mol, qat, -2.1546508620217987E-3_wp)

end subroutine test_e_m01


subroutine test_g_m02(error)

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
   call test_g(error, mol, qat)

end subroutine test_g_m02


subroutine test_p_m03(error)

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
   call test_p(error, mol, qat)

end subroutine test_p_m03


subroutine test_fock_m01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(integral_type) :: ints  
   ! type(potential_type) :: pot
   type(wavefunction_type) :: wfn
   real(wp) :: energy

   call get_structure(mol, "MB16-43", "01")

   energy = 0.0_wp

   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   call eeq_guess(mol, calc, wfn)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)

   call test_fock(error, mol, wfn, calc)

end subroutine test_fock_m01



end module test_solvation_ddx
