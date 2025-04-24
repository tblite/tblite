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

!> @file tblite/xtb/singlepoint.f90
!> Provides main entry point for performing single point calculations with
!> a #xtb_calculator instance.

!> Implementation of the single point calculation for a self-consistent
!> extended tight-binding Hamiltonian.
module tblite_xtb_singlepoint
   use mctc_env, only : wp, error_type, fatal_error, get_variable
   use mctc_io, only : structure_type
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_basis_type, only : get_cutoff, basis_type
   use tblite_blas, only : gemv
   use tblite_container, only : container_cache
   use tblite_context, only : context_type, escape
   use tblite_cutoff, only : get_lattice_points
   use tblite_integral_type, only : integral_type, new_integral
   use tblite_lapack_sygvr, only : sygvr_solver
   use tblite_output_format, only : format_string
   use tblite_results, only : results_type
   use tblite_scf, only : mixer_type, new_mixer, scf_info, next_scf, &
      & get_mixer_dimension, potential_type, new_potential
   use tblite_scf_solver, only : solver_type
   use tblite_timer, only : timer_type, format_time
   use tblite_wavefunction, only : wavefunction_type, get_density_matrix, &
      & get_alpha_beta_occupation, &
      & magnet_to_updown, updown_to_magnet
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_h0, only : get_selfenergy, get_hamiltonian, get_occupation, &
      & get_hamiltonian_gradient
   use tblite_post_processing_type, only : collect_containers_caches
   use tblite_post_processing_list, only : post_processing_list
   implicit none
   private

   public :: xtb_singlepoint

   real(wp), parameter :: cn_cutoff = 25.0_wp


   character(len=*), parameter :: real_format = "(es20.13)"
   character(len=25), parameter :: &
      label_halogen = "halogen bonding energy", &
      label_repulsion = "repulsion energy", &
      label_dispersion = "dispersion energy", &
      label_other = "interaction energy", &
      label_electrons = "number of electrons", &
      label_cutoff = "integral cutoff", &
      label_electronic = "electronic energy", &
      label_total = "total energy"

contains


!> Entry point for performing single point calculation using the xTB calculator
subroutine xtb_singlepoint(ctx, mol, calc, wfn, accuracy, energy, gradient, sigma, &
      & verbosity, results, post_process)
   !> Calculation context
   type(context_type), intent(inout) :: ctx
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Accuracy for computation
   real(wp), intent(in) :: accuracy
   !> Total energy
   real(wp), intent(out) :: energy
   !> Gradient with respect to cartesian coordinates
   real(wp), contiguous, intent(out), optional :: gradient(:, :)
   !> Strain derivatives with respect to strain deformations
   real(wp), contiguous, intent(out), optional :: sigma(:, :)
   !> Verbosity level of output
   integer, intent(in), optional :: verbosity
   !> Container for storing additional results
   type(results_type), intent(out), optional :: results
   type(post_processing_list), intent(inout), optional :: post_process
   
   logical :: grad, converged, econverged, pconverged
   integer :: prlevel
   real(wp) :: econv, pconv, cutoff, elast, nel
   real(wp), allocatable :: energies(:), edisp(:), erep(:), exbond(:), eint(:), eelec(:)
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), dEdcn(:)
   real(wp), allocatable :: selfenergy(:), dsedcn(:), lattr(:, :), wdensity(:, :, :)
   type(integral_type) :: ints
   type(potential_type) :: pot
   type(container_cache), allocatable :: ccache, dcache, icache, hcache, rcache
   class(mixer_type), allocatable :: mixer
   type(timer_type) :: timer
   type(error_type), allocatable :: error
   type(scf_info) :: info
   class(solver_type), allocatable :: solver
   type(adjacency_list) :: list
   type(container_cache), allocatable :: cache_list(:)
   
   integer :: iscf, spin

   call timer%push("total")

   if (present(verbosity)) then
      prlevel = verbosity
   else
      prlevel = ctx%verbosity
   end if

   econv = 1.e-6_wp*accuracy
   pconv = 2.e-5_wp*accuracy

   grad = present(gradient) .and. present(sigma)

   allocate(energies(mol%nat), source=0.0_wp)
   allocate(erep(mol%nat), source=0.0_wp)
   allocate(edisp(mol%nat), source=0.0_wp)
   allocate(eint(mol%nat), source=0.0_wp)
   allocate(exbond(mol%nat), source=0.0_wp)
   if (grad) then
      gradient(:, :) = 0.0_wp
      sigma(:, :) = 0.0_wp
   end if

   call get_occupation(mol, calc%bas, calc%h0, wfn%nocc, wfn%n0at, wfn%n0sh)
   nel = sum(wfn%n0at) - mol%charge
   if (mod(mol%uhf, 2) == mod(nint(nel), 2)) then
      wfn%nuhf = mol%uhf
   else
      if (mol%uhf /= 0) then
         call fatal_error(error, "Total number of electrons ("//format_string(nint(nel), "(i0)")//") and "//&
            & "number unpaired electrons ("//format_string(mol%uhf, "(i0)")//") is not compatible")
         call ctx%set_error(error)
         return
      end if
      wfn%nuhf = mod(nint(nel), 2)
   end if
   call get_alpha_beta_occupation(wfn%nocc, wfn%nuhf, wfn%nel(1), wfn%nel(2))

   if (allocated(calc%halogen)) then
      call timer%push("halogen")
      allocate(hcache)
      call calc%halogen%update(mol, hcache)
      call calc%halogen%get_engrad(mol, hcache, exbond, gradient, sigma)
      if (prlevel > 1) &
         call ctx%message(label_halogen // format_string(sum(exbond), real_format) // " Eh")
      energies(:) = energies + exbond
      call timer%pop
   end if

   if (allocated(calc%repulsion)) then
      call timer%push("repulsion")
      allocate(rcache)
      call calc%repulsion%update(mol, rcache)
      call calc%repulsion%get_engrad(mol, rcache, erep, gradient, sigma)
      if (prlevel > 1) &
         call ctx%message(label_repulsion // format_string(sum(erep), real_format) // " Eh")
      energies(:) = energies + erep
      call timer%pop
   end if

   if (allocated(calc%dispersion)) then
      call timer%push("dispersion")
      allocate(dcache)
      call calc%dispersion%update(mol, dcache)
      call calc%dispersion%get_engrad(mol, dcache, edisp, gradient, sigma)
      if (prlevel > 1) &
         call ctx%message(label_dispersion // format_string(sum(edisp), real_format) // " Eh")
      energies(:) = energies + edisp
      call timer%pop
   end if

   if (allocated(calc%interactions)) then
      call timer%push("interactions")
      allocate(icache)
      call calc%interactions%update(mol, icache)
      call calc%interactions%get_engrad(mol, icache, eint, gradient, sigma)
      if (prlevel > 1) &
         call ctx%message(label_other // format_string(sum(eint), real_format) // " Eh")
      energies(:) = energies + eint
      call timer%pop
   end if

   call new_potential(pot, mol, calc%bas, wfn%nspin)
   if (allocated(calc%coulomb)) then
      allocate(ccache)
      call timer%push("coulomb")
      call calc%coulomb%update(mol, ccache)
      call timer%pop
   end if

   if (prlevel > 1) &
      call ctx%message(label_electrons // format_string(wfn%nocc, real_format) // " e")

   call timer%push("hamiltonian")
   if (allocated(calc%ncoord)) then
      allocate(cn(mol%nat))
      if (grad) then
         allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
      end if
      call calc%ncoord%get_cn(mol, cn, dcndr, dcndL)
   end if

   allocate(selfenergy(calc%bas%nsh), dsedcn(calc%bas%nsh))
   call get_selfenergy(calc%h0, mol%id, calc%bas%ish_at, calc%bas%nsh_id, cn=cn, &
      & selfenergy=selfenergy, dsedcn=dsedcn)

   cutoff = get_cutoff(calc%bas, accuracy)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call new_adjacency_list(list, mol, lattr, cutoff)

   if (prlevel > 1) then
      call ctx%message(label_cutoff // format_string(cutoff, real_format) // " bohr")
      call ctx%message("")
   end if

   call new_integral(ints, calc%bas%nao)
   call get_hamiltonian(mol, lattr, list, calc%bas, calc%h0, selfenergy, &
      & ints%overlap, ints%dipole, ints%quadrupole, ints%hamiltonian)

   call ctx%new_solver(solver, ints%overlap, wfn%nel, wfn%kt)
   call timer%pop

   call timer%push("scc")
   allocate(eelec(mol%nat), source=0.0_wp)
   iscf = 0
   converged = .false.
   info = calc%variable_info()
   call new_mixer(mixer, calc%max_iter, wfn%nspin*get_mixer_dimension(mol, calc%bas, info), &
      & calc%mixer_damping)
   if (prlevel > 0) then
      call ctx%message(repeat("-", 60))
      call ctx%message("  cycle        total energy    energy error   density error")
      call ctx%message(repeat("-", 60))
   end if
   do while(.not.converged .and. iscf < calc%max_iter)
      elast = sum(eelec)
      call next_scf(iscf, mol, calc%bas, wfn, solver, mixer, info, &
         & calc%coulomb, calc%dispersion, calc%interactions, ints, pot, &
         & ccache, dcache, icache, eelec, error)
      econverged = abs(sum(eelec) - elast) < econv
      pconverged = mixer%get_error() < pconv
      converged = econverged .and. pconverged
      if (prlevel > 0) then
         call ctx%message(format_string(iscf, "(i7)") // &
            & format_string(sum(eelec + energies), "(g24.13)") // &
            & escape(merge(ctx%terminal%green, ctx%terminal%red, econverged)) // &
            & format_string(sum(eelec) - elast, "(es16.7)") // &
            & escape(merge(ctx%terminal%green, ctx%terminal%red, pconverged)) // &
            & format_string(mixer%get_error(), "(es16.7)") // &
            & escape(ctx%terminal%reset))
      end if
      if (allocated(error)) then
         call ctx%set_error(error)
         exit
      end if
   end do
   if (prlevel > 0) then
      call ctx%message(repeat("-", 60))
      call ctx%message("")
   end if
   energies(:) = energies + eelec
   energy = sum(energies)
   if (present(results)) then
      results%energies = energies
   end if
   call timer%pop

   if (prlevel > 1) then
      call ctx%message(label_electronic // format_string(sum(eelec), real_format) // " Eh")
      call ctx%message(label_total // format_string(sum(energies), real_format) // " Eh")
      call ctx%message("")
   end if

   if (ctx%failed()) then
      call ctx%delete_solver(solver)
      return
   end if

   if (grad) then
      if (allocated(calc%coulomb)) then
         call timer%push("coulomb")
         call calc%coulomb%get_gradient(mol, ccache, wfn, gradient, sigma)
         call timer%pop
      end if

      if (allocated(calc%dispersion)) then
         call timer%push("dispersion")
         call calc%dispersion%get_gradient(mol, dcache, wfn, gradient, sigma)
         call timer%pop
      end if

      if (allocated(calc%interactions)) then
         call timer%push("interactions")
         call calc%interactions%get_gradient(mol, icache, wfn, gradient, sigma)
         call timer%pop
      end if

      call timer%push("hamiltonian")
      allocate(dEdcn(mol%nat))
      dEdcn(:) = 0.0_wp

      allocate(wdensity(calc%bas%nao, calc%bas%nao, wfn%nspin))
      call solver%get_wdensity(wfn%coeff, ints%overlap, wfn%emo, wfn%focc, wdensity, error)
      call updown_to_magnet(wfn%density)
      call updown_to_magnet(wdensity)
      call get_hamiltonian_gradient(mol, lattr, list, calc%bas, calc%h0, selfenergy, &
         & dsedcn, pot, wfn%density, wdensity, dEdcn, gradient, sigma)
      call magnet_to_updown(wfn%density)

      if (allocated(dcndr)) then
         call gemv(dcndr, dEdcn, gradient, beta=1.0_wp)
      end if
      if (allocated(dcndL)) then
         call gemv(dcndL, dEdcn, sigma, beta=1.0_wp)
      end if
      call timer%pop
   end if

   call ctx%delete_solver(solver)
   if (ctx%failed()) return

   if (present(post_process)) then
      call timer%push("post processing")
      call collect_containers_caches(rcache, ccache, hcache, dcache, icache, calc, cache_list)
      call post_process%compute(mol, wfn, ints, calc, cache_list, ctx, prlevel)
      call post_process%print_csv(mol)
      if (prlevel > 1) call ctx%message(post_process%info(prlevel, " | "))
      call post_process%print_timer(prlevel, ctx)
      deallocate(cache_list)
      call timer%pop()
   end if
   if (present(results)) then
      if (allocated(results%dict)) deallocate(results%dict)
      allocate(results%dict)
      if (present(post_process)) then 
         call post_process%pack_res(mol, results)
      end if
   end if

   if (calc%save_integrals .and. present(results)) then
      call move_alloc(ints%overlap, results%overlap)
      call move_alloc(ints%hamiltonian, results%hamiltonian)
   end if

   block
      integer :: it
      real(wp) :: ttime, stime
      character(len=*), parameter :: label(*) = [character(len=20):: &
         & "repulsion", "halogen", "dispersion", "interactions", &
         & "coulomb", "hamiltonian", "post processing", "scc"]
      if (prlevel > 0) then
         ttime = timer%get("total")
         call ctx%message(" total:"//repeat(" ", 16)//format_time(ttime))
      end if
      if (prlevel > 1) then
         do it = 1, size(label)
            stime = timer%get(label(it))
            if (stime <= epsilon(0.0_wp)) cycle
            call ctx%message(" - "//label(it)//format_time(stime) &
               & //" ("//format_string(int(stime/ttime*100), '(i3)')//"%)")
         end do
         call ctx%message("")
      end if
   end block

   if (.not.converged) then
      call fatal_error(error, "SCF not converged in "//format_string(iscf, '(i0)')//" cycles")
      call ctx%set_error(error)
   end if


end subroutine xtb_singlepoint

end module tblite_xtb_singlepoint
