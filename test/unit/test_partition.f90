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

module test_partition
   use dftd3_partition, only : owns_d3_pair => owns_pair
   use dftd4_partition, only : owns_d4_pair => owns_pair
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_basis_type, only : get_cutoff
   use tblite_container, only : container_cache, container_type
   use tblite_container_list, only : cache_list
   use tblite_context, only : context_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_external_field, only : electric_field, new_electric_field
   use tblite_features, only : get_tblite_feature, tblite_has_mpi
   use tblite_integral_type, only : integral_type, new_integral
   use tblite_mpi_utils, only : get_mpi_comm_world, mpi_allreduce_sum, mpi_sync_error
   use tblite_partition, only : work_partition, new_work_partition, &
      & owns_index, owns_pair, serial_work_partition
   use tblite_post_processing_list, only : post_processing_list, add_post_processing
   use tblite_results, only : results_type
   use tblite_scf_potential, only : potential_type, new_potential
   use tblite_solvation, only : solvation_input, solvation_type, alpb_input, cds_input, &
      & new_solvation, new_solvation_cds
   use tblite_timer, only : timer_type
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_h0, only : get_selfenergy, get_hamiltonian
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: collect_partition

   integer, parameter :: nparts = 3
   real(wp), parameter :: thr = 1.0e-9_wp

contains


!> Collect all exported unit tests
subroutine collect_partition(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("invalid", test_invalid), &
      new_unittest("context", test_context), &
      new_unittest("ownership", test_ownership), &
      new_unittest("absent", test_absent), &
      new_unittest("feature", test_feature), &
      new_unittest("mpi-unavailable", test_mpi_unavailable), &
      new_unittest("scf-partition", test_scf_partition), &
      new_unittest("post-processing", test_post_processing), &
      new_unittest("gfn1-mol", test_gfn1_mol), &
      new_unittest("gfn2-mol", test_gfn2_mol), &
      new_unittest("gfn1-pbc", test_gfn1_pbc), &
      new_unittest("gfn2-pbc", test_gfn2_pbc), &
      new_unittest("hamiltonian", test_hamiltonian), &
      new_unittest("solvation", test_solvation), &
      new_unittest("field", test_field) &
      ]

end subroutine collect_partition


!> Out of range parts have to be rejected
subroutine test_invalid(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(work_partition) :: partition
   type(error_type), allocatable :: partition_error
   integer :: icase
   integer, parameter :: invalid(2, 4) = reshape(&
      & [-1, 3, 3, 3, 4, 3, 0, 0], [2, 4])

   call new_work_partition(error, partition, 1, nparts)
   if (allocated(error)) return
   do icase = 1, size(invalid, 2)
      call new_work_partition(partition_error, partition, &
         & invalid(1, icase), invalid(2, icase))
      call check(error, allocated(partition_error), "Invalid partition was accepted")
      if (allocated(error)) return
      call check(error, partition%get_part(), 1)
      if (allocated(error)) return
      call check(error, partition%get_nparts(), nparts)
      if (allocated(error)) return
   end do

end subroutine test_invalid


!> The context carries the partition for a black box calculation
subroutine test_context(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(error_type), allocatable :: partition_error

   if (.not.owns_pair(ctx%partition, 3, 2)) then
      call test_failed(error, "Default context does not own the complete work")
      return
   end if

   call ctx%set_partition(1, nparts, partition_error)
   if (allocated(partition_error)) then
      call test_failed(error, "Valid work partition was rejected by the context")
      return
   end if

   if (owns_pair(ctx%partition, 1, 1)) then
      call test_failed(error, "Context partition was not applied")
      return
   end if

   call ctx%set_partition(nparts, nparts, partition_error)
   if (.not.allocated(partition_error)) then
      call test_failed(error, "Invalid work partition was accepted by the context")
      return
   end if

   call check(error, ctx%partition%get_part(), 1)
   if (allocated(error)) return
   call check(error, ctx%partition%get_nparts(), nparts)
   if (allocated(error)) return

   ! Setting an external partition must discard any former MPI reduction mode.
   ctx%comm = 0
   call ctx%set_partition(0, 1, error)
   if (allocated(error)) return
   call check(error, .not.allocated(ctx%comm))

end subroutine test_context


!> Every index/pair has one owner; the D3/D4 adapters select the same work.
subroutine test_ownership(error)
   type(error_type), allocatable, intent(out) :: error
   type(work_partition) :: partitions(nparts)
   integer :: iat, jat, part

   do part = 1, nparts
      call new_work_partition(error, partitions(part), part-1, nparts)
      if (allocated(error)) return
   end do
   do iat = 1, 12
      call check(error, count(owns_index(partitions, iat)), 1)
      if (allocated(error)) return
      do jat = 1, iat
         call check(error, count(owns_pair(partitions, iat, jat)), 1)
         if (allocated(error)) return
         do part = 1, nparts
            call check(error, owns_pair(partitions(part), iat, jat) .eqv. &
               & owns_d3_pair(partitions(part)%get_d3(), iat, jat))
            if (allocated(error)) return
            call check(error, owns_pair(partitions(part), iat, jat) .eqv. &
               & owns_d4_pair(partitions(part)%get_d4(), iat, jat))
            if (allocated(error)) return
         end do
      end do
   end do
   ! The triangular index must not overflow a default integer.
   call check(error, count(owns_pair(partitions, 100000, 100000)), 1)
end subroutine test_ownership


!> An absent partition selects the complete work
subroutine test_absent(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check(error, owns_index(idx=7) .and. owns_pair(iat=5, jat=3))
   if (allocated(error)) return
   call check(error, owns_index(serial_work_partition, 7) &
      & .and. owns_pair(serial_work_partition, 5, 3))

end subroutine test_absent


!> The MPI feature must be queryable by name and as a compile time constant
subroutine test_feature(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check(error, get_tblite_feature("mpi"), tblite_has_mpi)
   if (allocated(error)) return

   call check(error, .not.get_tblite_feature("not-a-feature"))

end subroutine test_feature


!> Without MPI support the MPI entry points have to report an error rather than
!> silently do nothing, the distributed calculation itself is covered in test/mpi
subroutine test_mpi_unavailable(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(error_type), allocatable :: mpi_error
   real(wp) :: val(1)

   ! Also safe in MPI builds: MPI_Initialized may be called before MPI_Init.
   call ctx%set_mpi(mpi_error)
   call check(error, allocated(mpi_error) .and. .not.allocated(ctx%comm))
   if (allocated(error)) return
   if (tblite_has_mpi) return

   deallocate(mpi_error)
   val = 1.0_wp
   call mpi_allreduce_sum(mpi_error, val, get_mpi_comm_world())
   call check(error, allocated(mpi_error))
   if (allocated(error)) return
   deallocate(mpi_error)
   call mpi_sync_error(mpi_error, get_mpi_comm_world())
   call check(error, allocated(mpi_error))

end subroutine test_mpi_unavailable


!> Reject both mismatched partitions and partitioned SCF without reductions.
subroutine test_scf_partition(error)
   type(error_type), allocatable, intent(out) :: error
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(error_type), allocatable :: failure
   real(wp) :: energy
   integer :: icase
   character(len=16), parameter :: message(*) = [character(len=16) :: &
      & "does not match", "requires an MPI"]

   call get_structure(mol, "MB16-43", "01")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return
   call ctx%set_partition(1, nparts, error)
   if (allocated(error)) return
   ctx%verbosity = 0
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, 300.0_wp)

   do icase = 1, 2
      if (icase == 2) call calc%set_partition(ctx%partition)
      call xtb_singlepoint(ctx, mol, calc, wfn, 1.0_wp, energy)
      call check(error, ctx%failed(), "Invalid SCF partition was accepted")
      if (allocated(error)) return
      call ctx%get_error(failure)
      call check(error, index(failure%message, trim(message(icase))) > 0)
      if (allocated(error)) return
   end do
end subroutine test_scf_partition


!> The xTB-ML features are evaluated from the partitioned caches and cannot be
!> summed afterwards, a distributed calculation has to reject them
subroutine test_post_processing(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(results_type) :: res
   type(post_processing_list) :: pproc
   type(integral_type) :: ints
   type(cache_list) :: caches
   type(timer_type) :: timer
   character(len=:), allocatable :: label

   call get_structure(mol, "MB16-43", "01")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return

   label = "xtbml"
   call add_post_processing(pproc, mol, label, error)
   if (allocated(error)) return

   call ctx%set_partition(1, nparts, error)
   if (allocated(error)) return
   ctx%verbosity = 0
   call calc%set_partition(ctx%partition)

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, 300.0_wp)
   ! Call post-processing directly: the SCF driver rejects an unreduced
   ! partition before it reaches the independent xTB-ML restriction.
   allocate(res%dict)
   call pproc%compute(mol, wfn, ints, calc, caches, 1.0_wp, ctx, timer, 0, res)

   if (.not.ctx%failed()) then
      call test_failed(error, "Post-processing of a distributed calculation was allowed")
      return
   end if
   call ctx%get_error(error)
   if (index(error%message, "xTB-ML") == 0) then
      deallocate(error)
      call test_failed(error, "Expected the xTB-ML partition restriction")
      return
   end if
   deallocate(error)

end subroutine test_post_processing


!> Summing the diatomic blocks of all parts has to reproduce the complete
!> integral and core Hamiltonian matrices
subroutine test_hamiltonian(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(work_partition) :: partition
   type(adjacency_list) :: list
   type(integral_type) :: full, part_ints, summed
   real(wp) :: cutoff
   real(wp), allocatable :: cn(:), selfenergy(:), lattr(:, :)
   integer :: part

   call get_structure(mol, "MB16-43", "01")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return

   allocate(cn(mol%nat), selfenergy(calc%bas%nsh))
   call calc%ncoord%get_cn(mol, cn)
   call get_selfenergy(calc%h0, mol%id, calc%bas%ish_at, calc%bas%nsh_id, cn=cn, &
      & selfenergy=selfenergy)

   cutoff = get_cutoff(calc%bas, 1.0_wp)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call new_adjacency_list(list, mol, lattr, cutoff)

   call new_integral(full, calc%bas%nao)
   call get_hamiltonian(mol, lattr, list, calc%bas, calc%h0, selfenergy, &
      & full%overlap, full%dipole, full%quadrupole, full%hamiltonian)

   call new_integral(summed, calc%bas%nao)
   call new_integral(part_ints, calc%bas%nao)
   summed%overlap(:, :) = 0.0_wp
   summed%hamiltonian(:, :) = 0.0_wp
   summed%dipole(:, :, :) = 0.0_wp
   summed%quadrupole(:, :, :) = 0.0_wp

   do part = 0, nparts - 1
      call new_work_partition(error, partition, part, nparts)
      if (allocated(error)) return
      call get_hamiltonian(mol, lattr, list, calc%bas, calc%h0, selfenergy, &
         & part_ints%overlap, part_ints%dipole, part_ints%quadrupole, &
         & part_ints%hamiltonian, partition)
      summed%overlap(:, :) = summed%overlap + part_ints%overlap
      summed%hamiltonian(:, :) = summed%hamiltonian + part_ints%hamiltonian
      summed%dipole(:, :, :) = summed%dipole + part_ints%dipole
      summed%quadrupole(:, :, :) = summed%quadrupole + part_ints%quadrupole
   end do

   if (sum(abs(full%overlap)) < 1.0e-6_wp) then
      call test_failed(error, "Overlap reference is empty")
      return
   end if

   if (any(abs(summed%overlap - full%overlap) > thr) &
      & .or. any(abs(summed%hamiltonian - full%hamiltonian) > thr) &
      & .or. any(abs(summed%dipole - full%dipole) > thr) &
      & .or. any(abs(summed%quadrupole - full%quadrupole) > thr)) then
      call test_failed(error, "Partitioned integrals do not match")
   end if

end subroutine test_hamiltonian


!> The Born interaction matrix and the solvent accessible surface are partitioned,
!> the Born radii themselves are evaluated for the full system on every part
subroutine test_solvation(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc

   call get_structure(mol, "MB16-43", "04")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return
   call add_solvation(calc, mol, error)
   if (allocated(error)) return

   call check_calculator(error, mol, calc)

end subroutine test_solvation


!> Attach an ALPB and a CDS container to the calculator
subroutine add_solvation(calc, mol, error)

   !> Single-point calculator
   type(xtb_calculator), intent(inout) :: calc

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(solvation_input) :: input
   class(solvation_type), allocatable :: solv
   class(container_type), allocatable :: cont

   input%alpb = alpb_input(80.2_wp, solvent="water", alpb=.true.)
   call new_solvation(solv, mol, input, error, "gfn2")
   if (allocated(error)) return
   call move_alloc(solv, cont)
   call calc%push_back(cont)

   input%cds = cds_input(alpb=.true., solvent="water")
   call new_solvation_cds(solv, mol, input, error, "gfn2")
   if (allocated(error)) return
   call move_alloc(solv, cont)
   call calc%push_back(cont)

end subroutine add_solvation


subroutine test_gfn1_mol(error)
   type(error_type), allocatable, intent(out) :: error
   call check_partitioned(error, "MB16-43", "02", .false.)
end subroutine test_gfn1_mol


subroutine test_gfn2_mol(error)
   type(error_type), allocatable, intent(out) :: error
   call check_partitioned(error, "MB16-43", "03", .true.)
end subroutine test_gfn2_mol


subroutine test_gfn1_pbc(error)
   type(error_type), allocatable, intent(out) :: error
   call check_partitioned(error, "X23", "urea", .false.)
end subroutine test_gfn1_pbc


subroutine test_gfn2_pbc(error)
   type(error_type), allocatable, intent(out) :: error
   call check_partitioned(error, "X23", "urea", .true.)
end subroutine test_gfn2_pbc


!> An external field is not an interaction loop, only the first part carries it
subroutine test_field(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(electric_field) :: efield
   class(container_type), allocatable :: cont

   call get_structure(mol, "MB16-43", "01")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return

   call new_electric_field(efield, [1.0e-2_wp, -2.0e-2_wp, 3.0e-2_wp])
   cont = efield
   call calc%push_back(cont)

   call check_calculator(error, mol, calc)

end subroutine test_field


!> Summing all parts has to reproduce the complete evaluation of all containers
subroutine check_partitioned(error, set, name, gfn2)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Structure set to test with
   character(len=*), intent(in) :: set

   !> Structure record to test with
   character(len=*), intent(in) :: name

   !> Whether to use GFN2-xTB instead of GFN1-xTB
   logical, intent(in) :: gfn2

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   call get_structure(mol, set, name)
   if (gfn2) then
      call new_gfn2_calculator(calc, mol, error)
   else
      call new_gfn1_calculator(calc, mol, error)
   end if
   if (allocated(error)) return

   call check_calculator(error, mol, calc)

end subroutine check_partitioned


!> One summation/check path for molecular, periodic, field and solvation cases.
subroutine check_calculator(error, mol, calc)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type), intent(in) :: mol
   type(xtb_calculator), intent(inout) :: calc
   type(work_partition) :: partition
   real(wp), allocatable :: full(:), actual(:), summed(:)
   integer :: part

   call evaluate(mol, calc, full)
   call check(error, maxval(abs(full)) > thr, "Empty reference")
   if (allocated(error)) return
   allocate(summed(size(full)), source=0.0_wp)
   do part = 0, nparts-1
      call new_work_partition(error, partition, part, nparts)
      if (allocated(error)) return
      call calc%set_partition(partition)
      call evaluate(mol, calc, actual)
      summed(:) = summed + actual
   end do
   call check(error, maxval(abs(summed-full)), 0.0_wp, thr=thr)
   if (allocated(error)) return

   call calc%set_partition(serial_work_partition)
   call evaluate(mol, calc, actual)
   call check(error, maxval(abs(actual-full)), 0.0_wp, thr=thr)
end subroutine check_calculator


!> Pack energies, derivatives and all potential components for one comparison.
subroutine evaluate(mol, calc, output)
   type(structure_type), intent(in) :: mol
   type(xtb_calculator), intent(in) :: calc
   real(wp), allocatable, intent(out) :: output(:)
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   real(wp) :: energies(mol%nat), gradient(3, mol%nat), sigma(3, 3)

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, 300.0_wp)
   call model_wavefunction(wfn, mol, calc)
   call new_potential(pot, mol, calc%bas, 1)
   call pot%reset()
   energies(:) = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   if (allocated(calc%repulsion)) call run(calc%repulsion)
   if (allocated(calc%halogen)) call run(calc%halogen)
   if (allocated(calc%dispersion)) call run(calc%dispersion)
   if (allocated(calc%coulomb)) call run(calc%coulomb)
   if (allocated(calc%interactions)) call run(calc%interactions)
   output = [energies, reshape(gradient, [size(gradient)]), reshape(sigma, [9]), &
      & reshape(pot%vat, [size(pot%vat)]), reshape(pot%vsh, [size(pot%vsh)]), &
      & reshape(pot%vdp, [size(pot%vdp)]), reshape(pot%vqp, [size(pot%vqp)])]

contains

   subroutine run(cont)
      class(container_type), intent(in) :: cont
      type(container_cache) :: cache

      call cont%update(mol, cache)
      call cont%get_engrad(mol, cache, energies, gradient, sigma)
      call cont%get_energy(mol, cache, wfn, energies)
      call cont%get_potential(mol, cache, wfn, pot)
      call cont%get_gradient(mol, cache, wfn, gradient, sigma)
   end subroutine run

end subroutine evaluate


!> Deterministic charge distribution to probe the selfconsistent contributions
subroutine model_wavefunction(wfn, mol, calc)

   !> Wavefunction data
   type(wavefunction_type), intent(inout) :: wfn

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc

   integer :: iat, ish, ii, k

   do iat = 1, mol%nat
      wfn%qat(iat, 1) = 0.1_wp * sin(real(iat, wp))
      do k = 1, 3
         wfn%dpat(k, iat, 1) = 0.05_wp * cos(real(iat + k, wp))
      end do
      do k = 1, 5
         wfn%qpat(k, iat, 1) = 0.02_wp * sin(real(2*iat + k, wp))
      end do
      ! the anisotropic terms assume traceless quadrupole moments
      wfn%qpat(6, iat, 1) = -wfn%qpat(1, iat, 1) - wfn%qpat(3, iat, 1)
      ii = calc%bas%ish_at(iat)
      do ish = 1, calc%bas%nsh_at(iat)
         wfn%qsh(ii+ish, 1) = wfn%qat(iat, 1) / real(calc%bas%nsh_at(iat), wp)
      end do
   end do

end subroutine model_wavefunction


end module test_partition
