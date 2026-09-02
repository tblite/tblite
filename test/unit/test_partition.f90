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
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_basis_type, only : get_cutoff
   use tblite_container, only : container_cache, container_type
   use tblite_context, only : context_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_external_field, only : electric_field, new_electric_field
   use tblite_features, only : get_tblite_feature, tblite_has_mpi
   use tblite_integral_type, only : integral_type, new_integral
   use tblite_mpi_utils, only : get_mpi_comm_world, new_mpi_work_partition, &
      & mpi_allreduce_sum
   use tblite_partition, only : work_partition, new_work_partition, &
      & owns_index, owns_pair
   use tblite_post_processing_list, only : post_processing_list, add_post_processing
   use tblite_results, only : results_type
   use tblite_scf_potential, only : potential_type, new_potential
   use tblite_solvation, only : solvation_input, solvation_type, alpb_input, cds_input, &
      & new_solvation, new_solvation_cds
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

   !> Collected output of all containers of a calculator
   type :: container_output
      real(wp), allocatable :: energies(:)
      real(wp), allocatable :: gradient(:, :)
      real(wp), allocatable :: sigma(:, :)
      real(wp), allocatable :: vat(:, :)
      real(wp), allocatable :: vsh(:, :)
      real(wp), allocatable :: vdp(:, :, :)
      real(wp), allocatable :: vqp(:, :, :)
   end type container_output

contains


!> Collect all exported unit tests
subroutine collect_partition(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("invalid", test_invalid), &
      new_unittest("context", test_context), &
      new_unittest("disjoint-index", test_disjoint_index), &
      new_unittest("disjoint-pair", test_disjoint_pair), &
      new_unittest("absent", test_absent), &
      new_unittest("feature", test_feature), &
      new_unittest("mpi-unavailable", test_mpi_unavailable), &
      new_unittest("mpi-mismatch", test_mpi_mismatch), &
      new_unittest("unreduced", test_unreduced), &
      new_unittest("post-processing", test_post_processing), &
      new_unittest("serial", test_serial), &
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

   do icase = 1, size(invalid, 2)
      call new_work_partition(partition_error, partition, &
         & invalid(1, icase), invalid(2, icase))
      if (.not.allocated(partition_error)) then
         call test_failed(error, "Invalid work partition was accepted")
         return
      end if
      deallocate(partition_error)
   end do

   call new_work_partition(partition_error, partition, 0, 1)
   if (allocated(partition_error)) then
      call test_failed(error, "Serial work partition was rejected")
   end if

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
   end if

end subroutine test_context


!> Every one-dimensional unit of work belongs to exactly one part
subroutine test_disjoint_index(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(work_partition) :: partition
   integer :: idx, part, owner

   do idx = 1, 50
      owner = 0
      do part = 0, nparts - 1
         call new_work_partition(error, partition, part, nparts)
         if (allocated(error)) return
         if (owns_index(partition, idx)) owner = owner + 1
      end do
      call check(error, owner, 1)
      if (allocated(error)) return
   end do

end subroutine test_disjoint_index


!> Every atom pair belongs to exactly one part
subroutine test_disjoint_pair(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(work_partition) :: partition
   integer :: iat, jat, part, owner

   do iat = 1, 12
      do jat = 1, iat
         owner = 0
         do part = 0, nparts - 1
            call new_work_partition(error, partition, part, nparts)
            if (allocated(error)) return
            if (owns_pair(partition, iat, jat)) owner = owner + 1
         end do
         call check(error, owner, 1)
         if (allocated(error)) return
      end do
   end do

end subroutine test_disjoint_pair


!> An absent partition selects the complete work
subroutine test_absent(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   if (.not.owns_index(idx=7)) then
      call test_failed(error, "Absent partition does not own an index")
      return
   end if

   if (.not.owns_pair(iat=5, jat=3)) then
      call test_failed(error, "Absent partition does not own a pair")
      return
   end if

   if (.not.owns_pair(work_partition(), 5, 3)) then
      call test_failed(error, "Serial partition does not own a pair")
   end if

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
   type(work_partition) :: partition
   type(error_type), allocatable :: mpi_error
   real(wp) :: val(1)

   ! this tester is not an MPI program and must not call into a live library
   if (tblite_has_mpi) return

   call ctx%set_mpi(mpi_error)
   if (.not.allocated(mpi_error)) then
      call test_failed(error, "Uninitialized MPI was accepted by the context")
      return
   end if

   if (allocated(ctx%comm)) then
      call test_failed(error, "Context enabled MPI without a usable library")
      return
   end if

   deallocate(mpi_error)
   call new_mpi_work_partition(mpi_error, partition, get_mpi_comm_world())
   if (.not.allocated(mpi_error)) then
      call test_failed(error, "Uninitialized MPI produced a work partition")
      return
   end if

   deallocate(mpi_error)
   val = 1.0_wp
   call mpi_allreduce_sum(mpi_error, val, get_mpi_comm_world())
   if (.not.allocated(mpi_error)) then
      call test_failed(error, "Uninitialized MPI performed a reduction")
   end if

end subroutine test_mpi_unavailable


!> Reducing results of a calculator that does not share the partition of the
!> context would double count every contribution
subroutine test_mpi_mismatch(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   real(wp) :: energy

   call get_structure(mol, "MB16-43", "01")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return

   call ctx%set_partition(1, nparts, error)
   if (allocated(error)) return
   ctx%comm = 0
   ctx%verbosity = 0

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, 300.0_wp)
   call xtb_singlepoint(ctx, mol, calc, wfn, 1.0_wp, energy)

   if (.not.ctx%failed()) then
      call test_failed(error, "Mismatched work partition was not reported")
      return
   end if
   call ctx%get_error(error)
   deallocate(error)

end subroutine test_mpi_mismatch


!> A partitioned calculator without a context reducing the partial results
!> would converge the SCF against an incomplete potential
subroutine test_unreduced(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(work_partition) :: partition
   type(wavefunction_type) :: wfn
   real(wp) :: energy

   call get_structure(mol, "MB16-43", "01")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return

   call new_work_partition(error, partition, 1, nparts)
   if (allocated(error)) return
   call calc%set_partition(partition)
   ctx%verbosity = 0

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, 300.0_wp)
   call xtb_singlepoint(ctx, mol, calc, wfn, 1.0_wp, energy)

   if (.not.ctx%failed()) then
      call test_failed(error, "Unreduced work partition was not reported")
      return
   end if
   call ctx%get_error(error)
   deallocate(error)

end subroutine test_unreduced


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
   real(wp) :: energy
   character(len=:), allocatable :: label

   call get_structure(mol, "MB16-43", "01")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return

   label = "xtbml"
   call add_post_processing(pproc, mol, label, error)
   if (allocated(error)) return

   ! the context reduces, so the calculator has to carry the same partition
   call ctx%set_partition(1, nparts, error)
   if (allocated(error)) return
   ctx%comm = 0
   ctx%verbosity = 0
   call calc%set_partition(ctx%partition)

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, 300.0_wp)
   call xtb_singlepoint(ctx, mol, calc, wfn, 1.0_wp, energy, results=res, &
      & post_process=pproc)

   if (.not.ctx%failed()) then
      call test_failed(error, "Post-processing of a distributed calculation was allowed")
      return
   end if
   call ctx%get_error(error)
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
   type(work_partition) :: partition
   type(container_output) :: full, summed, part_result
   integer :: part

   call get_structure(mol, "MB16-43", "04")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return
   call add_solvation(calc, mol, error)
   if (allocated(error)) return

   call evaluate(mol, calc, full)

   do part = 0, nparts - 1
      call new_work_partition(error, partition, part, nparts)
      if (allocated(error)) return
      call calc%set_partition(partition)
      call evaluate(mol, calc, part_result)
      if (part == 0) then
         summed = part_result
      else
         call accumulate(summed, part_result)
      end if
   end do

   call compare(error, summed, full, "Solvation")

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


!> Passing the serial partition must be identical to leaving it at its default
subroutine test_serial(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(container_output) :: full, serial

   call get_structure(mol, "MB16-43", "01")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return

   call evaluate(mol, calc, full)

   call calc%set_partition(work_partition())
   call evaluate(mol, calc, serial)

   call compare(error, serial, full, "Serial partition")

end subroutine test_serial


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
   type(work_partition) :: partition
   type(electric_field) :: efield
   class(container_type), allocatable :: cont
   type(container_output) :: full, summed, part_result
   integer :: part

   call get_structure(mol, "MB16-43", "01")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return

   call new_electric_field(efield, [1.0e-2_wp, -2.0e-2_wp, 3.0e-2_wp])
   cont = efield
   call calc%push_back(cont)

   call evaluate(mol, calc, full)

   do part = 0, nparts - 1
      call new_work_partition(error, partition, part, nparts)
      if (allocated(error)) return
      call calc%set_partition(partition)
      call evaluate(mol, calc, part_result)
      if (part == 0) then
         summed = part_result
      else
         call accumulate(summed, part_result)
      end if
   end do

   call compare(error, summed, full, "Electric field")

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
   type(work_partition) :: partition
   type(container_output) :: full, summed, part_result
   integer :: part

   call get_structure(mol, set, name)
   if (gfn2) then
      call new_gfn2_calculator(calc, mol, error)
   else
      call new_gfn1_calculator(calc, mol, error)
   end if
   if (allocated(error)) return

   call evaluate(mol, calc, full)

   do part = 0, nparts - 1
      call new_work_partition(error, partition, part, nparts)
      if (allocated(error)) return
      call calc%set_partition(partition)
      call evaluate(mol, calc, part_result)
      if (part == 0) then
         summed = part_result
      else
         call accumulate(summed, part_result)
      end if
   end do

   call compare(error, summed, full, "Partitioned containers")

end subroutine check_partitioned


!> Evaluate every container of the calculator for a fixed model wavefunction
subroutine evaluate(mol, calc, output)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc

   !> Collected container contributions
   type(container_output), intent(out) :: output

   type(wavefunction_type) :: wfn
   type(potential_type) :: pot

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, 300.0_wp)
   call model_wavefunction(wfn, mol, calc)

   call new_potential(pot, mol, calc%bas, 1)
   call pot%reset()

   allocate(output%energies(mol%nat), source=0.0_wp)
   allocate(output%gradient(3, mol%nat), source=0.0_wp)
   allocate(output%sigma(3, 3), source=0.0_wp)

   if (allocated(calc%repulsion)) call run(calc%repulsion, mol, wfn, pot, output)
   if (allocated(calc%halogen)) call run(calc%halogen, mol, wfn, pot, output)
   if (allocated(calc%dispersion)) call run(calc%dispersion, mol, wfn, pot, output)
   if (allocated(calc%coulomb)) call run(calc%coulomb, mol, wfn, pot, output)
   if (allocated(calc%interactions)) call run(calc%interactions, mol, wfn, pot, output)

   output%vat = pot%vat
   output%vsh = pot%vsh
   output%vdp = pot%vdp
   output%vqp = pot%vqp

end subroutine evaluate


!> Accumulate all contributions of a single container
subroutine run(cont, mol, wfn, pot, output)

   !> Interaction container to evaluate
   class(container_type), intent(in) :: cont

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn

   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   !> Collected container contributions
   type(container_output), intent(inout) :: output

   type(container_cache) :: cache

   call cont%update(mol, cache)
   call cont%get_engrad(mol, cache, output%energies, output%gradient, output%sigma)
   call cont%get_energy(mol, cache, wfn, output%energies)
   call cont%get_potential(mol, cache, wfn, pot)
   call cont%get_gradient(mol, cache, wfn, output%gradient, output%sigma)

end subroutine run


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


!> Add the contributions of one part to the running sum
subroutine accumulate(summed, part_result)

   !> Running sum over all parts
   type(container_output), intent(inout) :: summed

   !> Contributions of a single part
   type(container_output), intent(in) :: part_result

   summed%energies(:) = summed%energies + part_result%energies
   summed%gradient(:, :) = summed%gradient + part_result%gradient
   summed%sigma(:, :) = summed%sigma + part_result%sigma
   summed%vat(:, :) = summed%vat + part_result%vat
   summed%vsh(:, :) = summed%vsh + part_result%vsh
   summed%vdp(:, :, :) = summed%vdp + part_result%vdp
   summed%vqp(:, :, :) = summed%vqp + part_result%vqp

end subroutine accumulate


!> Compare the summed parts against the complete evaluation
subroutine compare(error, actual, expected, label)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Summed contributions of all parts
   type(container_output), intent(in) :: actual

   !> Contributions of the complete calculation
   type(container_output), intent(in) :: expected

   !> Label of the compared quantity
   character(len=*), intent(in) :: label

   ! a vanishing reference would make every comparison below pass trivially
   if (sum(abs(expected%energies)) < 1.0e-6_wp) then
      call test_failed(error, label//" reference is empty")
      return
   end if

   call check(error, sum(actual%energies), sum(expected%energies), thr=thr)
   if (allocated(error)) return

   if (any(abs(actual%energies - expected%energies) > thr)) then
      call test_failed(error, label//" energies do not match")
      return
   end if

   if (any(abs(actual%gradient - expected%gradient) > thr)) then
      call test_failed(error, label//" gradient does not match")
      return
   end if

   if (any(abs(actual%sigma - expected%sigma) > thr)) then
      call test_failed(error, label//" virial does not match")
      return
   end if

   if (any(abs(actual%vat - expected%vat) > thr) &
      & .or. any(abs(actual%vsh - expected%vsh) > thr) &
      & .or. any(abs(actual%vdp - expected%vdp) > thr) &
      & .or. any(abs(actual%vqp - expected%vqp) > thr)) then
      call test_failed(error, label//" potential does not match")
   end if

end subroutine compare


end module test_partition
