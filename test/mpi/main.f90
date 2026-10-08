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

!> Distributing a single point calculation over MPI ranks must reproduce the
!> serial result on every rank.
module column_test_solver
   use mctc_env, only : wp
   use tblite_context_solver, only : context_solver
   use tblite_lapack_scalapack, only : psygvd_solver, new_psygvd
   use tblite_lapack_solver, only : lapack_solver
   use tblite_scf_solver, only : solver_type
   implicit none
   private
   public :: column_factory

   type, extends(lapack_solver) :: column_factory
   contains
      procedure :: new => new_column_solver
      procedure :: distributed_columns => use_columns
   end type column_factory
contains
   logical function use_columns(self, n, comm) result(distribute)
      class(column_factory), intent(in) :: self
      integer, intent(in) :: n
      integer, intent(in), optional :: comm
      distribute = present(comm)
   end function use_columns
   subroutine new_column_solver(self, solver, overlap, nel, kt, comm)
      class(column_factory), intent(inout) :: self
      class(solver_type), allocatable, intent(out) :: solver
      real(wp), intent(in) :: overlap(:, :), nel(:), kt
      integer, intent(in), optional :: comm
      type(psygvd_solver), allocatable :: tmp
      allocate(tmp)
      call new_psygvd(tmp, overlap, nel, kt, comm)
      call move_alloc(tmp, solver)
   end subroutine new_column_solver
end module column_test_solver

program test_mpi_singlepoint
   use column_test_solver, only : column_factory
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type, new
   use mpi_f08, only : MPI_COMM_WORLD, MPI_Abort, MPI_Comm_rank, MPI_Comm_size, &
      & MPI_Finalize, MPI_Init
   use mstore, only : get_structure
   use tblite_ceh_ceh, only : new_ceh_calculator
   use tblite_ceh_singlepoint, only : ceh_singlepoint
   use tblite_container, only : container_type
   use tblite_context, only : context_type
   use tblite_features, only : tblite_has_scalapack
   use tblite_integral_type, only : integral_type, new_integral
   use tblite_lapack_scalapack, only : psygvd_solver, new_psygvd, &
      & distribute_diagonalization
   use tblite_lapack_sygvd, only : sygvd_solver, new_sygvd
   use tblite_mpi_utils, only : mpi_allreduce_sum, mpi_density_columns, &
      & mpi_transpose_columns, mpi_gather_columns, mpi_expand_matrix, mpi_multiply_columns, mpi_expand_integrals
   use tblite_partition, only : work_partition, serial_work_partition, operator(==), column_range
   use tblite_post_processing_list, only : post_processing_list, add_post_processing
   use tblite_results, only : results_type
   use tblite_scf_potential, only : potential_type, new_potential, add_pot_to_h1, &
      & mpi_add_pot_to_h1, reduce_potential
   use tblite_solvation, only : solvation_input, solvation_type, alpb_input, cds_input, &
      & new_solvation, new_solvation_cds
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none

   real(wp), parameter :: thr = 1.0e-8_wp

   type(structure_type) :: mol
   type(error_type), allocatable :: error
   real(wp) :: serial_energy, mpi_energy
   real(wp) :: serial_sigma(3, 3), mpi_sigma(3, 3)
   real(wp), allocatable :: serial_gradient(:, :), mpi_gradient(:, :)
   real(wp), allocatable :: serial_wbo(:, :, :), mpi_wbo(:, :, :)
   integer :: stat, nspin

   call MPI_Init(stat)

   call get_structure(mol, "MB16-43", "01")
   allocate(serial_gradient(3, mol%nat), mpi_gradient(3, mol%nat))

   do nspin = 1, 2
      call run(mol, .false., nspin, serial_energy, serial_gradient, serial_sigma, serial_wbo, error)
      call check_error(error)
      call run(mol, .true., nspin, mpi_energy, mpi_gradient, mpi_sigma, mpi_wbo, error)
      call check_error(error)

      call assert(abs(mpi_energy - serial_energy) < thr, "energy")
      call assert(all(abs(mpi_gradient - serial_gradient) < thr), "gradient")
      call assert(all(abs(mpi_sigma - serial_sigma) < thr), "virial")
      call assert(allocated(serial_wbo) .and. allocated(mpi_wbo), "bond order results")
      call assert(all(abs(mpi_wbo - serial_wbo) < thr), "bond orders")
      if (nspin == 1) then
         call run(mol, .true., nspin, mpi_energy, mpi_gradient, mpi_sigma, mpi_wbo, error, retain=.false.)
         call check_error(error)
         call assert(all(abs(mpi_wbo-serial_wbo) < thr), "post-processing before matrix release")
      end if
      if (tblite_has_scalapack) then
         call run(mol, .true., nspin, mpi_energy, mpi_gradient, mpi_sigma, mpi_wbo, error, .true.)
         call check_error(error)
         call assert(abs(mpi_energy-serial_energy) < thr, "column SCF energy")
         call assert(all(abs(mpi_gradient-serial_gradient) < thr), "column SCF gradient")
         call assert(all(abs(mpi_sigma-serial_sigma) < thr), "column SCF virial")
         call assert(all(abs(mpi_wbo-serial_wbo) < thr), "column SCF bond orders")
         if (nspin == 1) then
            call run(mol, .true., nspin, mpi_energy, mpi_gradient, mpi_sigma, mpi_wbo, error, .true., .true.)
            call check_error(error)
            call assert(all(abs(mpi_wbo-serial_wbo) < thr), "dense post-processing compatibility")
            call run(mol, .true., nspin, mpi_energy, mpi_gradient, mpi_sigma, mpi_wbo, error, &
               & columns=.true., retain=.false.)
            call check_error(error)
            call assert(all(abs(mpi_wbo-serial_wbo) < thr), "column post-processing before matrix release")
            call run(mol, .true., nspin, mpi_energy, mpi_gradient, mpi_sigma, mpi_wbo, error, &
               & columns=.true., full_post=.true., retain=.false.)
            call check_error(error)
            call assert(all(abs(mpi_wbo-serial_wbo) < thr), "dense post-processing before matrix release")
         end if
      end if
   end do

   call check_sync_error()
   call check_ceh()
   call check_eigensolver(37)
   call check_eigensolver(193)
   call check_column_product(5)
   call check_column_product(1)
   call check_column_product(0)
   call check_potential_reduction()
   call check_large_reduction()
   call check_hamiltonian(mol, 1)
   call check_hamiltonian(mol, 2)
   block
      type(structure_type) :: atom
      call new(atom, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call check_hamiltonian(atom, 2)
   end block
   call check_density(37, 1)
   call check_density(37, 2)

   call MPI_Finalize(stat)

contains

   !> Multiple communication chunks, including an incomplete last chunk.
   subroutine check_large_reduction()
      real(wp), allocatable :: values(:)
      integer :: i, rank, nranks

      call MPI_Comm_rank(MPI_COMM_WORLD, rank, stat)
      call MPI_Comm_size(MPI_COMM_WORLD, nranks, stat)
      allocate(values(2**22 + 17))
      do i = 1, size(values)
         values(i) = real(modulo(i, 97) + rank, wp)
      end do
      call mpi_allreduce_sum(error, values, MPI_COMM_WORLD%MPI_VAL)
      call check_error(error)
      do i = 1, size(values)
         if (values(i) /= real(nranks*modulo(i, 97) + nranks*(nranks-1)/2, wp)) then
            call assert(.false., "large real reduction")
         end if
      end do
      deallocate(values)
      allocate(values(0))
      call mpi_allreduce_sum(error, values, MPI_COMM_WORLD%MPI_VAL)
      call check_error(error)
   end subroutine check_large_reduction

   subroutine run(mol, distributed, nspin, energy, gradient, sigma, wbo, error, columns, full_post, retain)
      type(structure_type), intent(in) :: mol
      logical, intent(in) :: distributed
      logical, intent(in), optional :: columns, full_post, retain
      integer, intent(in) :: nspin
      real(wp), allocatable, intent(out) :: wbo(:, :, :)
      real(wp), intent(out) :: energy, gradient(:, :), sigma(:, :)
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
      type(post_processing_list) :: post
      type(results_type) :: results
      character(len=:), allocatable :: label
      real(wp), allocatable :: csc(:, :)
      integer :: spin, i

      ctx%verbosity = 0
      if (present(columns)) then
         if (columns) ctx%solver = column_factory()
      end if
      call new_gfn2_calculator(calc, mol, error)
      if (allocated(error)) return
      call add_solvation(calc, mol, error)
      if (allocated(error)) return

      if (distributed) then
         call ctx%set_mpi(error)
         if (allocated(error)) return
         call calc%set_partition(ctx%partition)
      end if

      calc%save_integrals = .true.
      if (present(retain)) calc%retain_matrices = retain
      label = "bond-orders"
      call add_post_processing(post, mol, label, error)
      if (allocated(error)) return
      if (present(full_post)) then
         ! Emulate an existing consumer that requires full AO matrices.
         if (full_post) post%list(1)%pproc%local_matrices = .false.
      end if
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, nspin, 300.0_wp)
      call xtb_singlepoint(ctx, mol, calc, wfn, 1.0_wp, energy, gradient, sigma, &
         & results=results, post_process=post)
      if (ctx%failed()) call ctx%get_error(error)
      if (allocated(error)) return
      call results%dict%get_entry("bond-orders", wbo)
      call assert(size(results%hamiltonian, 2) == calc%bas%nao, "saved full Hamiltonian")
      if (.not.calc%retain_matrices) then
         call assert(.not.allocated(wfn%coeff) .and. .not.allocated(wfn%density), "released AO matrices")
         return
      end if
      call assert(all(shape(wfn%coeff) == [calc%bas%nao, calc%bas%nao, nspin]), "full returned coefficients")
      call assert(all(shape(wfn%density) == shape(wfn%coeff)), "full returned density")
      do spin = 1, nspin
         csc = matmul(transpose(wfn%coeff(:, :, spin)), matmul(results%overlap, wfn%coeff(:, :, spin)))
         do i = 1, calc%bas%nao
            csc(i, i) = csc(i, i) - 1.0_wp
         end do
         call assert(all(abs(csc) < thr), "returned coefficient orthonormality")
      end do
   end subroutine run

   subroutine add_solvation(calc, mol, error)
      type(xtb_calculator), intent(inout) :: calc
      type(structure_type), intent(in) :: mol
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

   !> A failure on the last rank alone has to become visible to every rank,
   !> this call hangs instead of returning if the synchronization is missing
   subroutine check_sync_error()
      type(error_type), allocatable :: error
      type(context_type) :: ctx
      type(work_partition) :: partition
      real(wp) :: val(1)
      integer :: rank, nranks

      call MPI_Comm_rank(MPI_COMM_WORLD, rank, stat)
      call MPI_Comm_size(MPI_COMM_WORLD, nranks, stat)
      call ctx%set_mpi(error, MPI_COMM_WORLD%MPI_VAL)
      call check_error(error)
      partition = ctx%partition
      call ctx%set_partition(nranks, nranks, error)
      call assert(allocated(error), "invalid partition")
      call assert(allocated(ctx%comm) .and. partition == ctx%partition, "preserved context")

      ! One mismatched calculator must reject the SCF collectively.
      if (nranks > 1) then
         if (rank == nranks-1) partition = serial_work_partition
         call ctx%check_partition(partition, error)
         call assert(allocated(error), "partition mismatch synchronization")
         deallocate(error)
      end if

      if (rank == nranks-1) then
         call fatal_error(error, "Failure on the last rank")
         call ctx%set_error(error)
      end if
      call ctx%sync_error()
      call assert(ctx%failed(), "error synchronization")
      call ctx%get_error(error)
      if (rank == nranks-1) call assert(error%message == "Failure on the last rank", &
         & "original error message")
      val = 1.0_wp
      call mpi_allreduce_sum(error, val, ctx%comm)
      call assert(all(val == 1.0_wp), "pending error skips reduction")
      call ctx%sync_error()
      call assert(.not.ctx%failed(), "error-free synchronization")

      call ctx%set_partition(0, 1, error)
      call check_error(error)
      call assert(.not.allocated(ctx%comm), "external partition clears MPI")
   end subroutine check_sync_error

   !> The CEH charges of a distributed run have to match the serial ones
   subroutine check_ceh()
      real(wp), allocatable :: serial_qat(:), mpi_qat(:)

      call run_ceh(mol, .false., serial_qat, error)
      call check_error(error)
      call run_ceh(mol, .true., mpi_qat, error)
      call check_error(error)

      call assert(all(abs(mpi_qat - serial_qat) < thr), "CEH charges")
   end subroutine check_ceh

   subroutine run_ceh(mol, distributed, qat, error)
      type(structure_type), intent(in) :: mol
      logical, intent(in) :: distributed
      real(wp), allocatable, intent(out) :: qat(:)
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn

      ctx%verbosity = 0
      call new_ceh_calculator(calc, mol, error)
      if (allocated(error)) return

      if (distributed) then
         call ctx%set_mpi(error)
         if (allocated(error)) return
         call calc%set_partition(ctx%partition)
      end if

      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, 4000.0_wp)
      call ceh_singlepoint(ctx, calc, mol, wfn, 1.0_wp)
      if (ctx%failed()) then
         call ctx%get_error(error)
         return
      end if
      qat = wfn%qat(:, 1)
   end subroutine run_ceh

   !> The distributed eigensolver has to reproduce the eigenvalues of the
   !> replicated one for the block sizes and process grids it picks
   subroutine check_eigensolver(n)
      integer, intent(in) :: n
      real(wp) :: amat(n, n), bmat(n, n), work(n, n)
      real(wp) :: hmat(n, n), eval(n), reference(n), nel(2)
      type(sygvd_solver) :: serial
      type(psygvd_solver) :: distributed
      integer :: i, j, iteration

      if (.not.tblite_has_scalapack) return

      ! a diagonally dominant pair keeps the generalized problem well conditioned
      do i = 1, n
         do j = 1, n
            work(j, i) = sin(real(i*j, wp))
         end do
      end do
      amat(:, :) = 0.5_wp*(work + transpose(work))
      bmat(:, :) = 0.0_wp
      do i = 1, n
         amat(i, i) = amat(i, i) + real(i, wp)
         bmat(i, i) = 1.0_wp + 0.1_wp*real(modulo(i, 3), wp)
      end do

      nel = [real(n, wp), real(n, wp)]
      call new_sygvd(serial, bmat, nel, 0.0_wp)
      call new_psygvd(distributed, bmat, nel, 0.0_wp, MPI_COMM_WORLD%MPI_VAL)
      do iteration = 1, 3
         ! Reuse the grid for a different Hamiltonian, as in successive SCF steps.
         amat(1, 1) = amat(1, 1) + 0.1_wp
         ! Only one rank owns this entry; all ranks must refactor together.
         if (iteration == 3) bmat(n, n) = bmat(n, n) + 0.05_wp
         hmat(:, :) = amat
         call serial%solve(hmat, bmat, reference, error)
         call check_error(error)
         hmat(:, :) = amat
         call distributed%solve(hmat, bmat, eval, error)
         call check_error(error)
         call assert(all(abs(eval - reference) < thr), "eigenvalues")
         work(:, :) = matmul(amat, hmat) - matmul(bmat, hmat)*spread(eval, 1, n)
         call assert(maxval(abs(work)) < thr, "eigenvector residual")
         work(:, :) = matmul(transpose(hmat), matmul(bmat, hmat))
         do i = 1, n
            work(i, i) = work(i, i) - 1.0_wp
         end do
         call assert(maxval(abs(work)) < thr, "eigenvector normalization")
      end do
      call distributed%delete()

      ! a matrix this small is faster on a single rank
      call assert(.not.distribute_diagonalization(4, MPI_COMM_WORLD%MPI_VAL), &
         & "small matrix fallback")
   end subroutine check_eigensolver

   subroutine check_column_product(n)
      integer, intent(in) :: n
      real(wp) :: a(n, n), b(n, n), reference(n, n), gathered(n, n), weights(n)
      real(wp), allocatable :: local(:, :), product(:, :)
      integer :: i, j, rank, nranks, columns(2)
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, stat)
      call MPI_Comm_size(MPI_COMM_WORLD, nranks, stat)
      columns = column_range(n, rank, nranks)
      do j = 1, n
         do i = 1, n
            a(i, j) = sin(real(i+2*j, wp))
            b(i, j) = cos(real(3*i+j, wp))
         end do
      end do
      local = a(:, columns(1):columns(2))
      allocate(product(n, size(local, 2)))
      call mpi_transpose_columns(error, local, product, MPI_COMM_WORLD%MPI_VAL)
      call check_error(error)
      call mpi_gather_columns(error, product, gathered, n, MPI_COMM_WORLD%MPI_VAL)
      call check_error(error)
      call assert(all(abs(gathered-transpose(a)) < thr), "column transpose including empty parts")
      call mpi_multiply_columns(error, local, b(:, columns(1):columns(2)), product, MPI_COMM_WORLD%MPI_VAL)
      call mpi_gather_columns(error, product, gathered, n, MPI_COMM_WORLD%MPI_VAL)
      call check_error(error)
      reference = matmul(a, b)
      call assert(all(abs(gathered-reference) < thr), "column product including empty parts")
      weights = [(cos(real(i, wp)), i=1, n)]
      call mpi_density_columns(error, weights, local, product, MPI_COMM_WORLD%MPI_VAL)
      call mpi_gather_columns(error, product, gathered, n, MPI_COMM_WORLD%MPI_VAL)
      call check_error(error)
      reference = matmul(a*spread(weights, 1, n), transpose(a))
      call assert(all(abs(gathered-reference) < thr), "weighted column density including empty parts")
   end subroutine check_column_product

   !> Pack/unpack every potential component and spin, reusing the buffer.
   subroutine check_potential_reduction()
      type(xtb_calculator) :: calc
      type(potential_type) :: pot
      integer :: rank, nranks, step
      real(wp) :: value, total

      call new_gfn2_calculator(calc, mol, error)
      call check_error(error)
      call new_potential(pot, mol, calc%bas, 2)
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, stat)
      call MPI_Comm_size(MPI_COMM_WORLD, nranks, stat)
      do step = 1, 2
         value = real((rank+1)*step, wp)
         total = real(nranks*(nranks+1)*step, wp)/2.0_wp
         pot%vat = value
         pot%vsh = -2*value
         pot%vdp = 3*value
         pot%vqp = -4*value
         call reduce_potential(error, pot, MPI_COMM_WORLD%MPI_VAL)
         call check_error(error)
         call assert(all(abs(pot%vat-total) < thr), "packed atom potential")
         call assert(all(abs(pot%vsh+2*total) < thr), "packed shell potential")
         call assert(all(abs(pot%vdp-3*total) < thr), "packed dipole potential")
         call assert(all(abs(pot%vqp+4*total) < thr), "packed quadrupole potential")
      end do
   end subroutine check_potential_reduction

   !> Uneven/empty column ranges and both spin representations must match the
   !> serial symmetric expansion, including non-symmetric multipole integrals.
   subroutine check_hamiltonian(mol, nspin)
      type(structure_type), intent(in) :: mol
      integer, intent(in) :: nspin
      type(xtb_calculator) :: calc
      type(potential_type) :: pot, reference_pot
      type(integral_type) :: ints, local_ints
      type(potential_type) :: saved_pot
      real(wp), allocatable :: h(:, :, :), reference(:, :, :)
      integer :: i, j, k, spin, n, rank, nranks, columns(2)

      call new_gfn2_calculator(calc, mol, error)
      call check_error(error)
      n = calc%bas%nao
      call new_integral(ints, n)
      call new_potential(pot, mol, calc%bas, nspin)
      call pot%reset()
      do j = 1, n
         do i = 1, n
            ints%hamiltonian(i, j) = cos(real(i+j, wp))
            ints%overlap(i, j) = sin(real(i*j, wp))
            do k = 1, 3
               ints%dipole(k, i, j) = cos(real(i+2*j+3*k, wp))
            end do
            do k = 1, 6
               ints%quadrupole(k, i, j) = sin(real(3*i+j+2*k, wp))
            end do
         end do
      end do
      do spin = 1, nspin
         do i = 1, mol%nat
            pot%vat(i, spin) = 0.1_wp*real(i+spin, wp)
            pot%vdp(:, i, spin) = [0.1_wp, -0.3_wp, 0.2_wp]*real(i+spin, wp)
            pot%vqp(:, i, spin) = [(0.02_wp*real(i+spin+k, wp), k=1, 6)]
         end do
         pot%vsh(:, spin) = [(0.01_wp*real(i-spin, wp), i=1, calc%bas%nsh)]
      end do
      reference_pot = pot
      saved_pot = pot
      allocate(h(n, n, nspin), reference(n, n, nspin))
      call add_pot_to_h1(calc%bas, ints, reference_pot, reference)
      h = -123.0_wp
      call mpi_add_pot_to_h1(error, calc%bas, ints, pot, h, MPI_COMM_WORLD%MPI_VAL)
      call check_error(error)
      call assert(all(abs(h-reference) < thr), "Hamiltonian column expansion")
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, stat)
      call MPI_Comm_size(MPI_COMM_WORLD, nranks, stat)
      columns = column_range(n, rank, nranks)
      call new_integral(local_ints, n, columns)
      local_ints%overlap = ints%overlap(:, columns(1):columns(2))
      local_ints%hamiltonian = ints%hamiltonian(:, columns(1):columns(2))
      local_ints%dipole = ints%dipole(:, :, columns(1):columns(2))
      local_ints%quadrupole = ints%quadrupole(:, :, columns(1):columns(2))
      call mpi_allreduce_sum(error, local_ints, MPI_COMM_WORLD%MPI_VAL)
      call assert(allocated(error), "reject reduction of disjoint integral columns")
      deallocate(error)
      deallocate(h)
      allocate(h(n, columns(2)-columns(1)+1, nspin))
      call mpi_add_pot_to_h1(error, calc%bas, local_ints, saved_pot, h, MPI_COMM_WORLD%MPI_VAL)
      call check_error(error)
      call assert(all(abs(h-reference(:, columns(1):columns(2), :)) < thr), "local Hamiltonian expansion")
      call mpi_expand_integrals(error, local_ints, MPI_COMM_WORLD%MPI_VAL)
      call check_error(error)
      call assert(.not.local_ints%local, "expanded integral layout")
      call assert(all(local_ints%overlap == ints%overlap), "expanded overlap")
      call assert(all(local_ints%hamiltonian == ints%hamiltonian), "expanded Hamiltonian")
      call assert(all(local_ints%dipole == ints%dipole), "expanded dipoles")
      call assert(all(local_ints%quadrupole == ints%quadrupole), "expanded quadrupoles")
   end subroutine check_hamiltonian

   !> Occupied/partially occupied orbitals, both spins and signed energy weights.
   subroutine check_density(n, nspin)
      integer, intent(in) :: n, nspin
      type(sygvd_solver) :: serial
      type(psygvd_solver) :: distributed
      real(wp) :: hs(n, n, nspin), hd(n, n, nspin), smat(n, n), nel(2)
      real(wp) :: es(n, nspin), ed(n, nspin), fs(n, 2), fd(n, 2), saved(n, 2)
      real(wp) :: ps(n, n, nspin), pd(n, n, nspin)
      integer :: i, j, spin, rank, nranks, columns(2)
      real(wp), allocatable :: hc(:, :, :), pc(:, :, :)

      if (.not.tblite_has_scalapack) return
      smat = 0.0_wp
      do i = 1, n
         smat(i, i) = 1.0_wp + 0.01_wp*real(i, wp)
         do j = 1, n
            hs(j, i, 1) = 0.01_wp*cos(real(i*j, wp))
         end do
         hs(i, i, 1) = hs(i, i, 1) + real(i, wp) - 0.6_wp*n
      end do
      do spin = 2, nspin
         hs(:, :, spin) = 0.9_wp*hs(:, :, 1)
      end do
      hd = hs
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, stat)
      call MPI_Comm_size(MPI_COMM_WORLD, nranks, stat)
      columns = column_range(n, rank, nranks)
      hc = hs(:, columns(1):columns(2), :)
      nel = [0.45_wp*n, 0.35_wp*n]
      call new_sygvd(serial, smat, nel, 0.3_wp)
      call new_psygvd(distributed, smat, nel, 0.3_wp, MPI_COMM_WORLD%MPI_VAL)
      call serial%get_density(hs, smat, es, fs, ps, error)
      call check_error(error)
      call distributed%get_density(hd, smat, ed, fd, pd, error)
      call check_error(error)
      call assert(all(abs(ps-pd) < thr), "density matrix")
      call assert(all(abs(fs-fd) < thr), "orbital occupations")
      call distributed%delete()
      call new_psygvd(distributed, smat(:, columns(1):columns(2)), nel, 0.3_wp, MPI_COMM_WORLD%MPI_VAL)
      allocate(pc(n, size(hc, 2), nspin))
      call distributed%get_density(hc, smat(:, columns(1):columns(2)), ed, fd, pc, error)
      call check_error(error)
      call assert(all(abs(pc-ps(:, columns(1):columns(2), :)) < thr), "PBLAS column density")
      call assert(all(abs(fs-fd) < thr), "column orbital occupations")
      saved = fd
      call serial%get_wdensity(hs, smat, es, fs, ps, error)
      call check_error(error)
      call distributed%get_wdensity(hc, smat(:, columns(1):columns(2)), ed, fd, pc, error)
      call mpi_expand_matrix(error, pc, MPI_COMM_WORLD%MPI_VAL)
      pd = pc
      call check_error(error)
      call assert(all(abs(ps-pd) < thr), "energy weighted density matrix")
      call assert(all(abs(saved-fd) < thr), "restored orbital occupations")
      call distributed%delete()
   end subroutine check_density

   subroutine check_error(error)
      type(error_type), allocatable, intent(in) :: error
      if (allocated(error)) then
         write(*, "(2a)") "[Fatal] ", error%message
         call MPI_Abort(MPI_COMM_WORLD, 1, stat)
      end if
   end subroutine check_error

   subroutine assert(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not.condition) then
         write(*, "(3a)") "[Fatal] Distributed ", label, " does not match the serial result"
         call MPI_Abort(MPI_COMM_WORLD, 1, stat)
      end if
   end subroutine assert

end program test_mpi_singlepoint
