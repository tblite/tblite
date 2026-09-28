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
program test_mpi_singlepoint
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use mpi_f08, only : MPI_COMM_WORLD, MPI_Abort, MPI_Comm_rank, MPI_Comm_size, &
      & MPI_Finalize, MPI_Init
   use mstore, only : get_structure
   use tblite_ceh_ceh, only : new_ceh_calculator
   use tblite_ceh_singlepoint, only : ceh_singlepoint
   use tblite_container, only : container_type
   use tblite_context, only : context_type
   use tblite_features, only : tblite_has_scalapack
   use tblite_lapack_scalapack, only : psygvd_solver, new_psygvd, &
      & distribute_diagonalization
   use tblite_lapack_sygvd, only : sygvd_solver, new_sygvd
   use tblite_mpi_utils, only : mpi_allreduce_sum
   use tblite_partition, only : work_partition, serial_work_partition, operator(==)
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
   integer :: stat

   call MPI_Init(stat)

   call get_structure(mol, "MB16-43", "01")
   allocate(serial_gradient(3, mol%nat), mpi_gradient(3, mol%nat))

   call run(mol, .false., serial_energy, serial_gradient, serial_sigma, error)
   call check_error(error)
   call run(mol, .true., mpi_energy, mpi_gradient, mpi_sigma, error)
   call check_error(error)

   call assert(abs(mpi_energy - serial_energy) < thr, "energy")
   call assert(all(abs(mpi_gradient - serial_gradient) < thr), "gradient")
   call assert(all(abs(mpi_sigma - serial_sigma) < thr), "virial")

   call check_sync_error()
   call check_ceh()
   call check_eigensolver(37)
   call check_eigensolver(193)

   call MPI_Finalize(stat)

contains

   subroutine run(mol, distributed, energy, gradient, sigma, error)
      type(structure_type), intent(in) :: mol
      logical, intent(in) :: distributed
      real(wp), intent(out) :: energy, gradient(:, :), sigma(:, :)
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn

      ctx%verbosity = 0
      call new_gfn2_calculator(calc, mol, error)
      if (allocated(error)) return
      call add_solvation(calc, mol, error)
      if (allocated(error)) return

      if (distributed) then
         call ctx%set_mpi(error)
         if (allocated(error)) return
         call calc%set_partition(ctx%partition)
      end if

      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, 300.0_wp)
      call xtb_singlepoint(ctx, mol, calc, wfn, 1.0_wp, energy, gradient, sigma)
      if (ctx%failed()) call ctx%get_error(error)
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
      if (rank == nranks-1) partition = serial_work_partition
      call ctx%check_partition(partition, error)
      call assert(allocated(error), "partition mismatch synchronization")
      deallocate(error)

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
      do iteration = 1, 2
         ! Reuse the grid for a different Hamiltonian, as in successive SCF steps.
         amat(1, 1) = amat(1, 1) + 0.1_wp
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
