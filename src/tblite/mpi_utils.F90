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

#ifndef TBLITE_HAS_MPI
#define TBLITE_HAS_MPI 0
#endif

!> @file tblite/mpi_utils.F90
!> Provides the MPI primitives used to distribute a partitioned calculation

!> Thin wrappers around the MPI calls needed for a partitioned calculation.
!>
!> An absent communicator leaves serial results untouched. An explicit
!> communicator requires MPI support.
!>
!> Communicators are passed as plain integer handles to keep the MPI types out
!> of the rest of the library, ``MPI_Comm%MPI_VAL`` converts an mpi_f08 handle.
module tblite_mpi_utils
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_blas, only : gemm
   use tblite_integral_type, only : integral_type
   use tblite_partition, only : work_partition, new_work_partition, column_range
#if TBLITE_HAS_MPI
   use mpi_f08, only : MPI_Comm, MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, &
      & MPI_LOGICAL, MPI_LOR, MPI_SUM, MPI_Allreduce, MPI_Comm_rank, MPI_Comm_size, &
      & MPI_Finalize, MPI_Init, MPI_Initialized, MPI_Allgatherv, MPI_Alltoallv, MPI_Bcast
#endif
   implicit none
   private

   public :: get_mpi_comm_world, new_mpi_work_partition, mpi_allreduce_sum
   public :: mpi_sync_error, mpi_startup, mpi_shutdown
   public :: mpi_transpose_columns, mpi_gather_columns, mpi_expand_matrix
   public :: mpi_multiply_columns, mpi_density_columns, mpi_expand_integrals

   !> Bound MPI-library temporary storage to 16 MiB payloads for real reductions.
   integer, parameter :: reduce_chunk_size = 2**21

   !> Sum a partitioned result over all ranks of a communicator, in place.
   !> An absent communicator or pending error skips the reduction.
   interface mpi_allreduce_sum
      module procedure :: allreduce_sum_r1
      module procedure :: allreduce_sum_r2
      module procedure :: allreduce_sum_r3
      module procedure :: allreduce_sum_integrals
   end interface mpi_allreduce_sum

   interface mpi_expand_matrix
      module procedure :: expand_matrix_r2, expand_matrix_r3
   end interface

   character(len=*), parameter :: no_mpi = "tblite was built without MPI support"


contains


!> Enter the MPI environment, a no-op without MPI support
subroutine mpi_startup()

#if TBLITE_HAS_MPI
   integer :: stat

   call MPI_Init(stat)
#endif

end subroutine mpi_startup


!> Leave the MPI environment, a no-op without MPI support
subroutine mpi_shutdown()

#if TBLITE_HAS_MPI
   integer :: stat

   call MPI_Finalize(stat)
#endif

end subroutine mpi_shutdown


!> Handle of the global communicator, meaningless without MPI support
function get_mpi_comm_world() result(comm)

   !> Handle of the global communicator
   integer :: comm

#if TBLITE_HAS_MPI
   comm = MPI_COMM_WORLD%MPI_VAL
#else
   comm = 0
#endif

end function get_mpi_comm_world


!> Derive the work partition of this rank from the size of a communicator
subroutine new_mpi_work_partition(error, partition, comm)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Share of the interaction loops evaluated by this rank
   type(work_partition), intent(out) :: partition

   !> Communicator to distribute over
   integer, intent(in) :: comm

   integer :: rank, nranks, stat

#if TBLITE_HAS_MPI
   logical :: ready

   ! calling into MPI before MPI_Init aborts the whole program
   call MPI_Initialized(ready, stat)
   if (stat /= 0 .or. .not.ready) then
      call fatal_error(error, "MPI is not initialized")
      return
   end if

   call MPI_Comm_rank(MPI_Comm(comm), rank, stat)
   if (stat /= 0) then
      call fatal_error(error, "Could not determine the rank of this process")
      return
   end if
   call MPI_Comm_size(MPI_Comm(comm), nranks, stat)
   if (stat /= 0) then
      call fatal_error(error, "Could not determine the size of the communicator")
      return
   end if
   call new_work_partition(error, partition, rank, nranks)
#else
   call fatal_error(error, no_mpi)
#endif

end subroutine new_mpi_work_partition


!> Make a failure on any rank of the communicator visible to all of them, a rank
!> leaving a collective on its own would deadlock the remaining ones
subroutine mpi_sync_error(error, comm)

   !> Error handling, allocated on every rank if any rank failed
   type(error_type), allocatable, intent(inout) :: error

   !> Communicator to synchronize over
   integer, intent(in), optional :: comm

   logical :: failed
   integer :: stat

   if (.not.present(comm)) return

#if TBLITE_HAS_MPI
   failed = allocated(error)
   call MPI_Allreduce(MPI_IN_PLACE, failed, 1, MPI_LOGICAL, MPI_LOR, MPI_Comm(comm), stat)
   if (stat /= 0) failed = .true.

   if (failed .and. .not.allocated(error)) then
      call fatal_error(error, "Calculation failed on another rank")
   end if
#else
   if (.not.allocated(error)) call fatal_error(error, no_mpi)
#endif

end subroutine mpi_sync_error


subroutine allreduce_sum_r1(error, array, comm)
   !> Error handling
   type(error_type), allocatable, intent(inout) :: error
   !> Partial result of this rank, replaced by the sum over all ranks
   real(wp), contiguous, intent(inout) :: array(:)
   !> Communicator to reduce over
   integer, intent(in), optional :: comm

   integer :: stat, first, count

   if (.not.present(comm) .or. allocated(error)) return

#if TBLITE_HAS_MPI
   ! An in-place collective can still allocate a full-size temporary inside MPI.
   ! Stream large arrays so these buffers do not scale with the full AO tensor.
   first = 1
   do
      count = min(reduce_chunk_size, size(array)-first+1)
      call MPI_Allreduce(MPI_IN_PLACE, array(first:first+count-1), count, MPI_DOUBLE_PRECISION, &
         & MPI_SUM, MPI_Comm(comm), stat)
      if (stat /= 0) then
         call fatal_error(error, "Could not reduce partitioned results")
         return
      end if
      if (count == size(array)-first+1) exit
      first = first + count
   end do
#else
   call fatal_error(error, no_mpi)
#endif
end subroutine allreduce_sum_r1


subroutine allreduce_sum_r2(error, array, comm)
   !> Error handling
   type(error_type), allocatable, intent(inout) :: error
   !> Partial result of this rank, replaced by the sum over all ranks
   real(wp), contiguous, intent(inout), target :: array(:, :)
   !> Communicator to reduce over
   integer, intent(in), optional :: comm

   ! Preserve contiguity through rank remapping to avoid a full copy-in/out.
   real(wp), contiguous, pointer :: flat(:)

   flat(1:size(array)) => array
   call allreduce_sum_r1(error, flat, comm)
end subroutine allreduce_sum_r2


subroutine allreduce_sum_r3(error, array, comm)
   !> Error handling
   type(error_type), allocatable, intent(inout) :: error
   !> Partial result of this rank, replaced by the sum over all ranks
   real(wp), contiguous, intent(inout), target :: array(:, :, :)
   !> Communicator to reduce over
   integer, intent(in), optional :: comm

   real(wp), contiguous, pointer :: flat(:)

   flat(1:size(array)) => array
   call allreduce_sum_r1(error, flat, comm)
end subroutine allreduce_sum_r3


!> Assemble the one-electron integrals before constructing the Hamiltonian.
subroutine allreduce_sum_integrals(error, ints, comm)
   type(error_type), allocatable, intent(inout) :: error
   type(integral_type), intent(inout) :: ints
   integer, intent(in), optional :: comm

   if (.not.present(comm)) return
   if (allocated(error)) return
   if (ints%local) then
      call fatal_error(error, "Complete local integral columns must not be sum-reduced")
      return
   end if
   call mpi_allreduce_sum(error, ints%overlap, comm)
   call mpi_allreduce_sum(error, ints%hamiltonian, comm)
   call mpi_allreduce_sum(error, ints%dipole, comm)
   call mpi_allreduce_sum(error, ints%quadrupole, comm)
end subroutine allreduce_sum_integrals


!> Transpose a square matrix stored as complete local columns. Only the local
!> block and two equally sized communication buffers are needed.
subroutine mpi_transpose_columns(error, a, b, comm)
   type(error_type), allocatable, intent(inout) :: error
   real(wp), contiguous, intent(in) :: a(:, :)
   real(wp), contiguous, intent(out) :: b(:, :)
   integer, intent(in), optional :: comm
   type(work_partition) :: partition
   integer :: n, nlocal, part, bounds(2), i, j, k, stat
   integer, allocatable :: counts(:), displs(:)
   real(wp), allocatable :: sendbuf(:), recvbuf(:)

   if (allocated(error)) return
   if (present(comm)) call new_mpi_work_partition(error, partition, comm)
   if (allocated(error)) return
   if (partition%get_nparts() == 1) then
      b = transpose(a)
      return
   end if
#if TBLITE_HAS_MPI
   n = size(a, 1)
   nlocal = size(a, 2)
   allocate(counts(partition%get_nparts()), displs(partition%get_nparts()))
   allocate(sendbuf(size(a)), recvbuf(size(a)))
   k = 0
   do part = 0, partition%get_nparts()-1
      bounds = column_range(n, part, partition%get_nparts())
      displs(part+1) = k
      counts(part+1) = nlocal*(bounds(2)-bounds(1)+1)
      do j = 1, nlocal
         do i = bounds(1), bounds(2)
            k = k + 1
            sendbuf(k) = a(i, j)
         end do
      end do
   end do
   call MPI_Alltoallv(sendbuf, counts, displs, MPI_DOUBLE_PRECISION, &
      & recvbuf, counts, displs, MPI_DOUBLE_PRECISION, MPI_Comm(comm), stat)
   if (stat /= 0) then
      call fatal_error(error, "Could not transpose distributed matrix columns")
      return
   end if
   k = 0
   do part = 0, partition%get_nparts()-1
      bounds = column_range(n, part, partition%get_nparts())
      do j = bounds(1), bounds(2)
         do i = 1, nlocal
            k = k + 1
            b(j, i) = recvbuf(k)
         end do
      end do
   end do
#endif
end subroutine mpi_transpose_columns

!> Collect columns into a complete matrix, also for flattened tensors. When a
!> is absent, b already holds this rank's columns and the gather is in place.
subroutine mpi_gather_columns(error, a, b, n, comm)
   type(error_type), allocatable, intent(inout) :: error
   real(wp), contiguous, intent(in), optional :: a(:, :)
   real(wp), contiguous, intent(inout) :: b(:, :)
   integer, intent(in) :: n
   integer, intent(in), optional :: comm
   type(work_partition) :: partition
   integer :: part, bounds(2), stat
   integer, allocatable :: counts(:), displs(:)

   if (allocated(error)) return
   if (present(comm)) call new_mpi_work_partition(error, partition, comm)
   if (allocated(error)) return
   if (partition%get_nparts() == 1) then
      if (present(a)) b = a
      return
   end if
#if TBLITE_HAS_MPI
   allocate(counts(partition%get_nparts()), displs(partition%get_nparts()))
   do part = 0, partition%get_nparts()-1
      bounds = column_range(n, part, partition%get_nparts())
      counts(part+1) = size(b, 1)*(bounds(2)-bounds(1)+1)
      displs(part+1) = size(b, 1)*(bounds(1)-1)
   end do
   if (present(a)) then
      call MPI_Allgatherv(a, size(a), MPI_DOUBLE_PRECISION, b, counts, displs, &
         & MPI_DOUBLE_PRECISION, MPI_Comm(comm), stat)
   else
      call MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE_PRECISION, b, counts, displs, &
         & MPI_DOUBLE_PRECISION, MPI_Comm(comm), stat)
   end if
   if (stat /= 0) call fatal_error(error, "Could not collect distributed matrix columns")
#endif
end subroutine mpi_gather_columns

subroutine expand_matrix_r2(error, a, comm)
   type(error_type), allocatable, intent(inout) :: error
   real(wp), allocatable, intent(inout) :: a(:, :)
   integer, intent(in), optional :: comm
   real(wp), allocatable :: full(:, :)
   if (allocated(error)) return
   allocate(full(size(a, 1), size(a, 1)))
   call mpi_gather_columns(error, a, full, size(a, 1), comm)
   if (.not.allocated(error)) call move_alloc(full, a)
end subroutine expand_matrix_r2

subroutine expand_matrix_r3(error, a, comm)
   type(error_type), allocatable, intent(inout) :: error
   real(wp), allocatable, intent(inout) :: a(:, :, :)
   integer, intent(in), optional :: comm
   real(wp), allocatable :: full(:, :, :)
   integer :: spin
   if (allocated(error)) return
   allocate(full(size(a, 1), size(a, 1), size(a, 3)))
   do spin = 1, size(a, 3)
      call mpi_gather_columns(error, a(:, :, spin), full(:, :, spin), size(a, 1), comm)
   end do
   if (.not.allocated(error)) call move_alloc(full, a)
end subroutine expand_matrix_r3

!> C = A B, with all three matrices column distributed. Broadcast one local
!> A panel at a time, retaining O(N^2 / ranks) storage throughout.
subroutine mpi_multiply_columns(error, a, b, c, comm)
   type(error_type), allocatable, intent(inout) :: error
   real(wp), contiguous, intent(in) :: a(:, :), b(:, :)
   real(wp), contiguous, intent(out) :: c(:, :)
   integer, intent(in), optional :: comm
   type(work_partition) :: partition
   real(wp), allocatable :: panel(:, :)
   integer :: n, part, bounds(2), nc, stat

   if (allocated(error)) return
   if (present(comm)) call new_mpi_work_partition(error, partition, comm)
   if (allocated(error)) return
   if (partition%get_nparts() == 1) then
      call gemm(a, b, c)
      return
   end if
#if TBLITE_HAS_MPI
   n = size(a, 1)
   bounds = column_range(n, 0, partition%get_nparts())
   allocate(panel(n, bounds(2)))
   c = 0.0_wp
   do part = 0, partition%get_nparts()-1
      bounds = column_range(n, part, partition%get_nparts())
      nc = bounds(2)-bounds(1)+1
      if (nc == 0) cycle
      if (part == partition%get_part()) panel(:, :nc) = a
      call MPI_Bcast(panel, n*nc, MPI_DOUBLE_PRECISION, part, MPI_Comm(comm), stat)
      if (stat /= 0) then
         call fatal_error(error, "Could not broadcast a distributed matrix panel")
         return
      end if
      if (size(c, 2) > 0) call gemm(panel(:, :nc), b(bounds(1):bounds(2), :), c, beta=1.0_wp)
   end do
#endif
end subroutine mpi_multiply_columns

!> Density from MO-column-distributed coefficients. Used by spin-resolved
!> post-processing; ScaLAPACK uses PBLAS for its repeated SCF density builds.
subroutine mpi_density_columns(error, focc, coeff, pmat, comm)
   type(error_type), allocatable, intent(inout) :: error
   real(wp), intent(in) :: focc(:)
   real(wp), contiguous, intent(in) :: coeff(:, :)
   real(wp), contiguous, intent(out) :: pmat(:, :)
   integer, intent(in), optional :: comm
   type(work_partition) :: partition
   real(wp), allocatable :: weighted(:, :), transposed(:, :)
   integer :: bounds(2), j

   if (allocated(error)) return
   if (present(comm)) call new_mpi_work_partition(error, partition, comm)
   if (allocated(error)) return
   bounds = partition%get_columns(size(coeff, 1))
   allocate(weighted(size(coeff, 1), size(coeff, 2)), transposed(size(coeff, 1), size(coeff, 2)))
   do j = 1, size(coeff, 2)
      weighted(:, j) = coeff(:, j)*focc(bounds(1)+j-1)
   end do
   call mpi_transpose_columns(error, weighted, transposed, comm)
   deallocate(weighted)
   call mpi_multiply_columns(error, coeff, transposed, pmat, comm)
end subroutine mpi_density_columns

!> Compatibility path for consumers explicitly requesting full AO integrals.
subroutine mpi_expand_integrals(error, ints, comm)
   type(error_type), allocatable, intent(inout) :: error
   type(integral_type), intent(inout) :: ints
   integer, intent(in), optional :: comm
   if (.not.ints%local .or. allocated(error)) return
   call mpi_expand_matrix(error, ints%overlap, comm)
   call mpi_expand_matrix(error, ints%hamiltonian, comm)
   call expand_multipoles(error, ints%dipole, comm)
   call expand_multipoles(error, ints%quadrupole, comm)
   if (allocated(error)) return
   ints%local = .false.
   ints%columns = [1, size(ints%overlap, 1)]
end subroutine mpi_expand_integrals

subroutine expand_multipoles(error, a, comm)
   type(error_type), allocatable, intent(inout) :: error
   real(wp), allocatable, target, intent(inout) :: a(:, :, :)
   integer, intent(in), optional :: comm
   real(wp), allocatable, target :: full(:, :, :)
   real(wp), contiguous, pointer :: source(:, :), dest(:, :)
   integer :: n, nrow
   if (allocated(error)) return
   n = size(a, 2)
   nrow = size(a, 1)*n
   allocate(full(size(a, 1), n, n))
   source(1:nrow, 1:size(a, 3)) => a
   dest(1:nrow, 1:n) => full
   call mpi_gather_columns(error, source, dest, n, comm)
   if (.not.allocated(error)) call move_alloc(full, a)
end subroutine expand_multipoles

end module tblite_mpi_utils
