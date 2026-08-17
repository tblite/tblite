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
!> Without MPI support every entry point reports an error, the communicator
!> handles are meaningless in that case.
!>
!> Communicators are passed as plain integer handles to keep the MPI types out
!> of the rest of the library, ``MPI_Comm%MPI_VAL`` converts an mpi_f08 handle.
module tblite_mpi_utils
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_partition, only : work_partition, new_work_partition
#if TBLITE_HAS_MPI
   use mpi_f08, only : MPI_Comm, MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, &
      & MPI_LOGICAL, MPI_LOR, MPI_SUM, MPI_Allreduce, MPI_Comm_rank, MPI_Comm_size, &
      & MPI_Initialized
#endif
   implicit none
   private

   public :: get_mpi_comm_world, new_mpi_work_partition, mpi_allreduce_sum
   public :: mpi_sync_error

   !> Sum a partitioned result over all ranks of a communicator, in place
   interface mpi_allreduce_sum
      module procedure :: allreduce_sum_r0
      module procedure :: allreduce_sum_r1
      module procedure :: allreduce_sum_r2
      module procedure :: allreduce_sum_r3
   end interface mpi_allreduce_sum

   character(len=*), parameter :: no_mpi = "tblite was built without MPI support"


contains


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
   integer, intent(in) :: comm

   logical :: failed
   integer :: stat

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


subroutine allreduce_sum_r0(error, val, comm)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Partial result of this rank, replaced by the sum over all ranks
   real(wp), intent(inout) :: val
   !> Communicator to reduce over
   integer, intent(in) :: comm

   real(wp) :: buffer(1)

   buffer(1) = val
   call allreduce_sum_r1(error, buffer, comm)
   val = buffer(1)
end subroutine allreduce_sum_r0


subroutine allreduce_sum_r1(error, array, comm)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Partial result of this rank, replaced by the sum over all ranks
   real(wp), contiguous, intent(inout) :: array(:)
   !> Communicator to reduce over
   integer, intent(in) :: comm

   integer :: stat

#if TBLITE_HAS_MPI
   call MPI_Allreduce(MPI_IN_PLACE, array, size(array), MPI_DOUBLE_PRECISION, &
      & MPI_SUM, MPI_Comm(comm), stat)
   if (stat /= 0) call fatal_error(error, "Could not reduce partitioned results")
#else
   call fatal_error(error, no_mpi)
#endif
end subroutine allreduce_sum_r1


subroutine allreduce_sum_r2(error, array, comm)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Partial result of this rank, replaced by the sum over all ranks
   real(wp), contiguous, intent(inout), target :: array(:, :)
   !> Communicator to reduce over
   integer, intent(in) :: comm

   real(wp), pointer :: flat(:)

   flat(1:size(array)) => array
   call allreduce_sum_r1(error, flat, comm)
end subroutine allreduce_sum_r2


subroutine allreduce_sum_r3(error, array, comm)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Partial result of this rank, replaced by the sum over all ranks
   real(wp), contiguous, intent(inout), target :: array(:, :, :)
   !> Communicator to reduce over
   integer, intent(in) :: comm

   real(wp), pointer :: flat(:)

   flat(1:size(array)) => array
   call allreduce_sum_r1(error, flat, comm)
end subroutine allreduce_sum_r3


end module tblite_mpi_utils
