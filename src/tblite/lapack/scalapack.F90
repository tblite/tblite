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

#ifndef TBLITE_HAS_SCALAPACK
#define TBLITE_HAS_SCALAPACK 0
#endif

!> @file tblite/lapack/scalapack.F90
!> Provides a distributed solver for the general symmetric eigenvalue problem

!> Distributed divide-and-conquer solver based on ScaLAPACK.
!>
!> The Hamiltonian and overlap matrices are replicated on every rank, so the
!> block-cyclic distribution is a local copy without communication and only the
!> eigenvectors have to be collected again. This distributes the cubically
!> scaling diagonalization, the memory still holds the full matrices everywhere.
!>
!> Without ScaLAPACK support the solver reports an error instead of falling back
!> to a replicated diagonalization.
module tblite_lapack_scalapack
   use mctc_env, only : sp, dp, wp, error_type, fatal_error
   use tblite_features, only : tblite_has_scalapack
   use tblite_mpi_utils, only : get_mpi_comm_world, mpi_allreduce_sum, &
      & new_mpi_work_partition
   use tblite_output_format, only : format_string
   use tblite_partition, only : work_partition
   use tblite_scf_diag, only : diag_solver_type
   implicit none
   private

   public :: psygvd_solver, new_psygvd, distribute_diagonalization

#if TBLITE_HAS_SCALAPACK
   interface
      !> Number of rows or columns of a block-cyclic matrix owned by a process
      pure integer function numroc(n, nb, iproc, isrcproc, nprocs)
         integer, intent(in) :: n, nb, iproc, isrcproc, nprocs
      end function numroc

      !> Initialize the descriptor of a block-cyclic matrix
      pure subroutine descinit(desc, m, n, mb, nb, irsrc, icsrc, ctxt, lld, info)
         integer, intent(out) :: desc(*)
         integer, intent(in) :: m, n, mb, nb, irsrc, icsrc, ctxt, lld
         integer, intent(out) :: info
      end subroutine descinit

      !> Cholesky factorization of a distributed symmetric positive definite matrix
      subroutine pdpotrf(uplo, n, a, ia, ja, desca, info)
         import :: dp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n, ia, ja, desca(*)
         real(dp), intent(inout) :: a(*)
         integer, intent(out) :: info
      end subroutine pdpotrf

      !> Reduce a distributed general eigenvalue problem to standard form
      subroutine pdsygst(ibtype, uplo, n, a, ia, ja, desca, b, ib, jb, descb, &
            & scale, info)
         import :: dp
         integer, intent(in) :: ibtype
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n, ia, ja, desca(*), ib, jb, descb(*)
         real(dp), intent(inout) :: a(*)
         real(dp), intent(in) :: b(*)
         real(dp), intent(out) :: scale
         integer, intent(out) :: info
      end subroutine pdsygst

      !> Divide-and-conquer solver for a distributed standard eigenvalue problem
      subroutine pdsyevd(jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, &
            & work, lwork, iwork, liwork, info)
         import :: dp
         character(len=1), intent(in) :: jobz, uplo
         integer, intent(in) :: n, ia, ja, desca(*), iz, jz, descz(*)
         real(dp), intent(inout) :: a(*)
         real(dp), intent(out) :: w(*), z(*)
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork, liwork
         integer, intent(inout) :: iwork(*)
         integer, intent(out) :: info
      end subroutine pdsyevd

      !> Solve a distributed triangular system with multiple right hand sides
      subroutine pdtrsm(side, uplo, transa, diag, m, n, alpha, a, ia, ja, desca, &
            & b, ib, jb, descb)
         import :: dp
         character(len=1), intent(in) :: side, uplo, transa, diag
         integer, intent(in) :: m, n, ia, ja, desca(*), ib, jb, descb(*)
         real(dp), intent(in) :: alpha, a(*)
         real(dp), intent(inout) :: b(*)
      end subroutine pdtrsm

      subroutine blacs_get(ctxt, what, val)
         integer, intent(in) :: ctxt, what
         integer, intent(out) :: val
      end subroutine blacs_get

      subroutine blacs_pinfo(mypnum, nprocs)
         integer, intent(out) :: mypnum, nprocs
      end subroutine blacs_pinfo

      subroutine blacs_gridinit(ctxt, order, nprow, npcol)
         integer, intent(inout) :: ctxt
         character(len=1), intent(in) :: order
         integer, intent(in) :: nprow, npcol
      end subroutine blacs_gridinit

      subroutine blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
         integer, intent(in) :: ctxt
         integer, intent(out) :: nprow, npcol, myrow, mycol
      end subroutine blacs_gridinfo

      subroutine blacs_gridexit(ctxt)
         integer, intent(in) :: ctxt
      end subroutine blacs_gridexit
   end interface
#endif

   !> Default block size of the block-cyclic distribution
   integer, parameter :: default_block_size = 64

   character(len=*), parameter :: no_scalapack = &
      & "tblite was built without ScaLAPACK support"


   !> Distributed solver for the general symmetric eigenvalue problem
   type, extends(diag_solver_type) :: psygvd_solver
      !> Communicator the diagonalization is distributed over
      integer :: comm = 0
      !> BLACS context of the process grid, negative until the grid is created
      integer :: ctxt = -1
      !> Dimension of the eigenvalue problem
      integer :: n = 0
      !> Block size of the block-cyclic distribution
      integer :: nb = default_block_size
      !> Shape of the process grid
      integer :: nprow = 1, npcol = 1
      !> Position of this rank in the process grid
      integer :: myrow = 0, mycol = 0
      !> Local shape of the distributed matrices
      integer :: mloc = 0, nloc = 0
      !> Descriptor shared by all distributed matrices
      integer :: desc(9) = 0
   contains
      procedure :: solve_sp
      procedure :: solve_dp
      procedure :: delete
   end type psygvd_solver


contains


!> Whether distributing the diagonalization over a communicator pays off, a
!> single rank or a matrix too small to fill one block per process row is
!> faster with a replicated diagonalization
function distribute_diagonalization(n, comm) result(distribute)

   !> Dimension of the eigenvalue problem
   integer, intent(in) :: n

   !> Communicator the calculation is distributed over
   integer, intent(in) :: comm

   !> Whether to use the distributed solver
   logical :: distribute

   type(error_type), allocatable :: error
   type(work_partition) :: partition
   integer :: nprow, npcol

   distribute = .false.
   if (.not.tblite_has_scalapack) return

   call new_mpi_work_partition(error, partition, comm)
   if (allocated(error) .or. partition%nparts <= 1) return

   call grid_shape(partition%nparts, nprow, npcol)
   distribute = n >= default_block_size*max(nprow, npcol)

end function distribute_diagonalization


!> Create a distributed eigenvalue solver, the process grid is set up on the
!> first solve where errors can be reported
subroutine new_psygvd(self, overlap, nel, kt, comm)
   !> Instance of the distributed solver
   type(psygvd_solver), intent(out) :: self
   !> Overlap matrix
   real(wp), intent(in) :: overlap(:, :)
   !> Number of electrons per spin channel
   real(wp), intent(in) :: nel(:)
   !> Electronic temperature
   real(wp), intent(in) :: kt
   !> Communicator to distribute the diagonalization over
   integer, intent(in) :: comm

   self%n = size(overlap, 1)
   self%nel = nel
   self%kt = kt
   self%comm = comm
end subroutine new_psygvd


!> Create the process grid and the matrix descriptor
subroutine setup_grid(self, error)
   !> Instance of the distributed solver
   class(psygvd_solver), intent(inout) :: self
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat

#if TBLITE_HAS_SCALAPACK
   integer :: iam, nprocs

   ! BLACS derives the grid from the default system context, which is the global
   ! communicator, a partitioned subgroup would map to the wrong processes
   if (self%comm /= get_mpi_comm_world()) then
      call fatal_error(error, "ScaLAPACK solver requires the global communicator")
      return
   end if

   call blacs_pinfo(iam, nprocs)
   call grid_shape(nprocs, self%nprow, self%npcol)

   call blacs_get(0, 0, self%ctxt)
   call blacs_gridinit(self%ctxt, "R", self%nprow, self%npcol)
   call blacs_gridinfo(self%ctxt, self%nprow, self%npcol, self%myrow, self%mycol)

   ! a block larger than the share of a process row leaves ranks without work
   self%nb = max(1, min(default_block_size, self%n/max(self%nprow, self%npcol)))

   self%mloc = numroc(self%n, self%nb, self%myrow, 0, self%nprow)
   self%nloc = numroc(self%n, self%nb, self%mycol, 0, self%npcol)

   call descinit(self%desc, self%n, self%n, self%nb, self%nb, 0, 0, self%ctxt, &
      & max(1, self%mloc), stat)
   if (stat /= 0) then
      call fatal_error(error, "(descinit) failed to describe the distributed matrix.&
         & info="//format_string(stat, "(i0)"))
   end if
#else
   call fatal_error(error, no_scalapack)
#endif
end subroutine setup_grid


!> Solve the general symmetric eigenvalue problem in double precision
subroutine solve_dp(self, hmat, smat, eval, error)
   !> Instance of the distributed solver
   class(psygvd_solver), intent(inout) :: self
   !> Hamiltonian matrix, contains the eigenvectors on output
   real(dp), contiguous, intent(inout) :: hmat(:, :)
   !> Overlap matrix
   real(dp), contiguous, intent(in) :: smat(:, :)
   !> Eigenvalues
   real(dp), contiguous, intent(inout) :: eval(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(dp), allocatable :: aloc(:, :), bloc(:, :), zloc(:, :), work(:)
   integer, allocatable :: iwork(:)
   real(dp) :: scale, wquery(1)
   integer :: info, lwork, liwork, trilwmin, ormlwmin, iquery(1)

#if TBLITE_HAS_SCALAPACK
   if (self%ctxt < 0) then
      call setup_grid(self, error)
      if (allocated(error)) return
   end if

   allocate(aloc(max(1, self%mloc), max(1, self%nloc)))
   allocate(bloc(max(1, self%mloc), max(1, self%nloc)))
   allocate(zloc(max(1, self%mloc), max(1, self%nloc)))
   call scatter(self, hmat, aloc)
   call scatter(self, smat, bloc)

   call pdpotrf("u", self%n, bloc, 1, 1, self%desc, info)
   if (info /= 0) then
      call handle_info(error, "pdpotrf", info)
      return
   end if

   call pdsygst(1, "u", self%n, aloc, 1, 1, self%desc, bloc, 1, 1, self%desc, &
      & scale, info)
   if (info /= 0) then
      call handle_info(error, "pdsygst", info)
      return
   end if

   call pdsyevd("v", "u", self%n, aloc, 1, 1, self%desc, eval, zloc, 1, 1, &
      & self%desc, wquery, -1, iquery, -1, info)
   if (info /= 0) then
      call handle_info(error, "pdsyevd", info)
      return
   end if
   ! the query under-reports and the documented pdsyevd minimum does not cover
   ! the pdormtr call it feeds its work array to, so both are used as a floor
   trilwmin = 3*self%n + max(self%nb*(self%mloc + 1), 3*self%nb)
   ormlwmin = (self%mloc + self%nloc + 2*self%nb)*self%nb + self%nb**2
   lwork = max(nint(wquery(1)), ormlwmin, &
      & max(1 + 6*self%n + 2*self%mloc*self%nloc, trilwmin) + 2*self%n)
   liwork = max(iquery(1), 7*self%n + 8*self%npcol + 2)
   allocate(work(lwork), iwork(liwork))

   call pdsyevd("v", "u", self%n, aloc, 1, 1, self%desc, eval, zloc, 1, 1, &
      & self%desc, work, lwork, iwork, liwork, info)
   if (info /= 0) then
      call handle_info(error, "pdsyevd", info)
      return
   end if

   ! back-transform the eigenvectors of the standard problem, x = U⁻¹ y,
   ! pdsygst documents scale as always 1 so the eigenvalues need no correction
   call pdtrsm("l", "u", "n", "n", self%n, self%n, 1.0_dp, bloc, 1, 1, &
      & self%desc, zloc, 1, 1, self%desc)

   call gather(self, zloc, hmat)
   call mpi_allreduce_sum(error, hmat, self%comm)
#else
   call fatal_error(error, no_scalapack)
#endif
end subroutine solve_dp


!> Single precision is not provided, the solver is only used for the
!> double precision self-consistent iterations
subroutine solve_sp(self, hmat, smat, eval, error)
   !> Instance of the distributed solver
   class(psygvd_solver), intent(inout) :: self
   !> Hamiltonian matrix, contains the eigenvectors on output
   real(sp), contiguous, intent(inout) :: hmat(:, :)
   !> Overlap matrix
   real(sp), contiguous, intent(in) :: smat(:, :)
   !> Eigenvalues
   real(sp), contiguous, intent(inout) :: eval(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call fatal_error(error, "ScaLAPACK solver does not support single precision")
end subroutine solve_sp


!> Release the process grid
subroutine delete(self)
   !> Instance of the distributed solver
   class(psygvd_solver), intent(inout) :: self

#if TBLITE_HAS_SCALAPACK
   if (self%ctxt >= 0) call blacs_gridexit(self%ctxt)
#endif
   self%ctxt = -1
end subroutine delete


!> Copy the share of this rank out of the replicated matrix, the whole matrix is
!> available everywhere so this needs no communication
pure subroutine scatter(self, glob, loc)
   !> Instance of the distributed solver
   class(psygvd_solver), intent(in) :: self
   !> Replicated matrix
   real(dp), intent(in) :: glob(:, :)
   !> Block-cyclic share of this rank
   real(dp), intent(out) :: loc(:, :)

   integer :: iloc, jloc

   do jloc = 1, self%nloc
      do iloc = 1, self%mloc
         loc(iloc, jloc) = glob(local_to_global(self, iloc, self%myrow, self%nprow), &
            & local_to_global(self, jloc, self%mycol, self%npcol))
      end do
   end do
end subroutine scatter


!> Place the share of this rank into the replicated matrix, the remaining
!> entries are left zero for the reduction to fill in
pure subroutine gather(self, loc, glob)
   !> Instance of the distributed solver
   class(psygvd_solver), intent(in) :: self
   !> Block-cyclic share of this rank
   real(dp), intent(in) :: loc(:, :)
   !> Replicated matrix
   real(dp), intent(out) :: glob(:, :)

   integer :: iloc, jloc

   glob(:, :) = 0.0_dp
   do jloc = 1, self%nloc
      do iloc = 1, self%mloc
         glob(local_to_global(self, iloc, self%myrow, self%nprow), &
            & local_to_global(self, jloc, self%mycol, self%npcol)) = loc(iloc, jloc)
      end do
   end do
end subroutine gather


!> Global index of a local index of the block-cyclic distribution
pure function local_to_global(self, idx, iproc, nprocs) result(glob)
   !> Instance of the distributed solver
   class(psygvd_solver), intent(in) :: self
   !> Local index
   integer, intent(in) :: idx
   !> Position of this rank along the distributed dimension
   integer, intent(in) :: iproc
   !> Number of ranks along the distributed dimension
   integer, intent(in) :: nprocs
   !> Global index
   integer :: glob

   glob = (((idx - 1)/self%nb)*nprocs + iproc)*self%nb + modulo(idx - 1, self%nb) + 1
end function local_to_global


!> Squarest process grid for a number of ranks
pure subroutine grid_shape(nranks, nprow, npcol)
   !> Number of ranks to arrange
   integer, intent(in) :: nranks
   !> Shape of the process grid
   integer, intent(out) :: nprow, npcol

   nprow = int(sqrt(real(nranks)))
   do while (nprow > 1 .and. modulo(nranks, nprow) /= 0)
      nprow = nprow - 1
   end do
   npcol = nranks/max(1, nprow)
end subroutine grid_shape


subroutine handle_info(error, name, info)
   type(error_type), allocatable, intent(out) :: error
   character(len=*), intent(in) :: name
   integer, intent(in) :: info

   call fatal_error(error, "("//name//") failed to solve eigenvalue problem.&
      & info="//format_string(info, "(i0)"))
end subroutine handle_info


end module tblite_lapack_scalapack
