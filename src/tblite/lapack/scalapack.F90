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
!> Accepts either replicated matrices or complete local columns. The latter
!> are redistributed directly to the 2-D BLACS grid, including eigenvectors
!> and density matrices; no rank needs a full AO matrix during SCF.
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
   use tblite_scf_diag, only : diag_solver_type, delete_diag_solver
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

      subroutine pdgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ctxt)
         import :: dp
         integer, intent(in) :: m, n, ia, ja, desca(*), ib, jb, descb(*), ctxt
         real(dp), intent(in) :: a(*)
         real(dp), intent(out) :: b(*)
      end subroutine pdgemr2d

      subroutine pdgemm(transa, transb, m, n, k, alpha, a, ia, ja, desca, &
            & b, ib, jb, descb, beta, c, ic, jc, descc)
         import :: dp
         character(len=1), intent(in) :: transa, transb
         integer, intent(in) :: m, n, k, ia, ja, desca(*), ib, jb, descb(*), ic, jc, descc(*)
         real(dp), intent(in) :: alpha, beta, a(*), b(*)
         real(dp), intent(inout) :: c(*)
      end subroutine pdgemm

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
      integer :: column_ctxt = -1, column_desc(9) = 0
      logical :: local_columns = .false.
      !> Dimension of the eigenvalue problem
      integer :: n = 0
      !> Number of process columns, used for the workspace bound
      integer :: npcol = 1
      !> Global indices of the rows and columns owned by this rank
      integer, allocatable :: rows(:), cols(:)
      !> Descriptor shared by all distributed matrices
      integer :: desc(9) = 0
      !> Local matrices and eigensolver workspace reused across SCF iterations
      real(dp), allocatable :: aloc(:, :), bloc(:, :), zloc(:, :), work(:)
      integer, allocatable :: iwork(:)
      !> Unfactored overlap snapshot; invalidated on any rank's local change
      real(dp), allocatable :: overlap(:, :), ploc(:, :)
      logical :: overlap_valid = .false.
   contains
      procedure :: solve_sp
      procedure :: solve_dp
      procedure :: get_density_matrix => density_matrix
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
   if (allocated(error) .or. partition%get_nparts() <= 1) return

   call grid_shape(partition%get_nparts(), nprow, npcol)
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
   self%local_columns = size(overlap, 2) /= self%n
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
   integer :: iam, nprocs, nprow, myrow, mycol, nb, i

   ! BLACS derives the grid from the default system context, which is the global
   ! communicator, a partitioned subgroup would map to the wrong processes
   if (self%comm /= get_mpi_comm_world()) then
      call fatal_error(error, "ScaLAPACK solver requires the global communicator")
      return
   end if

   call blacs_pinfo(iam, nprocs)
   call grid_shape(nprocs, nprow, self%npcol)

   call blacs_get(0, 0, self%ctxt)
   call blacs_gridinit(self%ctxt, "R", nprow, self%npcol)
   call blacs_gridinfo(self%ctxt, nprow, self%npcol, myrow, mycol)

   if (self%local_columns) then
      call blacs_get(0, 0, self%column_ctxt)
      call blacs_gridinit(self%column_ctxt, "R", 1, nprocs)
      call descinit(self%column_desc, self%n, self%n, self%n, (self%n+nprocs-1)/nprocs, &
         & 0, 0, self%column_ctxt, max(1, self%n), stat)
      if (stat /= 0) then
         call handle_info(error, "column descinit", stat)
         return
      end if
   end if

   ! a block larger than the share of a process row leaves ranks without work
   nb = max(1, min(default_block_size, self%n/max(nprow, self%npcol)))
   ! Build the block-cyclic index maps once, then scatter/gather by indexing.
   self%rows = [(((i-1)/nb*nprow + myrow)*nb + modulo(i-1, nb) + 1, &
      & i=1, numroc(self%n, nb, myrow, 0, nprow))]
   self%cols = [(((i-1)/nb*self%npcol + mycol)*nb + modulo(i-1, nb) + 1, &
      & i=1, numroc(self%n, nb, mycol, 0, self%npcol))]

   call descinit(self%desc, self%n, self%n, nb, nb, 0, 0, self%ctxt, &
      & max(1, size(self%rows)), stat)
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

   real(dp) :: scale, wquery(1), changed(1)
   integer :: info, lwork, liwork, trilwmin, ormlwmin, iquery(1), mloc, nloc, nb

#if TBLITE_HAS_SCALAPACK
   if (self%ctxt < 0) then
      call setup_grid(self, error)
      if (allocated(error)) return
   end if

   mloc = size(self%rows)
   nloc = size(self%cols)
   nb = self%desc(5)
   if (.not.allocated(self%aloc)) then
      allocate(self%aloc(max(1, mloc), max(1, nloc)), &
         & self%bloc(max(1, mloc), max(1, nloc)), self%zloc(max(1, mloc), max(1, nloc)), &
         & self%ploc(max(1, mloc), max(1, nloc)))
      if (self%local_columns) then
         allocate(self%overlap(size(smat, 1), size(smat, 2)))
      else
         allocate(self%overlap(mloc, nloc))
      end if
   end if
   if (self%local_columns) then
      call pdgemr2d(self%n, self%n, hmat, 1, 1, self%column_desc, self%aloc, &
         & 1, 1, self%desc, self%ctxt)
   else
      self%aloc(:mloc, :nloc) = hmat(self%rows, self%cols)
   end if

   changed = 1.0_wp
   if (self%overlap_valid) then
      if (self%local_columns) then
         changed = merge(1.0_wp, 0.0_wp, any(self%overlap /= smat))
      else
         changed = merge(1.0_wp, 0.0_wp, any(self%overlap /= smat(self%rows, self%cols)))
      end if
   end if
   call mpi_allreduce_sum(error, changed, self%comm)
   if (allocated(error)) return
   if (changed(1) > 0.0_wp) then
      self%overlap_valid = .false.
      if (self%local_columns) then
         self%overlap = smat
         call pdgemr2d(self%n, self%n, smat, 1, 1, self%column_desc, self%bloc, &
            & 1, 1, self%desc, self%ctxt)
      else
         self%overlap = smat(self%rows, self%cols)
         self%bloc(:mloc, :nloc) = self%overlap
      end if
      call pdpotrf("u", self%n, self%bloc, 1, 1, self%desc, info)
      if (info /= 0) then
         call handle_info(error, "pdpotrf", info)
         return
      end if
      self%overlap_valid = .true.
   end if

   call pdsygst(1, "u", self%n, self%aloc, 1, 1, self%desc, self%bloc, 1, 1, self%desc, &
      & scale, info)
   if (info /= 0) then
      call handle_info(error, "pdsygst", info)
      return
   end if

   if (.not.allocated(self%work)) then
      call pdsyevd("v", "u", self%n, self%aloc, 1, 1, self%desc, eval, self%zloc, 1, 1, &
         & self%desc, wquery, -1, iquery, -1, info)
      if (info /= 0) then
         call handle_info(error, "pdsyevd", info)
         return
      end if
      ! the query under-reports and the documented pdsyevd minimum does not cover
      ! the pdormtr call it feeds its work array to, so both are used as a floor
      trilwmin = 3*self%n + max(nb*(mloc + 1), 3*nb)
      ormlwmin = (mloc + nloc + 2*nb)*nb + nb**2
      lwork = max(nint(wquery(1)), ormlwmin, &
         & max(1 + 6*self%n + 2*mloc*nloc, trilwmin) + 2*self%n)
      liwork = max(iquery(1), 7*self%n + 8*self%npcol + 2)
      allocate(self%work(lwork), self%iwork(liwork))
   end if

   call pdsyevd("v", "u", self%n, self%aloc, 1, 1, self%desc, eval, self%zloc, 1, 1, &
      & self%desc, self%work, size(self%work), self%iwork, size(self%iwork), info)
   if (info /= 0) then
      call handle_info(error, "pdsyevd", info)
      return
   end if

   ! back-transform the eigenvectors of the standard problem, x = U⁻¹ y,
   ! pdsygst documents scale as always 1 so the eigenvalues need no correction
   call pdtrsm("l", "u", "n", "n", self%n, self%n, 1.0_dp, self%bloc, 1, 1, &
      & self%desc, self%zloc, 1, 1, self%desc)

   if (self%local_columns) then
      call pdgemr2d(self%n, self%n, self%zloc, 1, 1, self%desc, hmat, &
         & 1, 1, self%column_desc, self%ctxt)
   else
      hmat(:, :) = 0.0_dp
      hmat(self%rows, self%cols) = self%zloc(:mloc, :nloc)
      call mpi_allreduce_sum(error, hmat, self%comm)
   end if
#else
   call fatal_error(error, no_scalapack)
#endif
end subroutine solve_dp


!> Build P (or energy-weighted P) entirely on the 2-D grid. The public
!> arrays are redistributed according to the caller's storage layout.
subroutine density_matrix(self, focc, coeff, pmat, error)
   class(psygvd_solver), intent(inout) :: self
   real(wp), intent(in) :: focc(:)
   real(wp), contiguous, intent(in) :: coeff(:, :)
   real(wp), contiguous, intent(inout) :: pmat(:, :)
   type(error_type), allocatable, intent(inout) :: error
   integer :: j

   if (allocated(error)) return
#if TBLITE_HAS_SCALAPACK
   if (self%local_columns) then
      call pdgemr2d(self%n, self%n, coeff, 1, 1, self%column_desc, self%zloc, &
         & 1, 1, self%desc, self%ctxt)
   else
      self%zloc(:size(self%rows), :size(self%cols)) = coeff(self%rows, self%cols)
   end if
   do j = 1, size(self%cols)
      self%aloc(:, j) = self%zloc(:, j)*focc(self%cols(j))
   end do
   call pdgemm("n", "t", self%n, self%n, self%n, 1.0_wp, self%aloc, 1, 1, self%desc, &
      & self%zloc, 1, 1, self%desc, 0.0_wp, self%ploc, 1, 1, self%desc)
   if (self%local_columns) then
      call pdgemr2d(self%n, self%n, self%ploc, 1, 1, self%desc, pmat, &
         & 1, 1, self%column_desc, self%ctxt)
   else
      pmat = 0.0_wp
      pmat(self%rows, self%cols) = self%ploc(:size(self%rows), :size(self%cols))
      call mpi_allreduce_sum(error, pmat, self%comm)
   end if
#else
   call fatal_error(error, no_scalapack)
#endif
end subroutine density_matrix


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
   if (self%column_ctxt >= 0) call blacs_gridexit(self%column_ctxt)
   if (self%ctxt >= 0) call blacs_gridexit(self%ctxt)
#endif
   self%ctxt = -1
   self%column_ctxt = -1
   if (allocated(self%rows)) deallocate(self%rows, self%cols)
   if (allocated(self%aloc)) deallocate(self%aloc, self%bloc, self%zloc, self%overlap)
   if (allocated(self%ploc)) deallocate(self%ploc)
   if (allocated(self%work)) deallocate(self%work, self%iwork)
   self%overlap_valid = .false.
   call delete_diag_solver(self)
end subroutine delete


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
