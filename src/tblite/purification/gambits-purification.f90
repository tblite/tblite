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

!> @file tblite/purification/gambits-purification.f90
!> Wrapper for the Gambits purification library
module tblite_purification_solver
    use mctc_env, only : sp, dp, error_type, wp, fatal_error
    use tblite_scf_solver, only : solver_type
    use tblite_lapack_solver, only : lapack_solver
    use tblite_cusolver_sygvd, only : new_sygvd_gpu, sygvd_cusolver
    use tblite_lapack_sygvd, only : new_sygvd, sygvd_solver
    use tblite_lapack_solver, only : solver_type
    use tblite_timer, only : timer_type, format_time
    use tblite_wavefunction_type, only : wavefunction_type
    use gambits, only : gambits_context_type
    use iso_c_binding
    implicit none
    private

    public :: new_purification

    public :: purification_precision, purification_type, purification_runmode

    type enum_purification_type
        !> no puridfication
        integer(c_size_t) :: none = 0
        !> TC2 purification
        integer(c_size_t) :: tc2 = 2
        !> accelerated SP2
        integer(c_size_t) :: tc2accel = 12
        !> Mc Weeney Purification
        integer(c_size_t) :: mcweeney = 3
        !> TRS4 purification
        integer(c_size_t) :: trs4 = 4
    end type

    type enum_precision
        !> mixed precision
        integer(c_size_t) :: mixed = 0
        !> full double precision
        integer(c_size_t) :: double = 1
        !> full single precision
        integer(c_size_t) :: single = 2
        !> approx FP64
        integer(c_size_t) :: approx = 3
    end type

    type enum_runmode
        !> default switching CPU and GPU
        integer(c_size_t) :: default = 0
        !> all on the CPU
        integer(c_size_t) :: cpu = 1
        !> all GPU
        integer(c_size_t) :: gpu = 2
    end type
    
   type(enum_precision), parameter :: purification_precision = enum_precision()
   type(enum_purification_type), parameter :: purification_type = enum_purification_type()
   type(enum_runmode), parameter :: purification_runmode = enum_runmode()


    type, public, extends(solver_type) :: purification_solver
      integer(c_size_t) :: precision = purification_precision%mixed
      integer(c_size_t) :: type = purification_type%tc2
      integer(c_size_t) :: runmode = purification_runmode%default
      real(c_double) :: thresh = 5.0e-07_c_double
      type(c_ptr) :: solver_ptr
      integer(c_size_t) :: maxiter = 100
      integer :: iscf = 0
      class(solver_type), allocatable :: lapack_solv
      logical :: trans = .false.
      logical :: purification_success = .true.
      type(gambits_context_type) :: ctx 
   contains
      procedure :: solve_sp
      procedure :: solve_dp
      procedure :: purify_dp
      procedure :: got_transform
      procedure :: delete
      procedure :: get_density_matrix
      procedure :: get_energy_w_density_matrix
      procedure :: reset
   end type

   interface
      type(c_ptr) function PurifyFockSetup(ctx, nmo, type, run, prec) bind(C, name="PurifyFockSetup")
         use iso_c_binding
         type(c_ptr), value :: ctx
         integer(c_size_t), value, intent(in) :: nmo
         integer(c_size_t), value, intent(in) :: type, run, prec
      end function
      subroutine SetTransformationMatrix(ctx, ptr, X) bind(C, name="SetTransformationMatrix")
         use iso_c_binding
         type(c_ptr), value :: ctx
         type(c_ptr), value :: ptr
         real(c_double),  intent(in) :: X(*)
      end subroutine
      subroutine GetDensityAO(ctx, ptr, Fock, Dens, elnum, info) bind(C, name="GetDensityAO")
         use iso_c_binding
         type(c_ptr), value :: ctx
         type(c_ptr), value :: ptr
         real(c_double) :: Fock(*)
         real(c_double) :: Dens(*)
         real(c_double), value :: elnum
         integer(kind=4) :: info
      end subroutine
      subroutine GetEnergyWDensityAO(ctx, ptr, Fock, Dens, WDens) bind(C, name="GetEnergyWDensityAO")
         use iso_c_binding
         type(c_ptr), value :: ctx
         type(c_ptr), value :: ptr
         real(c_double) :: Fock(*)
         real(c_double) :: Dens(*)
         real(c_double) :: WDens(*)
      end subroutine
      subroutine DeletePointer(ptr) bind(C, name="DeletePointer")
         use iso_c_binding
         type(c_ptr), value :: ptr
      end subroutine
      subroutine ResetLib(ptr) bind(C, name="Reset")
         use iso_c_binding
         type(c_ptr), value :: ptr
      end subroutine
   end interface

   contains

   subroutine reset(self)
      class(purification_solver) :: self
      self%trans = .false.
      self%purification_success = .true.
      call ResetLib(self%solver_ptr)
   end subroutine
   

   function got_transform(self) result(trans)
      class(purification_solver) :: self
      logical :: trans
      trans = self%trans;
   end function

   subroutine new_purification(self, type_, run_, prec_, ndim, maxiter, thresh)
      type(purification_solver) :: self
      integer(c_size_t) :: type_, run_, prec_
      integer :: ndim
      integer(c_size_t), optional :: maxiter
      real(c_double), optional :: thresh
      self%type = type_
      self%precision = prec_
      self%runmode = run_

      !use LAPACK for molecules smaller than 750 basis functions
      select case(run_)
      case(purification_runmode%cpu)
         block
            type(sygvd_solver), allocatable :: tmp
            allocate(tmp) 
            call new_sygvd(tmp, ndim)
            call move_alloc(tmp, self%lapack_solv)
         end block
      case(purification_runmode%gpu)
         block
            type(sygvd_cusolver), allocatable :: tmp
            allocate(tmp) 
            call new_sygvd_gpu(tmp, ndim, .true.)
            call move_alloc(tmp, self%lapack_solv)
         end block
      case(purification_runmode%default)
         if (ndim < 750) then
            block
               type(sygvd_solver), allocatable :: tmp
               allocate(tmp) 
               call new_sygvd(tmp, ndim)
               call move_alloc(tmp, self%lapack_solv)
            end block
         else
            block
               type(sygvd_cusolver), allocatable :: tmp
               allocate(tmp) 
               call new_sygvd_gpu(tmp, ndim, .true.)
               call move_alloc(tmp, self%lapack_solv)
            end block
            end if
      end select


      
      
      if (present(maxiter)) self%maxiter = maxiter
      if (present(thresh)) self%thresh = thresh
      call self%ctx%setup(int(1, kind=c_size_t), self%maxiter, self%thresh)
      self%solver_ptr =  PurifyFockSetup(self%ctx%ptr, int(ndim, kind=c_size_t), &
      self%type, self%runmode, self%precision)
   end subroutine

   subroutine solve_sp(self, hmat, smat, eval, error)
      class(purification_solver), intent(inout) :: self
      real(sp), contiguous, intent(inout) :: hmat(:, :)
      real(sp), contiguous, intent(in) :: smat(:, :)
      real(sp), contiguous, intent(inout) :: eval(:)
      type(error_type), allocatable, intent(out) :: error
   end subroutine



   subroutine solve_dp(self, hmat, smat, eval, error)
      class(purification_solver), intent(inout) :: self
      real(dp), contiguous, intent(inout) :: hmat(:, :)
      real(dp), contiguous, intent(in) :: smat(:, :)
      real(dp), contiguous, intent(inout) :: eval(:)
      real(c_double), allocatable, target :: tmp(:, :)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: error_msg
      if (self%ctx%failed()) then
         call self%ctx%get_error(error_msg)
         call fatal_error(error, error_msg)
      end if
      call self%lapack_solv%solve(hmat, smat, eval, error)
      allocate(tmp(size(hmat, dim=1), size(hmat, dim=2)), source=0.0_c_double)
      tmp = hmat
      if (self%purification_success) call SetTransformationMatrix(self%ctx%ptr, self%solver_ptr, tmp)
      self%trans = .true.
   end subroutine solve_dp

   subroutine purify_dp(self, hmat, smat, dens, nel, error)
        class(purification_solver), intent(inout) :: self
        real(c_double), contiguous, intent(inout)  :: hmat(:, :)
        real(c_double), contiguous, intent(in) :: smat(:, :)
        real(c_double), contiguous, intent(inout) :: dens(:, :)
        real(c_double) :: nel
        integer(kind=4) :: info
        type(error_type), allocatable, intent(out) :: error
        character(len=:), allocatable :: error_msg
        character(len=:), allocatable :: msg

         self%iscf = self%iscf +1
         info = 42
         call GetDensityAO(self%ctx%ptr, self%solver_ptr, hmat, dens, nel/2, info)
         if (info /= 0) self%purification_success = .false.
         call self%ctx%get_message(msg)
         write(*,*) msg
         if (self%ctx%failed()) then
            call self%ctx%get_error(error_msg)
            call fatal_error(error, error_msg)
         end if
         
    end subroutine purify_dp

    subroutine delete(self)
      class(purification_solver) :: self
      call self%lapack_solv%delete() 
      deallocate(self%lapack_solv)
      call self%ctx%delete()
      call DeletePointer(self%solver_ptr)
      self%solver_ptr = c_null_ptr
    end subroutine

    subroutine get_density_matrix(self, focc, coeff, pmat)
      class(purification_solver) :: self
      real(wp), intent(in) :: focc(:)
      real(wp), contiguous, intent(in) :: coeff(:, :)
      real(wp), contiguous, intent(out) :: pmat(:, :)
   
      real(wp), allocatable :: scratch(:, :)
      integer :: iao, jao
      call self%lapack_solv%get_density_matrix(focc, coeff, pmat)
   end subroutine get_density_matrix


   subroutine get_energy_w_density_matrix(self, wfn, wdensity)
      class(purification_solver) :: self
      type(wavefunction_type), intent(inout) :: wfn
      real(wp) :: wdensity(:,:,:)
      real(wp), allocatable :: tmp(:)
      integer :: spin
      if (self%purification_success) then
         do spin = 1, wfn%nspin
            call GetEnergyWDensityAO(self%ctx%ptr, self%solver_ptr, wfn%coeff(:, :, spin), wfn%density(:, :, spin), wdensity(:, :, spin))         
         end do
      else
         call self%lapack_solv%get_energy_w_density_matrix(wfn, wdensity)
      end if
   end subroutine


end module
