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
    use tblite_cusolver_sygvd, only : new_sygvd_gpu, sygvd_cusolver
    use tblite_lapack_sygvd, only : new_sygvd, sygvd_solver
    use tblite_timer, only : timer_type, format_time
    use tblite_wavefunction_type, only : wavefunction_type
#if TBLITE_HAS_GAMBITS
    use gambits, only : gambits_context_type
#endif
    use tblite_features, only : tblite_has_cuda
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
        !> mixed precision with Ozaki-type FP16 splits
        integer(c_size_t) :: mixed_fp16 = 3
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

   type, public :: dmp_input
      !> Type/Algorithm of purification
      integer(c_size_t) :: type = purification_type%tc2 
      !> Numerical precision of the DMP solver
      integer(c_size_t) :: precision = purification_precision%mixed
      !> Run mode of the solver, default is to switch to GPU for ndim > 750
      integer(c_size_t) :: runmode = purification_runmode%default
   end type

   type, public, extends(solver_type) :: purification_solver
      type(dmp_input) :: input
      real(c_double) :: thresh = 5.0e-07_c_double
      type(c_ptr) :: solver_ptr = c_null_ptr
      integer :: ndim = 0
      integer(c_size_t) :: maxiter = 100
      integer :: iscf = 0
      class(solver_type), allocatable :: lapack_solv
      logical :: transform = .false.
      logical :: purification_success = .true.
      logical :: needs_orbitals = .false.
#if TBLITE_HAS_GAMBITS
      type(gambits_context_type) :: ctx = gambits_context_type(ptr=c_null_ptr)
#endif
   contains
      procedure :: get_density
      procedure :: finish_wavefunction
      procedure :: get_wdensity
      procedure :: got_transform
      procedure :: delete
      procedure :: reset
      final :: finalize
   end type

   interface
      type(c_ptr) function SetupPurification(ctx, nmo, type, run, prec) bind(C, name="SetupPurification")
         use iso_c_binding
         type(c_ptr), value :: ctx
         integer(c_size_t), value, intent(in) :: nmo
         integer(c_int), value, intent(in) :: type, run, prec
      end function
      subroutine SetTransformationMatrix(ctx, ptr, X) bind(C, name="SetTransformationMatrix")
         use iso_c_binding
         type(c_ptr), value :: ctx
         type(c_ptr), value :: ptr
         real(c_double),  intent(in) :: X(*)
      end subroutine
      subroutine GetDensityAO(ctx, ptr, Fock, Dens, elnum) bind(C, name="GetDensityAO")
         use iso_c_binding
         type(c_ptr), value :: ctx
         type(c_ptr), value :: ptr
         real(c_double) :: Fock(*)
         real(c_double) :: Dens(*)
         real(c_double), value :: elnum
      end subroutine
      subroutine DeletePurification(ptr) bind(C, name="DeletePurification")
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
   self%transform = .false.
   self%purification_success = .true.
#if TBLITE_HAS_GAMBITS
   call ResetLib(self%solver_ptr)
#endif
end subroutine
   

function got_transform(self) result(trans)
   class(purification_solver) :: self
   logical :: trans
   trans = self%transform;
end function

subroutine new_purification(self, overlap, nel, kt, dmp_inp)
   type(purification_solver), intent(inout) :: self
   real(wp), intent(in) :: overlap(:, :), nel(:), kt
   type(dmp_input), intent(in) :: dmp_inp
   logical :: use_gpu
   integer(c_int) :: runmode

   self%ndim = size(overlap, 1)
   self%input = dmp_inp
   self%nel = nel
   self%kt = kt
   self%needs_orbitals = .false.

#if TBLITE_HAS_GAMBITS
   if (kt /= 0.0_wp) return
   if (dmp_inp%runmode == purification_runmode%gpu .and. .not.tblite_has_cuda) return
   ! A new geometry needs a fresh overlap factorization, even at fixed AO count.
   if (allocated(self%lapack_solv)) deallocate(self%lapack_solv)
   use_gpu = dmp_inp%runmode == purification_runmode%gpu .or. &
      & (dmp_inp%runmode == purification_runmode%default .and. tblite_has_cuda .and. self%ndim >= 750)
   runmode = int(merge(purification_runmode%gpu, purification_runmode%cpu, use_gpu), c_int)
   if (use_gpu) then
      allocate(sygvd_cusolver :: self%lapack_solv)
      select type (diag => self%lapack_solv)
      type is (sygvd_cusolver)
         call new_sygvd_gpu(diag, overlap, nel, 0.0_wp)
      end select
   else
      allocate(sygvd_solver :: self%lapack_solv)
      select type (diag => self%lapack_solv)
      type is (sygvd_solver)
         call new_sygvd(diag, overlap, nel, 0.0_wp)
      end select
   end if
   if (c_associated(self%solver_ptr)) then
      call self%reset()
   else
      call self%ctx%setup(1_c_size_t, self%maxiter, self%thresh)
      self%solver_ptr = SetupPurification(self%ctx%ptr, int(self%ndim,c_size_t), &
         & int(dmp_inp%type,c_int), runmode, int(dmp_inp%precision,c_int))
   end if
#endif
end subroutine

subroutine get_density(self, hmat, smat, eval, focc, density, error)
   class(purification_solver), intent(inout) :: self
   !> Overlap matrix
   real(wp), contiguous, intent(in) :: smat(:, :)
   !> Hamiltonian matrix, can contains eigenvectors on output
   real(wp), contiguous, intent(inout) :: hmat(:, :, :)
   !> Eigenvalues
   real(wp), contiguous, intent(inout) :: eval(:, :)
   !> Occupation numbers
   real(wp), contiguous, intent(inout) :: focc(:, :)
   !> Density matrix
   real(wp), contiguous, intent(inout) :: density(:, :, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   character(len=:), allocatable :: error_msg
   real(c_double), allocatable, target :: tmp(:, :, :)
   real(c_double) :: nel, nel1
   integer :: nspin, spin

#if TBLITE_HAS_GAMBITS
   if (self%kt /= 0.0_wp) then
      call fatal_error(error, "Density matrix purification requires zero electronic temperature")
      return
   end if
   if (self%input%runmode == purification_runmode%gpu .and. .not.tblite_has_cuda) then
      call fatal_error(error, "GPU purification requires CUDA support")
      return
   end if
   nspin = size(hmat, dim=3)

   if (.not. self%transform .or. (.not. self%purification_success)) then
      if (self%ctx%failed()) then
         call self%ctx%get_error(error_msg)
         call fatal_error(error, error_msg)
         return
      end if
      call self%lapack_solv%get_density(hmat, smat, eval, focc, density, error)
      if (allocated(error)) return
      if (self%purification_success) call SetTransformationMatrix(self%ctx%ptr, self%solver_ptr, hmat(:, :, 1))
      self%transform = .true.
   else
      self%needs_orbitals = .true.
      select case(nspin)
      case(1)
         nel = sum(self%nel)
         if (mod(int(nel), 2) /= 0) then
            nel1 = ceiling(nel/2)
            allocate(tmp(size(density, dim=1), size(density, dim=2), nspin), source=0.0_c_double)
            call GetDensityAO(self%ctx%ptr, self%solver_ptr, hmat, tmp, nel1)
            density = tmp
            call GetDensityAO(self%ctx%ptr, self%solver_ptr, hmat, tmp, nel-nel1)
            density = density + tmp
         else
            call GetDensityAO(self%ctx%ptr, self%solver_ptr, hmat, density(:, :, 1), nel/2)
            density = 2 * density 
         end if
      case(2)
         do spin = 1, nspin
            nel = self%nel(spin)
            call GetDensityAO(self%ctx%ptr, self%solver_ptr, hmat(:, :, spin), density(:, :, spin), nel)
         end do
      end select   
   end if

   if (self%ctx%failed()) then
      call self%ctx%get_error(error_msg)
      write(*,*) error_msg
      self%purification_success = .false.
      self%needs_orbitals = .false.
      call self%lapack_solv%get_density(hmat, smat, eval, focc, density, error)
   end if
#else
   call fatal_error(error, "Density matrix purification requires a build with GAMBITS")
#endif
end subroutine

!> Recover eigenvectors and occupations once from the final SCF Hamiltonian.
subroutine finish_wavefunction(self, wfn, smat, error)
   class(purification_solver), intent(inout) :: self
   type(wavefunction_type), intent(inout) :: wfn
   real(wp), contiguous, intent(in) :: smat(:, :)
   type(error_type), allocatable, intent(out) :: error

   if (.not.self%needs_orbitals) return
   call self%lapack_solv%get_density(wfn%coeff, smat, wfn%emo, wfn%focc, wfn%density, error)
   self%needs_orbitals = .false.
end subroutine finish_wavefunction

subroutine delete(self)
   class(purification_solver), intent(inout) :: self
#if TBLITE_HAS_GAMBITS
   if (c_associated(self%solver_ptr)) then
      call DeletePurification(self%solver_ptr)
      self%solver_ptr = c_null_ptr
   end if
   if (c_associated(self%ctx%ptr)) then
      call self%ctx%delete()
      self%ctx%ptr = c_null_ptr
   end if
#endif
   if (allocated(self%lapack_solv)) deallocate(self%lapack_solv)
end subroutine

subroutine finalize(self)
   type(purification_solver), intent(inout) :: self
   call self%delete()
end subroutine

subroutine get_wdensity(self, hmat, smat, eval, focc, density, error)
   class(purification_solver), intent(inout) :: self
   !> Overlap matrix
   real(wp), contiguous, intent(in) :: smat(:, :)
   !> Hamiltonian matrix, can contains eigenvectors on output
   real(wp), contiguous, intent(inout) :: hmat(:, :, :)
   !> Eigenvalues
   real(wp), contiguous, intent(inout) :: eval(:, :)
   !> Occupation numbers
   real(wp), contiguous, intent(inout) :: focc(:, :)
   !> Density matrix
   real(wp), contiguous, intent(inout) :: density(:, :, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call self%lapack_solv%get_wdensity(hmat, smat, eval, focc, density, error)
end subroutine


end module
