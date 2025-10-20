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

    public :: new_purification, gambits_context_type, DeletePurification

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
      type(c_ptr) :: solver_ptr
      integer(c_size_t) :: maxiter = 100
      integer :: iscf = 0
      class(solver_type), allocatable :: lapack_solv
      logical :: transform = .false.
      logical :: purification_success = .true.
      type(gambits_context_type) :: ctx 
      !> Pointer to the C++ SYGVD instance
      type(c_ptr) :: sygvd_ptr = c_null_ptr
   contains
      procedure :: get_density
      procedure :: get_wdensity
      procedure :: got_transform
      procedure :: delete
      procedure :: reset
   end type

   interface
      type(c_ptr) function SetupPurification(ctx, nmo, type, run, prec) bind(C, name="SetupPurification")
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
      subroutine GetDensityAO(ctx, ptr, Fock, Dens, elnum) bind(C, name="GetDensityAO")
         use iso_c_binding
         type(c_ptr), value :: ctx
         type(c_ptr), value :: ptr
         real(c_double) :: Fock(*)
         real(c_double) :: Dens(*)
         real(c_double), value :: elnum
      end subroutine
      subroutine GetEnergyWDensityAO(ctx, ptr, Fock, Dens, WDens) bind(C, name="GetEnergyWDensityAO")
         use iso_c_binding
         type(c_ptr), value :: ctx
         type(c_ptr), value :: ptr
         real(c_double) :: Fock(*)
         real(c_double) :: Dens(*)
         real(c_double) :: WDens(*)
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
   call ResetLib(self%solver_ptr)
end subroutine
   

function got_transform(self) result(trans)
   class(purification_solver) :: self
   logical :: trans
   trans = self%transform;
end function

subroutine new_purification(self, overlap, nel, kt, dmp_inp, dmp_ptr, gvd_ptr)
   !> Create a new purification solver instance
   type(purification_solver), intent(inout) :: self
   !> Overlap matrix
   real(wp), contiguous, intent(in) :: overlap(:, :)
   !> Number of electrons
   real(wp), contiguous, intent(in) :: nel(:)
   !> kt Term
   real(wp), intent(in) :: kt
   !> Type of purification, runmode and precision
   type(dmp_input), intent(in) :: dmp_inp
   !> Pointer to the C++ DMP instance
   type(c_ptr), intent(inout) :: dmp_ptr
   !> Pointer to the C++ SGVD instance
   type(c_ptr), intent(inout) :: gvd_ptr
   integer :: ndim
   ndim = size(overlap, dim=1)
   self%input = dmp_inp
   !use LAPACK for molecules smaller than 750 basis functions
   if (.not. allocated(self%lapack_solv)) then
   select case(dmp_inp%runmode)
   case(purification_runmode%cpu)
      block
         type(sygvd_solver), allocatable :: tmp
         allocate(tmp) 
         call new_sygvd(tmp, overlap, nel, kt)
         call move_alloc(tmp, self%lapack_solv)
      end block
   case(purification_runmode%gpu)
      block
         type(sygvd_cusolver), allocatable :: tmp
         allocate(tmp) 
         call new_sygvd_gpu(tmp, overlap, nel, kt, gvd_ptr, .true.)
         call move_alloc(tmp, self%lapack_solv)
      end block
   case(purification_runmode%default)
      if (ndim < 750) then
         block
            type(sygvd_solver), allocatable :: tmp
            allocate(tmp) 
            call new_sygvd(tmp, overlap, nel, kt)
            call move_alloc(tmp, self%lapack_solv)
         end block
      else
         block
            type(sygvd_cusolver), allocatable :: tmp
            allocate(tmp) 
            call new_sygvd_gpu(tmp, overlap, nel, kt, gvd_ptr, .true.)
            call move_alloc(tmp, self%lapack_solv)
         end block
         end if
   end select
   call self%ctx%setup(int(1, kind=c_size_t), self%maxiter, self%thresh)
   end if
   
   if (.not. c_associated(dmp_ptr)) then
      dmp_ptr = SetupPurification(self%ctx%ptr, int(ndim, kind=c_size_t), self%input%type, self%input%runmode, self%input%precision)
   else
      call ResetLib(dmp_ptr)
   end if
   self%nel = nel
   self%solver_ptr = dmp_ptr
   self%kt = 0.0_wp
   if (c_associated(gvd_ptr)) self%sygvd_ptr = gvd_ptr
   
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

   nspin = size(hmat, dim=3)

   if (.not. self%transform .or. (.not. self%purification_success)) then
      if (self%ctx%failed()) then
         call self%ctx%get_error(error_msg)
         call fatal_error(error, error_msg)
      end if
      call self%lapack_solv%get_density(hmat, smat, eval, focc, density, error)
      allocate(tmp(size(hmat, dim=1), size(hmat, dim=2), 1), source=0.0_c_double)
      tmp = hmat
      if (self%purification_success) call SetTransformationMatrix(self%ctx%ptr, self%solver_ptr, tmp)
      self%transform = .true.
   else
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
            call GetDensityAO(self%ctx%ptr, self%solver_ptr, hmat, density(:, :, spin), nel)
         end do
      end select   
   end if

   if (self%ctx%failed()) then
      call self%ctx%get_error(error_msg)
      write(*,*) error_msg
      self%purification_success = .false.
      call self%lapack_solv%get_density(hmat, smat, eval, focc, density, error)
   end if

end subroutine

subroutine delete(self, ptr)
   class(purification_solver), intent(inout) :: self
   !> Pointer to the C++ solver instance
   type(c_ptr), intent(inout), optional :: ptr

   if (present(ptr)) then
      if (.not.c_associated(self%solver_ptr)) return
      if (c_associated(ptr)) call DeletePurification(ptr)
      if (c_associated(self%sygvd_ptr)) call self%lapack_solv%delete(self%sygvd_ptr)
      self%sygvd_ptr = c_null_ptr
      self%solver_ptr = c_null_ptr
      if (allocated(self%lapack_solv)) deallocate(self%lapack_solv)
   else
      call self%reset()
      
   end if

   
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

   real(c_double), allocatable, target :: tmp(:, :, :)
   
   integer :: spin

   tmp = density

   if (self%purification_success) then
      do spin = 1, size(hmat, dim=3)
         call GetEnergyWDensityAO(self%ctx%ptr, self%solver_ptr, hmat(:, :, spin), tmp(:, :, spin),&
          density(:, :, spin))         
      end do
   else
      call self%lapack_solv%get_wdensity(hmat, smat, eval, focc, density, error)
   end if
end subroutine


end module
