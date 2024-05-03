module tblite_purification_solver
    use mctc_env, only : sp, dp, error_type
    use tblite_scf_solver, only : solver_type
    use iso_c_binding
    implicit none
    private

    public :: new_purification

    public :: purification_precision, purification_type, purification_runmode

    type enum_purification_type
        !> no puridfication
        integer(c_size_t) :: none = 0
        !> TC2 purification
        integer(c_size_t) :: tc2 = 1
        !> TRS4 purification
        integer(c_size_t) :: trs4 = 2
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
      real(c_double) :: thresh = 1.0e-05
      integer(c_size_t) :: maxiter = 200
      integer :: iscf = 0
   contains
      procedure :: solve_sp
      procedure :: solve_dp
      procedure :: purify_dp
   end type

   interface
      subroutine GetBoundsGershgorin(n, H, cpot, hmax, hmin) bind(C, name="GetBoundsGershgorin")
         use iso_c_binding
         integer(c_size_t), value, intent(in) :: n
         type(c_ptr), value, intent(in) :: H
         real(c_double), intent(out) :: cpot, hmax, hmin
      end subroutine
      subroutine GetDensityGuessAO(nao, nel, Fock, Overlap, Dens) bind(C, name="GetDensityGuessAO")
         use iso_c_binding
         integer(c_size_t), value, intent(in) :: nao
         real(c_double), value, intent(in) :: nel
         type(c_ptr), value, intent(in) :: Fock
         type(c_ptr), value, intent(in) :: Overlap
         type(c_ptr), value :: Dens
      end subroutine
      subroutine PurifyDensity(nmo, nel, threshPP, maxit, Dens, type, run, prec) bind(C, name="PurifyDensity")
         use iso_c_binding
         integer(c_size_t), value, intent(in) :: nmo
         real(c_double), value, intent(in) :: nel
         real(c_double), value, intent(in) :: threshPP
         integer(c_size_t), value, intent(in) :: maxit
         type(c_ptr), value :: Dens
         integer(c_size_t), value, intent(in) :: type, run, prec
      end subroutine
   end interface

   contains

   subroutine new_purification(self, type_, run_, prec_, maxiter, thresh)
      type(purification_solver) :: self
      integer(c_size_t) :: type_, run_, prec_
      integer(c_size_t), optional :: maxiter
      real(c_double), optional :: thresh

      self%type = type_
      self%precision = prec_
      self%runmode = run_

      if (present(maxiter)) self%maxiter = maxiter
      if (present(thresh)) self%thresh = thresh

   end subroutine


    subroutine solve_sp(self, hmat, smat, eval, error)
        class(purification_solver), intent(inout) :: self
        real(sp), contiguous, intent(inout) :: hmat(:, :)
        real(sp), contiguous, intent(in) :: smat(:, :)
        real(sp), contiguous, intent(inout) :: eval(:)
        type(error_type), allocatable, intent(out) :: error


    end subroutine solve_sp

    subroutine solve_dp(self, hmat, smat, eval, error)
      class(purification_solver), intent(inout) :: self
      real(dp), contiguous, intent(inout) :: hmat(:, :)
      real(dp), contiguous, intent(in) :: smat(:, :)
      real(dp), contiguous, intent(inout) :: eval(:)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: cpot, hmax, hmin 

   end subroutine solve_dp

   subroutine purify_dp(self, hmat, smat, dens, nel, error)
        class(purification_solver), intent(inout) :: self
        real(c_double), contiguous, target, intent(inout)  :: hmat(:, :)
        real(c_double), contiguous, target, intent(in) :: smat(:, :)
        real(c_double), contiguous, target, intent(inout) :: dens(:, :)
        real(c_double), contiguous, pointer  :: Fptr(:, :)
        real(c_double), contiguous, pointer :: Sptr(:, :)
        real(c_double), contiguous, pointer :: Dptr(:, :)
        integer(c_size_t) :: nmo
        real(c_double) :: nel
        type(error_type), allocatable, intent(out) :: error
        real(c_double) :: cpot, hmax, hmin 
        type(c_ptr) :: F, S, D

         Fptr => hmat
         Sptr => smat
         Dptr => dens

         F = c_loc(Fptr)
         S = c_loc(Sptr)
         D = c_loc(Dptr)

         nmo = size(smat,dim=1)

        call GetBoundsGershgorin(nmo, F, cpot, hmax, hmin)
         call GetDensityGuessAO(nmo, nel/2, F, S, D)
         self%iscf = self%iscf +1 
        call PurifyDensity(nmo, nel/2, self%thresh, self%maxiter, D, self%type, self%runmode, self%precision)

    end subroutine purify_dp

end module
