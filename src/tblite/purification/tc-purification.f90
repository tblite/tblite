module tblite_purification_solver
    use mctc_env, only : sp, dp, error_type, wp
    use tblite_scf_solver, only : solver_type
    use tblite_lapack_solver, only : lapack_solver
    use tblite_lapack_sygvd, only : new_sygvd, sygvd_solver
    use tblite_blas, only : gemm
    use tblite_timer, only : timer_type, format_time
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
      type(sygvd_solver), allocatable :: lapack_solv
      real(c_double), allocatable :: X(:, :)
      type(timer_type) :: timer
   contains
      procedure :: solve_sp
      procedure :: solve_dp
      procedure :: purify_dp
      procedure :: got_transform
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
      subroutine GetDensityGuess(nao, efermi, hmax, hmin, nel, Fock, Dens) bind(C, name="GetDensityGuess")
         use iso_c_binding
         integer(c_size_t), value, intent(in) :: nao
         real(c_double), value, intent(in) :: nel, efermi, hmax, hmin
         type(c_ptr), value, intent(in) :: Fock
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

   function got_transform(self) result(trans)
      class(purification_solver) :: self
      logical :: trans
      trans = .false.
      if (allocated(self%X)) trans = .true.
   end function

   subroutine new_purification(self, type_, run_, prec_, ndim, maxiter, thresh)
      type(purification_solver) :: self
      integer(c_size_t) :: type_, run_, prec_
      integer :: ndim
      integer(c_size_t), optional :: maxiter
      real(c_double), optional :: thresh
      type(sygvd_solver), allocatable :: tmp
      
      allocate(self%lapack_solv)
      call new_sygvd(self%lapack_solv, ndim)
      
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
      call self%timer%push("Diagonalize")
      call self%lapack_solv%solve(hmat, smat, eval, error)
      allocate(self%X(size(hmat, dim=1), size(hmat, dim=2)), source=0.0_dp)
      self%X = hmat
      call self%timer%pop()
   end subroutine solve_dp

   subroutine purify_dp(self, hmat, smat, dens, nel, error)
        class(purification_solver), intent(inout) :: self
        real(c_double), contiguous, target, intent(inout)  :: hmat(:, :)
        real(c_double), contiguous, target, intent(in) :: smat(:, :)
        real(c_double), contiguous, target, intent(inout) :: dens(:, :)
        real(c_double), allocatable :: tmp(:,:)
        real(c_double), contiguous, pointer  :: Fptr(:, :)
        real(c_double), contiguous, pointer :: Sptr(:, :)
        real(c_double), contiguous, pointer :: Dptr(:, :)
        integer(c_size_t) :: nmo
        real(c_double) :: nel
        type(error_type), allocatable, intent(out) :: error
        real(c_double) :: cpot, hmax, hmin 
        type(c_ptr) :: F, S, D
         call self%timer%push("TransformF")
         allocate(tmp(size(hmat, dim=1), size(hmat, dim=2)))
         call gemm(hmat, self%X, tmp)
         call gemm(self%X, tmp, hmat, transa='t')
         call self%timer%pop()
         Fptr => hmat
         Sptr => smat
         Dptr => dens

         F = c_loc(Fptr)
         S = c_loc(Sptr)
         D = c_loc(Dptr)

         nmo = size(smat,dim=1)
         call self%timer%push("Guess")
         call GetBoundsGershgorin(nmo, F, cpot, hmax, hmin)
         call GetDensityGuess(nmo,cpot, hmax, hmin, nel/2, F, D)
         call self%timer%pop()
         self%iscf = self%iscf +1
         call self%timer%push("Purification") 
         call PurifyDensity(nmo, nel/2, self%thresh, self%maxiter, D, self%type, self%runmode, self%precision)
         call self%timer%pop()
         call self%timer%push("TransformD")
         call gemm(dens, self%X, tmp, alpha=1.0_dp, beta=0.0_dp, transb='t')
         call gemm(self%X, tmp, dens, alpha=2.0_dp)
         call self%timer%pop()
         block
         integer :: it
         real(wp) :: stime
         character(len=*), parameter :: label(*) = [character(len=20):: &
         & "TransformF", "Guess", "Purification", "TransformD", "Diagonalize"]
         do it = 1, size(label)
            stime = self%timer%get(label(it))
            if (stime <= epsilon(0.0_wp)) cycle
            write(*,*) " - "//label(it)//format_time(stime)
         end do            

         end block
    end subroutine purify_dp

end module
