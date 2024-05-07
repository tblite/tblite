module tblite_purification_solver_context
   use mctc_env, only : sp, dp, error_type, fatal_error
    use tblite_context_solver, only: context_solver
    use tblite_scf_solver, only: solver_type
    use tblite_purification_solver, only: purification_solver, new_purification, purification_runmode, purification_precision, purification_type
    use iso_c_binding
    implicit none
    private

    public :: purification_precision, purification_type, purification_runmode


    !> Generator for purification based electronic solvers
   type, public, extends(context_solver) :: purification_solver_context
      !> Selected electronic solver algorithm
      integer(c_size_t) :: type = purification_type%tc2
      integer(c_size_t) :: precision = purification_precision%mixed
      integer(c_size_t) :: runmode = purification_runmode%default
   contains
      !> Create new instance of electronic solver
      procedure :: new
      !> Delete an electronic solver instance
      procedure :: delete
   end type purification_solver_context

contains
   !> Create new electronic solver
subroutine new(self, solver, ndim)
   !> Instance of the solver factory
   class(purification_solver_context), intent(inout) :: self
   !> New electronic solver
   class(solver_type), allocatable, intent(out) :: solver
   !> Dimension of the eigenvalue problem
   integer, intent(in) :: ndim
   type(purification_solver), allocatable :: tmp

   allocate(tmp)
   call new_purification(tmp, self%type, self%runmode, self%precision, ndim)
   call move_alloc(tmp, solver)
   
   
end subroutine new

!> Delete electronic solver instance
subroutine delete(self, solver)
   !> Instance of the solver factory
   class(purification_solver_context), intent(inout) :: self
   !> Electronic solver instance
   class(solver_type), allocatable, intent(inout) :: solver
   call solver%delete()
   if (allocated(solver)) deallocate(solver)
end subroutine delete
    

end module