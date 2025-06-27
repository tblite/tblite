module tblite_purification_solver_context
   use mctc_env, only : sp, dp, error_type, fatal_error, wp
    use tblite_context_solver, only: context_solver
    use tblite_scf_solver, only: solver_type
    use tblite_purification_solver, only: purification_solver, new_purification
    use tblite_purification_solver, only: dmp_input, purification_precision, purification_type, purification_runmode
    use iso_c_binding
    implicit none
    private

    public :: purification_precision, purification_type, purification_runmode, dmp_input


    !> Generator for purification based electronic solvers
   type, public, extends(context_solver) :: purification_solver_context
      type(dmp_input) :: dmp_input
      !> Pointer to store the C++ DMP instance
      type(c_ptr) :: dmp_ptr = c_null_ptr
      !> Pointer to store the C++ SYGVD instance
      type(c_ptr) :: sgvd_ptr = c_null_ptr
      !> Reuse the solver instance
      logical :: reuse = .false.
   contains
      !> Create new instance of electronic solver
      procedure :: new
      !> Delete an electronic solver instance
      procedure :: delete
   end type purification_solver_context

contains
   !> Create new electronic solver
subroutine new(self, solver, overlap, nel, kt)
   !> Instance of the solver factory
   class(purification_solver_context), intent(inout) :: self
   !> New electronic solver
   class(solver_type), allocatable, intent(out) :: solver
   !> Overlap matrix
   real(wp), intent(in) :: overlap(:, :)
   !> Number of electrons per spin channel
   real(wp), intent(in) :: nel(:)
   !> Electronic temperature
   real(wp), intent(in) :: kt
   type(purification_solver), allocatable :: tmp
   
   allocate(tmp)
   call new_purification(tmp, overlap, nel,  kt,  self%dmp_input, self%dmp_ptr, self%sgvd_ptr)
   call move_alloc(tmp, solver)
   
   
end subroutine new

!> Delete electronic solver instance
subroutine delete(self, solver)
   !> Instance of the solver factory
   class(purification_solver_context), intent(inout) :: self
   !> Electronic solver instance
   class(solver_type), allocatable, intent(inout) :: solver
   
   
   if (allocated(solver)) then
      if (self%reuse .and. c_associated(self%dmp_ptr)) then
         call solver%delete()
      else
         write(*, '(a)') "Deleting solver instance"
         call solver%delete(self%dmp_ptr)
         self%dmp_ptr = c_null_ptr
         self%sgvd_ptr = c_null_ptr
      end if
      deallocate(solver)
   end if
   
end subroutine delete
    

end module