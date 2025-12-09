module tblite_purification_solver_context
   use mctc_env, only : sp, dp, error_type, fatal_error, wp
    use tblite_context_solver, only: context_solver
    use tblite_scf_solver, only: solver_type
    use tblite_purification_solver, only: purification_solver, new_purification, gambits_context_type, DeletePurification
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
   contains
      !> Create new instance of electronic solver
      procedure :: new
      !> Delete an electronic solver instance
      procedure :: delete
   end type purification_solver_context

   interface purification_solver_context
      module procedure new_purification_solver_context
   end interface purification_solver_context

contains

type(purification_solver_context) function new_purification_solver_context(dmp_inp, reuse)
   !> Input parameters for the DMP solver
   type(dmp_input), intent(in) :: dmp_inp
   !> Reuse existing solver instance
   logical, intent(in), optional :: reuse

   new_purification_solver_context%dmp_input = dmp_inp
   if (present(reuse)) new_purification_solver_context%reuse = reuse
end function new_purification_solver_context

!> Create new electronic solver
subroutine new(self, overlap, nel, kt)
   !> Instance of the solver factory
   class(purification_solver_context), intent(inout) :: self
   !> Overlap matrix
   real(wp), intent(in) :: overlap(:, :)
   !> Number of electrons per spin channel
   real(wp), intent(in) :: nel(:)
   !> Electronic temperature
   real(wp), intent(in) :: kt
   type(purification_solver), allocatable :: tmp

   if (self%ndim /= size(overlap, 1)) then
      self%ndim = size(overlap, 1)
      if (allocated(self%solver)) then
         call self%solver%delete(self%dmp_ptr)
      end if
      allocate(tmp)
      call new_purification(tmp, overlap, nel,  kt,  self%dmp_input, self%dmp_ptr, self%sgvd_ptr)
      call move_alloc(tmp, self%solver)
   end if

end subroutine new

!> Delete electronic solver instance
subroutine delete(self)
   !> Instance of the solver factory
   class(purification_solver_context), intent(inout) :: self

   if (self%reuse) then
      if (allocated(self%solver)) call self%solver%delete()
   else
      if (c_associated(self%dmp_ptr)) then
         call self%solver%delete(self%dmp_ptr)
         self%dmp_ptr = c_null_ptr
         self%sgvd_ptr = c_null_ptr
         self%ndim = 0
      end if
      if (allocated(self%solver)) deallocate(self%solver) 
   end if
      
   
end subroutine delete
    

end module