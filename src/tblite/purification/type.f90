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
      logical :: reuse = .false.
      class(solver_type), allocatable :: cached
   contains
      procedure :: new
      procedure :: delete
   end type purification_solver_context

   interface purification_solver_context
      module procedure new_purification_solver_context
   end interface

contains

type(purification_solver_context) function new_purification_solver_context(dmp_inp, reuse)
   type(dmp_input), intent(in) :: dmp_inp
   logical, intent(in), optional :: reuse

   new_purification_solver_context%dmp_input = dmp_inp
   if (present(reuse)) new_purification_solver_context%reuse = reuse
end function

subroutine new(self, solver, overlap, nel, kt, comm)
   class(purification_solver_context), intent(inout) :: self
   class(solver_type), allocatable, intent(out) :: solver
   real(wp), intent(in) :: overlap(:, :), nel(:), kt
   integer, intent(in), optional :: comm
   logical :: compatible

   if (allocated(self%cached)) then
      select type (cached => self%cached)
      type is (purification_solver)
         compatible = cached%ndim == size(overlap, 1) .and. cached%input%type == self%dmp_input%type &
            & .and. cached%input%runmode == self%dmp_input%runmode .and. cached%input%precision == self%dmp_input%precision
         if (compatible) call move_alloc(self%cached, solver)
      end select
      if (allocated(self%cached)) deallocate(self%cached)
   end if
   if (.not.allocated(solver)) allocate(purification_solver :: solver)
   select type (solver)
   type is (purification_solver)
      call new_purification(solver, overlap, nel, kt, self%dmp_input)
   end select
end subroutine

subroutine delete(self, solver)
   class(purification_solver_context), intent(inout) :: self
   class(solver_type), allocatable, intent(inout) :: solver

   if (self%reuse) then
      call move_alloc(solver, self%cached)
   else if (allocated(solver)) then
      deallocate(solver)
   end if
end subroutine

end module tblite_purification_solver_context
