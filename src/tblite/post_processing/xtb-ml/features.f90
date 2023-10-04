module tblite_xtbml_feature_type
    use mctc_env, only : wp
   use tblite_wavefunction_type, only : wavefunction_type
   use mctc_io, only : structure_type
   use tblite_integral_type, only : integral_type
   use tblite_container, only : container_cache
   use tblite_context , only : context_type
   use tblite_xtbml_convolution, only : xtbml_convolution_type
   use tblite_double_dictionary, only : double_dictionary_type
  use tblite_container, only : container_list
    use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

type, public, abstract :: xtbml_feature_type
    integer :: n_features
    character(len=:), allocatable :: label
    type(double_dictionary_type), allocatable :: dict, dict_ext
contains
    procedure(compute_features), deferred :: compute_features
    procedure(compute_extended), deferred :: compute_extended
    procedure :: get_n_features
    procedure :: info
    procedure :: setup
end type xtbml_feature_type

character(len=*), parameter :: label = "General feature class"

abstract interface
    subroutine compute_features(self, mol, wfn, integrals, calc, cache_list, prlevel, ctx)
        import :: wp, wavefunction_type, structure_type, integral_type, xtb_calculator,&
        & container_cache, context_type, xtbml_feature_type, container_list
        class(xtbml_feature_type), intent(inout) :: self
        !> Molecular structure data
        type(structure_type), intent(in) :: mol
        !> Wavefunction strcuture data
        type(wavefunction_type), intent(in) :: wfn
        !> Integral container
        type(integral_type) :: integrals
        !> Single-point calculator
        type(xtb_calculator), intent(in) :: calc
        !> List of containers 
        type(container_cache), intent(inout) :: cache_list(:)
        !> Context
        type(context_type),intent(inout) :: ctx
        !> Print Level
        integer, intent(in) :: prlevel
    end subroutine
    subroutine compute_extended(self, mol, wfn, integrals, calc, cache_list, prlevel, ctx, convolution)
        import :: wp, wavefunction_type, structure_type, integral_type, xtb_calculator,&
        & container_cache, context_type, xtbml_feature_type, xtbml_convolution_type, container_list
        class(xtbml_feature_type), intent(inout) :: self
        !> Molecular structure data
        type(structure_type), intent(in) :: mol
        !> Wavefunction strcuture data
        type(wavefunction_type), intent(in) :: wfn
        !> Integral container
        type(integral_type) :: integrals
        !> Single-point calculator
        type(xtb_calculator), intent(in) :: calc
        !> List of containers 
        type(container_cache), intent(inout) :: cache_list(:)
        !> Context type
        type(context_type),intent(inout) :: ctx
        !> Print Level
        integer, intent(in) :: prlevel
        !> Convolution container
        type(xtbml_convolution_type) :: convolution 
    end subroutine

end interface

contains

subroutine setup(self)
    class(xtbml_feature_type) :: self
    self%label = label
end subroutine

function get_n_features(self) result(n)
    class(xtbml_feature_type) :: self
    integer :: n 

    n = self%dict%get_n_entries()
    n = n + self%dict_ext%get_n_entries()
end function

pure function info(self, verbosity, indent) result(str)
   !> Instance of the interaction container
   class(xtbml_feature_type), intent(in) :: self
   !> Verbosity level
   integer, intent(in) :: verbosity
   !> Indentation level
   character(len=*), intent(in) :: indent
   !> Information on the container
   character(len=:), allocatable :: str

   if (allocated(self%label)) then
      str = indent // self%label 
   else
      str = "Unknown"
   end if
end function info


end module  