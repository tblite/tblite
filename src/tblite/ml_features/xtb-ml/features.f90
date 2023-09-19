module tblite_xtbml_feature_type
    use mctc_env, only : wp
   use tblite_wavefunction_type, only : wavefunction_type
   use mctc_io, only : structure_type
   use tblite_integral_type, only : integral_type
   use tblite_basis_type, only : basis_type
   use tblite_container, only : container_cache
   use tblite_context , only : context_type
   use tblite_xtbml_convolution, only : xtbml_convolution_type
   use tblite_double_dictionary, only : double_dictionary_type
  use tblite_container, only : container_list
   implicit none
   private

type, public, abstract :: xtbml_feature_type
    integer :: n_features
    type(double_dictionary_type), allocatable :: dict, dict_ext
contains
    procedure(compute_features), deferred :: compute_features
    procedure(compute_extended), deferred :: compute_extended
    procedure :: get_n_features
end type xtbml_feature_type

abstract interface
    subroutine compute_features(self, mol, wfn, integrals, bas, contain_list, cache_list, prlevel, ctx)
        import :: wp, wavefunction_type, structure_type, integral_type, basis_type,&
        & container_cache, context_type, xtbml_feature_type, container_list
        class(xtbml_feature_type), intent(inout) :: self
        !> Molecular structure data
        type(structure_type), intent(in) :: mol
        !> Wavefunction strcuture data
        type(wavefunction_type), intent(in) :: wfn
        !> Integral container
        type(integral_type) :: integrals
        !> Single-point calculator
        type(basis_type), intent(in) :: bas
        !> List of containers 
        type(container_list), intent(inout) :: contain_list
        !> Container
        type(container_cache), intent(inout) :: cache_list(:)
        !> Context type
        type(context_type),intent(inout) :: ctx
        !> Print Level
        integer, intent(in) :: prlevel
    end subroutine
    subroutine compute_extended(self, mol, wfn, integrals, bas, contain_list, cache_list, prlevel, ctx, convolution)
        import :: wp, wavefunction_type, structure_type, integral_type, basis_type,&
        & container_cache, context_type, xtbml_feature_type, xtbml_convolution_type, container_list
        class(xtbml_feature_type), intent(inout) :: self
        !> Molecular structure data
        type(structure_type), intent(in) :: mol
        !> Wavefunction strcuture data
        type(wavefunction_type), intent(in) :: wfn
        !> Integral container
        type(integral_type) :: integrals
        !> Single-point calculator
        type(basis_type), intent(in) :: bas
        !> List of containers 
        type(container_list), intent(inout) :: contain_list
        !> Container
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

function get_n_features(self) result(n)
    class(xtbml_feature_type) :: self
    integer :: n 

    n = self%dict%get_n_entries()
    n = n + self%dict_ext%get_n_entries()
end function


end module  