module tblite_xtbml_feature_type
    use mctc_env, only : wp
   use tblite_wavefunction, only : wavefunction_type
   use mctc_io, only : structure_type
   use tblite_integral_type, only : integral_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_container, only : container_cache
   use tblite_context , only : context_type
   use tblite_xtbml_convolution, only : xtbml_convolution_type
   implicit none
   private

type, public, abstract :: xtbml_feature_type

contains
    procedure(compute_features), deferred :: compute_features
    procedure(compute_extended), deferred :: compute_extended
    procedure(get_n_features), deferred :: get_n_features
end type xtbml_feature_type

abstract interface
    subroutine compute_features(self, mol, wfn, integrals, calc, cache, prlevel, ctx)
        import :: wp, wavefunction_type, structure_type, wavefunction_type, integral_type, xtb_calculator,&
        & container_cache, context_type, xtbml_feature_type
        class(xtbml_feature_type), intent(inout) :: self
        !> Molecular structure data
        type(structure_type), intent(in) :: mol
        !> Wavefunction strcuture data
        type(wavefunction_type), intent(in) :: wfn
        !> Integral container
        type(integral_type) :: integrals
        !> Single-point calculator
        type(xtb_calculator), intent(in) :: calc
        !> Container
        type(container_cache), intent(inout) :: cache
        !> Context type
        type(context_type),intent(inout) :: ctx
        !> Print Level
        integer, intent(in) :: prlevel
    end subroutine
    subroutine compute_extended(self, mol, wfn, integrals, calc, cache, prlevel, ctx, convolution)
        import :: wp, wavefunction_type, structure_type, wavefunction_type, integral_type, xtb_calculator,&
        & container_cache, context_type, xtbml_feature_type, xtbml_convolution_type
        class(xtbml_feature_type), intent(inout) :: self
        !> Molecular structure data
        type(structure_type), intent(in) :: mol
        !> Wavefunction strcuture data
        type(wavefunction_type), intent(in) :: wfn
        !> Integral container
        type(integral_type) :: integrals
        !> Single-point calculator
        type(xtb_calculator), intent(in) :: calc
        !> Container
        type(container_cache), intent(inout) :: cache
        !> Context type
        type(context_type),intent(inout) :: ctx
        !> Print Level
        integer, intent(in) :: prlevel
        !> Convolution container
        type(xtbml_convolution_type) :: convolution 
    end subroutine
    subroutine get_n_features(self,n)
        import :: xtbml_feature_type
        class(xtbml_feature_type), intent(inout) :: self
        !> Number of features in feature type
        integer :: n
    end subroutine
end interface


end module  