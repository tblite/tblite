module tblite_xtbml_features
    use mctc_env, only : wp
    use tblite_ml_features_type, only : ml_features_type
    use tblite_xtbml_convolution, only : xtbml_convolution_type
    use tblite_xtbml_geometry_based, only : xtbml_geometry_features_type
    use tblite_xtbml_density_based, only : xtbml_density_features_type
    use tblite_wavefunction_type, only : wavefunction_type
    use mctc_io, only : structure_type
    use tblite_basis_type, only : basis_type
    use tblite_results, only : results_type
    use tblite_integral_type, only : integral_type
    use tblite_basis_type, only : basis_type
    use tblite_container, only : container_cache
    use tblite_results, only : results_type
    use tblite_context, only : context_type
    use tblite_double_dictionary, only : double_dictionary_type
    implicit none
    private
    public :: xtbml_type, new_xtbml_features
    type, extends(ml_features_type) :: xtbml_type
        type(xtbml_geometry_features_type), allocatable :: geom
        type(xtbml_density_features_type), allocatable :: dens
        type(xtbml_convolution_type), allocatable :: conv
    contains
        procedure :: compute
    end type

contains

    subroutine new_xtbml_features(param, new_xtbml_model)
        use tblite_param_ml_features, only : ml_features_record
        type(ml_features_record) :: param
        logical :: xtbml_geometry, xtbml_density, xtbml_convolution, xtbml_tensor
        
        type(xtbml_type), intent(out) :: new_xtbml_model

        if (param%xtbml_geometry) allocate(new_xtbml_model%geom)
        if (param%xtbml_density) then
            allocate(new_xtbml_model%dens)
            new_xtbml_model%dens%return_xyz = param%xtbml_tensor
        end if
        if (param%xtbml_convolution) then
            allocate(new_xtbml_model%conv)
            new_xtbml_model%conv%a = param%xtbml_a
            new_xtbml_model%conv%n_a = size(param%xtbml_a)
        end if
        !if (param%xtbml_orbital_energy)

    end subroutine

    subroutine compute(self, mol, wfn, integrals, bas, ccache, dcache, rcache, ctx, prlevel, dict)
        class(xtbml_type),intent(in) :: self
        !> Molecular structure data
        type(structure_type), intent(in) :: mol
        !> Wavefunction strcuture data
        type(wavefunction_type), intent(in) :: wfn
        !> integral container for dipole and quadrupole integrals for CAMMs
        type(integral_type) :: integrals
        !> Single-point calculator conatiner
        type(basis_type), intent(in) :: bas
        !> Context container for writing to stdout
        type(context_type), intent(inout) :: ctx
        !> Compute cache containers
        type(container_cache),intent(inout) :: ccache, dcache, rcache
        type(double_dictionary_type), intent(inout) :: dict
        type(container_cache) :: void_cache
        integer :: prlevel
        type(xtbml_type) :: ml_model

        ml_model = self

        if (allocated(ml_model%geom)) then
            associate(category => ml_model%geom)
                call category%compute_features(mol, wfn, integrals, bas, void_cache, prlevel, ctx)
                dict = dict + category%dict
            end associate
        end if

        if (allocated(ml_model%dens)) then
            associate(category => ml_model%dens)
                call category%compute_features(mol, wfn, integrals, bas, void_cache, prlevel, ctx)
                dict = dict + category%dict
            end associate
        end if

        if (allocated(ml_model%conv)) then 
            if (allocated(ml_model%geom)) then
                associate(category => ml_model%geom)
                    call category%compute_extended(mol, wfn, integrals, bas, void_cache, prlevel, ctx, ml_model%conv)
                    dict = dict + category%dict_ext
                end associate
            end if
    
            if (allocated(ml_model%dens)) then
                associate(category => ml_model%dens)
                    call category%compute_extended(mol, wfn, integrals, bas, void_cache, prlevel, ctx, ml_model%conv)
                    dict = dict + category%dict_ext
                end associate
            end if
        end if            

    end subroutine
    
end module tblite_xtbml_features