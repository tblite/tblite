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
    use tblite_xtbml_orbital_energy, only : xtbml_orbital_features_type
    use tblite_container, only : container_list
    use tblite_xtbml_energy_features, only : xtbml_energy_features_type
    implicit none
    private
    public :: xtbml_type, new_xtbml_features
    type, extends(ml_features_type) :: xtbml_type
        type(xtbml_geometry_features_type), allocatable :: geom
        type(xtbml_density_features_type), allocatable :: dens
        type(xtbml_orbital_features_type), allocatable :: orb
        type(xtbml_energy_features_type), allocatable :: energy
        type(xtbml_convolution_type), allocatable :: conv
    contains
        procedure :: compute
        procedure :: pack_res
        procedure :: info
    end type
    character(len=*), parameter :: label = "xtbml features"
contains

subroutine pack_res(self, mol, dict, res)
    class(xtbml_type),intent(in) :: self
    type(structure_type), intent(in) :: mol
    type(results_type), intent(inout) :: res
    type(double_dictionary_type) :: dict 
    integer :: i, n
    real(wp), allocatable :: tmp_array(:)
    character(len=:), allocatable :: tmp_label
    !the partitioning weights are also included allthough they are stored seperately 
    n = dict%get_n_entries()
    res%n_features = n

    allocate(res%ml_features(mol%nat, n))
    allocate(res%ml_labels(res%n_features))

    do i = 1, n
        call dict%get_label(i, tmp_label)
        res%ml_labels(i) = tmp_label
        
        call dict%get_entry(i, tmp_array)
        if (trim(tmp_label) == "w_tot") then
            allocate(res%w_xtbml(mol%nat))
            res%w_xtbml = tmp_array
            cycle
        end if
        res%ml_features(:, i) = tmp_array 
    end do
    
end subroutine

    subroutine new_xtbml_features(param, new_xtbml_model)
        use tblite_param_ml_features, only : ml_features_record
        type(ml_features_record) :: param 
        type(xtbml_type), intent(inout) :: new_xtbml_model
        
        new_xtbml_model%label = label

        if (param%xtbml_geometry) then
            allocate(new_xtbml_model%geom)
            call new_xtbml_model%geom%setup()
        end if

        if (param%xtbml_density) then
            allocate(new_xtbml_model%dens)
            call new_xtbml_model%dens%setup()
            new_xtbml_model%dens%return_xyz = param%xtbml_tensor
        end if
        if (param%xtbml_orbital_energy) then 
            allocate(new_xtbml_model%orb)
            call new_xtbml_model%orb%setup()
        end if
        if (param%xtbml_energy) then 
            allocate(new_xtbml_model%energy)
            call new_xtbml_model%energy%setup()
        end if
        if (param%xtbml_convolution) then
            allocate(new_xtbml_model%conv)
            new_xtbml_model%conv%a = param%xtbml_a
            new_xtbml_model%conv%n_a = size(param%xtbml_a)
            call new_xtbml_model%conv%setup()
        end if

    end subroutine

    subroutine compute(self, mol, wfn, integrals, bas, contain_list, ctx, prlevel, dict)
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
        type(container_list), intent(inout) :: contain_list
        type(double_dictionary_type), intent(inout) :: dict
        integer :: prlevel
        type(xtbml_type) :: ml_model

        ml_model = self

        if (allocated(ml_model%geom)) then
            associate(category => ml_model%geom)
                call category%compute_features(mol, wfn, integrals, bas, contain_list, prlevel, ctx)
                dict = dict + category%dict
            end associate
        end if

        if (allocated(ml_model%dens)) then
            associate(category => ml_model%dens)
                call category%compute_features(mol, wfn, integrals, bas, contain_list, prlevel, ctx)
                dict = dict + category%dict
            end associate
        end if

        if (allocated(ml_model%orb)) then
            associate(category => ml_model%orb)
                call category%compute_features(mol, wfn, integrals, bas, contain_list, prlevel, ctx)
                dict = dict + category%dict
            end associate
        end if

        if (allocated(ml_model%energy)) then
            associate(category => ml_model%energy)
                call category%compute_features(mol, wfn, integrals, bas, contain_list, prlevel, ctx)
                dict = dict + category%dict
            end associate
        end if

        if (allocated(ml_model%conv)) then 

            if (allocated(ml_model%geom)) then
                ml_model%conv%cn = ml_model%geom%cn_atom
                call ml_model%conv%compute_kernel(mol)
                associate(category => ml_model%geom)
                    call category%compute_extended(mol, wfn, integrals, bas, contain_list, prlevel, ctx, ml_model%conv)
                    dict = dict + category%dict_ext
                end associate
            else
                call ml_model%conv%compute_kernel(mol) 
            end if
    
            if (allocated(ml_model%dens)) then
                associate(category => ml_model%dens)
                    call category%compute_extended(mol, wfn, integrals, bas, contain_list, prlevel, ctx, ml_model%conv)
                    dict = dict + category%dict_ext
                end associate
            end if

            if (allocated(ml_model%orb)) then
                associate(category => ml_model%orb)
                    call category%compute_extended(mol, wfn, integrals, bas, contain_list, prlevel, ctx, ml_model%conv)
                    dict = dict + category%dict_ext
                end associate
            end if
        end if            

    end subroutine

pure function info(self, verbosity, indent) result(str)
    !> Instance of the interaction container
    class(xtbml_type), intent(in) :: self
    !> Verbosity level
    integer, intent(in) :: verbosity
    !> Indentation level
    character(len=*), intent(in) :: indent
    !> Information on the container
    character(len=:), allocatable :: str
    character(len=*), parameter :: nl = new_line('a')
    character(len=:), allocatable :: category_indent 
    if (allocated(self%label)) then
       str = self%label
    else
       str = "Unknown"
    end if

    category_indent = indent // " * "
    
    if (allocated(self%geom)) then 
        str = str // nl // self%geom%info(verbosity, category_indent)
    end if 

    if (allocated(self%dens)) then 
        str = str // nl // self%dens%info(verbosity, category_indent)
    end if
    
    if (allocated(self%orb)) then 
        str = str // nl // self%orb%info(verbosity, category_indent)
    end if

    if (allocated(self%energy)) then 
        str = str // nl // self%energy%info(verbosity, category_indent)
    end if

    if (allocated(self%conv)) then
        str = str // nl // self%conv%info(verbosity, category_indent)
    end if
 end function info
    
end module tblite_xtbml_features