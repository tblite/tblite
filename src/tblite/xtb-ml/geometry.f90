module tblite_xtbml_geometry_based
    use mctc_env, only : wp
    use tblite_xtbml_feature_type, only : xtbml_feature_type
    use mctc_env, only : wp
   use tblite_wavefunction, only : wavefunction_type
   use mctc_io, only : structure_type
   use tblite_integral_type, only : integral_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_container, only : container_cache
   use tblite_context , only : context_type
   use tblite_xtbml_convolution, only : xtbml_convolution_type
   use tblite_double_dictionary, only : double_dictionary_type
    implicit none
    private

    integer, parameter ::features = 1
    integer, parameter :: ext_features = 1
    

type, public, extends(xtbml_feature_type) :: xtbml_geometry_features_type
    character(len=*) :: label = "geometry-based features"
    real(wp), allocatable ::  cn_atom(:)
    real(wp), allocatable ::  delta_cn(:, :)
    integer :: n_features
    type(enum_geometry_features) :: labels = enum_geometry_features()
    type(enum_ext_geometry_features) :: ext_labels = enum_geometry_features()
    type(double_dictionary_type) :: dict

contains
    procedure :: compute_features
    procedure :: compute_extended
    procedure :: get_n_features
    procedure :: get_feature_labels
end type

contains

subroutine compute_features(self, mol, wfn, integrals, calc, cache, prlevel, ctx)
    use tblite_ncoord_exp, only : new_exp_ncoord, exp_ncoord_type
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
    self%n_features = self%n_features + features

    allocate(self%cn_atom(mol%nat))
    
    call new_exp_ncoord(ncoord_exp, mol)
    call ncoord_exp%get_cn(mol, self%cn_atom) 

    call self%dict%add_entry("CN", self%cn_atom)
    
end subroutine

subroutine compute_extended(self, mol, wfn, integrals, calc, cache, prlevel, ctx, convolution)
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
    
    allocate(self%delta_cn(mol%nat, convolution%n_a))
    
    call get_delta_cn(mol%nat, n, self%cn_atom, mol%id, mol%xyz, self%delta_cn, convolution)
    self%n_features = self%n_features + ext_features

    call self%dict%add_entry("delta_CN", self%delta_cn) 
end subroutine

subroutine get_delta_cn(nat, n_a, cn, at, xyz, delta_cn, conv)
    use tblite_timer, only : timer_type, format_time
    integer, intent(in) :: nat, at(nat), n_a
    real(wp), intent(in) :: cn(nat), xyz(3, nat)
    real(wp), intent(out) :: delta_cn(nat, n_a)
    real(wp):: delta_cn_tmp(nat, n_a), stime
    integer :: i, j, k
    !> Convolution container
    type(xtbml_convolution_type) :: conv

    
    delta_cn = 0.0_wp
    !$omp parallel do default(none) collapse(2)&
    !$omp shared(nat, conv%inv_cn_a, n_a)&
    !$omp private(delta_cn,i,j,k)
    do k = 1, n_a
       do i = 1, nat
          do j = 1, nat
             if (i == j) cycle
             !$omp atomic
             delta_cn(i, k) = delta_cn(i, k) + cn(j)/conv%inv_cn_a(i, j, k)
 
          end do
       end do
    end do
    !$omp end parallel do
 end subroutine





end module tblite_xtbml_geometry_based   