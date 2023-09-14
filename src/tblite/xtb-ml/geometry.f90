module tblite_xtbml_geometry_based
    use mctc_env, only : wp
    use tblite_xtbml_feature_type, only : xtbml_feature_type
    use mctc_env, only : wp
   use tblite_wavefunction_type, only : wavefunction_type
   use mctc_io, only : structure_type
   use tblite_integral_type, only : integral_type
   use tblite_basis_type, only : basis_type
   use tblite_container, only : container_cache
   use tblite_context , only : context_type
   use tblite_xtbml_convolution, only : xtbml_convolution_type
    implicit none
    private

    integer, parameter ::features = 1
    integer, parameter :: ext_features = 1
    

type, public, extends(xtbml_feature_type) :: xtbml_geometry_features_type
    character(len=23) :: label = "geometry-based features"
    real(wp), allocatable ::  cn_atom(:)
    real(wp), allocatable ::  delta_cn(:, :)

contains
    procedure :: compute_features
    procedure :: compute_extended
end type

contains



subroutine compute_features(self, mol, wfn, integrals, bas, cache, prlevel, ctx)
    use tblite_ncoord_exp, only : new_exp_ncoord, exp_ncoord_type
    class(xtbml_geometry_features_type), intent(inout) :: self
    !> Molecular structure data
    type(structure_type), intent(in) :: mol
    !> Wavefunction strcuture data
    type(wavefunction_type), intent(in) :: wfn
    !> Integral container
    type(integral_type) :: integrals
    !> Single-point calculator
    type(basis_type), intent(in) :: bas
    !> Container
    type(container_cache), intent(inout) :: cache
    !> Context type
    type(context_type),intent(inout) :: ctx
    !> Print Level
    integer, intent(in) :: prlevel
    type(exp_ncoord_type) :: ncoord_exp
    self%n_features = self%n_features + features
    allocate(self%dict)
    allocate(self%cn_atom(mol%nat))
    
    call new_exp_ncoord(ncoord_exp, mol)
    call ncoord_exp%get_cn(mol, self%cn_atom) 

    call self%dict%add_entry("CN", self%cn_atom)
    
end subroutine

subroutine compute_extended(self, mol, wfn, integrals, bas, cache, prlevel, ctx, convolution)
    class(xtbml_geometry_features_type), intent(inout) :: self
    !> Molecular structure data
    type(structure_type), intent(in) :: mol
    !> Wavefunction strcuture data
    type(wavefunction_type), intent(in) :: wfn
    !> Integral container
    type(integral_type) :: integrals
    !> Single-point calculator
    type(basis_type), intent(in) :: bas
    !> Container
    type(container_cache), intent(inout) :: cache
    !> Context type
    type(context_type),intent(inout) :: ctx
    !> Print Level
    integer, intent(in) :: prlevel
    !> Convolution container
    type(xtbml_convolution_type) :: convolution
    allocate(self%dict_ext)
    allocate(self%delta_cn(mol%nat, convolution%n_a))
    convolution%cn = self%cn_atom 
    call get_delta_cn(mol%nat, self%cn_atom, mol%id, mol%xyz, self%delta_cn, convolution)
    self%n_features = self%n_features + ext_features

    call self%dict_ext%add_entry("delta_CN", self%delta_cn) 
end subroutine

subroutine get_delta_cn(nat, cn, at, xyz, delta_cn, conv)
    use tblite_timer, only : timer_type, format_time
    integer, intent(in) :: nat, at(nat)
    type(xtbml_convolution_type) :: conv
    real(wp), intent(in) :: cn(nat), xyz(3, nat)
    real(wp), intent(out) :: delta_cn(nat, conv%n_a)
    real(wp):: delta_cn_tmp(nat, conv%n_a), stime
    integer :: i, j, k
    !> Convolution container
    

    delta_cn = 0.0_wp
    !$omp parallel do default(none) collapse(2)&
    !$omp shared(nat, conv, cn)&
    !$omp private(delta_cn,i,j,k)
    do k = 1, conv%n_a
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