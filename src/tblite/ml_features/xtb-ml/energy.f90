module tblite_xtbml_energy_features
    use mctc_env, only : wp
   use mctc_io_convert, only : autoev
   use tblite_xtbml_feature_type, only : xtbml_feature_type
   use tblite_wavefunction_type, only : wavefunction_type
   use mctc_io, only : structure_type
   use tblite_integral_type, only : integral_type
   use tblite_basis_type, only : basis_type 
   use tblite_container, only : container_cache
   use tblite_context , only : context_type
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_xtbml_convolution, only : xtbml_convolution_type
   use tblite_xtbml_atomic_frontier, only : atomic_frontier_orbitals
  use tblite_container, only : container_list, container_type
   implicit none
   private

   type, public, extends(xtbml_feature_type) :: xtbml_energy_features_type
      character(len=21) :: label = "energy-based features"
    contains
      procedure :: compute_features
      procedure :: compute_extended
    end type

contains

subroutine compute_features(self, mol, wfn, integrals, bas, contain_list, cache_list, prlevel, ctx)
    use tblite_repulsion, only : tb_repulsion
    use tblite_xtb_coulomb, only : tb_coulomb
  use tblite_scf_iterator, only : get_electronic_energy, reduce
  use tblite_disp_d3, only : d3_dispersion, new_d3_dispersion
  use tblite_classical_halogen, only : halogen_correction
  use tblite_disp_d4, only : d4_dispersion, new_d4_dispersion
  class(xtbml_energy_features_type), intent(inout) :: self
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
  type(d3_dispersion), allocatable :: d3
  type(d4_dispersion), allocatable :: d4
  class(container_type), allocatable :: cont
  real(wp), allocatable :: tmp_energy(:), e_ao(:), e_disp_tot(:), e_disp_ATM(:)
  integer :: i

  allocate(self%dict)
  
  allocate(e_ao(bas%nao), source=0.0_wp)
  allocate(tmp_energy(mol%nat), source=0.0_wp)
  call get_electronic_energy(integrals%hamiltonian, wfn%density, e_ao)
  call reduce(tmp_energy, e_ao, bas%ao2at)

  call self%dict%add_entry("E_EHT", tmp_energy)

  do i = 1, size(cache_list)
    tmp_energy = 0.0_wp
    call contain_list%pop(cont, 1)
    select type(cont)
    type is (tb_repulsion)
        call cont%update(mol, cache_list(i))
        call cont%get_engrad(mol, cache_list(i), tmp_energy)
        call self%dict%add_entry("E_rep", tmp_energy)
    type is (tb_coulomb)
        call cont%update(mol, cache_list(i))
        if (allocated(cont%es2)) then
            call cont%es2%update(mol, cache_list(i))
            call cont%es2%get_energy(mol, cache_list(i), wfn, tmp_energy)
        end if
        if (allocated(cont%es3)) then
            call cont%es3%update(mol, cache_list(i))
            call cont%es3%get_energy(mol, cache_list(i), wfn, tmp_energy)
        end if
        call self%dict%add_entry("E_ies_ixc", tmp_energy)
        if (allocated(cont%aes2)) then
            tmp_energy = 0.0_wp
            call cont%aes2%get_AXC(mol, wfn, tmp_energy)
            call self%dict%add_entry("E_AXC", tmp_energy)
            tmp_energy = 0.0_wp 
            call cont%aes2%get_energy_aes_xtb(mol, cache_list(i), wfn, tmp_energy)
            call self%dict%add_entry("E_AES", tmp_energy) 
        end if
    type is (halogen_correction)
        call cont%update(mol, cache_list(i))
        call cont%get_engrad(mol, cache_list(i), tmp_energy)
        call self%dict%add_entry("E_HX", tmp_energy)
    type is (d3_dispersion)
        allocate(e_disp_tot(mol%nat), e_disp_ATM(mol%nat), source=0.0_wp)
        call cont%update(mol, cache_list(i))
        call cont%get_engrad(mol, cache_list(i), e_disp_tot)
        write(*,*) cont%param%s9
        allocate(d3)
        call new_d3_dispersion(d3, mol, s6=0.0_wp, s8=0.0_wp, a1=cont%param%a1, a2=cont%param%a2, s9=cont%param%s9)
        call d3%update(mol, cache_list(i))
        call d3%get_engrad(mol, cache_list(i), e_disp_ATM)
        call self%dict%add_entry("E_disp2", e_disp_tot-e_disp_ATM)
        call self%dict%add_entry("E_disp3", e_disp_ATM) 
    type is (d4_dispersion)
        allocate(e_disp_tot(mol%nat), e_disp_ATM(mol%nat), source=0.0_wp)
        call cont%update(mol, cache_list(i))
        call cont%get_engrad(mol, cache_list(i), e_disp_tot)
        call cont%get_energy(mol, cache_list(i), wfn, e_disp_tot)
            
        allocate(d4)
        call new_d4_dispersion(d4, mol, s6=0.0_wp, s8=0.0_wp, a1=cont%param%a1, a2=cont%param%a2, s9=cont%param%s9)
        call d4%update(mol, cache_list(i))
        call d4%get_engrad(mol, cache_list(i), e_disp_ATM)
        call self%dict%add_entry("E_disp2", e_disp_tot-e_disp_ATM)
        call self%dict%add_entry("E_disp3", e_disp_ATM)
    class default
        call cont%update(mol, cache_list(i))
        call cont%get_engrad(mol, cache_list(i), tmp_energy)
        call cont%get_energy(mol, cache_list(i), wfn, tmp_energy)
        call self%dict%add_entry(cont%info(0, ""), tmp_energy)
    end select
  end do

end subroutine

subroutine compute_extended(self, mol, wfn, integrals, bas, contain_list, cache_list, prlevel, ctx, convolution)
  use tblite_output_format, only : format_string 
  class(xtbml_energy_features_type), intent(inout) :: self
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
end module