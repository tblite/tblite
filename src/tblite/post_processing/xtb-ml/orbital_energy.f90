module tblite_xtbml_orbital_energy
   use mctc_env, only : wp
   use mctc_io_convert, only : autoev
   use tblite_xtbml_feature_type, only : xtbml_feature_type
   use tblite_wavefunction_type, only : wavefunction_type
   use mctc_io, only : structure_type
   use tblite_integral_type, only : integral_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_container, only : container_cache
   use tblite_context , only : context_type
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_xtbml_convolution, only : xtbml_convolution_type
   use tblite_xtbml_atomic_frontier, only : atomic_frontier_orbitals
  use tblite_container, only : container_list
   implicit none
   private
  character(len=*), parameter :: label = "orbital energy-based features"
   type, public, extends(xtbml_feature_type) :: xtbml_orbital_features_type
      
      real(wp),allocatable ::  response(:)
      real(wp),allocatable ::  egap(:)
      real(wp),allocatable ::  chempot(:)
      real(wp),allocatable ::  ehoao_a(:)
      real(wp),allocatable ::  eluao_a(:)
      real(wp),allocatable ::  ehoao_b(:)
      real(wp),allocatable ::  eluao_b(:)
      real(wp),allocatable ::  delta_chempot(:,:)
      real(wp),allocatable ::  delta_egap(:,:)
      real(wp),allocatable ::  delta_eluao(:,:)
      real(wp),allocatable ::  delta_ehoao(:,:)
    contains
      procedure :: compute_features
      procedure :: compute_extended
      procedure, private :: allocate
      procedure, private :: allocate_extended
      procedure :: setup 
    end type

contains

subroutine setup(self)
  class(xtbml_orbital_features_type) :: self
  self%label = label
  if (allocated(self%dict)) deallocate(self%dict)
  allocate(self%dict)
  if (allocated(self%dict_ext)) deallocate(self%dict_ext)
  allocate(self%dict_ext)
end subroutine

subroutine allocate(self, nat)
  class(xtbml_orbital_features_type) :: self
  integer :: nat

  allocate(self%response(nat), source=0.0_wp)
  allocate(self%egap(nat), source=0.0_wp)
  allocate(self%chempot(nat), source=0.0_wp)
  allocate(self%ehoao_a(nat), source=0.0_wp)
  allocate(self%eluao_a(nat), source=0.0_wp)
  allocate(self%ehoao_b(nat), source=0.0_wp)
  allocate(self%eluao_b(nat), source=0.0_wp)

end subroutine

subroutine allocate_extended(self, nat, n_a)
  class(xtbml_orbital_features_type) :: self
  integer :: nat, n_a
  
  allocate(self%delta_chempot(nat, n_a), source=0.0_wp)
  allocate(self%delta_egap(nat, n_a), source=0.0_wp)
  allocate(self%delta_eluao(nat, n_a), source=0.0_wp)
  allocate(self%delta_ehoao(nat, n_a), source=0.0_wp)

end subroutine

subroutine compute_features(self, mol, wfn, integrals, calc, cache_list, prlevel, ctx)
  class(xtbml_orbital_features_type), intent(inout) :: self
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
  logical :: print_afo = .false.
  real(wp) :: focc_(2,size(wfn%focc, dim=1))
  integer :: i, j
  real(wp) :: nel_
  
  call self%allocate(mol%nat)
  self%label = label

  focc_ = 0.0_wp

  if (prlevel > 2) print_afo = .true.

  if (size(wfn%nel(:)) > 1) then
    do j = 1,2
        nel_ = wfn%nel(j)
        do i = 1, size(wfn%focc)
          if (nel_ > 1.0_wp) then
              focc_(j,i) = 1.0_wp
              nel_ = nel_ - 1.0_wp
          else
              focc_(j,i) = nel_
              exit
          end if
        end do
    end do
  else
    focc_(1,:) = wfn%focc(:,1)/2.0_wp
    focc_(2,:) = wfn%focc(:,1)/2.0_wp
  end if

  call atomic_frontier_orbitals(mol%nat, calc%bas%nao, focc_(1,:), focc_(2,:), wfn%emo(:, 1)*autoev, &
    calc%bas%ao2at, wfn%coeff(:, :, 1), integrals%overlap(:, :), &
    self%response, self%egap, self%chempot, self%ehoao_a, &
    self%eluao_a, self%ehoao_b, self%eluao_b, print_afo, ctx)
  associate(dict => self%dict)
    call dict%add_entry("response", self%response)
    call dict%add_entry("gap", self%egap)
    call dict%add_entry("chem_pot", self%chempot)
    call dict%add_entry("HOAO_a", self%ehoao_a)
    call dict%add_entry("LUAO_a", self%eluao_a)
    call dict%add_entry("HOAO_b", self%ehoao_b)
    call dict%add_entry("LUAO_b", self%eluao_b)
  end associate
end subroutine

subroutine compute_extended(self, mol, wfn, integrals, calc, cache_list, prlevel, ctx, convolution)
  use tblite_output_format, only : format_string 
  class(xtbml_orbital_features_type), intent(inout) :: self
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

  integer :: n, i
  character(len=:), allocatable :: a_label
  real(wp), allocatable :: beta(:, :, :)
  real(wp) :: hl_gap
  
  call self%allocate_extended(mol%nat, convolution%n_a)
  

  n=convolution%n_a
  allocate(beta(mol%nat, mol%nat, n))
  call get_beta(mol%nat, n, mol%id, mol%xyz, beta, convolution%kernel)

  call get_chem_pot_ext(mol%nat, n, beta, self%chempot, self%delta_chempot)
  hl_gap=(wfn%emo(wfn%homo(1), 1)-wfn%emo(wfn%homo(1)+1, 1)) * autoev

  call get_e_gap_ext(mol%nat, n, hl_gap, beta, self%egap, self%delta_egap)

  call get_ehoao_ext(mol%nat, n, self%delta_chempot, self%delta_egap, self%delta_ehoao)

  call get_eluao_ext(mol%nat, n, self%delta_chempot, self%delta_egap, self%delta_eluao)
  
  associate( dict => self%dict_ext)
  do i = 1, n
    a_label = "_"//trim(adjustl(format_string(convolution%a(i), '(f12.2)')) )
    if (a_label .eq. "_1.00") a_label = ""
    call dict%add_entry("delta_gap"//a_label, self%delta_egap(:, i))
    call dict%add_entry("delta_chem_pot"//a_label, self%delta_chempot(:, i))
    call dict%add_entry("delta_HOAO"//a_label, self%delta_ehoao(:, i)) 
    call dict%add_entry("delta_LUAO"//a_label, self%delta_eluao(:, i))
  end do
  end associate
end subroutine

subroutine get_beta(nat, n_a, at, xyz, beta, kernel)
  intrinsic :: sum
  integer, intent(in) :: nat, at(nat), n_a
  real(wp), intent(in) :: xyz(3, nat), kernel(nat, nat, n_a)
  real(wp) :: beta(nat, nat, n_a)
  real(wp) :: sigma_tot, sigma(nat, nat, n_a)
  real(wp) :: damp_func

  integer :: a, b, k
  do k = 1, n_a
     do a = 1, nat
        do b = 1, nat
           damp_func = kernel(a, b, k)
           sigma(a, b, k) = 1/damp_func
        end do
     end do
  end do

  do k = 1, n_a
     do a = 1, nat
        do b = 1, nat
           beta(a, b, k) = sigma(a, b, k)/sum(sigma(a, :, k))
        end do
     end do
  end do
end subroutine

subroutine get_chem_pot_ext(nat, n_a, beta, chempot, chempot_ext)
  integer, intent(in) :: nat, n_a
  real(wp), intent(in) :: beta(nat, nat, n_a)
  real(wp), intent(in) :: chempot(nat)
  real(wp) :: chempot_ext(nat, n_a)
  integer :: a, b, k
  do k = 1, n_a
     do a = 1, nat
        do b = 1, nat
          chempot_ext(a, k) = chempot_ext(a, k) + beta(a, b, k) * chempot(b)
        end do
     end do
  end do
end subroutine

subroutine get_e_gap_ext(nat, n_a, hl_gap, beta, e_gap, e_gap_ext)
  integer, intent(in) :: nat, n_a
  real(wp), intent(in) :: hl_gap
  real(wp), intent(in) :: beta(nat, nat, n_a)
  real(wp), intent(in) :: e_gap(nat)
  real(wp) :: e_gap_ext(nat, n_a), e_gap_tot
  integer :: a, b, k

  e_gap_tot = 0.0_wp
  do k = 1, n_a
     do a = 1, nat
        e_gap_tot = e_gap_tot + e_gap(a)
        do b = 1, nat
           e_gap_ext(a, k) = e_gap_ext(a, k) + beta(a, b, k) * e_gap(b)
        end do
     end do
  end do
  ! correction factor for mol. properties
  !e_gap_ext = e_gap_ext * (hl_gap * nat/e_gap_tot)
end subroutine

subroutine get_ehoao_ext(nat, n_a, chempot_ext, e_gap_ext, ehoao_ext)
  integer, intent(in) :: nat, n_a
  real(wp), intent(in) :: chempot_ext(nat, n_a)
  real(wp), intent(in) :: e_gap_ext(nat, n_a)
  real(wp), intent(out) :: ehoao_ext(nat, n_a)
  integer :: a, k
  do k = 1, n_a
     do a = 1, nat
        ehoao_ext(a, k) = chempot_ext(a, k) - e_gap_ext(a, k) / 2
     end do
  end do
end subroutine

subroutine get_eluao_ext(nat, n_a, chempot_ext, e_gap_ext, eluao_ext)
  integer, intent(in) :: nat, n_a
  real(wp), intent(in) :: chempot_ext(nat, n_a)
  real(wp), intent(in) :: e_gap_ext(nat, n_a)
  real(wp), intent(out) :: eluao_ext(nat, n_a)
  integer :: a, k
  do k = 1, n_a
     do a = 1, nat
        eluao_ext(a, k) = chempot_ext(a, k) + e_gap_ext(a, k) / 2
     end do
  end do
end subroutine

end module 