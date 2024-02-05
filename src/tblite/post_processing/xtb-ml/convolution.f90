module tblite_xtbml_convolution
    use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_ncoord_xtbml, only : xtbml_ncoord_type, new_xtbml_ncoord
  
    implicit none
    private
    real(wp) :: k1 = 16.0_wp
    public :: xtbml_convolution_type
    type :: xtbml_convolution_type
        real(wp), allocatable :: rcov(:)
        real(wp), allocatable :: a(:)
        real(wp), allocatable :: cn(:, :)
        integer :: n_a
        real(wp), allocatable :: kernel(:, :, :)
        character(len=:), allocatable :: label
    contains
        procedure :: setup
        procedure :: compute_kernel
        procedure, private :: populate_kernel
        procedure, private :: get_rcov
        procedure, private :: compute_cn
        procedure :: info
    end type
    character(len=*), parameter :: label = "CN-based convolution"


    contains

    subroutine setup(self)
        class(xtbml_convolution_type), intent(inout) :: self
        self%label = label
    end subroutine    

    subroutine compute_kernel(self, mol)
        class(xtbml_convolution_type), intent(inout) :: self
        type(structure_type), intent(in) :: mol
 
        call self%get_rcov(mol)
        call self%populate_kernel(mol%nat, mol%id, mol%xyz)
        call self%compute_cn(mol)
        
    end subroutine

    subroutine compute_cn(self, mol)
        
        class(xtbml_convolution_type) :: self
        type(structure_type) :: mol
        type(xtbml_ncoord_type) :: ncoord_xtbml
        integer :: i, n_a
        n_a = size(self%a)
        allocate(self%cn(mol%nat, n_a), source=0.0_wp)
        do i = 1, n_a
            call new_xtbml_ncoord(ncoord_xtbml, mol, a=self%a(i))
            call ncoord_xtbml%get_cn(mol, self%cn(:, i))
        end do 
    end subroutine

    subroutine get_rcov(self,mol)
        use tblite_data_covrad, only : get_covalent_rad
        class(xtbml_convolution_type) :: self
        type(structure_type), intent(in) :: mol
        if (allocated(self%rcov)) then
           deallocate (self%rcov)
        end if
        allocate (self%rcov(mol%nid), source=0.0_wp)
        self%rcov(:) = get_covalent_rad(mol%num)
    end subroutine
    
    subroutine populate_kernel(self, nat, at, xyz)
        class(xtbml_convolution_type) :: self
        integer, intent(in) :: nat, at(nat)
        real(wp), intent(in) :: xyz(:, :)
        real(wp) :: result
        integer :: i, j, k, n_a
        n_a = size(self%a)
    
        if (allocated(self%kernel)) then
           deallocate (self%kernel)
        end if
        allocate (self%kernel(nat, nat, n_a), source=0.0_wp)
        !$omp parallel do default(none) collapse(2)&
        !$omp shared(self, nat, at, xyz, n_a)&
        !$omp private(result,i,j,k)
        do k = 1, n_a
           do i = 1, nat
              do j = 1, nat
                 !if (i == j) cycle
                 call inv_cn(self, nat, i, j, at, xyz, self%a(k), result)
                 self%kernel(i, j, k) = result
              end do
           end do
        end do
        
     
     end subroutine populate_kernel
    
    subroutine inv_cn(self, nat, a, b, at, xyz, dampening_fact, result)
        type(xtbml_convolution_type) :: self
        integer, intent(in) :: a, b, nat
        integer, intent(in) :: at(:)
        real(wp), intent(in)  :: xyz(3, nat), dampening_fact
        real(wp), intent(out) :: result
        real(wp) :: rab(3), r, rco, r2
     
        result = 0.0_wp
     
        rab = xyz(1:3, a) - xyz(1:3, b)
        r2 = sum(rab**2)
        r = sqrt(r2)
     
        rco = dampening_fact*(self%rcov(at(a)) + self%rcov(at(b)))
     
        result = 1.0_wp/exp_count(k1, r, rco)
     
    end subroutine
    
pure elemental function exp_count(k, r, r0) result(count)
       real(wp), intent(in) :: k
       real(wp), intent(in) :: r
       real(wp), intent(in) :: r0
       real(wp) :: count
       count = 1.0_wp/(1.0_wp + exp(-k*(r0/r - 1.0_wp)))
    end function exp_count
    
    pure function info(self, verbosity, indent) result(str)
    !> Instance of the interaction container
    class(xtbml_convolution_type), intent(in) :: self
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

end module tblite_xtbml_convolution