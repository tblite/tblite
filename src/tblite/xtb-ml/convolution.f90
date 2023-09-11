module tblite_xtbml_convolution
    use mctc_env, only : wp
   use mctc_io, only : structure_type
   
  
    implicit none
    private

    type xtbml_convolution_type
        real(wp), allocatable :: rcov(:)
        real(wp), allocatable :: inv_cn_a(:, :, :)
        real(wp), allocatable :: a(:)
        real(wp) :: k1 = 16.0_wp
        integer :: n_a
    contains
        procedure :: setup
        procedure, private :: populate_inv_cn_array
        procedure, private :: get_rcov
    end type

    contains

    subroutine setup(self, mol, a_array)
        class(xtbml_convolution_type), intent(in) :: self
        type(structure_type), intent(in) :: mol
        real(wp), intent(in) :: a_array(:)

        self%a = a_array
        self%n_a = size(a_array)

        call self%populate_inv_cn_array(mol%nat, mol%at, mol%xyz)
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
    
    subroutine populate_inv_cn_array(self, nat, at, xyz)
        class(xtbml_convolution_type) :: self
        integer, intent(in) :: nat, at(nat)
        real(wp), intent(in) :: xyz(:, :)
        real(wp) :: result
        integer :: i, j, k, n_a
        n_a = size(self%a)
        if (allocated(self%rcov)) then
            deallocate (self%rcov)
        end if
        allocate (self%inv_rcov(nat), source=0.0_wp)
        call get_rcov(mol)
        if (allocated(self%inv_cn_a)) then
           deallocate (slef%inv_cn_a)
        end if
        allocate (self%inv_cn_a(nat, nat, n_a), source=0.0_wp)
        !$omp parallel do default(none) collapse(2)&
        !$omp shared(self%a,nat, at,xyz,self%inv_cn_a,n_a)&
        !$omp private(result,i,j,k)
        do k = 1, n_a
           do i = 1, nat
              do j = 1, nat
                 !if (i == j) cycle
                 call inv_cn(nat, i, j, at, xyz, self%a(k), result)
                 self%inv_cn_a(i, j, k) = result
              end do
           end do
        end do
        !$omp end parallel do
     
     end subroutine populate_inv_cn_array
    
    subroutine inv_cn(nat, a, b, at, xyz, dampening_fact, result)
       
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
    

end module tblite_xtbml_convolution