module tblite_ml_features_collect_containers
implicit none
private
public :: collect_containers_caches
contains
subroutine collect_containers_caches(rep, rcache, coulomb, ccache, hal, hcache, disp, dcache, inter, icache,&
    & contain_list, cach_list)
    use tblite_container_list, only : taint, cache_list
    use tblite_container_cache, only : resize
    use tblite_container, only : container_list, container_type, container_cache
    use tblite_repulsion_effective, only : tb_repulsion
    use tblite_xtb_coulomb, only : tb_coulomb
    use tblite_disp, only : dispersion_type
    use tblite_classical_halogen, only : halogen_correction
    type(cache_list), pointer, intent(inout) :: cach_list
    type(container_cache) :: rcache, ccache, hcache, dcache, icache
    type(container_list), allocatable, intent(inout) :: contain_list
    type(container_list), allocatable :: inter 
    class(container_type), allocatable :: tmp_container
    type(tb_repulsion), allocatable :: rep
    !> Collection of all Coulombic interactions
    type(tb_coulomb), allocatable :: coulomb
    !> Halogen bonding correction
    type(halogen_correction), allocatable :: hal
    !> London-dispersion interaction
    class(dispersion_type), allocatable :: disp

    integer :: index, i
    allocate(contain_list)
  
    index = 1
    if (allocated(rep)) then
        tmp_container = rep
        call contain_list%push_back(tmp_container)
        call taint(rcache, cach_list)
        call resize(cach_list%list, index)
        index = index + 1 
    end if
    if (allocated(coulomb)) then
        tmp_container = coulomb
        call contain_list%push_back(tmp_container)
        call taint(ccache, cach_list)
        call resize(cach_list%list, index)
        index = index + 1
    end if
    if (allocated(hal)) then         
        tmp_container = hal
        call contain_list%push_back(tmp_container)
        call taint(hcache, cach_list)
        call resize(cach_list%list, index)
        index = index +1
    end if
    if (allocated(disp)) then
        tmp_container = disp
        call contain_list%push_back(tmp_container)
        call taint(dcache, cach_list)
        call resize(cach_list%list, index)
        index = index + 1
    end if
    if (allocated(inter)) then   
      do i =1, inter%get_n_containers()
        call inter%pop(tmp_container)
        call contain_list%push_back(tmp_container)
        call taint(icache, cach_list)
        call resize(cach_list%list, index) 
        index = index + 1
      end do
    end if
    end subroutine
  
end module