module tblite_ml_features_collect_containers
implicit none
private
public :: collect_containers_caches, distribute_cache

contains
subroutine collect_containers_caches(rep, rcache, coulomb, ccache, hal, hcache, disp, dcache, inter, icache,&
    & contain_list)
    use tblite_container_list, only : taint, cache_list
    use tblite_container_cache, only : resize
    use tblite_container, only : container_list, container_type, container_cache
    use tblite_repulsion_effective, only : tb_repulsion
    use tblite_xtb_coulomb, only : tb_coulomb
    use tblite_disp, only : dispersion_type
    use tblite_classical_halogen, only : halogen_correction
    type(container_cache), allocatable, intent(in) :: rcache, ccache, hcache, dcache, icache
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
    index = 1
    allocate(contain_list)
    if (allocated(rep)) then
        tmp_container = rep
        call contain_list%push_back(tmp_container, rcache)
        index = index + 1 
    end if
    if (allocated(coulomb)) then
        tmp_container = coulomb
        call contain_list%push_back(tmp_container, ccache)
        index = index + 1
    end if
    if (allocated(hal)) then         
        tmp_container = hal
        call contain_list%push_back(tmp_container, hcache)
        index = index +1
    end if
    if (allocated(disp)) then
        tmp_container = disp
        call contain_list%push_back(tmp_container, dcache)
        index = index + 1
    end if
    if (allocated(inter)) then   
      do i =1, inter%get_n_containers()
        call inter%pop(tmp_container)
        call contain_list%push_back(tmp_container, icache)
        index = index + 1
      end do
    end if
end subroutine

subroutine distribute_cache(rcache, ccache, hcache, dcache, icache, contain_list)
    use tblite_container, only : container_cache, container_type
    use tblite_container_list, only : container_list
    use tblite_repulsion, only : tb_repulsion
    use tblite_xtb_coulomb, only : tb_coulomb
    use tblite_scf_iterator, only : get_electronic_energy, reduce
    use tblite_disp_d3, only : d3_dispersion
    use tblite_classical_halogen, only : halogen_correction
    use tblite_disp_d4, only : d4_dispersion
    type(container_cache), allocatable, intent(inout) :: rcache, ccache, hcache, dcache, icache
    type(container_list), allocatable, intent(inout) :: contain_list
    type(container_cache), allocatable :: cache
    class(container_type), allocatable :: cont
    integer :: i
    do i = 1, contain_list%get_n_containers()
        call contain_list%pop(cont, cache=cache)
        select type(cont)
        type is (tb_repulsion)
            rcache = cache
        type is (tb_coulomb)
            ccache = cache
        type is (halogen_correction)
            hcache = cache
        type is (d3_dispersion)
            dcache = cache
        type is (d4_dispersion)
            dcache = cache
        class default
            icache = cache
        end select
    end do
end subroutine
end module