! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/xtb-ml/type.f90
module xtbml_class
    use mctc_env, only : wp
    use tblite_wavefunction, only : wavefunction_type
    use mctc_io, only : structure_type
    use tblite_basis_type, only : basis_type
    use tblite_results, only : results_type
    use tblite_integral_type, only : integral_type
    use tblite_xtb_calculator, only : xtb_calculator
    use tblite_container, only : container_cache
    use tblite_results, only: results_type
    use tblite_timer, only: timer_type, format_time
    implicit none
    character(len=12), parameter :: toml_file = "a_array.toml"
    private
    type(timer_type) :: timer

    type,public,abstract :: xtbml_type
        integer :: n_features
        character(len=30),allocatable :: feature_labels(:) 
        character(len=30),allocatable :: delta_labels(:) 
        real(wp),allocatable :: a(:)
        real(wp),ALLOCATABLE ::  w_tot  (:)
        !
        !
        real(wp),ALLOCATABLE ::  cn_atom  (:) 
        real(wp),ALLOCATABLE ::  delta_cn  (:,:)
        !
        !
        real(wp),ALLOCATABLE ::  mulliken_shell (:)
        real(wp),ALLOCATABLE ::  dipm_shell (:)
        real(wp),ALLOCATABLE ::  qm_shell (:)
        real(wp),ALLOCATABLE ::  partial_charge_atom  (:)
        real(wp),ALLOCATABLE ::  delta_partial_charge  (:,:)
        real(wp),ALLOCATABLE ::  dipm_atom  (:)
        real(wp),ALLOCATABLE ::  delta_dipm  (:,:)
        real(wp),ALLOCATABLE ::  qm_atom  (:)
        real(wp),ALLOCATABLE ::  delta_qm  (:,:)
        real(wp),ALLOCATABLE ::  delta_dipm_e  (:,:)
        real(wp),ALLOCATABLE ::  delta_qm_e  (:,:)
        real(wp),ALLOCATABLE ::  delta_dipm_Z  (:,:)
        real(wp),ALLOCATABLE ::  delta_qm_Z  (:,:)
        !
        !> shell dipm xyz
        real(wp),ALLOCATABLE ::  dipm_shell_xyz (:,:)
        !> dipm xyz
        real(wp),ALLOCATABLE ::  dipm_atom_xyz  (:,:)
        !> delta dipm xyz
        real(wp),ALLOCATABLE ::  delta_dipm_xyz  (:,:,:)
        !> shell qm xyz
        real(wp),ALLOCATABLE ::  qm_shell_xyz (:,:)
        !> qm xyz
        real(wp),ALLOCATABLE ::  qm_atom_xyz  (:,:)
        !> delta qm xyz
        real(wp),ALLOCATABLE ::  delta_qm_xyz  (:,:,:)
        !
        !> delta dipm only electron effect
        real(wp),ALLOCATABLE ::  delta_dipm_e_xyz (:,:,:)
        !> delta qm only electron effect
        real(wp),ALLOCATABLE ::  delta_qm_e_xyz (:,:,:)
        !> delta dipm only nuclear effect
        real(wp),ALLOCATABLE ::  delta_dipm_Z_xyz (:,:,:)
        !> delta qm only nuclear effect
        real(wp),ALLOCATABLE ::  delta_qm_Z_xyz (:,:,:)
        !
        !
        !seperate alpha and beta
        real(wp),ALLOCATABLE ::  response  (:)
        real(wp),ALLOCATABLE ::  egap  (:)
        real(wp),ALLOCATABLE ::  chempot  (:)
        real(wp),ALLOCATABLE ::  ehoao_a  (:)
        real(wp),ALLOCATABLE ::  eluao_a  (:)
        real(wp),ALLOCATABLE ::  ehoao_b  (:)
        real(wp),ALLOCATABLE ::  eluao_b  (:)
        real(wp),ALLOCATABLE ::  e_rep_atom  (:)
        real(wp),ALLOCATABLE ::  e_EHT  (:)
        real(wp),ALLOCATABLE ::  e_disp_2  (:)
        real(wp),ALLOCATABLE ::  e_disp_3  (:)
        real(wp),ALLOCATABLE ::  e_ies_ixc  (:)
        real(wp),ALLOCATABLE ::  e_aes  (:)
        real(wp),ALLOCATABLE ::  e_axc (:)
        !energy based features; extensions
        real(wp),ALLOCATABLE ::  delta_chempot(:,:)
        real(wp),ALLOCATABLE ::  delta_egap(:,:)
        real(wp),ALLOCATABLE ::  delta_eluao(:,:)
        real(wp),ALLOCATABLE ::  delta_ehoao(:,:)
    contains
        procedure,private :: allocate => allocate_ml
        procedure :: get_geometry_density_based
        procedure :: get_energy_based
        procedure :: get_extended_frontier
        !> Generate the xtbml features
        procedure(get_xtbml),deferred :: get_xtbml
        procedure(pack_res),deferred :: pack_res
        procedure :: compute_extended
        procedure :: print_out
        procedure :: pop_a
        procedure :: print_timer
    end type xtbml_type
        
    !<
    abstract interface
        !> Routine that computes the xtbml features
        subroutine get_xtbml(self,mol,wfn,integrals,erep,calc,ccache,dcache,prlevel,a_array,res)
            import :: wp, structure_type, wavefunction_type,integral_type,&
                    xtb_calculator,container_cache,results_type,xtbml_type
            class(xtbml_type),intent(inout) :: self
            !> Molecular structure data
            type(structure_type), intent(in) :: mol
            !> Wavefunction strcuture data
            type(wavefunction_type), intent(in) :: wfn
            type(integral_type) :: integrals
            !> Single-point calculator
            type(xtb_calculator), intent(in) :: calc
            type(container_cache),intent(inout) :: ccache,dcache
            type(results_type), intent(inout) :: res
            real(wp), INTENT(IN) ::  erep(mol%nat)
            integer, intent(in) :: prlevel
            real(wp), intent(in),allocatable  :: a_array(:)
       
        end subroutine get_xtbml
        !> Routine to pack the xtbml features into the result container
        subroutine pack_res(self,nat,nsh_tot,at2nsh,e_tot,labels,res)
            import :: wp,results_type, xtbml_type
            class(xtbml_type),intent(inout) :: self
            integer, intent(in) :: nat,nsh_tot,at2nsh(nat)
            real(wp), intent(in) :: e_tot
            type(results_type),intent(inout) :: res
            character(len=30),intent(in) :: labels(:)
        end subroutine pack_res
            
    end interface
    contains

    subroutine allocate_ml(self,nat,nshell,n_a)
        implicit none
        class(xtbml_type) :: self
        integer,intent(in)  :: nat,nshell,n_a
        
        allocate( self%w_tot (nat),     source = 0.0_wp )
        !
        allocate( self%cn_atom (nat),   source = 0.0_wp )
        allocate( self%delta_cn (nat,n_a),  source = 0.0_wp )
        !
        allocate(self%mulliken_shell(nshell) ,source = 0.0_wp)
        allocate(self%dipm_shell(nshell) ,source = 0.0_wp)
        allocate(self%qm_shell(nshell) ,source = 0.0_wp)
        allocate(self%partial_charge_atom(nat) ,source = 0.0_wp)
        allocate(self%delta_partial_charge(nat,n_a) ,source = 0.0_wp)
        allocate(self%dipm_atom(nat) ,source = 0.0_wp)
        allocate(self%delta_dipm(nat,n_a) ,source = 0.0_wp)
        allocate(self%qm_atom(nat) ,source = 0.0_wp)
        allocate(self%delta_qm(nat,n_a) ,source = 0.0_wp)
        allocate(self%delta_dipm_e(nat,n_a) ,source = 0.0_wp)
        allocate(self%delta_qm_e(nat,n_a) ,source = 0.0_wp)
        allocate(self%delta_dipm_Z(nat,n_a) ,source = 0.0_wp)
        allocate(self%delta_qm_Z(nat,n_a) ,source = 0.0_wp)
        !
        allocate(self%dipm_shell_xyz(3,nshell),source = 0.0_wp)
        allocate(self%dipm_atom_xyz(3,nat),source = 0.0_wp)
        allocate(self%delta_dipm_xyz(3,nat,n_a),source = 0.0_wp)
        allocate(self%qm_shell_xyz(6,nshell),source = 0.0_wp)
        allocate(self%qm_atom_xyz(6,nat),source = 0.0_wp)
        allocate(self%delta_qm_xyz(6,nat,n_a),source = 0.0_wp)
        !
        allocate(self%delta_dipm_e_xyz(3,nat,n_a),source = 0.0_wp)
        allocate(self%delta_qm_e_xyz(6,nat,n_a),source = 0.0_wp)
        allocate(self%delta_dipm_Z_xyz(3,nat,n_a),source = 0.0_wp)
        allocate(self%delta_qm_Z_xyz(6,nat,n_a),source = 0.0_wp)
        !
        !
        !seperate alpha and beta
        allocate(self%response(nat),source = 0.0_wp)
        allocate(self%egap(nat),source = 0.0_wp)
        allocate(self%chempot(nat),source = 0.0_wp)
        allocate(self%ehoao_a(nat),source = 0.0_wp)
        allocate(self%eluao_a(nat),source = 0.0_wp)
        allocate(self%ehoao_b(nat),source = 0.0_wp)
        allocate(self%eluao_b(nat),source = 0.0_wp)
        !
        allocate(self%e_rep_atom(nat),source = 0.0_wp)
        allocate(self%e_EHT(nat),source = 0.0_wp)
        allocate(self%e_disp_2(nat),source = 0.0_wp)
        allocate(self%e_disp_3(nat),source = 0.0_wp)
        allocate(self%e_ies_ixc(nat),source = 0.0_wp)
        allocate(self%e_aes(nat),source = 0.0_wp)
        allocate(self%e_axc(nat),source = 0.0_wp)
        !Extensions
        allocate( self%delta_chempot(nat,n_a),       source = 0.0_wp )
        allocate( self%delta_egap(nat,n_a),       source = 0.0_wp )
        allocate( self%delta_eluao(nat,n_a),       source = 0.0_wp )
        allocate( self%delta_ehoao(nat,n_a),       source = 0.0_wp )
    end subroutine allocate_ml

    subroutine get_geometry_density_based(self,mol,wfn,integrals,calc)
        use tblite_data_covrad, only: get_covalent_rad
        use tblite_ncoord_exp, only:  new_exp_ncoord, exp_ncoord_type
        use xtbml_functions
        use tblite_wavefunction_mulliken, only: get_mulliken_shell_multipoles
        class(xtbml_type) :: self
        !> Molecular structure data
        type(structure_type), intent(in) :: mol
        !> Wavefunction strcuture data
        type(wavefunction_type), intent(in) :: wfn
        type(integral_type) :: integrals
        !> Single-point calculator
        type(xtb_calculator), intent(in) :: calc
        type(exp_ncoord_type) :: ncoord_exp
        integer :: mu
        real(wp) :: z(mol%nat),dipm_shell_tmp(3,calc%bas%nsh,1),qm_shell_tmp(6,calc%bas%nsh,1)
        !allocate self type
        call timer%push("total")
        call timer%push("geometric")
        call self%allocate(mol%nat,calc%bas%nsh,size(self%a))
        call get_rcov(mol)
        call new_exp_ncoord(ncoord_exp,mol)
        call ncoord_exp%get_cn(mol,self%cn_atom)
        call timer%pop()
        
        call timer%push("density") 
        !shellwise mulliken charges
        call mulliken_shellwise(calc%bas%nao,calc%bas%nsh,calc%bas%ao2sh,wfn%density(:,:,wfn%nspin),&
        integrals%overlap,self%mulliken_shell)
        !call sum_up_spin(wfn%qat,self%partial_charge_atom)
        call mol_set_nuclear_charge(mol%nat,mol%num,mol%id,z)

        do mu = 1, calc%bas%nsh
            self%partial_charge_atom(calc%bas%sh2at(mu)) = self%partial_charge_atom(calc%bas%sh2at(mu)) +  self%mulliken_shell(mu)
        enddo
        self%partial_charge_atom = -self%partial_charge_atom + z
        !multipole moments shellwise und then atomwise
        
        call get_mulliken_shell_multipoles(calc%bas,integrals%dipole, wfn%density, &
        & dipm_shell_tmp)
        self%dipm_shell_xyz = dipm_shell_tmp(:,:,1)
        call get_mulliken_shell_multipoles(calc%bas,integrals%quadrupole, wfn%density, &
        & qm_shell_tmp)
        self%qm_shell_xyz = qm_shell_tmp(:,:,1)
        
        self%dipm_atom_xyz = wfn%dpat(:,:,1)
        self%qm_atom_xyz = wfn%qpat(:,:,1)

        call comp_norm(calc%bas%nsh,dipm_shell_tmp,qm_shell_tmp,self%dipm_shell,self%qm_shell)
        call comp_norm(mol%nat,wfn%dpat,wfn%qpat,self%dipm_atom,self%qm_atom)
        call timer%pop()
        !call comp_norm(mol%nat,self%delta_dipm_xyz,self%delta_qm_xyz,self%delta_dipm,self%delta_qm)
        
    end subroutine get_geometry_density_based

    subroutine get_energy_based(self,mol,wfn,calc,integrals,ccache,dcache,erep,e_gfn2_tot)
        use tblite_container, only : container_cache
        use tblite_scf_iterator, only : get_electronic_energy,reduce
        use xtbml_functions, only: get_total_xtb_weights
        class(xtbml_type) :: self
        !> Molecular structure data
        type(structure_type), intent(in) :: mol
        !> Wavefunction strcuture data
        type(wavefunction_type), intent(in) :: wfn
        type(integral_type) :: integrals
        !> Single-point calculator
        type(xtb_calculator), intent(in) :: calc
        type(container_cache),intent(inout) :: ccache,dcache
        type(container_cache) :: dcache2
        real(wp),intent(in) ::erep(mol%nat)
        real(wp),intent(inout) :: e_gfn2_tot
        real(wp), allocatable :: e_ao(:),e_disp(:)
        call timer%push("energy")
        call calc%coulomb%aes2%get_AXC(mol,wfn,self%e_axc)
        call calc%coulomb%aes2%get_energy_aes_xtb(mol,ccache,wfn,self%e_aes)
        
        call calc%coulomb%es2%get_energy(mol,ccache,wfn,self%e_ies_ixc)
        call calc%coulomb%es3%get_energy(mol,ccache,wfn,self%e_ies_ixc)
        call calc%dispersion_3body%update(mol,dcache2)
        call calc%dispersion_3body%get_engrad(mol,dcache2,self%e_disp_3)
        
        allocate (e_disp(mol%nat),source=0.0_wp)
        call calc%dispersion%update(mol,dcache)
        call calc%dispersion%get_energy(mol,dcache,wfn,e_disp)

        self%e_disp_2(:) = e_disp(:) - self%e_disp_3(:)

        self%e_rep_atom = erep
        
        !Compute E_EHT
        allocate(e_ao(calc%bas%nao),source=0.0_wp)
        call get_electronic_energy(integrals%hamiltonian,wfn%density,e_ao)
        call reduce(self%e_EHT,e_ao,calc%bas%ao2at)
        !Compute partition weights based on total energy expression
        call get_total_xtb_weights(mol%nat,self%e_EHT,self%e_rep_atom,self%e_disp_2,self%e_disp_3,&
                                self%e_ies_ixc,self%e_aes,self%e_axc,self%w_tot,e_gfn2_tot)
        call timer%pop()
    end subroutine get_energy_based

    subroutine print_out(self,out,nat,at,id2at,res)
        class(xtbml_type),intent(inout) :: self
        integer, INTENT(IN) :: out,nat,at(nat),id2at(nat)
        type(results_type), intent(in) :: res
        integer :: i,j
        write(out,'(a)', advance="no") "Atom,w_xtb_tot"//','
        do i=1, self%n_features-1
            write(out,'(a)', advance="no") trim(res%xtbml_labels(i))//','
        end do
        write(out,'(a)') trim(res%xtbml_labels(self%n_features))

        do i = 1, nat
            write(out,'(i2,a)', advance="no") at(id2at(i)),',' 
            write(out,'(f12.8,a)', advance="no") res%w_xtbml(i),','
            do j = 1, self%n_features-1
                write(out,'(f14.8,a)', advance="no") res%ml_features(i,j),',' 
            end do
            write(out,'(f14.8)') res%ml_features(i,self%n_features)
        enddo
    end subroutine print_out


    subroutine get_extended_frontier(self,mol,wfn)
        use xtbml_functions, only : get_beta,get_chem_pot_ext,get_e_gap_ext,&
        get_ehoao_ext,get_eluao_ext
        use mctc_io_convert, only : autoev
        class(xtbml_type),intent(inout) :: self
        !> Molecular structure data
        type(structure_type), intent(in) :: mol
        !> Wavefunction strcuture data
        type(wavefunction_type), intent(in) :: wfn
        real(wp) :: beta(mol%nat,mol%nat,size(self%a)),hl_gap
        integer :: n
        call timer%push("extended frontier")
        n =size(self%a)
        call get_beta(mol%nat,n,mol%id,mol%xyz,beta)

        call get_chem_pot_ext(mol%nat,n,beta,self%chempot,self%delta_chempot)
        hl_gap = (wfn%emo(wfn%homo(1),1) - wfn%emo(wfn%homo(1)+1,1))*autoev
        
        call get_e_gap_ext(mol%nat,n,hl_gap,beta,self%egap,self%delta_egap)

        call get_ehoao_ext(mol%nat,n,self%delta_chempot,self%delta_egap,self%delta_ehoao)

        call get_eluao_ext(mol%nat,n,self%delta_chempot,self%delta_egap,self%delta_eluao)
        call timer%pop()
        call timer%push("print out")
    end subroutine get_extended_frontier

    subroutine compute_extended(self,mol,wfn,calc)
        use xtbml_functions
        class(xtbml_type),intent(inout) :: self
        !> Molecular structure data
        type(structure_type), intent(in) :: mol
        !> Wavefunction strcuture data
        type(wavefunction_type), intent(in) :: wfn
        !> array of a values to used
        !> Single-point calculator
        type(xtb_calculator), intent(in) :: calc
        integer :: n
        real(wp) :: mull_charge_atomic(mol%nat)
        real(wp) :: z(mol%nat)
        call timer%push("extended")
        n = size(self%a)
        call populate_inv_cn_array(mol%nat,mol%id,mol%xyz,self%a)
        call mol_set_nuclear_charge(mol%nat,mol%num,mol%id,z)
        !compute delta CN
        call get_delta_cn(mol%nat,n,self%cn_atom,mol%id,mol%xyz,self%delta_cn)
        !delta partial charge
        call get_delta_partial(mol%nat,n,self%partial_charge_atom,mol%id,mol%xyz,&
        self%cn_atom,self%delta_partial_charge)
        !delta multipole moments
        call get_delta_mm(mol%nat,n,self%partial_charge_atom,wfn%dpat,wfn%qpat,mol%id,mol%xyz,&
        self%cn_atom,self%delta_dipm_xyz,self%delta_qm_xyz)

        call comp_norm_3(mol%nat,n,self%delta_dipm_xyz,self%delta_qm_xyz,self%delta_dipm,self%delta_qm)
        call get_delta_mm_Z(mol%nat,n,z,wfn%dpat,wfn%qpat,mol%id,mol%xyz,self%cn_atom,self%delta_dipm_Z_xyz,&
        self%delta_qm_Z_xyz)
        call sum_up_mulliken(mol%nat,calc%bas%nsh,calc%bas%sh2at,calc%bas%sh2at,self%mulliken_shell,mull_charge_atomic)
        call get_delta_mm_p(mol%nat,n,mull_charge_atomic,wfn%dpat,wfn%qpat,mol%id,mol%xyz,self%cn_atom,&
        self%delta_dipm_e_xyz,self%delta_qm_e_xyz)

        call comp_norm_3(mol%nat,n,self%delta_dipm_e_xyz,self%delta_qm_e_xyz,self%delta_dipm_e,self%delta_qm_e)
        call comp_norm_3(mol%nat,n,self%delta_dipm_Z_xyz,self%delta_qm_Z_xyz,self%delta_dipm_Z,self%delta_qm_Z)
        call timer%pop()
    end subroutine compute_extended

    subroutine pop_a(self)
        use tblite_os, only: file_exists
        use tblite_toml
        
        class(xtbml_type),intent(inout) :: self
        type(toml_table),allocatable :: table
        type(toml_table), pointer :: child
        type(toml_array), pointer :: array
        type(toml_error),allocatable :: error
        integer :: io,i
        
        if (file_exists(toml_file)) then
            open(file=toml_file,newunit=io,status="old")
            call toml_parse(table,io,error)
            close(io)
            if (allocated(error)) then
                print '(a)', "Error: "//error%message
                stop 1
            end if
            call get_value(table,"xtbml",child)
            call get_value(child, "a", array)
            
            allocate(self%a(len(array)))
            do i = 1, size(self%a)
                call get_value(array,i,self%a(i))
            enddo
        else
            allocate(self%a(1))
            self%a = [1.0_wp]
        end if

    end subroutine pop_a

    subroutine print_timer(self,ftime)
        use tblite_output_format, only : format_string
        class(xtbml_type) :: self
        real(wp),intent(in) :: ftime
        integer :: it
        real(wp) :: ttime, stime
        character(len=*), parameter :: label(*) = [character(len=20):: &
        & "geometric","density", "energy", "frontier", "extended", "extended frontier","print out"]
                call timer%pop()
                   ttime = timer%get("total")
                   write(*,*) (" total:"//repeat(" ", 16)//format_time(ttime))
                
                
                   do it = 1, size(label)
                      stime = timer%get(label(it))
                      if (stime <= epsilon(0.0_wp)) cycle
                      write(*,*)(" - "//label(it)//format_time(stime) &
                         & //" ("//format_string(int(stime/ttime*100), '(i3)')//"%)")
                   end do
                   write(*,*)(" - "//"frontier"//format_time(ftime) &
                   & //" ("//format_string(int(ftime/ttime*100), '(i3)')//"%)") 
                  write(*,*)
    end subroutine print_timer


end module xtbml_class
 