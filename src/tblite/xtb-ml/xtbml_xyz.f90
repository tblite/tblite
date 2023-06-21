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

!> @file tblite/xtb-ml/xtbml_xyz.f90
module xtbml_xyz

    use mctc_env, only : wp
    use tblite_wavefunction, only : wavefunction_type
    use mctc_io, only : structure_type
    use tblite_basis_type, only : basis_type
    use tblite_results, only : results_type
    use tblite_integral_type, only : integral_type
    use tblite_xtb_calculator, only : xtb_calculator
    use tblite_container, only : container_cache
    use tblite_scf_iterator, only : get_electronic_energy,reduce
    use tblite_results, only: results_type
    use mctc_io_convert, only : autoev
    
    use xtbml_class, only: xtbml_type
    integer,parameter :: n_delta = 39
    logical, parameter :: debug = .false.
type, public, extends(xtbml_type) :: xtbml_xyz_type
        
contains
    procedure :: get_xtbml
    procedure :: pack_res
end type xtbml_xyz_type
contains
    subroutine get_xtbml(self,mol,wfn,integrals,erep,calc,ccache,dcache,prlevel,a_array,res)
        use xtbml_functions
        use tblite_timer, only: timer_type
        implicit none
        class(xtbml_xyz_type),intent(inout) :: self
        type(structure_type), intent(in) :: mol
        type(integral_type) :: integrals
        !> Single-point calculator
        type(xtb_calculator), intent(in) :: calc
        type(wavefunction_type), intent(in) :: wfn
        type(container_cache),intent(inout) :: ccache,dcache
        type(results_type),intent(inout) :: res
        real(wp), intent(in), allocatable :: a_array(:)
        real(wp), INTENT(IN) ::  erep(mol%nat)
        integer, intent(in) :: prlevel
        real(wp) :: e_gfn2_tot,ftime
        integer :: ml_out
        logical :: print_afo
        type(timer_type) :: timer
        character(len=30), allocatable :: tmp_labels(:)
        if (allocated(a_array)) then
            allocate(self%a(size(a_array)))
            self%a = a_array
        else
            call self%pop_a
        endif
        self%n_features = 103 + ((size(self%a)-1)*n_delta)
        self%a = 1.0_wp
        allocate(self%feature_labels(self%n_features))
        tmp_labels = [ character(len=30) :: "CN",&
        &"q_s","q_p","q_d",&
        &"dipm_s","dipm_p","dipm_d",&
        &"dipm_s_x","dipm_s_y","dipm_s_z",&
        &"dipm_p_x","dipm_p_y","dipm_p_z",&
        &"dipm_d_x","dipm_d_y","dipm_d_z",&
        &"qm_s","qm_p","qm_d",& ! 19
        &"qm_s_xx","qm_s_xy","qm_s_yy","qm_s_xz","qm_s_yz","qm_s_zz",&
        &"qm_p_xx","qm_p_xy","qm_p_yy","qm_p_xz","qm_p_yz","qm_p_zz",&
        &"qm_d_xx","qm_d_xy","qm_d_yy","qm_d_xz","qm_d_yz","qm_d_zz",& !37
        &"p_A","dipm_A",& !38,39
        &"dipm_A_x","dipm_A_y","dipm_A_z",& !40 41 42
        &"qm_A",& !43
        &"qm_A_xx","qm_A_xy","qm_A_yy","qm_A_xz","qm_A_yz","qm_A_zz",& !44 45 46 47 48 49
        &"response","gap","chem.pot","HOAO_a","LUAO_a","HOAO_b","LUAO_b",& ! 50 51 52 53 54 55 56
        &"E_repulsion","E_EHT","E_disp_2","E_disp_3","E_ies_ixc","E_aes","E_axc","E_tot"]

        allocate(self%delta_labels(n_delta))
        self%delta_labels = [ character(len=30) :: "delta_CN","delta_p_A",&! 1 2
        &"delta_dipm_A",&! 3
        &"delta_dipm_A_x","delta_dipm_A_y","delta_dipm_A_z",& ! 4 5 6 
        &"delta_qm",& ! 7
        &"delta_qm_A_xx","delta_qm_A_xy","delta_qm_A_yy","delta_qm_A_xz","delta_qm_A_yz","delta_qm_A_zz",& !8 9 10 11 12 13
        &"delta_dipm_e",& ! 14
        &"delta_dipm_e_x","delta_dipm_e_y","delta_dipm_e_z",& ! 15 16 17
        &"delta_qm_e",& ! 18
        &"delta_qm_e_xx","delta_qm_e_xy","delta_qm_e_yy","delta_qm_e_xz","delta_qm_e_yz","delta_qm_e_zz",& !19 20 21 22 23 24
        &"delta_dipm_Z",& !25
        &"delta_dipm_Z_x","delta_dipm_Z_y","delta_dipm_Z_z",&
        &"delta_qm_Z",&
        &"delta_qm_Z_xx","delta_qm_Z_xy","delta_qm_Z_yy","delta_qm_Z_xz","delta_qm_Z_yz","delta_qm_Z_zz",&
        &"delta_gap","delta_chem_pot","delta_HOAO","delta_LUAO"]
        
        !get individual coulombic energy contributions in an atomwise vector
        call self%get_geometry_density_based(mol,wfn,integrals,calc)
        call self%get_energy_based(mol,wfn,calc,integrals,ccache,dcache,erep,e_gfn2_tot)

        print_afo = .false.
        if (prlevel > 1) then
            print_afo = .true.
        end if
        call timer%push("frontier")
        call atomic_frontier_orbitals(mol%nat,calc%bas%nao,wfn%focca,wfn%foccb,wfn%emo(:,1)*autoev,calc%bas%ao2at,wfn%coeff(:,:,1),&
        integrals%overlap(:,:),self%response,self%egap,self%chempot,self%ehoao_a,self%eluao_a,self%ehoao_b,self%eluao_b,print_afo)
        call timer%pop()
        ftime = timer%get("frontier")
        call self%compute_extended(mol,wfn,calc)
        call self%get_extended_frontier(mol,wfn)
        call self%pack_res(mol%nat,calc%bas%nsh,calc%bas%nsh_at,e_gfn2_tot,tmp_labels,res)
        !res%w_xtbml = self%w_tot
        if (prlevel > 1) then
            ml_out = 42
            open(file='ml_feature_tblite.csv', newunit=ml_out)
            call self%print_out(ml_out,mol%nat,mol%num,mol%id,res)
        endif
        if (debug) then 
            call self%print_timer(ftime)
         endif
    end subroutine get_xtbml

    subroutine pack_res(self,nat,nsh_tot,at2nsh,e_tot,labels,res)
        use xtbml_functions, only : pack_mult_xyz, pack_mult_xyz_shell,pack_shellwise
        use tblite_output_format, only: format_string
        implicit none
        integer, intent(in) :: nat,nsh_tot,at2nsh(nat)
        real(wp), intent(in) :: e_tot
        type(results_type),intent(inout) :: res
        class(xtbml_xyz_type), intent(inout) :: self
        integer :: i, nsh, k, n_other, offset
        character(len=30),intent(in) :: labels(:)

        allocate(res%ml_features(nat,self%n_features),source=0.0_wp)

        res%ml_features(:,1) = self%cn_atom(:)
        call pack_shellwise(self%mulliken_shell,res,2,at2nsh,nat) !2-4
        call pack_shellwise(self%dipm_shell,res,5,at2nsh,nat) ! 5-7
        call pack_mult_xyz_shell(self%dipm_shell_xyz,res,8,nat,at2nsh) !packs xyz for s to d shell 8-16
        call pack_shellwise(self%qm_shell,res,17,at2nsh,nat) ! 17-19
        call pack_mult_xyz_shell(self%qm_shell_xyz,res,20,nat,at2nsh) !20-37
        res%ml_features(:,38) = self%partial_charge_atom(:)
        res%ml_features(:,39) = self%dipm_atom(:)
        call pack_mult_xyz(self%dipm_atom_xyz,res,40,nat) !41-43
        res%ml_features(:,43) = self%qm_atom(:)
        call pack_mult_xyz(self%qm_atom_xyz,res,44,nat) !45-50
        res%ml_features(:,50) = self%response(:)
        res%ml_features(:,51) = self%egap(:)
        res%ml_features(:,52) = self%chempot(:)
        res%ml_features(:,53) = self%ehoao_a(:)
        res%ml_features(:,54) = self%eluao_a(:)
        res%ml_features(:,55) = self%ehoao_b(:)
        res%ml_features(:,56) = self%eluao_b(:)
        res%ml_features(:,57) = self%e_rep_atom(:)
        res%ml_features(:,58) = self%e_EHT(:)
        res%ml_features(:,59) = self%e_disp_2(:)
        res%ml_features(:,60) = self%e_disp_3(:)
        res%ml_features(:,61) = self%e_ies_ixc(:)
        res%ml_features(:,62) = self%e_aes(:)
        res%ml_features(:,63) = self%e_axc(:)
        res%ml_features(:,64) = e_tot
        n_other = 64
        
        do i =1, n_other
            self%feature_labels(i) = trim(labels(i))
        enddo
        
        do k =1 , size(self%a)
            offset = n_other + (k-1)*n_delta
            res%ml_features(:,offset+1) = self%delta_cn(:,k)
            res%ml_features(:,offset+2) = self%delta_partial_charge(:,k)
            res%ml_features(:,offset+3) = self%delta_dipm(:,k)
            call pack_mult_xyz(self%delta_dipm_xyz(:,:,k),res,offset+4,nat) 
            res%ml_features(:,offset+7) = self%delta_qm(:,k)
            call pack_mult_xyz(self%delta_qm_xyz(:,:,k),res,offset+8,nat) 
            res%ml_features(:,offset+14) = self%delta_dipm_e(:,k)
            call pack_mult_xyz(self%delta_dipm_e_xyz(:,:,k),res,offset+15,nat)  
            res%ml_features(:,offset+18) = self%delta_qm_e(:,k)
            call pack_mult_xyz(self%delta_qm_e_xyz(:,:,k),res,offset+19,nat) 
            res%ml_features(:,offset+25) = self%delta_dipm_Z(:,k)
            call pack_mult_xyz(self%delta_dipm_Z_xyz(:,:,k),res,offset+26,nat)   
            res%ml_features(:,offset+29) = self%delta_qm_Z(:,k)
            call pack_mult_xyz(self%delta_qm_Z_xyz(:,:,k),res,offset+30,nat)  
            res%ml_features(:,offset+36) = self%delta_egap(:,k)
            res%ml_features(:,offset+37) = self%delta_chempot(:,k)
            res%ml_features(:,offset+38) = self%delta_ehoao(:,k)
            res%ml_features(:,offset+39) = self%delta_eluao(:,k)
            if (self%a(1) .ne. 1.0_wp) then 
            do i = 1, n_delta
                self%feature_labels(offset+i) = trim(self%delta_labels(i)) //'_' //adjustl(format_string(self%a(k),'(f12.2)'))
            enddo
            else
                do i = 1, n_delta
                    self%feature_labels(offset+i) = trim(self%delta_labels(i))
                enddo
            endif
        enddo

        res%n_features = self%n_features
        call move_alloc(self%feature_labels,res%xtbml_labels)
        call move_alloc(self%w_tot,res%w_xtbml)
        
    end subroutine


end module xtbml_xyz
