module test_xtbml
    use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_context_type, only : context_type
   use tblite_basis_type
   use tblite_wavefunction , only : wavefunction_type, new_wavefunction
   use xtbml_base
   use xtbml_functions
   use tblite_integral_type, only : integral_type, new_integral
   use tblite_xtb_h0, only : get_selfenergy, get_hamiltonian, get_occupation, &
      & get_hamiltonian_gradient
   use tblite_wavefunction_mulliken, only : get_mulliken_shell_multipoles
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   use mctc_io_convert, only : aatoau
    implicit none
    private
    real(wp), parameter :: thr = 100*epsilon(1.0_wp)
    real(wp), parameter :: acc = 0.001_wp
    real(wp), parameter :: kt = 300.0_wp * 3.166808578545117e-06_wp
    real(wp), parameter :: thr2 = 10e-4
    public :: collect_xtbml
    contains
    
    subroutine collect_xtbml(testsuite)
    !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("xtbml-mulliken-charges", test_mulliken_charges_shell_h2p),&
      new_unittest("xtbml-dipm-mol-sum-up-h2+", test_dipm_shell_h2p),&
      new_unittest("xtbml-dipm-mol-sum-up-co2", test_dipm_shell_co2),&
      new_unittest("xtbml-qp-sum-up-benzene", test_qp_shell_benz),&
      new_unittest("xtbml-qp-dipm-benzene-high-a", test_qp_shell_benz_high_a),&
      new_unittest("xtbml-rot-translation", test_translation_rotation_co2)&
      ]

end subroutine collect_xtbml

    subroutine test_mulliken_charges_shell_h2p(error)

        use mctc_io_structure, only : new
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        type(context_type) :: ctx
        type(structure_type) :: mol
        type(xtb_calculator) :: calc
        type(wavefunction_type) :: wfn
        type(results_type) :: res
        real(wp) :: mulliken_shell(2)
        real(wp) :: energy 
        real(wp), parameter :: xyz(3, 2) = reshape((/0.0_wp,0.0_wp,0.35_wp,&
                                                            &0.0_wp,0.0_wp,-0.35_wp/),shape=(/3,2/))
        integer, parameter :: num(2) = (/1,1/)

        call new(mol,num,xyz*aatoau,uhf=1,charge=1.0_wp)
        call new_gfn2_calculator(calc,mol)
        call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
        calc%save_integrals = .true.
        call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0,results=res)
        
        call mulliken_shellwise(calc%bas%nao,calc%bas%nsh,calc%bas%ao2sh,wfn%density(:,:,wfn%nspin),&
        res%overlap,mulliken_shell)
        if (sum(mulliken_shell) - 1.0_wp > thr) then
            call test_failed(error, "Charge is not summing up to 1 electron for H2+")
            print'(3es21.14)', mulliken_shell
        end if

    end subroutine test_mulliken_charges_shell_h2p


    subroutine test_dipm_shell_h2p(error)

        use mctc_io_structure, only : new
        use tblite_adjlist, only : adjacency_list, new_adjacency_list
        use tblite_ncoord_exp, only:  new_exp_ncoord, exp_ncoord_type
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        type(context_type) :: ctx
        type(structure_type) :: mol
        type(xtb_calculator) :: calc
        type(wavefunction_type) :: wfn
        type(results_type) :: res
        
        integer, parameter :: nsh = 2, nat=2
        real(wp) :: mol_dipm(3), mol_dipm_delta(3),energy
        real(wp) :: partial(nat),dipm_xyz(3,nat),qm_xyz(6,nat)
        real(wp) :: delta_partial(nat),delta_dipm_xyz(3,nat),delta_qm_xyz(6,nat)
        
        
        real(wp), parameter :: xyz(3, nat) = reshape((/0.0_wp,0.0_wp,0.35_wp,&
                                                            &0.0_wp,0.0_wp,-0.35_wp/),shape=(/3,2/))
        integer, parameter :: num(nat) = (/1,1/)
        integer :: i
        

        call new(mol,num,xyz*aatoau,uhf=1,charge=1.0_wp)
        call new_gfn2_calculator(calc,mol)
        call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
        
        call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)

        call compute_mult_ext_mult(mol,wfn,calc,partial,delta_partial,dipm_xyz,delta_dipm_xyz,qm_xyz,delta_qm_xyz)
        mol_dipm = 0.0_wp
        do i= 1,2
            mol_dipm = mol_dipm+ dipm_xyz(:,i)+xyz(:,i)*partial(i)
        enddo 

        if (sum(mol_dipm)  > thr) then
            call test_failed(error, "Molecular dipole moment is non zero, for dipm")
            print'(3es21.14)', mol_dipm
        end if
        mol_dipm_delta = 0.0_wp
        do i= 1,2
            mol_dipm_delta = mol_dipm_delta + delta_dipm_xyz(:,i)+xyz(:,i)*delta_partial(i)
        enddo 

        if (sum(mol_dipm_delta)  > thr) then
            call test_failed(error, "Molecular dipole moment is non zero, for dipm_delta")
            print'(3es21.14)', mol_dipm_delta
        end if

        if (sum(mol_dipm_delta-mol_dipm)  > thr) then
            call test_failed(error, "Molecular dipole moment of extended and non-extended are not equal")
            print'(3es21.14)', mol_dipm_delta-mol_dipm
        end if
        
        

    end subroutine test_dipm_shell_h2p

    subroutine test_dipm_shell_co2(error)

        use mctc_io_structure, only : new
        use tblite_adjlist, only : adjacency_list, new_adjacency_list
        use tblite_ncoord_exp, only:  new_exp_ncoord, exp_ncoord_type
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        type(context_type) :: ctx
        type(structure_type) :: mol
        type(xtb_calculator) :: calc
        type(wavefunction_type) :: wfn
        type(integral_type) :: ints
        type(results_type) :: res
        type(adjacency_list) :: list
        integer, parameter :: nat=3
        real(wp) :: mol_dipm(3), mol_dipm_delta(3), energy
        real(wp) :: partial(nat),dipm_xyz(3,nat),qm_xyz(6,nat)
        real(wp) :: delta_partial(nat),delta_dipm_xyz(3,nat),delta_qm_xyz(6,nat)
        integer :: i,j

        real(wp), parameter :: xyz(3, nat) = reshape((/0.00000,0.00000,1.0000000,&
                                                            &0.0000,0.00000,-1.0000000,&
                                                            &0.000000,0.000000,0.000000/),shape=(/3,nat/))
        integer, parameter :: num(nat) = (/8,8,6/)
        mol_dipm = 0.0_wp
        mol_dipm_delta = 0.0_wp
        
        call new(mol,num,xyz*aatoau,uhf=0,charge=0.0_wp)
        call new_gfn2_calculator(calc,mol)
        call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
        calc%save_integrals = .true.
        call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0,results=res)

        call compute_mult_ext_mult(mol,wfn,calc,partial,delta_partial,dipm_xyz,delta_dipm_xyz,qm_xyz,delta_qm_xyz)
        
        do i= 1,nat
            do j = 1,3
                mol_dipm(j) = mol_dipm(j)+ dipm_xyz(j,i)+mol%xyz(j,i)*partial(i)
            enddo
        enddo 

        if (norm2(mol_dipm)  > thr) then
            call test_failed(error, "Molecular dipole moment is non zero, for dipm")
            print'(3es21.14)', mol_dipm
        end if

        do i= 1,nat
            do j = 1, 3
                mol_dipm_delta(j) = mol_dipm_delta(j) + delta_dipm_xyz(j,i)+mol%xyz(j,i)*delta_partial(i)
            enddo
        enddo 

        if (norm2(mol_dipm_delta)  > thr) then
            call test_failed(error, "Molecular dipole moment is non zero, for dipm_delta")
            print'(3es21.14)', mol_dipm_delta
        end if

        if (norm2(mol_dipm_delta-mol_dipm)  > thr) then
            call test_failed(error, "Molecular dipole moment of extended and non-extended are not equal")
            print'(3es21.14)', mol_dipm_delta-mol_dipm
        end if

    end subroutine test_dipm_shell_co2

    subroutine test_qp_shell_benz(error)

        use mctc_io_structure, only : new
        use tblite_adjlist, only : adjacency_list, new_adjacency_list
        use tblite_ncoord_exp, only:  new_exp_ncoord, exp_ncoord_type
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        type(context_type) :: ctx
        type(structure_type) :: mol
        type(xtb_calculator) :: calc
        type(wavefunction_type) :: wfn
        type(integral_type) :: ints
        type(results_type) :: res
        type(adjacency_list) :: list
        integer, parameter :: nat=12
        real(wp) :: mol_dipm(3), mol_dipm_delta(3), energy
        real(wp) :: partial(nat),dipm_xyz(3,nat),qm_xyz(6,nat)
        real(wp) :: delta_partial(nat),delta_dipm_xyz(3,nat),delta_qm_xyz(6,nat),mol_qm(6),delta_mol_qm(6)
        integer :: i,j

        real(wp), parameter :: xyz(3, nat) = reshape((/&
        &1.06880660023529,       -0.46478030005927,        0.00009471732781,&
        &2.45331325533661,       -0.46484444142679,        0.00042084312846,&
        &3.14559996373466,        0.73409173401801,        0.00008724271521,&
        &2.45340597395333,        1.93314591537421,       -0.00086301044874,&
        &1.06895328716368,        1.93321130458200,       -0.00141386731889,&
        &0.37663958202671,        0.73422879200405,       -0.00090929198808,&
        &0.52864276199175,       -1.40035288735680,        0.00049380958906,&
        &2.99337903563419,       -1.40047547903112,        0.00121759015506,&
        &4.22591259698321,        0.73408789322206,        0.00025398801936,&
        &2.99365942822711,        2.86866756976346,       -0.00166131228918,&
        &0.52879830433456,        2.86879139255056,       -0.00224874122149,&
        &-0.70367266962110,        0.73433126635962,       -0.00138296766859&
        &/),shape=(/3,nat/))
        integer, parameter :: num(nat) = (/6,6,6,6,6,6,1,1,1,1,1,1/)
        
        mol_dipm = 0.0_wp
        mol_dipm_delta = 0.0_wp
        call new(mol,num,xyz*aatoau,uhf=0,charge=0.0_wp)
        call new_gfn2_calculator(calc,mol)
        call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
        call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)

        call compute_mult_ext_mult(mol,wfn,calc,partial,delta_partial,dipm_xyz,delta_dipm_xyz,qm_xyz,delta_qm_xyz)
        
        do i= 1,nat
            do j = 1,3
                mol_dipm(j) = mol_dipm(j)+dipm_xyz(j,i)+ mol%xyz(j,i)*partial(i)
            enddo
        enddo
        
        do i= 1,nat
            do j =1,3
                mol_dipm_delta(j) = mol_dipm_delta(j) + delta_dipm_xyz(j,i)+mol%xyz(j,i)*delta_partial(i)
            enddo 
        enddo 

        if (norm2(mol_dipm_delta-mol_dipm)  > thr) then
            call test_failed(error, "Molecular dipole moment of extended and non-extended are not equal")
            print'(3es21.14)', norm2(mol_dipm_delta-mol_dipm)
        end if


        call compute_traceless_mol_qm(mol%nat,mol%xyz,partial,dipm_xyz,qm_xyz,mol_qm)


        call compute_traceless_mol_qm(mol%nat,mol%xyz,delta_partial,delta_dipm_xyz,delta_qm_xyz,delta_mol_qm)


        if (norm2(mol_qm-delta_mol_qm)  > thr2) then
            call test_failed(error, "Molecular quadrupole moment of extended and non-extended are not equal")
            print'(3es21.14)', norm2(mol_dipm_delta-mol_dipm)
        end if
    end subroutine test_qp_shell_benz

    subroutine test_qp_shell_benz_high_a(error)

        use mctc_io_structure, only : new
        use tblite_adjlist, only : adjacency_list, new_adjacency_list
        use tblite_ncoord_exp, only:  new_exp_ncoord, exp_ncoord_type
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        type(context_type) :: ctx
        type(structure_type) :: mol
        type(xtb_calculator) :: calc
        type(wavefunction_type) :: wfn
        type(integral_type) :: ints
        type(results_type) :: res
        type(adjacency_list) :: list
        integer, parameter :: nat=12
        real(wp) :: mol_dipm(3), mol_dipm_delta(3), energy
        real(wp) :: partial(nat),dipm_xyz(3,nat),qm_xyz(6,nat)
        real(wp) :: delta_partial(nat),delta_dipm_xyz(3,nat),delta_qm_xyz(6,nat),mol_qm(6),delta_mol_qm(6)
        integer :: i,j

        real(wp), parameter :: xyz(3, nat) = reshape((/&
        &1.06880660023529,       -0.46478030005927,        0.00009471732781,&
        &2.45331325533661,       -0.46484444142679,        0.00042084312846,&
        &3.14559996373466,        0.73409173401801,        0.00008724271521,&
        &2.45340597395333,        1.93314591537421,       -0.00086301044874,&
        &1.06895328716368,        1.93321130458200,       -0.00141386731889,&
        &0.37663958202671,        0.73422879200405,       -0.00090929198808,&
        &0.52864276199175,       -1.40035288735680,        0.00049380958906,&
        &2.99337903563419,       -1.40047547903112,        0.00121759015506,&
        &4.22591259698321,        0.73408789322206,        0.00025398801936,&
        &2.99365942822711,        2.86866756976346,       -0.00166131228918,&
        &0.52879830433456,        2.86879139255056,       -0.00224874122149,&
        &-0.70367266962110,        0.73433126635962,       -0.00138296766859&
        &/),shape=(/3,nat/))
        integer, parameter :: num(nat) = (/6,6,6,6,6,6,1,1,1,1,1,1/)
        
        mol_dipm = 0.0_wp
        mol_dipm_delta = 0.0_wp
        call new(mol,num,xyz*aatoau,uhf=0,charge=0.0_wp)
        call new_gfn2_calculator(calc,mol)
        call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
        !calc%save_integrals = .true.
        call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)

        call compute_mult_ext_mult_high_a(mol,wfn,calc,partial,delta_partial,dipm_xyz,delta_dipm_xyz,qm_xyz,delta_qm_xyz)
        
        do i= 1,nat
            do j = 1,3
                mol_dipm(j) = mol_dipm(j)+delta_dipm_xyz(j,i)+ mol%xyz(j,i)*partial(i)
            enddo
        enddo

        mol_dipm_delta = sum(delta_dipm_xyz,dim=2)

        if (norm2(mol_dipm_delta-mol_dipm)  > thr2) then
            call test_failed(error, "Molecular dipole moment of extended and non-extended are not equal")
            print'(3es21.14)', norm2(mol_dipm_delta-mol_dipm)
        end if

        call compute_traceless_mol_qm(mol%nat,mol%xyz,partial,dipm_xyz,qm_xyz,mol_qm)

        delta_mol_qm = sum(delta_qm_xyz,dim=2)

        if (norm2(mol_qm-delta_mol_qm)  > thr2) then
            call test_failed(error, "Molecular quadrupole moment of extended and non-extended are not equal")
            print'(3es21.14)', norm2(mol_qm-delta_mol_qm)
        end if
    end subroutine test_qp_shell_benz_high_a


    subroutine test_translation_rotation_co2(error)
        use mctc_io_structure, only : new
        use tblite_blas, only : gemm
        type(context_type) :: ctx
        type(structure_type) :: mol
        type(xtb_calculator) :: calc
        type(wavefunction_type) :: wfn
        !> Error handling
        integer,parameter :: nat = 3 
        type(error_type), allocatable, intent(out) :: error
        real(wp), parameter :: xyz(3, nat) = reshape((/0.00000,0.00000,1.0000000,&
                                                            &0.0000,0.00000,-1.0000000,&
                                                            &0.000000,0.000000,0.000000/),shape=(/3,nat/))
    
        integer, parameter :: num(nat) = (/8,8,6/)
        type(results_type) :: res
        real(wp) :: rot_matrix(3,3),xyz_rot(3,nat),energy, xyz_trans(3,nat)
        real(wp), allocatable :: xtbml(:,:),xtbml_rot(:,:),xtbml_trans(:,:)
        integer :: i

        call new(mol,num,xyz*aatoau,uhf=0,charge=0.0_wp)
        call new_gfn2_calculator(calc,mol)
        call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
        calc%xtbml = 1
        !calc%save_integrals = .true.
        call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0,results= res)

        allocate(xtbml(nat,size(res%ml_features,dim=2)),source=res%ml_features)
        !rotation by 45 degrees
        rot_matrix = (reshape((/&
        &1.0_wp,0.0_wp,0.0_wp,&
        &0.0_wp,0.52532_wp,-0.85090352453_wp,&
        &0.0_wp,0.85090352453_wp,0.52532_wp/),shape=(/3,3/)))
        call gemm(rot_matrix,xyz,xyz_rot)
        
        call new(mol,num,xyz_rot*aatoau,uhf=0,charge=0.0_wp)
        call new_gfn2_calculator(calc,mol)
        call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
        calc%xtbml = 1
        !calc%save_integrals = .true.
        call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0,results= res)

        allocate(xtbml_rot(nat,size(res%ml_features,dim=2)),source=res%ml_features)
        
        if (any(xtbml-xtbml_rot  > thr2)) then
            call test_failed(error, "Rotational invariance is not respected.")
            print'(3es21.14)', xtbml-xtbml_rot
        end if

        xyz_trans = xyz(:,:)
        
        do i = 1,nat
            xyz_trans(1,i) = xyz_trans(1,i) +5.0_wp
        enddo

        call new(mol,num,xyz_trans*aatoau,uhf=0,charge=0.0_wp)
        call new_gfn2_calculator(calc,mol)
        call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
        calc%xtbml = 1
        !calc%save_integrals = .true.
        call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0,results= res)

        allocate(xtbml_trans(nat,size(res%ml_features,dim=2)),source=res%ml_features)
        if (any(xtbml-xtbml_trans  > thr2)) then
            call test_failed(error, "Molecular quadrupole moment of extended and non-extended are not equal")
            print'(3es21.14)', xtbml-xtbml_trans
        end if
    end subroutine

    subroutine compute_mult_ext_mult(mol,wfn,calc,partial_charge,delta_partial,dipm,delta_dipm,qm,delta_qm)
        use tblite_ncoord_exp, only:  new_exp_ncoord, exp_ncoord_type

        type(context_type) :: ctx
        type(structure_type) :: mol
        type(xtb_calculator) :: calc
        type(wavefunction_type) :: wfn
        type(exp_ncoord_type) :: ncoord_exp
        real(wp) :: dipm_shell_tmp(3,calc%bas%nsh,1),partial_charge(mol%nat),a(1)
        real(wp) :: energy, z(mol%nat), cutoff,cn(mol%nat),delta_partial(mol%nat),delta_dipm_xyz(3,mol%nat,1),&
        &delta_qm_xyz(6,mol%nat,1)
        real(wp) :: dipm(3,mol%nat),delta_dipm(3,mol%nat)
        real(wp) :: qm(6,mol%nat),delta_qm(6,mol%nat) 
        integer :: n
       
        call new_gfn2_calculator(calc,mol)
        call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
        calc%save_integrals = .true.
        call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)
        
        !call mulliken_shellwise(calc%bas%nao,calc%bas%nsh,calc%bas%ao2sh,wfn%density(:,:,wfn%nspin),&
        !res%overlap,mulliken_shell)
        
        call mol_set_nuclear_charge(mol%nat,mol%num,mol%id,z)

        partial_charge =  wfn%qat(:,1) 
        
        call get_rcov(mol)
        call new_exp_ncoord(ncoord_exp,mol)
        call ncoord_exp%get_cn(mol,cn)
        
        a = 1.0_wp
        n = size(a)
        call populate_inv_cn_array(mol%nat,mol%id,mol%xyz,a)
        !delta partial charge
        call get_delta_partial(mol%nat,n,partial_charge,mol%id,mol%xyz,&
        cn,delta_partial)
        call get_delta_mm(mol%nat,n,partial_charge,wfn%dpat,wfn%qpat,mol%id,mol%xyz,&
        cn,delta_dipm_xyz,delta_qm_xyz)

        dipm = wfn%dpat(:,:,1)
        qm = wfn%qpat(:,:,1)

        delta_dipm = delta_dipm_xyz(:,:,1)
        delta_qm = delta_qm_xyz(:,:,1)

    end subroutine

    subroutine compute_traceless_mol_qm(n,xyz,q,dipm,qp,mol_qm)
        integer :: n,i,l,k,j
        real(wp) :: xyz(3,n),dipm(3,n), qp(6,n),q(n)
        real(wp) :: mol_qm(6)
        real(wp) :: tma(6),tmb(6),tmc(6),dum
        tma = 0.0_wp
        tmb = 0.0_wp
        tmc = 0.0_wp
        do i = 1,n
            l = 0
            do j = 1,3
                do k = 1,j
                    l = lin(k,j)
                    tma(l) = tma(l)+xyz(j,i)*xyz(k,i)*q(i)
                    tmb(l) = tmb(l)+dipm(k,i)*xyz(j,i)+dipm(j,i)*xyz(k,i)
                    tmc(l) = tmc(l)+qp(l,i)
                enddo
            enddo
        enddo
        ! remove traces and multiply with 3/2 in q and dip parts
        dum = tma(1)+tma(3)+tma(6)
        dum = 0.50_wp*dum
        tma = 1.50_wp*tma
        l = 0
        do j = 1,3
            l = l+j
            tma(l) = tma(l)-dum
        enddo
        dum = tmb(1)+tmb(3)+tmb(6)
        dum = 0.50_wp*dum
        tmb = 1.50_wp*tmb
        l = 0
        do j = 1,3
            l = l+j
            tmb(l) = tmb(l)-dum
        enddo
        mol_qm = tma+tmb+tmc
   end subroutine

   pure elemental integer function lin(i1,i2)
    !$acc routine seq
    integer,intent(in) :: i1,i2
    integer :: idum1,idum2
    idum1=max(i1,i2)
    idum2=min(i1,i2)
    lin=idum2+idum1*(idum1-1)/2
    return
  end function lin

  subroutine compute_mult_ext_mult_high_a(mol,wfn,calc,partial_charge,delta_partial,dipm,delta_dipm,qm,delta_qm)
    use tblite_ncoord_exp, only:  new_exp_ncoord, exp_ncoord_type

    type(context_type) :: ctx
    type(structure_type) :: mol
    type(xtb_calculator) :: calc
    type(wavefunction_type) :: wfn
    type(exp_ncoord_type) :: ncoord_exp
    real(wp) :: dipm_shell_tmp(3,calc%bas%nsh,1),partial_charge(mol%nat),a(1)
    real(wp) :: energy, z(mol%nat), cutoff,cn(mol%nat),delta_partial(mol%nat),delta_dipm_xyz(3,mol%nat,1),&
    &delta_qm_xyz(6,mol%nat,1)
    real(wp) :: dipm(3,mol%nat),delta_dipm(3,mol%nat)
    real(wp) :: qm(6,mol%nat),delta_qm(6,mol%nat) 
    integer :: n
   
    call new_gfn2_calculator(calc,mol)
    call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
    calc%save_integrals = .true.
    call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)
    
    !call mulliken_shellwise(calc%bas%nao,calc%bas%nsh,calc%bas%ao2sh,wfn%density(:,:,wfn%nspin),&
    !res%overlap,mulliken_shell)
    
    call mol_set_nuclear_charge(mol%nat,mol%num,mol%id,z)

    partial_charge =  wfn%qat(:,1) 
    
    call get_rcov(mol)
    
    a = 1000.0_wp
    n = size(a)
    call populate_inv_cn_array(mol%nat,mol%id,mol%xyz,a)
    call compute_cn(mol%nat,cn)
   
    !delta partial charge
    call get_delta_partial(mol%nat,n,partial_charge,mol%id,mol%xyz,&
    cn,delta_partial)
    call get_delta_mm(mol%nat,n,partial_charge,wfn%dpat,wfn%qpat,mol%id,mol%xyz,&
    cn,delta_dipm_xyz,delta_qm_xyz)

    dipm = wfn%dpat(:,:,1)
    qm = wfn%qpat(:,:,1)

    delta_dipm = delta_dipm_xyz(:,:,1)
    delta_qm = delta_qm_xyz(:,:,1)

    end subroutine compute_mult_ext_mult_high_a

  

end module