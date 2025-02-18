module tblite_test_xtbml
    use mctc_env, only : wp
    use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
       & test_failed
    use mctc_io_structure, only : new
    use mctc_io, only : structure_type
    use mstore, only : get_structure
    use tblite_context_type, only : context_type
    use tblite_basis_type
    use tblite_wavefunction , only : wavefunction_type, new_wavefunction
    use tblite_integral_type, only : integral_type, new_integral
    use tblite_xtb_h0, only : get_selfenergy, get_hamiltonian, get_occupation, &
       & get_hamiltonian_gradient
    use tblite_wavefunction_mulliken, only : get_mulliken_shell_multipoles
    use tblite_xtb_calculator, only : xtb_calculator
    use tblite_xtb_gfn2, only : new_gfn2_calculator
    use tblite_xtb_gfn1, only : new_gfn1_calculator
    use tblite_xtb_singlepoint, only : xtb_singlepoint
    use mctc_io_convert, only : aatoau
    use tblite_post_processing_list, only : add_post_processing, post_processing_list
    use tblite_param_post_processing, only : post_processing_param_list
    use tblite_param_xtbml_features, only : xtbml_features_record
    use tblite_results, only : results_type
    use tblite_param_serde, only : serde_record
    use tblite_blas, only : gemm
    use tblite_double_dictionary, only : double_dictionary_type
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
      new_unittest("xtbml-energy-sum-up-gfn2", test_energy_sum_up_gfn2),&
      new_unittest("xtbml-energy-sum-up-gfn1", test_energy_sum_up_gfn1),& 
      new_unittest("xtbml-rot", test_rotation_co2),&
      new_unittest("xtbml-translation", test_translation_co2)&
      ]
 
end subroutine collect_xtbml
 
 subroutine test_mulliken_charges_shell_h2p(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error
    type(context_type) :: ctx
    type(structure_type) :: mol
    type(xtb_calculator) :: calc
    type(wavefunction_type) :: wfn
    type(post_processing_list), allocatable :: pproc
    type(xtbml_features_record), allocatable :: xtbml_param
    type(post_processing_param_list), allocatable :: pparam
    type(results_type) :: res
    class(serde_record), allocatable :: tmp_record
    real(wp), allocatable :: mulliken_shell(:)
    real(wp) :: energy
    real(wp), parameter :: xyz(3, 2) = reshape((/0.0_wp,0.0_wp,0.35_wp,&
       &0.0_wp,0.0_wp,-0.35_wp/),shape=(/3,2/))
    integer, parameter :: num(2) = (/1,1/)
 
    call new(mol, num, xyz*aatoau, uhf=1, charge=1.0_wp)
    call new_gfn2_calculator(calc, mol, error)
    call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
    calc%save_integrals = .true.
    allocate(pproc)
    allocate(xtbml_param)
    allocate(pparam)
    xtbml_param%xtbml_density = .true.
    call move_alloc(xtbml_param, tmp_record)
    call pparam%push(tmp_record)
    call add_post_processing(pproc, pparam)
    call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc)
 
    call res%dict%get_entry("q_A", mulliken_shell)
    if (sum(mulliken_shell) - 1.0_wp > thr) then
       call test_failed(error, "Charge is not summing up to 1 electron for H2+")
       print'(3es21.14)', mulliken_shell
    end if
 
 end subroutine test_mulliken_charges_shell_h2p
 
 
subroutine test_dipm_shell_h2p(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(results_type) :: res
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   class(serde_record), allocatable :: tmp_record
   integer, parameter :: nsh = 2, nat=2
   real(wp) :: mol_dipm(3), mol_dipm_delta(3), energy
   real(wp) :: dipm_xyz(3,nat),qm_xyz(6,nat)
   real(wp), allocatable :: tmp_array(:)
   real(wp), allocatable :: partial(:), delta_partial(:)
   real(wp) :: delta_dipm_xyz(3,nat),delta_qm_xyz(6,nat)
 
   real(wp), parameter :: xyz(3, nat) = reshape((/0.0_wp,0.0_wp,0.35_wp,&
      &0.0_wp,0.0_wp,-0.35_wp/),shape=(/3,2/))
   integer, parameter :: num(nat) = (/1,1/)
   integer :: i
 
 
   call new(mol,num,xyz*aatoau,uhf=1,charge=1.0_wp)
   call new_gfn2_calculator(calc,mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_density = .true.
   xtbml_param%xtbml_tensor = .true.
   xtbml_param%xtbml_convolution = .true.
   xtbml_param%xtbml_a = [1.0_wp]
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc)   
 
   call res%dict%get_entry("dipm_A_x", tmp_array)
   dipm_xyz(1, :) = tmp_array
   call res%dict%get_entry("dipm_A_y", tmp_array)
   dipm_xyz(2, :) = tmp_array
   call res%dict%get_entry("dipm_A_z", tmp_array)
   dipm_xyz(3, :) = tmp_array
   call res%dict%get_entry("q_A", partial)
   mol_dipm = 0.0_wp
   
   call res%dict%get_entry("delta_dipm_A_x", tmp_array)
   delta_dipm_xyz(1, :) = tmp_array
   call res%dict%get_entry("delta_dipm_A_y", tmp_array)
   delta_dipm_xyz(2, :) = tmp_array
   call res%dict%get_entry("delta_dipm_A_z", tmp_array)
   delta_dipm_xyz(3, :) = tmp_array
   call res%dict%get_entry("delta_q_A", delta_partial)
   
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
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(integral_type) :: ints
   type(results_type) :: res
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   class(serde_record), allocatable :: tmp_record
   integer, parameter :: nat=3
   real(wp) :: mol_dipm(3), mol_dipm_delta(3), energy
   real(wp) :: dipm_xyz(3,nat),qm_xyz(6,nat)
   real(wp), allocatable :: tmp_array(:)
   real(wp), allocatable :: partial(:), delta_partial(:)
   real(wp) :: delta_dipm_xyz(3,nat),delta_qm_xyz(6,nat)
   integer :: i,j
 
   real(wp), parameter :: xyz(3, nat) = reshape((/0.00000,0.00000,1.0000000,&
      &0.0000,0.00000,-1.0000000,&
      &0.000000,0.000000,0.000000/),shape=(/3,nat/))
   integer, parameter :: num(nat) = (/8,8,6/)
   mol_dipm = 0.0_wp
   mol_dipm_delta = 0.0_wp
 
   call new(mol,num,xyz*aatoau,uhf=0,charge=0.0_wp)
   call new_gfn2_calculator(calc,mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   calc%save_integrals = .true.
   
   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_density = .true.
   xtbml_param%xtbml_tensor = .true.
   xtbml_param%xtbml_convolution = .true.
   xtbml_param%xtbml_a = [1.0_wp]
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc) 
   
   call res%dict%get_entry("dipm_A_x", tmp_array)
   dipm_xyz(1, :) = tmp_array
   call res%dict%get_entry("dipm_A_y", tmp_array)
   dipm_xyz(2, :) = tmp_array
   call res%dict%get_entry("dipm_A_z", tmp_array)
   dipm_xyz(3, :) = tmp_array
   call res%dict%get_entry("q_A", partial)
   mol_dipm = 0.0_wp
   
   call res%dict%get_entry("delta_dipm_A_x", tmp_array)
   delta_dipm_xyz(1, :) = tmp_array
   call res%dict%get_entry("delta_dipm_A_y", tmp_array)
   delta_dipm_xyz(2, :) = tmp_array
   call res%dict%get_entry("delta_dipm_A_z", tmp_array)
   delta_dipm_xyz(3, :) = tmp_array
   call res%dict%get_entry("delta_q_A", delta_partial)

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
 
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(integral_type) :: ints
   type(results_type) :: res
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   class(serde_record), allocatable :: tmp_record 
   integer, parameter :: nat=12
   real(wp) :: mol_dipm(3), mol_dipm_delta(3), energy
   real(wp) :: dipm_xyz(3,nat),qm_xyz(6,nat)
   real(wp), allocatable :: tmp_array(:)
   real(wp), allocatable :: partial(:), delta_partial(:)
   real(wp) :: delta_dipm_xyz(3,nat),delta_qm_xyz(6,nat),mol_qm(6),delta_mol_qm(6)
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
   mol_qm = 0.0_wp
   delta_mol_qm = 0.0_wp
   delta_qm_xyz = 0.0_wp
   delta_dipm_xyz = 0.0_wp

   dipm_xyz = 0.0_wp
   qm_xyz = 0.0_wp


   call new(mol,num,xyz*aatoau,uhf=0,charge=0.0_wp)
   call new_gfn2_calculator(calc,mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_density = .true.
   xtbml_param%xtbml_tensor = .true.
   xtbml_param%xtbml_convolution = .true.
   xtbml_param%xtbml_a = [1.0_wp]
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc) 

   call res%dict%get_entry("dipm_A_x", tmp_array)
   dipm_xyz(1, :) = tmp_array
   call res%dict%get_entry("dipm_A_y", tmp_array)
   dipm_xyz(2, :) = tmp_array
   call res%dict%get_entry("dipm_A_z", tmp_array)
   dipm_xyz(3, :) = tmp_array
   call res%dict%get_entry("q_A", partial)

   
   call res%dict%get_entry("delta_dipm_A_x", tmp_array)
   delta_dipm_xyz(1, :) = tmp_array
   call res%dict%get_entry("delta_dipm_A_y", tmp_array)
   delta_dipm_xyz(2, :) = tmp_array
   call res%dict%get_entry("delta_dipm_A_z", tmp_array)
   delta_dipm_xyz(3, :) = tmp_array
   call res%dict%get_entry("delta_q_A", delta_partial)
   
   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_xx", tmp_array)
   qm_xyz(1, :) = tmp_array
   deallocate(tmp_array) 
   call res%dict%get_entry("qm_A_xy", tmp_array)
   qm_xyz(2, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_yy", tmp_array)
   qm_xyz(3, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_xz", tmp_array)
   qm_xyz(4, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_yz", tmp_array)
   qm_xyz(5, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_zz", tmp_array)
   qm_xyz(6, :) = tmp_array

   deallocate(tmp_array)
   call res%dict%get_entry("delta_qm_A_xx", tmp_array)
   delta_qm_xyz(1, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("delta_qm_A_xy", tmp_array)
   delta_qm_xyz(2, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("delta_qm_A_yy", tmp_array)
   delta_qm_xyz(3, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("delta_qm_A_xz", tmp_array)
   delta_qm_xyz(4, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("delta_qm_A_yz", tmp_array)
   delta_qm_xyz(5, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("delta_qm_A_zz", tmp_array)
   delta_qm_xyz(6, :) = tmp_array

   do i= 1,nat
      do j = 1,3
         mol_dipm(j) = mol_dipm(j)+dipm_xyz(j,i)+ mol%xyz(j,i)*partial(i)
      enddo
   enddo
   mol_dipm_delta = 0.0_wp
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
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(integral_type) :: ints
   type(results_type) :: res
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   class(serde_record), allocatable :: tmp_record
   integer, parameter :: nat=12
   real(wp) :: mol_dipm(3), mol_dipm_delta(3), energy
   real(wp) :: dipm_xyz(3,nat),qm_xyz(6,nat)
   real(wp), allocatable :: tmp_array(:)
   real(wp), allocatable :: partial(:), delta_partial(:)
   real(wp) :: delta_dipm_xyz(3,nat),delta_qm_xyz(6,nat),mol_qm(6),delta_mol_qm(6)
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
   mol_qm = 0.0_wp
   delta_mol_qm = 0.0_wp
   delta_qm_xyz = 0.0_wp
   delta_dipm_xyz = 0.0_wp

   dipm_xyz = 0.0_wp
   qm_xyz = 0.0_wp
   
   call new(mol,num,xyz*aatoau,uhf=0,charge=0.0_wp)
   call new_gfn2_calculator(calc,mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_density = .true.
   xtbml_param%xtbml_tensor = .true.
   xtbml_param%xtbml_convolution = .true.
   xtbml_param%xtbml_a = [1000.0_wp]
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc) 
   
   call res%dict%get_entry("dipm_A_x", tmp_array)
   dipm_xyz(1, :) = tmp_array
   call res%dict%get_entry("dipm_A_y", tmp_array)
   dipm_xyz(2, :) = tmp_array
   call res%dict%get_entry("dipm_A_z", tmp_array)
   dipm_xyz(3, :) = tmp_array
   call res%dict%get_entry("q_A", partial)
   
   deallocate(tmp_array)
   call res%dict%get_entry("delta_dipm_A_x_1000.00", tmp_array)
   delta_dipm_xyz(1, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("delta_dipm_A_y_1000.00", tmp_array)
   delta_dipm_xyz(2, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("delta_dipm_A_z_1000.00", tmp_array)
   delta_dipm_xyz(3, :) = tmp_array
   call res%dict%get_entry("delta_q_A_1000.00", delta_partial)

   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_xx", tmp_array)
   qm_xyz(1, :) = tmp_array
   deallocate(tmp_array) 
   call res%dict%get_entry("qm_A_xy", tmp_array)
   qm_xyz(2, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_yy", tmp_array)
   qm_xyz(3, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_xz", tmp_array)
   qm_xyz(4, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_yz", tmp_array)
   qm_xyz(5, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("qm_A_zz", tmp_array)
   qm_xyz(6, :) = tmp_array

   deallocate(tmp_array)
   call res%dict%get_entry("delta_qm_A_xx_1000.00", tmp_array)
   delta_qm_xyz(1, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("delta_qm_A_xy_1000.00", tmp_array)
   delta_qm_xyz(2, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("delta_qm_A_yy_1000.00", tmp_array)
   delta_qm_xyz(3, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("delta_qm_A_xz_1000.00", tmp_array)
   delta_qm_xyz(4, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("delta_qm_A_yz_1000.00", tmp_array)
   delta_qm_xyz(5, :) = tmp_array
   deallocate(tmp_array)
   call res%dict%get_entry("delta_qm_A_zz_1000.00", tmp_array)
   delta_qm_xyz(6, :) = tmp_array
  
   do i= 1,nat
      do j = 1,3
         mol_dipm(j) = mol_dipm(j)+dipm_xyz(j,i)+ mol%xyz(j,i)*partial(i)
      enddo
   enddo
   mol_dipm_delta = sum(delta_dipm_xyz,dim=2)
   if (norm2(mol_dipm_delta-mol_dipm)  > thr2) then
      call test_failed(error, "Molecular dipole moment of extended and non-extended are not equal")
      print'(3es21.14)', norm2(mol_dipm_delta-mol_dipm)
   end if
 
   call compute_traceless_mol_qm(mol%nat,mol%xyz,partial,dipm_xyz,qm_xyz,mol_qm)
   
   delta_mol_qm = sum(delta_qm_xyz, dim=2)
   
   if (norm2(mol_qm-delta_mol_qm)  > thr2) then
      call test_failed(error, "Molecular quadrupole moment of extended and non-extended are not equal")
      print'(3es21.14)', norm2(mol_qm-delta_mol_qm)
   end if
end subroutine test_qp_shell_benz_high_a
 
 
subroutine test_rotation_co2(error)
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   class(serde_record), allocatable :: tmp_record
   !> Error handling
   integer,parameter :: nat = 3
   type(error_type), allocatable, intent(out) :: error
   real(wp), parameter :: xyz(3, nat) = reshape((/0.00000,0.00000,1.0000000,&
      &0.0000,0.00000,-1.0000000,&
      &0.000000,0.000000,0.000000/),shape=(/3,nat/))
 
   integer, parameter :: num(nat) = (/8,8,6/)
   type(results_type) :: res, res_
   real(wp) :: rot_matrix(3,3),xyz_rot(3,nat),energy, xyz_trans(3,nat)
   real(wp), allocatable :: xtbml(:,:),xtbml_rot(:,:),xtbml_trans(:,:)
   integer :: i
 
   call new(mol,num,xyz*aatoau,uhf=0,charge=0.0_wp)
   call new_gfn2_calculator(calc,mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_energy = .true.
   xtbml_param%xtbml_geometry = .true.
   xtbml_param%xtbml_density = .true.
   xtbml_param%xtbml_tensor = .false.
   xtbml_param%xtbml_convolution = .true.
   xtbml_param%xtbml_a = [1.0_wp]
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc) 
   
 
   rot_matrix = (reshape((/&
      &1.0_wp,0.0_wp,0.0_wp,&
      &0.0_wp,0.52532_wp,-0.85090352453_wp,&
      &0.0_wp,0.85090352453_wp,0.52532_wp/),shape=(/3,3/)))
   call gemm(rot_matrix,xyz,xyz_rot)
 
   call new(mol,num,xyz_rot*aatoau,uhf=0,charge=0.0_wp)
   call new_gfn2_calculator(calc,mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   deallocate(pproc)
   allocate(pproc)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res_, post_process=pproc) 
  
   if (.not.(compare_dict(res%dict, res_%dict, thr2))) then
      call test_failed(error, "Rotational invariance is not respected.")
   end if
 
end subroutine

subroutine test_translation_co2(error)
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   class(serde_record), allocatable :: tmp_record
   !> Error handling
   integer,parameter :: nat = 3
   type(error_type), allocatable, intent(out) :: error
   real(wp), parameter :: xyz(3, nat) = reshape((/0.00000,0.00000,1.0000000,&
      &0.0000,0.00000,-1.0000000,&
      &0.000000,0.000000,0.000000/),shape=(/3,nat/))
 
   integer, parameter :: num(nat) = (/8,8,6/)
   type(results_type) :: res, res_
   real(wp) :: rot_matrix(3,3),xyz_rot(3,nat),energy, xyz_trans(3,nat)
   real(wp), allocatable :: xtbml(:,:),xtbml_rot(:,:),xtbml_trans(:,:)
   integer :: i
 
   call new(mol,num,xyz*aatoau,uhf=0,charge=0.0_wp)
   call new_gfn2_calculator(calc,mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_energy = .true.
   xtbml_param%xtbml_geometry = .true.
   xtbml_param%xtbml_density = .true.
   xtbml_param%xtbml_tensor = .false.
   xtbml_param%xtbml_convolution = .true.
   xtbml_param%xtbml_a = [1.0_wp]
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc) 
   
   xyz_trans = xyz(:,:)
 
   do i = 1,nat
      xyz_trans(1,i) = xyz_trans(1,i) +5.0_wp
   enddo
 
   call new(mol,num,xyz_trans*aatoau,uhf=0,charge=0.0_wp)
   call new_gfn2_calculator(calc,mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   deallocate(pproc)
   allocate(pproc)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res_, post_process=pproc) 
  
   if (.not.(compare_dict(res%dict, res_%dict, thr2))) then
      call test_failed(error, "Translational invariance is not respected.")
   end if
 
end subroutine

function compare_dict(lhs, rhs, thr_) result(equal)
   type(double_dictionary_type), intent(in) :: lhs, rhs
   real(kind=wp), intent(in) :: thr_
   logical :: equal
   integer :: i
   real(wp), allocatable :: array1(:)
   real(wp), allocatable :: array2(:, :)
   real(wp), allocatable :: array3(:, :, :)
   real(wp), allocatable :: array1_(:)
   real(wp), allocatable :: array2_(:, :)
   real(wp), allocatable :: array3_(:, :, :)
   equal = .false.
   i = lhs%get_n_entries()
   if (i /= rhs%get_n_entries()) then
      return
   end if

   do i = 1, lhs%get_n_entries()
      
      call lhs%get_entry(i, array1)
      if (allocated(array1)) then
         call rhs%get_entry(i, array1_)
         if (allocated(array1_)) then
            
            if (any(array1-array1_  > thr_)) then
               return
            end if
            continue
         else 
            return
         end if
      end if

      call lhs%get_entry(i, array2)
      if (allocated(array2)) then
         call rhs%get_entry(i, array2_)
         if (allocated(array2_)) then
            if (any(array2-array2_  > thr_)) return
            continue
         else 
            return
         end if
      end if

      call lhs%get_entry(i, array3)
      if (allocated(array3)) then
         call rhs%get_entry(i, array3_)
         if (allocated(array3_)) then
            if (any(array3-array3_  > thr_)) return
            continue
         else 
            return
         end if
      end if
   end do

   equal = .true.

end function



subroutine test_energy_sum_up_gfn2(error)
 
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(integral_type) :: ints
   type(results_type) :: res
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   class(serde_record), allocatable :: tmp_record 
   integer, parameter :: nat=12
   real(wp), allocatable :: tmp_array(:)
   real(wp) :: energy, sum_energy
   integer :: i,j
   character(len=:), allocatable :: label1
 
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

   call new(mol,num,xyz*aatoau,uhf=0,charge=0.0_wp)
   call new_gfn2_calculator(calc, mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_energy = .true.
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc) 

   call res%dict%get_entry("E_tot", tmp_array) 
   sum_energy = sum_energy + sum(tmp_array)
   print'(3es21.14)', abs(sum_energy-energy)
   if (abs(sum_energy-energy)  > thr) then
      call test_failed(error, "GFN2: Energy features don't add up to total energy.")
      print'(3es21.14)', abs(sum_energy-energy)
   end if


   call res%dict%remove_entry("E_tot")
   call res%dict%remove_entry("w_tot")

   sum_energy = 0.0_wp
   do i = 1, res%dict%get_n_entries()
      call res%dict%get_entry(i, tmp_array)
      call res%dict%get_label(i, label1)
      write(*,*) label1
      write(*,*) tmp_array
      sum_energy = sum_energy + sum(tmp_array)
      deallocate(tmp_array)
   end do
   print'(3es21.14)', abs(sum_energy-energy)
   if (abs(sum_energy-energy)  > thr) then
      call test_failed(error, "GFN2: Energy features don't add up to total energy.")
      print'(3es21.14)', abs(sum_energy-energy)
   end if

end subroutine test_energy_sum_up_gfn2

subroutine test_energy_sum_up_gfn1(error)
 
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(integral_type) :: ints
   type(results_type) :: res
   type(post_processing_list), allocatable :: pproc
   type(xtbml_features_record), allocatable :: xtbml_param
   type(post_processing_param_list), allocatable :: pparam
   class(serde_record), allocatable :: tmp_record 
   integer, parameter :: nat=12
   real(wp), allocatable :: tmp_array(:)
   real(wp) :: energy, sum_energy
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

   call new(mol,num,xyz*aatoau,uhf=0,charge=0.0_wp)
   call new_gfn1_calculator(calc,mol,error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   allocate(pproc)
   allocate(xtbml_param)
   allocate(pparam)
   xtbml_param%xtbml_energy = .true.
   call move_alloc(xtbml_param, tmp_record)
   call pparam%push(tmp_record)
   call add_post_processing(pproc, pparam)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, results=res, post_process=pproc) 

   call res%dict%get_entry("E_tot", tmp_array) 
   sum_energy = sum_energy + sum(tmp_array)
   

   if (abs(sum_energy-energy)  > thr) then
      call test_failed(error, "GFN1: Energy features don't add up to total energy.")
      print'(3es21.14)', abs(sum_energy-energy)
   end if
   call res%dict%remove_entry("E_tot")
   call res%dict%remove_entry("w_tot")

   sum_energy = 0.0_wp
   do i = 1, res%dict%get_n_entries()
      call res%dict%get_entry(i, tmp_array) 
      sum_energy = sum_energy + sum(tmp_array)
      deallocate(tmp_array)
   end do
   
   
   if (abs(sum_energy-energy)  > thr) then
      call test_failed(error, "GFN1: Energy features don't add up to total energy.")
      print'(3es21.14)', abs(sum_energy-energy)
   end if

end subroutine test_energy_sum_up_gfn1
 
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

end module