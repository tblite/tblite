module test_post_processing
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_context_type, only : context_type
   use tblite_param_post_processing, only : post_processing_record_list
   use tblite_param_post_processing_molmom, only: molmom_record
   use tblite_post_processing_list, only : post_processing_list, add_post_processing
   use tblite_results, only : results_type
   use tblite_toml, only : toml_table, add_table, set_value, toml_key, get_value
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction, eeq_guess, &
      & get_density_matrix
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: collect_post_processing

   real(wp), parameter :: kt = 300.0_wp * 3.166808578545117e-06_wp
   real(wp), parameter :: acc = 0.01_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr1 = 1e5*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
contains

subroutine collect_post_processing(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [&
      new_unittest("check h2 wbo", test_h2_wbo), &
      new_unittest("test timer print", test_timer_print), &
      new_unittest("molmom param dump", test_molmom_dump_param), &
      new_unittest("post proc param load", test_pproc_load_param), &
      new_unittest("post proc param dump", test_pproc_dump_param), &
      new_unittest("molmom param dipm", test_molmom_dipm_param, should_fail=.true.), &
      new_unittest("molmom param qp", test_molmom_qp_param, should_fail=.true.), &
      new_unittest("check-m01-localization", test_m01_localization) &
   ]
end subroutine collect_post_processing


subroutine test_h2_wbo(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: n_atoms = 2
   type(structure_type) :: mol
   type(context_type) :: ctx
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(post_processing_list) :: pproc
   type(results_type) :: res
   real(kind=wp) :: energy
   real(kind=wp), allocatable :: wbo(:, :, :), wbo_exp(:, :, :)
   real(kind=wp), allocatable :: xyz(:, :)
   character(len=:), allocatable :: wbo_label
   integer, parameter :: atoms(2) =  [1, 1]
   allocate(xyz(3, n_atoms))
   xyz = reshape([&
      &+0.00000000_wp, +0.000000000_wp, +0.472429040_wp,&
      &+0.00000000_wp, +0.000000000_wp, -0.472429040_wp],&
      & shape(xyz))

   call new(mol, atoms, xyz, charge=+1.0_wp, uhf=1)
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   wbo_label = "bond-orders"
   call add_post_processing(pproc, mol, wbo_label, error)
   if (allocated(error)) return
   call eeq_guess(mol, calc, wfn, error)
   if (allocated(error)) return
   energy = 0.0_wp

   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, results=res, post_process=pproc, verbosity=0)
   call res%dict%get_entry("bond-orders", wbo)

   allocate(wbo_exp(2, 2, 2))
   wbo_exp = reshape([&
      &0.000000_wp, +0.50000000_wp, &
      &0.500000000_wp, 0.00000000_wp, &
      &0.000000_wp, +0.50000000_wp, &
      &0.500000000_wp, 0.00000000_wp],&
      & shape(wbo_exp))

   if (any(abs(wbo - wbo_exp) > thr)) then
         call test_failed(error, "Gradient of energy does not match")
         print"(3es21.14)", wbo
         print'("---")'
         print"(3es21.14)", wbo_exp
   end if
   if (allocated(error)) return

   deallocate(wbo, wbo_exp)
   mol%charge = 0.0_wp
   mol%uhf = 0
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, results=res, post_process=pproc, verbosity=0)
   call res%dict%get_entry("bond-orders", wbo)
   allocate(wbo_exp(2, 2, 1))
   wbo_exp = reshape([&
      &0.00000000_wp, +1.0000000_wp,&
      &+1.00000000_wp, 0.00000000_wp ],&
      & shape(wbo_exp))

   if (any(abs(wbo - wbo_exp) > thr)) then
      call test_failed(error, "Gradient of energy does not match")
      print"(3es21.14)", wbo
      print'("---")'
      print"(3es21.14)", wbo_exp
   end if
   if (allocated(error)) return

   deallocate(wbo, wbo_exp)
   mol%charge = -1.0_wp
   mol%uhf = 1
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, results=res, post_process=pproc, verbosity=0)
   call res%dict%get_entry("bond-orders", wbo)
   allocate(wbo_exp(2, 2, 2))
   wbo_exp = reshape([&
      &0.00000000_wp, +0.5000000_wp,&
      &+0.50000000_wp, 0.00000000_wp, &
      &0.00000000_wp, -0.5000000_wp,&
      &-0.50000000_wp, 0.00000000_wp ],&
      & shape(wbo_exp))

   if (any(abs(wbo - wbo_exp) > thr)) then
         call test_failed(error, "Gradient of energy does not match")
         print"(3es21.14)", wbo
         print'("---")'
         print"(3es21.14)", wbo_exp
   end if
   if (allocated(error)) return

   deallocate(wbo, wbo_exp)
   mol%charge = -2.0_wp
   mol%uhf = 0
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, results=res, post_process=pproc, verbosity=0)
   call res%dict%get_entry("bond-orders", wbo)
   allocate(wbo_exp(2, 2, 1))
   wbo_exp = reshape([&
      &0.00000000_wp, +0.0000000_wp,&
      &+0.00000000_wp, 0.00000000_wp ],&
      & shape(wbo_exp))

   if (any(abs(wbo - wbo_exp) > thr)) then
         call test_failed(error, "Gradient of energy does not match")
         print"(3es21.14)", wbo
         print'("---")'
         print"(3es21.14)", wbo_exp
   end if

end subroutine test_h2_wbo

subroutine test_timer_print(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: n_atoms = 2
   type(structure_type) :: mol
   type(context_type) :: ctx
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(post_processing_list) :: pproc
   type(results_type) :: res
   real(kind=wp) :: energy
   real(kind=wp), allocatable :: xyz(:, :)
   character(len=:), allocatable :: wbo_label, molmom_label, xtbml_label
   integer, parameter :: atoms(2) =  [1, 1]
   allocate(xyz(3, n_atoms))
   xyz = reshape([&
      &+0.00000000_wp, +0.000000000_wp, +0.472429040_wp,&
      &+0.00000000_wp, +0.000000000_wp, -0.472429040_wp],&
      & shape(xyz))

   call new(mol, atoms, xyz, charge=+1.0_wp, uhf=1)
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   wbo_label = "bond-orders"
   call add_post_processing(pproc, mol, wbo_label, error)
   if (allocated(error)) return
   molmom_label = "molmom"
   call add_post_processing(pproc, mol, molmom_label, error)
   if (allocated(error)) return
   xtbml_label = "xtbml"
   call add_post_processing(pproc, mol, xtbml_label, error)
   if (allocated(error)) return
   call eeq_guess(mol, calc, wfn, error)
   if (allocated(error)) return
   energy = 0.0_wp

   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, results=res, post_process=pproc, verbosity=3)

end subroutine test_timer_print

subroutine test_molmom_dipm_param(error)
   type(error_type), allocatable, intent(out) :: error
   type(toml_table), pointer :: table_multipole
   type(toml_table) :: table_post_proc
   type(molmom_record) :: param
   table_post_proc = toml_table()
   call add_table(table_post_proc, "molecular-multipole", table_multipole)
   call set_value(table_multipole, "dipole", 5)
   call set_value(table_multipole, "quadrupole", .true.)

   call param%load(table_post_proc, error)
   table_post_proc = toml_table()
   call add_table(table_post_proc, "molecular-multipole", table_multipole)
   call set_value(table_multipole, "dipole", .true.)
   call set_value(table_multipole, "quadrupole", 42)

   call param%load(table_post_proc, error)

end subroutine test_molmom_dipm_param

subroutine test_molmom_qp_param(error)
   type(error_type), allocatable, intent(out) :: error
   type(toml_table), pointer :: table_multipole
   type(toml_table) :: table_post_proc
   type(molmom_record) :: param

   table_post_proc = toml_table()
   call add_table(table_post_proc, "molecular-multipole", table_multipole)
   call set_value(table_multipole, "dipole", .true.)
   call set_value(table_multipole, "quadrupole", 42)

   call param%load(table_post_proc, error)

end subroutine test_molmom_qp_param

subroutine test_molmom_dump_param(error)
   type(error_type), allocatable, intent(out) :: error
   type(toml_table), pointer :: table_multipole, child
   type(toml_table) :: table_post_proc
   type(toml_table) :: new_table
   type(toml_key), allocatable :: list(:)
   type(molmom_record) :: param

   table_post_proc = toml_table()
   call add_table(table_post_proc, "molecular-multipole", table_multipole)
   call set_value(table_multipole, "dipole", .true.)
   call set_value(table_multipole, "quadrupole", .false.)

   call param%load(table_post_proc, error)
   call check(error, param%moldipm, .true.)
   if (allocated(error)) return
   call check(error, param%molqp, .false.)
   if (allocated(error)) return

   new_table = toml_table()
   call param%dump(new_table, error)
   call param%load(new_table, error)
   call new_table%get_keys(list)

   call check(error, size(list), 1)
   if (allocated(error)) return
   call get_value(new_table, list(1)%key, child)
   call child%get_keys(list)
   call check(error, size(list), 2)
   if (allocated(error)) return

end subroutine test_molmom_dump_param

subroutine test_pproc_load_param(error)
   type(error_type), allocatable, intent(out) :: error
   type(toml_table), pointer ::  table_entries
   type(toml_table) :: table_multipole
   type(post_processing_record_list) :: param

   table_multipole = toml_table()
   call add_table(table_multipole, "molecular-multipole", table_entries)
   call set_value(table_entries, "dipole", .true.)
   call set_value(table_entries, "quadrupole", .false.)
   call param%load(table_multipole, error)

end subroutine test_pproc_load_param

subroutine test_pproc_dump_param(error)
   type(error_type), allocatable, intent(out) :: error
   type(toml_table), pointer ::  table_entries, child
   type(toml_table) :: table_multipole
   type(toml_table) :: new_table
   type(toml_key), allocatable :: list(:)
   type(post_processing_record_list) :: param

   table_multipole = toml_table()
   call add_table(table_multipole, "molecular-multipole", table_entries)
   call set_value(table_entries, "dipole", .true.)
   call set_value(table_entries, "quadrupole", .false.)
   call param%load(table_multipole, error)
   new_table = toml_table()
   call param%dump(new_table, error)
   call param%load(new_table, error)
   call new_table%get_keys(list)

   call check(error, size(list), 1)
   if (allocated(error)) return
   call get_value(new_table, list(1)%key, child)
   call child%get_keys(list)
   call check(error, size(list), 2)
   if (allocated(error)) return
end subroutine test_pproc_dump_param


subroutine test_m01_localization(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(context_type) :: ctx
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(post_processing_list) :: pproc
   type(results_type) :: res
   real(wp) :: energy
   real(wp), allocatable :: lmo(:, :, :), pmat_can(:, :), pmat_loc(:, :), ortho(:, :)
   character(len=:), allocatable :: label
   integer :: nao, nocc, i

   call get_structure(mol, "MB16-43", "01")

   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return
   calc%save_integrals = .true.
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)

   label = "lmo-foster-boys"
   call add_post_processing(pproc, mol, label, error)
   if (allocated(error)) return

   call eeq_guess(mol, calc, wfn, error)
   if (allocated(error)) return
   energy = 0.0_wp
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, results=res, &
      & post_process=pproc, verbosity=0)

   call res%dict%get_entry("localized-orbitals", lmo)
   if (.not.allocated(lmo)) then
      call test_failed(error, "No localized orbitals in the result dictionary")
      return
   end if

   nao = calc%bas%nao
   nocc = nint(wfn%nel(1))

   call check(error, size(lmo, 1), nao)
   if (allocated(error)) return
   call check(error, size(lmo, 2), nao)
   if (allocated(error)) return
   call check(error, size(lmo, 3), wfn%nspin)
   if (allocated(error)) return

   ! The occupied block must actually have been rotated, not just copied
   if (maxval(abs(lmo(:, :nocc, 1) - wfn%coeff(:, :nocc, 1))) < thr2) then
      call test_failed(error, "Localization did not modify the occupied orbitals")
   end if
   if (allocated(error)) return

   ! Virtual orbitals must remain identical to the canonical coefficients
   if (any(abs(lmo(:, nocc+1:, 1) - wfn%coeff(:, nocc+1:, 1)) > thr)) then
      call test_failed(error, "Virtual orbitals were changed by the localization")
   end if
   if (allocated(error)) return

   ! Localized occupied orbitals must remain orthonormal with respect to the overlap
   allocate(ortho(nocc, nocc))
   ortho = matmul(transpose(lmo(:, :nocc, 1)), matmul(res%overlap, lmo(:, :nocc, 1)))
   do i = 1, nocc
      ortho(i, i) = ortho(i, i) - 1.0_wp
   end do
   if (maxval(abs(ortho)) > thr1) then
      call test_failed(error, "Localized occupied orbitals are not orthonormal")
   end if
   if (allocated(error)) return

   ! The occupied subspace, and hence the density matrix, must be unchanged.
   allocate(pmat_can(nao, nao), pmat_loc(nao, nao))
   pmat_can = wfn%density(:, :, 1)
   call get_density_matrix(wfn%focc(:, 1) + wfn%focc(:, 2), lmo(:, :, 1), pmat_loc)
   if (any(abs(pmat_can - pmat_loc) > thr1)) then
      call test_failed(error, "Localization changed the occupied-space density matrix")
      print "(3es21.14)", pmat_can
      print '("---")'
      print "(3es21.14)", pmat_loc
   end if
   if (allocated(error)) return

end subroutine test_m01_localization

end module test_post_processing
