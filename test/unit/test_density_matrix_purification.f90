module test_density_matrix_purification
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
   & test_failed
   use tblite_double_dictionary, only : double_dictionary_type
   use mstore, only : get_structure
   use tblite_context_type, only : context_type
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction, eeq_guess
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   use tblite_lapack_solver, only : lapack_solver
   use tblite_purification_solver_context, only : purification_solver_context, purification_type, &
      purification_runmode, purification_precision, dmp_input
   use mctc_io, only : structure_type

   implicit none
   private

   public :: collect_DMP

   real(wp), parameter :: thr = 5e-6_wp
   real(wp), parameter :: thr2 = 10*sqrt(epsilon(1.0_wp))
   real(wp), parameter :: etemp = 0.0_wp
   real(wp), parameter :: kt = etemp * 3.166808578545117e-06_wp
   real(wp), parameter :: acc = 1.0_wp

contains

   subroutine collect_DMP(testsuite)

      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("SP2 GFN2 energy", test_sp2_gfn2_energy), &
         new_unittest("SP2 GFN2 gradient", test_sp2_gfn2_egrad), &
         new_unittest("SP2 GFN1 energy", test_sp2_gfn1_energy), &
         new_unittest("SP2 GFN1 gradient", test_sp2_gfn1_egrad), &
         new_unittest("TRS4 GFN2 energy", test_trs4_gfn2_energy), &
         new_unittest("TRS4 GFN2 gradient", test_trs4_gfn2_egrad), &
         new_unittest("TRS4 GFN1 energy", test_trs4_gfn1_energy), &
         new_unittest("TRS4 GFN1 gradient", test_trs4_gfn1_egrad), &
         new_unittest("SP2-accel GFN2 energy", test_sp2_accel_gfn2_energy), &
         new_unittest("SP2-accel GFN2 gradient", test_sp2_accel_gfn2_egrad), &
         new_unittest("SP2-accel GFN1 energy", test_sp2_accel_gfn1_energy), &
         new_unittest("SP2-accel GFN1 gradient", test_sp2_accel_gfn1_egrad), &
         new_unittest("DMP failover", test_sp2_failover) &
         ]
   end subroutine collect_DMP

   subroutine test_sp2_gfn2_energy(error)
      type(error_type),allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn,  wfn1
      type(context_type) :: ctx, ctx1
      type(dmp_input) :: inp

      real(wp) :: gvd_energy = 0.0_wp
      real(wp) :: dmp_energy = 0.0_wp

      call get_structure(mol, "MB16-43", "01")
      call new_gfn2_calculator(calc, mol, error)
      ctx%solver = lapack_solver(1)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, gvd_energy, verbosity=0)
      inp%type = purification_type%tc2
      inp%runmode = purification_runmode%cpu
      inp%precision = purification_precision%double
      ctx1%solver = purification_solver_context(inp)
      call new_wavefunction(wfn1, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx1, mol, calc, wfn1, acc, dmp_energy, verbosity=0)

      call check(error, abs(gvd_energy - dmp_energy) < thr, "SP2 GFN2 energy")

   end subroutine

   subroutine test_sp2_gfn1_energy(error)
      type(error_type),allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn, wfn1
      type(context_type) :: ctx, ctx1
      type(dmp_input) :: inp

      real(wp) :: gvd_energy = 0.0_wp
      real(wp) :: dmp_energy = 0.0_wp

      call get_structure(mol, "MB16-43", "01")
      call new_gfn1_calculator(calc, mol, error)
      ctx%solver = lapack_solver(1)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, gvd_energy, verbosity=0)

      inp%type = purification_type%tc2
      inp%runmode = purification_runmode%cpu
      inp%precision = purification_precision%double
      ctx1%solver = purification_solver_context(inp)
      call new_wavefunction(wfn1, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx1, mol, calc, wfn1, acc, dmp_energy, verbosity=0)

      call check(error, abs(gvd_energy - dmp_energy) < thr, "SP2 GFN2 energy")

   end subroutine

   subroutine test_sp2_gfn2_egrad(error)
      type(error_type),allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn, wfn1
      type(context_type) :: ctx, ctx1
      type(dmp_input) :: inp

      real(wp), allocatable :: gvd_grad(:,:), dmp_grad(:,:), sigma(:, :)
      real(wp) :: gvd_energy = 0.0_wp
      real(wp) :: dmp_energy = 0.0_wp

      call get_structure(mol, "MB16-43", "01")
      allocate(gvd_grad(3, mol%nat), source=0.0_wp)
      allocate(dmp_grad(3, mol%nat), source=1.0_wp)
      allocate(sigma(3, 3), source=0.0_wp)
      call new_gfn2_calculator(calc, mol, error)
      ctx%solver = lapack_solver(1)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, gvd_energy, gvd_grad, sigma, verbosity=0)

      inp%type = purification_type%tc2
      inp%runmode = purification_runmode%cpu
      inp%precision = purification_precision%double
      ctx1%solver = purification_solver_context(inp)
      
      call new_wavefunction(wfn1, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx1, mol, calc, wfn1, acc, dmp_energy, dmp_grad, sigma, verbosity=0)

      call check(error, abs(gvd_energy - dmp_energy) < thr, "SP2 GFN2 energy")

      if (any(abs(dmp_grad - gvd_grad) > thr2)) then
         call test_failed(error, "Gradient of energy does not match")
         print'(3es21.14)', dmp_grad
         print'("---")'
         print'(3es21.14)', gvd_grad
         print'("---")'
         print'(3es21.14)', dmp_grad - gvd_grad
      end if

   end subroutine

   subroutine test_sp2_gfn1_egrad(error)
      type(error_type),allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn, wfn1
      type(context_type) :: ctx, ctx1
      type(dmp_input) :: inp

      real(wp), allocatable :: gvd_grad(:,:), dmp_grad(:,:), sigma(:, :)
      real(wp) :: gvd_energy = 0.0_wp
      real(wp) :: dmp_energy = 0.0_wp

      call get_structure(mol, "MB16-43", "01")
      allocate(gvd_grad(3, mol%nat), source=0.0_wp)
      allocate(dmp_grad(3, mol%nat), source=1.0_wp)
      allocate(sigma(3, 3), source=0.0_wp)
      call new_gfn1_calculator(calc, mol, error)
      ctx%solver = lapack_solver(1)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, gvd_energy, gvd_grad, sigma, verbosity=0)

      inp%type = purification_type%tc2
      inp%runmode = purification_runmode%cpu
      inp%precision = purification_precision%double
      ctx1%solver = purification_solver_context(inp)
      
      call new_wavefunction(wfn1, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx1, mol, calc, wfn1, acc, dmp_energy, dmp_grad, sigma, verbosity=0)

      call check(error, abs(gvd_energy - dmp_energy) < thr, "SP2 GFN2 energy")

      if (any(abs(dmp_grad - gvd_grad) > thr2)) then
         call test_failed(error, "Gradient of energy does not match")
         print'(3es21.14)', dmp_grad
         print'("---")'
         print'(3es21.14)', gvd_grad
         print'("---")'
         print'(3es21.14)', dmp_grad - gvd_grad
      end if

   end subroutine

   subroutine test_trs4_gfn2_energy(error)
      type(error_type),allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn, wfn1
      type(context_type) :: ctx, ctx1
      type(dmp_input) :: inp

      real(wp) :: gvd_energy = 0.0_wp
      real(wp) :: dmp_energy = 0.0_wp

      call get_structure(mol, "MB16-43", "01")
      call new_gfn2_calculator(calc, mol, error)
      ctx%solver = lapack_solver(1)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, gvd_energy, verbosity=0)

      inp%type = purification_type%trs4
      inp%runmode = purification_runmode%cpu
      inp%precision = purification_precision%double
      ctx1%solver = purification_solver_context(inp)
      call new_wavefunction(wfn1, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx1, mol, calc, wfn1, acc, dmp_energy, verbosity=0)

      call check(error, abs(gvd_energy - dmp_energy) < thr, "SP2 GFN2 energy")

   end subroutine

   subroutine test_trs4_gfn1_energy(error)
      type(error_type),allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn, wfn1
      type(context_type) :: ctx, ctx1
      type(dmp_input) :: inp

      real(wp) :: gvd_energy = 0.0_wp
      real(wp) :: dmp_energy = 0.0_wp

      call get_structure(mol, "MB16-43", "01")
      call new_gfn1_calculator(calc, mol, error)
      ctx%solver = lapack_solver(1)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, gvd_energy, verbosity=0)

      inp%type = purification_type%trs4
      inp%runmode = purification_runmode%cpu
      inp%precision = purification_precision%double
      ctx1%solver = purification_solver_context(inp)
      
      call new_wavefunction(wfn1, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx1, mol, calc, wfn1, acc, dmp_energy, verbosity=0)

      call check(error, abs(gvd_energy - dmp_energy) < thr, "SP2 GFN2 energy")

   end subroutine

   subroutine test_trs4_gfn2_egrad(error)
      type(error_type),allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn, wfn1
      type(context_type) :: ctx, ctx1
      type(dmp_input) :: inp

      real(wp), allocatable :: gvd_grad(:,:), dmp_grad(:,:), sigma(:, :)
      real(wp) :: gvd_energy = 0.0_wp
      real(wp) :: dmp_energy = 0.0_wp

      call get_structure(mol, "MB16-43", "01")
      allocate(gvd_grad(3, mol%nat), source=0.0_wp)
      allocate(dmp_grad(3, mol%nat), source=1.0_wp)
      allocate(sigma(3, 3), source=0.0_wp)
      call new_gfn2_calculator(calc, mol, error)
      ctx%solver = lapack_solver(1)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, gvd_energy, gvd_grad, sigma, verbosity=0)

      inp%type = purification_type%trs4
      inp%runmode = purification_runmode%cpu
      inp%precision = purification_precision%double
      ctx1%solver = purification_solver_context(inp)
      call new_wavefunction(wfn1, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx1, mol, calc, wfn1, acc, dmp_energy, dmp_grad, sigma, verbosity=0)

      call check(error, abs(gvd_energy - dmp_energy) < thr, "SP2 GFN2 energy")

      if (any(abs(dmp_grad - gvd_grad) > thr2)) then
         call test_failed(error, "Gradient of energy does not match")
         print'(3es21.14)', dmp_grad
         print'("---")'
         print'(3es21.14)', gvd_grad
         print'("---")'
         print'(3es21.14)', dmp_grad - gvd_grad
      end if

   end subroutine

   subroutine test_trs4_gfn1_egrad(error)
      type(error_type),allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn, wfn1
      type(context_type) :: ctx, ctx1
      type(dmp_input) :: inp

      real(wp), allocatable :: gvd_grad(:,:), dmp_grad(:,:), sigma(:, :)
      real(wp) :: gvd_energy = 0.0_wp
      real(wp) :: dmp_energy = 0.0_wp

      call get_structure(mol, "MB16-43", "01")
      allocate(gvd_grad(3, mol%nat), source=0.0_wp)
      allocate(dmp_grad(3, mol%nat), source=1.0_wp)
      allocate(sigma(3, 3), source=0.0_wp)
      call new_gfn1_calculator(calc, mol, error)
      ctx%solver = lapack_solver(1)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, gvd_energy, gvd_grad, sigma, verbosity=0)

      inp%type = purification_type%trs4
      inp%runmode = purification_runmode%cpu
      inp%precision = purification_precision%double
      ctx1%solver = purification_solver_context(inp)
      call new_wavefunction(wfn1, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx1, mol, calc, wfn1, acc, dmp_energy, dmp_grad, sigma, verbosity=0)

      call check(error, abs(gvd_energy - dmp_energy) < thr, "SP2 GFN2 energy")

      if (any(abs(dmp_grad - gvd_grad) > thr2)) then
         call test_failed(error, "Gradient of energy does not match")
         print'(3es21.14)', dmp_grad
         print'("---")'
         print'(3es21.14)', gvd_grad
         print'("---")'
         print'(3es21.14)', dmp_grad - gvd_grad
      end if

   end subroutine

   subroutine test_sp2_accel_gfn2_energy(error)
      type(error_type),allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn, wfn1
      type(context_type) :: ctx, ctx1
      type(dmp_input) :: inp

      real(wp) :: gvd_energy = 0.0_wp
      real(wp) :: dmp_energy = 0.0_wp

      call get_structure(mol, "MB16-43", "01")
      call new_gfn2_calculator(calc, mol, error)
      ctx%solver = lapack_solver(1)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, gvd_energy, verbosity=0)

      inp%type = purification_type%tc2accel
      inp%runmode = purification_runmode%cpu
      inp%precision = purification_precision%double
      ctx1%solver = purification_solver_context(inp)
      call new_wavefunction(wfn1, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx1, mol, calc, wfn1, acc, dmp_energy, verbosity=0)

      call check(error, abs(gvd_energy - dmp_energy) < thr, "SP2 GFN2 energy")

   end subroutine

   subroutine test_sp2_accel_gfn1_energy(error)
      type(error_type),allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn, wfn1
      type(context_type) :: ctx, ctx1
      type(dmp_input) :: inp

      real(wp) :: gvd_energy = 0.0_wp
      real(wp) :: dmp_energy = 0.0_wp

      call get_structure(mol, "MB16-43", "01")
      call new_gfn1_calculator(calc, mol, error)
      ctx%solver = lapack_solver(1)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, gvd_energy, verbosity=0)

      inp%type = purification_type%tc2accel
      inp%runmode = purification_runmode%cpu
      inp%precision = purification_precision%double
      ctx1%solver = purification_solver_context(inp)
      call new_wavefunction(wfn1, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx1, mol, calc, wfn1, acc, dmp_energy, verbosity=0)

      call check(error, abs(gvd_energy - dmp_energy) < thr, "SP2 GFN2 energy")

   end subroutine

   subroutine test_sp2_accel_gfn2_egrad(error)
      type(error_type),allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn, wfn1
      type(context_type) :: ctx, ctx1
      type(dmp_input) :: inp

      real(wp), allocatable :: gvd_grad(:,:), dmp_grad(:,:), sigma(:, :)
      real(wp) :: gvd_energy = 0.0_wp
      real(wp) :: dmp_energy = 0.0_wp

      call get_structure(mol, "MB16-43", "01")
      allocate(gvd_grad(3, mol%nat), source=0.0_wp)
      allocate(dmp_grad(3, mol%nat), source=1.0_wp)
      allocate(sigma(3, 3), source=0.0_wp)
      call new_gfn2_calculator(calc, mol, error)
      ctx%solver = lapack_solver(1)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, gvd_energy, gvd_grad, sigma, verbosity=0)

      inp%type = purification_type%tc2accel
      inp%runmode = purification_runmode%cpu
      inp%precision = purification_precision%double
      ctx1%solver = purification_solver_context(inp)
      call new_wavefunction(wfn1, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx1, mol, calc, wfn1, acc, dmp_energy, dmp_grad, sigma, verbosity=0)

      call check(error, abs(gvd_energy - dmp_energy) < thr, "SP2 GFN2 energy")

      if (any(abs(dmp_grad - gvd_grad) > thr2)) then
         call test_failed(error, "Gradient of energy does not match")
         print'(3es21.14)', dmp_grad
         print'("---")'
         print'(3es21.14)', gvd_grad
         print'("---")'
         print'(3es21.14)', dmp_grad - gvd_grad
      end if

   end subroutine

   subroutine test_sp2_accel_gfn1_egrad(error)
      type(error_type),allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn, wfn1
      type(context_type) :: ctx, ctx1
      type(dmp_input) :: inp

      real(wp), allocatable :: gvd_grad(:,:), dmp_grad(:,:), sigma(:, :)
      real(wp) :: gvd_energy = 0.0_wp
      real(wp) :: dmp_energy = 0.0_wp

      call get_structure(mol, "MB16-43", "01")
      allocate(gvd_grad(3, mol%nat), source=0.0_wp)
      allocate(dmp_grad(3, mol%nat), source=1.0_wp)
      allocate(sigma(3, 3), source=0.0_wp)
      call new_gfn1_calculator(calc, mol, error)
      ctx%solver = lapack_solver(1)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, gvd_energy, gvd_grad, sigma, verbosity=0)

      inp%type = purification_type%tc2accel
      inp%runmode = purification_runmode%cpu
      inp%precision = purification_precision%double
      ctx1%solver = purification_solver_context(inp)
      call new_wavefunction(wfn1, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx1, mol, calc, wfn1, acc, dmp_energy, dmp_grad, sigma, verbosity=0)

      call check(error, abs(gvd_energy - dmp_energy) < thr, "SP2 GFN2 energy")

      if (any(abs(dmp_grad - gvd_grad) > thr2)) then
         call test_failed(error, "Gradient of energy does not match")
         print'(3es21.14)', dmp_grad
         print'("---")'
         print'(3es21.14)', gvd_grad
         print'("---")'
         print'(3es21.14)', dmp_grad - gvd_grad
      end if

   end subroutine

   subroutine test_sp2_failover(error)
      type(error_type),allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn, wfn1
      type(context_type) :: ctx, ctx1
      type(dmp_input) :: inp

      real(wp) :: gvd_energy = 0.0_wp
      real(wp) :: dmp_energy = 0.0_wp

      call get_structure(mol, "IL16", "008A")
      !set false charge to trigger failover
      mol%charge = 0
      call new_gfn2_calculator(calc, mol, error)
      
      ctx%solver = lapack_solver(1)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, gvd_energy, verbosity=0)

      inp%type = purification_type%tc2accel
      inp%runmode = purification_runmode%cpu
      inp%precision = purification_precision%double
      ctx1%solver = purification_solver_context(inp)
      call new_wavefunction(wfn1, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      call xtb_singlepoint(ctx1, mol, calc, wfn1, acc, dmp_energy, verbosity=0)

      call check(error, abs(gvd_energy - dmp_energy) < thr, "SP2 GFN2 energy")


   end subroutine


end module test_density_matrix_purification

