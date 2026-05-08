module test_post_processing
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_blas, only : gemm
   use tblite_context_type, only : context_type
   use tblite_param_post_processing_molmom, only: molmom_record
   use tblite_param_post_processing, only : post_processing_record_list
   use tblite_post_processing_list, only : post_processing_list, add_post_processing
   use tblite_results, only : results_type
   use tblite_toml, only : toml_table, add_table, set_value, toml_key, get_value
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction, eeq_guess
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: collect_post_processing

   real(wp), parameter :: kt = 300.0_wp * 3.166808578545117e-06_wp
   real(wp), parameter :: acc = 0.01_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
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
      new_unittest("trafo-pcl", test_trafo_pcl) &
   ]
end subroutine


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
         print'(3es21.14)', wbo
         print'("---")'
         print'(3es21.14)', wbo_exp
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
      print'(3es21.14)', wbo
      print'("---")'
      print'(3es21.14)', wbo_exp
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
         print'(3es21.14)', wbo
         print'("---")'
         print'(3es21.14)', wbo_exp
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
         print'(3es21.14)', wbo
         print'("---")'
         print'(3es21.14)', wbo_exp
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
   character(len=:), allocatable :: wbo_label, molmom_label, trafo_label, xtbml_label
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
   trafo_label = "trafo"
   call add_post_processing(pproc, mol, trafo_label, error)
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


subroutine test_density_trafo(mol, calc, wfn, ref, error, thr_in)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Calculator instance
   type(xtb_calculator), intent(in) :: calc
   !> Wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Reference density matrix in Cartesian basis
   real(wp), intent(in) :: ref(:, :, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: spin
   real(wp) :: energy, thr_
   real(wp), allocatable :: coeff_cart(:, :, :), focc(:), density_cart(:, :, :)
   character(len=:), allocatable :: trafo_label
   type(context_type) :: ctx
   type(post_processing_list) :: pproc
   type(results_type) :: res
   
   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   ! Set up post-processing with spherical-to-cartesian transformation
   trafo_label = "trafo"
   call add_post_processing(pproc, mol, trafo_label, error)
   if (allocated(error)) return

   ! Perform single point calculation
   call eeq_guess(mol, calc, wfn, error)
   if (allocated(error)) return
   energy = 0.0_wp
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, &
      & results=res, post_process=pproc)

   ! Retrieve transformed MO coefficients
   call res%dict%get_entry("cartesian-mos", coeff_cart)
   
   ! Calculate density matrix in the cartesian basis
   allocate(focc(calc%bas%nao))
   allocate(density_cart(calc%bas%nao_cart, calc%bas%nao_cart, wfn%nspin), source=0.0_wp)
   do spin = 1, wfn%nspin
      if (wfn%nspin == 1) then
         focc = wfn%focc(:, 1) + wfn%focc(:, 2)
      else
         focc = wfn%focc(:, spin)
      end if

      call gemm(coeff_cart(:, :, spin) * spread(focc, 1, calc%bas%nao_cart), &
         & transpose(coeff_cart(:, :, spin)), density_cart(:, :, spin))
   end do

   if (any(abs(density_cart - ref) > thr_)) then
      call test_failed(error, "Cartesian density matrix does not match")
      print'(3es21.14)', density_cart
      print'("---")'
      print'(3es21.14)', ref
      print'("---")'
      print'(3es21.14)', density_cart-ref
      write(*,*) "Max abs. deviation:", maxval(abs(density_cart - ref))
   end if

end subroutine test_density_trafo

subroutine test_trafo_pcl(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: density_cart(20, 20, 1) = reshape([&
      &  1.96022219468921E+00_wp, -4.31628878279072E-16_wp, -3.69459754888647E-01_wp, &
      & -1.11098031904690E-16_wp,  4.71236654550718E-02_wp,  4.71236654550742E-02_wp, &
      & -9.42473309101460E-02_wp,  5.62965590903630E-16_wp, -4.06828497062294E-16_wp, &
      & -8.52830081768616E-16_wp, -1.89854934222313E-01_wp,  8.09624498998321E-16_wp, &
      &  7.25134761302635E-02_wp,  1.59916712486622E-16_wp, -2.34268009256941E-02_wp, &
      & -2.34268009256942E-02_wp,  4.68536018513883E-02_wp,  3.48317988599535E-16_wp, &
      &  1.88088687294335E-16_wp, -1.76463194084002E-16_wp, -4.31628878279072E-16_wp, &
      &  1.07791095451279E+00_wp, -8.31505858137497E-16_wp,  7.22200077518664E-14_wp, &
      & -6.20379008988276E-17_wp, -3.24781092909577E-16_wp,  3.86818993808405E-16_wp, &
      &  1.07556891827876E-16_wp, -8.61810622865278E-15_wp, -5.04737692105317E-02_wp, &
      &  1.13747411177022E-16_wp,  1.91867623619473E-01_wp, -1.02920847827749E-15_wp, &
      & -2.68951527715444E-14_wp, -9.11978402679842E-17_wp,  6.66905072726581E-17_wp, &
      &  2.45073329953261E-17_wp, -2.29032624584476E-16_wp, -4.62824223390612E-15_wp, &
      & -7.61574005556389E-02_wp, -3.69459754888647E-01_wp, -8.31505858137497E-16_wp, &
      &  5.06882651413875E-01_wp,  6.61199693299378E-16_wp,  5.03840047492308E-03_wp, &
      &  5.03840047492275E-03_wp, -1.00768009498458E-02_wp, -1.14256085434253E-16_wp, &
      &  3.31421482075486E-16_wp,  1.82617289715317E-16_wp, -1.39976031367795E-02_wp, &
      &  2.04803054630828E-18_wp, -7.33326590175914E-01_wp, -4.24957452481721E-16_wp, &
      & -1.19309309178673E-03_wp, -1.19309309178663E-03_wp,  2.38618618357336E-03_wp, &
      &  3.14132994583873E-18_wp, -1.88721574223030E-16_wp, -1.26943910123494E-17_wp, &
      & -1.11098031904690E-16_wp,  7.22200077518664E-14_wp,  6.61199693299378E-16_wp, &
      &  1.07791095451263E+00_wp,  1.42412783202672E-16_wp,  1.00685168645336E-16_wp, &
      & -2.43097951848008E-16_wp, -2.48680357460767E-16_wp, -5.04737692105127E-02_wp, &
      & -8.43769498715119E-15_wp, -4.77025776053928E-17_wp, -2.65620858641569E-14_wp, &
      &  5.25687101242358E-16_wp,  1.91867623619530E-01_wp, -1.45451846493814E-16_wp, &
      &  5.70778474664071E-17_wp,  8.83739990274073E-17_wp,  2.92526454109793E-16_wp, &
      & -7.61574005556284E-02_wp, -4.41660596983695E-15_wp,  4.71236654550718E-02_wp, &
      & -6.20379008988276E-17_wp,  5.03840047492308E-03_wp,  1.42412783202672E-16_wp, &
      &  1.97110953674444E-03_wp,  1.97110953674450E-03_wp, -3.94221907348894E-03_wp, &
      &  1.32961974553750E-17_wp, -1.55896748336360E-17_wp,  6.40735691835751E-18_wp, &
      & -3.38948815549140E-02_wp,  2.16655877069712E-16_wp, -2.24138588507156E-02_wp, &
      & -5.72903643784898E-17_wp, -8.09838349982601E-04_wp, -8.09838349982605E-04_wp, &
      &  1.61967669996521E-03_wp,  9.02129887592915E-18_wp, -8.57688494343901E-18_wp, &
      & -1.97876393374441E-17_wp,  4.71236654550742E-02_wp, -3.24781092909577E-16_wp, &
      &  5.03840047492275E-03_wp,  1.00685168645336E-16_wp,  1.97110953674450E-03_wp, &
      &  1.97110953674456E-03_wp, -3.94221907348906E-03_wp,  1.32961974553757E-17_wp, &
      & -2.40519587348739E-17_wp,  1.25681864548440E-17_wp, -3.38948815549142E-02_wp, &
      &  1.24259581954180E-16_wp, -2.24138588507157E-02_wp, -1.42095173197775E-16_wp, &
      & -8.09838349982632E-04_wp, -8.09838349982636E-04_wp,  1.61967669996527E-03_wp, &
      &  9.02129887592965E-18_wp, -4.90329933150083E-18_wp, -7.96338946182486E-19_wp, &
      & -9.42473309101460E-02_wp,  3.86818993808405E-16_wp, -1.00768009498458E-02_wp, &
      & -2.43097951848008E-16_wp, -3.94221907348894E-03_wp, -3.94221907348906E-03_wp, &
      &  7.88443814697800E-03_wp, -2.65923949107507E-17_wp,  3.96416335685099E-17_wp, &
      & -1.89755433732015E-17_wp,  6.77897631098282E-02_wp, -3.40915459023892E-16_wp, &
      &  4.48277177014313E-02_wp,  1.99385537576265E-16_wp,  1.61967669996523E-03_wp, &
      &  1.61967669996524E-03_wp, -3.23935339993047E-03_wp, -1.80425977518588E-17_wp, &
      &  1.34801842749398E-17_wp,  2.05839782836266E-17_wp,  5.62965590903630E-16_wp, &
      &  1.07556891827876E-16_wp, -1.14256085434253E-16_wp, -2.48680357460767E-16_wp, &
      &  1.32961974553750E-17_wp,  1.32961974553757E-17_wp, -2.65923949107507E-17_wp, &
      &  3.03868259254195E-31_wp,  2.95679184093213E-17_wp,  4.00370232499131E-17_wp, &
      & -5.51399477298318E-17_wp,  3.53975478135562E-16_wp,  3.41688307716610E-17_wp, &
      &  8.88793920882234E-17_wp, -6.62731219503248E-18_wp, -6.62731219503245E-18_wp, &
      &  1.32546243900649E-17_wp, -1.06242698042029E-31_wp,  1.63217263204427E-17_wp, &
      & -1.07382345176700E-17_wp, -4.06828497062294E-16_wp, -8.61810622865278E-15_wp, &
      &  3.31421482075486E-16_wp, -5.04737692105127E-02_wp, -1.55896748336360E-17_wp, &
      & -2.40519587348739E-17_wp,  3.96416335685099E-17_wp,  2.95679184093213E-17_wp, &
      &  3.41889525777604E-02_wp,  1.17007098454636E-15_wp, -9.57915379548323E-17_wp, &
      &  3.47118167542959E-15_wp, -6.67167431068387E-16_wp,  2.27433048909729E-01_wp, &
      & -2.28934973120123E-18_wp,  3.10411422805550E-17_wp, -2.87517925493538E-17_wp, &
      &  1.03182578415653E-17_wp,  1.34969454111339E-03_wp,  6.11490025281825E-16_wp, &
      & -8.52830081768616E-16_wp, -5.04737692105317E-02_wp,  1.82617289715317E-16_wp, &
      & -8.43769498715119E-15_wp,  6.40735691835751E-18_wp,  1.25681864548440E-17_wp, &
      & -1.89755433732015E-17_wp,  4.00370232499131E-17_wp,  1.17007098454636E-15_wp, &
      &  3.41889525777625E-02_wp,  1.92950814569396E-16_wp,  2.27433048909735E-01_wp, &
      &  1.23395282461191E-16_wp,  4.02802791121815E-15_wp, -2.27416693403459E-17_wp, &
      &  3.25960150188036E-17_wp, -9.85434567845777E-18_wp, -7.97202226222352E-17_wp, &
      &  5.55545193181572E-16_wp,  1.34969454111463E-03_wp, -1.89854934222313E-01_wp, &
      &  1.13747411177022E-16_wp, -1.39976031367795E-02_wp, -4.77025776053928E-17_wp, &
      & -3.38948815549140E-02_wp, -3.38948815549142E-02_wp,  6.77897631098282E-02_wp, &
      & -5.51399477298318E-17_wp, -9.57915379548323E-17_wp,  1.92950814569396E-16_wp, &
      &  1.97254097945181E+00_wp,  6.32122311057884E-17_wp,  1.62414000110393E-01_wp, &
      &  2.18292545659950E-18_wp,  7.69029197968299E-03_wp,  7.69029197968316E-03_wp, &
      & -1.53805839593662E-02_wp,  6.67419266931980E-17_wp,  3.54773941002931E-18_wp, &
      &  9.94892415398064E-16_wp,  8.09624498998321E-16_wp,  1.91867623619473E-01_wp, &
      &  2.04803054630830E-18_wp, -2.65620858641569E-14_wp,  2.16655877069712E-16_wp, &
      &  1.24259581954180E-16_wp, -3.40915459023892E-16_wp,  3.53975478135562E-16_wp, &
      &  3.47118167542959E-15_wp,  2.27433048909735E-01_wp,  6.32122311057884E-17_wp, &
      &  1.79039127284740E+00_wp, -1.58698189670615E-16_wp,  9.72139035937403E-15_wp, &
      & -3.12677515642898E-16_wp,  1.81426312876657E-16_wp,  1.31251202766241E-16_wp, &
      & -7.12641717691887E-16_wp,  1.87697080100691E-15_wp, -3.00207544958856E-02_wp, &
      &  7.25134761302635E-02_wp, -1.02920847827749E-15_wp, -7.33326590175914E-01_wp, &
      &  5.25687101242358E-16_wp, -2.24138588507156E-02_wp, -2.24138588507157E-02_wp, &
      &  4.48277177014313E-02_wp,  3.41688307716610E-17_wp, -6.67167431068387E-16_wp, &
      &  1.23395282461191E-16_wp,  1.62414000110393E-01_wp, -1.58698189670615E-16_wp, &
      &  1.19108922406209E+00_wp, -2.07339616780469E-16_wp,  8.57919007600030E-03_wp, &
      &  8.57919007600018E-03_wp, -1.71583801520005E-02_wp, -9.54742728903731E-17_wp, &
      &  1.58300091699251E-16_wp,  2.85677087660719E-16_wp,  1.59916712486622E-16_wp, &
      & -2.68951527715444E-14_wp, -4.24957452481721E-16_wp,  1.91867623619530E-01_wp, &
      & -5.72903643784898E-17_wp, -1.42095173197775E-16_wp,  1.99385537576265E-16_wp, &
      &  8.88793920882234E-17_wp,  2.27433048909729E-01_wp,  4.02629318774217E-15_wp, &
      &  2.18292545659951E-18_wp,  9.72139035937403E-15_wp, -2.07339616780469E-16_wp, &
      &  1.79039127284737E+00_wp, -9.67122129575804E-17_wp,  2.57384468356189E-16_wp, &
      & -1.60672255398608E-16_wp,  2.30473471113603E-16_wp, -3.00207544958891E-02_wp, &
      &  2.29937596740726E-15_wp, -2.34268009256941E-02_wp, -9.11978402679842E-17_wp, &
      & -1.19309309178673E-03_wp, -1.45451846493814E-16_wp, -8.09838349982601E-04_wp, &
      & -8.09838349982632E-04_wp,  1.61967669996523E-03_wp, -6.62731219503248E-18_wp, &
      & -2.28934973120123E-18_wp, -2.27416693403459E-17_wp,  7.69029197968299E-03_wp, &
      & -3.12677515642898E-16_wp,  8.57919007600030E-03_wp, -9.67122129575804E-17_wp, &
      &  3.63656149333574E-04_wp,  3.63656149333575E-04_wp, -7.27312298667149E-04_wp, &
      & -4.77930967387482E-18_wp,  1.02256863982123E-17_wp,  1.53607452191457E-17_wp, &
      & -2.34268009256942E-02_wp,  6.66905072726581E-17_wp, -1.19309309178663E-03_wp, &
      &  5.70778474664071E-17_wp, -8.09838349982605E-04_wp, -8.09838349982636E-04_wp, &
      &  1.61967669996524E-03_wp, -6.62731219503245E-18_wp,  3.10411422805550E-17_wp, &
      &  3.25960150188036E-17_wp,  7.69029197968316E-03_wp,  1.81426312876657E-16_wp, &
      &  8.57919007600018E-03_wp,  2.57384468356189E-16_wp,  3.63656149333575E-04_wp, &
      &  3.63656149333576E-04_wp, -7.27312298667151E-04_wp, -4.77930967387496E-18_wp, &
      & -7.06529069837249E-18_wp, -1.63261802835476E-19_wp,  4.68536018513883E-02_wp, &
      &  2.45073329953261E-17_wp,  2.38618618357336E-03_wp,  8.83739990274073E-17_wp, &
      &  1.61967669996521E-03_wp,  1.61967669996527E-03_wp, -3.23935339993047E-03_wp, &
      &  1.32546243900649E-17_wp, -2.87517925493538E-17_wp, -9.85434567845777E-18_wp, &
      & -1.53805839593662E-02_wp,  1.31251202766241E-16_wp, -1.71583801520005E-02_wp, &
      & -1.60672255398608E-16_wp, -7.27312298667149E-04_wp, -7.27312298667151E-04_wp, &
      &  1.45462459733430E-03_wp,  9.55861934774978E-18_wp, -3.16039569983986E-18_wp, &
      & -1.51974834163102E-17_wp,  3.48317988599535E-16_wp, -2.29032624584476E-16_wp, &
      &  3.14132994583873E-18_wp,  2.92526454109793E-16_wp,  9.02129887592915E-18_wp, &
      &  9.02129887592965E-18_wp, -1.80425977518588E-17_wp, -1.06242698042029E-31_wp, &
      &  1.03182578415653E-17_wp, -7.97202226222352E-17_wp,  6.67419266931980E-17_wp, &
      & -7.12641717691887E-16_wp, -9.54742728903731E-17_wp,  2.30473471113603E-16_wp, &
      & -4.77930967387482E-18_wp, -4.77930967387496E-18_wp,  9.55861934774978E-18_wp, &
      &  4.81946035529280E-31_wp, -2.23403459632380E-17_wp,  2.24806227067167E-17_wp, &
      &  1.88088687294335E-16_wp, -4.62824223390612E-15_wp, -1.88721574223030E-16_wp, &
      & -7.61574005556284E-02_wp, -8.57688494343901E-18_wp, -4.90329933150083E-18_wp, &
      &  1.34801842749398E-17_wp,  1.63217263204427E-17_wp,  1.34969454111339E-03_wp, &
      &  5.55545193181572E-16_wp,  3.54773941002931E-18_wp,  1.87697080100691E-15_wp, &
      &  1.58300091699251E-16_wp, -3.00207544958891E-02_wp,  1.02256863982123E-17_wp, &
      & -7.06529069837249E-18_wp, -3.16039569983986E-18_wp, -2.23403459632380E-17_wp, &
      &  5.53508929793042E-03_wp,  2.74953670942324E-16_wp, -1.76463194084002E-16_wp, &
      & -7.61574005556389E-02_wp, -1.26943910123494E-17_wp, -4.42007541678890E-15_wp, &
      & -1.97876393374441E-17_wp, -7.96338946182485E-19_wp,  2.05839782836266E-17_wp, &
      & -1.07382345176700E-17_wp,  6.11490025281825E-16_wp,  1.34969454111463E-03_wp, &
      &  9.94892415398064E-16_wp, -3.00207544958856E-02_wp,  2.85677087660719E-16_wp, &
      &  2.29937596740726E-15_wp,  1.53607452191457E-17_wp, -1.63261802835476E-19_wp, &
      & -1.51974834163102E-17_wp,  2.24806227067167E-17_wp,  2.74953670942324E-16_wp, &
      &  5.53508929793114E-03_wp], shape(density_cart))

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn

   call get_structure(mol, "MB16-43", "PCl")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   call test_density_trafo(mol, calc, wfn, density_cart, error, thr_in=thr*100.0_wp)

end subroutine test_trafo_pcl

end module test_post_processing