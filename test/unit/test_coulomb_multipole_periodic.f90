! This file is part of tblite.
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

module test_coulomb_multipole_periodic
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_container_cache, only : container_cache
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_multipole, only : damped_multipole, new_damped_multipole
   use tblite_cutoff, only : get_lattice_points
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: collect_coulomb_multipole_periodic

   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   abstract interface
      subroutine multipole_maker(multipole, mol)
         import :: damped_multipole, structure_type
         type(damped_multipole), intent(out) :: multipole
         type(structure_type), intent(in) :: mol
      end subroutine multipole_maker
   end interface

contains


!> Collect all exported unit tests
subroutine collect_coulomb_multipole_periodic(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("energy-co2-realspace", test_e_effective_co2_realspace), &
      new_unittest("energy-urea-realspace", test_e_effective_urea_realspace), &
      new_unittest("energy-co2-supercell", test_e_effective_co2), &
      new_unittest("energy-urea-supercell", test_e_effective_urea) &
      ]

end subroutine collect_coulomb_multipole_periodic


subroutine test_energy(error, mol, qat, dpat, qpat, make_multipole, energy)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Expected multipole moments
   real(wp), intent(in) :: qat(:)

   !> Expected dipole moments
   real(wp), intent(in) :: dpat(:, :)

   !> Expected quadrupole moments
   real(wp), intent(in) :: qpat(:, :)

   !> Factory to create electrostatic objects
   procedure(multipole_maker) :: make_multipole

   !> Energy of the multipole interaction
   real(wp), intent(out) :: energy

   type(damped_multipole) :: multipole
   type(container_cache) :: cache
   type(coulomb_cache), pointer :: ccache
   type(wavefunction_type) :: wfn

   real(wp), allocatable :: e(:)

   allocate(e(mol%nat), source=0.0_wp)

   wfn%qat = reshape(qat, [size(qat), 1])
   wfn%dpat = reshape(dpat, [shape(dpat), 1])
   wfn%qpat = reshape(qpat, [shape(qpat), 1])
   call taint(cache, ccache)
   call ccache%update(mol)
   call make_multipole(multipole, mol)
   call multipole%update(mol, cache)

   call multipole%get_energy(mol, cache, wfn, e)
   energy = sum(e)

end subroutine test_energy


subroutine test_energy_realspace(error, mol, qat, dpat, qpat, make_multipole, energy)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Expected multipole moments
   real(wp), intent(in) :: qat(:)

   !> Expected dipole moments
   real(wp), intent(in) :: dpat(:, :)

   !> Expected quadrupole moments
   real(wp), intent(in) :: qpat(:, :)

   !> Factory to create electrostatic objects
   procedure(multipole_maker) :: make_multipole

   !> Energy of the multipole interaction
   real(wp), intent(out) :: energy

   type(damped_multipole) :: multipole
   type(container_cache) :: cache
   type(coulomb_cache), pointer :: ccache

   real(wp) :: e_sd, e_dd, e_sq
   real(wp), allocatable :: e(:)

   allocate(e(mol%nat), source=0.0_wp)

   call taint(cache, ccache)
   call ccache%update(mol)
   call make_multipole(multipole, mol)
   call multipole%update(mol, cache)

   call get_multipole_energy_realspace(mol, ccache%mrad, multipole%kdmp3, multipole%kdmp5, &
      & qat, dpat, qpat, e_sd, e_dd, e_sq)
   print *, "e_sd", e_sd, "e_dd", e_dd, "e_sq", e_sq
   energy = e_sd + e_dd + e_sq

end subroutine test_energy_realspace


subroutine test_generic(error, mol, qat, dpat, qpat, make_multipole)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Expected multipole moments
   real(wp), intent(in) :: qat(:)

   !> Expected dipole moments
   real(wp), intent(in) :: dpat(:, :)

   !> Expected quadrupole moments
   real(wp), intent(in) :: qpat(:, :)

   !> Factory to create electrostatic objects
   procedure(multipole_maker) :: make_multipole

   integer, parameter :: supercell(3) = [2, 2, 2]

   type(structure_type) :: mol_sc
   real(wp) :: energy, energy_sc
   real(wp), allocatable :: qat_sc(:), dpat_sc(:, :), qpat_sc(:, :)

   call test_energy(error, mol, qat, dpat, qpat, make_multipole, energy)

   call make_supercell(mol, mol_sc, supercell)
   qat_sc = [qat, qat, qat, qat, qat, qat, qat, qat]
   dpat_sc = reshape([dpat, dpat, dpat, dpat, dpat, dpat, dpat, dpat], &
      & [size(dpat, 1), product(supercell)*size(dpat, 2)])
   qpat_sc = reshape([qpat, qpat, qpat, qpat, qpat, qpat, qpat, qpat], &
      & [size(qpat, 1), product(supercell)*size(qpat, 2)])

   call test_energy(error, mol_sc, qat_sc, dpat_sc, qpat_sc, make_multipole, energy_sc)

   call check(error, energy, energy_sc/product(supercell))
   if (allocated(error)) then
      print *, energy, energy_sc/product(supercell), energy - energy_sc/product(supercell)
   end if
end subroutine test_generic

subroutine test_realspace(error, mol, qat, dpat, qpat, make_multipole)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Expected multipole moments
   real(wp), intent(in) :: qat(:)

   !> Expected dipole moments
   real(wp), intent(in) :: dpat(:, :)

   !> Expected quadrupole moments
   real(wp), intent(in) :: qpat(:, :)

   !> Factory to create electrostatic objects
   procedure(multipole_maker) :: make_multipole

   integer, parameter :: supercell(3) = [2, 2, 2]

   type(structure_type) :: mol_sc
   real(wp) :: energy, energy_real, energy_sc, energy_real_sc
   real(wp), allocatable :: qat_sc(:), dpat_sc(:, :), qpat_sc(:, :)

   call test_energy(error, mol, qat, dpat, qpat, make_multipole, energy)
   call test_energy_realspace(error, mol, qat, dpat, qpat, make_multipole, energy_real)

   call make_supercell(mol, mol_sc, supercell)
   qat_sc = [qat, qat, qat, qat, qat, qat, qat, qat]
   dpat_sc = reshape([dpat, dpat, dpat, dpat, dpat, dpat, dpat, dpat], &
      & [size(dpat, 1), product(supercell)*size(dpat, 2)])
   qpat_sc = reshape([qpat, qpat, qpat, qpat, qpat, qpat, qpat, qpat], &
      & [size(qpat, 1), product(supercell)*size(qpat, 2)])

   call test_energy(error, mol_sc, qat_sc, dpat_sc, qpat_sc, make_multipole, energy_sc)
   call test_energy_realspace(error, mol_sc, qat_sc, dpat_sc, qpat_sc, make_multipole, energy_real_sc)

   call check(error, energy_real, energy_real_sc/product(supercell), thr=sqrt(epsilon(1.0_wp)))
   if (allocated(error)) then
      print *, "realspace/realspace(sc)", energy_real, energy_real_sc/product(supercell), &
         & energy_real - energy_real_sc/product(supercell)
      return
   end if

   call check(error, energy_real, energy)
   if (allocated(error)) then
      print *, "realspace/ewald", energy_real, energy, &
         & energy_real - energy
   end if

   call check(error, energy_real, energy_sc/product(supercell))
   if (allocated(error)) then
      print *, "realspace/ewald(sc)", energy_real, energy_sc/product(supercell), &
         & energy_real - energy_sc/product(supercell)
   end if
end subroutine test_realspace


!> Factory to create electrostatic objects based on GFN2-xTB values
subroutine make_multipole2(multipole, mol)

   !> New electrostatic object
   type(damped_multipole), intent(out) :: multipole

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   real(wp), parameter :: kdmp3 = 3.0_wp, kdmp5 = 3.0_wp
   real(wp), parameter :: shift = 1.2_wp, kexp = 4.0_wp, rmax = 5.0_wp
   !> Dipole exchange-correlation kernel
   real(wp), parameter :: p_dkernel(20) = 0.01_wp * [&
      & 5.563889_wp,-1.000000_wp,-0.500000_wp,-0.613341_wp,-0.481186_wp, &
      &-0.411674_wp, 3.521273_wp,-4.935670_wp,-8.339183_wp,10.000000_wp, &
      & 0.000000_wp,-0.082005_wp, 2.633341_wp,-0.025750_wp, 2.110225_wp, &
      &-0.151117_wp,-2.536958_wp,-2.077329_wp,-0.103383_wp,-0.236675_wp]
   !> Quadrupole exchange-correlation kernel
   real(wp), parameter :: p_qkernel(20) = 0.01_wp * [&
      & 0.027431_wp,-0.337528_wp, 0.020000_wp,-0.058586_wp,-0.058228_wp, &
      & 0.213583_wp, 2.026786_wp,-0.310828_wp,-0.245955_wp,-0.500000_wp, &
      & 0.020000_wp,-0.005516_wp,-0.021887_wp,-0.080000_wp, 0.028679_wp, &
      & 0.442859_wp, 0.122783_wp,-1.083404_wp, 0.025000_wp, 0.010000_wp]
   real(wp), parameter :: p_rad(20) = [&
      & 1.4_wp, 3.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 3.0_wp, 1.9_wp, 1.8_wp, 2.4_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 3.9_wp, 2.1_wp, 3.1_wp, 2.5_wp, 5.0_wp, 5.0_wp, 5.0_wp]
   real(wp), parameter :: p_vcn(20) = [&
      & 1.0_wp, 1.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 2.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 2.0_wp]
   real(wp), allocatable :: dkernel(:), qkernel(:), rad(:), vcn(:)

   dkernel = p_dkernel(mol%num)
   qkernel = p_qkernel(mol%num)
   rad = p_rad(mol%num)
   vcn = p_vcn(mol%num)

   call new_damped_multipole(multipole, mol, kdmp3, kdmp5, dkernel, qkernel, &
      & shift, kexp, rmax, rad, vcn)

end subroutine make_multipole2


!> Factory to create electrostatic objects based on GFN2-xTB values
subroutine make_multipole2_without_kernel(multipole, mol)

   !> New electrostatic object
   type(damped_multipole), intent(out) :: multipole

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   call make_multipole2(multipole, mol)
   multipole%dkernel(:) = 0.0_wp
   multipole%qkernel(:) = 0.0_wp
end subroutine make_multipole2_without_kernel


!> Inspect container cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(coulomb_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(coulomb_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

!> Return reference to container cache after resolving its type
subroutine view(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(coulomb_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(coulomb_cache)
      ptr => target
   end select
end subroutine view


subroutine test_e_effective_co2(error)
   
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(12) = [&
      &  4.95105332967126E-01_wp,  4.95110445149787E-01_wp,  4.95109208803526E-01_wp, &
      &  4.95110553060372E-01_wp, -2.47570208775520E-01_wp, -2.47372219387123E-01_wp, &
      & -2.47367867558844E-01_wp, -2.47541478966424E-01_wp, -2.47535153228645E-01_wp, &
      & -2.47738212772185E-01_wp, -2.47741337427010E-01_wp, -2.47569061865040E-01_wp]
   real(wp), parameter :: dpat(3, 12) = reshape([&
      & -3.26654706213550E-04_wp, -6.13403589170031E-05_wp,  2.83929462708564E-04_wp, &
      & -2.68534201975523E-04_wp, -1.57447758402670E-05_wp,  2.68727495879782E-04_wp, &
      & -2.74196985286096E-04_wp, -4.56547709680684E-05_wp,  2.33174149569738E-04_wp, &
      & -3.32314146822759E-04_wp, -6.08326819539977E-08_wp,  3.19473368568051E-04_wp, &
      & -5.42513865726882E-02_wp, -5.43318111661449E-02_wp, -5.44179905961818E-02_wp, &
      & -5.43660117536740E-02_wp, -5.44385268257028E-02_wp,  5.43650907998299E-02_wp, &
      & -5.43539636657199E-02_wp,  5.44469674651156E-02_wp,  5.43590516972356E-02_wp, &
      &  5.44387502680685E-02_wp, -5.43324554152907E-02_wp,  5.42581516232861E-02_wp, &
      &  5.44165886089117E-02_wp,  5.43323030762035E-02_wp,  5.42418925720472E-02_wp, &
      &  5.43211854223397E-02_wp,  5.42458343314803E-02_wp, -5.43133771609431E-02_wp, &
      &  5.43327026728458E-02_wp, -5.42359656240572E-02_wp, -5.43211640102117E-02_wp, &
      & -5.42499802018258E-02_wp,  5.43502803322248E-02_wp, -5.44204207635553E-02_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 12) = reshape([&
      &  1.18066656167315E-05_wp, -4.34798463050927E-01_wp,  1.44462097284581E-07_wp, &
      & -4.34798514781465E-01_wp, -4.34797971420059E-01_wp, -1.19511277154594E-05_wp, &
      &  2.25596069447498E-05_wp, -4.34783690153835E-01_wp, -2.11985471754161E-05_wp, &
      &  4.34768022081199E-01_wp,  4.34785018224248E-01_wp, -1.36105976999978E-06_wp, &
      &  2.25103495342660E-05_wp,  4.34767817555497E-01_wp,  1.07631391967900E-05_wp, &
      &  4.34783902623050E-01_wp, -4.34783120715119E-01_wp, -3.32734887309449E-05_wp, &
      & -9.48473787443227E-06_wp,  4.34783953044164E-01_wp,  1.08281396449250E-05_wp, &
      & -4.34783789399495E-01_wp,  4.34768955922080E-01_wp, -1.34340177049275E-06_wp, &
      & -9.66531291224371E-05_wp, -1.18702391374322E-01_wp, -1.37547479104327E-05_wp, &
      & -1.18645267151075E-01_wp, -1.18600873572414E-01_wp,  1.10407877032648E-04_wp, &
      & -2.26931493848559E-05_wp, -1.18637258774011E-01_wp,  5.06129802909649E-05_wp, &
      &  1.18680813186751E-01_wp,  1.18647169151941E-01_wp, -2.79198309067752E-05_wp, &
      & -3.70191853311663E-05_wp,  1.18641364188384E-01_wp,  6.65002198458886E-05_wp, &
      &  1.18692897340129E-01_wp, -1.18644737187210E-01_wp, -2.94810345136121E-05_wp, &
      &  1.06241518233796E-04_wp,  1.18591745614562E-01_wp,  9.27256724547743E-06_wp, &
      & -1.18653931140339E-01_wp,  1.18709157591422E-01_wp, -1.15514085479163E-04_wp, &
      &  9.75875617358346E-05_wp, -1.18606840689317E-01_wp,  1.39412337278877E-05_wp, &
      & -1.18668677913009E-01_wp, -1.18717811830031E-01_wp, -1.11528795463611E-04_wp, &
      &  3.00943076451121E-05_wp, -1.18658461646011E-01_wp, -6.19964990261623E-05_wp, &
      &  1.18611269561002E-01_wp,  1.18657681228184E-01_wp,  3.19021913810502E-05_wp, &
      &  4.34882210389453E-05_wp,  1.18645813671023E-01_wp, -6.06140176001579E-05_wp, &
      &  1.18607408054076E-01_wp, -1.18660243635834E-01_wp,  1.71257965608795E-05_wp, &
      & -1.15717831159601E-04_wp,  1.18704043515894E-01_wp, -3.90833585128814E-06_wp, &
      & -1.18646463864190E-01_wp,  1.18587233297690E-01_wp,  1.19626167010667E-04_wp],&
      & shape(qpat))
   real(wp), parameter :: qat0(12) = 0.0_wp, dpat0(3, 12) = 0.0_wp, qpat0(6, 12) = 0.0_wp

   call get_structure(mol, "X23", "CO2")
   call test_generic(error, mol, qat, dpat, qpat, make_multipole2)
   call test_generic(error, mol, qat0, dpat, qpat0, make_multipole2)
   call test_generic(error, mol, qat, dpat0, qpat, make_multipole2)

end subroutine test_e_effective_co2


subroutine test_e_effective_co2_realspace(error)
   
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(12) = [&
      &  4.95105332967126E-01_wp,  4.95110445149787E-01_wp,  4.95109208803526E-01_wp, &
      &  4.95110553060372E-01_wp, -2.47570208775520E-01_wp, -2.47372219387123E-01_wp, &
      & -2.47367867558844E-01_wp, -2.47541478966424E-01_wp, -2.47535153228645E-01_wp, &
      & -2.47738212772185E-01_wp, -2.47741337427010E-01_wp, -2.47569061865040E-01_wp]
   real(wp), parameter :: dpat(3, 12) = reshape([&
      & -3.26654706213550E-04_wp, -6.13403589170031E-05_wp,  2.83929462708564E-04_wp, &
      & -2.68534201975523E-04_wp, -1.57447758402670E-05_wp,  2.68727495879782E-04_wp, &
      & -2.74196985286096E-04_wp, -4.56547709680684E-05_wp,  2.33174149569738E-04_wp, &
      & -3.32314146822759E-04_wp, -6.08326819539977E-08_wp,  3.19473368568051E-04_wp, &
      & -5.42513865726882E-02_wp, -5.43318111661449E-02_wp, -5.44179905961818E-02_wp, &
      & -5.43660117536740E-02_wp, -5.44385268257028E-02_wp,  5.43650907998299E-02_wp, &
      & -5.43539636657199E-02_wp,  5.44469674651156E-02_wp,  5.43590516972356E-02_wp, &
      &  5.44387502680685E-02_wp, -5.43324554152907E-02_wp,  5.42581516232861E-02_wp, &
      &  5.44165886089117E-02_wp,  5.43323030762035E-02_wp,  5.42418925720472E-02_wp, &
      &  5.43211854223397E-02_wp,  5.42458343314803E-02_wp, -5.43133771609431E-02_wp, &
      &  5.43327026728458E-02_wp, -5.42359656240572E-02_wp, -5.43211640102117E-02_wp, &
      & -5.42499802018258E-02_wp,  5.43502803322248E-02_wp, -5.44204207635553E-02_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 12) = reshape([&
      &  1.18066656167315E-05_wp, -4.34798463050927E-01_wp,  1.44462097284581E-07_wp, &
      & -4.34798514781465E-01_wp, -4.34797971420059E-01_wp, -1.19511277154594E-05_wp, &
      &  2.25596069447498E-05_wp, -4.34783690153835E-01_wp, -2.11985471754161E-05_wp, &
      &  4.34768022081199E-01_wp,  4.34785018224248E-01_wp, -1.36105976999978E-06_wp, &
      &  2.25103495342660E-05_wp,  4.34767817555497E-01_wp,  1.07631391967900E-05_wp, &
      &  4.34783902623050E-01_wp, -4.34783120715119E-01_wp, -3.32734887309449E-05_wp, &
      & -9.48473787443227E-06_wp,  4.34783953044164E-01_wp,  1.08281396449250E-05_wp, &
      & -4.34783789399495E-01_wp,  4.34768955922080E-01_wp, -1.34340177049275E-06_wp, &
      & -9.66531291224371E-05_wp, -1.18702391374322E-01_wp, -1.37547479104327E-05_wp, &
      & -1.18645267151075E-01_wp, -1.18600873572414E-01_wp,  1.10407877032648E-04_wp, &
      & -2.26931493848559E-05_wp, -1.18637258774011E-01_wp,  5.06129802909649E-05_wp, &
      &  1.18680813186751E-01_wp,  1.18647169151941E-01_wp, -2.79198309067752E-05_wp, &
      & -3.70191853311663E-05_wp,  1.18641364188384E-01_wp,  6.65002198458886E-05_wp, &
      &  1.18692897340129E-01_wp, -1.18644737187210E-01_wp, -2.94810345136121E-05_wp, &
      &  1.06241518233796E-04_wp,  1.18591745614562E-01_wp,  9.27256724547743E-06_wp, &
      & -1.18653931140339E-01_wp,  1.18709157591422E-01_wp, -1.15514085479163E-04_wp, &
      &  9.75875617358346E-05_wp, -1.18606840689317E-01_wp,  1.39412337278877E-05_wp, &
      & -1.18668677913009E-01_wp, -1.18717811830031E-01_wp, -1.11528795463611E-04_wp, &
      &  3.00943076451121E-05_wp, -1.18658461646011E-01_wp, -6.19964990261623E-05_wp, &
      &  1.18611269561002E-01_wp,  1.18657681228184E-01_wp,  3.19021913810502E-05_wp, &
      &  4.34882210389453E-05_wp,  1.18645813671023E-01_wp, -6.06140176001579E-05_wp, &
      &  1.18607408054076E-01_wp, -1.18660243635834E-01_wp,  1.71257965608795E-05_wp, &
      & -1.15717831159601E-04_wp,  1.18704043515894E-01_wp, -3.90833585128814E-06_wp, &
      & -1.18646463864190E-01_wp,  1.18587233297690E-01_wp,  1.19626167010667E-04_wp],&
      & shape(qpat))
   real(wp), parameter :: qat0(12) = 0.0_wp, dpat0(3, 12) = 0.0_wp, qpat0(6, 12) = 0.0_wp

   call get_structure(mol, "X23", "CO2")
   call test_realspace(error, mol, qat0, dpat, qpat0, make_multipole2_without_kernel)
   call test_realspace(error, mol, qat, dpat0, qpat, make_multipole2_without_kernel)

end subroutine test_e_effective_co2_realspace

subroutine test_e_effective_urea(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      & 2.14480738661017E-1_wp, 2.14527495272549E-1_wp, 2.14471598997661E-1_wp, &
      & 2.14516468108033E-1_wp, 2.45821990520972E-1_wp, 2.45784832080333E-1_wp, &
      & 2.45806360746891E-1_wp, 2.45774007559891E-1_wp, 3.19474245000812E-1_wp, &
      & 3.19579496746895E-1_wp,-2.94626505330352E-1_wp,-2.94543754718096E-1_wp, &
      &-2.94659848864979E-1_wp,-2.94567600483171E-1_wp,-6.50780911474268E-1_wp, &
      &-6.51058612824167E-1_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      & 8.94827368009637E-2_wp, 8.94995909146529E-2_wp, 7.30634961756909E-2_wp, &
      & 8.95014203078334E-2_wp,-8.94859532646144E-2_wp,-7.30378298851327E-2_wp, &
      &-8.94873352015873E-2_wp,-8.95044727986339E-2_wp, 7.30601341120621E-2_wp, &
      &-8.95090716572443E-2_wp, 8.94905477532190E-2_wp,-7.30362173059338E-2_wp, &
      &-9.01917250610402E-3_wp,-9.00201247838739E-3_wp,-1.47850060857521E-1_wp, &
      &-9.00995133125046E-3_wp, 9.02740339798291E-3_wp, 1.47848725725542E-1_wp, &
      & 9.01739468689571E-3_wp, 9.00092652137544E-3_wp,-1.47851490297087E-1_wp, &
      & 9.00924018813641E-3_wp,-9.02542526993692E-3_wp, 1.47850396429899E-1_wp, &
      & 3.58533790400338E-6_wp, 1.76804223012015E-5_wp, 7.90811791429655E-2_wp, &
      & 1.90955806626401E-5_wp, 1.99177597125356E-6_wp,-7.93788543388031E-2_wp, &
      &-4.56180855347671E-2_wp,-4.55680071125238E-2_wp, 5.92409368467603E-2_wp, &
      &-4.56101826239023E-2_wp, 4.56639049054392E-2_wp,-5.91934296227654E-2_wp, &
      & 4.55978078230863E-2_wp, 4.55693703731601E-2_wp, 5.92330972112299E-2_wp, &
      & 4.56164302678149E-2_wp,-4.56412184067195E-2_wp,-5.91891814963370E-2_wp, &
      & 6.11781066870311E-6_wp, 9.15928265940359E-7_wp,-2.04395992132628E-1_wp, &
      & 1.78950595055862E-6_wp,-4.23972317385840E-6_wp, 2.04360605943730E-1_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      &-1.23698956475140E-2_wp, 1.79751870775559E-2_wp,-1.23691904599371E-2_wp, &
      & 4.67729267242994E-3_wp, 4.68004166312949E-3_wp, 2.47390861074511E-2_wp, &
      &-1.23654488813711E-2_wp,-1.79817963662729E-2_wp,-1.23706010436822E-2_wp, &
      &-4.69179866811513E-3_wp, 4.68829261745078E-3_wp, 2.47360499250533E-2_wp, &
      &-1.23743195999797E-2_wp, 1.79683161428587E-2_wp,-1.23718141043582E-2_wp, &
      &-4.68144974224393E-3_wp,-4.68203929764922E-3_wp, 2.47461337043379E-2_wp, &
      &-1.23715411335551E-2_wp,-1.79777626390457E-2_wp,-1.23696001364202E-2_wp, &
      & 4.69398264928409E-3_wp,-4.69414370569245E-3_wp, 2.47411412699752E-2_wp, &
      &-7.32193778340168E-3_wp, 3.45599884203983E-2_wp,-7.32434931893968E-3_wp, &
      &-6.43942937553289E-3_wp,-6.44505606234306E-3_wp, 1.46462871023413E-2_wp, &
      &-7.32508525229408E-3_wp,-3.45618871061652E-2_wp,-7.32335313431533E-3_wp, &
      & 6.45162772813276E-3_wp,-6.44644470074471E-3_wp, 1.46484383866092E-2_wp, &
      &-7.32905564328019E-3_wp, 3.45586433766952E-2_wp,-7.32543310737930E-3_wp, &
      & 6.44113915814194E-3_wp, 6.44697435437880E-3_wp, 1.46544887506594E-2_wp, &
      &-7.32605090318348E-3_wp,-3.45610747417528E-2_wp,-7.32899993430305E-3_wp, &
      &-6.45372642877205E-3_wp, 6.44746289062807E-3_wp, 1.46550508374864E-2_wp, &
      & 1.11951790334487E-1_wp,-4.02377811299945E-1_wp, 1.12018181924718E-1_wp, &
      &-1.04989650284058E-6_wp,-1.61793643693583E-5_wp,-2.23969972259205E-1_wp, &
      & 1.11996170927313E-1_wp, 4.02383573647660E-1_wp, 1.11929776104217E-1_wp, &
      & 1.69716696800307E-5_wp, 3.87072152745758E-6_wp,-2.23925947031530E-1_wp, &
      & 3.55892547851532E-2_wp,-1.40674101590403E-3_wp, 3.55796422243525E-2_wp, &
      &-4.51437869371043E-2_wp,-4.51403945941909E-2_wp,-7.11688970095042E-2_wp, &
      & 3.56007063955610E-2_wp, 1.39242785502690E-3_wp, 3.56047482812895E-2_wp, &
      & 4.51288756967405E-2_wp,-4.51327168181670E-2_wp,-7.12054546768490E-2_wp, &
      & 3.55791799999949E-2_wp,-1.38863474463774E-3_wp, 3.55709732649653E-2_wp, &
      & 4.51424181512024E-2_wp, 4.51358771366434E-2_wp,-7.11501532649592E-2_wp, &
      & 3.55872643357090E-2_wp, 1.37913694676448E-3_wp, 3.56010570423592E-2_wp, &
      &-4.51277501739264E-2_wp, 4.51338445551804E-2_wp,-7.11883213780670E-2_wp, &
      &-2.34031908350780E-2_wp,-5.97048165007334E-2_wp,-2.34051674127262E-2_wp, &
      &-1.53933429016999E-6_wp, 9.26688498685427E-7_wp, 4.68083582478040E-2_wp, &
      &-2.34383117547964E-2_wp, 5.96710054013645E-2_wp,-2.34363320491785E-2_wp, &
      &-2.13796886588954E-6_wp,-2.36716531036422E-6_wp, 4.68746438039748E-2_wp],&
      & shape(qpat))
   real(wp), parameter :: qat0(16) = 0.0_wp, dpat0(3, 16) = 0.0_wp, qpat0(6, 16) = 0.0_wp

   call get_structure(mol, "X23", "urea")
   call test_generic(error, mol, qat, dpat, qpat, make_multipole2)
   call test_generic(error, mol, qat0, dpat, qpat0, make_multipole2)
   call test_generic(error, mol, qat, dpat0, qpat, make_multipole2)

end subroutine test_e_effective_urea

subroutine test_e_effective_urea_realspace(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      & 2.14480738661017E-1_wp, 2.14527495272549E-1_wp, 2.14471598997661E-1_wp, &
      & 2.14516468108033E-1_wp, 2.45821990520972E-1_wp, 2.45784832080333E-1_wp, &
      & 2.45806360746891E-1_wp, 2.45774007559891E-1_wp, 3.19474245000812E-1_wp, &
      & 3.19579496746895E-1_wp,-2.94626505330352E-1_wp,-2.94543754718096E-1_wp, &
      &-2.94659848864979E-1_wp,-2.94567600483171E-1_wp,-6.50780911474268E-1_wp, &
      &-6.51058612824167E-1_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      & 8.94827368009637E-2_wp, 8.94995909146529E-2_wp, 7.30634961756909E-2_wp, &
      & 8.95014203078334E-2_wp,-8.94859532646144E-2_wp,-7.30378298851327E-2_wp, &
      &-8.94873352015873E-2_wp,-8.95044727986339E-2_wp, 7.30601341120621E-2_wp, &
      &-8.95090716572443E-2_wp, 8.94905477532190E-2_wp,-7.30362173059338E-2_wp, &
      &-9.01917250610402E-3_wp,-9.00201247838739E-3_wp,-1.47850060857521E-1_wp, &
      &-9.00995133125046E-3_wp, 9.02740339798291E-3_wp, 1.47848725725542E-1_wp, &
      & 9.01739468689571E-3_wp, 9.00092652137544E-3_wp,-1.47851490297087E-1_wp, &
      & 9.00924018813641E-3_wp,-9.02542526993692E-3_wp, 1.47850396429899E-1_wp, &
      & 3.58533790400338E-6_wp, 1.76804223012015E-5_wp, 7.90811791429655E-2_wp, &
      & 1.90955806626401E-5_wp, 1.99177597125356E-6_wp,-7.93788543388031E-2_wp, &
      &-4.56180855347671E-2_wp,-4.55680071125238E-2_wp, 5.92409368467603E-2_wp, &
      &-4.56101826239023E-2_wp, 4.56639049054392E-2_wp,-5.91934296227654E-2_wp, &
      & 4.55978078230863E-2_wp, 4.55693703731601E-2_wp, 5.92330972112299E-2_wp, &
      & 4.56164302678149E-2_wp,-4.56412184067195E-2_wp,-5.91891814963370E-2_wp, &
      & 6.11781066870311E-6_wp, 9.15928265940359E-7_wp,-2.04395992132628E-1_wp, &
      & 1.78950595055862E-6_wp,-4.23972317385840E-6_wp, 2.04360605943730E-1_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      &-1.23698956475140E-2_wp, 1.79751870775559E-2_wp,-1.23691904599371E-2_wp, &
      & 4.67729267242994E-3_wp, 4.68004166312949E-3_wp, 2.47390861074511E-2_wp, &
      &-1.23654488813711E-2_wp,-1.79817963662729E-2_wp,-1.23706010436822E-2_wp, &
      &-4.69179866811513E-3_wp, 4.68829261745078E-3_wp, 2.47360499250533E-2_wp, &
      &-1.23743195999797E-2_wp, 1.79683161428587E-2_wp,-1.23718141043582E-2_wp, &
      &-4.68144974224393E-3_wp,-4.68203929764922E-3_wp, 2.47461337043379E-2_wp, &
      &-1.23715411335551E-2_wp,-1.79777626390457E-2_wp,-1.23696001364202E-2_wp, &
      & 4.69398264928409E-3_wp,-4.69414370569245E-3_wp, 2.47411412699752E-2_wp, &
      &-7.32193778340168E-3_wp, 3.45599884203983E-2_wp,-7.32434931893968E-3_wp, &
      &-6.43942937553289E-3_wp,-6.44505606234306E-3_wp, 1.46462871023413E-2_wp, &
      &-7.32508525229408E-3_wp,-3.45618871061652E-2_wp,-7.32335313431533E-3_wp, &
      & 6.45162772813276E-3_wp,-6.44644470074471E-3_wp, 1.46484383866092E-2_wp, &
      &-7.32905564328019E-3_wp, 3.45586433766952E-2_wp,-7.32543310737930E-3_wp, &
      & 6.44113915814194E-3_wp, 6.44697435437880E-3_wp, 1.46544887506594E-2_wp, &
      &-7.32605090318348E-3_wp,-3.45610747417528E-2_wp,-7.32899993430305E-3_wp, &
      &-6.45372642877205E-3_wp, 6.44746289062807E-3_wp, 1.46550508374864E-2_wp, &
      & 1.11951790334487E-1_wp,-4.02377811299945E-1_wp, 1.12018181924718E-1_wp, &
      &-1.04989650284058E-6_wp,-1.61793643693583E-5_wp,-2.23969972259205E-1_wp, &
      & 1.11996170927313E-1_wp, 4.02383573647660E-1_wp, 1.11929776104217E-1_wp, &
      & 1.69716696800307E-5_wp, 3.87072152745758E-6_wp,-2.23925947031530E-1_wp, &
      & 3.55892547851532E-2_wp,-1.40674101590403E-3_wp, 3.55796422243525E-2_wp, &
      &-4.51437869371043E-2_wp,-4.51403945941909E-2_wp,-7.11688970095042E-2_wp, &
      & 3.56007063955610E-2_wp, 1.39242785502690E-3_wp, 3.56047482812895E-2_wp, &
      & 4.51288756967405E-2_wp,-4.51327168181670E-2_wp,-7.12054546768490E-2_wp, &
      & 3.55791799999949E-2_wp,-1.38863474463774E-3_wp, 3.55709732649653E-2_wp, &
      & 4.51424181512024E-2_wp, 4.51358771366434E-2_wp,-7.11501532649592E-2_wp, &
      & 3.55872643357090E-2_wp, 1.37913694676448E-3_wp, 3.56010570423592E-2_wp, &
      &-4.51277501739264E-2_wp, 4.51338445551804E-2_wp,-7.11883213780670E-2_wp, &
      &-2.34031908350780E-2_wp,-5.97048165007334E-2_wp,-2.34051674127262E-2_wp, &
      &-1.53933429016999E-6_wp, 9.26688498685427E-7_wp, 4.68083582478040E-2_wp, &
      &-2.34383117547964E-2_wp, 5.96710054013645E-2_wp,-2.34363320491785E-2_wp, &
      &-2.13796886588954E-6_wp,-2.36716531036422E-6_wp, 4.68746438039748E-2_wp],&
      & shape(qpat))
   real(wp), parameter :: qat0(16) = 0.0_wp, dpat0(3, 16) = 0.0_wp, qpat0(6, 16) = 0.0_wp

   call get_structure(mol, "X23", "urea")
   call test_realspace(error, mol, qat0, dpat, qpat0, make_multipole2_without_kernel)
   call test_realspace(error, mol, qat, dpat0, qpat, make_multipole2_without_kernel)

end subroutine test_e_effective_urea_realspace


subroutine test_s_effective_oxacb(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      & 3.83148728601709E-1_wp, 3.83118715934305E-1_wp, 3.83080191213008E-1_wp, &
      & 3.82184417033339E-1_wp, 3.47311885772574E-1_wp, 3.45281176165792E-1_wp, &
      & 3.47791244554054E-1_wp, 3.47952947802927E-1_wp,-3.92870149599810E-1_wp, &
      &-3.92136115450962E-1_wp,-3.93531741241129E-1_wp,-3.92928999213039E-1_wp, &
      &-3.37416498051037E-1_wp,-3.36550419758369E-1_wp,-3.37512222809576E-1_wp, &
      &-3.36923160953788E-1_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      & 4.09688027943853E-2_wp, 8.33483477302867E-3_wp,-1.09648576234461E-2_wp, &
      &-4.09474996408295E-2_wp, 8.30233124676478E-3_wp, 1.09874908009843E-2_wp, &
      &-4.09278272213102E-2_wp,-8.32793028126641E-3_wp, 1.09615760238990E-2_wp, &
      & 4.11086171124463E-2_wp,-8.39441030068160E-3_wp,-1.05971938523367E-2_wp, &
      & 5.96300641349022E-2_wp,-4.34579903309867E-2_wp, 6.26025852647896E-2_wp, &
      &-5.98049356163273E-2_wp,-4.34306642465664E-2_wp,-6.22704459500676E-2_wp, &
      &-6.00922329406000E-2_wp, 4.36215944530855E-2_wp,-6.21183058238615E-2_wp, &
      & 5.93589247389314E-2_wp, 4.37917586457069E-2_wp, 6.24175529675174E-2_wp, &
      &-9.47993007577910E-2_wp,-9.21933679240633E-2_wp, 1.41162420287618E-1_wp, &
      & 9.47847920056643E-2_wp,-9.22795166591021E-2_wp,-1.41330245385699E-1_wp, &
      & 9.47957155116519E-2_wp, 9.20217761595336E-2_wp,-1.41075167719715E-1_wp, &
      &-9.47612703050034E-2_wp, 9.20970946041004E-2_wp, 1.41174078870679E-1_wp, &
      & 3.07186908919029E-3_wp, 9.78091090173267E-2_wp,-1.59405890565398E-1_wp, &
      &-2.89501267712337E-3_wp, 9.81266227352185E-2_wp, 1.59671635941258E-1_wp, &
      &-3.07984597314511E-3_wp,-9.78272060157212E-2_wp, 1.59285689458381E-1_wp, &
      & 3.45622260398958E-3_wp,-9.78261483786836E-2_wp,-1.59383710436135E-1_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      & 6.56379864949194E-2_wp,-5.73520679256482E-4_wp,-3.52084956587255E-2_wp, &
      & 7.71130683289638E-3_wp,-1.49139358369174E-2_wp,-3.04294908361939E-2_wp, &
      & 6.54814038446620E-2_wp, 5.17836513451199E-4_wp,-3.51478340536698E-2_wp, &
      & 7.70830546882409E-3_wp, 1.49356522560081E-2_wp,-3.03335697909922E-2_wp, &
      & 6.55759132522838E-2_wp,-5.57374828410120E-4_wp,-3.51772230507483E-2_wp, &
      & 7.71886247591620E-3_wp,-1.49101894592223E-2_wp,-3.03986902015354E-2_wp, &
      & 6.44858926966245E-2_wp, 6.21105139273002E-4_wp,-3.60657223374486E-2_wp, &
      & 8.53829013403492E-3_wp, 1.44562000883513E-2_wp,-2.84201703591760E-2_wp, &
      &-2.16992043497083E-1_wp,-6.14525993351003E-2_wp, 2.61869112357575E-1_wp, &
      & 7.43403756029076E-2_wp, 3.29057424398825E-1_wp,-4.48770688604918E-2_wp, &
      &-2.17253885199018E-1_wp, 6.16371296152827E-2_wp, 2.61441601565784E-1_wp, &
      & 7.36098779877098E-2_wp,-3.29050085539625E-1_wp,-4.41877163667672E-2_wp, &
      &-2.17095729603948E-1_wp,-6.13633493424920E-2_wp, 2.61949523181417E-1_wp, &
      & 7.41593953178476E-2_wp, 3.29100144825087E-1_wp,-4.48537935774700E-2_wp, &
      &-2.16773170522421E-1_wp, 6.14170456026348E-2_wp, 2.61763772584928E-1_wp, &
      & 7.44868464063092E-2_wp,-3.29037954066144E-1_wp,-4.49906020625066E-2_wp, &
      &-7.84671583301996E-2_wp,-2.93177074689629E-3_wp, 9.45239916587787E-2_wp, &
      &-1.35552573154701E-3_wp, 1.17958224427535E-1_wp,-1.60568333285798E-2_wp, &
      &-7.86660963425362E-2_wp, 3.27173353036481E-3_wp, 9.46888561547714E-2_wp, &
      &-1.70110485836640E-3_wp,-1.18041188489037E-1_wp,-1.60227598122358E-2_wp, &
      &-7.83224933313823E-2_wp,-2.95614679715344E-3_wp, 9.43219272902749E-2_wp, &
      &-1.39062937509669E-3_wp, 1.17927803395099E-1_wp,-1.59994339588925E-2_wp, &
      &-7.84941989805875E-2_wp, 2.98202100370017E-3_wp, 9.44901433808538E-2_wp, &
      &-1.32436849481419E-3_wp,-1.18022315802278E-1_wp,-1.59959444002657E-2_wp, &
      &-4.89935283276599E-3_wp,-1.35938726283863E-2_wp, 8.75281455685861E-3_wp, &
      & 1.70009780794440E-2_wp, 1.14577178485397E-2_wp,-3.85346172409284E-3_wp, &
      &-4.97349323995833E-3_wp, 1.37313962587995E-2_wp, 8.85061916230923E-3_wp, &
      & 1.71435461297655E-2_wp,-1.12771254648400E-2_wp,-3.87712592235134E-3_wp, &
      &-4.80233859626600E-3_wp,-1.36263939310827E-2_wp, 8.80915005635918E-3_wp, &
      & 1.70027332793215E-2_wp, 1.14775835735027E-2_wp,-4.00681146009285E-3_wp, &
      &-4.69609623172096E-3_wp, 1.35220039889045E-2_wp, 8.69541313382693E-3_wp, &
      & 1.69993038311265E-2_wp,-1.13773019767008E-2_wp,-3.99931690210575E-3_wp],&
      & shape(qpat))

   call get_structure(mol, "X23", "oxacb")
   ! call test_numsigma(error, mol, qat, dpat, qpat, make_multipole2)

end subroutine test_s_effective_oxacb

subroutine make_supercell(mol, mol_sc, rep)
   type(structure_type), intent(in) :: mol
   type(structure_type), intent(out) :: mol_sc
   integer, intent(in) :: rep(3)

   real(wp), allocatable :: xyz(:, :), lattice(:, :)
   integer, allocatable :: num(:)
   integer :: i, j, k, c

   num = reshape(spread([mol%num(mol%id)], 2, product(rep)), [product(rep)*mol%nat])
   lattice = reshape(&
      [rep(1)*mol%lattice(:, 1), rep(2)*mol%lattice(:, 2), rep(3)*mol%lattice(:, 3)], &
      shape(mol%lattice))
   allocate(xyz(3, product(rep)*mol%nat))
   c = 0
   do i = 0, rep(1)-1
      do j = 0, rep(2)-1
         do k = 0, rep(3)-1
            xyz(:, c+1:c+mol%nat) = mol%xyz &
               & + spread(matmul(mol%lattice, [real(wp):: i, j, k]), 2, mol%nat)
            c = c + mol%nat
         end do
      end do
   end do

   call new(mol_sc, num, xyz, lattice=lattice)
end subroutine make_supercell

subroutine get_multipole_energy_realspace(mol, rad, kdmp3, kdmp5, qat, dpat, qpat, e_sd, e_dd, e_sq)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multipole damping radii for all atoms
   real(wp), intent(in) :: rad(:)
   !> Damping function for inverse quadratic contributions
   real(wp), intent(in) :: kdmp3
   !> Damping function for inverse cubic contributions
   real(wp), intent(in) :: kdmp5
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Atomic dipole moments
   real(wp), intent(in) :: dpat(:, :)
   !> Atomic quadrupole moments
   real(wp), intent(in) :: qpat(:, :)
   !> Energy for charges and dipoles
   real(wp), intent(out) :: e_sd
   !> Energy for dipoles and dipoles
   real(wp), intent(out) :: e_dd
   !> Energy for charges and quadrupoles
   real(wp), intent(out) :: e_sq

   integer :: iat, jat, itr
   real(wp) :: r1, vec(3), g1, g3, fdmp3, fdmp5, tc(6), rr, dr(3, 3), t3, t5
   real(wp) :: sidj(3), sjdi(3), siqj(6), sjqi(6), didj(3, 3)
   real(wp), allocatable :: trans(:, :), e_kernel(:)

   real(wp), parameter :: unity(3, 3) = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   real(wp), parameter :: cutoff = 120.0_wp

   call get_lattice_points([.true.], mol%lattice, cutoff, trans)

   e_sd = 0.0_wp
   e_dd = 0.0_wp
   e_sq = 0.0_wp
   do iat = 1, mol%nat
      do jat = 1, mol%nat
         sidj = qat(iat) * dpat(:, jat)
         sjdi = qat(jat) * dpat(:, iat)
         didj = spread(dpat(:, iat), 1, 3) * spread(dpat(:, jat), 2, 3)
         siqj = qat(iat) * qpat(:, jat)
         sjqi = qat(jat) * qpat(:, iat)
         do itr = 1, size(trans, 2)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)

            r1 = norm2(vec)
            if (r1 < epsilon(1.0_wp) .or. r1 > cutoff) cycle
            g1 = 1.0_wp / r1
            g3 = g1 * g1 * g1

            dr(:, :) = spread(vec, 1, 3) * spread(vec, 2, 3) * g1*g1

            rr = 0.5_wp * (rad(jat) + rad(iat)) * g1
            fdmp3 = 1.0_wp / (1.0_wp + 6.0_wp * rr**kdmp3)
            fdmp5 = 1.0_wp / (1.0_wp + 6.0_wp * rr**kdmp5)

            t3 = g3 * fdmp3
            t5 = g3 * fdmp5

            e_sd = e_sd + t3 * dot_product(sidj + sjdi, vec) * 0.5_wp
            e_dd = e_dd + sum((unity - 3 * dr) * t5 * didj) * 0.5_wp
            tc(2) = 2*dr(1, 2)*t5
            tc(4) = 2*dr(1, 3)*t5
            tc(5) = 2*dr(2, 3)*t5
            tc(1) = dr(1, 1)*t5
            tc(3) = dr(2, 2)*t5
            tc(6) = dr(3, 3)*t5
            e_sq = e_sq + dot_product(siqj + sjqi, tc) * 0.5_wp
         end do
      end do
   end do
end subroutine get_multipole_energy_realspace

end module test_coulomb_multipole_periodic