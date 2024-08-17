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

module test_solvation_cds
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_container, only : container_cache
   use tblite_scf_potential, only : potential_type
   use tblite_solvation_alpb
   use tblite_solvation_cds
   use tblite_data_cds
   use tblite_solvation_data
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: collect_solvation_cds


   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains


!> Collect all exported unit tests
subroutine collect_solvation_cds(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("energy-neutral", test_e_neutral), &
      new_unittest("energy-charged", test_e_charged), &
      new_unittest("gradient-scf", test_g_cds), &
      new_unittest("gradient-nonscf", test_g_cdsnonscf), &
      new_unittest("potential", test_p_cds) &
      ]

end subroutine collect_solvation_cds


subroutine test_e(error, mol, input, qat, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Solvation model input
   type(cds_input), intent(in) :: input

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Reference energy
   real(wp), intent(in) :: ref


   type(cds_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp) :: energy(mol%nat)
   type(cds_input) :: tmpinput

   energy = 0.0_wp
   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1))

   tmpinput = input
   call get_cds_param(tmpinput, mol, error)
   if(allocated(error))then
     call test_failed(error, "Failed to get solvation CDS parameters")
   endif
   solv = cds_solvation(mol, tmpinput)

   call solv%update(mol, cache)
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)
   !> CDS has some non-selfconsistent part as well:
   call solv%get_engrad(mol, cache, energy)

   if (abs(sum(energy) - ref) > thr) then
      call test_failed(error, "Energy does not match reference")
      print *, sum(energy),'reference:',ref
   end if
end subroutine test_e


subroutine test_g(error, mol, input, qat)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Solvation model input
   type(cds_input), intent(in) :: input

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   type(cds_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp), allocatable :: gradient(:, :), numg(:, :)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3)
   integer :: ii, ic
   type(cds_input) :: tmpinput

   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1))

   tmpinput=input
   call get_cds_param(tmpinput, mol, error)
   if(allocated(error))then
     call test_failed(error, "Failed to get solvation CDS parameters")
   endif
   solv = cds_solvation(mol, tmpinput)

   allocate(numg(3, mol%nat), gradient(3, mol%nat))
   do ii = 1, mol%nat
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(ic, ii) = mol%xyz(ic, ii) + step
         call solv%update(mol, cache)
         call solv%get_potential(mol, cache, wfn, pot)
         call solv%get_energy(mol, cache, wfn, er)

         mol%xyz(ic, ii) = mol%xyz(ic, ii) - 2*step
         call solv%update(mol, cache)
         call solv%get_potential(mol, cache, wfn, pot)
         call solv%get_energy(mol, cache, wfn, el)

         mol%xyz(ic, ii) = mol%xyz(ic, ii) + step
         numg(ic, ii) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp

   call solv%update(mol, cache)
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)
   call solv%get_gradient(mol, cache, wfn, gradient, sigma)

   if (any(abs(gradient - numg) > thr)) then
      call test_failed(error, "Gradient does not match")
      print '(3es20.13)', gradient
      print '(a)', "---"
      print '(3es20.13)', numg
      print '(a)', "---"
      print '(3es20.13)', gradient - numg
   end if
end subroutine test_g


subroutine test_g_nonscf(error, mol, input, qat)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Solvation model input
   type(cds_input), intent(in) :: input

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   type(cds_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp), allocatable :: gradient(:, :), numg(:, :)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3)
   integer :: ii, ic
   type(cds_input) :: tmpinput

   tmpinput = input 
   call get_cds_param(tmpinput, mol, error)
   if(allocated(error))then
     call test_failed(error, "Failed to get solvation CDS parameters")
   endif
   solv = cds_solvation(mol, tmpinput)

   allocate(numg(3, mol%nat), gradient(3, mol%nat))
   do ii = 1, mol%nat
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(ic, ii) = mol%xyz(ic, ii) + step
         call solv%update(mol, cache)
         call solv%get_engrad(mol, cache, er)

         mol%xyz(ic, ii) = mol%xyz(ic, ii) - 2*step
         call solv%update(mol, cache)
         call solv%get_engrad(mol, cache, el)

         mol%xyz(ic, ii) = mol%xyz(ic, ii) + step
         numg(ic, ii) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp

   call solv%update(mol, cache)
   call solv%get_engrad(mol, cache, energy, gradient)

   if (any(abs(gradient - numg) > thr)) then
      call test_failed(error, "Gradient does not match")
      print '(3es20.13)', gradient
      print '(a)', "---"
      print '(3es20.13)', numg
      print '(a)', "---"
      print '(3es20.13)', gradient - numg
   end if
end subroutine test_g_nonscf

subroutine test_p(error, mol, input, qat)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Solvation model input
   type(cds_input), intent(in) :: input

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   type(cds_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), parameter :: thr = 1e+3_wp*sqrt(epsilon(1.0_wp))
   real(wp), allocatable :: vat(:)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat)
   integer :: ii
   type(cds_input) :: tmpinput


   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1))

   tmpinput = input
   call get_cds_param(tmpinput, mol, error)
   if(allocated(error))then
     call test_failed(error, "Failed to get solvation CDS parameters")
   endif
   solv = cds_solvation(mol, tmpinput)

   call solv%update(mol, cache)

   allocate(vat(mol%nat))
   do ii = 1, mol%nat
      er = 0.0_wp
      el = 0.0_wp
      wfn%qat(ii, 1) = wfn%qat(ii, 1) + step
      call solv%get_potential(mol, cache, wfn, pot)
      call solv%get_energy(mol, cache, wfn, er)

      wfn%qat(ii, 1) = wfn%qat(ii, 1) - 2*step
      call solv%get_potential(mol, cache, wfn, pot)
      call solv%get_energy(mol, cache, wfn, el)

      wfn%qat(ii, 1) = wfn%qat(ii, 1) + step
      vat(ii) = 0.5_wp*(sum(er) - sum(el))/step
   end do

   energy = 0.0_wp
   pot%vat(:, :) = 0.0_wp
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)

   if (any(abs([pot%vat] - vat) > thr)) then
      call test_failed(error, "Potential does not match")
      print '(3es20.13)', pot%vat
      print '(a)', "---"
      print '(3es20.13)', vat
      print '(a)', "---"
      print '(3es20.13)', [pot%vat] - vat
   end if
end subroutine test_p


subroutine test_e_neutral(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-8.99890404486076E-2_wp, 9.42168087556583E-2_wp,-1.49387631509499E-1_wp, &
      &-2.99114121895542E-1_wp, 4.85527734875224E-1_wp,-6.83156326406137E-2_wp, &
      & 1.50011293889337E-2_wp, 2.79368544459465E-1_wp,-1.24072452878322E-1_wp, &
      &-9.36760994051244E-2_wp,-2.19062123031622E-1_wp, 2.14538817685587E-1_wp, &
      & 3.06156726072831E-1_wp,-3.86105514712244E-1_wp,-1.51265171389388E-3_wp, &
      & 3.64255069977693E-2_wp]
   real(wp), parameter :: feps = 80.0_wp

   call get_structure(mol, "MB16-43", "04")

   call test_e(error, mol, cds_input(alpb=.true., solvent='aniline',      method='gfn1'), qat, -7.5989750654943374d-003)
   call test_e(error, mol, cds_input(alpb=.true., solvent='cs2',          method='gfn1'), qat, -1.1405775310930399d-002)
   call test_e(error, mol, cds_input(alpb=.true., solvent='methanol',     method='gfn1'), qat, -9.4336600684243963d-003)
   call test_e(error, mol, cds_input(alpb=.true., solvent='dmf',          method='gfn1'), qat, -1.4452113472146064d-002)
   if(allocated(error)) return

   call test_e(error, mol, cds_input(alpb=.false., solvent='ch2cl2',       method='gfn1'), qat,  -1.2674622292821234d-002)
   call test_e(error, mol, cds_input(alpb=.false., solvent='dmso',         method='gfn1'), qat,  -1.0251673911959865d-002)
   call test_e(error, mol, cds_input(alpb=.false., solvent='water',        method='gfn1'), qat,   2.8274080367559997d-003)
   call test_e(error, mol, cds_input(alpb=.false., solvent='thf',          method='gfn1'), qat,  -9.1185514342106630d-003)
   if(allocated(error)) return

   call test_e(error, mol, cds_input(alpb=.true., solvent='acetone',      method='gfn2'), qat, -1.5608029177684101d-002)
   call test_e(error, mol, cds_input(alpb=.true., solvent='benzene',      method='gfn2'), qat, -1.3607685212954617d-002)
   call test_e(error, mol, cds_input(alpb=.true., solvent='hexadecane',   method='gfn2'), qat, -1.2559137678418444d-002)
   call test_e(error, mol, cds_input(alpb=.true., solvent='water',        method='gfn2'), qat,  8.2141194896751066d-004)
   call test_e(error, mol, cds_input(alpb=.true., solvent='nhexane',      method='gfn2'), qat, -1.0916108842590913d-002)
   if(allocated(error)) return

   call test_e(error, mol, cds_input(alpb=.false., solvent='acetonitrile', method='gfn2'), qat,  -5.2242386726497692d-003)
   call test_e(error, mol, cds_input(alpb=.false., solvent='chcl3',        method='gfn2'), qat,  -1.3510126366481480d-002)
   call test_e(error, mol, cds_input(alpb=.false., solvent='ether',        method='gfn2'), qat,  -1.2799880638529441d-002)
   call test_e(error, mol, cds_input(alpb=.false., solvent='toluene',      method='gfn2'), qat,  -1.4112628162725245d-002)
   if(allocated(error)) return

end subroutine test_e_neutral

subroutine test_e_charged(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 2.91737167991831E-1_wp,-4.44921281646692E-1_wp, 7.03886898899541E-2_wp, &
      & 4.32458130745259E-2_wp, 3.90393832017139E-2_wp, 9.56766032864156E-2_wp, &
      & 6.01866392558221E-2_wp,-3.31751380042108E-1_wp, 1.56707781479786E-1_wp, &
      & 8.89443731456522E-2_wp,-6.03526873354187E-2_wp, 3.83058968439732E-2_wp, &
      & 9.51289165540499E-2_wp,-9.20860756307965E-2_wp, 7.12122392663975E-2_wp, &
      & 2.76742553014572E-1_wp,-4.38083248713717E-1_wp,-1.76798298145349E-1_wp, &
      & 2.14382387079637E-1_wp, 3.25215078814299E-1_wp,-4.30846325150764E-1_wp, &
      & 1.12908046155543E-1_wp, 6.61648434849400E-2_wp, 9.89171412174742E-2_wp, &
      & 5.80323177221752E-2_wp,-4.43388471029542E-1_wp, 3.38493755177292E-1_wp, &
      &-3.52969753592199E-1_wp, 4.59552645375000E-1_wp,-5.33442813398395E-1_wp, &
      &-5.23460418980519E-1_wp,-3.14158345734236E-1_wp, 8.63777453496141E-2_wp, &
      & 5.51896920278464E-2_wp, 5.60421076086023E-2_wp, 1.02875052304244E-1_wp, &
      & 5.43513415521492E-2_wp,-3.13791816210054E-1_wp, 1.68862983166011E-1_wp, &
      & 8.93759201291419E-2_wp,-8.90121909290432E-2_wp, 4.15150867426933E-2_wp, &
      & 1.22221651251480E-1_wp,-8.26447904844349E-2_wp, 1.00154406589010E-1_wp, &
      & 2.77572579099978E-1_wp,-4.29147369583175E-1_wp,-1.78581481555413E-1_wp, &
      & 2.09487890121871E-1_wp, 3.17109649407645E-1_wp,-4.62948099575570E-1_wp, &
      & 9.81620738022878E-2_wp, 5.14984224707990E-2_wp, 9.63222020737258E-2_wp, &
      & 3.80443799704811E-2_wp,-4.41189291092377E-1_wp, 3.13549746324888E-1_wp, &
      &-4.41335746902051E-1_wp, 3.01219329594079E-1_wp]

   call get_structure(mol, "UPU23", "0a")
   call test_e(error, mol, cds_input(alpb=.true., solvent='water', method='gfn2'), qat, -7.4055525062074719d-003)
   if(allocated(error)) return
   call test_e(error, mol, cds_input(alpb=.false., solvent='water', method='gfn2'), qat, -1.6875062934239744d-002)
   if(allocated(error)) return
   call test_e(error, mol, cds_input(alpb=.true., solvent='water', method='gfn1'), qat, -1.7447404840496523d-002)
   if(allocated(error)) return
   call test_e(error, mol, cds_input(alpb=.false., solvent='water', method='gfn1'), qat, 8.7872092832611444d-003)
   if(allocated(error)) return

end subroutine test_e_charged


subroutine test_g_cds(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 2.08159387594211E-1_wp,-3.78010519998818E-1_wp, 3.36498247356244E-2_wp, &
      &-4.11556158912895E-1_wp, 8.14928196660512E-2_wp,-2.00886649303053E-1_wp, &
      & 2.44756994282684E-1_wp, 2.54580499189089E-2_wp, 2.59835128092562E-1_wp, &
      & 4.21683321877209E-1_wp, 1.37097163086023E-1_wp, 4.06951664942900E-2_wp, &
      &-1.10955378625897E-1_wp,-6.44033540918074E-2_wp,-1.91525919028143E-1_wp, &
      &-9.54898757869102E-2_wp]

   call get_structure(mol, "MB16-43", "06")
   call test_g(error, mol, cds_input(alpb=.true., solvent='water', method='gfn1'), qat)
   if(allocated(error)) return
   call test_g(error, mol, cds_input(alpb=.true., solvent='water', method='gfn2'), qat)
   if(allocated(error)) return

end subroutine test_g_cds

subroutine test_g_cdsnonscf(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 2.08159387594211E-1_wp,-3.78010519998818E-1_wp, 3.36498247356244E-2_wp, &
      &-4.11556158912895E-1_wp, 8.14928196660512E-2_wp,-2.00886649303053E-1_wp, &
      & 2.44756994282684E-1_wp, 2.54580499189089E-2_wp, 2.59835128092562E-1_wp, &
      & 4.21683321877209E-1_wp, 1.37097163086023E-1_wp, 4.06951664942900E-2_wp, &
      &-1.10955378625897E-1_wp,-6.44033540918074E-2_wp,-1.91525919028143E-1_wp, &
      &-9.54898757869102E-2_wp]

   call get_structure(mol, "MB16-43", "06")
   call test_g_nonscf(error, mol, cds_input(alpb=.true., solvent='water', method='gfn2'), qat)
   if(allocated(error)) return
   call test_g_nonscf(error, mol, cds_input(alpb=.false., solvent='water', method='gfn1'), qat)
   if(allocated(error)) return     

end subroutine test_g_cdsnonscf


subroutine test_p_cds(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-2.05668345919710E-1_wp,-3.99553123071811E-1_wp, 3.29242774348191E-1_wp, &
      &-3.11737933844111E-1_wp, 3.58851882478133E-2_wp, 3.21886835736497E-1_wp, &
      & 4.14743455841314E-2_wp, 2.95727359478547E-2_wp,-5.06347224522431E-1_wp, &
      & 3.43067182413129E-1_wp, 6.88373767679515E-1_wp, 7.03357390141253E-2_wp, &
      &-9.62424552888750E-2_wp,-1.32209348056625E-1_wp, 9.78998441832186E-2_wp, &
      &-3.05979982450903E-1_wp]

   call get_structure(mol, "MB16-43", "08")
   call test_p(error, mol, cds_input(alpb=.true., solvent='water', method='gfn2'), qat)
   if(allocated(error)) return
   call test_p(error, mol, cds_input(alpb=.false., solvent='water', method='gfn2'), qat)
   if(allocated(error)) return
   call test_p(error, mol, cds_input(alpb=.true., solvent='water', method='gfn1'), qat)
   if(allocated(error)) return
   call test_p(error, mol, cds_input(alpb=.false., solvent='water', method='gfn1'), qat)
   if(allocated(error)) return

end subroutine test_p_cds


end module test_solvation_cds
