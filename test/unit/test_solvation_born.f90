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

module test_solvation_born
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_container, only : container_cache
   use tblite_scf_potential, only : potential_type
   use tblite_solvation_alpb
   use tblite_solvation_born
   use tblite_solvation_data
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: collect_solvation_born


   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains


!> Collect all exported unit tests
subroutine collect_solvation_born(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("born-1", test_mb01), &
      new_unittest("born-2", test_mb02), &
      new_unittest("born-3", test_mb03), &
      new_unittest("energy-p16", test_e_p16), &
      new_unittest("energy-still", test_e_still), &
      new_unittest("energy-charged", test_e_charged), &
      new_unittest("gradient-p16", test_g_p16), &
      new_unittest("gradient-still", test_g_still), &
      new_unittest("potential-p16", test_p_p16), &
      new_unittest("potential-still", test_p_still) &
      ]

end subroutine collect_solvation_born


subroutine test_numg(error, gbobc, mol)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Born radii integrator
   type(born_integrator), intent(inout) :: gbobc
   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   integer :: iat, ic
   real(wp), allocatable :: rad(:), sr(:), sl(:), draddr(:, :, :), numg(:, :, :)
   real(wp), parameter :: step = 1.0e-5_wp

   allocate(rad(mol%nat), sr(mol%nat), sl(mol%nat))
   allocate(draddr(3, mol%nat, mol%nat), numg(3, mol%nat, mol%nat))

   call gbobc%get_rad(mol, rad, draddr)

   do iat = 1, mol%nat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call gbobc%get_rad(mol, sr)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call gbobc%get_rad(mol, sl)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step

         numg(ic, iat, :) = 0.5_wp * (sr - sl) / step
      end do
   end do

   if (any(abs(numg - draddr) > thr2)) then
      call test_failed(error, "Born radii derivative does not much finite difference solution")
   end if
end subroutine test_numg


subroutine test_mb01(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   real(wp), allocatable :: rvdw(:), rad(:), draddr(:, :, :)
   real(wp), parameter :: ref(16) = [&
      & 4.07331798531438E+0_wp, 2.56451650429167E+0_wp, 2.97676954448882E+0_wp, &
      & 2.88833112759592E+0_wp, 2.59127476011008E+0_wp, 2.63750279425510E+0_wp, &
      & 3.56149571036025E+0_wp, 2.89090958281373E+0_wp, 4.53815592283277E+0_wp, &
      & 2.46342847720303E+0_wp, 2.72707461251522E+0_wp, 3.52933932564532E+0_wp, &
      & 3.66934919146868E+0_wp, 3.66697827876019E+0_wp, 3.37809284756764E+0_wp, &
      & 4.57932013403544E+0_wp]

   call get_structure(mol, "MB16-43", "01")

   allocate(rad(mol%nat), draddr(3, mol%nat, mol%nat))
   rvdw = get_vdw_rad_d3(mol%num)

   call new_born_integrator(gbobc, mol, rvdw)
   call gbobc%get_rad(mol, rad, draddr)

   if (any(abs(rad - ref) > thr2)) then
      call test_failed(error, "Born radii area values do not match")
      print '(es20.14e1)', rad
      return
   end if

   call test_numg(error, gbobc, mol)

end subroutine test_mb01


subroutine test_mb02(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   real(wp), allocatable :: rvdw(:), rad(:), draddr(:, :, :)
   real(wp), parameter :: ref(16) = [&
      & 3.43501192886252E+0_wp, 5.02404916999371E+0_wp, 4.72253865039692E+0_wp, &
      & 3.52096661104217E+0_wp, 4.76023330956437E+0_wp, 3.46119195863261E+0_wp, &
      & 3.17361370475619E+0_wp, 2.90775065608382E+0_wp, 4.94595287355805E+0_wp, &
      & 3.32592657749444E+0_wp, 4.54348353409109E+0_wp, 4.30924297002105E+0_wp, &
      & 3.47420343563851E+0_wp, 2.82302349370343E+0_wp, 6.67552064739394E+0_wp, &
      & 4.23898159491675E+0_wp]

   call get_structure(mol, "MB16-43", "02")

   allocate(rad(mol%nat), draddr(3, mol%nat, mol%nat))
   rvdw = get_vdw_rad_bondi(mol%num)

   call new_born_integrator(gbobc, mol, rvdw)
   call gbobc%get_rad(mol, rad, draddr)

   if (any(abs(rad - ref) > thr2)) then
      call test_failed(error, "Born radii area values do not match")
      print '(es20.14e1)', rad
      return
   end if

   call test_numg(error, gbobc, mol)

end subroutine test_mb02


subroutine test_mb03(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   real(wp), allocatable :: rvdw(:), rad(:), draddr(:, :, :)
   real(wp), parameter :: ref(16) = [&
      & 4.94764986698701E+0_wp, 3.95122747438791E+0_wp, 4.57075556245289E+0_wp, &
      & 5.46368225070994E+0_wp, 8.24269261139398E+0_wp, 5.68405471762112E+0_wp, &
      & 5.51002325309604E+0_wp, 4.75597020148093E+0_wp, 4.21190195089894E+0_wp, &
      & 4.32836770082885E+0_wp, 4.12684869499911E+0_wp, 5.15171226623248E+0_wp, &
      & 4.83223856996055E+0_wp, 3.02638025720185E+0_wp, 4.05683426506167E+0_wp, &
      & 4.63569783992096E+0_wp]

   call get_structure(mol, "MB16-43", "03")

   allocate(rad(mol%nat), draddr(3, mol%nat, mol%nat))
   rvdw = get_vdw_rad_cosmo(mol%num)

   call new_born_integrator(gbobc, mol, rvdw)
   call gbobc%get_rad(mol, rad, draddr)

   if (any(abs(rad - ref) > thr2)) then
      call test_failed(error, "Born radii area values do not match")
      print '(es20.14e1)', rad
      return
   end if

   call test_numg(error, gbobc, mol)

end subroutine test_mb03


subroutine test_e(error, mol, input, qat, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Solvation model input
   type(alpb_input), intent(in) :: input

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Reference energy
   real(wp), intent(in) :: ref

   type(alpb_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp) :: energy(mol%nat)

   energy = 0.0_wp
   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1))

   solv = alpb_solvation(mol, input)

   call solv%update(mol, cache)
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)

   if (abs(sum(energy) - ref) > thr) then
      call test_failed(error, "Energy does not match reference")
      print *, sum(energy)
   end if
end subroutine test_e


subroutine test_g(error, mol, input, qat)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Solvation model input
   type(alpb_input), intent(in) :: input

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   type(alpb_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp), allocatable :: gradient(:, :), numg(:, :)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3)
   integer :: ii, ic

   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1))

   solv = alpb_solvation(mol, input)

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


subroutine test_p(error, mol, input, qat)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Solvation model input
   type(alpb_input), intent(in) :: input

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   type(alpb_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), parameter :: thr = 1e+3_wp*sqrt(epsilon(1.0_wp))
   real(wp), allocatable :: vat(:)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat)
   integer :: ii

   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1))

   solv = alpb_solvation(mol, input)

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


subroutine test_e_p16(error)

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
   call test_e(error, mol, alpb_input(feps, kernel=born_kernel%p16), qat, &
      & -7.2620663020537416E-3_wp)

end subroutine test_e_p16


subroutine test_e_still(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 2.29208822115185E-1_wp, 4.70816658242009E-2_wp,-3.52834459119718E-2_wp, &
      & 3.09012847396269E-2_wp, 3.28898468854920E-1_wp,-2.01747405019535E-1_wp, &
      & 5.46554362391008E-2_wp,-1.09681283064574E-1_wp,-3.47340505091849E-1_wp, &
      & 1.84567817865267E-1_wp,-2.07552337277185E-1_wp, 4.67140380802351E-1_wp, &
      &-1.84261319200178E-2_wp,-1.05015595324833E-1_wp, 6.52511545312054E-2_wp, &
      &-3.82658324740237E-1_wp]
   real(wp), parameter :: feps = 80.0_wp

   call get_structure(mol, "MB16-43", "05")
   call test_e(error, mol, alpb_input(feps, kernel=born_kernel%still), qat, &
      & -5.7758191311385130E-3_wp)

end subroutine test_e_still


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
   real(wp), parameter :: feps = 80.0_wp

   call get_structure(mol, "UPU23", "0a")
   call test_e(error, mol, alpb_input(feps, kernel=born_kernel%p16, alpb=.true.), qat, &
      & -6.2968900200158134E-2_wp)
   if (allocated(error)) return

   call test_e(error, mol, alpb_input(feps, kernel=born_kernel%p16, alpb=.false.), qat, &
      & -6.3099611865806121E-2_wp)

end subroutine test_e_charged


subroutine test_g_p16(error)

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
   real(wp), parameter :: feps = 80.0_wp

   call get_structure(mol, "MB16-43", "06")
   call test_g(error, mol, alpb_input(feps, kernel=born_kernel%p16), qat)

end subroutine test_g_p16


subroutine test_g_still(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-1.57321098180703E-1_wp, 1.65233008998668E-1_wp, 3.22320267782066E-1_wp, &
      & 3.63564544135336E-2_wp, 4.85639267214320E-2_wp,-3.59203277893926E-1_wp, &
      &-1.93841260011383E-1_wp,-3.86495230324447E-1_wp, 3.10104147485353E-1_wp, &
      & 8.34907519580185E-2_wp,-3.62672063405622E-1_wp, 3.64143595819311E-1_wp, &
      & 3.34640678947868E-1_wp,-4.69881543486815E-1_wp,-1.89222615863620E-1_wp, &
      & 4.53784257040286E-1_wp]
   real(wp), parameter :: feps = 80.0_wp

   call get_structure(mol, "MB16-43", "07")
   call test_g(error, mol, alpb_input(feps, kernel=born_kernel%still), qat)

end subroutine test_g_still


subroutine test_p_p16(error)

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
   real(wp), parameter :: feps = 80.0_wp

   call get_structure(mol, "MB16-43", "08")
   call test_p(error, mol, alpb_input(feps, kernel=born_kernel%p16), qat)

end subroutine test_p_p16


subroutine test_p_still(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-5.25247525508526E-2_wp,-1.97060288484673E-2_wp,-1.36432079012167E-1_wp, &
      &-9.52347486058915E-2_wp, 4.66468767398605E-1_wp, 2.70593750695380E-2_wp, &
      &-1.82819518140298E-1_wp, 1.59893322865668E-1_wp,-1.04573523319519E-1_wp, &
      & 1.16709722834927E-2_wp, 4.11411958411730E-1_wp, 2.87449331948580E-2_wp, &
      &-4.12833246638830E-1_wp,-3.57145165960098E-1_wp,-9.54533576632005E-2_wp, &
      & 3.51473091515422E-1_wp]
   real(wp), parameter :: feps = 80.0_wp

   call get_structure(mol, "MB16-43", "09")
   call test_p(error, mol, alpb_input(feps, kernel=born_kernel%still), qat)

end subroutine test_p_still


end module test_solvation_born
