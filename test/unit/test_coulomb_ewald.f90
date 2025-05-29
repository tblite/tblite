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

module test_coulomb_ewald
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   use mctc_io_constants, only : pi
   use mstore, only : get_structure
   use tblite_blas, only : dot, symv, gemv
   use tblite_cutoff, only : get_lattice_points
   use tblite_coulomb_ewald, only : get_amat_ewald_3d, get_amat_ewald_mp_3d, &
      & get_energy_ewald_3d, get_energy_ewald_mp_3d, get_potential_ewald_3d, &
      & get_alpha, get_rec_cutoff
   implicit none
   private

   public :: collect_coulomb_ewald

   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr_sc = 1e-6_wp
   real(wp), parameter :: conv = sqrt(epsilon(0.0_wp))

contains


!> Collect all exported unit tests
subroutine collect_coulomb_ewald(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("ewald-oxacb", test_ewald_oxacb), &
      new_unittest("ewald-co2", test_ewald_co2), &
      new_unittest("ewald-ammonia", test_ewald_ammonia), &
      new_unittest("ewald-urea", test_ewald_urea), &
      new_unittest("ewald-pyrazine", test_ewald_pyrazine), &
      new_unittest("ewald-mp-co2", test_ewald_mp_co2) &
      ]

end subroutine collect_coulomb_ewald


subroutine test_ewald(error, mol, alpha, qat, e1, e2)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure
   type(structure_type), intent(in) :: mol
   !> Charges
   real(wp), intent(in) :: qat(:)
   !> Ewald energy
   real(wp), intent(out) :: e1, e2

   integer :: ind, iat
   real(wp) :: alpha, vol, rec_lat(3, 3)
   real(wp), allocatable :: rtrans(:, :), amat(:, :), energies(:), vat1(:), vat2(:)
   integer, allocatable :: nshell(:), offset(:)

   allocate(nshell(mol%nat), offset(mol%nat), source=0)
   ind = 0
   do iat = 1, mol%nat
      nshell(iat) = 1
      offset(iat) = ind
      ind = ind + nshell(iat)
   end do

   vol = abs(matdet_3x3(mol%lattice))
   rec_lat = (2*pi)*transpose(matinv_3x3(mol%lattice))
   call get_lattice_points([.true.], rec_lat, get_rec_cutoff(alpha, vol, conv), rtrans)
   rtrans = rtrans(:, 2:)

   allocate(amat(mol%nat, mol%nat), source=0.0_wp)
   call get_amat_ewald_3d(mol, nshell, offset, alpha, vol, rtrans, amat)

   allocate(energies(mol%nat), vat1(mol%nat), vat2(mol%nat), source=0.0_wp)
   call get_energy_ewald_3d(mol, alpha, vol, rtrans, qat, energies)
   call get_potential_ewald_3d(mol, alpha, vol, rtrans, qat, vat1)

   e1 = sum(energies)

   call symv(amat, qat, vat2)
   e2 = 0.5_wp * dot(qat, vat2)
end subroutine test_ewald

subroutine test_ewald_potential(error, mol, alpha, qat)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure
   type(structure_type), intent(in) :: mol
   !> Charges
   real(wp), intent(in) :: qat(:)

   integer :: ind, iat
   real(wp) :: alpha, vol, rec_lat(3, 3), e1, e2
   real(wp), allocatable :: rtrans(:, :), amat(:, :), energies(:), vat1(:), vat2(:)
   real(wp), allocatable :: el(:), er(:), qat_(:)
   integer, allocatable :: nshell(:), offset(:)
   real(wp), parameter :: step = 1e-6_wp

   allocate(nshell(mol%nat), offset(mol%nat), source=0)
   ind = 0
   do iat = 1, mol%nat
      nshell(iat) = 1
      offset(iat) = ind
      ind = ind + nshell(iat)
   end do

   vol = abs(matdet_3x3(mol%lattice))
   rec_lat = (2*pi)*transpose(matinv_3x3(mol%lattice))
   call get_lattice_points([.true.], rec_lat, get_rec_cutoff(alpha, vol, conv), rtrans)
   rtrans = rtrans(:, 2:)

   allocate(amat(mol%nat, mol%nat), source=0.0_wp)
   call get_amat_ewald_3d(mol, nshell, offset, alpha, vol, rtrans, amat)

   allocate(energies(mol%nat), vat1(mol%nat), vat2(mol%nat), source=0.0_wp)
   call get_energy_ewald_3d(mol, alpha, vol, rtrans, qat, energies)
   call get_potential_ewald_3d(mol, alpha, vol, rtrans, qat, vat1)

   e1 = sum(energies)

   call symv(amat, qat, vat2)
   e2 = 0.5_wp * dot(qat, vat2)

   ! do iat = 1, mol%nat
   !    call check(error, vat1(iat), vat2(iat), thr=thr)
   !    if (allocated(error)) then
   !       print *, "Ewald potential test failed for atom ", iat, &
   !          & " with value ", vat1(iat), " and ", vat2(iat), &
   !          & " difference ", vat1(iat) - vat2(iat)
   !       return
   !    end if
   ! end do

   allocate(er(mol%nat), el(mol%nat), source=0.0_wp)
   qat_ = qat
   do iat = 1, mol%nat
      qat_(iat) = qat(iat) + step
      er(:) = 0.0_wp
      call get_energy_ewald_3d(mol, alpha, vol, rtrans, qat_, er)
      qat_(iat) = qat(iat) - step
      el(:) = 0.0_wp
      call get_energy_ewald_3d(mol, alpha, vol, rtrans, qat_, el)
      vat2(iat) = (er(iat) - el(iat)) / (2 * step)
   end do

   do iat = 1, mol%nat
      call check(error, vat1(iat), vat2(iat), thr=thr)
      if (allocated(error)) then
         print *, "Ewald potential numdiff failed for atom ", iat, &
            & " with value ", vat1(iat), " and ", vat2(iat), &
            & " difference ", vat1(iat) - vat2(iat)
         return
      end if
   end do
end subroutine test_ewald_potential

subroutine test_generic(error, mol, qat, reference)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure
   type(structure_type), intent(in) :: mol
   !> Charges
   real(wp), intent(in) :: qat(:)
   !> Reference value
   real(wp), intent(in) :: reference

   integer, parameter :: supercell(3) = [2, 2, 2]
   type(structure_type) :: mol_sc
   real(wp), allocatable :: qat_sc(:)
   real(wp) :: e1, e2, e1_sc, e2_sc, alpha

   call get_alpha(mol%lattice, alpha, .false.)

   call test_ewald(error, mol, alpha, qat, e1, e2)
   ! print *, "Ewald energy: ", e1, e2, e1 - e2
   call check(error, e1, e2, thr=thr)
   if (allocated(error)) return

   qat_sc = [spread(qat, 2, product(supercell))]
   call make_supercell(mol, mol_sc, supercell)
   call test_ewald(error, mol_sc, alpha, qat_sc, e1_sc, e2_sc)
   ! print *, "Ewald energy in supercell: ", e1_sc / product(supercell), &
   !    & e2_sc / product(supercell), e1_sc / product(supercell) - e2_sc / product(supercell)
   call check(error, e1_sc, e2_sc, thr=thr * product(supercell))
   if (allocated(error)) return

   ! print *, "Super cell scaling: ", e1, e1_sc / product(supercell), &
   !    & e1 - e1_sc / product(supercell)
   call check(error, e1, e1_sc / product(supercell), thr=thr_sc)
   if (allocated(error)) return
   ! print *, "Super cell scaling: ", e2, e2_sc / product(supercell), &
   !    & e2 - e2_sc / product(supercell)
   call check(error, e2, e2_sc / product(supercell), thr=thr_sc)
   if (allocated(error)) return

   call check(error, e1, reference, thr=thr)
   if (allocated(error)) return
   ! print *, "Ewald energy: ", e1, reference, e1 - reference

   ! call test_ewald_potential(error, mol, alpha, qat)
end subroutine test_generic

subroutine test_ewald_oxacb(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 3.41731844312030E-1_wp, 3.41716020106239E-1_wp, 3.41730526585671E-1_wp,&
      & 3.41714427217954E-1_wp, 3.80996046757999E-1_wp, 3.80989821246195E-1_wp,&
      & 3.81000747720282E-1_wp, 3.80990494183703E-1_wp,-3.70406587264474E-1_wp,&
      &-3.70407565207006E-1_wp,-3.70417590212352E-1_wp,-3.70399716470705E-1_wp,&
      &-3.52322260586075E-1_wp,-3.52304269439196E-1_wp,-3.52313440903261E-1_wp,&
      &-3.52298498047004E-1_wp]

   call get_structure(mol, "X23", "oxacb")
   call test_generic(error, mol, qat, -0.25622014078408595_wp)

end subroutine test_ewald_oxacb

subroutine test_ewald_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 4.56275672862067E-1_wp, 4.56284770386671E-1_wp, 4.56284770386671E-1_wp,&
      & 4.56284770386671E-1_wp,-2.28127680925611E-1_wp,-2.28138283131909E-1_wp,&
      &-2.28145770512561E-1_wp,-2.28145770512561E-1_wp,-2.28150142163058E-1_wp,&
      &-2.28145770512561E-1_wp,-2.28138283131909E-1_wp,-2.28138283131909E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "CO2")
   call test_generic(error, mol, qat, -0.20685735924532045_wp)

end subroutine test_ewald_co2


subroutine test_ewald_urea(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 5.55723890858218E-1_wp, 5.55765354442035E-1_wp, 2.50200231242017E-1_wp,&
      & 2.50282053284422E-1_wp, 2.39786980460652E-1_wp, 2.39895142481200E-1_wp,&
      & 2.50103678240412E-1_wp, 2.50425041601730E-1_wp, 2.39464477136495E-1_wp,&
      & 2.40360053062669E-1_wp,-4.38369096728919E-1_wp,-4.38451412936599E-1_wp,&
      &-4.38310020776279E-1_wp,-4.38617373848238E-1_wp,-6.59141030224988E-1_wp,&
      &-6.59117968294813E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "urea")
   call test_generic(error, mol, qat, -0.43822886730428179_wp)

end subroutine test_ewald_urea

subroutine test_ewald_ammonia(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 2.95376975876519E-1_wp, 2.95376975876519E-1_wp, 2.95376975876519E-1_wp,&
      & 2.95329109335847E-1_wp, 2.95332441468412E-1_wp, 2.95347202855778E-1_wp,&
      & 2.95347202855779E-1_wp, 2.95329109335848E-1_wp, 2.95332441468411E-1_wp,&
      & 2.95347202855777E-1_wp, 2.95329109335847E-1_wp, 2.95332441468412E-1_wp,&
      &-8.86118742099358E-1_wp,-8.86012815503436E-1_wp,-8.86012815503437E-1_wp,&
      &-8.86012815503434E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "ammonia")
   call test_generic(error, mol, qat, -0.71791185165457716_wp)

end subroutine test_ewald_ammonia


subroutine test_ewald_pyrazine(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 1.23705026512437E-1_wp, 1.22537765989959E-1_wp, 1.23932626831231E-1_wp,&
      & 1.22418959326822E-1_wp, 1.23788033684569E-1_wp, 1.23643389058068E-1_wp,&
      & 1.22389811551880E-1_wp, 1.22402837718155E-1_wp, 3.17174594622166E-2_wp,&
      & 3.15585817789390E-2_wp, 3.14290005061543E-2_wp, 3.10358506314526E-2_wp,&
      & 3.11084225356749E-2_wp, 3.11187325528499E-2_wp, 3.10707965296438E-2_wp,&
      & 3.13508824430484E-2_wp,-3.09107070033859E-1_wp,-3.09144123845085E-1_wp,&
      &-3.08473365847409E-1_wp,-3.08483617386745E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "pyrazine")
   call test_generic(error, mol, qat, -5.7642313568676787E-2_wp)

end subroutine test_ewald_pyrazine

subroutine test_ewald_mp(error, mol, alpha, qat, dpat, qpat, e1, e2)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure
   type(structure_type), intent(in) :: mol
   !> Charges
   real(wp), intent(in) :: qat(:)
   !> Dipole moments
   real(wp), intent(in) :: dpat(:, :)
   !> Quadrupole moments
   real(wp), intent(in) :: qpat(:, :)
   !> Ewald energy
   real(wp), intent(out) :: e1, e2

   integer :: ind, iat
   real(wp) :: alpha, vol, rec_lat(3, 3)
   real(wp), allocatable :: amat_sd(:, :, :), amat_dd(:, :, :, :), amat_sq(:, :, :)
   real(wp), allocatable :: rtrans(:, :), energies(:)
   real(wp), allocatable :: vat1(:), vat2(:), vdpat1(:, :), vdpat2(:, :), vqpat1(:, :), vqpat2(:, :)

   vol = abs(matdet_3x3(mol%lattice))
   rec_lat = (2*pi)*transpose(matinv_3x3(mol%lattice))
   call get_lattice_points([.true.], rec_lat, get_rec_cutoff(alpha, vol, conv), rtrans)
   rtrans = rtrans(:, 2:)

   allocate(amat_sd(3, mol%nat, mol%nat), source=0.0_wp)
   allocate(amat_dd(3, mol%nat, 3, mol%nat), source=0.0_wp)
   allocate(amat_sq(6, mol%nat, mol%nat), source=0.0_wp)
   call get_amat_ewald_mp_3d(mol, alpha, vol, rtrans, amat_sd, amat_dd, amat_sq)

   allocate(energies(mol%nat), vat1(mol%nat), vat2(mol%nat), vdpat1(3, mol%nat), vdpat2(3, mol%nat), &
      & vqpat1(6, mol%nat), vqpat2(6, mol%nat), source=0.0_wp)
   call get_energy_ewald_mp_3d(mol, alpha, vol, rtrans, qat, dpat, qpat, energies)
   ! call get_potential_ewald_3d(mol, alpha, vol, rtrans, qat, vat1)

   e1 = sum(energies)

   call gemv(amat_sd, qat, vdpat2)
   call gemv(amat_dd, dpat, vdpat2, beta=1.0_wp, alpha=0.5_wp)
   call gemv(amat_sq, qat, vqpat2)

   e2 = dot(dpat, vdpat2) + dot(qpat, vqpat2)
end subroutine test_ewald_mp

subroutine test_generic_mp(error, mol, qat, dpat, qpat, reference)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure
   type(structure_type), intent(in) :: mol
   !> Charges
   real(wp), intent(in) :: qat(:)
   !> Dipole moments
   real(wp), intent(in) :: dpat(:, :)
   !> Quadrupole moments
   real(wp), intent(in) :: qpat(:, :)
   !> Reference value
   real(wp), intent(in) :: reference

   integer, parameter :: supercell(3) = [2, 2, 2]
   type(structure_type) :: mol_sc
   real(wp) :: e1, e2, e1_sc, e2_sc, alpha
   real(wp), allocatable :: qat0(:), dpat0(:, :), qpat0(:, :)
   real(wp), allocatable :: qat_sc(:), dpat_sc(:, :), qpat_sc(:, :)

   allocate(qat0(size(qat)), dpat0(3, size(dpat, 2)), qpat0(6, size(qpat, 2)), source=0.0_wp)

   call get_alpha(mol%lattice, alpha, .false.)

   call test_ewald_mp(error, mol, alpha, qat0, dpat, qpat0, e1, e2)
   call check(error, e1, e2, thr=thr)
   if (allocated(error)) then
      print *, "Ewald energy (dd): ", e1, e2, e1 - e2
      return
   end if

   call test_ewald_mp(error, mol, alpha, qat, dpat0, qpat, e1, e2)
   call check(error, e1, e2, thr=thr)
   if (allocated(error)) then
      print *, "Ewald energy (sq): ", e1, e2, e1 - e2
      return
   end if

   call test_ewald_mp(error, mol, alpha, qat, dpat, qpat0, e1, e2)
   call check(error, e1, e2, thr=thr)
   if (allocated(error)) then
      print *, "Ewald energy (sd): ", e1, e2, e1 - e2
      return
   end if

   call test_ewald_mp(error, mol, alpha, qat, dpat, qpat, e1, e2)
   call check(error, e1, e2, thr=thr)
   if (allocated(error)) then
      print *, "Ewald energy (all): ", e1, e2, e1 - e2
      return
   end if

   deallocate(qat0, dpat0, qpat0)

   qat_sc = [spread(qat, 2, product(supercell))]
   dpat_sc = reshape([dpat, dpat, dpat, dpat, dpat, dpat, dpat, dpat, dpat], &
      & [3, product(supercell)*mol%nat])
   qpat_sc = reshape([qpat, qpat, qpat, qpat, qpat, qpat, qpat, qpat, qpat], &
      & [6, product(supercell)*mol%nat])
   allocate(qat0(size(qat_sc)), dpat0(3, size(dpat_sc, 2)), qpat0(6, size(qpat_sc, 2)), source=0.0_wp)
   call make_supercell(mol, mol_sc, supercell)

   call test_ewald_mp(error, mol_sc, alpha, qat_sc, dpat_sc, qpat_sc, e1_sc, e2_sc)
   call check(error, e1_sc, e2_sc, thr=thr * product(supercell))
   if (allocated(error)) then
      print *, "Ewald energy in supercell: ", e1_sc, e2_sc, e1_sc - e2_sc
      return
   end if

   call check(error, e1, e1_sc / product(supercell), thr=thr_sc)
   if (allocated(error)) then
      print *, "Super cell scaling: ", e1, e1_sc / product(supercell), &
         & e1 - e1_sc / product(supercell)
      return
   end if
   call check(error, e2, e2_sc / product(supercell), thr=thr_sc)
   if (allocated(error)) then
      print *, "Super cell scaling: ", e2, e2_sc / product(supercell), &
         & e2 - e2_sc / product(supercell)
      return
   end if

   call check(error, e1, reference, thr=thr)
   if (allocated(error)) then
      print *, "Ewald energy: ", e1, reference, e1 - reference
      return
   end if
end subroutine test_generic_mp

subroutine test_ewald_mp_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 12
   type(structure_type) :: mol
   real(wp), parameter :: qat(nat) = [&
      &  4.95105332967126E-01_wp,  4.95110445149787E-01_wp,  4.95109208803526E-01_wp, &
      &  4.95110553060372E-01_wp, -2.47570208775520E-01_wp, -2.47372219387123E-01_wp, &
      & -2.47367867558844E-01_wp, -2.47541478966424E-01_wp, -2.47535153228645E-01_wp, &
      & -2.47738212772185E-01_wp, -2.47741337427010E-01_wp, -2.47569061865040E-01_wp]
   real(wp), parameter :: dpat(3, nat) = reshape([&
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
      & [3, nat])
   real(wp), parameter :: qpat(6, nat) = reshape([&
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
      & [6, nat])

   call get_structure(mol, "X23", "CO2")
   call test_generic_mp(error, mol, qat, dpat, qpat, 7.4194473613609238E-3_wp)

end subroutine test_ewald_mp_co2

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


end module test_coulomb_ewald
