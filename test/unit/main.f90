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

!> Driver for unit testing
program tester
   use, intrinsic :: iso_fortran_env, only : error_unit
   use mctc_env, only : get_argument
   use mctc_env_testing, only : run_testsuite, new_testsuite, testsuite_type, &
      & select_suite, run_selected
   use test_cgto_ortho, only : collect_cgto_ortho
   use test_coulomb_charge, only : collect_coulomb_charge
   use test_coulomb_multipole, only : collect_coulomb_multipole
   use test_fit, only : collect_fit
   use test_gfn1_xtb, only : collect_gfn1_xtb
   use test_gfn2_xtb, only : collect_gfn2_xtb
   use test_hamiltonian, only : collect_hamiltonian
   use test_halogen, only : collect_halogen
   use test_integral_multipole, only : collect_integral_multipole
   use test_integral_overlap, only : collect_integral_overlap
   use test_ipea1_xtb, only : collect_ipea1_xtb
   use test_ncoord_gfn, only : collect_ncoord_gfn
   use test_repulsion, only : collect_repulsion
   use test_slater_expansion, only : collect_slater_expansion
   use test_spin, only : collect_spin
   use test_solvation_born, only : collect_solvation_born
   use test_solvation_cpcm, only : collect_solvation_cpcm
   use test_solvation_surface, only : collect_solvation_surface
   use test_tagged_io, only : collect_tagged_io
   use test_xtb_external, only : collect_xtb_external
   use test_xtb_param, only : collect_xtb_param
   implicit none
   integer :: stat, is
   character(len=:), allocatable :: suite_name, test_name
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
      new_testsuite("tagged-io", collect_tagged_io), &
      new_testsuite("fit", collect_fit), &
      new_testsuite("repulsion", collect_repulsion), &
      new_testsuite("ncoord-gfn", collect_ncoord_gfn), &
      new_testsuite("solvation-born", collect_solvation_born), &
      new_testsuite("solvation-cpcm", collect_solvation_cpcm), &
      new_testsuite("solvation-surface", collect_solvation_surface), &
      new_testsuite("coulomb-charge", collect_coulomb_charge), &
      new_testsuite("coulomb-multipole", collect_coulomb_multipole), &
      new_testsuite("slater-expansion", collect_slater_expansion), &
      new_testsuite("cgto-ortho", collect_cgto_ortho), &
      new_testsuite("integral-overlap", collect_integral_overlap), &
      new_testsuite("integral-multipole", collect_integral_multipole), &
      new_testsuite("hamiltonian", collect_hamiltonian), &
      new_testsuite("halogen", collect_halogen), &
      new_testsuite("gfn1-xtb", collect_gfn1_xtb), &
      new_testsuite("ipea1-xtb", collect_ipea1_xtb), &
      new_testsuite("gfn2-xtb", collect_gfn2_xtb), &
      new_testsuite("xtb-external", collect_xtb_external), &
      new_testsuite("spin", collect_spin), &
      new_testsuite("xtb-param", collect_xtb_param) &
      ]

   call get_argument(1, suite_name)
   call get_argument(2, test_name)

   if (allocated(suite_name)) then
      is = select_suite(testsuites, suite_name)
      if (is > 0 .and. is <= size(testsuites)) then
         if (allocated(test_name)) then
            write(error_unit, fmt) "Suite:", testsuites(is)%name
            call run_selected(testsuites(is)%collect, test_name, error_unit, stat)
            if (stat < 0) then
               error stop 1
            end if
         else
            write(error_unit, fmt) "Testing:", testsuites(is)%name
            call run_testsuite(testsuites(is)%collect, error_unit, stat)
         end if
      else
         write(error_unit, fmt) "Available testsuites"
         do is = 1, size(testsuites)
            write(error_unit, fmt) "-", testsuites(is)%name
         end do
         error stop 1
      end if
   else
      do is = 1, size(testsuites)
         write(error_unit, fmt) "Testing:", testsuites(is)%name
         call run_testsuite(testsuites(is)%collect, error_unit, stat)
      end do
   end if

   if (stat > 0) then
      write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop 1
   end if


end program tester
