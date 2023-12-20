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

!> @dir tblite/ceh/ceh.f90
!> Contains the implementation of the Charge Extended Hückel (CEH) method.

module tblite_ceh_ceh
   !> mctc-lib
   use mctc_env, only : error_type, wp
   use mctc_io, only: structure_type
   !> Basis set
   use tblite_basis_ortho, only : orthogonalize
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_basis_type, only : cgto_type, new_basis, get_cutoff, basis_type
   !> Coordination number
   use tblite_ncoord, only : new_ncoord
   !> Calculation context
   use tblite_context, only : context_type
   use tblite_output_format, only: format_string
   !> Integrals
   use tblite_integral_type, only : integral_type, new_integral
   use tblite_integral_dipole, only: get_dipole_integrals, dipole_cgto, &
   & dipole_cgto_diat_scal, maxl, msao
   use tblite_integral_diat_trafo, only: relvec
   !> Wavefunction
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction, &
   & get_alpha_beta_occupation
   use tblite_wavefunction_mulliken, only: get_mulliken_shell_charges, &
   & get_mulliken_atomic_multipoles
   use tblite_scf_iterator, only: get_density, get_qat_from_qsh
   !> Additional potentials (external field)
   use tblite_scf, only: new_potential, potential_type
   use tblite_external_field, only : electric_field
   use tblite_container, only : container_type, container_cache
   use tblite_scf_potential, only: add_pot_to_h1
   !> Electronic solver
   use tblite_lapack_solver, only: lapack_solver
   use tblite_scf_solver, only : solver_type
   !> BLAS
   use tblite_blas, only: gemv
   !> CEH specific
   use tblite_ceh_calculator, only : ceh_calculator
   use tblite_ceh_h0, only : ceh_hamiltonian
   !> Miscelaneous
   use tblite_timer, only : timer_type, format_time

   implicit none
   private

   public :: ceh_guess, new_ceh_calculator

   integer, parameter, private :: max_elem = 86
   integer, parameter, private :: max_shell = 3

   !> Number of shells # MM, August 01, 2023
   integer, parameter :: nshell(max_elem) = [&
   & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 2, 3, & ! 1-20
   & 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, & ! 21-40
   & 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, & ! 41-60
   & 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, & ! 61-80
   & 3, 3, 3, 3, 3, 3]                                             ! 81-86

   !> Angular momentum of each shell # MM, August 01, 2023
   ! 0 = s, 1 = p, 2 = d
   ! CAUTION: Ordering from original CEH model is taken for consistency with the parameterization
   ! I.e., the ordering of the shells is always: "s", "p", "d"
   integer, parameter :: ang_shell(max_shell, max_elem) = reshape([&
   & 0, 0, 0,  0, 0, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0, & ! 1-7
   & 0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 2,  0, 1, 2, & ! 8-14
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 0,  0, 1, 2,  0, 1, 2, & ! 15-21
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 22-28
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 29-35
   & 0, 1, 2,  0, 1, 0,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 36-42
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 43-49
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 0,  0, 1, 2, & ! 50-56
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 57-63
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 64-70
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 71-77
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 78-84
   & 0, 1, 2,  0, 1, 2], shape(ang_shell))                                  ! 85-86

   !> Principal quantum number of each shell (see commeent above regarding angular momentum)
   integer, parameter :: principal_quantum_number(max_shell, max_elem) = reshape([&
   & 1, 2, 0,  1, 0, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0, & ! 1-7
   & 2, 2, 0,  2, 2, 0,  2, 2, 3,  3, 3, 0,  3, 3, 0,  3, 3, 3,  3, 3, 3, & ! 8-14
   & 3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  4, 4, 0,  4, 4, 3,  4, 4, 3, & ! 15-21
   & 4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3, & ! 22-28
   & 4, 4, 3,  4, 4, 3,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4, & ! 29-35
   & 4, 4, 4,  5, 5, 0,  5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 4, & ! 36-42
   & 5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 5, & ! 43-49
   & 5, 5, 5,  5, 5, 5,  5, 5, 5,  5, 5, 5,  5, 5, 5,  6, 6, 0,  6, 6, 5, & ! 50-56
   & 6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5, & ! 57-63
   & 6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5, & ! 64-70
   & 6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5, & ! 71-77
   & 6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5, & ! 78-84
   & 6, 6, 5,  6, 6, 5], shape(principal_quantum_number))                   ! 85-86

   !> Number of primitive gaussians per shell # MM, August 01, 2023
   integer, parameter :: number_of_primitives(max_shell, max_elem) = reshape([&
   & 4, 0, 0,  4, 0, 0,  4, 4, 0,  4, 4, 0,  4, 4, 0,  4, 4, 0,  4, 4, 0, & ! 1-7
   & 4, 4, 0,  4, 4, 0,  4, 4, 0,  4, 4, 0,  4, 4, 0,  4, 4, 4,  4, 4, 4, & ! 8-14
   & 4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 0,  4, 4, 4,  4, 4, 4, & ! 15-21
   & 4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4, & ! 22-28
   & 4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4, & ! 29-35
   & 4, 4, 4,  4, 4, 0,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4, & ! 36-42
   & 4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4, & ! 43-49
   & 4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  6, 6, 0,  6, 6, 4, & ! 50-56
   & 6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4, & ! 57-63
   & 6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4, & ! 64-70
   & 6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4, & ! 71-77
   & 6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4, & ! 78-84
   & 6, 6, 4,  6, 6, 4], shape(number_of_primitives))                       ! 85-86

   !> Reference occupation of the atom
   real(wp), parameter :: reference_occ(max_shell, max_elem) = reshape([&
   & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  1.0_wp, 0.0_wp, 0.0_wp, & ! 1-3
   & 2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 1.0_wp, 0.0_wp,  2.0_wp, 2.0_wp, 0.0_wp, & ! 4-6
   & 2.0_wp, 3.0_wp, 0.0_wp,  2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp, & ! 7-9
   & 2.0_wp, 6.0_wp, 0.0_wp,  1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, & ! 10-12
   & 2.0_wp, 1.0_wp, 0.0_wp,  2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp, & ! 13-15
   & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, & ! 16-18
   & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, & ! 19-21
   & 2.0_wp, 0.0_wp, 2.0_wp,  2.0_wp, 0.0_wp, 3.0_wp,  2.0_wp, 0.0_wp, 4.0_wp, & ! 22-24
   & 2.0_wp, 0.0_wp, 5.0_wp,  2.0_wp, 0.0_wp, 6.0_wp,  2.0_wp, 0.0_wp, 7.0_wp, & ! 25-27
   & 2.0_wp, 0.0_wp, 8.0_wp,  2.0_wp, 0.0_wp, 9.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, & ! 28-30
   & 2.0_wp, 1.0_wp, 0.0_wp,  2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp, & ! 31-33
   & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, & ! 34-36
   & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, & ! 37-39
   & 2.0_wp, 0.0_wp, 2.0_wp,  2.0_wp, 0.0_wp, 3.0_wp,  2.0_wp, 0.0_wp, 4.0_wp, & ! 40-42
   & 2.0_wp, 0.0_wp, 5.0_wp,  2.0_wp, 0.0_wp, 6.0_wp,  2.0_wp, 0.0_wp, 7.0_wp, & ! 43-45
   & 2.0_wp, 0.0_wp, 8.0_wp,  2.0_wp, 0.0_wp, 9.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, & ! 46-48
   & 2.0_wp, 1.0_wp, 0.0_wp,  2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp, & ! 49-51
   & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, & ! 52-54
   & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, & ! 55-57
   & 2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, & ! 58-60
   & 2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, & ! 61-63
   & 2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, & ! 64-66
   & 2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, & ! 67-69
   & 2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 2.0_wp, & ! 70-72
   & 2.0_wp, 0.0_wp, 3.0_wp,  2.0_wp, 0.0_wp, 4.0_wp,  2.0_wp, 0.0_wp, 5.0_wp, & ! 73-75
   & 2.0_wp, 0.0_wp, 6.0_wp,  2.0_wp, 0.0_wp, 7.0_wp,  2.0_wp, 0.0_wp, 8.0_wp, & ! 76-78
   & 2.0_wp, 0.0_wp, 9.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 1.0_wp, 0.0_wp, & ! 79-81
   & 2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp,  2.0_wp, 4.0_wp, 0.0_wp, & ! 82-84
   & 2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp], shape(reference_occ))     ! 85-86

   !> Exponent of the Slater function # MM, August 00, 2023
   real(wp), parameter :: slater_exponent(max_shell, max_elem) = reshape([&
   &  1.13056314_wp,  0.00000000_wp,  0.00000000_wp,  1.50065094_wp,  0.00000000_wp,  0.00000000_wp, &
   &  1.19413208_wp,  1.00435986_wp,  0.00000000_wp,  1.04160365_wp,  1.58776712_wp,  0.00000000_wp, &
   &  1.80615826_wp,  1.68397018_wp,  0.00000000_wp,  1.92578557_wp,  1.64086778_wp,  0.00000000_wp, &
   &  2.16223922_wp,  1.79669019_wp,  0.00000000_wp,  2.31036294_wp,  2.07689830_wp,  0.00000000_wp, &
   &  2.63689506_wp,  2.39472851_wp,  0.00000000_wp,  3.02968365_wp,  2.22048798_wp,  0.00000000_wp, &
   &  1.75555631_wp,  1.33838074_wp,  0.00000000_wp,  1.70435065_wp,  1.39174584_wp,  0.00000000_wp, &
   &  1.80712664_wp,  1.62636355_wp,  2.28258151_wp,  2.02930097_wp,  1.63990927_wp,  1.85910992_wp, &
   &  2.24859618_wp,  1.77685356_wp,  1.76755789_wp,  2.72470838_wp,  2.01797698_wp,  2.04253078_wp, &
   &  2.79738729_wp,  2.29520297_wp,  1.61932628_wp,  2.53272473_wp,  1.81071325_wp,  1.69070808_wp, &
   &  1.98719413_wp,  1.22379806_wp,  0.00000000_wp,  1.62021149_wp,  1.54648435_wp,  2.66580857_wp, &
   &  2.14684072_wp,  1.36068696_wp,  1.74878015_wp,  0.69639491_wp,  1.03767289_wp,  1.86366428_wp, &
   &  1.69851374_wp,  2.17779533_wp,  1.83191727_wp,  1.12031062_wp,  1.41500509_wp,  2.16118253_wp, &
   &  1.76160033_wp,  1.98316214_wp,  2.06245601_wp,  1.92205103_wp,  1.30492118_wp,  2.50946021_wp, &
   &  2.07573975_wp,  1.46306532_wp,  2.63099323_wp,  2.44652113_wp,  1.68710423_wp,  3.04820511_wp, &
   &  1.80368223_wp,  1.71760742_wp,  3.02263698_wp,  1.92342980_wp,  1.82453676_wp,  3.49908671_wp, &
   &  2.24151618_wp,  2.08921540_wp,  2.12567962_wp,  2.23205133_wp,  1.89430927_wp,  1.97513779_wp, &
   &  2.69088900_wp,  2.12267272_wp,  1.94145540_wp,  2.94651193_wp,  2.35079302_wp,  1.79682240_wp, &
   &  2.62665945_wp,  2.50566210_wp,  1.75613478_wp,  2.72057897_wp,  2.16290669_wp,  1.88632737_wp, &
   &  1.70116961_wp,  1.64545292_wp,  0.00000000_wp,  1.58807973_wp,  1.80661603_wp,  2.99671810_wp, &
   &  1.08395599_wp,  1.42679219_wp,  2.19662852_wp,  1.47634306_wp,  1.29489380_wp,  2.37913368_wp, &
   &  2.69077953_wp,  1.95585989_wp,  2.11756166_wp,  2.44347621_wp,  1.73442204_wp,  2.02065708_wp, &
   &  2.23908639_wp,  1.95016590_wp,  2.18967157_wp,  2.40381905_wp,  1.44282665_wp,  2.59962263_wp, &
   &  1.96295898_wp,  1.58471984_wp,  2.95895238_wp,  2.38179585_wp,  2.23204400_wp,  3.22819136_wp, &
   &  2.51271048_wp,  1.85519318_wp,  3.24389973_wp,  2.11457529_wp,  1.96275823_wp,  3.51477490_wp, &
   &  2.38640874_wp,  2.12255171_wp,  2.27363232_wp,  2.62089386_wp,  2.06375224_wp,  2.75185414_wp, &
   &  2.79382121_wp,  2.28410073_wp,  2.04070104_wp,  2.97205755_wp,  2.32287719_wp,  1.74734495_wp, &
   &  2.72692370_wp,  2.65571946_wp,  1.85299047_wp,  2.73000797_wp,  2.51650283_wp,  2.24007230_wp, &
   &  1.68663815_wp,  2.16354335_wp,  0.00000000_wp,  1.71010093_wp,  1.58444852_wp,  3.02035513_wp, &
   &  1.11292580_wp,  1.67649145_wp,  2.19394682_wp,  2.04931361_wp,  2.02737301_wp,  2.57681508_wp, &
   &  2.01436702_wp,  2.02975950_wp,  2.59448588_wp,  1.97942043_wp,  2.03214600_wp,  2.61215668_wp, &
   &  1.94447383_wp,  2.03453249_wp,  2.62982748_wp,  1.90952724_wp,  2.03691898_wp,  2.64749828_wp, &
   &  1.87458065_wp,  2.03930547_wp,  2.66516908_wp,  1.83963406_wp,  2.04169196_wp,  2.68283988_wp, &
   &  1.80468747_wp,  2.04407845_wp,  2.70051068_wp,  1.76974088_wp,  2.04646494_wp,  2.71818148_wp, &
   &  1.73479429_wp,  2.04885143_wp,  2.73585228_wp,  1.69984770_wp,  2.05123792_wp,  2.75352308_wp, &
   &  1.66490111_wp,  2.05362441_wp,  2.77119388_wp,  1.62995451_wp,  2.05601091_wp,  2.78886467_wp, &
   &  1.59500792_wp,  2.05839740_wp,  2.80653547_wp,  1.23998866_wp,  1.55902021_wp,  2.56002100_wp, &
   &  2.21476417_wp,  1.21067246_wp,  2.27168445_wp,  2.35183881_wp,  2.17165636_wp,  2.38705146_wp, &
   &  2.29642436_wp,  2.38947110_wp,  2.49180058_wp,  2.78850456_wp,  1.82182727_wp,  3.11161203_wp, &
   &  2.18906814_wp,  1.84759137_wp,  3.45359549_wp,  2.08875168_wp,  2.56532973_wp,  3.59828518_wp, &
   &  2.39586959_wp,  2.60766526_wp,  3.87307779_wp,  2.39636916_wp,  2.53950059_wp,  3.95389569_wp, &
   &  2.72134635_wp,  2.34645674_wp,  2.67736087_wp,  2.99323675_wp,  2.17128412_wp,  5.38232746_wp, &
   &  3.05653080_wp,  2.49075153_wp,  2.67638930_wp,  3.02905756_wp,  2.63479560_wp,  1.94913437_wp, &
   &  2.54746694_wp,  2.83550170_wp,  1.88029428_wp,  2.26386287_wp,  2.46706218_wp,  2.09966650_wp],&
   & shape(slater_exponent))

   !> Atomic orbital shell energy level # MM, August 00, 2023
   real(wp), parameter :: ceh_level(max_shell, max_elem) = reshape([&
   & -0.50000000_wp,  0.00000000_wp,  0.00000000_wp, -0.57125723_wp,  0.00000000_wp,  0.00000000_wp, &
   & -0.39565609_wp, -0.16907309_wp,  0.00000000_wp, -0.46895051_wp, -0.34955238_wp,  0.00000000_wp, &
   & -0.59747672_wp, -0.37636672_wp,  0.00000000_wp, -0.61483866_wp, -0.39087888_wp,  0.00000000_wp, &
   & -0.54872469_wp, -0.43093840_wp,  0.00000000_wp, -0.54699351_wp, -0.42145348_wp,  0.00000000_wp, &
   & -0.58434026_wp, -0.46052654_wp,  0.00000000_wp, -0.61062817_wp, -0.55100942_wp,  0.00000000_wp, &
   & -0.41788970_wp, -0.17987754_wp,  0.00000000_wp, -0.47886929_wp, -0.22957173_wp,  0.00000000_wp, &
   & -0.55132058_wp, -0.35311397_wp, -0.12176082_wp, -0.56321172_wp, -0.40581247_wp, -0.11309926_wp, &
   & -0.66917833_wp, -0.42673915_wp, -0.10120352_wp, -0.53720655_wp, -0.44594494_wp, -0.14447503_wp, &
   & -0.69824334_wp, -0.46955037_wp, -0.08507092_wp, -0.45558942_wp, -0.58263724_wp, -0.15380014_wp, &
   & -0.45067129_wp, -0.21425082_wp,  0.00000000_wp, -0.45071787_wp, -0.24715593_wp, -0.36116434_wp, &
   & -0.41411705_wp, -0.09023659_wp, -0.39582033_wp, -0.49489455_wp, -0.12171921_wp, -0.42870973_wp, &
   & -0.26775746_wp,  0.02131362_wp, -0.47682296_wp, -0.30698882_wp, -0.08166219_wp, -0.42886348_wp, &
   & -0.35893667_wp, -0.15681968_wp, -0.46796329_wp, -0.42189488_wp, -0.15158595_wp, -0.44515883_wp, &
   & -0.26304867_wp, -0.28077502_wp, -0.44505112_wp, -0.40887861_wp, -0.09411084_wp, -0.50469491_wp, &
   & -0.45774194_wp, -0.19313602_wp, -0.58014645_wp, -0.50481797_wp, -0.24203976_wp, -0.48025570_wp, &
   & -0.56565639_wp, -0.37817956_wp, -0.06106926_wp, -0.58170224_wp, -0.42023981_wp, -0.21861265_wp, &
   & -0.59773334_wp, -0.44273372_wp, -0.15372635_wp, -0.60341440_wp, -0.44377998_wp, -0.13702676_wp, &
   & -0.75228843_wp, -0.46070782_wp, -0.16311997_wp, -0.45400718_wp, -0.55540080_wp, -0.29627321_wp, &
   & -0.43480220_wp, -0.25268935_wp,  0.00000000_wp, -0.44963509_wp, -0.28690023_wp, -0.37294037_wp, &
   & -0.43635740_wp, -0.15301650_wp, -0.42985795_wp, -0.49887749_wp, -0.12526225_wp, -0.41994890_wp, &
   & -0.33637242_wp,  0.11506145_wp, -0.47063236_wp, -0.30236908_wp, -0.21307280_wp, -0.44059807_wp, &
   & -0.41053189_wp, -0.08929977_wp, -0.48059019_wp, -0.26112615_wp, -0.21074503_wp, -0.46642678_wp, &
   & -0.37000129_wp, -0.23755777_wp, -0.43199592_wp, -0.38974216_wp, -0.06650347_wp, -0.51904336_wp, &
   & -0.46812094_wp, -0.20754727_wp, -0.55921606_wp, -0.50134195_wp, -0.25122862_wp, -0.80255009_wp, &
   & -0.55190032_wp, -0.36279281_wp,  0.01756359_wp, -0.58038815_wp, -0.40263314_wp, -0.18395644_wp, &
   & -0.65306924_wp, -0.43960193_wp, -0.13809004_wp, -0.54509849_wp, -0.44924670_wp, -0.16182853_wp, &
   & -0.95593684_wp, -0.45693729_wp, -0.16969031_wp, -0.62532339_wp, -0.50865282_wp, -0.34806027_wp, &
   & -0.39616051_wp, -0.26832346_wp,  0.00000000_wp, -0.44826735_wp, -0.21976536_wp, -0.37260694_wp, &
   & -0.43176645_wp, -0.16248487_wp, -0.40748316_wp, -0.25060286_wp, -0.27387264_wp, -0.45115761_wp, &
   & -0.23930086_wp, -0.27602116_wp, -0.45173181_wp, -0.22799886_wp, -0.27816969_wp, -0.45230601_wp, &
   & -0.21669686_wp, -0.28031822_wp, -0.45288021_wp, -0.20539486_wp, -0.28246674_wp, -0.45345441_wp, &
   & -0.19409286_wp, -0.28461527_wp, -0.45402861_wp, -0.18279086_wp, -0.28676380_wp, -0.45460281_wp, &
   & -0.17148887_wp, -0.28891232_wp, -0.45517701_wp, -0.16018687_wp, -0.29106085_wp, -0.45575121_wp, &
   & -0.14888487_wp, -0.29320937_wp, -0.45632542_wp, -0.13758287_wp, -0.29535790_wp, -0.45689962_wp, &
   & -0.12628087_wp, -0.29750643_wp, -0.45747382_wp, -0.11497887_wp, -0.29965495_wp, -0.45804802_wp, &
   & -0.10367688_wp, -0.30180348_wp, -0.45862222_wp, -0.51413540_wp, -0.06435519_wp, -0.44273659_wp, &
   & -0.34859712_wp,  0.07959224_wp, -0.50110969_wp, -0.44821046_wp, -0.16269840_wp, -0.42026844_wp, &
   & -0.49513987_wp, -0.08699315_wp, -0.47219426_wp, -0.27378562_wp, -0.18370746_wp, -0.48313060_wp, &
   & -0.40293719_wp, -0.25066440_wp, -0.42516439_wp, -0.42858490_wp, -0.18782375_wp, -0.49166382_wp, &
   & -0.49533468_wp, -0.19144863_wp, -0.50644911_wp, -0.51587320_wp, -0.29145452_wp, -0.53800387_wp, &
   & -0.55257311_wp, -0.34898594_wp,  0.02606851_wp, -0.51728236_wp, -0.40709467_wp, -0.41799964_wp, &
   & -0.65355387_wp, -0.44221646_wp, -0.14897940_wp, -0.57437177_wp, -0.44267376_wp, -0.18531561_wp, &
   & -0.79336611_wp, -0.45416764_wp, -0.18618812_wp, -0.98227018_wp, -0.49955181_wp, -0.28755078_wp],&
   & shape(ceh_level))

   !> Dependence of orbital shell energy level on standard CN (shell-resolved) # MM, August 00, 2023
   real(wp), parameter :: ceh_kcn(max_shell, max_elem) = reshape([&
   & -0.01626957_wp,  0.00000000_wp,  0.00000000_wp, -0.14712486_wp,  0.00000000_wp,  0.00000000_wp, &
   &  0.00863841_wp,  0.00494496_wp,  0.00000000_wp, -0.08288884_wp,  0.00649581_wp,  0.00000000_wp, &
   &  0.01514747_wp, -0.00629044_wp,  0.00000000_wp, -0.00205353_wp, -0.04733032_wp,  0.00000000_wp, &
   & -0.05954040_wp, -0.08236997_wp,  0.00000000_wp, -0.07734635_wp, -0.16068330_wp,  0.00000000_wp, &
   & -0.01491119_wp, -0.27326542_wp,  0.00000000_wp, -0.05286409_wp, -0.40977373_wp,  0.00000000_wp, &
   &  0.10509451_wp, -0.01346266_wp,  0.00000000_wp, -0.02900627_wp, -0.01793791_wp,  0.00000000_wp, &
   &  0.00189897_wp, -0.00948821_wp, -0.02317406_wp, -0.01335877_wp, -0.01156999_wp, -0.00182405_wp, &
   &  0.03475586_wp, -0.03701144_wp, -0.02765386_wp, -0.01870669_wp, -0.07543686_wp, -0.01303772_wp, &
   &  0.04301284_wp, -0.11599174_wp, -0.03729682_wp, -0.20417822_wp, -0.07409862_wp, -0.02246543_wp, &
   &  0.34322188_wp,  0.01318320_wp,  0.00000000_wp, -0.01229019_wp,  0.03743403_wp,  0.01151337_wp, &
   & -0.12157369_wp, -0.06487047_wp,  0.12081987_wp,  0.05591202_wp, -0.03578924_wp,  0.04248427_wp, &
   & -0.04786277_wp, -0.02467660_wp,  0.00842692_wp, -0.01433692_wp,  0.03090773_wp, -0.01130523_wp, &
   & -0.05867434_wp, -0.02278037_wp, -0.01229249_wp,  0.00479160_wp, -0.03674844_wp, -0.02044685_wp, &
   & -0.09953767_wp,  0.02054226_wp, -0.03114463_wp,  0.11728275_wp, -0.05821580_wp, -0.03465944_wp, &
   &  0.03156397_wp,  0.00033800_wp, -0.05822776_wp, -0.04575552_wp, -0.01278541_wp, -0.11231628_wp, &
   & -0.01576923_wp,  0.02005206_wp, -0.05292968_wp, -0.03253403_wp,  0.00567502_wp,  0.00447648_wp, &
   & -0.01163885_wp, -0.02613921_wp,  0.00434336_wp,  0.00570942_wp, -0.06347916_wp, -0.00699032_wp, &
   &  0.05551107_wp, -0.09877832_wp, -0.01657677_wp, -0.24065971_wp, -0.05569917_wp,  0.02076511_wp, &
   &  0.37289705_wp,  0.00600915_wp,  0.00000000_wp,  0.00669819_wp,  0.07812375_wp, -0.00299273_wp, &
   &  0.08506897_wp, -0.05475700_wp,  0.03494791_wp, -0.01889109_wp, -0.04616713_wp,  0.02275678_wp, &
   & -0.00759356_wp, -0.05795806_wp,  0.01355074_wp, -0.04241473_wp,  0.02035114_wp, -0.00933872_wp, &
   & -0.04182419_wp, -0.04857672_wp,  0.00336497_wp, -0.09422350_wp, -0.01501399_wp, -0.01515242_wp, &
   & -0.10124578_wp,  0.00201696_wp, -0.04131672_wp,  0.03033805_wp, -0.08593120_wp, -0.02564891_wp, &
   &  0.08778117_wp, -0.01713825_wp, -0.07598911_wp, -0.02900415_wp, -0.02705660_wp,  0.15303569_wp, &
   & -0.04211022_wp,  0.02701719_wp, -0.08508553_wp, -0.00348280_wp, -0.01167965_wp,  0.01715785_wp, &
   &  0.03121942_wp, -0.02388064_wp, -0.00882082_wp, -0.00908386_wp, -0.05053057_wp, -0.00175911_wp, &
   &  0.28073776_wp, -0.07497037_wp, -0.06496945_wp, -0.04236284_wp, -0.05477997_wp,  0.02800537_wp, &
   &  0.16635470_wp,  0.07094546_wp,  0.00000000_wp,  0.03038949_wp,  0.06317391_wp,  0.01584630_wp, &
   &  0.07546793_wp, -0.04035467_wp,  0.04007654_wp, -0.10920538_wp,  0.01342832_wp,  0.03455530_wp, &
   & -0.11693057_wp,  0.01732656_wp,  0.03367679_wp, -0.12465575_wp,  0.02122480_wp,  0.03279829_wp, &
   & -0.13238094_wp,  0.02512304_wp,  0.03191978_wp, -0.14010612_wp,  0.02902128_wp,  0.03104127_wp, &
   & -0.14783131_wp,  0.03291952_wp,  0.03016276_wp, -0.15555649_wp,  0.03681775_wp,  0.02928425_wp, &
   & -0.16328167_wp,  0.04071599_wp,  0.02840574_wp, -0.17100686_wp,  0.04461423_wp,  0.02752723_wp, &
   & -0.17873204_wp,  0.04851247_wp,  0.02664872_wp, -0.18645723_wp,  0.05241071_wp,  0.02577022_wp, &
   & -0.19418241_wp,  0.05630895_wp,  0.02489171_wp, -0.20190759_wp,  0.06020718_wp,  0.02401320_wp, &
   & -0.20963278_wp,  0.06410542_wp,  0.02313469_wp,  0.02596189_wp, -0.06122644_wp,  0.03375121_wp, &
   & -0.00323387_wp, -0.09185345_wp,  0.03248178_wp, -0.01857132_wp,  0.00174848_wp, -0.01183479_wp, &
   & -0.02028725_wp, -0.05422920_wp,  0.00195873_wp, -0.07001977_wp, -0.02376318_wp, -0.01485681_wp, &
   & -0.10137711_wp,  0.00872310_wp, -0.04240400_wp, -0.03299086_wp, -0.01794242_wp, -0.05587321_wp, &
   &  0.04479228_wp, -0.03116888_wp, -0.09950389_wp, -0.07151097_wp,  0.01595180_wp, -0.09129392_wp, &
   & -0.05856389_wp,  0.02497592_wp, -0.05393627_wp, -0.09750939_wp,  0.01204451_wp, -0.01448104_wp, &
   &  0.01891168_wp, -0.01620237_wp, -0.00996082_wp, -0.04354041_wp, -0.04464475_wp,  0.02445958_wp, &
   &  0.13451910_wp, -0.04984156_wp, -0.03343742_wp,  0.08695188_wp, -0.06402913_wp,  0.02033733_wp],&
   & shape(ceh_kcn))

   !> Interaction type- and atom-wise resolved scal. fact. for overlap mat. elements # MM, August 02, 2023
   real(wp), parameter :: ceh_h0k(max_shell, max_elem) = reshape([&
   &  2.28423009_wp,  0.00000000_wp,  0.00000000_wp,  1.58277987_wp,  0.00000000_wp,  0.00000000_wp, &
   &  2.98891162_wp,  2.17097413_wp,  0.00000000_wp,  2.34504409_wp,  3.28893897_wp,  0.00000000_wp, &
   &  2.07564539_wp,  2.25048626_wp,  0.00000000_wp,  1.82191778_wp,  1.93197992_wp,  0.00000000_wp, &
   &  1.72051995_wp,  1.76469990_wp,  0.00000000_wp,  1.80961394_wp,  1.74570773_wp,  0.00000000_wp, &
   &  1.93384804_wp,  2.02401936_wp,  0.00000000_wp,  1.74025410_wp,  2.58304406_wp,  0.00000000_wp, &
   &  3.39014673_wp,  3.07156916_wp,  0.00000000_wp,  3.03055621_wp,  3.25403512_wp,  0.00000000_wp, &
   &  2.23746469_wp,  2.56462815_wp,  6.00000000_wp,  2.07861475_wp,  2.36094979_wp,  6.00000000_wp, &
   &  1.90209827_wp,  2.04051236_wp,  3.00000000_wp,  1.77730462_wp,  2.07830363_wp,  6.00000000_wp, &
   &  1.94529850_wp,  1.90586743_wp,  5.00000000_wp,  1.37890037_wp,  1.28294676_wp,  1.00000000_wp, &
   &  2.87448296_wp,  3.72226571_wp,  0.00000000_wp,  3.14084639_wp,  3.31703669_wp,  1.00000000_wp, &
   &  3.29568138_wp,  3.08465648_wp, 17.10675555_wp,  3.21154387_wp,  2.45464310_wp, 10.14857434_wp, &
   &  2.11408163_wp,  2.04530906_wp,  3.46167896_wp,  3.56167295_wp,  2.61312876_wp,  6.85137167_wp, &
   &  1.62878077_wp,  1.78203925_wp,  1.68990121_wp,  2.01115532_wp,  1.96976253_wp,  4.10611296_wp, &
   &  1.91973592_wp,  2.35569299_wp,  3.44672144_wp,  3.32694934_wp,  3.90587781_wp,  1.04934608_wp, &
   &  2.65648639_wp,  3.89784381_wp,  2.00000000_wp,  2.54632810_wp,  3.90115946_wp,  8.00000000_wp, &
   &  2.24804717_wp,  2.74423549_wp,  3.00000000_wp,  2.01811310_wp,  1.96329259_wp, 15.00000000_wp, &
   &  1.90311332_wp,  2.47030816_wp,  4.00000000_wp,  1.91826569_wp,  2.17391957_wp,  3.00000000_wp, &
   &  1.82793051_wp,  1.89810448_wp,  6.29065538_wp,  1.41919289_wp,  1.62689377_wp,  3.00000000_wp, &
   &  2.80939855_wp,  4.42883642_wp,  0.00000000_wp,  2.89996755_wp,  3.16966601_wp,  4.00000000_wp, &
   &  2.69813399_wp,  1.92787149_wp, 15.00000000_wp,  2.70415780_wp,  3.09701004_wp,  5.51720196_wp, &
   &  2.45919058_wp,  1.82695877_wp,  4.87327399_wp,  1.88861390_wp,  1.72703820_wp,  3.73801225_wp, &
   &  1.50572656_wp,  1.69657202_wp,  2.46044697_wp,  1.85064790_wp,  1.81317683_wp,  4.24377940_wp, &
   &  1.83220394_wp,  2.66380112_wp,  3.54924053_wp,  2.39574698_wp,  4.53044987_wp,  1.69378761_wp, &
   &  2.53183487_wp,  4.28938171_wp,  1.00000000_wp,  2.36735880_wp,  3.58658849_wp,  1.00000000_wp, &
   &  2.53514087_wp,  3.63918043_wp,  1.00000000_wp,  2.10857559_wp,  2.52249464_wp, 15.00000000_wp, &
   &  1.88894137_wp,  2.27206643_wp, 12.59823870_wp,  1.70796352_wp,  2.09146688_wp,  2.92888465_wp, &
   &  1.79341685_wp,  1.83721664_wp,  3.32908903_wp,  1.45123275_wp,  1.34649450_wp,  5.86632201_wp, &
   &  3.30475963_wp,  8.04774260_wp,  0.00000000_wp,  3.06227316_wp,  3.11860250_wp,  1.00000000_wp, &
   &  2.56722698_wp,  1.89071299_wp, 13.00000000_wp,  2.39955831_wp,  1.92593873_wp,  1.00000000_wp, &
   &  2.39122235_wp,  1.90607569_wp,  1.00000000_wp,  2.38288640_wp,  1.88621265_wp,  1.00000000_wp, &
   &  2.37455045_wp,  1.86634961_wp,  1.00000000_wp,  2.36621450_wp,  1.84648657_wp,  1.00000000_wp, &
   &  2.35787855_wp,  1.82662353_wp,  1.00000000_wp,  2.34954260_wp,  1.80676049_wp,  1.00000000_wp, &
   &  2.34120665_wp,  1.78689745_wp,  1.00000000_wp,  2.33287070_wp,  1.76703441_wp,  1.00000000_wp, &
   &  2.32453475_wp,  1.74717137_wp,  1.00000000_wp,  2.31619880_wp,  1.72730833_wp,  1.00000000_wp, &
   &  2.30786285_wp,  1.70744529_wp,  1.00000000_wp,  2.29952689_wp,  1.68758225_wp,  1.00000000_wp, &
   &  2.29119094_wp,  1.66771921_wp,  1.00000000_wp,  2.79331545_wp,  2.50791914_wp,  6.95602975_wp, &
   &  2.53172597_wp,  1.90212322_wp,  2.95767860_wp,  1.76773228_wp,  1.83393645_wp,  3.70487961_wp, &
   &  1.51411461_wp,  1.73091670_wp,  2.67692862_wp,  1.97219223_wp,  2.09071992_wp,  1.49630847_wp, &
   &  1.89488596_wp,  2.71429561_wp,  2.26016227_wp,  2.19525073_wp,  3.60152743_wp,  0.15145647_wp, &
   &  2.38603366_wp,  5.08498430_wp,  1.00000000_wp,  2.29577007_wp,  4.98265552_wp,  1.00000000_wp, &
   &  2.55949715_wp,  3.79876876_wp,  1.00000000_wp,  2.41361463_wp,  2.59977715_wp,  2.00000000_wp, &
   &  1.84332653_wp,  2.51356296_wp, 15.00000000_wp,  1.73599763_wp,  2.15552900_wp,  2.50000000_wp, &
   &  1.71327975_wp,  1.72568101_wp,  6.36853611_wp,  1.51399613_wp,  1.41998519_wp, 15.00000000_wp],&
   & shape(ceh_h0k))

   !> Dependence of orbital shell energy level on EN-weighted CN (atom-resolved) # MM, August 00, 2023
   real(wp), parameter :: ceh_kcnen(max_elem) = [&
   & -0.12768184_wp, -0.30481937_wp,  0.00302231_wp,  0.07200407_wp, -0.04294759_wp, -0.07955338_wp, &
   & -0.09493707_wp, -0.07198264_wp,  0.03080274_wp, -0.01999000_wp, -0.00232660_wp,  0.04746226_wp, &
   &  0.01262132_wp, -0.00849344_wp, -0.00867003_wp, -0.00730045_wp,  0.01104749_wp, -0.00040887_wp, &
   & -0.08045470_wp,  0.01026750_wp,  0.04205488_wp, -0.01175355_wp, -0.02151850_wp, -0.04419584_wp, &
   & -0.02732939_wp, -0.04970992_wp, -0.04216317_wp, -0.03601579_wp, -0.02613562_wp,  0.05601951_wp, &
   &  0.02220564_wp,  0.00675828_wp, -0.01083386_wp, -0.01866043_wp,  0.00431215_wp, -0.01972074_wp, &
   & -0.04276521_wp, -0.00344024_wp,  0.00887900_wp,  0.02360741_wp, -0.02038550_wp, -0.02639773_wp, &
   & -0.02954938_wp, -0.03654088_wp, -0.02867754_wp, -0.00282480_wp, -0.00734717_wp,  0.04537368_wp, &
   &  0.04130496_wp, -0.02390599_wp, -0.01759813_wp, -0.00914366_wp,  0.01045523_wp, -0.03440882_wp, &
   & -0.07875308_wp, -0.02476772_wp,  0.00043117_wp,  0.04086382_wp,  0.03949769_wp,  0.03813156_wp, &
   &  0.03676543_wp,  0.03539930_wp,  0.03403318_wp,  0.03266705_wp,  0.03130092_wp,  0.02993479_wp, &
   &  0.02856866_wp,  0.02720253_wp,  0.02583640_wp,  0.02447027_wp,  0.02310414_wp,  0.00251351_wp, &
   & -0.01970384_wp, -0.02447100_wp, -0.02845688_wp, -0.02733732_wp, -0.02856101_wp,  0.01266183_wp, &
   &  0.02207229_wp,  0.06732818_wp,  0.03457416_wp, -0.00347199_wp, -0.01934747_wp, -0.00166422_wp, &
   & -0.02282661_wp, -0.01354102_wp]

   !> Angular momentum-specific scaling factors for H0
   real(wp), parameter   :: kll(1:3) = [0.6366_wp, 0.9584_wp, 1.2320_wp]
   !> Conversion constant
   real(wp), parameter   :: kt = 3.166808578545117e-06_wp

   character(len=*), parameter :: real_format = "(es18.6)"
   character(len=*), parameter :: &
   &  label_charges = "CEH atomic charges", &
   &  label_dipole = "CEH molecular dipole moment / a.u."

contains

   !> Run the CEH calculation
   subroutine ceh_guess(ctx, calc, mol, error, wfn,verbosity)
      !> Calculation context
      type(context_type), intent(inout) :: ctx
      !> CEH calculator
      type(ceh_calculator), intent(inout) :: calc
      !> Molecular structure data
      type(structure_type), intent(in)  :: mol
      !> Error container
      type(error_type), allocatable, intent(out) :: error
      !> Wavefunction data
      type(wavefunction_type), intent(inout) :: wfn
      !> Verbosity level of output
      integer, intent(in), optional :: verbosity
      !> Molecular dipole moment
      real(wp) :: dipole(3) = 0.0_wp
      !> Integral container
      type(integral_type) :: ints
      !> Electronic solver
      class(solver_type), allocatable :: solver
      !> Potential type
      type(potential_type) :: pot
      !> Restart data for interaction containers
      type(container_cache) :: icache
      !> Timer
      type(timer_type) :: timer
      real(wp) :: ttime

      logical :: grad = .false.

      real(wp) :: elec_entropy
      real(wp) :: nel = 0.0_wp
      real(wp), allocatable :: tmp(:)

      integer :: i, prlevel

      call timer%push("wall time CEH")

      if (present(verbosity)) then
         prlevel = verbosity
      else
         prlevel = ctx%verbosity
      end if

      if (prlevel > 2) then
         call header(ctx)
      elseif (prlevel > 1) then
         call ctx%message("CEH guess")
      endif
      !> Gradient logical as future starting point (not implemented yet)
      !> Entry point could either be (i) modified wavefunction type (including derivatives),
      !> (iii) additional wavefunction derivative type (see old commits) or (ii) optional 
      !> dqdR variable in this routine
      grad = .false.

      !> Reference occupation for number of electrons and formal charges
      call get_reference_occ(mol, calc%bas, calc%hamiltonian%refocc)
      !> Define occupation
      call get_occupation(mol, calc%bas, calc%hamiltonian, wfn%nocc, wfn%n0at, wfn%n0sh)
      nel = sum(wfn%n0at) - mol%charge
      if (mod(mol%uhf, 2) == mod(nint(nel), 2)) then
         wfn%nuhf = mol%uhf
      else
         wfn%nuhf = mod(nint(nel), 2)
      end if
      call get_alpha_beta_occupation(wfn%nocc, wfn%nuhf, wfn%nel(1), wfn%nel(2))

      !> Initialize integrals
      call new_integral(ints, calc%bas%nao)
      ints%quadrupole = 0.0_wp
      !> Get Hamiltonian and integrals
      call get_hamiltonian(calc, mol, ints%overlap, ints%dipole, ints%hamiltonian)

      !> Get initial potential
      call new_potential(pot, mol, calc%bas, wfn%nspin)
      !> Set potential to zero
      call pot%reset
      if (allocated(calc%interactions)) then
         call calc%interactions%update(mol, icache)
         call calc%interactions%get_potential(mol, icache, wfn, pot)
      endif

      !> Add effective Hamiltonian to wavefunction
      call add_pot_to_h1(calc%bas, ints, pot, wfn%coeff)

      !> Solve the effective Hamiltonian
      call ctx%new_solver(solver, calc%bas%nao)

      !> Get the density matrix
      call get_density(wfn, solver, ints, elec_entropy, error)
      if (allocated(error)) then
         call ctx%set_error(error)
      end if

      !> Get charges and dipole moment from density and integrals
      call get_mulliken_shell_charges(calc%bas, ints%overlap, wfn%density, wfn%n0sh, &
      & wfn%qsh)
      call get_qat_from_qsh(calc%bas, wfn%qsh, wfn%qat)
      call get_mulliken_atomic_multipoles(calc%bas, ints%dipole, wfn%density, &
      & wfn%dpat)
      allocate(tmp(3), source = 0.0_wp)
      call gemv(mol%xyz, wfn%qat(:, 1), tmp)
      dipole(:) = tmp + sum(wfn%dpat(:, :, 1), 2)

      call timer%pop
      ttime = timer%get("wall time CEH")

      !> Printout of results
      if (prlevel > 2) then
         call ctx%message(label_charges)
         call ctx%message("Atom index      Charge / a.u.")
         do i = 1, mol%nat
            call ctx%message(format_string(i, "(i7)") // &
            & "   " // format_string(wfn%qat(i,1), real_format))
         end do
         call ctx%message(repeat("-", 60))
         call ctx%message(label_dipole)
         call ctx%message("     x           y           z")
         call ctx%message(format_string(dipole(1), "(f12.5)") // &
         & format_string(dipole(2), "(f12.5)") // &
         & format_string(dipole(3), "(f12.5)"))
         call ctx%message(repeat("-", 60))
      endif
      if (prlevel > 0) then
         call ctx%message(" - CEH single point"//repeat(" ", 4)//format_time(ttime))
         call ctx%message("")
      endif

   end subroutine ceh_guess

   subroutine new_ceh_calculator(calc,mol)
      !> Instance of the CEH evaluator
      type(ceh_calculator), intent(out) :: calc
      type(structure_type), intent(in)  :: mol

      call add_ceh_basis(calc, mol)
      call add_ncoord(calc, mol)


   end subroutine new_ceh_calculator

   subroutine add_ncoord(calc, mol)
      !> Instance of the xTB evaluator
      type(ceh_calculator), intent(inout) :: calc
      !> Molecular structure data
      type(structure_type), intent(in) :: mol

      call new_ncoord(calc%ncoordstd, mol, cn_type="ceh_std")
      call new_ncoord(calc%ncoorden, mol, cn_type="ceh_en")
   end subroutine add_ncoord

   subroutine add_ceh_basis(calc, mol)
      !> Instance of the CEH evaluator
      type(ceh_calculator), intent(inout) :: calc
      !> Molecular structure data
      type(structure_type), intent(in) :: mol

      integer :: isp, izp, ish, stat, ng, il
      integer, allocatable :: nsh_id(:)
      integer :: ang_idx(0:2), ortho(max_shell)
      type(cgto_type), allocatable :: cgto(:, :)

      nsh_id = nshell(mol%num)
      allocate(cgto(maxval(nsh_id), mol%nid))
      do isp = 1, mol%nid
         ang_idx = 0
         ortho = 0
         izp = mol%num(isp)
         do ish = 1, nsh_id(isp)
            il = ang_shell(ish, izp)
            ng = number_of_primitives(ish, izp)
            if (ang_idx(il) > 0) then
               ortho(ish) = ang_idx(il)
            else
               ang_idx(il) = ish
            end if
            call slater_to_gauss(ng, principal_quantum_number(ish, izp), il, &
            & slater_exponent(ish, izp), cgto(ish, isp), .true., stat)
         end do

         do ish = 1, nsh_id(isp)
            if (ortho(ish) > 0) then
               call orthogonalize(cgto(ortho(ish), isp), cgto(ish, isp))
            end if
         end do
      end do

      call new_basis(calc%bas, mol, nsh_id, cgto, 1.0_wp)

   end subroutine add_ceh_basis

   subroutine get_reference_occ(mol, bas, refocc)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set information
      type(basis_type), intent(in) :: bas
      !> Reference occupation numbers
      real(wp), allocatable, intent(out) :: refocc(:, :)
      logical, allocatable  :: valence(:,:)

      integer :: isp, izp, ish, il, mshell
      integer :: ang_idx(0:2)


      allocate(valence(3, mol%nid))
      do isp = 1, mol%nid
         ang_idx = 0
         izp = mol%num(isp)
         do ish = 1, nshell(izp)
            il = ang_shell(ish, izp)
            valence(ish, isp) = ang_idx(il) == 0
            if (valence(ish, isp)) ang_idx(il) = ish
         end do
      end do

      mshell = maxval(bas%nsh_id)
      allocate(refocc(mshell, mol%nid))
      refocc(:, :) = 0.0_wp
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, bas%nsh_id(isp)
            if (valence(ish, isp)) then
               refocc(ish, isp) = reference_occ(bas%cgto(ish, isp)%ang+1, izp)
            else
               refocc(ish, isp) = 0.0_wp
            end if
         end do
      end do
   end subroutine get_reference_occ

   subroutine get_occupation(mol, bas, hamiltonian, nocc, n0at, n0sh)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set information
      type(basis_type), intent(in) :: bas
      !> Hamiltonian interaction data
      type(ceh_hamiltonian), intent(in) :: hamiltonian
      !> Occupation number
      real(wp), intent(out) :: nocc
      !> Reference occupation for each atom
      real(wp), intent(out) :: n0at(:)
      !> Reference occupation for each shell
      real(wp), intent(out) :: n0sh(:)

      integer :: iat, ish, izp, ii

      nocc = -mol%charge
      n0at(:) = 0.0_wp
      n0sh(:) = 0.0_wp
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ii = bas%ish_at(iat)
         do ish = 1, bas%nsh_id(izp)
            nocc = nocc + hamiltonian%refocc(ish, izp)
            n0at(iat) = n0at(iat) + hamiltonian%refocc(ish, izp)
            n0sh(ii+ish) = n0sh(ii+ish) + hamiltonian%refocc(ish, izp)
         end do
      end do

   end subroutine get_occupation

   subroutine get_hamiltonian(calc, mol, overlap, dipole, hamiltonian)
      !> CEH Hamiltonian type
      type(ceh_hamiltonian)               :: self
      !> CEH calculator
      type(ceh_calculator), intent(inout) :: calc
      !> Molecular structure type
      type(structure_type), intent(in)    :: mol
      !> Overlap integral matrix
      real(wp), allocatable, intent(out)  :: overlap(:,:)
      !> Dipole moment integral matrix
      real(wp), allocatable, intent(out)  :: dipole(:, :, :)
      !> Full Hamiltonian matrix
      real(wp), allocatable, intent(out)  :: hamiltonian(:,:)
      !> Scaling factors for the diatomic frame for the three differnt bonding motifs
      !> (sigma, pi, delta)
      real(wp) :: ksig, kpi, kdel

      real(wp), allocatable   :: cn(:), cn_en(:), & 
      & overlap_diat(:,:)
      real(wp), allocatable :: overlap_tmp(:), dipole_tmp(:, :), &
      & overlap_scaled_tmp(:)
      real(wp) :: cutoff, cutoff2, felem, r2, vec(3), vec_diat_trafo(3)
      real(wp) :: dipole_tmp_tf(3)

      integer                 :: ii, offset_iat, offset_jat, jzp, izp, nao
      integer                 :: iat, ish, jsh, k, l, iao, jao, jat, kl

      !> allocate CEH diagonal elements
      allocate(self%hlevel(calc%bas%nsh), source=0.0_wp)

      !> allocate full CEH Hamiltonian
      allocate(hamiltonian(calc%bas%nao, calc%bas%nao), source=0.0_wp)

      !> calculate coordination number (CN) and CN-weighted energy
      if (allocated(calc%ncoordstd)) then
         allocate(cn(mol%nat))
         call calc%ncoordstd%get_cn(mol, cn)
      end if
      if (allocated(calc%ncoorden)) then
         allocate(cn_en(mol%nat))
         call calc%ncoorden%get_cn(mol, cn_en)
      end if

      !> define diagonal elements of CEH Hamiltonian
      !> shell-resolved, not AO resolved
      do iat = 1, mol%nat
         ii = calc%bas%ish_at(iat)
         do ish = 1, calc%bas%nsh_at(iat)
            self%hlevel(ii+ish) = ceh_level(ish, mol%num(mol%id(iat))) + &
            &             ceh_kcn(ish, mol%num(mol%id(iat))) * cn(iat) + &
            &                ceh_kcnen(mol%num(mol%id(iat))) * cn_en(iat)
         end do
      end do

      cutoff = get_cutoff(calc%bas)
      cutoff2 = cutoff**2

      !> Allocate matrices for overlap, scaled overlap and dipole moment integrals
      allocate(overlap(calc%bas%nao, calc%bas%nao), &
      & overlap_diat(calc%bas%nao, calc%bas%nao), &
      & dipole(3, calc%bas%nao, calc%bas%nao), source = 0.0_wp)
      !> Allocate temporary matrices for overlap, scaled overlap and dipole moment integrals
      allocate(overlap_tmp(msao(calc%bas%maxl)**2), &
      & overlap_scaled_tmp(msao(calc%bas%maxl)**2), & 
      & dipole_tmp(3, msao(calc%bas%maxl)**2))

      !> define effective CEH Hamiltonian
      do iat = 1, mol%nat
      izp = mol%id(iat)
         offset_iat = calc%bas%ish_at(iat)
         do ish = 1, calc%bas%nsh_at(iat)

            !> loop over all AOs in atoms before current atom
            !> -- AO loop over ish-AOs (-> iao) is done in the individual loops
            !> -- to simplifiy overlap integral calculation
            do jat = 1,iat-1
               jzp = mol%id(jat)
               offset_jat = calc%bas%ish_at(jat)
               !> Calculate relative integral aufpoint vectors
               vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
               r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
               if (r2 > cutoff2) cycle
               call relvec(vec, sqrt(r2), vec_diat_trafo)

               do jsh = 1, calc%bas%nsh_at(jat)
                  felem = ceh_h0_entry_od(mol%num(mol%id(iat)), mol%num(mol%id(jat)), ish, jsh, &
                  & self%hlevel(offset_iat + ish), self%hlevel(offset_jat + jsh))

                  !> Determine factors for diatomic overlap scaling from atom parameters
                  ksig = 2.0_wp / (1.0_wp / ceh_h0k(1,mol%num(mol%id(iat))) &
                  & + 1.0_wp / ceh_h0k(1,mol%num(mol%id(jat))) )
                  kpi = 2.0_wp / (1.0_wp / ceh_h0k(2,mol%num(mol%id(iat))) &
                  & + 1.0_wp / ceh_h0k(2,mol%num(mol%id(jat))) )
                  kdel = 2.0_wp / (1.0_wp / ceh_h0k(3,mol%num(mol%id(iat))) &
                  & + 1.0_wp / ceh_h0k(3,mol%num(mol%id(jat))) )
                  
                  !> Integral call for different atom and different shell (-> with diatomic scaling)
                  overlap_tmp = 0.0_wp
                  overlap_scaled_tmp = 0.0_wp
                  dipole_tmp = 0.0_wp
                  call dipole_cgto_diat_scal(calc%bas%cgto(jsh,jzp), calc%bas%cgto(ish,izp), &
                  & r2, vec, calc%bas%intcut, vec_diat_trafo, ksig, kpi, kdel, & 
                  & overlap_tmp, overlap_scaled_tmp, dipole_tmp)

                  nao = msao(calc%bas%cgto(jsh, jzp)%ang)
                  do iao = 1, calc%bas%nao_sh(ish + offset_iat)
                     !> AO iterator of i 
                     k = calc%bas%iao_sh(ish + offset_iat) + iao
                     do jao = 1, calc%bas%nao_sh(jsh + offset_jat)
                        l = calc%bas%iao_sh(jsh + offset_jat) + jao
                        kl = jao + nao*(iao-1)

                        overlap_diat(k, l) = overlap_scaled_tmp(kl)
                        overlap_diat(l, k) = overlap_scaled_tmp(kl)
                        overlap(k, l) = overlap_tmp(kl)
                        overlap(l, k) = overlap_tmp(kl)

                        hamiltonian(k, l) = overlap_diat(k, l) * felem
                        hamiltonian(l, k) = hamiltonian(k, l)

                        !> Shift dipole operator from Ket function (center i) 
                        !> to Bra function (center j) to save the redundant calculation 
                        call shift_operator(vec, overlap_tmp(kl), &
                           & dipole_tmp(:, kl), dipole_tmp_tf)
                        !> Order of l and k is not relevant for the result 
                        !> of the dipole moment but in this ordering it corresponds exactly to 
                        !> the result by the 'call get_dipole_integrals'
                        dipole(:, l, k) = dipole_tmp(:, kl)
                        dipole(:, k, l) = dipole_tmp_tf(:)
                     end do
                  enddo
               end do
            end do

            !> reset vectors to zero because we are only considering 1-center terms from now on
            vec(:) = 0.0_wp 
            r2 = 0.0_wp 

            !> loop over all AOs in shells before current shell (same atom)
            do jsh = 1, ish - 1
               felem = ceh_h0_entry_od(mol%num(mol%id(iat)), mol%num(mol%id(iat)), ish, jsh, &
               & self%hlevel(offset_iat + ish), self%hlevel(offset_iat + jsh))

               !> Integral call for same atom and different shell
               overlap_tmp = 0.0_wp
               dipole_tmp = 0.0_wp
               call dipole_cgto(calc%bas%cgto(jsh,izp), calc%bas%cgto(ish,izp), &
               & r2, vec, calc%bas%intcut, overlap_tmp,  dipole_tmp)

               nao = msao(calc%bas%cgto(jsh, izp)%ang)
               do iao = 1, calc%bas%nao_sh(ish + offset_iat)
                  !> AO iterator of i 
                  k = calc%bas%iao_sh(ish + offset_iat) + iao
                  do jao = 1, calc%bas%nao_sh(jsh + offset_iat)
                     l = calc%bas%iao_sh(jsh + offset_iat) + jao
                     kl = jao + nao*(iao-1)

                     overlap_diat(k, l) = overlap_tmp(kl)
                     overlap_diat(l, k) = overlap_tmp(kl)
                     overlap(k, l) = overlap_tmp(kl)
                     overlap(l, k) = overlap_tmp(kl)

                     hamiltonian(k, l) = overlap_diat(k, l) * felem
                     hamiltonian(l, k) = hamiltonian(k, l)

                     dipole(:, k, l) = dipole_tmp(:, kl)
                     dipole(:, l, k) = dipole_tmp(:, kl)
                  end do
               end do
            end do
            
            !> Integral call for same atom and same shell
            overlap_tmp = 0.0_wp
            dipole_tmp = 0.0_wp
            call dipole_cgto(calc%bas%cgto(ish,izp), calc%bas%cgto(ish,izp), &
            & r2, vec, calc%bas%intcut, overlap_tmp, dipole_tmp)
            
            !> loop over all AOs before current AO (same atom and shell)
            felem = ceh_h0_entry_od(mol%num(mol%id(iat)), mol%num(mol%id(iat)), ish, ish, &
            & self%hlevel(offset_iat + ish), self%hlevel(offset_iat + ish))
            nao = msao(calc%bas%cgto(ish, izp)%ang)
            do iao = 1, calc%bas%nao_sh(ish + offset_iat)
               !> AO iterator of i 
               k = calc%bas%iao_sh(ish + offset_iat) + iao
               do jao = 1, iao - 1
                  l = calc%bas%iao_sh(ish + offset_iat) + jao
                  kl = jao + nao*(iao-1)

                  overlap_diat(k, l) = overlap_tmp(kl)
                  overlap_diat(l, k) = overlap_tmp(kl)
                  overlap(k, l) = overlap_tmp(kl)
                  overlap(l, k) = overlap_tmp(kl)

                  hamiltonian(k, l) = overlap_diat(k, l) * felem
                  hamiltonian(l, k) = hamiltonian(k, l)

                  dipole(:, k, l) = dipole_tmp(:, jao + nao*(iao-1))
                  dipole(:, l, k) = dipole_tmp(:, jao + nao*(iao-1))
               end do
               !> diagonal term (AO(i) == AO(j))
               kl = iao + nao*(iao-1)

               overlap_diat(k, k) = overlap_tmp(kl)
               overlap(k, k) = overlap_tmp(kl)
               hamiltonian(k, k) = self%hlevel(offset_iat + ish)

               dipole(:, k, k) = dipole_tmp(:, kl)
            end do
         enddo
      enddo
      if (k /= calc%bas%nao) then
         error stop "ERROR: k /= calc%bas%nao"
      end if

      calc%hamiltonian = self

   end subroutine get_hamiltonian

   !> Shift dipole operator from Ket function (center i) to Bra function (center j),
   !> the dipole operator on the Bra function can be assembled from the lower moments
   !> on the Ket function and the displacement vector using horizontal shift rules.
   pure subroutine shift_operator(vec, s, di, dj)
      !> Displacement vector of center i and j
      real(wp),intent(in) :: vec(:)
      !> Overlap integral between basis functions
      real(wp),intent(in) :: s
      !> Dipole integral with operator on Ket function (center i)
      real(wp),intent(in) :: di(:)
      !> Dipole integral with operator on Bra function (center j)
      real(wp),intent(out) :: dj(:)

      ! Create dipole operator on Bra function from Ket function and shift contribution
      ! due to monopol displacement
      dj(1) = di(1) + vec(1)*s
      dj(2) = di(2) + vec(2)*s
      dj(3) = di(3) + vec(3)*s

   end subroutine shift_operator

   function ceh_h0_entry_od(ati,atj,ish,jsh, hleveli, hlevelj) result(level)
      integer, intent(in) :: ish, ati, jsh, atj
      real(wp), intent(in) :: hleveli, hlevelj
      real(wp) :: level

      integer :: ang_i, ang_j

      ang_i = ang_shell(ish, ati)
      ang_j = ang_shell(jsh, atj)
      level = 0.25_wp * (hleveli + hlevelj) * &
      & ( kll(ang_i + 1) + kll(ang_j + 1) )

   end function ceh_h0_entry_od

   subroutine header(ctx)
      !> Calculation context
      type(context_type), intent(inout) :: ctx
      call ctx%message(repeat("-", 60))
      call ctx%message('       C H A R G E    E X T E N D E D    H U C K E L (CEH) ')
      call ctx%message('                       SG, MM, AH, TF, May 2023            ')
      call ctx%message(repeat("-", 60))
   end subroutine header

end module tblite_ceh_ceh
