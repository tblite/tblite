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

module test_ceh
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
   & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_context_type, only : context_type
   use tblite_lapack_solver, only : lapack_solver, lapack_algorithm
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_cutoff, only : get_lattice_points
   use tblite_basis_type
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_xtb_h0, only : tb_hamiltonian, new_hamiltonian
   use tblite_data_covrad, only : get_covalent_rad
   use tblite_data_paulingen, only : get_pauling_en
   use tblite_ncoord_erf
   use tblite_ncoord_erf_en
   use tblite_ncoord_type, only : get_coordination_number

   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_ceh_singlepoint, only : ceh_singlepoint
   use tblite_ceh_ceh, only : ceh_h0spec, new_ceh_calculator
   use tblite_ceh_h0, only : get_scaled_selfenergy, get_hamiltonian

   use tblite_blas, only: gemv

   use tblite_container, only : container_type, container_cache
   use tblite_external_field, only : electric_field
   implicit none
   private

   public :: collect_ceh

   real(wp), parameter :: kt = 4000.0_wp * 3.166808578545117e-06_wp
   real(wp), parameter :: thr2 = 1.0e2_wp*sqrt(epsilon(1.0_wp))

contains

   !> Collect all exported unit tests
   subroutine collect_ceh(testsuite)

      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("scaled-selfenergy-H2", test_scaled_selfenergy_h2), &
         new_unittest("scaled-selfenergy-LiH", test_scaled_selfenergy_lih), &
         new_unittest("scaled-selfenergy-S2", test_scaled_selfenergy_s2), &
         new_unittest("scaled-selfenergy-SiH4", test_scaled_selfenergy_sih4), &
         new_unittest("scaled-selfenergy-AcCl6", test_scaled_selfenergy_accl6), &
         new_unittest("hamiltonian-H2", test_hamiltonian_h2), &
         new_unittest("hamiltonian-LiH", test_hamiltonian_lih), &
         new_unittest("hamiltonian-S2", test_hamiltonian_s2), &
         new_unittest("hamiltonian-SiH4", test_hamiltonian_sih4), &
         new_unittest("overlap_diat-H2", test_overlap_diat_h2), &
         new_unittest("overlap_diat-LiH", test_overlap_diat_lih), &
         new_unittest("overlap_diat-S2", test_overlap_diat_s2), &
         new_unittest("overlap_diat-SiH4", test_overlap_diat_sih4), &
         new_unittest("q-mol-h2", test_q_h2), &
         new_unittest("q-mol-lih", test_q_lih), &
         new_unittest("q-mol-sih4", test_q_sih4), &
         new_unittest("q-mol-cecl3", test_q_cecl3), &
         new_unittest("q-mol-accl6", test_q_accl6), &
         new_unittest("q-mol-panp", test_q_panp), &
         new_unittest("q-mol-mb01", test_q_mb01), &
         new_unittest("q-mol-mb02", test_q_mb02), &
         new_unittest("q-mol-mb03", test_q_mb03), &
         new_unittest("q-mol-mb04", test_q_mb04), &
         new_unittest("q-chrgd-efield-mol", test_q_ef_chrg_mb01), &
         new_unittest("d-mol", test_d_mb01), &
         new_unittest("d-field-mol", test_d_field_mb04), &
         new_unittest("d-field-change-mol", test_d_hcn) &
         ]

   end subroutine collect_ceh

   !> Testing on the CEH basis
   subroutine make_basis(bas, mol, ng)
      type(basis_type), intent(out) :: bas
      type(structure_type), intent(in) :: mol
      integer, intent(in) :: ng

      integer, parameter :: nsh(103) = [&
      & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, & ! 1-20
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, & ! 21-40
      & 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, & ! 41-60
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, & ! 61-80
      & 3, 3, 3, 3, 3, 3, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, & ! 81-100
      & 4, 4, 4]
      integer, parameter :: lsh(4, 103) = reshape([&
      & 0, 0, 0, 0,  0, 0, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0, & ! 1-6
      & 0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 2, 0, & ! 7-12
      & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 13-18
      & 0, 1, 0, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 19-24
      & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 0, 0, & ! 25-30
      & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 31-36
      & 0, 1, 0, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 37-42
      & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 0, 0, & ! 43-48
      & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 49-54
      & 0, 1, 0, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 55-60
      & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 61-66
      & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 67-72
      & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 73-78
      & 0, 1, 2, 0,  0, 1, 0, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 79-84
      & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 0, 0,  0, 1, 2, 0,  0, 1, 2, 3,  0, 1, 2, 3, & ! 85-90
      & 0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3, & ! 91-96
      & 0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3, & ! 97-102
      & 0, 1, 2, 3], shape(lsh))

      integer, parameter :: pqn(4, 103) = reshape([&
      & 1, 0, 0, 0,  1, 0, 0, 0,  2, 2, 0, 0,  2, 2, 0, 0,  2, 2, 0, 0,  2, 2, 0, 0, & ! 1-6
      & 2, 2, 0, 0,  2, 2, 0, 0,  2, 2, 0, 0,  2, 2, 0, 0,  3, 3, 0, 0,  3, 3, 3, 0, & ! 7-12
      & 3, 3, 3, 0,  3, 3, 3, 0,  3, 3, 3, 0,  3, 3, 3, 0,  3, 3, 3, 0,  3, 3, 3, 0, & ! 13-18
      & 4, 4, 0, 0,  4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 3, 0, & ! 19-24
      & 4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 0, 0, & ! 25-30
      & 4, 4, 4, 0,  4, 4, 4, 0,  4, 4, 4, 0,  4, 4, 4, 0,  4, 4, 4, 0,  4, 4, 4, 0, & ! 31-36
      & 5, 5, 0, 0,  5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 4, 0, & ! 37-42
      & 5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 0, 0, & ! 43-48
      & 5, 5, 5, 0,  5, 5, 5, 0,  5, 5, 5, 0,  5, 5, 5, 0,  5, 5, 5, 0,  5, 5, 5, 0, & ! 49-54
      & 6, 6, 0, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0, & ! 55-60
      & 6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0, & ! 61-66
      & 6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0, & ! 67-72
      & 6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0, & ! 73-78
      & 6, 6, 5, 0,  6, 6, 0, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0, & ! 79-84
      & 6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 0, 0,  6, 6, 5, 0,  6, 6, 5, 5,  6, 6, 5, 5, & ! 85-90
      & 6, 6, 5, 5,  6, 6, 5, 5,  6, 6, 5, 5,  6, 6, 5, 5,  6, 6, 5, 5,  6, 6, 5, 5, & ! 91-96
      & 6, 6, 5, 5,  6, 6, 5, 5,  6, 6, 5, 5,  6, 6, 5, 5,  6, 6, 5, 5,  6, 6, 5, 5, & ! 97-102
      & 6, 6, 5, 5], shape(pqn))

      real(wp), parameter :: zeta(4, 103) = reshape([&
      & 1.23363166_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, & ! 1
      & 2.27004605_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, & ! 2
      & 0.86185456_wp, 1.42017184_wp, 0.00000000_wp, 0.00000000_wp, & ! 3
      & 1.76817995_wp, 1.44095844_wp, 0.00000000_wp, 0.00000000_wp, & ! 4
      & 2.06339837_wp, 1.52051807_wp, 0.00000000_wp, 0.00000000_wp, & ! 5
      & 2.56058582_wp, 1.86484737_wp, 0.00000000_wp, 0.00000000_wp, & ! 6
      & 2.71233631_wp, 2.19848968_wp, 0.00000000_wp, 0.00000000_wp, & ! 7
      & 3.21585650_wp, 2.41309737_wp, 0.00000000_wp, 0.00000000_wp, & ! 8
      & 3.82146807_wp, 2.63063636_wp, 0.00000000_wp, 0.00000000_wp, & ! 9
      & 4.62721228_wp, 2.53599954_wp, 0.00000000_wp, 0.00000000_wp, & ! 10
      & 0.93221172_wp, 1.55333839_wp, 0.00000000_wp, 0.00000000_wp, & ! 11
      & 1.77220557_wp, 1.59942632_wp, 2.98596647_wp, 0.00000000_wp, & ! 12
      & 2.26040231_wp, 1.78718151_wp, 2.00990188_wp, 0.00000000_wp, & ! 13
      & 1.85259089_wp, 1.81733349_wp, 1.65269988_wp, 0.00000000_wp, & ! 14
      & 2.65701241_wp, 2.03189759_wp, 2.03883661_wp, 0.00000000_wp, & ! 15
      & 2.60609998_wp, 2.16530440_wp, 2.41888232_wp, 0.00000000_wp, & ! 16
      & 2.78818934_wp, 2.24732894_wp, 1.99081182_wp, 0.00000000_wp, & ! 17
      & 2.55424399_wp, 2.20946190_wp, 1.93619550_wp, 0.00000000_wp, & ! 18
      & 1.73713827_wp, 1.33788617_wp, 0.00000000_wp, 0.00000000_wp, & ! 19
      & 2.47982574_wp, 1.07250770_wp, 2.11920764_wp, 0.00000000_wp, & ! 20
      & 2.22449249_wp, 1.55418319_wp, 2.00953578_wp, 0.00000000_wp, & ! 21
      & 2.58879616_wp, 0.99441077_wp, 1.88561781_wp, 0.00000000_wp, & ! 22
      & 3.04370654_wp, 4.03007600_wp, 1.66329169_wp, 0.00000000_wp, & ! 23
      & 2.25012727_wp, 2.70681556_wp, 1.67501904_wp, 0.00000000_wp, & ! 24
      & 2.20605319_wp, 2.82019792_wp, 1.86102254_wp, 0.00000000_wp, & ! 25
      & 1.57297015_wp, 1.98621494_wp, 2.83790684_wp, 0.00000000_wp, & ! 26
      & 1.80826602_wp, 1.73675835_wp, 2.79767448_wp, 0.00000000_wp, & ! 27
      & 2.00758945_wp, 2.25075692_wp, 2.98291663_wp, 0.00000000_wp, & ! 28
      & 2.18159986_wp, 2.38459096_wp, 3.09502522_wp, 0.00000000_wp, & ! 29
      & 2.26376756_wp, 2.20362977_wp, 0.00000000_wp, 0.00000000_wp, & ! 30
      & 2.63822153_wp, 2.06752328_wp, 2.11361643_wp, 0.00000000_wp, & ! 31
      & 2.52891955_wp, 2.19441794_wp, 1.77661998_wp, 0.00000000_wp, & ! 32
      & 3.55667605_wp, 2.42075463_wp, 1.46579772_wp, 0.00000000_wp, & ! 33
      & 2.89652631_wp, 2.45421858_wp, 2.27883625_wp, 0.00000000_wp, & ! 34
      & 3.28921099_wp, 2.56526915_wp, 1.64501640_wp, 0.00000000_wp, & ! 35
      & 5.20988189_wp, 2.84336725_wp, 2.75838814_wp, 0.00000000_wp, & ! 36
      & 1.26972917_wp, 1.88730596_wp, 0.00000000_wp, 0.00000000_wp, & ! 37
      & 1.86880714_wp, 1.78546342_wp, 2.16012236_wp, 0.00000000_wp, & ! 38
      & 0.92001877_wp, 1.45732462_wp, 2.22901395_wp, 0.00000000_wp, & ! 39
      & 6.50647305_wp, 1.43202338_wp, 2.11971490_wp, 0.00000000_wp, & ! 40
      & 2.10973371_wp, 2.79944781_wp, 2.01897369_wp, 0.00000000_wp, & ! 41
      & 2.58413333_wp, 3.02795359_wp, 2.08733665_wp, 0.00000000_wp, & ! 42
      & 2.62141555_wp, 3.13487625_wp, 2.13259872_wp, 0.00000000_wp, & ! 43
      & 2.73984499_wp, 2.18167834_wp, 2.54609647_wp, 0.00000000_wp, & ! 44
      & 1.84057176_wp, 2.97482636_wp, 3.10693700_wp, 0.00000000_wp, & ! 45
      & 1.75622839_wp, 3.39424756_wp, 3.20265306_wp, 0.00000000_wp, & ! 46
      & 3.05018811_wp, 2.34951987_wp, 3.35332952_wp, 0.00000000_wp, & ! 47
      & 2.41999128_wp, 2.28892954_wp, 0.00000000_wp, 0.00000000_wp, & ! 48
      & 2.87813961_wp, 2.44659724_wp, 2.75773502_wp, 0.00000000_wp, & ! 49
      & 3.03823214_wp, 2.32082155_wp, 1.77513328_wp, 0.00000000_wp, & ! 50
      & 2.68750711_wp, 2.38565373_wp, 2.12596190_wp, 0.00000000_wp, & ! 51
      & 2.81071790_wp, 2.45274786_wp, 2.01871821_wp, 0.00000000_wp, & ! 52
      & 2.90686956_wp, 2.49377102_wp, 1.90073732_wp, 0.00000000_wp, & ! 53
      & 4.17531340_wp, 2.86937955_wp, 2.96894812_wp, 0.00000000_wp, & ! 54
      & 1.24299361_wp, 1.99142040_wp, 0.00000000_wp, 0.00000000_wp, & ! 55
      & 1.31400366_wp, 1.16438481_wp, 2.12759606_wp, 0.00000000_wp, & ! 56
      & 2.81737350_wp, 1.69863323_wp, 2.27369715_wp, 0.00000000_wp, & ! 57
      & 2.84503901_wp, 1.46018192_wp, 2.53498936_wp, 0.00000000_wp, & ! 58
      & 2.81697107_wp, 1.47545307_wp, 2.54350275_wp, 0.00000000_wp, & ! 59
      & 2.78890313_wp, 1.49072422_wp, 2.55201615_wp, 0.00000000_wp, & ! 60
      & 2.76083520_wp, 1.50599537_wp, 2.56052955_wp, 0.00000000_wp, & ! 61
      & 2.73276726_wp, 1.52126652_wp, 2.56904294_wp, 0.00000000_wp, & ! 62
      & 2.70469932_wp, 1.53653767_wp, 2.57755634_wp, 0.00000000_wp, & ! 63
      & 2.67663138_wp, 1.55180881_wp, 2.58606974_wp, 0.00000000_wp, & ! 64
      & 2.64856345_wp, 1.56707996_wp, 2.59458313_wp, 0.00000000_wp, & ! 65
      & 2.62049551_wp, 1.58235111_wp, 2.60309653_wp, 0.00000000_wp, & ! 66
      & 2.59242757_wp, 1.59762226_wp, 2.61160992_wp, 0.00000000_wp, & ! 67
      & 2.56435964_wp, 1.61289341_wp, 2.62012332_wp, 0.00000000_wp, & ! 68
      & 2.53629170_wp, 1.62816456_wp, 2.62863672_wp, 0.00000000_wp, & ! 69
      & 2.50822376_wp, 1.64343571_wp, 2.63715011_wp, 0.00000000_wp, & ! 70
      & 2.48015583_wp, 1.65870685_wp, 2.64566351_wp, 0.00000000_wp, & ! 71
      & 3.19537752_wp, 2.24853837_wp, 2.41492177_wp, 0.00000000_wp, & ! 72
      & 3.14122020_wp, 2.48723489_wp, 2.21933576_wp, 0.00000000_wp, & ! 73
      & 3.17661283_wp, 3.39538568_wp, 2.37502789_wp, 0.00000000_wp, & ! 74
      & 3.14538352_wp, 2.58361113_wp, 2.47139347_wp, 0.00000000_wp, & ! 75
      & 1.81565647_wp, 2.48106221_wp, 3.18585355_wp, 0.00000000_wp, & ! 76
      & 2.11798490_wp, 2.85857032_wp, 3.47048400_wp, 0.00000000_wp, & ! 77
      & 2.71241232_wp, 3.37886078_wp, 3.64124964_wp, 0.00000000_wp, & ! 78
      & 2.80572458_wp, 2.82570220_wp, 3.72064445_wp, 0.00000000_wp, & ! 79
      & 2.61951362_wp, 2.69607886_wp, 0.00000000_wp, 0.00000000_wp, & ! 80
      & 3.05383193_wp, 2.61683803_wp, 3.32179612_wp, 0.00000000_wp, & ! 81
      & 3.02135073_wp, 2.59250246_wp, 4.24674489_wp, 0.00000000_wp, & ! 82
      & 3.16405210_wp, 2.63238785_wp, 3.04625573_wp, 0.00000000_wp, & ! 83
      & 2.96133467_wp, 2.71388453_wp, 2.31022562_wp, 0.00000000_wp, & ! 84
      & 2.98240599_wp, 2.95960758_wp, 2.43778345_wp, 0.00000000_wp, & ! 85
      & 3.07936232_wp, 2.68589775_wp, 2.10311395_wp, 0.00000000_wp, & ! 86
      & 1.81913220_wp, 3.23064408_wp, 0.00000000_wp, 0.00000000_wp, & ! 87
      & 2.43263729_wp, 2.47485608_wp, 2.09113715_wp, 0.00000000_wp, & ! 88
      & 3.65108887_wp, 3.45440279_wp, 1.97314608_wp, 1.98901892_wp, & ! 89
      & 3.35816295_wp, 2.81245896_wp, 2.05947820_wp, 2.04247660_wp, & ! 90
      & 3.08262439_wp, 2.24936413_wp, 2.13696560_wp, 2.09705269_wp, & ! 91
      & 2.82447317_wp, 1.76511830_wp, 2.20560830_wp, 2.15274719_wp, & ! 92
      & 2.58370931_wp, 1.35972146_wp, 2.26540630_wp, 2.20956009_wp, & ! 93
      & 2.36033280_wp, 1.03317362_wp, 2.31635959_wp, 2.26749140_wp, & ! 94
      & 2.15434364_wp, 0.78547477_wp, 2.35846817_wp, 2.32654112_wp, & ! 95
      & 1.96574183_wp, 0.61662492_wp, 2.39173204_wp, 2.38670925_wp, & ! 96
      & 1.79452738_wp, 0.52662407_wp, 2.41615121_wp, 2.44799578_wp, & ! 97
      & 1.64070027_wp, 0.51547221_wp, 2.43172568_wp, 2.51040072_wp, & ! 98
      & 1.50426052_wp, 0.58316935_wp, 2.43845543_wp, 2.57392407_wp, & ! 99
      & 1.38520812_wp, 0.72971549_wp, 2.43634048_wp, 2.63856583_wp, & ! 100
      & 1.28354307_wp, 0.95511062_wp, 2.42538083_wp, 2.70432599_wp, & ! 101
      & 1.19926537_wp, 1.25935475_wp, 2.40557646_wp, 2.77120456_wp, & ! 102
      & 1.13237503_wp, 1.64244787_wp, 2.37692740_wp, 2.83920154_wp],& ! 103
      & shape(zeta))

      integer :: isp, izp, ish, stat
      integer, allocatable :: nshell(:)
      type(cgto_type), allocatable :: cgto(:, :)

      nshell = nsh(mol%num)
      allocate(cgto(maxval(nshell), mol%nid))
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, nshell(isp)
            call slater_to_gauss(ng, pqn(ish, izp), lsh(ish, izp), zeta(ish, izp), &
            & cgto(ish, isp), .true., stat)
         end do
      end do

      call new_basis(bas, mol, nshell, cgto, 1.0_wp)

   end subroutine make_basis

   subroutine test_scaled_selfenergy_mol(error, mol, ref)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: ref(:)

      type(basis_type) :: bas
      type(tb_hamiltonian) :: h0
      type(erf_ncoord_type) :: ncoord
      type(erf_en_ncoord_type) :: ncoord_en
      type(adjacency_list) :: list
      real(wp), parameter :: cn_cutoff = 30.0_wp
      real(wp), allocatable :: lattr(:, :), cn(:), cn_en(:), rcov(:), en(:)
      real(wp), allocatable :: selfenergy(:)
      real(wp) :: cutoff
      integer :: ii

      call make_basis(bas, mol, 6)

      call check(error, bas%nsh, size(ref, 1))
      if (allocated(error)) return

      call new_hamiltonian(h0, mol, bas, ceh_h0spec(mol))

      allocate(cn(mol%nat), cn_en(mol%nat), rcov(mol%nid), en(mol%nid))
      ! test with the standard Pyykko radii and Pauling EN (not as in CEH parametrization)
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)
      call new_erf_ncoord(ncoord, mol, cutoff=cn_cutoff, rcov=rcov)
      call new_erf_en_ncoord(ncoord_en, mol, cutoff=cn_cutoff, rcov=rcov)
      call get_lattice_points(mol%periodic, mol%lattice, cn_cutoff, lattr)
      call get_coordination_number(ncoord , mol, lattr, cn_cutoff, cn)
      call get_coordination_number(ncoord_en, mol, lattr, cn_cutoff, cn_en)

      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      call new_adjacency_list(list, mol, lattr, cutoff)

      allocate(selfenergy(bas%nsh))
      call get_scaled_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, cn=cn, cn_en=cn_en, &
      & selfenergy=selfenergy)

      do ii = 1, size(selfenergy, 1)
        call check(error, selfenergy(ii), ref(ii), thr=thr2)
        if (allocated(error)) then
           print '(2es20.13)', selfenergy(ii), ref(ii)
           return
        end if
      end do

   end subroutine test_scaled_selfenergy_mol


   subroutine test_hamiltonian_mol(error, mol, ref)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: ref(:, :)

      type(basis_type) :: bas
      type(tb_hamiltonian) :: h0
      type(erf_ncoord_type) :: ncoord
      type(erf_en_ncoord_type) :: ncoord_en
      type(adjacency_list) :: list
      real(wp), parameter :: cn_cutoff = 30.0_wp
      real(wp), allocatable :: lattr(:, :), cn(:), cn_en(:), rcov(:), en(:)
      real(wp), allocatable :: overlap(:, :), overlap_diat(:, :), dpint(:, :, :)
      real(wp), allocatable :: hamiltonian(:, :), selfenergy(:)
      real(wp) :: cutoff
      integer :: ii, jj

      call make_basis(bas, mol, 6)

      call check(error, bas%nao, size(ref, 1))
      if (allocated(error)) return

      call new_hamiltonian(h0, mol, bas, ceh_h0spec(mol))

      allocate(cn(mol%nat), cn_en(mol%nat), rcov(mol%nid), en(mol%nid))
      ! test with the standard Pyykko radii and Pauling EN (not as in CEH parametrization)
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)
      call new_erf_ncoord(ncoord, mol, cutoff=cn_cutoff, rcov=rcov)
      call new_erf_en_ncoord(ncoord_en, mol, cutoff=cn_cutoff, rcov=rcov)
      call get_lattice_points(mol%periodic, mol%lattice, cn_cutoff, lattr)
      call get_coordination_number(ncoord , mol, lattr, cn_cutoff, cn)
      call get_coordination_number(ncoord_en, mol, lattr, cn_cutoff, cn_en)

      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      call new_adjacency_list(list, mol, lattr, cutoff)

      allocate(selfenergy(bas%nsh))
      call get_scaled_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, cn=cn, cn_en=cn_en, &
      & selfenergy=selfenergy)

      allocate(overlap(bas%nao, bas%nao), overlap_diat(bas%nao, bas%nao), &
      & dpint(3, bas%nao, bas%nao), hamiltonian(bas%nao, bas%nao))

      call get_hamiltonian(mol, lattr, list, bas, h0, selfenergy, &
      & overlap, overlap_diat, dpint, hamiltonian)

      do ii = 1, size(hamiltonian, 2)
         do jj = 1, size(hamiltonian, 1)
            call check(error, hamiltonian(jj, ii), ref(jj, ii), thr=thr2)
            if (allocated(error)) then
               print '(3es21.13)', hamiltonian(jj, ii), ref(jj, ii), &
               & hamiltonian(jj, ii) - ref(jj, ii)
               return
            end if
         end do
      end do

   end subroutine test_hamiltonian_mol

   subroutine test_overlap_diat_mol(error, mol, ref)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: ref(:, :)

      type(basis_type) :: bas
      type(tb_hamiltonian) :: h0
      type(erf_ncoord_type) :: ncoord
      type(erf_en_ncoord_type) :: ncoord_en
      type(adjacency_list) :: list
      real(wp), parameter :: cn_cutoff = 30.0_wp
      real(wp), allocatable :: lattr(:, :), cn(:), cn_en(:), rcov(:), en(:)
      real(wp), allocatable :: overlap(:, :), overlap_diat(:, :), dpint(:, :, :)
      real(wp), allocatable :: hamiltonian(:, :), selfenergy(:)
      real(wp) :: cutoff
      integer :: ii, jj

      call make_basis(bas, mol, 6)

      call check(error, bas%nao, size(ref, 1))
      if (allocated(error)) return

      call new_hamiltonian(h0, mol, bas, ceh_h0spec(mol))

      allocate(cn(mol%nat), cn_en(mol%nat), rcov(mol%nid), en(mol%nid))
      ! test with the standard Pyykko radii and Pauling EN (not as in CEH parametrization)
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)
      call new_erf_ncoord(ncoord, mol, cutoff=cn_cutoff, rcov=rcov)
      call new_erf_en_ncoord(ncoord_en, mol, cutoff=cn_cutoff, rcov=rcov)
      call get_lattice_points(mol%periodic, mol%lattice, cn_cutoff, lattr)
      call get_coordination_number(ncoord , mol, lattr, cn_cutoff, cn)
      call get_coordination_number(ncoord_en, mol, lattr, cn_cutoff, cn_en)

      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      call new_adjacency_list(list, mol, lattr, cutoff)

      allocate(selfenergy(bas%nsh))
      call get_scaled_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, cn=cn, cn_en=cn_en, &
      & selfenergy=selfenergy)

      allocate(overlap(bas%nao, bas%nao), overlap_diat(bas%nao, bas%nao), &
      & dpint(3, bas%nao, bas%nao), hamiltonian(bas%nao, bas%nao))

      call get_hamiltonian(mol, lattr, list, bas, h0, selfenergy, &
      & overlap, overlap_diat, dpint, hamiltonian)

      do ii = 1, size(overlap_diat, 2)
         do jj = 1, size(overlap_diat, 1)
            call check(error, overlap_diat(jj, ii), ref(jj, ii), thr=thr2)
            if (allocated(error)) then
               print '(3es21.13)', overlap_diat(jj, ii), ref(jj, ii), &
               & overlap_diat(jj, ii) - ref(jj, ii)
               return
            end if
         end do
      end do

   end subroutine test_overlap_diat_mol

   subroutine test_q_gen(error, mol, ref)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type), intent(in) :: mol

      !> Reference CEH charges
      real(wp), intent(in) :: ref(:)

      type(context_type) :: ctx
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
      real(wp), allocatable :: cn(:)
      real(wp), parameter :: accuracy = 1e-8_wp
      integer :: i
      allocate(cn(mol%nat))

      call new_ceh_calculator(calc, mol, error)
      if (allocated(error)) return
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      ctx%verbosity = 0
      call ceh_singlepoint(ctx, calc, mol, wfn, accuracy)
      if (ctx%failed()) then
         call ctx%get_error(error)
         return
      end if

      do i = 1, mol%nat
         call check(error, wfn%qat(i,1), ref(i), thr=1e-6_wp)
         if (allocated(error)) then
            print '(3es21.13)',  wfn%qat(i,1), ref(i), &
            & wfn%qat(i,1) - ref(i)
            return
         end if
      enddo

   end subroutine test_q_gen


   subroutine test_scaled_selfenergy_h2(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nsh = 2
      real(wp), parameter :: scaled_selfenergy(nsh) = reshape([&
      & -5.2057326046758E-01_wp, -5.2057326046758E-01_wp & 
      &],shape(scaled_selfenergy))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "H2")
      call test_scaled_selfenergy_mol(error, mol, scaled_selfenergy)

   end subroutine test_scaled_selfenergy_h2

   subroutine test_scaled_selfenergy_lih(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nsh = 3
      real(wp), parameter :: scaled_selfenergy(nsh) = reshape([&
      & -5.7614182696741E-02_wp, -1.3057703854461E-01_wp, -3.6985761349230E-01_wp &
      &],shape(scaled_selfenergy))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "LiH")
      call test_scaled_selfenergy_mol(error, mol, scaled_selfenergy)

   end subroutine test_scaled_selfenergy_lih

   subroutine test_scaled_selfenergy_s2(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nsh = 6
      real(wp), parameter :: scaled_selfenergy(nsh) = reshape([&
      & -6.9008304496671E-01_wp, -5.6274208401578E-01_wp, -5.7343694597688E-02_wp, & 
      & -6.9008304496671E-01_wp, -5.6274208401578E-01_wp, -5.7343694597688E-02_wp &
      &], shape(scaled_selfenergy))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "S2")
      call test_scaled_selfenergy_mol(error, mol, scaled_selfenergy)

   end subroutine test_scaled_selfenergy_s2

   subroutine test_scaled_selfenergy_sih4(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nsh = 7
      real(wp), parameter :: scaled_selfenergy(nsh) = reshape([&
      & -7.0849504464403E-01_wp, -4.7605638972741E-01_wp, -1.8541704653682E-01_wp, & 
      & -4.8652644697198E-01_wp, -4.8652644697198E-01_wp, -4.8652644697198E-01_wp, &
      & -4.8652644697198E-01_wp], shape(scaled_selfenergy))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "SiH4")
      call test_scaled_selfenergy_mol(error, mol, scaled_selfenergy)

   end subroutine test_scaled_selfenergy_sih4

   subroutine test_scaled_selfenergy_accl6(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nsh = 22
      real(wp), parameter :: scaled_selfenergy(nsh) = reshape([&
      & 2.73904391613233_wp    , -4.2731759819762E-01_wp, -5.4115371441002E-01_wp, &
      &-1.67143635735189E-02_wp, -4.3819259611367E-01_wp, -2.9082233002044E-01_wp, &
      & 7.02743423215177E-02_wp, -4.3819425191345E-01_wp, -2.9084728314198E-01_wp, &
      & 7.02681834728156E-02_wp, -4.3822006265352E-01_wp, -2.9084720659706E-01_wp, &
      & 7.02683075671165E-02_wp, -4.3822174122876E-01_wp, -2.9085944552609E-01_wp, &
      & 7.02652903250950E-02_wp, -4.3819504378374E-01_wp, -2.9083806237387E-01_wp, &
      & 7.02704650332720E-02_wp, -4.3818833398193E-01_wp, -2.9082921173350E-01_wp, &
      & 7.02726245699396E-02_wp], shape(scaled_selfenergy))

      type(structure_type) :: mol

      call get_structure(mol, "f-block", "AcCl6")
      call test_scaled_selfenergy_mol(error, mol, scaled_selfenergy)

   end subroutine test_scaled_selfenergy_accl6


   subroutine test_hamiltonian_h2(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 2
      real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -5.20573260467584E-01_wp,  -5.12556505508839E-01_wp,  -5.12556505508839E-01_wp, &
      & -5.20573260467584E-01_wp],shape(hamiltonian))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "H2")
      call test_hamiltonian_mol(error, mol, hamiltonian)

   end subroutine test_hamiltonian_h2

   subroutine test_hamiltonian_lih(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 5
      real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -5.76141826967410E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.64597683251088E-01_wp,  0.00000000000000E+00_wp, &
      & -1.30577038544615E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.30577038544615E-01_wp,  0.00000000000000E+00_wp, -2.17667285277957E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.30577038544615E-01_wp,  0.00000000000000E+00_wp, -1.64597683251088E-01_wp, &
      &  0.00000000000000E+00_wp, -2.17667285277957E-01_wp,  0.00000000000000E+00_wp, &
      & -3.69857613492308E-01_wp],shape(hamiltonian))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "LiH")
      call test_hamiltonian_mol(error, mol, hamiltonian)

   end subroutine test_hamiltonian_lih

   subroutine test_hamiltonian_s2(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 18
      real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -6.90083044966714E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.08003200438262E-01_wp,  0.00000000000000E+00_wp,  2.54418424031753E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.19476852196354E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -5.62742084015788E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.36072924446956E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  9.36704132556010E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.62742084015788E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.54418424031753E-01_wp,  0.00000000000000E+00_wp,  4.02200104061939E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.56606945484210E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -5.62742084015788E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.36072924446956E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  9.36704132556010E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -5.73436945976887E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -5.31497588937989E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.73436945976887E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -9.36704132556010E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  2.34814881206664E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -5.73436945976887E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.19476852196354E-01_wp,  0.00000000000000E+00_wp,  1.56606945484210E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.50747597557354E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -5.73436945976887E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -9.36704132556010E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.34814881206664E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.73436945976887E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.31497588937989E-03_wp, &
      & -1.08003200438262E-01_wp,  0.00000000000000E+00_wp, -2.54418424031753E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.19476852196354E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -6.90083044966714E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.36072924446956E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -9.36704132556010E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -5.62742084015788E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  2.54418424031753E-01_wp,  0.00000000000000E+00_wp,  4.02200104061939E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  1.56606945484210E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.62742084015788E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.36072924446956E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -9.36704132556010E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -5.62742084015788E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -5.31497588937989E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -5.73436945976887E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  9.36704132556010E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  2.34814881206664E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.73436945976887E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.19476852196354E-01_wp,  0.00000000000000E+00_wp, -1.56606945484210E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.50747597557354E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -5.73436945976887E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  9.36704132556010E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.34814881206664E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -5.73436945976887E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.31497588937989E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.73436945976887E-02_wp],&
      & shape(hamiltonian))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "S2")
      call test_hamiltonian_mol(error, mol, hamiltonian)

   end subroutine test_hamiltonian_s2

   subroutine test_hamiltonian_sih4(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 13
      real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -7.08495044644038E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.38605666294581E-01_wp, -3.38605666294581E-01_wp, -3.38605666294581E-01_wp, &
      & -3.38605666294581E-01_wp,  0.00000000000000E+00_wp, -4.76056389727416E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.31521047829169E-01_wp,  2.31521047829169E-01_wp, &
      &  2.31521047829169E-01_wp, -2.31521047829169E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -4.76056389727416E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  2.31521047829169E-01_wp, &
      &  2.31521047829169E-01_wp, -2.31521047829169E-01_wp, -2.31521047829169E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -4.76056389727416E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.31521047829169E-01_wp,  2.31521047829169E-01_wp, -2.31521047829169E-01_wp, &
      &  2.31521047829169E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -1.85417046536822E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.59872557670323E-01_wp, -1.59872557670323E-01_wp, &
      &  1.59872557670323E-01_wp,  1.59872557670323E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.85417046536822E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.59872557670323E-01_wp, &
      & -1.59872557670323E-01_wp,  1.59872557670323E-01_wp, -1.59872557670323E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.85417046536822E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.72825675197766E-17_wp, -2.72825675197766E-17_wp, -2.72825675197766E-17_wp, &
      & -2.72825675197766E-17_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -1.85417046536822E-01_wp, &
      &  0.00000000000000E+00_wp,  1.59872557670323E-01_wp, -1.59872557670323E-01_wp, &
      & -1.59872557670323E-01_wp,  1.59872557670323E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.85417046536822E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.38605666294581E-01_wp, -2.31521047829169E-01_wp,  2.31521047829169E-01_wp, &
      & -2.31521047829169E-01_wp, -1.59872557670323E-01_wp,  1.59872557670323E-01_wp, &
      & -2.72825675197766E-17_wp,  1.59872557670323E-01_wp,  0.00000000000000E+00_wp, &
      & -4.86526446971988E-01_wp, -4.43220553944223E-02_wp, -4.43220553944223E-02_wp, &
      & -4.43220553944223E-02_wp, -3.38605666294581E-01_wp,  2.31521047829169E-01_wp, &
      &  2.31521047829169E-01_wp,  2.31521047829169E-01_wp, -1.59872557670323E-01_wp, &
      & -1.59872557670323E-01_wp, -2.72825675197766E-17_wp, -1.59872557670323E-01_wp, &
      &  0.00000000000000E+00_wp, -4.43220553944223E-02_wp, -4.86526446971988E-01_wp, &
      & -4.43220553944223E-02_wp, -4.43220553944223E-02_wp, -3.38605666294581E-01_wp, &
      &  2.31521047829169E-01_wp, -2.31521047829169E-01_wp, -2.31521047829169E-01_wp, &
      &  1.59872557670323E-01_wp,  1.59872557670323E-01_wp, -2.72825675197766E-17_wp, &
      & -1.59872557670323E-01_wp,  0.00000000000000E+00_wp, -4.43220553944223E-02_wp, &
      & -4.43220553944223E-02_wp, -4.86526446971988E-01_wp, -4.43220553944223E-02_wp, &
      & -3.38605666294581E-01_wp, -2.31521047829169E-01_wp, -2.31521047829169E-01_wp, &
      &  2.31521047829169E-01_wp,  1.59872557670323E-01_wp, -1.59872557670323E-01_wp, &
      & -2.72825675197766E-17_wp,  1.59872557670323E-01_wp,  0.00000000000000E+00_wp, &
      & -4.43220553944223E-02_wp, -4.43220553944223E-02_wp, -4.43220553944223E-02_wp, &
      & -4.86526446971988E-01_wp],shape(hamiltonian))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "SiH4")
      call test_hamiltonian_mol(error, mol, hamiltonian)

   end subroutine test_hamiltonian_sih4

   subroutine test_overlap_diat_h2(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 2
      real(wp), parameter :: overlap_diat(nao, nao) = reshape([&
      &  9.9999999988149E-01_wp,  1.5435023369101E+00_wp,  1.5435023369101E+00_wp,&
      &  9.9999999988149E-01_wp],shape(overlap_diat))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "H2")
      call test_overlap_diat_mol(error, mol, overlap_diat)

   end subroutine test_overlap_diat_h2

   subroutine test_overlap_diat_lih(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 5
      real(wp), parameter :: overlap_diat(nao, nao) = reshape([&
      &  1.0000000000060E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  1.2072400803924E+00_wp,  0.0000000000000E+00_wp,&
      &  9.9999999992569E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  9.9999999992569E-01_wp,  0.0000000000000E+00_wp,  1.0945054381789E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  9.9999999992569E-01_wp,  0.0000000000000E+00_wp,  1.2072400803924E+00_wp,&
      &  0.0000000000000E+00_wp,  1.0945054381789E+00_wp,  0.0000000000000E+00_wp,&
      &  9.9999999988149E-01_wp],shape(overlap_diat))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "LiH")
      call test_overlap_diat_mol(error, mol, overlap_diat)

   end subroutine test_overlap_diat_lih

   subroutine test_overlap_diat_s2(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 18
      real(wp), parameter :: overlap_diat(nao, nao) = reshape([&
      &  9.9999999986933E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  2.4534808116926E-01_wp,  0.0000000000000E+00_wp, -5.1101098902195E-01_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  3.5172654233339E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  9.9999999999806E-01_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  2.5407518750216E-01_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp, -2.8345529314904E-01_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  9.9999999999806E-01_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  5.1101098902195E-01_wp,  0.0000000000000E+00_wp, -7.5098751105891E-01_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  4.7390703316608E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  9.9999999999806E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  2.5407518750216E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp, -2.8345529314904E-01_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  9.9999999983021E-01_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  7.8547726772879E-02_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  9.9999999983021E-01_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  2.8345529314904E-01_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp, -3.4702274319026E-01_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  9.9999999983021E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  3.5172654233339E-01_wp,  0.0000000000000E+00_wp, -4.7390703316608E-01_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  3.7056901464493E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  9.9999999983021E-01_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  2.8345529314904E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp, -3.4702274319026E-01_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  9.9999999983021E-01_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  7.8547726772879E-02_wp,&
      &  2.4534808116926E-01_wp,  0.0000000000000E+00_wp,  5.1101098902195E-01_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  3.5172654233339E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  9.9999999986933E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  2.5407518750216E-01_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  2.8345529314904E-01_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  9.9999999999806E-01_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      & -5.1101098902195E-01_wp,  0.0000000000000E+00_wp, -7.5098751105891E-01_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      & -4.7390703316608E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  9.9999999999806E-01_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  2.5407518750216E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  2.8345529314904E-01_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  9.9999999999806E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  7.8547726772879E-02_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  9.9999999983021E-01_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp, -2.8345529314904E-01_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp, -3.4702274319026E-01_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  9.9999999983021E-01_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  3.5172654233339E-01_wp,  0.0000000000000E+00_wp,  4.7390703316608E-01_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  3.7056901464493E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  9.9999999983021E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      & -2.8345529314904E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp, -3.4702274319026E-01_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  9.9999999983021E-01_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  7.8547726772879E-02_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  9.9999999983021E-01_wp],&
      & shape(overlap_diat))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "S2")
      call test_overlap_diat_mol(error, mol, overlap_diat)

   end subroutine test_overlap_diat_s2

   subroutine test_overlap_diat_sih4(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 13
      real(wp), parameter :: overlap_diat(nao, nao) = reshape([&
      &  9.9999999986933E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  8.8837412170621E-01_wp,  8.8837412170621E-01_wp,  8.8837412170621E-01_wp,&
      &  8.8837412170621E-01_wp,  0.0000000000000E+00_wp,  9.9999999999806E-01_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  6.0523565223965E-01_wp, -6.0523565223965E-01_wp,&
      & -6.0523565223965E-01_wp,  6.0523565223965E-01_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  9.9999999999806E-01_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp, -6.0523565223965E-01_wp,&
      & -6.0523565223965E-01_wp,  6.0523565223965E-01_wp,  6.0523565223965E-01_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  9.9999999999806E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  6.0523565223965E-01_wp, -6.0523565223965E-01_wp,  6.0523565223965E-01_wp,&
      & -6.0523565223965E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  9.9999999983021E-01_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  5.2351746131409E-01_wp,  5.2351746131409E-01_wp,&
      & -5.2351746131409E-01_wp, -5.2351746131409E-01_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  9.9999999983021E-01_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp, -5.2351746131409E-01_wp,&
      &  5.2351746131409E-01_wp, -5.2351746131409E-01_wp,  5.2351746131409E-01_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  9.9999999983021E-01_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  8.9339288081803E-17_wp,  8.9339288081803E-17_wp,  8.9339288081803E-17_wp,&
      &  8.9339288081803E-17_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  9.9999999983021E-01_wp,&
      &  0.0000000000000E+00_wp, -5.2351746131409E-01_wp,  5.2351746131409E-01_wp,&
      &  5.2351746131409E-01_wp, -5.2351746131409E-01_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  9.9999999983021E-01_wp,  0.0000000000000E+00_wp,&
      &  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,  0.0000000000000E+00_wp,&
      &  8.8837412170621E-01_wp,  6.0523565223965E-01_wp, -6.0523565223965E-01_wp,&
      &  6.0523565223965E-01_wp,  5.2351746131409E-01_wp, -5.2351746131409E-01_wp,&
      &  8.9339288081803E-17_wp, -5.2351746131409E-01_wp,  0.0000000000000E+00_wp,&
      &  9.9999999988149E-01_wp,  1.4281072933032E-01_wp,  1.4281072933032E-01_wp,&
      &  1.4281072933032E-01_wp,  8.8837412170621E-01_wp, -6.0523565223965E-01_wp,&
      & -6.0523565223965E-01_wp, -6.0523565223965E-01_wp,  5.2351746131409E-01_wp,&
      &  5.2351746131409E-01_wp,  8.9339288081803E-17_wp,  5.2351746131409E-01_wp,&
      &  0.0000000000000E+00_wp,  1.4281072933032E-01_wp,  9.9999999988149E-01_wp,&
      &  1.4281072933032E-01_wp,  1.4281072933032E-01_wp,  8.8837412170621E-01_wp,&
      & -6.0523565223965E-01_wp,  6.0523565223965E-01_wp,  6.0523565223965E-01_wp,&
      & -5.2351746131409E-01_wp, -5.2351746131409E-01_wp,  8.9339288081803E-17_wp,&
      &  5.2351746131409E-01_wp,  0.0000000000000E+00_wp,  1.4281072933032E-01_wp,&
      &  1.4281072933032E-01_wp,  9.9999999988149E-01_wp,  1.4281072933032E-01_wp,&
      &  8.8837412170621E-01_wp,  6.0523565223965E-01_wp,  6.0523565223965E-01_wp,&
      & -6.0523565223965E-01_wp, -5.2351746131409E-01_wp,  5.2351746131409E-01_wp,&
      &  8.9339288081803E-17_wp, -5.2351746131409E-01_wp,  0.0000000000000E+00_wp,&
      &  1.4281072933032E-01_wp,  1.4281072933032E-01_wp,  1.4281072933032E-01_wp,&
      &  9.9999999988149E-01_wp],shape(overlap_diat))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "SiH4")
      call test_overlap_diat_mol(error, mol, overlap_diat)

   end subroutine test_overlap_diat_sih4

   subroutine test_q_h2(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(2) = reshape([ & 
      & 0.0000000000000_wp, 0.00000000000000_wp &
      &], shape(charges))

      call get_structure(mol, "MB16-43", "H2")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_h2

   subroutine test_q_lih(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(2) = reshape([ &
      & 0.452928818506870_wp, -0.452928818506872_wp &
      &], shape(charges))

      call get_structure(mol, "MB16-43", "LiH")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_lih

   subroutine test_q_sih4(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(5) = reshape([ &
      & 0.268915553151106_wp, -0.067228888686111_wp, -0.067228888686110_wp, &
      &-0.067228888686110_wp, -0.067228888686110_wp], shape(charges))

      call get_structure(mol, "MB16-43", "SiH4")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_sih4

   subroutine test_q_cecl3(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(4) = reshape([ &
      & 0.941257219125013_wp, -0.312389885339237_wp, -0.316670447603892_wp, &
      &-0.312196886181877_wp], shape(charges))

      call get_structure(mol, "f-block", "CeCl3")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_cecl3

   subroutine test_q_accl6(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(7) = reshape([ &
      &  0.288810442976471_wp, -0.04820498331048_wp, -0.048312400345871_wp, &
      & -0.047844815375011_wp, -0.04788094096742_wp, -0.048243946234083_wp, &
      & -0.048323356740988_wp], shape(charges))

      call get_structure(mol, "f-block", "AcCl6")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_accl6

   subroutine test_q_panp(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(2) = reshape([ &
      & -0.490095851470747_wp, 0.490095852181341_wp], shape(charges))

      call get_structure(mol, "f-block", "PaNp")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_panp

   subroutine test_q_mb01(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(16) = reshape([ &
      &  0.50933743182523_wp, -0.063950757122457_wp, -0.447396210062547_wp, &     
      & -0.06007626073194_wp, -0.228995347028063_wp,  0.081936572631240_wp, &
      & -0.04029707489635_wp, -0.384822906853029_wp, -0.214508333206973_wp, &     
      &  0.14648324095015_wp,  0.090840217217610_wp,  0.034875957186194_wp, &
      & -0.05930815144452_wp,  0.133798380818110_wp, -0.063944989141738_wp, &
      &  0.56602822987298_wp], shape(charges))

      call get_structure(mol, "MB16-43", "01")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_mb01

   subroutine test_q_mb02(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(16) = reshape([ &
      &-0.085528757393776_wp, -0.069997806147367_wp, -0.234351594245719_wp, &
      &-0.167766349660902_wp,  0.432632782697259_wp,  0.171300109658293_wp, &
      &-0.088439445405627_wp, -0.045209869368921_wp,  0.391614741682698_wp, &
      & 0.148323670914553_wp, -0.084649073627366_wp,  0.371028349885743_wp, &
      &-0.347496611788866_wp, -0.086873903658569_wp, -0.001912803785713_wp, &
      &-0.302673439755712_wp], shape(charges))

      call get_structure(mol, "MB16-43", "02")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_mb02

   subroutine test_q_mb03(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(16) = reshape([ &
      & 0.079458459812547_wp, -0.514984986958578_wp,  0.027119814956294_wp, &
      & 0.294936546309030_wp,  0.396051282334921_wp,  0.033859363521586_wp, &
      &-0.261951753885935_wp,  0.026314869404998_wp,  0.037212079428442_wp, &
      &-0.005505860141304_wp, -0.364181487304452_wp, -0.142392172388313_wp, &     
      &-0.286416123354516_wp,  0.100899051623118_wp,  0.558735071392839_wp, &     
      & 0.020845845249307_wp], shape(charges))

      call get_structure(mol, "MB16-43", "03")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_mb03

   subroutine test_q_mb04(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(16) = reshape([ &
      &-0.018039040440414_wp, -0.194798993283821_wp,  -0.076478151155481_wp, &
      &-0.167853095449230_wp,  0.304619635473503_wp,  -0.022689763513518_wp, &
      & 0.019714590105775_wp,  0.000713111871502_wp,  -0.047985458556249_wp, &
      &-0.104183443069051_wp, -0.161265309374659_wp,   0.538802940610920_wp, &
      & 0.190213172949517_wp, -0.317990923047223_wp,   0.034501918582603_wp, &
      & 0.022718808299253_wp], shape(charges))

      call get_structure(mol, "MB16-43", "04")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_mb04

   subroutine test_q_ef_chrg_mb01(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
      real(wp), parameter :: accuracy = 1e-8_wp
      class(container_type), allocatable :: cont      
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: ref(16) = reshape([ &
      &-5.42237346896788_wp, -0.77304500586496_wp,   2.5895850175165_wp, &
      &-0.92233780581096_wp,  6.99267602990832_wp,   0.4742366118103_wp, &
      &-0.11849846722517_wp,  4.22307140408149_wp,   1.5873873640455_wp, &
      & 0.31672778030780_wp,  0.99906183347402_wp, -10.5405405662106_wp, &
      &-3.80217066006454_wp,  1.92138378495190_wp,   3.8481284909192_wp, &
      & 0.62670765712891_wp], shape(ref))

      real(wp) :: efield(3)
      integer :: i

      efield = 0.0_wp
      efield(3) = 0.2_wp

      call get_structure(mol, "MB16-43", "01")
      mol%charge = 2.0_wp
      call new_ceh_calculator(calc, mol, error)
      if (allocated(error)) return
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      cont = electric_field(efield)
      call calc%push_back(cont)
      ctx%verbosity = 0
      call ceh_singlepoint(ctx, calc, mol, wfn, accuracy)
      if (ctx%failed()) then
         call ctx%get_error(error)
         return
      end if

      do i = 1, mol%nat
         call check(error, wfn%qat(i,1), ref(i), thr=5e-6_wp, message="Calculated charge&
         & does not match reference")
         if (allocated(error)) return
      enddo

   end subroutine test_q_ef_chrg_mb01

   subroutine test_d_mb01(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
      real(wp) :: dipole(3), tmp(3)
      real(wp), parameter :: accuracy = 1e-8_wp
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: ref(3) = reshape([ &
      & 0.201892497728508_wp, -1.15519893399684_wp, -1.91423938957019_wp &
      &], shape(ref))

      call get_structure(mol, "MB16-43", "01")

      call new_ceh_calculator(calc, mol, error)
      if (allocated(error)) return
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      ctx%verbosity = 0
      call ceh_singlepoint(ctx, calc, mol, wfn, accuracy)
      if (ctx%failed()) then
         call ctx%get_error(error)
         return
      end if
      tmp = 0.0_wp
      dipole = 0.0_wp
      call gemv(mol%xyz, wfn%qat(:, 1), tmp)
      dipole(:) = tmp + sum(wfn%dpat(:, :, 1), 2)

      if (any(abs(dipole - ref) > 1e-5_wp)) then
         call test_failed(error, "Numerical dipole moment does not match")
         print '(3es21.14)', dipole
         print '("---")'
         print '(3es21.14)', ref
         print '("---")'
         print '(3es21.14)', dipole - ref
      end if

   end subroutine test_d_mb01

   subroutine test_d_field_mb04(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
      class(container_type), allocatable :: cont
      real(wp) :: energy, efield(3), dipole(3), tmp(3)
      real(wp), parameter :: accuracy = 1e-8_wp
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: ref(3) = reshape([ & 
         & -7.6402587223855_wp,  83.5065044491344_wp,  0.55047274934631_wp &
         & ], shape(ref))

      call get_structure(mol, "MB16-43", "04")
      energy = 0.0_wp
      efield = 0.0_wp
      efield(2) = 0.2_wp

      call new_ceh_calculator(calc, mol, error)
      if (allocated(error)) return
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)

      cont = electric_field(efield)
      call calc%push_back(cont)

      ctx%verbosity = 0
      call ceh_singlepoint(ctx, calc, mol, wfn, accuracy)
      if (ctx%failed()) then
         call ctx%get_error(error)
         return
      end if
      tmp = 0.0_wp
      dipole = 0.0_wp
      call gemv(mol%xyz, wfn%qat(:, 1), tmp)
      dipole(:) = tmp + sum(wfn%dpat(:, :, 1), 2)

      if (any(abs(dipole - ref) > 1e-5_wp)) then
         call test_failed(error, "Numerical dipole moment does not match")
         print '(3es21.14)', dipole
         print '("---")'
         print '(3es21.14)', ref
         print '("---")'
         print '(3es21.14)', dipole - ref
      end if

   end subroutine test_d_field_mb04

   subroutine test_d_hcn(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(structure_type) :: mol1,mol2
      type(xtb_calculator) :: calc1,calc2
      type(wavefunction_type) :: wfn1,wfn2
      class(container_type), allocatable :: cont1,cont2
      real(wp), parameter :: accuracy = 1e-8_wp
      real(wp) :: efield(3), dip1(3), dip2(3), tmp(3)
      integer, parameter :: num(3) = reshape([ &
         7, &
         6, &
         1], shape(num))
      integer, parameter :: nat = 3
      real(wp) :: xyz(3, nat) = reshape([ &
      & -0.09604091224796_wp,  0.0_wp, 0.0_wp, &
      &  2.09604091224796_wp,  0.0_wp, 0.0_wp, &
      &  4.10859879422050_wp,  0.0_wp, 0.0_wp], &
      & shape(xyz))

      ctx%verbosity = 0
      call new(mol1, num, xyz)
      efield = 0.0_wp
      efield(1) = -0.1_wp
      call new_ceh_calculator(calc1, mol1, error)
      if (allocated(error)) return
      call new_wavefunction(wfn1, mol1%nat, calc1%bas%nsh, calc1%bas%nao, 1, kt)
      cont1 = electric_field(efield)
      call calc1%push_back(cont1)
      call ceh_singlepoint(ctx, calc1, mol1, wfn1, accuracy)
      if (ctx%failed()) then
         call ctx%get_error(error)
         return
      end if
      tmp = 0.0_wp
      dip1 = 0.0_wp
      call gemv(mol1%xyz, wfn1%qat(:, 1), tmp)
      dip1(:) = tmp + sum(wfn1%dpat(:, :, 1), 2)

      xyz(1, :) = xyz(1, :) - 1.0_wp
      call new(mol2, num, xyz)
      call new_ceh_calculator(calc2, mol2, error)
      if (allocated(error)) return
      call new_wavefunction(wfn2, mol2%nat, calc2%bas%nsh, calc2%bas%nao, 1, kt)
      cont2 = electric_field(efield)
      call calc2%push_back(cont2)
      call ceh_singlepoint(ctx, calc2, mol2, wfn2, accuracy)
      if (ctx%failed()) then
         call ctx%get_error(error)
         return
      end if
      tmp = 0.0_wp
      dip2 = 0.0_wp
      call gemv(mol2%xyz, wfn2%qat(:, 1), tmp)
      dip2(:) = tmp + sum(wfn2%dpat(:, :, 1), 2)

      if (any(abs(dip1 - dip2) > 1e-7_wp)) then
         call test_failed(error, "Numerical dipole moment does not match")
         print '(3es21.14)', dip1
         print '("---")'
         print '(3es21.14)', dip2
         print '("---")'
         print '(3es21.14)', dip1 - dip2
      end if

   end subroutine test_d_hcn
end module test_ceh
