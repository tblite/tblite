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

!> @file tblite/data/alpb.f90
!> Provides ALPB/GBSA cds parameters

!> ALPB/GBSA parameters
module tblite_data_cds
   use mctc_env, only : wp
   use mctc_io_symbols, only : to_number
   use tblite_solvation_cds, only: cds_input
   use mctc_io, only : structure_type
   use mctc_io_convert, only : aatoau, kcaltoau
   implicit none
   private

   public :: get_cds_param

   real(wp) :: probe = 1.13409020_wp * aatoau

  real(wp) :: tension(1:94) = 1.0e-5_wp * [ &
 -0.08533368_wp,   1.13711777_wp,  -5.83236188_wp, -17.01350889_wp,  -0.82891231_wp, &
 -0.51385188_wp,  -3.24932941_wp,   2.03601297_wp,   1.26697909_wp,   1.00726547_wp, &
 -3.05732760_wp,  -1.61595951_wp,  -0.10749537_wp,   1.73814581_wp,   1.01416364_wp, &
  1.20810520_wp,  -0.06835879_wp,   0.43292409_wp,  -2.69545448_wp,  -9.68055177_wp, &
 -6.86171080_wp,  -2.32183024_wp,  -2.10041202_wp,  -4.15076452_wp,  -2.32406972_wp, &
  9.24779587_wp,   4.48206277_wp,  -1.15972411_wp,  -0.21281688_wp,  -2.64500225_wp, &
 -2.36086956_wp,  -0.46303904_wp,  -0.88149455_wp,   0.23523157_wp,  -0.18620262_wp, &
  0.35105321_wp,  -2.88851792_wp, -11.04377179_wp,  -7.95128133_wp,   3.90876499_wp, &
 -2.48432528_wp,  -3.66936332_wp,  -4.43332314_wp,  -2.94937418_wp,   0.00028790_wp, &
 -0.93137790_wp,  -0.79778296_wp,  -0.92747581_wp,  -2.70394304_wp,  -0.43878679_wp, &
 -0.70393148_wp,  -0.77480977_wp,  -0.76873446_wp,  -0.06431749_wp,  -3.15995511_wp, &
 -5.92642054_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp, &
  0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp, &
  0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp, &
  0.00000000_wp,  -0.91732017_wp,  -0.76038638_wp,   4.55364802_wp,  -5.19397805_wp, &
 -0.97455175_wp,  -0.19095469_wp,   0.37263783_wp,   0.41221465_wp,  -1.77134941_wp, &
 -0.89456867_wp,   0.24690462_wp,   0.62621722_wp,   0.00000000_wp,   0.00000000_wp, &
  0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp, &
  0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp]


  real(wp) :: hbond(1:94) = -kcaltoau * [&
  8.09964704_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp, &
  2.23420044_wp,   1.80633279_wp,   2.22319193_wp,   3.56359195_wp,   0.00000000_wp, &
  0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   3.62926820_wp, &
  0.26212102_wp,   0.15334756_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp, &
  0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp, &
  0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp, &
  0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   2.37162766_wp,   2.02275702_wp, &
  0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp, &
  0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp, &
  0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp, &
  0.00000000_wp,   0.00000000_wp,   2.61716906_wp,   0.00000000_wp,   0.00000000_wp, &
  0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp, &
  0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp, &
  0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp, &
  0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp, &
  0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp, &
  0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp, &
  0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp, &
  0.00000000_wp,   0.00000000_wp,   0.00000000_wp,   0.00000000_wp]**2


contains


!> Get CDS parameters
subroutine get_cds_param(input, mol)
   !> Input of cds
   type(cds_input), intent(inout) :: input
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> iterator
   integer :: i


   if (.not. allocated(input%tension)) then
      allocate(input%tension(mol%nid))
   end if

   if (.not. allocated(input%hbond)) then
      allocate(input%hbond(mol%nid))
   end if

   !> set probe radius for water GFN2 alpb
   input%probe = probe

   !> set tension 
   input%tension = tension(mol%num)

   !> set hbond
   input%hbond = hbond(mol%num)

   return

end subroutine get_cds_param

end module tblite_data_cds
