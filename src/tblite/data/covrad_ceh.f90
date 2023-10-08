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

!> @file tblite/data/covrad_ceh.f90
!> Provides covalent radii for all elements (up to 86)
!> Covalent radii were fitted together with other empirical parameters # MM, August 02, 2023

module tblite_data_covrad_ceh
   use mctc_env, only : wp
   use mctc_io_convert, only : aatoau
   use mctc_io_symbols, only : to_number
   implicit none
   private

   public :: get_covalent_cehrad


   !> Covalent radii for DFT-D3 coordination number
   interface get_covalent_cehrad
      module procedure :: get_covalent_rad_num
      module procedure :: get_covalent_rad_sym
   end interface get_covalent_cehrad


   integer, parameter :: max_elem = 86

   !> Empirical radii for the CEH CN # MM, August 02, 2023
   !> already in a.u.
   real(wp), parameter :: ceh_cov_radii(max_elem) = [&
   &  0.57035470_wp,  0.63663000_wp,  1.61183561_wp,  1.47896082_wp,  1.43072572_wp,  1.43413810_wp, &
&  1.35274709_wp,  1.25722907_wp,  1.13811103_wp,  1.58691253_wp,  2.11082024_wp,  2.03286435_wp, &
&  1.99093514_wp,  2.12405887_wp,  2.06378172_wp,  1.95815921_wp,  1.80110747_wp,  1.98297456_wp, &
&  2.53436208_wp,  2.38331220_wp,  2.32183367_wp,  2.27983350_wp,  2.36227384_wp,  2.28992230_wp, &
&  2.08884537_wp,  1.95721844_wp,  2.00494356_wp,  1.90496208_wp,  1.95121459_wp,  1.97995506_wp, &
&  2.03020326_wp,  2.12653262_wp,  2.21608183_wp,  2.30941815_wp,  2.04703550_wp,  2.27500482_wp, &
&  2.73713344_wp,  2.52355454_wp,  2.51484190_wp,  2.51383609_wp,  2.52894575_wp,  2.71783178_wp, &
&  2.47383671_wp,  2.18217476_wp,  2.06022477_wp,  2.08877605_wp,  2.19362364_wp,  2.24348356_wp, &
&  2.32088123_wp,  2.41681804_wp,  2.52967611_wp,  2.51654972_wp,  2.37162867_wp,  2.25373961_wp, &
&  2.86931103_wp,  2.78052368_wp,  2.77183103_wp,  2.58463810_wp,  2.57071497_wp,  2.55679184_wp, &
&  2.54286872_wp,  2.52894559_wp,  2.51502246_wp,  2.50109933_wp,  2.48717620_wp,  2.47325307_wp, &
&  2.45932995_wp,  2.44540682_wp,  2.43148369_wp,  2.41756056_wp,  2.40363743_wp,  2.67061205_wp, &
&  2.61590595_wp,  2.69225275_wp,  2.50456586_wp,  2.18330586_wp,  2.08963987_wp,  2.16181230_wp, &
&  2.15789163_wp,  2.30990980_wp,  2.42191712_wp,  2.54016038_wp,  2.61105824_wp,  2.62500484_wp, &
&  2.58044220_wp,  2.55638791_wp ]

contains


!> Get covalent radius for a given element symbol
elemental function get_covalent_rad_sym(sym) result(rad)

   !> Element symbol
   character(len=*), intent(in) :: sym

   !> Covalent radius
   real(wp) :: rad

   rad = get_covalent_cehrad(to_number(sym))

end function get_covalent_rad_sym


!> Get covalent radius for a given atomic number
elemental function get_covalent_rad_num(num) result(rad)

   !> Atomic number
   integer, intent(in) :: num

   !> Covalent radius
   real(wp) :: rad

   if (num > 0 .and. num <= size(ceh_cov_radii)) then
      rad = ceh_cov_radii(num)
   else
      rad = 0.0_wp
   end if

end function get_covalent_rad_num


end module tblite_data_covrad_ceh
