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

!> @file tblite/data/paulingen.f90
!> Provides electronegativities for all elements used in the CEH model

!> Pauling electronegativities for the CEH model
module tblite_data_paulingen_ceh
   use mctc_env, only : wp
   use mctc_io_symbols, only : to_number
   implicit none
   private

   public :: get_pauling_en_ceh, pauling_en_ceh


   !> Get electronegativity for a species
   interface get_pauling_en_ceh
      module procedure :: get_pauling_en_ceh_symbol
      module procedure :: get_pauling_en_ceh_number
   end interface get_pauling_en_ceh


   integer, parameter :: max_elem = 86

   ! Pauling EN normalized to EN(F)=1
   ! TM and group 1/2 (from K on) hand optimized
   ! also adjusted: Rn,Xe,Kr,He,LNs
   ! Used for EN-scaled Coordination number in CEH
   real(wp), parameter :: pauling_en_ceh(max_elem) = (1.0_wp/3.98_wp) * [ &
   & 2.200_wp,3.100_wp, & ! H,He
   & 0.980_wp,1.570_wp,2.040_wp,2.550_wp,3.040_wp,3.440_wp,3.980_wp,4.500_wp, & ! Li-Ne
   & 0.930_wp,1.310_wp,1.610_wp,1.900_wp,2.190_wp,2.580_wp,3.160_wp,3.500_wp, & ! Na-Ar
   & 0.700_wp,1.050_wp, & ! K,Ca
   &                   1.280_wp,1.350_wp,1.620_wp,1.710_wp,1.800_wp, & ! Sc-
   &                   1.850_wp,1.930_wp,1.870_wp,1.870_wp,1.600_wp, & ! -Zn
   &                   1.810_wp,2.010_wp,2.180_wp,2.550_wp,2.960_wp,3.200_wp, & ! Ga-Kr
   & 0.700_wp,0.900_wp, & ! Rb,Sr
   &                   1.320_wp,1.380_wp,1.570_wp,1.800_wp,1.900_wp, & ! Y-
   &                   2.180_wp,2.300_wp,2.100_wp,1.800_wp,1.600_wp, & ! -Cd
   &                   1.780_wp,1.960_wp,2.050_wp,2.100_wp,2.660_wp,2.750_wp, & ! In-Xe
   & 0.700_wp,0.800_wp, & ! Cs,Ba
   &          1.050_wp,1.350_wp,1.350_wp,1.350_wp,1.350_wp,1.350_wp,1.350_wp, & ! La-Eu
   &          1.350_wp,1.350_wp,1.350_wp,1.350_wp,1.350_wp,1.350_wp,1.350_wp, & ! Gd-Yb
   &                   1.350_wp,1.350_wp,1.530_wp,1.610_wp,1.730_wp, & ! Lu
   &                   1.920_wp,2.150_wp,2.010_wp,2.000_wp,1.600_wp, & ! -Hg
   &                   1.620_wp,2.330_wp,2.020_wp,2.000_wp,2.200_wp,2.600_wp ] ! Tl-Rn

contains


   !> Get electronegativity for species with a given symbol
   elemental function get_pauling_en_ceh_symbol(symbol) result(en)

      !> Element symbol
      character(len=*), intent(in) :: symbol

      !> atomic EN
      real(wp) :: en

      en = get_pauling_en_ceh(to_number(symbol))

   end function get_pauling_en_ceh_symbol


   !> Get electronegativity for species with a given atomic number
   elemental function get_pauling_en_ceh_number(number) result(en)

      !> Atomic number
      integer, intent(in) :: number

      !> atomic EN
      real(wp) :: en

      if (number > 0 .and. number <= size(pauling_en_ceh, dim=1)) then
         en = pauling_en_ceh(number)
      else
         en = -1.0_wp
      end if

   end function get_pauling_en_ceh_number


end module tblite_data_paulingen_ceh
