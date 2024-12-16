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

!> @file tblite/data/spin.f90
!> Provides spin constants for all elements

!> Spin constants
module tblite_data_spin
   use mctc_env, only : wp
   use mctc_io_symbols, only : to_number
   implicit none
   private

   public :: get_spin_constant


   !> Get spin constant for species
   interface get_spin_constant
      module procedure :: get_spin_constant_symbol
      module procedure :: get_spin_constant_number
   end interface get_spin_constant

   integer, parameter :: lidx(0:2, 0:2) = reshape(& ! ss sp sd  sp pp pd  sd pd dd
      & [1, 2, 4,  2, 3, 5,  4, 5, 6], shape(lidx))

   real(wp), parameter :: spin_constants(6, 86) = reshape([&  ! ss, sp, pp, sd, pd, dd
      -0.0716250_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0865500_wp,-0.0386630_wp,-0.0674250_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0178000_wp,-0.0139500_wp,-0.0180500_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0229750_wp,-0.0186250_wp,-0.0175750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0272500_wp,-0.0219370_wp,-0.0195750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0305000_wp,-0.0250250_wp,-0.0226750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0330750_wp,-0.0275000_wp,-0.0254500_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0350750_wp,-0.0295380_wp,-0.0278500_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0369000_wp,-0.0311870_wp,-0.0299250_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0383000_wp,-0.0326250_wp,-0.0317250_wp,-0.0141250_wp,-0.0152500_wp,-0.0413500_wp, &
      -0.0150750_wp,-0.0133370_wp,-0.0229250_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0165000_wp,-0.0130750_wp,-0.0175000_wp,-0.0093750_wp,-0.0179630_wp,-0.0223500_wp, &
      -0.0182500_wp,-0.0138380_wp,-0.0139750_wp,-0.0082250_wp,-0.0117000_wp,-0.0129000_wp, &
      -0.0195250_wp,-0.0150000_wp,-0.0143750_wp,-0.0084500_wp,-0.0116120_wp,-0.0140000_wp, &
      -0.0205750_wp,-0.0161250_wp,-0.0149000_wp,-0.0093000_wp,-0.0119870_wp,-0.0148250_wp, &
      -0.0213250_wp,-0.0170130_wp,-0.0155000_wp,-0.0100370_wp,-0.0121750_wp,-0.0149500_wp, &
      -0.0217500_wp,-0.0177130_wp,-0.0160500_wp,-0.0109750_wp,-0.0126620_wp,-0.0150750_wp, &
      -0.0221500_wp,-0.0183630_wp,-0.0165500_wp,-0.0118870_wp,-0.0131130_wp,-0.0153000_wp, &
      -0.0106500_wp,-0.0109000_wp,-0.0164750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0118000_wp,-0.0104500_wp,-0.0134500_wp,-0.0055000_wp,-0.0035130_wp,-0.0101750_wp, &
      -0.0127250_wp,-0.0108620_wp,-0.0138500_wp,-0.0047880_wp,-0.0024130_wp,-0.0125250_wp, &
      -0.0134250_wp,-0.0112380_wp,-0.0146500_wp,-0.0043380_wp,-0.0019750_wp,-0.0138750_wp, &
      -0.0140750_wp,-0.0114630_wp,-0.0152750_wp,-0.0040500_wp,-0.0017250_wp,-0.0149250_wp, &
      -0.0144750_wp,-0.0116120_wp,-0.0160000_wp,-0.0037250_wp,-0.0014630_wp,-0.0157750_wp, &
      -0.0149000_wp,-0.0118000_wp,-0.0167250_wp,-0.0034870_wp,-0.0013120_wp,-0.0165000_wp, &
      -0.0154000_wp,-0.0120250_wp,-0.0177500_wp,-0.0032880_wp,-0.0011630_wp,-0.0171250_wp, &
      -0.0157000_wp,-0.0120250_wp,-0.0187000_wp,-0.0031500_wp,-0.0010250_wp,-0.0177500_wp, &
      -0.0161500_wp,-0.0122000_wp,-0.0197000_wp,-0.0030370_wp,-0.0009130_wp,-0.0183000_wp, &
      -0.0166500_wp,-0.0123500_wp,-0.0203000_wp,-0.0028250_wp,-0.0008620_wp,-0.0188250_wp, &
      -0.0168500_wp,-0.0123250_wp,-0.0214500_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0172250_wp,-0.0128120_wp,-0.0134000_wp,-0.0085250_wp,-0.0130000_wp,-0.0157750_wp, &
      -0.0174500_wp,-0.0133500_wp,-0.0135500_wp,-0.0081120_wp,-0.0128130_wp,-0.0175250_wp, &
      -0.0178750_wp,-0.0137630_wp,-0.0135750_wp,-0.0080500_wp,-0.0123250_wp,-0.0176500_wp, &
      -0.0180000_wp,-0.0141250_wp,-0.0136250_wp,-0.0081500_wp,-0.0120130_wp,-0.0172000_wp, &
      -0.0181000_wp,-0.0143750_wp,-0.0136750_wp,-0.0082750_wp,-0.0117500_wp,-0.0167750_wp, &
      -0.0181250_wp,-0.0145500_wp,-0.0137000_wp,-0.0086880_wp,-0.0118000_wp,-0.0164250_wp, &
      -0.0095500_wp,-0.0096000_wp,-0.0167250_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0105750_wp,-0.0092870_wp,-0.0125500_wp,-0.0074000_wp,-0.0059000_wp,-0.0079500_wp, &
      -0.0115000_wp,-0.0098250_wp,-0.0136000_wp,-0.0072880_wp,-0.0046250_wp,-0.0090000_wp, &
      -0.0121500_wp,-0.0099880_wp,-0.0162250_wp,-0.0066250_wp,-0.0036250_wp,-0.0098750_wp, &
      -0.0125750_wp,-0.0102620_wp,-0.0191750_wp,-0.0060630_wp,-0.0029250_wp,-0.0104750_wp, &
      -0.0129000_wp,-0.0105000_wp,-0.0222250_wp,-0.0055750_wp,-0.0024250_wp,-0.0109000_wp, &
      -0.0131250_wp,-0.0106250_wp,-0.0247250_wp,-0.0051250_wp,-0.0020120_wp,-0.0113000_wp, &
      -0.0133500_wp,-0.0107620_wp,-0.0276000_wp,-0.0047370_wp,-0.0016750_wp,-0.0116000_wp, &
      -0.0135500_wp,-0.0108380_wp,-0.0320500_wp,-0.0043750_wp,-0.0014250_wp,-0.0118750_wp, &
      -0.0136500_wp,-0.0112120_wp,-0.0287000_wp,-0.0041250_wp,-0.0012880_wp,-0.0121250_wp, &
      -0.0139250_wp,-0.0110500_wp,-0.0397750_wp,-0.0038870_wp,-0.0009630_wp,-0.0124000_wp, &
      -0.0138500_wp,-0.0105000_wp,-0.0196500_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0142250_wp,-0.0105500_wp,-0.0115750_wp,-0.0050370_wp,-0.0093750_wp,-0.0100000_wp, &
      -0.0143000_wp,-0.0108750_wp,-0.0116750_wp,-0.0046880_wp,-0.0090750_wp,-0.0118750_wp, &
      -0.0145250_wp,-0.0111250_wp,-0.0116250_wp,-0.0043750_wp,-0.0087130_wp,-0.0124250_wp, &
      -0.0145250_wp,-0.0112500_wp,-0.0115750_wp,-0.0041870_wp,-0.0081750_wp,-0.0121750_wp, &
      -0.0146500_wp,-0.0113870_wp,-0.0114750_wp,-0.0044620_wp,-0.0083620_wp,-0.0128250_wp, &
      -0.0146500_wp,-0.0114250_wp,-0.0114500_wp,-0.0048750_wp,-0.0085750_wp,-0.0132000_wp, &
      -0.0082000_wp,-0.0085880_wp,-0.0153000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0092500_wp,-0.0083000_wp,-0.0113750_wp,-0.0063870_wp,-0.0042250_wp,-0.0079250_wp, &
      -0.0099000_wp,-0.0084250_wp,-0.0114000_wp,-0.0059370_wp,-0.0033750_wp,-0.0090250_wp, &
      -0.0100000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0100000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0100000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0100000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0100000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0100000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0100000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0100000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0100000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0100000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0100000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0100000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0100000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0100000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0121750_wp,-0.0096750_wp,-0.0126250_wp,-0.0076250_wp,-0.0041130_wp,-0.0104250_wp, &
      -0.0123000_wp,-0.0095750_wp,-0.0134000_wp,-0.0071380_wp,-0.0034630_wp,-0.0109250_wp, &
      -0.0125000_wp,-0.0094620_wp,-0.0144500_wp,-0.0066880_wp,-0.0029130_wp,-0.0112500_wp, &
      -0.0126000_wp,-0.0096620_wp,-0.0148000_wp,-0.0063000_wp,-0.0026130_wp,-0.0114500_wp, &
      -0.0127000_wp,-0.0092000_wp,-0.0205750_wp,-0.0059380_wp,-0.0021120_wp,-0.0115500_wp, &
      -0.0127500_wp,-0.0092750_wp,-0.0209250_wp,-0.0056880_wp,-0.0019120_wp,-0.0116000_wp, &
      -0.0127500_wp,-0.0092250_wp,-0.0222500_wp,-0.0054370_wp,-0.0017870_wp,-0.0117000_wp, &
      -0.0129000_wp,-0.0089380_wp,-0.0306750_wp,-0.0052500_wp,-0.0015000_wp,-0.0117750_wp, &
      -0.0129250_wp,-0.0091870_wp,-0.0292750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0133500_wp,-0.0091120_wp,-0.0107250_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0135750_wp,-0.0094250_wp,-0.0110000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0136750_wp,-0.0095380_wp,-0.0109500_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0137500_wp,-0.0096380_wp,-0.0108500_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0137750_wp,-0.0096750_wp,-0.0107250_wp,-0.0026000_wp,-0.0073630_wp,-0.0119000_wp, &
      -0.0139000_wp,-0.0097380_wp,-0.0106500_wp,-0.0028750_wp,-0.0078120_wp,-0.0130000_wp], &
      shape(spin_constants))

contains


!> Get spin constant for species with a given symbol
elemental function get_spin_constant_symbol(jang, iang, symbol) result(wll)

   !> Angular momentum of shells
   integer, intent(in) :: jang, iang

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> Spin constant
   real(wp) :: wll

   wll = get_spin_constant(jang, iang, to_number(symbol))

end function get_spin_constant_symbol


!> Get spin constant for species with a given atomic number
elemental function get_spin_constant_number(jang, iang, number) result(wll)

   !> Angular momentum of shells
   integer, intent(in) :: jang, iang

   !> Atomic number
   integer, intent(in) :: number

   !> Spin constant
   real(wp) :: wll

   if (number > 0 .and. number <= size(spin_constants, dim=2)) then
      wll = spin_constants(lidx(jang, iang), number)
   else
      wll = -1.0_wp
   end if

end function get_spin_constant_number


end module tblite_data_spin