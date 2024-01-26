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

   public :: get_spin_constant, read_spin_constants

   !> Get spin constant for species
   interface get_spin_constant
      module procedure :: get_spin_constant_symbol
      module procedure :: get_spin_constant_number
   end interface get_spin_constant

   integer, parameter :: lidx(0:2, 0:2) = reshape(& ! ss sp sd  sp pp pd  sd pd dd
      & [1, 2, 4,  2, 3, 5,  4, 5, 6], shape(lidx))

   real(wp) :: spin_constants(6, 86) = reshape([&  ! ss, sp, pp, sd, pd, dd
      -0.0716750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0866000_wp,-0.0386380_wp,-0.0673250_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0177750_wp,-0.0139380_wp,-0.0180500_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0229250_wp,-0.0186130_wp,-0.0175750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0270750_wp,-0.0218880_wp,-0.0195750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0305750_wp,-0.0249870_wp,-0.0226250_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0330750_wp,-0.0274880_wp,-0.0254750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0352000_wp,-0.0295380_wp,-0.0278500_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0369000_wp,-0.0312000_wp,-0.0299000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0384000_wp,-0.0326750_wp,-0.0317250_wp,-0.0142250_wp,-0.0153000_wp,-0.0413750_wp, &
      -0.0151250_wp,-0.0133500_wp,-0.0229250_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0165750_wp,-0.0131380_wp,-0.0174000_wp,-0.0094000_wp,-0.0179750_wp,-0.0223750_wp, &
      -0.0181500_wp,-0.0138870_wp,-0.0139750_wp,-0.0081250_wp,-0.0116870_wp,-0.0128000_wp, &
      -0.0195250_wp,-0.0150500_wp,-0.0143500_wp,-0.0083870_wp,-0.0116250_wp,-0.0140000_wp, &
      -0.0205500_wp,-0.0161250_wp,-0.0149250_wp,-0.0093000_wp,-0.0119620_wp,-0.0147250_wp, &
      -0.0213250_wp,-0.0170000_wp,-0.0155250_wp,-0.0100500_wp,-0.0121750_wp,-0.0148750_wp, &
      -0.0216500_wp,-0.0176620_wp,-0.0160500_wp,-0.0109630_wp,-0.0126500_wp,-0.0150750_wp, &
      -0.0221500_wp,-0.0184250_wp,-0.0166500_wp,-0.0118370_wp,-0.0131630_wp,-0.0153000_wp, &
      -0.0106750_wp,-0.0108870_wp,-0.0164750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0118000_wp,-0.0104000_wp,-0.0134500_wp,-0.0055000_wp,-0.0034750_wp,-0.0102000_wp, &
      -0.0128500_wp,-0.0109250_wp,-0.0139500_wp,-0.0047870_wp,-0.0024750_wp,-0.0125000_wp, &
      -0.0134250_wp,-0.0112250_wp,-0.0146500_wp,-0.0043380_wp,-0.0020000_wp,-0.0138750_wp, &
      -0.0139500_wp,-0.0115120_wp,-0.0152500_wp,-0.0039870_wp,-0.0017120_wp,-0.0149250_wp, &
      -0.0145000_wp,-0.0116500_wp,-0.0158000_wp,-0.0037380_wp,-0.0016130_wp,-0.0157250_wp, &
      -0.0150000_wp,-0.0117880_wp,-0.0167250_wp,-0.0034870_wp,-0.0013250_wp,-0.0165250_wp, &
      -0.0154000_wp,-0.0119250_wp,-0.0178500_wp,-0.0033380_wp,-0.0011630_wp,-0.0171750_wp, &
      -0.0158250_wp,-0.0120500_wp,-0.0186000_wp,-0.0031380_wp,-0.0010000_wp,-0.0177750_wp, &
      -0.0161500_wp,-0.0122000_wp,-0.0196000_wp,-0.0029870_wp,-0.0009000_wp,-0.0182750_wp, &
      -0.0165500_wp,-0.0123120_wp,-0.0202250_wp,-0.0028370_wp,-0.0008500_wp,-0.0187750_wp, &
      -0.0168500_wp,-0.0123750_wp,-0.0213500_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0172250_wp,-0.0128000_wp,-0.0134500_wp,-0.0085750_wp,-0.0129880_wp,-0.0158500_wp, &
      -0.0175500_wp,-0.0133500_wp,-0.0135750_wp,-0.0081130_wp,-0.0128130_wp,-0.0175500_wp, &
      -0.0178500_wp,-0.0137750_wp,-0.0136000_wp,-0.0080500_wp,-0.0123620_wp,-0.0175250_wp, &
      -0.0180000_wp,-0.0141250_wp,-0.0136500_wp,-0.0081130_wp,-0.0119620_wp,-0.0173000_wp, &
      -0.0181000_wp,-0.0144120_wp,-0.0136750_wp,-0.0082750_wp,-0.0117630_wp,-0.0167750_wp, &
      -0.0181500_wp,-0.0146630_wp,-0.0138000_wp,-0.0087380_wp,-0.0118000_wp,-0.0165250_wp, &
      -0.0095750_wp,-0.0096000_wp,-0.0167250_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0106750_wp,-0.0093250_wp,-0.0124500_wp,-0.0074870_wp,-0.0059130_wp,-0.0079500_wp, &
      -0.0114500_wp,-0.0097370_wp,-0.0134750_wp,-0.0072750_wp,-0.0047250_wp,-0.0089750_wp, &
      -0.0120000_wp,-0.0098880_wp,-0.0160250_wp,-0.0066500_wp,-0.0036370_wp,-0.0098750_wp, &
      -0.0126000_wp,-0.0102750_wp,-0.0190750_wp,-0.0060500_wp,-0.0029370_wp,-0.0104500_wp, &
      -0.0128750_wp,-0.0105120_wp,-0.0222500_wp,-0.0055500_wp,-0.0024250_wp,-0.0109000_wp, &
      -0.0131500_wp,-0.0106870_wp,-0.0246250_wp,-0.0051250_wp,-0.0020250_wp,-0.0113000_wp, &
      -0.0133750_wp,-0.0107500_wp,-0.0275750_wp,-0.0047630_wp,-0.0016750_wp,-0.0116500_wp, &
      -0.0135500_wp,-0.0108370_wp,-0.0320500_wp,-0.0043870_wp,-0.0014250_wp,-0.0118500_wp, &
      -0.0133500_wp,-0.0111500_wp,-0.0285000_wp,-0.0042380_wp,-0.0013250_wp,-0.0121500_wp, &
      -0.0139250_wp,-0.0110370_wp,-0.0397500_wp,-0.0038880_wp,-0.0010250_wp,-0.0124000_wp, &
      -0.0139750_wp,-0.0105000_wp,-0.0196750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0141250_wp,-0.0105620_wp,-0.0115750_wp,-0.0049370_wp,-0.0093250_wp,-0.0100250_wp, &
      -0.0142750_wp,-0.0109000_wp,-0.0116750_wp,-0.0046880_wp,-0.0091250_wp,-0.0118750_wp, &
      -0.0145250_wp,-0.0110750_wp,-0.0116250_wp,-0.0043620_wp,-0.0087130_wp,-0.0125250_wp, &
      -0.0145250_wp,-0.0112750_wp,-0.0115500_wp,-0.0041880_wp,-0.0082120_wp,-0.0122500_wp, &
      -0.0146500_wp,-0.0113380_wp,-0.0114750_wp,-0.0044500_wp,-0.0083120_wp,-0.0128250_wp, &
      -0.0146250_wp,-0.0114750_wp,-0.0114500_wp,-0.0049130_wp,-0.0086250_wp,-0.0131750_wp, &
      -0.0081750_wp,-0.0085370_wp,-0.0151750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0092500_wp,-0.0082880_wp,-0.0113500_wp,-0.0062880_wp,-0.0042750_wp,-0.0079500_wp, &
      -0.0099000_wp,-0.0084250_wp,-0.0114000_wp,-0.0059250_wp,-0.0033620_wp,-0.0090250_wp, &
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
      -0.0121500_wp,-0.0096250_wp,-0.0125000_wp,-0.0076750_wp,-0.0041750_wp,-0.0103750_wp, &
      -0.0124000_wp,-0.0095750_wp,-0.0134000_wp,-0.0071380_wp,-0.0034750_wp,-0.0109250_wp, &
      -0.0125000_wp,-0.0095370_wp,-0.0144500_wp,-0.0066880_wp,-0.0029250_wp,-0.0112250_wp, &
      -0.0126000_wp,-0.0096120_wp,-0.0148250_wp,-0.0062880_wp,-0.0026120_wp,-0.0114500_wp, &
      -0.0126000_wp,-0.0092000_wp,-0.0206000_wp,-0.0059250_wp,-0.0021250_wp,-0.0115750_wp, &
      -0.0126000_wp,-0.0092750_wp,-0.0209000_wp,-0.0057060_wp,-0.0019750_wp,-0.0116000_wp, &
      -0.0127250_wp,-0.0092120_wp,-0.0222000_wp,-0.0054250_wp,-0.0017500_wp,-0.0117000_wp, &
      -0.0127500_wp,-0.0088880_wp,-0.0306000_wp,-0.0052880_wp,-0.0015370_wp,-0.0117500_wp, &
      -0.0129250_wp,-0.0091880_wp,-0.0291750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0132500_wp,-0.0091120_wp,-0.0107000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0134750_wp,-0.0094120_wp,-0.0109750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0135750_wp,-0.0095250_wp,-0.0109500_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0138500_wp,-0.0096380_wp,-0.0108500_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0137750_wp,-0.0097380_wp,-0.0107250_wp,-0.0025500_wp,-0.0073620_wp,-0.0119250_wp, &
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

!> Read spin constants from external file
subroutine read_spin_constants(filename)

   !> Name of the inputfile
   character(len=*),intent(in) :: filename
   integer :: input 
   integer :: i

   open(newunit=input, file=filename, status='old')
   do i=1,86
     read(input, *) spin_constants(1,i), spin_constants(2,i), spin_constants(3,i), spin_constants(4,i), & 
        & spin_constants(5,i), spin_constants(6,i)
   end do
    

   close(input)

   !> Read from spin_param.txt and modify spin_constants accordingly 

end subroutine read_spin_constants




end module tblite_data_spin
