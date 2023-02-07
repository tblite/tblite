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
      -0.6146750_wp,-0.0335870_wp,-0.1258000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0177750_wp,-0.0139370_wp,-0.0180250_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0228750_wp,-0.0186250_wp,-0.0175750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0273250_wp,-0.0220370_wp,-0.0196000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0305000_wp,-0.0250370_wp,-0.0226750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0330750_wp,-0.0274870_wp,-0.0254500_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0351000_wp,-0.0295500_wp,-0.0278500_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0368000_wp,-0.0312120_wp,-0.0299000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.5500750_wp,-0.1283000_wp,-0.2260000_wp,-0.0119250_wp,-0.0167370_wp,-0.0807250_wp, &
      -0.0151000_wp,-0.0133370_wp,-0.0229250_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0165000_wp,-0.0131750_wp,-0.0175000_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0181500_wp,-0.0138370_wp,-0.0140000_wp,-0.0081870_wp,-0.0116370_wp,-0.0128750_wp, &
      -0.0195250_wp,-0.0150000_wp,-0.0143750_wp,-0.0085000_wp,-0.0116250_wp,-0.0141000_wp, &
      -0.0205750_wp,-0.0160620_wp,-0.0149000_wp,-0.0093120_wp,-0.0119250_wp,-0.0148250_wp, &
      -0.0213250_wp,-0.0170120_wp,-0.0154750_wp,-0.0099870_wp,-0.0121250_wp,-0.0149500_wp, &
      -0.0218500_wp,-0.0177120_wp,-0.0160750_wp,-0.0109870_wp,-0.0126120_wp,-0.0150750_wp, &
      -0.3424750_wp,-0.0778250_wp,-0.1206750_wp,-0.0155000_wp,-0.0230750_wp,-0.0519250_wp, &
      -0.0106750_wp,-0.0109000_wp,-0.0164750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0118000_wp,-0.0104500_wp,-0.0134500_wp,-0.0055250_wp,-0.0035620_wp,-0.0101000_wp, &
      -0.0127250_wp,-0.0108750_wp,-0.0138500_wp,-0.0047870_wp,-0.0024750_wp,-0.0125000_wp, &
      -0.0134250_wp,-0.0111750_wp,-0.0145500_wp,-0.0043500_wp,-0.0019620_wp,-0.0139000_wp, &
      -0.0139750_wp,-0.0114620_wp,-0.0152750_wp,-0.0039870_wp,-0.0017120_wp,-0.0149250_wp, &
      -0.0151750_wp,-0.0124370_wp,-0.0212250_wp,-0.0041750_wp,-0.0016750_wp,-0.0138500_wp, &
      -0.0150000_wp,-0.0117870_wp,-0.0168500_wp,-0.0035000_wp,-0.0013250_wp,-0.0164750_wp, &
      -0.0154000_wp,-0.0119120_wp,-0.0178500_wp,-0.0032870_wp,-0.0011620_wp,-0.0171250_wp, &
      -0.0158250_wp,-0.0120870_wp,-0.0187000_wp,-0.0031370_wp,-0.0010620_wp,-0.0177250_wp, &
      -0.0161500_wp,-0.0121870_wp,-0.0197000_wp,-0.0030370_wp,-0.0009000_wp,-0.0183000_wp, &
      -0.0171250_wp,-0.0131870_wp,-0.0303750_wp,-0.0028370_wp,-0.0006500_wp,-0.0173750_wp, &
      -0.0167750_wp,-0.0123870_wp,-0.0214500_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0172250_wp,-0.0128500_wp,-0.0134250_wp,-0.0085250_wp,-0.0129750_wp,-0.0158500_wp, &
      -0.0175500_wp,-0.0133500_wp,-0.0135750_wp,-0.0080750_wp,-0.0127750_wp,-0.0174250_wp, &
      -0.0178750_wp,-0.0137750_wp,-0.0135500_wp,-0.0080620_wp,-0.0123000_wp,-0.0175250_wp, &
      -0.0180000_wp,-0.0141250_wp,-0.0136250_wp,-0.0081370_wp,-0.0120370_wp,-0.0172000_wp, &
      -0.0181000_wp,-0.0143750_wp,-0.0137250_wp,-0.0082750_wp,-0.0117750_wp,-0.0166750_wp, &
      -0.2990250_wp,-0.0665870_wp,-0.1018750_wp,-0.0125750_wp,-0.0212870_wp,-0.0483000_wp, &
      -0.0095500_wp,-0.0095370_wp,-0.0168250_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0105750_wp,-0.0092870_wp,-0.0124250_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0114250_wp,-0.0094500_wp,-0.0123250_wp,-0.0067620_wp,-0.0039250_wp,-0.0097750_wp, &
      -0.0118500_wp,-0.0096000_wp,-0.0134250_wp,-0.0061250_wp,-0.0030750_wp,-0.0107750_wp, &
      -0.0126000_wp,-0.0102620_wp,-0.0190750_wp,-0.0060500_wp,-0.0028870_wp,-0.0104750_wp, &
      -0.0129250_wp,-0.0104870_wp,-0.0223250_wp,-0.0055500_wp,-0.0024250_wp,-0.0109250_wp, &
      -0.0131250_wp,-0.0106750_wp,-0.0246250_wp,-0.0051120_wp,-0.0020250_wp,-0.0113000_wp, &
      -0.0133500_wp,-0.0107500_wp,-0.0276000_wp,-0.0047500_wp,-0.0017250_wp,-0.0116000_wp, &
      -0.0135250_wp,-0.0108750_wp,-0.0320500_wp,-0.0044000_wp,-0.0014120_wp,-0.0118500_wp, &
      -0.0189750_wp,-0.0239370_wp,-0.1802000_wp,-0.0020870_wp,-0.0014750_wp,-0.0113250_wp, &
      -0.0139000_wp,-0.0110500_wp,-0.0398000_wp,-0.0039000_wp,-0.0010250_wp,-0.0124250_wp, &
      -0.0138500_wp,-0.0105000_wp,-0.0196500_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0141250_wp,-0.0105500_wp,-0.0115750_wp,-0.0050620_wp,-0.0093750_wp,-0.0101000_wp, &
      -0.0143000_wp,-0.0109370_wp,-0.0116750_wp,-0.0046870_wp,-0.0091250_wp,-0.0118750_wp, &
      -0.0144250_wp,-0.0110620_wp,-0.0116250_wp,-0.0043750_wp,-0.0087250_wp,-0.0124250_wp, &
      -0.0145250_wp,-0.0112370_wp,-0.0115250_wp,-0.0041870_wp,-0.0081620_wp,-0.0122500_wp, &
      -0.0145750_wp,-0.0113370_wp,-0.0114500_wp,-0.0044500_wp,-0.0083120_wp,-0.0128250_wp, &
      -0.2557500_wp,-0.0555870_wp,-0.0856250_wp,-0.0046620_wp,-0.0133370_wp,-0.0373500_wp, &
      -0.0081750_wp,-0.0085370_wp,-0.0151750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0093500_wp,-0.0082000_wp,-0.0113750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0099000_wp,-0.0084250_wp,-0.0113750_wp,-0.0059250_wp,-0.0033120_wp,-0.0090250_wp, &
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
      -0.0121750_wp,-0.0096120_wp,-0.0126000_wp,-0.0076500_wp,-0.0041870_wp,-0.0104250_wp, &
      -0.0123250_wp,-0.0095870_wp,-0.0134000_wp,-0.0071500_wp,-0.0034750_wp,-0.0109250_wp, &
      -0.0125000_wp,-0.0095620_wp,-0.0144500_wp,-0.0066750_wp,-0.0029750_wp,-0.0112250_wp, &
      -0.0126000_wp,-0.0096620_wp,-0.0148000_wp,-0.0062620_wp,-0.0026120_wp,-0.0114250_wp, &
      -0.0125750_wp,-0.0092500_wp,-0.0205750_wp,-0.0059370_wp,-0.0021120_wp,-0.0115750_wp, &
      -0.0127250_wp,-0.0092750_wp,-0.0210250_wp,-0.0056750_wp,-0.0019000_wp,-0.0116500_wp, &
      -0.0130750_wp,-0.0102620_wp,-0.0335250_wp,-0.0055620_wp,-0.0017870_wp,-0.0111250_wp, &
      -0.0131750_wp,-0.0100250_wp,-0.0530000_wp,-0.0052750_wp,-0.0015250_wp,-0.0111750_wp, &
      -0.0129250_wp,-0.0092370_wp,-0.0292750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0132750_wp,-0.0091120_wp,-0.0107250_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0134750_wp,-0.0093500_wp,-0.0109750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0136500_wp,-0.0095370_wp,-0.0109750_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0137250_wp,-0.0096370_wp,-0.0108500_wp, 0.0000000_wp, 0.0000000_wp, 0.0000000_wp, &
      -0.0137750_wp,-0.0097370_wp,-0.0107250_wp,-0.0026250_wp,-0.0073620_wp,-0.0119250_wp, &
      -0.2544000_wp,-0.0504000_wp,-0.0806250_wp,-0.0010870_wp,-0.0110500_wp,-0.0351750_wp], &
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
