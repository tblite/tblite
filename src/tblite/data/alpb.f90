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
!> Provides ALPB/GBSA parameters

!> ALPB/GBSA parameters
module tblite_data_alpb
   use mctc_env, only : wp
   use mctc_io_symbols, only : to_number
   use tblite_solvation_alpb, only: alpb_input
   use mctc_io, only : structure_type
   implicit none
   private

   public :: get_alpb_param


   !> Get spin constant for species
   !interface get_alpb_param
   !   module procedure :: get_alpb_param
   !end interface get_alpb_param

 
   real(wp) :: sx(1:94) = [0.18678116_wp, 1.99854836_wp, 0.50934487_wp, 0.30000000_wp, 0.93372749_wp, &
 0.73948749_wp, 0.77003311_wp, 0.30000000_wp, 0.62227524_wp, 1.22892076_wp, &
 0.65065466_wp, 1.15301804_wp, 1.00517744_wp, 0.93204996_wp, 0.80000000_wp, &
 0.30000000_wp, 0.80000000_wp, 0.90930774_wp, 0.30000000_wp, 1.13114016_wp, &
 1.06981655_wp, 1.42134411_wp, 1.09146204_wp, 1.14487481_wp, 0.93238901_wp, &
 1.66030803_wp, 1.03734066_wp, 1.21328994_wp, 0.96459170_wp, 0.70218231_wp, &
 0.94458397_wp, 0.87849331_wp, 1.04104293_wp, 1.01890919_wp, 1.42471942_wp, &
 1.31471665_wp, 0.30000000_wp, 0.99506730_wp, 1.31750068_wp, 0.91790577_wp, &
 1.13352069_wp, 0.90365194_wp, 0.97192416_wp, 1.13693166_wp, 1.13908614_wp, &
 1.03723586_wp, 0.99373444_wp, 0.57959099_wp, 0.80138068_wp, 0.45705909_wp, &
 1.12709787_wp, 0.80213876_wp, 1.37998347_wp, 1.10387200_wp, 0.30000000_wp, &
 0.30000000_wp, 0.80000000_wp, 0.80000000_wp, 0.80000000_wp, 0.80000000_wp, &
 0.80000000_wp, 0.80000000_wp, 0.80000000_wp, 0.80000000_wp, 0.80000000_wp, &
 0.80000000_wp, 0.80000000_wp, 0.80000000_wp, 0.80000000_wp, 0.80000000_wp, &
 0.80000000_wp, 1.11415978_wp, 0.95102453_wp, 1.31024463_wp, 1.21914682_wp, &
 0.91559959_wp, 0.99579666_wp, 0.90140483_wp, 1.09420509_wp, 0.02609846_wp, &
 0.49510247_wp, 0.69445869_wp, 0.54304437_wp, 0.80000000_wp, 0.80000000_wp, &
 0.80000000_wp, 0.80000000_wp, 0.80000000_wp, 0.80000000_wp, 0.80000000_wp, &
 0.80000000_wp, 0.80000000_wp, 0.80000000_wp, 0.80000000_wp]


contains


!> Get ALPB/GBSA parameters
subroutine get_alpb_param(input, mol)
   !> Input of ALPB
   type(alpb_input), intent(inout) :: input
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> iterator
   integer :: i

   !> born scale for GFN2/water
   input%born_scale = 1.47438678_wp

   !> born offset for GFN2/water 
   input%born_offset = 0.00000000_wp

   !> set descreening -> loop over nat and assign sx
   do i = 1, mol%nat
      input%descreening(i) = sx(mol%num(i))
   end do

   return

end subroutine get_alpb_param

end module tblite_data_alpb
