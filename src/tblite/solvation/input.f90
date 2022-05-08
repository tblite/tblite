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

module tblite_solvation_input
   use tblite_solvation_alpb, only : alpb_input
   use tblite_solvation_cpcm, only : cpcm_input
   implicit none
   private

   public :: solvation_input


   !> Collection of possible solvation models
   type :: solvation_input
      !> Input for CPCM solvation model
      type(cpcm_input), allocatable :: cpcm
      !> Input for ALPB solvation model
      type(alpb_input), allocatable :: alpb
   end type solvation_input

end module tblite_solvation_input
