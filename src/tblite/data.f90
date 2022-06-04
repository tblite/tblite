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

!> @dir tblite/data
!> Contains element data used for defining interactions.

!> @file tblite/data.f90
!> Reexports access to element-specific data.

!> Proxy module for providing access to element data.
module tblite_data
   use tblite_data_atomicrad, only : get_atomic_rad
   use tblite_data_covrad, only : get_covalent_rad
   use tblite_data_paulingen, only : get_pauling_en
   use tblite_data_spin, only : get_spin_constant
   implicit none

   public :: get_atomic_rad, get_covalent_rad, get_pauling_en, get_spin_constant
end module tblite_data
