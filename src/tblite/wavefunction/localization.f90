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

!> @dir tblite/wavefunction/localization
!> Contains orbital localization method implementations

!> @file tblite/wavefunction/localization.f90
!> Provides reexports of orbital localization related types and procedures

!> Proxy module for wavefunction orbital localization types and procedures
module tblite_wavefunction_localization
   use tblite_wavefunction_localization_fosterboys, only : &
      & fosterboys_localization_type, new_fosterboys_localization
   use tblite_wavefunction_localization_type, only : localization_type, &
      & localization_method, get_orbital_centers, get_localization_id
   implicit none
   public

end module tblite_wavefunction_localization
