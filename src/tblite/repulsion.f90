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

!> @dir tblite/repulsion
!> Contains the repulsive interaction implementations.

!> @file tblite/repulsion.f90
!> Reexport for the repulsion container implementations.

!> Proxy module for repulsion interactions.
module tblite_repulsion
   use tblite_repulsion_effective, only : tb_repulsion, new_repulsion
   use tblite_repulsion_type, only : repulsion_type
   implicit none
   private

   public :: tb_repulsion, new_repulsion

end module tblite_repulsion
