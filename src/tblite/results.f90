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

!> Container for holding results produced by a calculation
module tblite_results
   use mctc_env, only : wp
   implicit none
   private

   public :: results_type


   !> Results container
   type results_type
      !> Atom-resolved energies
      real(wp), allocatable :: energies(:)
      !> Overlap integrals
      real(wp), allocatable :: overlap(:, :)
      !> (Core) Hamiltonian integrals
      real(wp), allocatable :: hamiltonian(:, :)
   end type results_type

end module tblite_results
