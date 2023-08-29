
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

!> @file tblite/ceh/h0.f90
!> Provides the Hamiltonian type for CEH.

module tblite_ceh_h0
   use mctc_env, only : wp
   implicit none
   private

   type, public :: ceh_hamiltonian
      !> Atomic level information
      real(wp), allocatable :: hlevel(:)
      !> Reference occupation
      real(wp), allocatable :: refocc(:,:)
   end type ceh_hamiltonian

end module tblite_ceh_h0