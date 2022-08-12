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

!> @dir tblite/scf
!> Contains the implementation for the self-consistent field iterations

!> @file tblite/scf.f90
!> Reexports the self-consistent field iteration related functionality.

!> Proxy module for rexports from SCF related modules
module tblite_scf
   use tblite_scf_mixer, only : mixer_type, new_mixer
   use tblite_scf_info, only : scf_info
   use tblite_scf_iterator, only : next_scf, get_mixer_dimension
   use tblite_scf_potential, only : potential_type, new_potential
   implicit none
   private

   public :: mixer_type, new_mixer, scf_info, next_scf, get_mixer_dimension, &
      & potential_type, new_potential

end module tblite_scf
