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

!> @file tblite/results.f90
!> Provides a container for storing additional calculation results.

!> Container for holding results produced by a calculation.
module tblite_results
   use mctc_env, only : wp
   use tblite_double_dictionary, only : double_dictionary_type
   implicit none
   private


   !> Results container
   type, public :: results_type
      !> Atom-resolved energies
      real(wp), allocatable :: energies(:)
      !> Overlap integrals
      real(wp), allocatable :: overlap(:, :)
      !> (Core) Hamiltonian integrals
      real(wp), allocatable :: hamiltonian(:, :)
      !> Wiberg/Mayer bond orders
      real(wp), allocatable :: bond_orders(:, :, :)
      !> post processing values (nat, n_post_proc_labels)
      real(wp), allocatable :: post_proc_values(:, :)
      !> number of ml features
      integer :: n_post_proc_labels = 0
      !> labels of the ml features
      character(len=30),allocatable :: post_proc_labels(:)
      type(double_dictionary_type), allocatable :: dict
   end type results_type

end module tblite_results
