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

!> Definition of the abstract base class for classical interactions.
module tblite_repulsion_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none
   private

   public :: repulsion_type

   !> Abstract base class for classical interactions, like repulsion interactions.
   !>
   !> This class provides a method to retrieve the contributions to the energy,
   !> gradient and virial within a given cutoff.
   type, abstract :: repulsion_type
   contains
      !> Obtain energy, and optionally gradient and virial, for classical contribution
      procedure(get_engrad), deferred :: get_engrad
   end type repulsion_type

   abstract interface
      !> Interface for the evaluation of the classical contribution
      subroutine get_engrad(self, mol, trans, cutoff, energy, gradient, sigma)
         import :: repulsion_type, structure_type, wp
         !> Instance of the repulsion container
         class(repulsion_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Lattice points
         real(wp), intent(in) :: trans(:, :)
         !> Real space cutoff
         real(wp), intent(in) :: cutoff
         !> Repulsion energy
         real(wp), intent(inout) :: energy
         !> Molecular gradient of the repulsion energy
         real(wp), intent(inout), optional :: gradient(:, :)
         !> Strain derivatives of the repulsion energy
         real(wp), intent(inout), optional :: sigma(:, :)
      end subroutine get_engrad
   end interface

end module tblite_repulsion_type
