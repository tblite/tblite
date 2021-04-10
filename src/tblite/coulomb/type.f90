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

!> Definition of the abstract base class for electrostatic and coulombic interactions
module tblite_coulomb_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_scf_info, only : scf_info
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: coulomb_type

   type, abstract :: coulomb_type
   contains
      procedure(update), deferred :: update
      procedure(get_energy), deferred :: get_energy
      procedure(get_potential), deferred :: get_potential
      procedure(get_gradient), deferred :: get_gradient
   end type coulomb_type

   abstract interface
      subroutine update(self, mol, cache)
         import :: coulomb_type, structure_type, coulomb_cache
         !> Instance of the electrostatic container
         class(coulomb_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Reusable data container
         type(coulomb_cache), intent(inout) :: cache
      end subroutine update

      subroutine get_energy(self, mol, cache, wfn, energy)
         import :: coulomb_type, structure_type, coulomb_cache, wavefunction_type, wp
         !> Instance of the electrostatic container
         class(coulomb_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Reusable data container
         type(coulomb_cache), intent(inout) :: cache
         !> Wavefunction data
         type(wavefunction_type), intent(in) :: wfn
         !> Electrostatic energy
         real(wp), intent(inout) :: energy
      end subroutine get_energy

      subroutine get_potential(self, mol, cache, wfn, pot)
         import :: coulomb_type, structure_type, coulomb_cache, wavefunction_type, &
            & potential_type
         !> Instance of the electrostatic container
         class(coulomb_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Reusable data container
         type(coulomb_cache), intent(inout) :: cache
         !> Wavefunction data
         type(wavefunction_type), intent(in) :: wfn
         !> Density dependent potential
         type(potential_type), intent(inout) :: pot
      end subroutine get_potential

      subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
         import :: coulomb_type, structure_type, coulomb_cache, wavefunction_type, wp
         !> Instance of the electrostatic container
         class(coulomb_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Reusable data container
         type(coulomb_cache), intent(inout) :: cache
         !> Wavefunction data
         type(wavefunction_type), intent(in) :: wfn
         !> Molecular gradient of the repulsion energy
         real(wp), contiguous, intent(inout) :: gradient(:, :)
         !> Strain derivatives of the repulsion energy
         real(wp), contiguous, intent(inout) :: sigma(:, :)
      end subroutine get_gradient
   end interface

end module tblite_coulomb_type
