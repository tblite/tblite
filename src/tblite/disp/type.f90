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

!> Definition of abstract base class for dispersion interactions
module tblite_disp_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_disp_cache, only : dispersion_cache
   use tblite_scf_info, only : scf_info
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: dispersion_type

   !> Abstract base class for dispersion interactions
   type, abstract :: dispersion_type
   contains
      !> Update dispersion cache
      procedure :: update
      !> Get information about density dependent quantities used in the energy
      procedure :: variable_info
      !> Evaluate non-selfconsistent part of the dispersion correction
      procedure :: get_engrad
      !> Evaluate selfconsistent energy of the dispersion correction
      procedure :: get_energy
      !> Evaluate charge dependent potential shift from the dispersion correction
      procedure :: get_potential
      !> Evaluate gradient contributions from the selfconsistent dispersion correction
      procedure :: get_gradient
   end type dispersion_type


contains


!> Update dispersion cache
subroutine update(self, mol, cache)
   !> Instance of the dispersion correction
   class(dispersion_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(dispersion_cache), intent(inout) :: cache
end subroutine update


!> Evaluate non-selfconsistent part of the dispersion correction
subroutine get_engrad(self, mol, cache, energy, gradient, sigma)
   !> Instance of the dispersion correction
   class(dispersion_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(dispersion_cache), intent(inout) :: cache
   !> Dispersion energy
   real(wp), intent(inout) :: energy
   !> Dispersion gradient
   real(wp), contiguous, intent(inout), optional :: gradient(:, :)
   !> Dispersion virial
   real(wp), contiguous, intent(inout), optional :: sigma(:, :)
end subroutine get_engrad


!> Evaluate selfconsistent energy of the dispersion correction
subroutine get_energy(self, mol, cache, wfn, energy)
   !> Instance of the dispersion correction
   class(dispersion_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(dispersion_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Dispersion energy
   real(wp), intent(inout) :: energy
end subroutine get_energy


!> Evaluate charge dependent potential shift from the dispersion correction
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the dispersion correction
   class(dispersion_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(dispersion_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot
end subroutine get_potential


!> Evaluate gradient contributions from the selfconsistent dispersion correction
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the dispersion correction
   class(dispersion_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(dispersion_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Dispersion gradient
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Dispersion virial
   real(wp), contiguous, intent(inout) :: sigma(:, :)
end subroutine get_gradient


!> Get information about density dependent quantities used in the energy
pure function variable_info(self) result(info)
   !> Instance of the electrostatic container
   class(dispersion_type), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info()
end function variable_info


end module tblite_disp_type
