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

!> Semiclassical DFT-D3 dispersion correction
module tblite_disp_d3
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_disp_cache, only : dispersion_cache
   use tblite_disp_type, only : dispersion_type
   use dftd3, only : d3_model, new_d3_model, rational_damping_param, realspace_cutoff, &
      & get_dispersion
   implicit none
   private

   public :: d3_dispersion, new_d3_dispersion


   !> Container for DFT-D3 type dispersion correction
   type, extends(dispersion_type) :: d3_dispersion
      type(d3_model) :: model
      type(rational_damping_param) :: param
      type(realspace_cutoff) :: cutoff
   contains
      procedure :: get_engrad
   end type d3_dispersion

contains

subroutine new_d3_dispersion(self, mol, s6, s8, a1, a2, s9)
   !> Instance of the dispersion correction
   type(d3_dispersion), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: s6, s8, a1, a2, s9

   call new_d3_model(self%model, mol)
   self%param = rational_damping_param(s6=s6, s8=s8, s9=s9, a1=a1, a2=a2, alp=14.0_wp)
   self%cutoff = realspace_cutoff(cn=25.0_wp, disp3=25.0_wp, disp2=50.0_wp)
end subroutine new_d3_dispersion

!> Evaluate non-selfconsistent part of the dispersion correction
subroutine get_engrad(self, mol, cache, energy, gradient, sigma)
   !> Instance of the dispersion correction
   class(d3_dispersion), intent(in) :: self
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

   real(wp), allocatable :: tgradient(:, :), tsigma(:, :)

   if (present(gradient) .and. present(sigma)) then
      allocate(tgradient(3, mol%nat), tsigma(3, 3))
   end if

   call get_dispersion(mol, self%model, self%param, self%cutoff, &
      & energy, tgradient, tsigma)

   if (present(gradient) .and. present(sigma)) then
      gradient(:, :) = gradient + tgradient
      sigma(:, :) = sigma + tsigma
   end if
end subroutine get_engrad


end module tblite_disp_d3
