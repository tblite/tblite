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

module tblite_xtb_calculator
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_classical_halogen, only : halogen_correction
   use tblite_disp, only : dispersion_type
   use tblite_ncoord, only : ncoord_type
   use tblite_repulsion_effective, only : tb_repulsion
   use tblite_xtb_coulomb, only : tb_coulomb
   use tblite_xtb_h0, only : tb_hamiltonian
   implicit none
   private

   public :: xtb_calculator

   real(wp), parameter :: mixer_damping_default = 0.4_wp
   integer, parameter :: max_iter_default = 250

   type :: xtb_calculator
      type(basis_type) :: bas
      type(tb_hamiltonian) :: h0
      class(ncoord_type), allocatable :: ncoord
      type(tb_repulsion), allocatable :: repulsion
      type(tb_coulomb), allocatable :: coulomb
      type(halogen_correction), allocatable :: halogen
      class(dispersion_type), allocatable :: dispersion
      real(wp) :: mixer_damping = mixer_damping_default
      integer :: max_iter = max_iter_default
   contains
      procedure :: variable_info
   end type xtb_calculator


contains


subroutine update(self, mol)
   class(xtb_calculator), intent(inout) :: self
   type(structure_type), intent(in) :: mol
end subroutine update


pure function variable_info(self) result(info)
   use tblite_scf_info, only : scf_info, max
   !> Instance of the electrostatic container
   class(xtb_calculator), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info()

   if (allocated(self%coulomb)) then
      info = max(info, self%coulomb%variable_info())
   end if

   if (allocated(self%dispersion)) then
      info = max(info, self%dispersion%variable_info())
   end if

end function variable_info


end module tblite_xtb_calculator
