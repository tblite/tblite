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

!> @file tblite/xtb/calculator.f90
!> Provides the calculator type for holding xTB Hamiltonian parametrization.

!> Implementation of calculator type for the extended-tight binding Hamiltonian.
!> The #tblite_xtb_calculator::xtb_calculator collects the basic interactions
!> required to perform a tight-binding calculation.
module tblite_ceh_calculator
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_basis_ortho, only : orthogonalize
   use tblite_basis_type, only : basis_type, new_basis, cgto_type
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_classical_halogen, only : halogen_correction, new_halogen_correction
   use tblite_container, only : container_type, container_list
   use tblite_coulomb_charge, only : coulomb_kernel, new_gamma_coulomb, gamma_coulomb, &
   & new_effective_coulomb, effective_coulomb, average_interface, &
   & harmonic_average, arithmetic_average, geometric_average
   use tblite_coulomb_multipole, only : new_damped_multipole
   use tblite_coulomb_thirdorder, only : new_onsite_thirdorder
   use tblite_disp, only : dispersion_type, d4_dispersion, new_d4_dispersion, &
   & d3_dispersion, new_d3_dispersion
   use tblite_ncoord_ceh, only : ncoord_type_ceh, new_ncoord
   use tblite_param, only : param_record
   use tblite_repulsion, only : new_repulsion
   use tblite_repulsion_effective, only : tb_repulsion
   use tblite_xtb_coulomb, only : tb_coulomb
   use tblite_ceh_h0, only : ceh_hamiltonian
   use tblite_xtb_spec, only : tb_h0spec
   implicit none
   private

   ! public :: new_xtb_calculator
   ! public :: param_h0spec

   !> Extended tight-binding calculator
   type, public :: ceh_calculator
      !> Basis set definition
      type(basis_type) :: bas
      !> Core Hamiltonian
      type(ceh_hamiltonian) :: hamiltonian
      !> Coordination number for modifying the self-energies
      class(ncoord_type_ceh), allocatable :: ncoord
      !> Repulsion energy interactions
      type(tb_repulsion), allocatable :: repulsion
      !> Collection of all Coulombic interactions
      type(tb_coulomb), allocatable :: coulomb
      !> Halogen bonding correction
      type(halogen_correction), allocatable :: halogen
      !> London-dispersion interaction
      class(dispersion_type), allocatable :: dispersion
      !> Store calculated integral intermediates
      logical :: save_integrals = .false.
      !> List of additional interaction containers
      type(container_list), allocatable :: interactions
   contains
      !> Get information about density dependent quantities used in the energy
      procedure :: variable_info
      !> Information on calculator
      procedure :: info
      !> Add an interaction container
      procedure :: push_back
      !> Remove an interaction container
      procedure :: pop
   end type ceh_calculator

contains

!> Add an interaction container
   subroutine push_back(self, cont)
      !> Instance of the tight-binding calculator
      class(ceh_calculator), intent(inout) :: self
      !> Container to be added
      class(container_type), allocatable, intent(inout) :: cont

      if (.not.allocated(self%interactions)) allocate(self%interactions)

      call self%interactions%push_back(cont)
   end subroutine push_back


!> Add a container
   subroutine pop(self, cont, idx)
      !> Instance of the tight-binding calculator
      class(ceh_calculator), intent(inout) :: self
      !> Container to be removed
      class(container_type), allocatable, intent(out) :: cont
      !> Index to remove container from
      integer, intent(in), optional :: idx

      if (.not.allocated(self%interactions)) return

      call self%interactions%pop(cont, idx)
   end subroutine pop


   pure function variable_info(self) result(info)
      use tblite_scf_info, only : scf_info, max
      !> Instance of the electrostatic container
      class(ceh_calculator), intent(in) :: self
      !> Information on the required potential data
      type(scf_info) :: info

      info = scf_info()

      if (allocated(self%coulomb)) then
         info = max(info, self%coulomb%variable_info())
      end if

      if (allocated(self%dispersion)) then
         info = max(info, self%dispersion%variable_info())
      end if

      if (allocated(self%interactions)) then
         info = max(info, self%interactions%variable_info())
      end if

   end function variable_info


!> Information on container
   pure function info(self, verbosity, indent) result(str)
      !> Instance of the interaction container
      class(ceh_calculator), intent(in) :: self
      !> Verbosity level
      integer, intent(in) :: verbosity
      !> Indentation level
      character(len=*), intent(in) :: indent
      !> Information on the container
      character(len=:), allocatable :: str

      character(len=*), parameter :: nl = new_line('a')

      str = "CEH calculator"

      if (allocated(self%repulsion)) then
         str = str // nl // indent // self%repulsion%info(verbosity, indent)
      end if

      if (allocated(self%coulomb)) then
         str = str // nl // indent // self%coulomb%info(verbosity, indent)
      end if

      if (allocated(self%dispersion)) then
         str = str // nl // indent // self%dispersion%info(verbosity, indent)
      end if

      if (allocated(self%interactions)) then
         str = str // nl // indent // self%interactions%info(verbosity, indent)
      end if
   end function info


end module tblite_ceh_calculator
