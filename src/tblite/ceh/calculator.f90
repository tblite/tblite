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

!> @file tblite/ceh/calculator.f90
!> Provides the calculator type for holding CEH Hamiltonian parametrization.

!> Implementation of calculator type for the Charge-Extended Hückel Hamiltonian.
module tblite_ceh_calculator
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_basis_ortho, only : orthogonalize
   use tblite_basis_type, only : basis_type, new_basis, cgto_type
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_container, only : container_type, container_list
   use tblite_ncoord, only : ncoord_type, new_ncoord
   use tblite_param, only : param_record
   use tblite_ceh_h0, only : ceh_hamiltonian
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
      class(ncoord_type), allocatable :: ncoordstd, ncoorden
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

      if (allocated(self%interactions)) then
         str = str // nl // indent // self%interactions%info(verbosity, indent)
      end if
   end function info


end module tblite_ceh_calculator
