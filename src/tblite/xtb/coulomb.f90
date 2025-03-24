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

!> @file tblite/xtb/coulomb.f90
!> Provides a collection for all Coulomb related interactions

!> Collects all Coulomb related interaction in as single interaction container.
module tblite_xtb_coulomb
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_container, only : container_cache, container_type
   use tblite_coulomb, only : coulomb_charge_type, damped_multipole, onsite_thirdorder
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   !> Collection of Coulombic interactions
   type, public, extends(container_type) :: tb_coulomb
      !> Isotroptic second-order electrostatics
      class(coulomb_charge_type), allocatable :: es2
      !> Anisotropic second-order electrostatic
      type(damped_multipole), allocatable :: aes2
      !> Onsite third-order electrostatic
      type(onsite_thirdorder), allocatable :: es3
   contains
      procedure :: update
      procedure :: variable_info
      procedure :: get_energy
      procedure :: get_potential
      procedure :: get_potential_gradient
      procedure :: get_gradient
      procedure :: info
   end type tb_coulomb

contains


subroutine update(self, mol, cache)
   !> Instance of the electrostatic container
   class(tb_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   if (allocated(self%es2)) then
      call self%es2%update(mol, cache)
   end if

   if (allocated(self%aes2)) then
      call self%aes2%update(mol, cache)
   end if

   if (allocated(self%es3)) then
      call self%es3%update(mol, cache)
   end if
end subroutine update


subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the electrostatic container
   class(tb_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Electrostatic energy
   real(wp), intent(inout) :: energies(:)

   if (allocated(self%es2)) then
      call self%es2%get_energy(mol, cache, wfn, energies)
   end if

   if (allocated(self%aes2)) then
      call self%aes2%get_energy(mol, cache, wfn, energies)
   end if

   if (allocated(self%es3)) then
      call self%es3%get_energy(mol, cache, wfn, energies)
   end if
end subroutine get_energy


subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the electrostatic container
   class(tb_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   if (allocated(self%es2)) then
      call self%es2%get_potential(mol, cache, wfn, pot)
   end if

   if (allocated(self%aes2)) then
      call self%aes2%get_potential(mol, cache, wfn, pot)
   end if

   if (allocated(self%es3)) then
      call self%es3%get_potential(mol, cache, wfn, pot)
   end if
end subroutine get_potential


subroutine get_potential_gradient(self, mol, cache, wfn, pot)
   !> Instance of the electrostatic container
   class(tb_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   if (allocated(self%es2)) then
      call self%es2%get_potential_gradient(mol, cache, wfn, pot)
   end if

   if (allocated(self%aes2)) then
      call self%aes2%get_potential_gradient(mol, cache, wfn, pot)
   end if

   if (allocated(self%es3)) then
      call self%es3%get_potential_gradient(mol, cache, wfn, pot)
   end if
end subroutine get_potential_gradient


subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the electrostatic container
   class(tb_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular gradient of the repulsion energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the repulsion energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   if (allocated(self%es2)) then
      call self%es2%get_gradient(mol, cache, wfn, gradient, sigma)
   end if

   if (allocated(self%aes2)) then
      call self%aes2%get_gradient(mol, cache, wfn, gradient, sigma)
   end if

   if (allocated(self%es3)) then
      call self%es3%get_gradient(mol, cache, wfn, gradient, sigma)
   end if
end subroutine get_gradient


pure function variable_info(self) result(info)
   use tblite_scf_info, only : scf_info, max
   !> Instance of the electrostatic container
   class(tb_coulomb), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info()

   if (allocated(self%es2)) then
      info = max(info, self%es2%variable_info())
   end if

   if (allocated(self%aes2)) then
      info = max(info, self%aes2%variable_info())
   end if

   if (allocated(self%es3)) then
      info = max(info, self%es3%variable_info())
   end if

end function variable_info


!> Information on container
pure function info(self, verbosity, indent) result(str)
   !> Instance of the interaction container
   class(tb_coulomb), intent(in) :: self
   !> Verbosity level
   integer, intent(in) :: verbosity
   !> Indentation level
   character(len=*), intent(in) :: indent
   !> Information on the container
   character(len=:), allocatable :: str

   character(len=*), parameter :: nl = new_line('a'), marker = " * "

   if (allocated(self%label)) then
      str = self%label
   else
      str = "Coulomb electrostatics"
   end if

   if (allocated(self%es2)) then
      str = str // nl // indent // marker // &
         & self%es2%info(verbosity, indent//marker)
   end if

   if (allocated(self%aes2)) then
      str = str // nl // indent // marker // &
         & self%aes2%info(verbosity, indent//marker)
   end if

   if (allocated(self%es3)) then
      str = str // nl // indent // marker // &
         & self%es3%info(verbosity, indent//marker)
   end if
end function info


end module tblite_xtb_coulomb
