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

!> @file tblite/solvation/shift.f90
!> Provides solvent and state specific shifts to the free energy.

module tblite_solvation_shift
    use mctc_env, only : wp, error_type, fatal_error
    use mctc_io, only : structure_type
    use tblite_container_cache, only : container_cache
    use tblite_scf_potential, only : potential_type
    use tblite_wavefunction_type, only : wavefunction_type
    use tblite_solvation_type, only : solvation_type
    implicit none
    private
 
    public :: shift_solvation, new_shift, shift_input

    type :: shift_input
        !> reference state
        integer :: state = 1
        !> solvent specific free energy shift
        real(wp) :: gshift = 0.0_wp
        !> molar mass of solvent
        real(wp), allocatable :: molar_mass 
        !> density of solvent
        real(wp), allocatable :: rho
        !> temperature in Kelvin
        real(wp) :: temperature = 298.15_wp
        !> Linearized Poisson-Boltzmann model for parameter selection
        logical :: alpb = .true.
        !> Method for parameter selection
        character(len=:), allocatable :: method
        !> Solvent for parameter selection
      character(len=:), allocatable :: solvent
    end type shift_input

    !> Possible reference states for the solution
    type :: TSolutionStateEnum

       !> 1 l of ideal gas and 1 l of liquid solution
       integer :: gsolv = 1

       !> 1 bar of ideal gas and 1 mol/L of liquid solution at infinite dilution
       integer :: reference = 2

       !> 1 bar of ideal gas and 1 mol/L of liquid solution
       integer :: mol1bar = 3

    end type TSolutionStateEnum

    !> Actual solvation state enumerator
    type(TSolutionStateEnum), parameter :: solutionState = TSolutionStateEnum()

    real(wp), parameter :: refDensity = 1.0e-3_wp ! kg__au/(1.0e10_wp*AA__Bohr)**3
    real(wp), parameter :: refMolecularMass = 1.0_wp ! amu__au
    real(wp), parameter :: idealGasMolVolume = 24.79_wp
    real(wp), parameter :: ambientTemperature = 298.15_wp ! * Boltzman

    type, extends(solvation_type) :: shift_solvation
       !> total shift
       real(wp) :: total_shift 
       contains   
       !> Update cache from container
       procedure :: update
       !> Return dependency on density
       !>procedure :: variable_info
       !> Get shift energy
       procedure :: get_engrad
       !> Get solvation energy
       procedure :: get_energy
       !> Get solvation potential
       procedure :: get_potential
       !> Get solvation gradient
       procedure :: get_gradient
    end type shift_solvation

    !> Provide constructor for shift
    interface shift_solvation
        module procedure :: create_shift
    end interface shift_solvation

    type :: shift_cache
         !> Total free energy shift
         real(wp) :: total_shift
    end type shift_cache

    !> Identifier for container 
    character(len=*), parameter :: label = "free energy shift for state and solvent"

contains

!> Calculate the solvent and state shift
subroutine new_shift(self, input)
   !> Instance of the solvation model
   type(shift_solvation) :: self
   !> Input for shift solvation
   type(shift_input), intent(in) :: input
   !> State shift
   real(wp) :: stateshift
      
   self%label = label
   stateshift = get_stateshift(input%state, input%temperature, input%rho, input%molar_mass)
   self%total_shift = input%gshift + stateshift 
end subroutine new_shift

!> Type constructor for shift
function create_shift(input) result(self)
   !> Instance of the solvation model
   type(shift_solvation) :: self
   !> Input for solvation shift
   type(shift_input), intent(in) :: input
   
   call new_shift(self, input)
end function create_shift

!> Update cache from container
subroutine update(self, mol, cache)
   !> Instance of the solvation model
   class(shift_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   type(shift_cache), pointer :: ptr

   call taint(cache, ptr)
   ptr%total_shift = self%total_shift
end subroutine update

!> Evaluate non-selfconsistent part of the interaction
subroutine get_engrad(self, mol, cache, energies, gradient, sigma)
   !> Instance of the solvation model
   class(shift_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Interaction energy
   real(wp), intent(inout) :: energies(:)
   !> Interaction gradient
   real(wp), contiguous, intent(inout), optional :: gradient(:, :)
   !> Interaction virial
   real(wp), contiguous, intent(inout), optional :: sigma(:, :)
   
   type(shift_cache), pointer :: ptr
   
   call view(cache, ptr)
    
   energies(:) = energies + ptr%total_shift/real(mol%nat)
end subroutine get_engrad

!> Get solvation energy
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the solvation model
   class(shift_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Solvation free energy
   real(wp), intent(inout) :: energies(:)

   type(shift_cache), pointer :: ptr

   call view(cache, ptr)

end subroutine get_energy

!> Get shift potential
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the solvation model
   class(shift_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   type(shift_cache), pointer :: ptr

   call view(cache, ptr)

end subroutine get_potential

!> Get shift gradient
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the solvation model
   class(shift_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular gradient of the solvation free energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the solvation free energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   type(shift_cache), pointer :: ptr

   call view(cache, ptr)

end subroutine get_gradient

function get_stateshift(state, temperature, density, molecularMass) &
   & result(stateshift)

   !> Reference state
   integer, intent(in) :: state

   !> Temperature of the solution
   real(wp), intent(in) :: temperature

   !> Mass density of the solvent
   real(wp), intent(in) :: density

   !> Molecular mass of the solvent
   real(wp), intent(in) :: molecularMass

   !> Resulting shift to the solvation free energy
   real(wp) :: stateshift

   !  Boltzmann constant in Eh/K
   real(wp),parameter :: boltzmann = 3.166808578545117e-06_wp

   select case(state)
   case default
      stateshift = 0.0_wp
   case(solutionState%reference)
      stateshift = temperature * boltzmann &
         & * (log(idealGasMolVolume * temperature / ambientTemperature) &
         & + log(density/refDensity * refMolecularMass/molecularMass))
   case(solutionState%mol1bar)
      stateshift = temperature * boltzmann &
         & * log(idealGasMolVolume * temperature / ambientTemperature)
   end select

end function get_stateshift

subroutine taint(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(shift_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(shift_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

subroutine view(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(shift_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(shift_cache)
      ptr => target
   end select
end subroutine view

end module tblite_solvation_shift
