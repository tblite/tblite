! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
! License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/scf/mixer/input.f90
!> Provides input types for electronic mixers
module tblite_scf_mixer_input
   use mctc_env, only : wp
   implicit none
   private

   public :: scf_version, trial_version
   public :: scf_version_to_string, trial_version_to_string

   !> Default value for self-consistent iteration mixing
   real(wp), parameter :: mixer_damping_default = 0.4_wp

   !> Default value for self-consistent iteration mixing memory, 0 uses max_iter
   integer, parameter :: mixer_memory_default = 0

   !> Default maximum number of self-consistent iterations
   integer, parameter :: max_iter_default = 250

   !> Default number of standard SCF cycles before trial strategies are enabled
   integer, parameter :: trial_start_default = 4

   type :: scf_enum
      integer :: broyden = 1
   end type scf_enum
   type(scf_enum), parameter :: scf_version = scf_enum()

   type :: trial_enum
      integer :: default = 10
      integer :: oda = 11
      integer :: mesa = 12
   end type trial_enum
   type(trial_enum), parameter :: trial_version = trial_enum()

   !> Input for selecting electronic mixer
   type, public :: mixer_input
      !> Annealing schedule
      type(anneal_input), allocatable :: anneal
      !> Maximum number of self-consistent iterations
      integer :: max_iter = max_iter_default
      !> Mixer memory
      integer :: memory = mixer_memory_default
      !> Damping parameter for mixing
      real(wp) :: damping = mixer_damping_default
      !> Outer SCF loop version to use
      integer :: scf = scf_version%broyden
      !> Trial inner SCF loop version to use
      integer :: trial = trial_version%default
      !> Number of standard SCF cycles before trial strategies are enabled
      integer :: trial_start = trial_start_default
   end type mixer_input

   type, public :: anneal_input
      !> Initial electronic temperature for SCF annealing in Eh.
      real(wp) :: initial_kt
      !> Number of SCF cycles to keep the initial electronic temperature before annealing (default: 50).
      integer :: hold = 50
      !> Number of SCF cycles used to lower the electronic temperature to the calculation temperature (default: 50).
      integer :: cycles = 50
   end type anneal_input

contains

pure function scf_version_to_string(scf) result(str)
   integer, intent(in) :: scf
   character(len=:), allocatable :: str

   select case(scf)
   case(scf_version%broyden)
      str = "broyden"
   case default
      str = "unknown"
   end select
end function scf_version_to_string

pure function trial_version_to_string(trial) result(str)
   integer, intent(in) :: trial
   character(len=:), allocatable :: str

   select case(trial)
   case(trial_version%default)
      str = "default"
   case(trial_version%oda)
      str = "oda"
   case(trial_version%mesa)
      str = "mesa"
   case default
      str = "unknown"
   end select
end function trial_version_to_string

end module tblite_scf_mixer_input