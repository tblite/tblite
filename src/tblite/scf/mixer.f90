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

!> @dir tblite/scf/mixer
!> Routines for implementing electronic mixing

!> @file tblite/scf/mixer.f90
!> Proxy module for electronic mixing routines

!> Provides an electronic mixer implementation
module tblite_scf_mixer
   use mctc_env, only : wp
   use tblite_scf_mixer_broyden, only : broyden_mixer, new_broyden
   use tblite_scf_mixer_type, only : mixer_type
   implicit none
   private

   public :: mixer_type, new_mixer

contains

subroutine new_mixer(self, memory, ndim, damp)
   class(mixer_type), allocatable, intent(out) :: self
   integer, intent(in) :: memory
   integer, intent(in) :: ndim
   real(wp), intent(in) :: damp

   block
      type(broyden_mixer), allocatable :: mixer
      allocate(mixer)
      call new_broyden(mixer, memory, ndim, damp)
      call move_alloc(mixer, self)
   end block
end subroutine new_mixer

end module tblite_scf_mixer
