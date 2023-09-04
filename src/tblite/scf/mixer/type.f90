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

!> @file tblite/scf/mixer/type.f90
!> Base class for electronic mixing

!> Provides base class for electronic mixing routines
module tblite_scf_mixer_type
   use mctc_env, only : error_type, wp
   implicit none
   private


   !> Abstract base class for electronic mixing
   type, public, abstract :: mixer_type
   contains
      !> Apply mixing to the density
      procedure(next), deferred :: next
      !> Set new density
      generic :: set => set_1d, set_2d, set_3d
      !> Set new density from 1D array
      procedure(set_1d), deferred :: set_1d
      !> Set new density from 2D array
      procedure :: set_2d
      !> Set new density from 3D array
      procedure :: set_3d
      !> Set difference between new and old density
      generic :: diff => diff_1d, diff_2d, diff_3d
      !> Set difference between new and old density from 1D array
      procedure(diff_1d), deferred :: diff_1d
      !> Set difference between new and old density from 2D array
      procedure :: diff_2d
      !> Set difference between new and old density from 3D array
      procedure :: diff_3d
      !> Get density
      generic :: get => get_1d, get_2d, get_3d
      !> Get density as 1D array
      procedure(get_1d), deferred :: get_1d
      !> Get density as 2D array
      procedure :: get_2d
      !> Get density as 3D array
      procedure :: get_3d
      !> Get error metric from mixing
      procedure(get_error), deferred :: get_error
   end type mixer_type

   abstract interface
      !> Set new density from 1D array
      subroutine set_1d(self, qvec)
         import :: mixer_type, wp
         !> Instance of the electronic mixer
         class(mixer_type), intent(inout) :: self
         !> Density vector
         real(wp), intent(in) :: qvec(:)
      end subroutine set_1d

      !> Set difference between new and old density from 1D array
      subroutine diff_1d(self, qvec)
         import :: mixer_type, wp
         !> Instance of the electronic mixer
         class(mixer_type), intent(inout) :: self
         !> Density vector
         real(wp), intent(in) :: qvec(:)
      end subroutine diff_1d

      !> Get density as 1D array
      subroutine get_1d(self, qvec)
         import :: mixer_type, wp
         !> Instance of the electronic mixer
         class(mixer_type), intent(inout) :: self
         !> Density vector
         real(wp), intent(out) :: qvec(:)
      end subroutine get_1d

      !> Apply mixing to the density
      subroutine next(self, error)
         import :: mixer_type, error_type
         !> Instance of the electronic mixer
         class(mixer_type), intent(inout) :: self
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine next

      !> Get error metric from mixing
      pure function get_error(self) result(error)
         import :: mixer_type, wp
         !> Instance of the electronic mixer
         class(mixer_type), intent(in) :: self
         !> Error metric
         real(wp) :: error
      end function get_error
   end interface

contains

!> Set new density from 2D array
subroutine set_2d(self, qvec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(inout) :: self
   !> Density vector
   real(wp), contiguous, intent(in), target :: qvec(:, :)
   real(wp), pointer :: qptr(:)
   qptr(1:size(qvec)) => qvec
   call self%set(qptr)
end subroutine set_2d

!> Set new density from 3D array
subroutine set_3d(self, qvec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(inout) :: self
   !> Density vector
   real(wp), contiguous, intent(in), target :: qvec(:, :, :)
   real(wp), pointer :: qptr(:)
   qptr(1:size(qvec)) => qvec
   call self%set(qptr)
end subroutine set_3d

!> Set difference between new and old density from 2D array
subroutine diff_2d(self, qvec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(inout) :: self
   !> Density vector
   real(wp), contiguous, intent(in), target :: qvec(:, :)
   real(wp), pointer :: qptr(:)
   qptr(1:size(qvec)) => qvec
   call self%diff(qptr)
end subroutine diff_2d

!> Set difference between new and old density from 3D array
subroutine diff_3d(self, qvec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(inout) :: self
   !> Density vector
   real(wp), contiguous, intent(in), target :: qvec(:, :, :)
   real(wp), pointer :: qptr(:)
   qptr(1:size(qvec)) => qvec
   call self%diff(qptr)
end subroutine diff_3d

!> Get density as 2D array
subroutine get_2d(self, qvec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(inout) :: self
   !> Density vector
   real(wp), contiguous, intent(out), target :: qvec(:, :)
   real(wp), pointer :: qptr(:)
   qptr(1:size(qvec)) => qvec
   call self%get(qptr)
end subroutine get_2d

!> Get density as 3D array
subroutine get_3d(self, qvec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(inout) :: self
   !> Density vector
   real(wp), contiguous, intent(out), target :: qvec(:, :, :)
   real(wp), pointer :: qptr(:)
   qptr(1:size(qvec)) => qvec
   call self%get(qptr)
end subroutine get_3d

end module tblite_scf_mixer_type
