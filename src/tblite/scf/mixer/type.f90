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

   public :: mixer_type

   type, public, abstract :: mixer_type
   contains
      procedure(next), deferred :: next
      generic :: set => set_1d, set_2d, set_3d
      procedure(set_1d), deferred :: set_1d
      procedure :: set_2d
      procedure :: set_3d
      generic :: diff => diff_1d, diff_2d, diff_3d
      procedure(diff_1d), deferred :: diff_1d
      procedure :: diff_2d
      procedure :: diff_3d
      generic :: get => get_1d, get_2d, get_3d
      procedure(get_1d), deferred :: get_1d
      procedure :: get_2d
      procedure :: get_3d
      procedure(get_error), deferred :: get_error
   end type mixer_type

   abstract interface
      subroutine set_1d(self, qvec)
         import :: mixer_type, wp
         class(mixer_type), intent(inout) :: self
         real(wp), intent(in) :: qvec(:)
      end subroutine set_1d

      subroutine diff_1d(self, qvec)
         import :: mixer_type, wp
         class(mixer_type), intent(inout) :: self
         real(wp), intent(in) :: qvec(:)
      end subroutine diff_1d

      subroutine get_1d(self, qvec)
         import :: mixer_type, wp
         class(mixer_type), intent(inout) :: self
         real(wp), intent(out) :: qvec(:)
      end subroutine get_1d

      subroutine next(self, error)
         import :: mixer_type, error_type
         class(mixer_type), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine next

      pure function get_error(self) result(error)
         import :: mixer_type, wp
         class(mixer_type), intent(in) :: self
         real(wp) :: error
      end function get_error
   end interface

contains

subroutine set_2d(self, qvec)
   class(mixer_type), intent(inout) :: self
   real(wp), contiguous, intent(in), target :: qvec(:, :)
   real(wp), pointer :: qptr(:)
   qptr(1:size(qvec)) => qvec
   call self%set(qptr)
end subroutine set_2d

subroutine set_3d(self, qvec)
   class(mixer_type), intent(inout) :: self
   real(wp), contiguous, intent(in), target :: qvec(:, :, :)
   real(wp), pointer :: qptr(:)
   qptr(1:size(qvec)) => qvec
   call self%set(qptr)
end subroutine set_3d

subroutine diff_2d(self, qvec)
   class(mixer_type), intent(inout) :: self
   real(wp), contiguous, intent(in), target :: qvec(:, :)
   real(wp), pointer :: qptr(:)
   qptr(1:size(qvec)) => qvec
   call self%diff(qptr)
end subroutine diff_2d

subroutine diff_3d(self, qvec)
   class(mixer_type), intent(inout) :: self
   real(wp), contiguous, intent(in), target :: qvec(:, :, :)
   real(wp), pointer :: qptr(:)
   qptr(1:size(qvec)) => qvec
   call self%diff(qptr)
end subroutine diff_3d

subroutine get_2d(self, qvec)
   class(mixer_type), intent(inout) :: self
   real(wp), contiguous, intent(out), target :: qvec(:, :)
   real(wp), pointer :: qptr(:)
   qptr(1:size(qvec)) => qvec
   call self%get(qptr)
end subroutine get_2d

subroutine get_3d(self, qvec)
   class(mixer_type), intent(inout) :: self
   real(wp), contiguous, intent(out), target :: qvec(:, :, :)
   real(wp), pointer :: qptr(:)
   qptr(1:size(qvec)) => qvec
   call self%get(qptr)
end subroutine get_3d

end module tblite_scf_mixer_type
