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

!> @file tblite/disp/cache.f90
!> Provides a cache for dispersion interactions

!> Reusable data container for dispersion related calculations
module tblite_disp_cache
   use mctc_env, only : wp
   implicit none
   private
   public :: dispersion_cache

   type :: dispersion_cache
      real(wp), allocatable :: dispmat(:, :, :, :)
      real(wp), allocatable :: gwvec(:, :)
      real(wp), allocatable :: dgwdq(:, :)
      real(wp), allocatable :: vvec(:, :)
      real(wp), allocatable :: cn(:)
      real(wp), allocatable :: dcndr(:, :, :)
      real(wp), allocatable :: dcndL(:, :, :)
   end type dispersion_cache


end module tblite_disp_cache
