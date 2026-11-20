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

!> @file tblite/post_processing/xtb-ml/cache.f90
!> Provides a global cache for all xtb-ML features

!> Data container for mutable data in the xtb-ML feature calculation
module tblite_post_processing_xtbml_cache
   use mctc_env, only : wp
   implicit none
   private

   public :: xtbml_cache

   !> Cache for xTB-ML features calculation
   type :: xtbml_cache
      ! Geometry-based features
      !> Atomic coordination numbers
      real(wp), allocatable :: cn(:)

      ! Density-based features
      !> Shell populations (updown representation)
      real(wp), allocatable :: psh(:, :)
      !> Atomic charges (updown representation)
      real(wp), allocatable :: qat(:, :)
      !> Shell dipole moments (updown representation)
      real(wp), allocatable :: dpsh(:, :, :)
      !> Shell quadrupole moments (updown representation)
      real(wp), allocatable :: qpsh(:, :, :)
      !> Atomic dipole moments (updown representation)
      real(wp), allocatable :: dpsh_norm(:, :)
      !> Atomic quadrupole moments (updown representation)
      real(wp), allocatable :: qpsh_norm(:, :)
      !> Atomic dipole moment norms (updown representation)
      real(wp), allocatable :: dpat_norm(:, :)
      !> Atomic quadrupole moment norms (updown representation)
      real(wp), allocatable :: qpat_norm(:, :)

      ! Energy-based features
      !> Effective atomic HOMO-LUMO gap
      real(wp), allocatable ::  egap(:, :)
      !> Effective Fermi level per atom
      real(wp), allocatable ::  chempot(:, :)

      ! Convolution
      !> Convolution kernel
      real(wp), allocatable :: conv_kernel(:, :, :)
      !> Convolution coordination number at different CN length scales
      real(wp), allocatable :: conv_cn(:, :)

   end type xtbml_cache

end module tblite_post_processing_xtbml_cache
