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

!> @file tblite/coulomb/cache.f90
!> Provides a cache specific for all Coulomb interactions

!> Data container for mutable data in electrostatic calculations
module tblite_coulomb_cache
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   use tblite_coulomb_ewald, only : get_alpha, get_rec_cutoff
   use tblite_cutoff, only : get_lattice_points
   use tblite_wignerseitz, only : wignerseitz_cell, new_wignerseitz_cell
   implicit none
   private

   public :: coulomb_cache


   type :: coulomb_cache
      real(wp) :: alpha
      real(wp) :: vol
      type(wignerseitz_cell) :: wsc
      real(wp), allocatable :: amat(:, :)
      real(wp), allocatable :: vvec(:)
      real(wp), allocatable :: rtrans(:, :)

      real(wp), allocatable :: cn(:)
      real(wp), allocatable :: dcndr(:, :, :)
      real(wp), allocatable :: dcndL(:, :, :)
      real(wp), allocatable :: mrad(:)
      real(wp), allocatable :: dmrdcn(:)

      real(wp), allocatable :: amat_sd(:, :, :)
      real(wp), allocatable :: amat_dd(:, :, :, :)
      real(wp), allocatable :: amat_sq(:, :, :)
   contains
      procedure :: update
   end type coulomb_cache

   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))
   real(wp), parameter :: conv = eps

contains


subroutine update(self, mol)
   !> Instance of the electrostatic container
   class(coulomb_cache), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   if (any(mol%periodic)) then
      call new_wignerseitz_cell(self%wsc, mol)
      call get_alpha(mol%lattice, self%alpha, .false.)

      self%vol = abs(matdet_3x3(mol%lattice))
      call get_rec_trans(mol%lattice, self%alpha, self%vol, conv, self%rtrans)
   end if

end subroutine update

!> Get reciprocal lattice translations
subroutine get_rec_trans(lattice, alpha, volume, conv, trans)
   !> Lattice parameters
   real(wp), intent(in) :: lattice(:, :)
   !> Parameter for Ewald summation
   real(wp), intent(in) :: alpha
   !> Cell volume
   real(wp), intent(in) :: volume
   !> Tolerance for Ewald summation
   real(wp), intent(in) :: conv
   !> Translation vectors
   real(wp), allocatable, intent(out) :: trans(:, :)

   real(wp) :: rec_lat(3, 3)

   rec_lat = 2*pi*transpose(matinv_3x3(lattice))
   call get_lattice_points([.true.], rec_lat, get_rec_cutoff(alpha, volume, conv), trans)
   trans = trans(:, 2:)

end subroutine get_rec_trans

end module tblite_coulomb_cache
