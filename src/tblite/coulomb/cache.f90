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
   use tblite_coulomb_ewald, only : get_alpha, ewald_cache
   use tblite_partition, only : pair_list
   use tblite_wignerseitz, only : wignerseitz_cell, new_wignerseitz_cell
   implicit none
   private

   public :: coulomb_cache


   type :: coulomb_cache
      type(pair_list) :: pairs
      type(ewald_cache) :: charge_ewald, multipole_ewald
      real(wp), allocatable, private :: xyz(:, :)
      real(wp), private :: lattice(3, 3)
      logical, private :: periodic(3)
      real(wp) :: alpha
      real(wp) :: alpha_multipole
      type(wignerseitz_cell) :: wsc
      !> Contiguous local blocks, packed once per geometry for SCF contractions
      real(wp), allocatable :: local_amat(:, :, :)
      !> Charge-dipole blocks (3, 2, npair); multipole on row atom (1) or neighbour (2)
      real(wp), allocatable :: local_sd(:, :, :)
      !> Dipole-dipole blocks (3, 3, npair), with row-atom components first
      real(wp), allocatable :: local_dd(:, :, :)
      !> Charge-quadrupole blocks (6, 2, npair), with the same directions as local_sd
      !> Components are xx, xy, yy, xz, yz, zz, with doubled mixed products
      real(wp), allocatable :: local_sq(:, :, :)
      real(wp), allocatable :: vvec(:)

      real(wp), allocatable :: cn(:)
      logical :: cn_derivs_valid = .false.
      real(wp), allocatable :: dcndr(:, :, :)
      real(wp), allocatable :: dcndL(:, :, :)
      real(wp), allocatable :: mrad(:)
      real(wp), allocatable :: dmrdcn(:)

   contains
      procedure :: update
   end type coulomb_cache


contains


subroutine update(self, mol)
   !> Instance of the electrostatic container
   class(coulomb_cache), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   if (any(mol%periodic)) then
      if (allocated(self%xyz)) then
         if (size(self%xyz, 2) == mol%nat) then
            if (all(self%xyz == mol%xyz) .and. all(self%lattice == mol%lattice) &
               & .and. all(self%periodic .eqv. mol%periodic)) return
         end if
      end if
      call new_wignerseitz_cell(self%wsc, mol)
      call get_alpha(mol%lattice, self%alpha, .false.)
      call get_alpha(mol%lattice, self%alpha_multipole, .true.)
      self%xyz = mol%xyz
      self%lattice = mol%lattice
      self%periodic = mol%periodic
   end if

end subroutine update

end module tblite_coulomb_cache
