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

module tblite_adjlist
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_resize, only : resize
   implicit none
   private

   public :: adjacency_list, new_adjacency_list

   !> Neighbourlist in CSR format
   type :: adjacency_list
      integer, allocatable :: inl(:)
      integer, allocatable :: nnl(:)
      integer, allocatable :: nlat(:)
      integer, allocatable :: nltr(:)
   end type adjacency_list

contains

   subroutine new_adjacency_list(self, mol, trans, cutoff)
      type(adjacency_list), intent(out) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: trans(:, :)
      real(wp), intent(in) :: cutoff

      allocate(self%inl(mol%nat), source=0)
      allocate(self%nnl(mol%nat), source=0)
      call generate(mol, trans, cutoff, self%inl, self%nnl, self%nlat, self%nltr)
   end subroutine new_adjacency_list

   subroutine generate(mol, trans, cutoff, inl, nnl, nlat, nltr)
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: trans(:, :)
      real(wp), intent(in) :: cutoff
      integer, intent(inout) :: inl(:)
      integer, intent(inout) :: nnl(:)
      integer, allocatable, intent(out) :: nlat(:)
      integer, allocatable, intent(out) :: nltr(:)

      integer :: iat, jat, itr, img
      real(wp) :: r2, vec(3), cutoff2

      img = 0
      cutoff2 = cutoff**2

      call resize(nlat, 10*mol%nat)
      call resize(nltr, 10*mol%nat)

      do iat = 1, mol%nat
         inl(iat) = img
         do jat = 1, iat
            do itr = 1, size(trans, 2)
               vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
               r2 = sum(vec**2)
               if (r2 < epsilon(cutoff2) .or. r2 > cutoff2) cycle
               img = img + 1
               if (size(nlat) < img) call resize(nlat)
               if (size(nltr) < img) call resize(nltr)
               nlat(img) = jat
               nltr(img) = itr
            end do
         end do
         nnl(iat) = img - inl(iat)
      end do

      call resize(nlat, img)
      call resize(nltr, img)

   end subroutine generate

end module tblite_adjlist
