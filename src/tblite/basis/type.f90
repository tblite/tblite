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

!> @file tblite/basis/type.f90
!> Provides data types for managing basis set information

!> Gaussian type basis set data
module tblite_basis_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none
   private

   public :: new_basis, get_cutoff

   !> Maximum contraction length of basis functions.
   !> The limit is chosen as twice the maximum size returned by the STO-NG expansion
   integer, parameter :: maxg = 12

   !> Contracted Gaussian type basis function
   type, public :: cgto_type
      !> Angular momentum of this basis function
      integer :: ang = -1
      !> Contraction length of this basis function
      integer :: nprim = 0
      !> Exponent of the primitive Gaussian functions
      real(wp) :: alpha(maxg) = 0.0_wp
      !> Contraction coefficients of the primitive Gaussian functions,
      !> might contain normalization
      real(wp) :: coeff(maxg) = 0.0_wp
   end type cgto_type

   !> Collection of information regarding the basis set of a system
   type, public :: basis_type
      !> Maximum angular momentum of all basis functions,
      !> used to determine scratch size in integral calculation
      integer :: maxl = 0
      !> Number of shells in this basis set
      integer :: nsh = 0
      !> Number of spherical atomic orbitals in this basis set
      integer :: nao = 0
      !> Integral cutoff as maximum exponent of Gaussian product theoreom to consider
      real(wp) :: intcut = 0.0_wp
      !> Smallest primitive exponent in the basis set
      real(wp) :: min_alpha = huge(0.0_wp)
      !> Number of shells for each species
      integer, allocatable :: nsh_id(:)
      !> Number of shells for each atom
      integer, allocatable :: nsh_at(:)
      !> Number of spherical atomic orbitals for each shell
      integer, allocatable :: nao_sh(:)
      !> Index offset for each shell in the atomic orbital space
      integer, allocatable :: iao_sh(:)
      !> Index offset for each atom in the shell space
      integer, allocatable :: ish_at(:)
      !> Mapping from spherical atomic orbitals to the respective atom
      integer, allocatable :: ao2at(:)
      !> Mapping from spherical atomic orbitals to the respective shell
      integer, allocatable :: ao2sh(:)
      !> Mapping from shells to the respective atom
      integer, allocatable :: sh2at(:)
      !> Contracted Gaussian basis functions forming the basis set
      type(cgto_type), allocatable :: cgto(:, :)
   end type basis_type

   !> Get optimal real space cutoff for integral evaluation
   interface get_cutoff
      module procedure :: get_cutoff
   end interface get_cutoff

contains

!> Create a new basis set
subroutine new_basis(self, mol, nshell, cgto, acc)
   !> Instance of the basis set data
   type(basis_type), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells per species
   integer, intent(in) :: nshell(:)
   !> Contracted Gaussian basis functions for each shell and species
   type(cgto_type), intent(in) :: cgto(:, :)
   !> Calculation accuracy
   real(wp), intent(in) :: acc

   integer :: iat, isp, ish, iao, ii
   real(wp) :: min_alpha

   self%nsh_id = nshell
   self%cgto = cgto
   self%intcut = integral_cutoff(acc)

   ! Make count of shells for each atom
   self%nsh_at = nshell(mol%id)

   ! Create mapping between atoms and shells
   self%nsh = sum(self%nsh_at)
   allocate(self%ish_at(mol%nat), self%sh2at(self%nsh))
   ii = 0
   do iat = 1, mol%nat
      self%ish_at(iat) = ii
      do ish = 1, self%nsh_at(iat)
         self%sh2at(ii+ish) = iat
      end do
      ii = ii + self%nsh_at(iat)
   end do

   ! Make count of spherical orbitals for each shell
   allocate(self%nao_sh(self%nsh))
   do iat = 1, mol%nat
      isp = mol%id(iat)
      ii = self%ish_at(iat)
      do ish = 1, self%nsh_at(iat)
         self%nao_sh(ii+ish) = 2*cgto(ish, isp)%ang + 1
      end do
   end do

   ! Create mapping between shells and spherical orbitals, also map directly back to atoms
   self%nao = sum(self%nao_sh)
   allocate(self%iao_sh(self%nsh), self%ao2sh(self%nao), self%ao2at(self%nao))
   ii = 0
   do ish = 1, self%nsh
      self%iao_sh(ish) = ii
      do iao = 1, self%nao_sh(ish)
         self%ao2sh(ii+iao) = ish
         self%ao2at(ii+iao) = self%sh2at(ish)
      end do
      ii = ii + self%nao_sh(ish)
   end do

   ii = 0
   do iat = 1, mol%nat
      isp = mol%id(iat)
      do ish = 1, nshell(isp)
         self%iao_sh(ish+self%ish_at(iat)) = ii
         ii = ii + 2*cgto(ish, isp)%ang + 1
      end do
   end do

   min_alpha = huge(acc)
   do isp = 1, size(nshell)
      do ish = 1, nshell(isp)
         self%maxl = max(self%maxl, cgto(ish, isp)%ang)
         min_alpha = min(min_alpha, minval(cgto(ish, isp)%alpha(:cgto(ish, isp)%nprim)))
      end do
   end do

   self%min_alpha = min_alpha

end subroutine new_basis

!> Determine required real space cutoff for the basis set
pure function get_cutoff(self, acc) result(cutoff)
   !> Instance of the basis set data
   type(basis_type), intent(in) :: self
   !> Accuracy for the integral cutoff
   real(wp), intent(in), optional :: acc
   !> Required realspace cutoff
   real(wp) :: cutoff

   real(wp) :: intcut
   real(wp), parameter :: max_cutoff = 40.0_wp

   if (present(acc)) then
      intcut = integral_cutoff(acc)
   else
      intcut = self%intcut
   end if
   ! ai * aj * cutoff2 / (ai + aj) == intcut
   cutoff = min(sqrt(2.0_wp*intcut/self%min_alpha), max_cutoff)

end function get_cutoff


!> Create integral cutoff from accuracy value
pure function integral_cutoff(acc) result(intcut)
   !> Accuracy for the integral cutoff
   real(wp), intent(in) :: acc
   !> Integral cutoff
   real(wp) :: intcut

   real(wp), parameter :: min_intcut = 5.0_wp, max_intcut = 25.0_wp, &
      & max_acc = 1.0e-4_wp, min_acc = 1.0e+3_wp

   intcut = clip(max_intcut - 10*log10(clip(acc, min_acc, max_acc)), min_intcut, max_intcut)
end function integral_cutoff


pure function clip(val, min_val, max_val) result(res)
   real(wp), intent(in) :: val, min_val, max_val
   real(wp) :: res
   res = min(max(val, min_val), max_val)
end function clip

end module tblite_basis_type
