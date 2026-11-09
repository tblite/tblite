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

!> @file tblite/coulomb/charge/effective.f90
!> Provides an effective Coulomb operator for isotropic electrostatic interactions

!> Isotropic second-order electrostatics using an effective Coulomb operator
module tblite_coulomb_charge_effective
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_math, only : outerprod
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_charge_type, only : coulomb_charge_type
   use tblite_coulomb_ewald, only : ewald_cache, get_amat_rec_3d, get_damat_rec_3d
   use tblite_partition, only : pair_list
   use tblite_wignerseitz, only : wignerseitz_cell, get_wignerseitz_weights
   implicit none
   private

   public :: new_effective_coulomb
   public :: average_interface, harmonic_average, arithmetic_average, geometric_average


   !> Effective, Klopman-Ohno-type, second-order electrostatics
   type, public, extends(coulomb_charge_type) :: effective_coulomb
      !> Hubbard parameter for each shell and species
      real(wp), allocatable :: hubbard(:, :, :, :)
      !> Exponent of Coulomb kernel
      real(wp) :: gexp
   contains
      !> Evaluate the Coulomb matrix blocks of the locally owned atom pairs
      procedure :: get_coulomb_matrix
      !> Evaluate the Coulomb matrix derivatives of one atom pair
      procedure :: get_pair_derivs
   end type effective_coulomb


   !> Geometry shared by all shell combinations of one periodic atom pair.
   type :: ko_distances
      !> Squared distances to the periodic images
      real(wp), allocatable :: r2(:)
      !> Inverse distances, zero for excluded zero-distance images
      real(wp), allocatable :: invr(:)
      !> Cubes of the inverse distances
      real(wp), allocatable :: invr3(:)
   contains
      !> Release the pair geometry
      procedure :: delete => delete_ko_distances
   end type ko_distances

   abstract interface
      !> Average Hubbard parameter for two shells
      pure function average_interface(gi, gj) result(gij)
         import :: wp
         !> Hubbard parameter of shell i
         real(wp), intent(in) :: gi
         !> Hubbard parameter of shell j
         real(wp), intent(in) :: gj
         !> Averaged Hubbard parameter
         real(wp) :: gij
      end function average_interface
   end interface

   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))
   real(wp), parameter :: conv = epsilon(0.0_wp)
   real(wp), parameter :: euler_gamma = 0.57721566490153286061_wp
   real(wp), parameter :: s3_constant = 0.5_wp*(log(pi) - 2.0_wp + euler_gamma &
      & + 2.0_wp*log(2.0_wp))
   character(len=*), parameter :: label = "isotropic Klopman-Ohno-Mataga-Nishimoto electrostatics"

contains

!> Construct new effective electrostatic interaction container
subroutine new_effective_coulomb(self, mol, gexp, hubbard, average, nshell)
   !> Instance of the electrostatic container
   type(effective_coulomb), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Exponent of Coulomb kernel
   real(wp), intent(in) :: gexp
   !> Averaging function for Hubbard parameter of a shell-pair
   procedure(average_interface) :: average
   !> Hubbard parameter for all shells and species
   real(wp), intent(in) :: hubbard(:, :)
   !> Number of shells for each species
   integer, intent(in), optional :: nshell(:)

   integer :: mshell
   integer :: isp, jsp, ish, jsh, ind, iat

   self%label = label

   self%shell_resolved = present(nshell)
   if (present(nshell)) then
      mshell = maxval(nshell)
      self%nshell = nshell(mol%id)
   else
      mshell = 1
      self%nshell = spread(1, 1, mol%nat)
   end if
   allocate(self%offset(mol%nat))
   ind = 0
   do iat = 1, mol%nat
      self%offset(iat) = ind
      ind = ind + self%nshell(iat)
   end do

   self%gexp = gexp

   if (present(nshell)) then
      allocate(self%hubbard(mshell, mshell, mol%nid, mol%nid))
      do isp = 1, mol%nid
         do jsp = 1, mol%nid
            self%hubbard(:, :, jsp, isp) = 0.0_wp
            do ish = 1, nshell(isp)
               do jsh = 1, nshell(jsp)
                  self%hubbard(jsh, ish, jsp, isp) = &
                     & average(hubbard(ish, isp), hubbard(jsh, jsp))
               end do
            end do
         end do
      end do
   else
      allocate(self%hubbard(1, 1, mol%nid, mol%nid))
      do isp = 1, mol%nid
         do jsp = 1, mol%nid
            self%hubbard(1, 1, jsp, isp) = average(hubbard(1, isp), hubbard(1, jsp))
         end do
      end do
   end if

end subroutine new_effective_coulomb


!> Harmonic averaging functions for hardnesses in GFN1-xTB
pure function harmonic_average(gi, gj) result(gij)
   !> Hubbard parameter of shell i
   real(wp), intent(in) :: gi
   !> Hubbard parameter of shell j
   real(wp), intent(in) :: gj
   !> Averaged Hubbard parameter
   real(wp) :: gij

   gij = 2.0_wp/(1.0_wp/gi+1.0_wp/gj)

end function harmonic_average


!> Arithmetic averaging functions for hardnesses in GFN2-xTB
pure function arithmetic_average(gi, gj) result(gij)
   !> Hubbard parameter of shell i
   real(wp), intent(in) :: gi
   !> Hubbard parameter of shell j
   real(wp), intent(in) :: gj
   !> Averaged Hubbard parameter
   real(wp) :: gij

   gij = 0.5_wp*(gi+gj)

end function arithmetic_average


!> Geometric averaging functions for hardnesses
pure function geometric_average(gi, gj) result(gij)
   !> Hubbard parameter of shell i
   real(wp), intent(in) :: gi
   !> Hubbard parameter of shell j
   real(wp), intent(in) :: gj
   !> Averaged Hubbard parameter
   real(wp) :: gij

   gij = sqrt(gi*gj)

end function geometric_average


!> Evaluate the Coulomb matrix blocks of the locally owned atom pairs
subroutine get_coulomb_matrix(self, mol, cache)
   !> Instance of the electrostatic container
   class(effective_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container, local pair list and Coulomb blocks to overwrite
   type(coulomb_cache), intent(inout) :: cache

   cache%local_amat(:, :, :) = 0.0_wp
   if (any(mol%periodic)) then
      call cache%charge_ewald%update(mol%lattice, cache%alpha, conv, abs(self%gexp-2.0_wp) < eps)
      call get_amat_3d(mol, self%nshell, self%hubbard, self%gexp, &
         & cache%wsc, cache%alpha, cache%charge_ewald, cache%local_amat, cache%pairs)
   else
      call get_amat_0d(mol, self%nshell, self%hubbard, self%gexp, cache%local_amat, cache%pairs)
   end if
end subroutine get_coulomb_matrix


!> Evaluate Coulomb matrix for finite systems
subroutine get_amat_0d(mol, nshell, hubbard, gexp, amat, pairs)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each atom
   integer, intent(in) :: nshell(:)
   !> Hubbard parameter parameter for each shell
   real(wp), intent(in) :: hubbard(:, :, :, :)
   !> Exponent of Coulomb kernel
   real(wp), intent(in) :: gexp
   !> Coulomb blocks in local pair order
   real(wp), intent(inout) :: amat(:, :, :)
   !> Locally owned atom pairs in row order
   type(pair_list), intent(in) :: pairs

   integer :: iat, jat, izp, jzp, ipair, jpair, ish, jsh
   real(wp) :: vec(3), r1, r1g, gam, tmp

   ! Cyclic rows balance the triangular pair loop without changing MPI ownership.
   !$omp parallel do default(none) schedule(static, 1) &
   !$omp shared(amat, mol, nshell, hubbard, gexp, pairs) &
   !$omp private(iat, izp, ipair, ish, jat, jzp, jpair, jsh, gam, vec, r1, r1g, tmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do ipair = pairs%offset(iat-1)+1, pairs%offset(iat)
         jat = pairs%neighbour(ipair)
         if (jat > iat) cycle
         jpair = pairs%find(jat, iat)
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         r1g = r1**gexp
         do ish = 1, nshell(iat)
            do jsh = 1, merge(ish, nshell(jat), iat == jat)
               gam = hubbard(jsh, ish, jzp, izp)
               tmp = gam
               if (iat /= jat) tmp = 1.0_wp/(r1g + gam**(-gexp))**(1.0_wp/gexp)
               if (iat /= jat .or. ish /= jsh) amat(jsh, ish, jpair) = amat(jsh, ish, jpair) + tmp
               amat(ish, jsh, ipair) = amat(ish, jsh, ipair) + tmp
            end do
         end do
      end do
   end do

end subroutine get_amat_0d


!> Evaluate the coulomb matrix for 3D systems
subroutine get_amat_3d(mol, nshell, hubbard, gexp, wsc, alpha, ewald, amat, pairs)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached Ewald lattice data
   type(ewald_cache), intent(in) :: ewald
   !> Number of shells per atom
   integer, intent(in) :: nshell(:)
   !> Hubbard parameter of the shells
   real(wp), intent(in) :: hubbard(:, :, :, :)
   !> Exponent of the interaction kernel
   real(wp), intent(in) :: gexp
   !> Wigner-Seitz cell
   type(wignerseitz_cell), intent(in) :: wsc
   !> Convergence factor for Ewald sum
   real(wp), intent(in) :: alpha
   !> Coulomb blocks in local pair order
   real(wp), intent(inout) :: amat(:, :, :)
   !> Locally owned atom pairs in row order
   type(pair_list), intent(in) :: pairs

   integer :: iat, jat, izp, jzp, img, ipair, jpair, ish, jsh
   real(wp) :: vec(3), gam, dtmp, rtmp, stmp, aval
   real(wp) :: weight(size(wsc%tridx, 1))

   if (abs(gexp - 2.0_wp) < eps) then
      call get_amat_ko_3d(mol, nshell, hubbard, alpha, ewald, amat, pairs)
      return
   end if

   !$omp parallel do default(none) schedule(static, 1) shared(amat) &
   !$omp shared(mol, nshell, hubbard, gexp, wsc, ewald, alpha, pairs) &
   !$omp private(iat, izp, jat, jzp, ipair, jpair, ish, jsh, gam, weight, vec, dtmp, rtmp, stmp, aval)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do ipair = pairs%offset(iat-1)+1, pairs%offset(iat)
         jat = pairs%neighbour(ipair)
         if (jat > iat) cycle
         jpair = pairs%find(jat, iat)
         jzp = mol%id(jat)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         call get_wignerseitz_weights(wsc, jat, iat, vec, weight)
         call get_amat_dir_3d(vec, alpha, ewald%dtrans, dtmp)
         call get_amat_rec_3d(vec, ewald, rtmp)
         if (iat == jat) rtmp = rtmp - 2 * alpha / sqrtpi
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
            do ish = 1, nshell(iat)
               do jsh = 1, merge(ish, nshell(jat), iat == jat)
                  gam = hubbard(jsh, ish, jzp, izp)
                  call get_amat_wsc_3d(vec, gam, gexp, stmp)
                  aval = dtmp + rtmp + stmp
                  if (iat == jat) aval = aval + gam
                  aval = aval * weight(img)
                  if (iat /= jat .or. ish /= jsh) amat(jsh, ish, jpair) = amat(jsh, ish, jpair) + aval
                  amat(ish, jsh, ipair) = amat(ish, jsh, ipair) + aval
               end do
            end do
         end do
      end do

   end do

end subroutine get_amat_3d


!> Evaluate the periodic Klopman-Ohno matrix with the generalized Ewald sum
subroutine get_amat_ko_3d(mol, nshell, hubbard, alpha, ewald, amat, pairs)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached Ewald lattice data
   type(ewald_cache), intent(in) :: ewald
   !> Number of shells for each atom
   integer, intent(in) :: nshell(:)
   !> Hubbard parameter for each shell and species
   real(wp), intent(in) :: hubbard(:, :, :, :)
   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha
   !> Coulomb matrix elements
   real(wp), intent(inout) :: amat(:, :, :)
   !> Locally owned atom pairs in row order
   type(pair_list), intent(in) :: pairs

   integer :: iat, jat, izp, jzp, ipair, jpair, ish, jsh
   real(wp) :: vec(3), gam, val, s1, s3, sr
   type(ko_distances) :: distances

   !$omp parallel default(none) shared(amat) &
   !$omp shared(mol, nshell, hubbard, alpha, ewald, pairs) &
   !$omp private(iat, izp, ipair, ish, jat, jzp, jpair, jsh, gam, vec, val, s1, s3, sr, distances)
   !$omp do schedule(static, 1)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do ipair = pairs%offset(iat-1)+1, pairs%offset(iat)
         jat = pairs%neighbour(ipair)
         if (jat > iat) cycle
         jpair = pairs%find(jat, iat)
         jzp = mol%id(jat)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         call get_ko_distances(vec, ewald%strans, distances)
         call get_amat_ko_ewald_3d(vec, alpha, ewald, distances, iat == jat, s1, s3)
         do ish = 1, nshell(iat)
            do jsh = 1, merge(ish, nshell(jat), iat == jat)
               gam = hubbard(jsh, ish, jzp, izp)
               call get_amat_ko_short_3d(gam, distances, sr)
               val = s1 + sr - 0.5_wp*s3/(gam*gam)
               if (iat == jat) val = val + gam
               if (iat /= jat .or. ish /= jsh) amat(jsh, ish, jpair) = amat(jsh, ish, jpair) + val
               amat(ish, jsh, ipair) = amat(ish, jsh, ipair) + val
            end do
         end do
      end do

   end do
   !$omp end do
   call distances%delete()
   !$omp end parallel
end subroutine get_amat_ko_3d


!> Shell-independent generalized Ewald sums for one atom pair
subroutine get_amat_ko_ewald_3d(rij, alpha, ewald, distances, onsite, s1, s3)
   !> Displacement between the atom pair
   real(wp), intent(in) :: rij(3)
   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha
   !> Pair geometry for the short-range translations
   type(ko_distances), intent(in) :: distances
   !> Cached Ewald lattice data
   type(ewald_cache), intent(in) :: ewald
   !> Whether both indices belong to the same atom
   logical, intent(in) :: onsite
   !> Generalized Ewald sum for 1/r
   real(wp), intent(out) :: s1
   !> Generalized Ewald sum for 1/r^3
   real(wp), intent(out) :: s3

   integer :: itr
   real(wp) :: r1, r2, k, tmp, phase

   call get_amat_dir_3d(rij, alpha, ewald%dtrans, s1)
   if (onsite) s1 = s1 - 2.0_wp*alpha/sqrtpi

   s3 = 0.0_wp
   do itr = 1, size(distances%r2)
      if (distances%invr(itr) == 0.0_wp) cycle
      r2 = distances%r2(itr)
      r1 = sqrt(r2)
      s3 = s3 + erfc(alpha*r1)/(r2*r1) &
         & + 2.0_wp*alpha*exp(-alpha*alpha*r2)/(sqrtpi*r2)
   end do

   tmp = 0.0_wp
   do itr = 1, size(ewald%weight)
      phase = cos(dot_product(rij, ewald%vec(:, itr)))
      tmp = tmp + ewald%weight(itr)*phase
      s3 = s3 + ewald%weight3(itr)*phase
   end do
   s1 = s1 + tmp
   k = alpha/sqrtpi
   s3 = s3 + 4.0_wp*pi/ewald%vol*(log(k) + s3_constant)
   if (onsite) s3 = s3 - 4.0_wp*pi/3.0_wp*k**3
end subroutine get_amat_ko_ewald_3d


!> Distances from an atom pair to all real-space translations
subroutine get_ko_distances(rij, trans, distances)
   !> Displacement between the atom pair
   real(wp), intent(in) :: rij(3)
   !> Lattice translation vectors
   real(wp), intent(in) :: trans(:, :)
   !> Pair geometry for the short-range translations
   type(ko_distances), intent(inout) :: distances
   integer :: itr, ntrans
   real(wp) :: vec(3), r1, g1

   ntrans = size(trans, 2)
   if (.not.allocated(distances%r2)) then
      allocate(distances%r2(ntrans), distances%invr(ntrans), distances%invr3(ntrans))
   end if
   do itr = 1, ntrans
      vec = rij + trans(:, itr)
      distances%r2(itr) = dot_product(vec, vec)
      r1 = sqrt(distances%r2(itr))
      g1 = 0.0_wp
      if (r1 >= eps) g1 = 1.0_wp/r1
      distances%invr(itr) = g1
      distances%invr3(itr) = g1*g1*g1
   end do
end subroutine get_ko_distances


!> Release each component independently
subroutine delete_ko_distances(self)
   !> Pair geometry to release
   class(ko_distances), intent(inout) :: self

   if (allocated(self%r2)) deallocate(self%r2)
   if (allocated(self%invr)) deallocate(self%invr)
   if (allocated(self%invr3)) deallocate(self%invr3)
end subroutine delete_ko_distances


!> Hardness-dependent short-range Klopman-Ohno residual
subroutine get_amat_ko_short_3d(gam, distances, val)
   !> Averaged Hubbard parameter
   real(wp), intent(in) :: gam
   !> Pair geometry for the short-range translations
   type(ko_distances), intent(in) :: distances
   !> Short-range interaction contribution
   real(wp), intent(out) :: val
   integer :: itr
   real(wp) :: gam2

   gam2 = 1.0_wp/(gam*gam)
   val = 0.0_wp
   do itr = 1, size(distances%r2)
      if (distances%invr(itr) == 0.0_wp) cycle
      val = val + 1.0_wp/sqrt(distances%r2(itr) + gam2) - distances%invr(itr) &
         & + 0.5_wp*gam2*distances%invr3(itr)
   end do
end subroutine get_amat_ko_short_3d


!> Calculate direct space Ewald contribution for a pair under 3D periodic boundary conditions
subroutine get_amat_dir_3d(rij, alp, trans, amat)
   !> Distance between pair
   real(wp), intent(in) :: rij(3)
   !> Convergence factor
   real(wp), intent(in) :: alp
   !> Translation vectors to consider
   real(wp), intent(in) :: trans(:, :)
   !> Interaction matrix element
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: vec(3), r1

   amat = 0.0_wp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      amat = amat + erfc(alp*r1)/r1
   end do

end subroutine get_amat_dir_3d


!> Calculate Wigner-Seitz averaged short range correction for a pair
subroutine get_amat_wsc_3d(rij, gam, gexp, amat)
   !> Distance between pair
   real(wp), intent(in) :: rij(3)
   !> Hubbard parameter
   real(wp), intent(in) :: gam
   !> Exponent for interaction kernel
   real(wp), intent(in) :: gexp
   !> Interaction matrix element
   real(wp), intent(out) :: amat

   real(wp) :: r1

   amat = 0.0_wp

   r1 = norm2(rij)
   if (r1 < eps) return
   amat = 1.0_wp/(r1**gexp + gam**(-gexp))**(1.0_wp/gexp) - 1.0_wp/r1

end subroutine get_amat_wsc_3d


!> Evaluate the Coulomb matrix derivatives of one atom pair. For iat == jat only the
!> strain derivatives of the shell pairs jsh <= ish are used, otherwise all shell pairs.
subroutine get_pair_derivs(self, mol, cache, iat, jat, dG, dS)
   !> Instance of the electrostatic container
   class(effective_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container, updated for this geometry
   type(coulomb_cache), intent(in) :: cache
   !> First atom
   integer, intent(in) :: iat
   !> Second atom
   integer, intent(in) :: jat
   !> Derivatives with respect to the position of the first atom, indexed by (:, jsh, ish)
   real(wp), intent(out) :: dG(:, :, :)
   !> Strain derivatives, indexed by (:, :, jsh, ish), only evaluated for periodic systems
   real(wp), intent(out) :: dS(:, :, :, :)

   if (any(mol%periodic)) then
      if (abs(self%gexp-2.0_wp) < eps) then
         call get_damat_pair_ko_3d(mol, self%nshell, self%hubbard, cache%alpha, &
            & cache%charge_ewald, iat, jat, dG, dS)
      else
         call get_damat_pair_3d(mol, self%nshell, self%hubbard, self%gexp, cache%wsc, &
            & cache%alpha, cache%charge_ewald, iat, jat, dG, dS)
      end if
   else
      call get_damat_pair_0d(mol, self%nshell, self%hubbard, self%gexp, iat, jat, dG)
   end if
end subroutine get_pair_derivs


!> Coulomb matrix derivatives of one atom pair in a finite system
subroutine get_damat_pair_0d(mol, nshell, hubbard, gexp, iat, jat, dG)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each atom
   integer, intent(in) :: nshell(:)
   !> Hubbard parameter for each shell and species
   real(wp), intent(in) :: hubbard(:, :, :, :)
   !> Exponent of Coulomb kernel
   real(wp), intent(in) :: gexp
   !> First atom
   integer, intent(in) :: iat
   !> Second atom
   integer, intent(in) :: jat
   !> Derivatives with respect to the position of the first atom
   real(wp), intent(out) :: dG(:, :, :)

   integer :: izp, jzp, ish, jsh
   real(wp) :: vec(3), r1, r1g, rfac, gam, dtmp

   izp = mol%id(iat)
   jzp = mol%id(jat)
   vec = mol%xyz(:, iat) - mol%xyz(:, jat)
   r1 = norm2(vec)
   r1g = r1**gexp
   rfac = -r1**(gexp-2.0_wp)
   do ish = 1, nshell(iat)
      do jsh = 1, nshell(jat)
         gam = hubbard(jsh, ish, jzp, izp)
         dtmp = 1.0_wp/(r1g + gam**(-gexp))
         dtmp = rfac*dtmp*dtmp**(1.0_wp/gexp)
         dG(:, jsh, ish) = dtmp*vec
      end do
   end do
end subroutine get_damat_pair_0d


!> Coulomb matrix derivatives of one atom pair in a 3D system using Wigner-Seitz images
subroutine get_damat_pair_3d(mol, nshell, hubbard, gexp, wsc, alpha, ewald, iat, jat, dG, dS)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each atom
   integer, intent(in) :: nshell(:)
   !> Hubbard parameter for each shell and species
   real(wp), intent(in) :: hubbard(:, :, :, :)
   !> Exponent of Coulomb kernel
   real(wp), intent(in) :: gexp
   !> Wigner-Seitz cell
   type(wignerseitz_cell), intent(in) :: wsc
   !> Convergence factor for Ewald sum
   real(wp), intent(in) :: alpha
   !> Cached Ewald lattice data
   type(ewald_cache), intent(in) :: ewald
   !> First atom
   integer, intent(in) :: iat
   !> Second atom
   integer, intent(in) :: jat
   !> Derivatives with respect to the position of the first atom
   real(wp), intent(out) :: dG(:, :, :)
   !> Strain derivatives
   real(wp), intent(out) :: dS(:, :, :, :)

   integer :: izp, jzp, ish, jsh, img
   logical :: need_weight_energy
   real(wp) :: vec(3), gam, stmp
   real(wp) :: dGd(3), dSd(3, 3), dGr(3), dSr(3, 3), dGsr(3), dSsr(3, 3)
   real(wp) :: weight(size(wsc%tridx, 1))
   real(wp) :: dwdr(3, size(wsc%tridx, 1)), dwdL(3, 3, size(wsc%tridx, 1))

   izp = mol%id(iat)
   jzp = mol%id(jat)
   vec = mol%xyz(:, iat) - mol%xyz(:, jat)
   call get_wignerseitz_weights(wsc, jat, iat, vec, weight, dwdr, dwdL)
   need_weight_energy = any(dwdr(:, :wsc%nimg(jat, iat)) /= 0.0_wp) &
      & .or. any(dwdL(:, :, :wsc%nimg(jat, iat)) /= 0.0_wp)
   call get_damat_dir_3d(vec, alpha, ewald%dtrans, dGd, dSd)
   call get_damat_rec_3d(vec, ewald, dGr, dSr)
   dG = 0.0_wp
   dS = 0.0_wp
   do img = 1, wsc%nimg(jat, iat)
      vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
      do ish = 1, nshell(iat)
         do jsh = 1, merge(ish, nshell(jat), iat == jat)
            gam = hubbard(jsh, ish, jzp, izp)
            stmp = 0.0_wp
            if (need_weight_energy) call get_amat_wsc_3d(vec, gam, gexp, stmp)
            call get_damat_wsc_3d(vec, gam, gexp, dGsr, dSsr)
            dG(:, jsh, ish) = dG(:, jsh, ish) + (dGd + dGr + dGsr)*weight(img) &
               & + stmp*dwdr(:, img)
            dS(:, :, jsh, ish) = dS(:, :, jsh, ish) + (dSd + dSr + dSsr)*weight(img) &
               & + stmp*dwdL(:, :, img)
         end do
      end do
   end do
end subroutine get_damat_pair_3d


!> Coulomb matrix derivatives of one atom pair with the generalized Klopman-Ohno Ewald sum
subroutine get_damat_pair_ko_3d(mol, nshell, hubbard, alpha, ewald, iat, jat, dG, dS)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each atom
   integer, intent(in) :: nshell(:)
   !> Hubbard parameter for each shell and species
   real(wp), intent(in) :: hubbard(:, :, :, :)
   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha
   !> Cached Ewald lattice data
   type(ewald_cache), intent(in) :: ewald
   !> First atom
   integer, intent(in) :: iat
   !> Second atom
   integer, intent(in) :: jat
   !> Derivatives with respect to the position of the first atom
   real(wp), intent(out) :: dG(:, :, :)
   !> Strain derivatives
   real(wp), intent(out) :: dS(:, :, :, :)

   integer :: izp, jzp, ish, jsh, itr
   real(wp) :: rij(3), vec(3), vv(3, 3), r1, r2, g1, invr3, invr5, gam, gam2, tmp, dtmp
   real(wp) :: dG1(3), dS1(3, 3), dG3(3), dS3(3, 3), dGk1(3), dSk1(3, 3), dGk3(3), dSk3(3, 3)

   izp = mol%id(iat)
   jzp = mol%id(jat)
   rij = mol%xyz(:, iat) - mol%xyz(:, jat)
   dG3 = 0.0_wp
   dS3 = 0.0_wp
   dG = 0.0_wp
   dS = 0.0_wp
   ! The real-space 1/r^3 sum and the hardness-dependent residual share the short-range images
   do itr = 1, size(ewald%strans, 2)
      vec = rij + ewald%strans(:, itr)
      r2 = dot_product(vec, vec)
      r1 = sqrt(r2)
      if (r1 < eps) cycle
      g1 = 1.0_wp/r1
      invr3 = g1*g1*g1
      invr5 = invr3*g1*g1
      vv = outerprod(vec, vec)
      dtmp = -3.0_wp*erfc(alpha*r1)/(r2*r2*r1) &
         & - 6.0_wp*alpha*exp(-alpha*alpha*r2)/(sqrtpi*r2*r2) &
         & - 4.0_wp*alpha**3*exp(-alpha*alpha*r2)/(sqrtpi*r2)
      dG3 = dG3 + dtmp*vec
      dS3 = dS3 + dtmp*vv
      do ish = 1, nshell(iat)
         do jsh = 1, merge(ish, nshell(jat), iat == jat)
            gam = hubbard(jsh, ish, jzp, izp)
            gam2 = 1.0_wp/(gam*gam)
            tmp = 1.0_wp/(r2 + gam2)
            dtmp = -tmp*sqrt(tmp) + invr3 - 1.5_wp*gam2*invr5
            dG(:, jsh, ish) = dG(:, jsh, ish) + dtmp*vec
            dS(:, :, jsh, ish) = dS(:, :, jsh, ish) + dtmp*vv
         end do
      end do
   end do

   call get_damat_dir_3d(rij, alpha, ewald%dtrans, dG1, dS1)
   call get_damat_ko_rec_3d(rij, alpha, ewald, dGk1, dSk1, dGk3, dSk3)
   dG1 = dG1 + dGk1
   dS1 = dS1 + dSk1
   dG3 = dG3 + dGk3
   dS3 = dS3 + dSk3
   do ish = 1, nshell(iat)
      do jsh = 1, merge(ish, nshell(jat), iat == jat)
         gam = hubbard(jsh, ish, jzp, izp)
         dG(:, jsh, ish) = dG1 + dG(:, jsh, ish) - 0.5_wp*dG3/(gam*gam)
         dS(:, :, jsh, ish) = dS1 + dS(:, :, jsh, ish) - 0.5_wp*dS3/(gam*gam)
      end do
   end do
end subroutine get_damat_pair_ko_3d


!> Reciprocal-space derivatives of the generalized Ewald sums for 1/r and 1/r^3
subroutine get_damat_ko_rec_3d(rij, alpha, ewald, dg1, ds1, dg3, ds3)
   !> Displacement between the atom pair
   real(wp), intent(in) :: rij(3)
   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha
   !> Cached Ewald lattice data
   type(ewald_cache), intent(in) :: ewald
   !> Cartesian derivatives of the 1/r sum
   real(wp), intent(out) :: dg1(3)
   !> Strain derivatives of the 1/r sum
   real(wp), intent(out) :: ds1(3, 3)
   !> Cartesian derivatives of the 1/r^3 sum
   real(wp), intent(out) :: dg3(3)
   !> Strain derivatives of the 1/r^3 sum
   real(wp), intent(out) :: ds3(3, 3)

   integer :: itr
   real(wp) :: vec(3), phase, sink, cosk, k0
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))

   dg1 = 0.0_wp
   ds1 = 0.0_wp
   dg3 = 0.0_wp
   ds3 = 0.0_wp
   do itr = 1, size(ewald%weight)
      vec = ewald%vec(:, itr)
      phase = dot_product(rij, vec)
      sink = sin(phase)
      cosk = cos(phase)
      dg1 = dg1 - ewald%weight(itr)*sink*vec
      ds1 = ds1 + cosk*ewald%strain(:, :, itr)
      dg3 = dg3 - ewald%weight3(itr)*sink*vec
      ds3 = ds3 + cosk*ewald%strain3(:, :, itr)
   end do
   k0 = 4.0_wp*pi/ewald%vol*(log(alpha/sqrtpi) + s3_constant)
   ds3 = ds3 - k0*unity
end subroutine get_damat_ko_rec_3d


!> Calculate direct space Ewald contribution for a pair under 3D periodic boundary conditions
subroutine get_damat_dir_3d(rij, alp, trans, dg, ds)
   !> Distance between pair
   real(wp), intent(in) :: rij(3)
   !> Convergence factor
   real(wp), intent(in) :: alp
   !> Translation vectors to consider
   real(wp), intent(in) :: trans(:, :)
   !> Derivative with respect to cartesian displacements
   real(wp), intent(out) :: dg(3)
   !> Derivative with respect to strain deformations
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: vec(3), r1, r2, dtmp, alp2

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   alp2 = alp*alp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      r2 = r1*r1
      dtmp = -erfc(alp*r1)/(r2*r1) - 2*alp*exp(-r2*alp2)/(sqrtpi*r2)
      dg(:) = dg + dtmp * vec
      ds = ds + dtmp*outerprod(vec, vec)
   end do

end subroutine get_damat_dir_3d


!> Calculate Wigner-Seitz averaged short range correction for a pair
subroutine get_damat_wsc_3d(rij, gam, gexp, dg, ds)
   !> Distance between pair
   real(wp), intent(in) :: rij(3)
   !> Hubbard parameter
   real(wp), intent(in) :: gam
   !> Exponent for interaction kernel
   real(wp), intent(in) :: gexp
   !> Derivative with respect to cartesian displacements
   real(wp), intent(out) :: dg(3)
   !> Derivative with respect to strain deformations
   real(wp), intent(out) :: ds(3, 3)

   real(wp) :: r1, r2, dtmp

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   r1 = norm2(rij)
   if (r1 < eps) return
   r2 = r1*r1
   dtmp = 1.0_wp / (r1**gexp + gam**(-gexp))
   dtmp = -r1**(gexp-2.0_wp) * dtmp * dtmp**(1.0_wp/gexp) &
      & + 1.0_wp/(r2*r1)
   dg(:) = dtmp * rij
   ds = dtmp*outerprod(rij, rij)

end subroutine get_damat_wsc_3d


end module tblite_coulomb_charge_effective
