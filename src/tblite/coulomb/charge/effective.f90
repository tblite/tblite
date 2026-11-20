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
   use mctc_io_math, only : matdet_3x3, outerprod
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_charge_type, only : coulomb_charge_type, get_coulomb_matrix_from_blocks
   use tblite_coulomb_ewald, only : get_dir_cutoff, ewald_cache
   use tblite_cutoff, only : get_lattice_points
   use tblite_partition, only : owns_pair, pair_list
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
      !> Evaluate Coulomb matrix
      procedure :: get_coulomb_matrix
      !> Evaluate the locally owned Coulomb blocks directly
      procedure :: get_coulomb_blocks
      !> Evaluate uncontracted derivatives of Coulomb matrix
      procedure :: get_coulomb_derivs
      !> Contract derivatives directly for gradients and strain
      procedure :: get_coulomb_derivs_contracted
   end type effective_coulomb


   !> Geometry shared by all shell combinations of one periodic atom pair.
   type :: ko_distances
      !> Displacements to each periodic image, with shape (3, ntrans)
      real(wp), allocatable :: vec(:, :)
      !> Squared distances to the periodic images
      real(wp), allocatable :: r2(:)
      !> Inverse distances, zero for excluded zero-distance images
      real(wp), allocatable :: invr(:)
      !> Cubes of the inverse distances
      real(wp), allocatable :: invr3(:)
      !> Fifth powers of the inverse distances
      real(wp), allocatable :: invr5(:)
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
   real(wp), parameter :: ko_cutoff = 40.0_wp
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


!> Evaluate coulomb matrix
subroutine get_coulomb_matrix(self, mol, cache, amat)
   !> Instance of the electrostatic container
   class(effective_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   !> Coulomb matrix
   real(wp), contiguous, intent(out) :: amat(:, :)

   call get_coulomb_matrix_from_blocks(self, mol, cache, self%partition, amat)
end subroutine get_coulomb_matrix


!> Evaluate Coulomb blocks without constructing a full shell matrix.
subroutine get_coulomb_blocks(self, mol, cache, pairs, blocks)
   !> Instance of the electrostatic container
   class(effective_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   !> Locally owned atom pairs in row order
   type(pair_list), intent(in) :: pairs
   !> Coulomb blocks indexed by row shell, neighbour shell, and local pair
   real(wp), contiguous, intent(out) :: blocks(:, :, :)

   blocks(:, :, :) = 0.0_wp
   if (any(mol%periodic)) then
      call cache%charge_ewald%update(mol%lattice, cache%alpha, conv, abs(self%gexp-2.0_wp) < eps)
      call get_amat_3d(mol, self%nshell, self%hubbard, self%gexp, &
         & cache%wsc, cache%alpha, cache%charge_ewald, blocks, pairs)
   else
      call get_amat_0d(mol, self%nshell, self%hubbard, self%gexp, blocks, pairs)
   end if
end subroutine get_coulomb_blocks


!> Get real lattice translations
subroutine get_dir_trans(lattice, alpha, conv, trans)
   !> Lattice parameters
   real(wp), intent(in) :: lattice(:, :)
   !> Parameter for Ewald summation
   real(wp), intent(in) :: alpha
   !> Tolerance for Ewald summation
   real(wp), intent(in) :: conv
   !> Translation vectors
   real(wp), allocatable, intent(out) :: trans(:, :)

   call get_lattice_points([.true.], lattice, get_dir_cutoff(alpha, conv), trans)

end subroutine get_dir_trans


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
   !> Cached reciprocal-space coefficients
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
   real(wp), allocatable :: dtrans(:, :)

   if (abs(gexp - 2.0_wp) < eps) then
      call get_amat_ko_3d(mol, nshell, hubbard, alpha, ewald, amat, pairs)
      return
   end if

   call get_dir_trans(mol%lattice, alpha, conv, dtrans)

   !$omp parallel do default(none) schedule(static, 1) shared(amat) &
   !$omp shared(mol, nshell, hubbard, gexp, wsc, dtrans, ewald, alpha, pairs) &
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
         call get_amat_dir_3d(vec, alpha, dtrans, dtmp)
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
   !> Cached reciprocal-space coefficients
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
   real(wp) :: vec(3), gam, val, vol, s1, s3, sr
   real(wp), allocatable :: dtrans(:, :), strans(:, :)
   type(ko_distances) :: distances

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, alpha, conv, dtrans)
   call get_lattice_points([.true.], mol%lattice, ko_cutoff, strans)

   !$omp parallel default(none) shared(amat) &
   !$omp shared(mol, nshell, hubbard, alpha, vol, dtrans, ewald, strans, pairs) &
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
         call get_ko_distances(vec, strans, distances)
         call get_amat_ko_ewald_3d(vec, alpha, vol, dtrans, ewald, distances, &
            & iat == jat, s1, s3)
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
subroutine get_amat_ko_ewald_3d(rij, alpha, vol, dtrans, ewald, distances, &
      & onsite, s1, s3)
   !> Displacement between the atom pair
   real(wp), intent(in) :: rij(3)
   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha
   !> Unit cell volume
   real(wp), intent(in) :: vol
   !> Real-space Ewald translation vectors
   real(wp), intent(in) :: dtrans(:, :)
   !> Pair geometry for the short-range translations
   type(ko_distances), intent(in) :: distances
   !> Cached reciprocal-space coefficients
   type(ewald_cache), intent(in) :: ewald
   !> Whether both indices belong to the same atom
   logical, intent(in) :: onsite
   !> Generalized Ewald sum for 1/r
   real(wp), intent(out) :: s1
   !> Generalized Ewald sum for 1/r^3
   real(wp), intent(out) :: s3

   integer :: itr
   real(wp) :: r1, r2, k, tmp, phase

   call get_amat_dir_3d(rij, alpha, dtrans, s1)
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
   s3 = s3 + 4.0_wp*pi/vol*(log(k) + s3_constant)
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
   real(wp) :: r1, g1

   ntrans = size(trans, 2)
   if (.not.allocated(distances%r2)) then
      allocate(distances%vec(3, ntrans), distances%r2(ntrans), distances%invr(ntrans), &
         & distances%invr3(ntrans), distances%invr5(ntrans))
   end if
   do itr = 1, ntrans
      distances%vec(:, itr) = rij + trans(:, itr)
      distances%r2(itr) = dot_product(distances%vec(:, itr), distances%vec(:, itr))
      r1 = sqrt(distances%r2(itr))
      g1 = 0.0_wp
      if (r1 >= eps) g1 = 1.0_wp/r1
      distances%invr(itr) = g1
      distances%invr3(itr) = g1*g1*g1
      distances%invr5(itr) = distances%invr3(itr)*g1*g1
   end do
end subroutine get_ko_distances

!> Release each component independently
subroutine delete_ko_distances(self)
   !> Pair geometry to release
   class(ko_distances), intent(inout) :: self

   if (allocated(self%vec)) deallocate(self%vec)
   if (allocated(self%r2)) deallocate(self%r2)
   if (allocated(self%invr)) deallocate(self%invr)
   if (allocated(self%invr3)) deallocate(self%invr3)
   if (allocated(self%invr5)) deallocate(self%invr5)
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

!> Calculate reciprocal space contributions for a pair under 3D periodic boundary conditions
subroutine get_amat_rec_3d(rij, ewald, amat)
   !> Displacement between the atom pair
   real(wp), intent(in) :: rij(3)
   !> Cached reciprocal-space coefficients
   type(ewald_cache), intent(in) :: ewald
   !> Coulomb matrix elements
   real(wp), intent(out) :: amat
   integer :: itr

   amat = 0.0_wp
   do itr = 1, size(ewald%weight)
      amat = amat + cos(dot_product(rij, ewald%vec(:, itr)))*ewald%weight(itr)
   end do
end subroutine get_amat_rec_3d


!> Evaluate uncontracted derivatives of Coulomb matrix
subroutine get_coulomb_derivs(self, mol, cache, qat, qsh, dadr, dadL, atrace)
   !> Instance of the electrostatic container
   class(effective_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Shell-resolved partial charges
   real(wp), intent(in) :: qsh(:)
   !> Derivative of interactions with respect to cartesian displacements
   real(wp), contiguous, intent(out) :: dadr(:, :, :)
   !> Derivative of interactions with respect to strain deformations
   real(wp), contiguous, intent(out) :: dadL(:, :, :)
   !> On-site derivatives with respect to cartesian displacements
   real(wp), contiguous, intent(out) :: atrace(:, :)

   if (any(mol%periodic)) then
      call cache%charge_ewald%update(mol%lattice, cache%alpha, conv, abs(self%gexp-2.0_wp) < eps)
   end if
   if (self%shell_resolved) then
      call get_damat(self, mol, cache, qsh, dadr, dadL, atrace)
   else
      call get_damat(self, mol, cache, qat, dadr, dadL, atrace)
   end if

end subroutine get_coulomb_derivs


!> Accumulate charge-contracted derivatives without a potential Jacobian
subroutine get_coulomb_derivs_contracted(self, mol, cache, qat, qsh, gradient, sigma)
   !> Instance of the electrostatic container
   class(effective_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Shell-resolved partial charges
   real(wp), intent(in) :: qsh(:)
   !> Molecular gradient to accumulate into
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives to accumulate into
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   if (any(mol%periodic)) then
      call cache%charge_ewald%update(mol%lattice, cache%alpha, conv, abs(self%gexp-2.0_wp) < eps)
   end if
   if (self%shell_resolved) then
      call get_damat_contracted(self, mol, cache, qsh, gradient, sigma)
   else
      call get_damat_contracted(self, mol, cache, qat, gradient, sigma)
   end if
end subroutine get_coulomb_derivs_contracted


!> Assemble the full potential derivatives from individual atom pairs
subroutine get_damat(self, mol, cache, qvec, dadr, dadL, atrace)
   !> Instance of the electrostatic container
   class(effective_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(in) :: cache
   !> Atomic or shell-resolved partial charges
   real(wp), intent(in) :: qvec(:)
   !> Derivative of interactions with respect to cartesian displacements
   real(wp), intent(out) :: dadr(:, :, :)
   !> Derivative of interactions with respect to strain deformations
   real(wp), intent(out) :: dadL(:, :, :)
   !> On-site derivatives with respect to cartesian displacements
   real(wp), intent(out) :: atrace(:, :)

   integer :: iat, jat, ii, jj, ish, jsh
   real(wp) :: dG(3, size(self%hubbard, 1), size(self%hubbard, 1))
   real(wp) :: dS(3, 3, size(self%hubbard, 1), size(self%hubbard, 1))
   real(wp), allocatable :: dtrans(:, :), strans(:, :)
   real(wp), allocatable :: itrace(:, :), didr(:, :, :), didL(:, :, :)
   type(ko_distances) :: distances

   atrace = 0.0_wp
   dadr = 0.0_wp
   dadL = 0.0_wp
   if (any(mol%periodic)) then
      call get_dir_trans(mol%lattice, cache%alpha, conv, dtrans)
      if (abs(self%gexp-2.0_wp) < eps) then
         call get_lattice_points([.true.], mol%lattice, ko_cutoff, strans)
      end if
   end if

   !$omp parallel default(none) shared(atrace, dadr, dadL, self, mol, cache, qvec, dtrans, strans) &
   !$omp private(iat, jat, ii, jj, ish, jsh, dG, dS, itrace, didr, didL, distances)
   allocate(itrace, source=atrace)
   allocate(didr, source=dadr)
   allocate(didL, source=dadL)
   !$omp do schedule(static, 1)
   do iat = 1, mol%nat
      ii = self%offset(iat)
      do jat = 1, iat-1
         if (.not.owns_pair(self%partition, iat, jat)) cycle
         jj = self%offset(jat)
         call get_damat_pair(self, mol, cache, iat, jat, dtrans, strans, distances, dG, dS)
         do ish = 1, self%nshell(iat)
            do jsh = 1, self%nshell(jat)
               itrace(:, ii+ish) = itrace(:, ii+ish) + dG(:, jsh, ish)*qvec(jj+jsh)
               itrace(:, jj+jsh) = itrace(:, jj+jsh) - dG(:, jsh, ish)*qvec(ii+ish)
               didr(:, iat, jj+jsh) = didr(:, iat, jj+jsh) + dG(:, jsh, ish)*qvec(ii+ish)
               didr(:, jat, ii+ish) = didr(:, jat, ii+ish) - dG(:, jsh, ish)*qvec(jj+jsh)
               didL(:, :, jj+jsh) = didL(:, :, jj+jsh) + dS(:, :, jsh, ish)*qvec(ii+ish)
               didL(:, :, ii+ish) = didL(:, :, ii+ish) + dS(:, :, jsh, ish)*qvec(jj+jsh)
            end do
         end do
      end do

      if (.not.any(mol%periodic) .or. .not.owns_pair(self%partition, iat, iat)) cycle
      call get_damat_pair(self, mol, cache, iat, iat, dtrans, strans, distances, dG, dS)
      do ish = 1, self%nshell(iat)
         do jsh = 1, ish-1
            didL(:, :, ii+jsh) = didL(:, :, ii+jsh) + dS(:, :, jsh, ish)*qvec(ii+ish)
            didL(:, :, ii+ish) = didL(:, :, ii+ish) + dS(:, :, jsh, ish)*qvec(ii+jsh)
         end do
         didL(:, :, ii+ish) = didL(:, :, ii+ish) + dS(:, :, ish, ish)*qvec(ii+ish)
      end do
   end do
   !$omp end do
   !$omp critical (get_damat_)
   atrace = atrace + itrace
   dadr = dadr + didr
   dadL = dadL + didL
   !$omp end critical (get_damat_)
   deallocate(itrace, didr, didL)
   call distances%delete()
   !$omp end parallel
end subroutine get_damat


!> Accumulate gradient and strain without allocating full potential derivatives
subroutine get_damat_contracted(self, mol, cache, qvec, gradient, sigma)
   !> Instance of the electrostatic container
   class(effective_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(in) :: cache
   !> Atomic or shell-resolved partial charges
   real(wp), intent(in) :: qvec(:)
   !> Molecular gradient to accumulate into
   real(wp), intent(inout) :: gradient(:, :)
   !> Strain derivatives to accumulate into
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, ii, jj, ish, jsh
   real(wp) :: qq
   real(wp) :: dG(3, size(self%hubbard, 1), size(self%hubbard, 1))
   real(wp) :: dS(3, 3, size(self%hubbard, 1), size(self%hubbard, 1))
   real(wp), allocatable :: dtrans(:, :), strans(:, :)
   type(ko_distances) :: distances

   if (any(mol%periodic)) then
      call get_dir_trans(mol%lattice, cache%alpha, conv, dtrans)
      if (abs(self%gexp-2.0_wp) < eps) then
         call get_lattice_points([.true.], mol%lattice, ko_cutoff, strans)
      end if
   end if

   !$omp parallel default(none) reduction(+:gradient, sigma) &
   !$omp shared(self, mol, cache, qvec, dtrans, strans) &
   !$omp private(iat, jat, ii, jj, ish, jsh, qq, dG, dS, distances)
   !$omp do schedule(static, 1)
   do iat = 1, mol%nat
      ii = self%offset(iat)
      do jat = 1, iat-1
         if (.not.owns_pair(self%partition, iat, jat)) cycle
         jj = self%offset(jat)
         call get_damat_pair(self, mol, cache, iat, jat, dtrans, strans, distances, dG, dS)
         do ish = 1, self%nshell(iat)
            do jsh = 1, self%nshell(jat)
               qq = qvec(ii+ish)*qvec(jj+jsh)
               gradient(:, iat) = gradient(:, iat) + dG(:, jsh, ish)*qq
               gradient(:, jat) = gradient(:, jat) - dG(:, jsh, ish)*qq
               sigma = sigma + dS(:, :, jsh, ish)*qq
            end do
         end do
      end do

      if (.not.any(mol%periodic) .or. .not.owns_pair(self%partition, iat, iat)) cycle
      call get_damat_pair(self, mol, cache, iat, iat, dtrans, strans, distances, dG, dS)
      do ish = 1, self%nshell(iat)
         do jsh = 1, ish-1
            sigma = sigma + dS(:, :, jsh, ish)*(qvec(ii+ish)*qvec(ii+jsh))
         end do
         sigma = sigma + dS(:, :, ish, ish)*(0.5_wp*qvec(ii+ish)**2)
      end do
   end do
   !$omp end do
   call distances%delete()
   !$omp end parallel
end subroutine get_damat_contracted


!> Evaluate derivatives for one atom pair, independently of their contraction.
!> Only jsh <= ish is evaluated for on-site pairs; their Cartesian part is unused.
subroutine get_damat_pair(self, mol, cache, iat, jat, dtrans, strans, distances, dG, dS)
   !> Instance of the electrostatic container
   class(effective_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(in) :: cache
   !> First atom
   integer, intent(in) :: iat
   !> Second atom
   integer, intent(in) :: jat
   !> Real-space Ewald translations, allocated only for periodic systems
   real(wp), allocatable, intent(in) :: dtrans(:, :)
   !> Short-range translations, allocated only for periodic Klopman-Ohno sums
   real(wp), allocatable, intent(in) :: strans(:, :)
   !> Reusable pair geometry for the Klopman-Ohno short-range sum
   type(ko_distances), intent(inout) :: distances
   !> Cartesian derivatives, indexed by component and shells (jsh, ish)
   real(wp), intent(out) :: dG(:, :, :)
   !> Strain derivatives, indexed by components and shells (jsh, ish)
   real(wp), intent(out) :: dS(:, :, :, :)

   integer :: izp, jzp, ish, jsh
   real(wp) :: vec(3), gam, r1, dtmp
   real(wp) :: dG1(3), dS1(3, 3), dG3(3), dS3(3, 3), dGsr(3), dSsr(3, 3)

   izp = mol%id(iat)
   jzp = mol%id(jat)
   vec = mol%xyz(:, iat) - mol%xyz(:, jat)
   dG = 0.0_wp
   dS = 0.0_wp
   if (.not.any(mol%periodic)) then
      r1 = norm2(vec)
      do ish = 1, self%nshell(iat)
         do jsh = 1, self%nshell(jat)
            gam = self%hubbard(jsh, ish, jzp, izp)
            dtmp = 1.0_wp/(r1**self%gexp + gam**(-self%gexp))
            dtmp = -r1**(self%gexp-2.0_wp)*dtmp*dtmp**(1.0_wp/self%gexp)
            dG(:, jsh, ish) = dtmp*vec
            dS(:, :, jsh, ish) = outerprod(dG(:, jsh, ish), vec)
         end do
      end do
   else if (abs(self%gexp-2.0_wp) < eps) then
      call get_ko_distances(vec, strans, distances)
      call get_damat_ko_ewald_3d(vec, cache%alpha, abs(matdet_3x3(mol%lattice)), dtrans, &
         & cache%charge_ewald, distances, dG1, dS1, dG3, dS3)
      do ish = 1, self%nshell(iat)
         do jsh = 1, merge(ish, self%nshell(jat), iat == jat)
            gam = self%hubbard(jsh, ish, jzp, izp)
            call get_damat_ko_short_3d(gam, distances, dGsr, dSsr)
            dG(:, jsh, ish) = dG1 + dGsr - 0.5_wp*dG3/(gam*gam)
            dS(:, :, jsh, ish) = dS1 + dSsr - 0.5_wp*dS3/(gam*gam)
         end do
      end do
   else
      block
         integer :: img
         logical :: need_weight_energy
         real(wp) :: stmp, dGr(3), dSr(3, 3)
         real(wp) :: weight(size(cache%wsc%tridx, 1))
         real(wp) :: dwdr(3, size(cache%wsc%tridx, 1)), dwdL(3, 3, size(cache%wsc%tridx, 1))

         call get_wignerseitz_weights(cache%wsc, jat, iat, vec, weight, dwdr, dwdL)
         need_weight_energy = any(dwdr(:, :cache%wsc%nimg(jat, iat)) /= 0.0_wp) &
            & .or. any(dwdL(:, :, :cache%wsc%nimg(jat, iat)) /= 0.0_wp)
         call get_damat_dir_3d(vec, cache%alpha, dtrans, dG1, dS1)
         call get_damat_rec_3d(vec, cache%charge_ewald, dGr, dSr)
         do img = 1, cache%wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - cache%wsc%trans(:, cache%wsc%tridx(img, jat, iat))
            do ish = 1, self%nshell(iat)
               do jsh = 1, merge(ish, self%nshell(jat), iat == jat)
                  gam = self%hubbard(jsh, ish, jzp, izp)
                  stmp = 0.0_wp
                  if (need_weight_energy) call get_amat_wsc_3d(vec, gam, self%gexp, stmp)
                  call get_damat_wsc_3d(vec, gam, self%gexp, dGsr, dSsr)
                  dG(:, jsh, ish) = dG(:, jsh, ish) + (dG1 + dGr + dGsr)*weight(img) + stmp*dwdr(:, img)
                  dS(:, :, jsh, ish) = dS(:, :, jsh, ish) + (dS1 + dSr + dSsr)*weight(img) + stmp*dwdL(:, :, img)
               end do
            end do
         end do
      end block
   end if
end subroutine get_damat_pair


!> Derivatives of the shell-independent generalized Ewald sums
subroutine get_damat_ko_ewald_3d(rij, alpha, vol, dtrans, ewald, distances, &
      & dg1, ds1, dg3, ds3)
   !> Displacement between the atom pair
   real(wp), intent(in) :: rij(3)
   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha
   !> Unit cell volume
   real(wp), intent(in) :: vol
   !> Real-space Ewald translation vectors
   real(wp), intent(in) :: dtrans(:, :)
   !> Pair geometry for the short-range translations
   type(ko_distances), intent(in) :: distances
   !> Cached reciprocal-space coefficients
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
   real(wp) :: vec(3), r1, r2, phase, sink, cosk, dtmp, k0
   real(wp) :: dgtmp(3), dstmp(3, 3)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))

   call get_damat_dir_3d(rij, alpha, dtrans, dg1, ds1)

   dg3 = 0.0_wp
   ds3 = 0.0_wp
   do itr = 1, size(distances%r2)
      if (distances%invr(itr) == 0.0_wp) cycle
      vec = distances%vec(:, itr)
      r2 = distances%r2(itr)
      r1 = sqrt(r2)
      dtmp = -3.0_wp*erfc(alpha*r1)/(r2*r2*r1) &
         & - 6.0_wp*alpha*exp(-alpha*alpha*r2)/(sqrtpi*r2*r2) &
         & - 4.0_wp*alpha**3*exp(-alpha*alpha*r2)/(sqrtpi*r2)
      dg3 = dg3 + dtmp*vec
      ds3 = ds3 + dtmp*outerprod(vec, vec)
   end do

   dgtmp = 0.0_wp
   dstmp = 0.0_wp
   do itr = 1, size(ewald%weight)
      vec = ewald%vec(:, itr)
      phase = dot_product(rij, vec)
      sink = sin(phase)
      cosk = cos(phase)
      dgtmp = dgtmp - ewald%weight(itr)*sink*vec
      dstmp = dstmp + cosk*ewald%strain(:, :, itr)
      dg3 = dg3 - ewald%weight3(itr)*sink*vec
      ds3 = ds3 + cosk*ewald%strain3(:, :, itr)
   end do
   dg1 = dg1 + dgtmp
   ds1 = ds1 + dstmp
   k0 = 4.0_wp*pi/vol*(log(alpha/sqrtpi) + s3_constant)
   ds3 = ds3 - k0*unity
end subroutine get_damat_ko_ewald_3d

!> Derivatives of the hardness-dependent short-range residual
subroutine get_damat_ko_short_3d(gam, distances, dg, ds)
   !> Averaged Hubbard parameter
   real(wp), intent(in) :: gam
   !> Pair geometry for the short-range translations
   type(ko_distances), intent(in) :: distances
   !> Derivative with respect to cartesian displacements
   real(wp), intent(out) :: dg(3)
   !> Derivative with respect to strain deformations
   real(wp), intent(out) :: ds(3, 3)
   integer :: itr
   real(wp) :: gam2, tmp, dtmp, vec(3)

   gam2 = 1.0_wp/(gam*gam)
   dg = 0.0_wp
   ds = 0.0_wp
   do itr = 1, size(distances%r2)
      if (distances%invr(itr) == 0.0_wp) cycle
      vec = distances%vec(:, itr)
      tmp = 1.0_wp/(distances%r2(itr) + gam2)
      dtmp = -tmp*sqrt(tmp) + distances%invr3(itr) - 1.5_wp*gam2*distances%invr5(itr)
      dg = dg + dtmp*vec
      ds = ds + dtmp*outerprod(vec, vec)
   end do
end subroutine get_damat_ko_short_3d

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

!> Calculate reciprocal space contributions for a pair under 3D periodic boundary conditions
subroutine get_damat_rec_3d(rij, ewald, dg, ds)
   !> Displacement between the atom pair
   real(wp), intent(in) :: rij(3)
   !> Cached reciprocal-space coefficients
   type(ewald_cache), intent(in) :: ewald
   !> Derivative with respect to cartesian displacements
   real(wp), intent(out) :: dg(3)
   !> Derivative with respect to strain deformations
   real(wp), intent(out) :: ds(3, 3)
   integer :: itr
   real(wp) :: phase

   dg = 0.0_wp
   ds = 0.0_wp
   do itr = 1, size(ewald%weight)
      phase = dot_product(rij, ewald%vec(:, itr))
      dg = dg - sin(phase)*ewald%weight(itr)*ewald%vec(:, itr)
      ds = ds + cos(phase)*ewald%strain(:, :, itr)
   end do
end subroutine get_damat_rec_3d


end module tblite_coulomb_charge_effective
