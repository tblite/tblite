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
   use mctc_io_math, only : matdet_3x3
   use tblite_container_cache, only : container_cache
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_charge_type, only : coulomb_charge_type
   use tblite_coulomb_ewald, only : get_dir_cutoff, ewald_cache
   use tblite_cutoff, only : get_lattice_points
   use tblite_partition, only : work_partition, owns_pair
   use tblite_wavefunction_type, only : wavefunction_type
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
      !> Evaluate uncontracted derivatives of Coulomb matrix
      procedure :: get_coulomb_derivs
      !> Contract derivatives directly for forces and strain
      procedure :: get_gradient
   end type effective_coulomb


   !> Thread-local output of the same pair derivative kernels: either forces
   !> and strain, or the full potential Jacobian required by response methods.
   type :: derivative_buffer
      real(wp), allocatable :: dr(:, :, :), dL(:, :, :), trace(:, :), gradient(:, :)
      real(wp) :: sigma(3, 3)
   end type derivative_buffer

   !> Geometry shared by all shell combinations of one periodic atom pair.
   type :: ko_distances
      real(wp), allocatable :: vec(:, :), r2(:), invr(:), invr3(:), invr5(:)
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
   amat(:, :) = 0.0_wp

    if (any(mol%periodic)) then
       call cache%charge_ewald%update(mol%lattice, cache%alpha, conv, abs(self%gexp-2.0_wp) < eps)
       call get_amat_3d(mol, self%nshell, self%offset, self%hubbard, self%gexp, &
          & cache%wsc, cache%alpha, cache%charge_ewald, amat, self%partition)
    else
      call get_amat_0d(mol, self%nshell, self%offset, self%hubbard, self%gexp, amat, &
         & self%partition)
   end if

end subroutine get_coulomb_matrix


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
subroutine get_amat_0d(mol, nshell, offset, hubbard, gexp, amat, partition)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each atom
   integer, intent(in) :: nshell(:)
   !> Index offset for each shell
   integer, intent(in) :: offset(:)
   !> Hubbard parameter parameter for each shell
   real(wp), intent(in) :: hubbard(:, :, :, :)
   !> Exponent of Coulomb kernel
   real(wp), intent(in) :: gexp
   !> Coulomb matrix
   real(wp), intent(inout) :: amat(:, :)
   !> Share of the atom pairs evaluated here, absent selects the complete work
   type(work_partition), intent(in), optional :: partition

   integer :: iat, jat, izp, jzp, ii, jj, ish, jsh
   real(wp) :: vec(3), r1, r1g, gam, tmp

   ! Cyclic rows balance the triangular pair loop without changing MPI ownership.
   !$omp parallel do default(none) schedule(static, 1) &
   !$omp shared(amat, mol, nshell, offset, hubbard, gexp, partition) &
   !$omp private(iat, izp, ii, ish, jat, jzp, jj, jsh, gam, vec, r1, r1g, tmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = offset(iat)
      do jat = 1, iat-1
         if (.not.owns_pair(partition, iat, jat)) cycle
         jzp = mol%id(jat)
         jj = offset(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         r1g = r1**gexp
         do ish = 1, nshell(iat)
            do jsh = 1, nshell(jat)
               gam = hubbard(jsh, ish, jzp, izp)
               tmp = 1.0_wp/(r1g + gam**(-gexp))**(1.0_wp/gexp)
               amat(jj+jsh, ii+ish) = amat(jj+jsh, ii+ish) + tmp
               amat(ii+ish, jj+jsh) = amat(ii+ish, jj+jsh) + tmp
            end do
         end do
      end do
      if (.not.owns_pair(partition, iat, iat)) cycle
      do ish = 1, nshell(iat)
         do jsh = 1, ish-1
            gam = hubbard(jsh, ish, izp, izp)
            amat(ii+jsh, ii+ish) = amat(ii+jsh, ii+ish) + gam
            amat(ii+ish, ii+jsh) = amat(ii+ish, ii+jsh) + gam
         end do
         amat(ii+ish, ii+ish) = amat(ii+ish, ii+ish) + hubbard(ish, ish, izp, izp)
      end do
   end do

end subroutine get_amat_0d

!> Evaluate the coulomb matrix for 3D systems
subroutine get_amat_3d(mol, nshell, offset, hubbard, gexp, wsc, alpha, ewald, amat, partition)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   type(ewald_cache), intent(in) :: ewald
   !> Number of shells per atom
   integer, intent(in) :: nshell(:)
   !> Index offset for each atom
   integer, intent(in) :: offset(:)
   !> Hubbard parameter of the shells
   real(wp), intent(in) :: hubbard(:, :, :, :)
   !> Exponent of the interaction kernel
   real(wp), intent(in) :: gexp
   !> Wigner-Seitz cell
   type(wignerseitz_cell), intent(in) :: wsc
   !> Convergence factor for Ewald sum
   real(wp), intent(in) :: alpha
   !> Coulomb matrix
   real(wp), intent(inout) :: amat(:, :)
   !> Share of the atom pairs evaluated here, absent selects the complete work
   type(work_partition), intent(in), optional :: partition

   integer :: iat, jat, izp, jzp, img, ii, jj, ish, jsh
   real(wp) :: vec(3), gam, dtmp, rtmp, stmp, aval
   real(wp) :: weight(size(wsc%tridx, 1))
   real(wp), allocatable :: dtrans(:, :)

   if (abs(gexp - 2.0_wp) < eps) then
      call get_amat_ko_3d(mol, nshell, offset, hubbard, alpha, ewald, amat, partition)
      return
   end if

   call get_dir_trans(mol%lattice, alpha, conv, dtrans)

   !$omp parallel do default(none) schedule(static, 1) shared(amat) &
   !$omp shared(mol, nshell, offset, hubbard, gexp, wsc, dtrans, ewald, alpha, partition) &
   !$omp private(iat, izp, jat, jzp, ii, jj, ish, jsh, gam, weight, vec, dtmp, rtmp, stmp, aval)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = offset(iat)
      do jat = 1, iat-1
         if (.not.owns_pair(partition, iat, jat)) cycle
         jzp = mol%id(jat)
         jj = offset(jat)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         call get_wignerseitz_weights(wsc, jat, iat, vec, weight)
         call get_amat_dir_3d(vec, alpha, dtrans, dtmp)
         call get_amat_rec_3d(vec, ewald, rtmp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
            do ish = 1, nshell(iat)
               do jsh = 1, nshell(jat)
                  gam = hubbard(jsh, ish, jzp, izp)
                  call get_amat_wsc_3d(vec, gam, gexp, stmp)
                  aval = (dtmp + rtmp + stmp) * weight(img)
                  amat(jj+jsh, ii+ish) = amat(jj+jsh, ii+ish) + aval
                  amat(ii+ish, jj+jsh) = amat(ii+ish, jj+jsh) + aval
               end do
            end do
         end do
      end do

      if (.not.owns_pair(partition, iat, iat)) cycle
      vec = 0.0_wp
      call get_wignerseitz_weights(wsc, iat, iat, vec, weight)
      call get_amat_dir_3d(vec, alpha, dtrans, dtmp)
      call get_amat_rec_3d(vec, ewald, rtmp)
      rtmp = rtmp - 2 * alpha / sqrtpi
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         do ish = 1, nshell(iat)
            do jsh = 1, ish-1
               gam = hubbard(jsh, ish, izp, izp)
               call get_amat_wsc_3d(vec, gam, gexp, stmp)
               aval = (dtmp + rtmp + stmp + gam) * weight(img)
               amat(ii+jsh, ii+ish) = amat(ii+jsh, ii+ish) + aval
               amat(ii+ish, ii+jsh) = amat(ii+ish, ii+jsh) + aval
            end do
            gam = hubbard(ish, ish, izp, izp)
            call get_amat_wsc_3d(vec, gam, gexp, stmp)
            aval = (dtmp + rtmp + stmp + gam) * weight(img)
            amat(ii+ish, ii+ish) = amat(ii+ish, ii+ish) + aval
         end do
      end do

   end do

end subroutine get_amat_3d

!> Evaluate the periodic Klopman-Ohno matrix with the generalized Ewald sum
subroutine get_amat_ko_3d(mol, nshell, offset, hubbard, alpha, ewald, amat, partition)
   type(structure_type), intent(in) :: mol
   type(ewald_cache), intent(in) :: ewald
   integer, intent(in) :: nshell(:), offset(:)
   real(wp), intent(in) :: hubbard(:, :, :, :), alpha
   real(wp), intent(inout) :: amat(:, :)
   type(work_partition), intent(in), optional :: partition

   integer :: iat, jat, izp, jzp, ii, jj, ish, jsh
   real(wp) :: vec(3), gam, val, vol, s1, s3, sr
   real(wp), allocatable :: dtrans(:, :), strans(:, :)
   type(ko_distances) :: distances

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, alpha, conv, dtrans)
   call get_lattice_points([.true.], mol%lattice, ko_cutoff, strans)

   !$omp parallel default(none) shared(amat) &
   !$omp shared(mol, nshell, offset, hubbard, alpha, vol, dtrans, ewald, strans, partition) &
   !$omp private(iat, izp, ii, ish, jat, jzp, jj, jsh, gam, vec, val, s1, s3, sr, distances)
   !$omp do schedule(static, 1)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = offset(iat)
      do jat = 1, iat-1
         if (.not.owns_pair(partition, iat, jat)) cycle
         jzp = mol%id(jat)
         jj = offset(jat)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         call get_ko_distances(vec, strans, distances)
         call get_amat_ko_ewald_3d(vec, alpha, vol, dtrans, ewald, distances, &
            & .false., s1, s3)
         do ish = 1, nshell(iat)
            do jsh = 1, nshell(jat)
               gam = hubbard(jsh, ish, jzp, izp)
               call get_amat_ko_short_3d(gam, distances, sr)
               val = s1 + sr - 0.5_wp*s3/(gam*gam)
               amat(jj+jsh, ii+ish) = amat(jj+jsh, ii+ish) + val
               amat(ii+ish, jj+jsh) = amat(ii+ish, jj+jsh) + val
            end do
         end do
      end do

      if (.not.owns_pair(partition, iat, iat)) cycle
      vec = 0.0_wp
      call get_ko_distances(vec, strans, distances)
      call get_amat_ko_ewald_3d(vec, alpha, vol, dtrans, ewald, distances, &
         & .true., s1, s3)
      do ish = 1, nshell(iat)
         do jsh = 1, ish-1
            gam = hubbard(jsh, ish, izp, izp)
            call get_amat_ko_short_3d(gam, distances, sr)
            val = s1 + sr - 0.5_wp*s3/(gam*gam) + gam
            amat(ii+jsh, ii+ish) = amat(ii+jsh, ii+ish) + val
            amat(ii+ish, ii+jsh) = amat(ii+ish, ii+jsh) + val
         end do
         gam = hubbard(ish, ish, izp, izp)
         call get_amat_ko_short_3d(gam, distances, sr)
         val = s1 + sr - 0.5_wp*s3/(gam*gam) + gam
         amat(ii+ish, ii+ish) = amat(ii+ish, ii+ish) + val
      end do
   end do
   !$omp end do
   if (allocated(distances%r2)) then
      deallocate(distances%vec, distances%r2, distances%invr, distances%invr3, distances%invr5)
   end if
   !$omp end parallel
end subroutine get_amat_ko_3d

!> Shell-independent generalized Ewald sums for one atom pair
subroutine get_amat_ko_ewald_3d(rij, alpha, vol, dtrans, ewald, distances, &
      & onsite, s1, s3)
   real(wp), intent(in) :: rij(3), alpha, vol
   real(wp), intent(in) :: dtrans(:, :)
   type(ko_distances), intent(in) :: distances
   type(ewald_cache), intent(in) :: ewald
   logical, intent(in) :: onsite
   real(wp), intent(out) :: s1, s3

   integer :: itr
   real(wp) :: vec(3), r1, r2, k, tmp, phase

   call get_amat_dir_3d(rij, alpha, dtrans, s1)
   if (onsite) s1 = s1 - 2.0_wp*alpha/sqrtpi

   s3 = 0.0_wp
   do itr = 1, size(distances%r2)
      if (distances%invr(itr) == 0.0_wp) cycle
      vec = distances%vec(:, itr)
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
   real(wp), intent(in) :: rij(3), trans(:, :)
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

!> Hardness-dependent short-range Klopman-Ohno residual
subroutine get_amat_ko_short_3d(gam, distances, val)
   real(wp), intent(in) :: gam
   type(ko_distances), intent(in) :: distances
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
   real(wp), intent(in) :: rij(3)
   type(ewald_cache), intent(in) :: ewald
   real(wp), intent(out) :: amat
   integer :: itr

   amat = 0.0_wp
   do itr = 1, size(ewald%weight)
      amat = amat + cos(dot_product(rij, ewald%vec(:, itr)))*ewald%weight(itr)
   end do
end subroutine get_amat_rec_3d


!> Accumulate forces without materializing a potential Jacobian.
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   class(effective_coulomb), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: cache
   type(wavefunction_type), intent(in) :: wfn
   real(wp), contiguous, intent(inout) :: gradient(:, :), sigma(:, :)

   select type(ptr => cache%raw)
   type is(coulomb_cache)
      call get_derivatives(self, mol, ptr, wfn%qat(:, 1), wfn%qsh(:, 1), &
         & gradient=gradient, sigma=sigma)
   end select
end subroutine get_gradient

!> Both consumers use identical pair derivatives and work ownership.
subroutine get_derivatives(self, mol, cache, qat, qsh, dadr, dadL, atrace, gradient, sigma)
   class(effective_coulomb), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(coulomb_cache), intent(inout) :: cache
   real(wp), intent(in) :: qat(:), qsh(:)
   real(wp), intent(out), optional :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
   real(wp), intent(inout), optional :: gradient(:, :), sigma(:, :)
   real(wp) :: qvec(sum(self%nshell))

   if (self%shell_resolved) then
      qvec(:) = qsh
   else
      qvec(:) = qat
   end if
   if (any(mol%periodic)) then
      call cache%charge_ewald%update(mol%lattice, cache%alpha, conv, abs(self%gexp-2.0_wp) < eps)
      call get_damat_3d(mol, self%nshell, self%offset, self%hubbard, self%gexp, &
         & cache%wsc, cache%alpha, cache%charge_ewald, qvec, dadr, dadL, atrace, &
         & self%partition, gradient, sigma)
   else
      call get_damat_0d(mol, self%nshell, self%offset, self%hubbard, self%gexp, qvec, &
         & dadr, dadL, atrace, self%partition, gradient, sigma)
   end if
end subroutine get_derivatives

subroutine new_derivative_buffer(buffer, nat, nsh, contracted)
   type(derivative_buffer), intent(out) :: buffer
   integer, intent(in) :: nat, nsh
   logical, intent(in) :: contracted

   if (contracted) then
      allocate(buffer%gradient(3, nat), source=0.0_wp)
      buffer%sigma = 0.0_wp
   else
      allocate(buffer%dr(3, nat, nsh), buffer%dL(3, 3, nsh), buffer%trace(3, nsh), source=0.0_wp)
   end if
end subroutine new_derivative_buffer

!> Add an unordered shell pair; diagonal pairs carry half the strain weight.
pure subroutine add_derivative(buffer, iat, jat, ii, jj, qvec, dg, ds)
   type(derivative_buffer), intent(inout) :: buffer
   integer, intent(in) :: iat, jat, ii, jj
   real(wp), intent(in) :: qvec(:), dg(3), ds(3, 3)
   real(wp) :: qq

   if (allocated(buffer%gradient)) then
      qq = qvec(ii)*qvec(jj)
      if (iat /= jat) then
         buffer%gradient(:, iat) = buffer%gradient(:, iat) + dg*qq
         buffer%gradient(:, jat) = buffer%gradient(:, jat) - dg*qq
      end if
      if (ii == jj) qq = 0.5_wp*qq
      buffer%sigma = buffer%sigma + ds*qq
   else
      if (iat /= jat) then
         buffer%trace(:, ii) = buffer%trace(:, ii) + dg*qvec(jj)
         buffer%trace(:, jj) = buffer%trace(:, jj) - dg*qvec(ii)
         buffer%dr(:, iat, jj) = buffer%dr(:, iat, jj) + dg*qvec(ii)
         buffer%dr(:, jat, ii) = buffer%dr(:, jat, ii) - dg*qvec(jj)
      end if
      buffer%dL(:, :, jj) = buffer%dL(:, :, jj) + ds*qvec(ii)
      if (ii /= jj) buffer%dL(:, :, ii) = buffer%dL(:, :, ii) + ds*qvec(jj)
   end if
end subroutine add_derivative

!> Only a thread reduction is performed here; MPI sums the local results later.
subroutine reduce_derivatives(buffer, dadr, dadL, atrace, gradient, sigma)
   type(derivative_buffer), intent(inout) :: buffer
   real(wp), intent(inout), optional :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
   real(wp), intent(inout), optional :: gradient(:, :), sigma(:, :)

   !$omp critical (coulomb_derivatives)
   if (allocated(buffer%gradient)) then
      gradient = gradient + buffer%gradient
      sigma = sigma + buffer%sigma
   else
      dadr = dadr + buffer%dr
      dadL = dadL + buffer%dL
      atrace = atrace + buffer%trace
   end if
   !$omp end critical (coulomb_derivatives)
   if (allocated(buffer%gradient)) then
      deallocate(buffer%gradient)
   else
      deallocate(buffer%dr, buffer%dL, buffer%trace)
   end if
end subroutine reduce_derivatives

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
   call get_derivatives(self, mol, cache, qat, qsh, dadr, dadL, atrace)

end subroutine get_coulomb_derivs


!> Evaluate uncontracted derivatives of Coulomb matrix for finite system
subroutine get_damat_0d(mol, nshell, offset, hubbard, gexp, qvec, dadr, dadL, atrace, &
      & partition, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each atom
   integer, intent(in) :: nshell(:)
   !> Index offset for each shell
   integer, intent(in) :: offset(:)
   !> Hubbard parameter for each shell and species
   real(wp), intent(in) :: hubbard(:, :, :, :)
   !> Exponent of Coulomb kernel
   real(wp), intent(in) :: gexp
   !> Partial charge vector
   real(wp), intent(in) :: qvec(:)
   !> Derivative of interactions with respect to cartesian displacements
   real(wp), intent(out), optional :: dadr(:, :, :)
   !> Derivative of interactions with respect to strain deformations
   real(wp), intent(out), optional :: dadL(:, :, :)
   !> On-site derivatives with respect to cartesian displacements
   real(wp), intent(out), optional :: atrace(:, :)
   !> Share of the atom pairs evaluated here, absent selects the complete work
   type(work_partition), intent(in), optional :: partition
   real(wp), intent(inout), optional :: gradient(:, :), sigma(:, :)

   integer :: iat, jat, izp, jzp, ii, jj, ish, jsh, b
   real(wp) :: vec(3), r1, gam, dtmp, dG(3), dS(3, 3)
   type(derivative_buffer) :: buffer

   if (present(dadr)) then
      atrace = 0.0_wp
      dadr = 0.0_wp
      dadL = 0.0_wp
   end if

   !$omp parallel default(none) &
   !$omp shared(atrace, dadr, dadL, gradient, sigma, mol, qvec, hubbard, nshell, offset, gexp, partition) &
   !$omp private(iat, izp, ii, ish, jat, jzp, jj, jsh, gam, r1, vec, dG, dS, dtmp, b) &
   !$omp private(buffer)
   call new_derivative_buffer(buffer, mol%nat, size(qvec), present(gradient))
   !$omp do schedule(static, 1)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = offset(iat)
      do jat = 1, iat-1
         if (.not.owns_pair(partition, iat, jat)) cycle
         jzp = mol%id(jat)
         jj = offset(jat)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         r1 = norm2(vec)
         do ish = 1, nshell(iat)
            do jsh = 1, nshell(jat)
               gam = hubbard(jsh, ish, jzp, izp)
               dtmp = 1.0_wp / (r1**gexp + gam**(-gexp))
               dtmp = -r1**(gExp-2.0_wp) * dtmp * dtmp**(1.0_wp/gExp)
               dG = dtmp*vec
               do b = 1, 3
                  dS(:, b) = dG*vec(b)
               end do
               call add_derivative(buffer, iat, jat, ii+ish, jj+jsh, qvec, dG, dS)
            end do
         end do
      end do
   end do
   call reduce_derivatives(buffer, dadr, dadL, atrace, gradient, sigma)
   !$omp end parallel

end subroutine get_damat_0d

!> Evaluate uncontracted derivatives of Coulomb matrix for 3D periodic system
subroutine get_damat_3d(mol, nshell, offset, hubbard, gexp, wsc, alpha, ewald, qvec, &
      & dadr, dadL, atrace, partition, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   type(ewald_cache), intent(in) :: ewald
   !> Number of shells for each atom
   integer, intent(in) :: nshell(:)
   !> Index offset for each shell
   integer, intent(in) :: offset(:)
   !> Hubbard parameter for each shell and species
   real(wp), intent(in) :: hubbard(:, :, :, :)
   !> Exponent of Coulomb kernel
   real(wp), intent(in) :: gexp
   !> Wigner-Seitz image information
   type(wignerseitz_cell), intent(in) :: wsc
   !> Convergence factor for Ewald sum
   real(wp), intent(in) :: alpha
   !> Partial charge vector
   real(wp), intent(in) :: qvec(:)
   !> Derivative of interactions with respect to cartesian displacements
   real(wp), intent(out), optional :: dadr(:, :, :)
   !> Derivative of interactions with respect to strain deformations
   real(wp), intent(out), optional :: dadL(:, :, :)
   !> On-site derivatives with respect to cartesian displacements
   real(wp), intent(out), optional :: atrace(:, :)
   !> Share of the atom pairs evaluated here, absent selects the complete work
   type(work_partition), intent(in), optional :: partition
   real(wp), intent(inout), optional :: gradient(:, :), sigma(:, :)

   integer :: iat, jat, izp, jzp, img, ii, jj, ish, jsh
   logical :: need_weight_energy
   real(wp) :: gam, stmp, vec(3), dG(3), dS(3, 3)
   real(wp) :: dGd(3), dSd(3, 3), dGr(3), dSr(3, 3), dGw(3), dSw(3, 3)
   real(wp) :: weight(size(wsc%tridx, 1))
   real(wp) :: dwdr(3, size(wsc%tridx, 1)), dwdL(3, 3, size(wsc%tridx, 1))
   type(derivative_buffer) :: buffer
   real(wp), allocatable :: dtrans(:, :)

   if (present(dadr)) then
      atrace = 0.0_wp
      dadr = 0.0_wp
      dadL = 0.0_wp
   end if

   if (abs(gexp - 2.0_wp) < eps) then
      call get_damat_ko_3d(mol, nshell, offset, hubbard, alpha, ewald, qvec, dadr, &
         & dadL, atrace, partition, gradient, sigma)
      return
   end if

   call get_dir_trans(mol%lattice, alpha, conv, dtrans)

   !$omp parallel default(none) shared(atrace, dadr, dadL, gradient, sigma) &
   !$omp shared(mol, wsc, alpha, dtrans, ewald, qvec, hubbard, nshell, offset, gexp) &
   !$omp shared(partition) &
   !$omp private(iat, izp, jat, jzp, img, ii, jj, ish, jsh, gam, stmp, need_weight_energy) &
   !$omp private(weight, dwdr, dwdL) &
   !$omp private(vec, dG, dS, dGr, dSr, dGd, dSd, dGw, dSw, buffer)
   call new_derivative_buffer(buffer, mol%nat, size(qvec), present(gradient))
   !$omp do schedule(static, 1)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = offset(iat)
      do jat = 1, iat-1
         if (.not.owns_pair(partition, iat, jat)) cycle
         jzp = mol%id(jat)
         jj = offset(jat)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         call get_wignerseitz_weights(wsc, jat, iat, vec, weight, dwdr, dwdL)
         need_weight_energy = any(dwdr(:, :wsc%nimg(jat, iat)) /= 0.0_wp) &
            & .or. any(dwdL(:, :, :wsc%nimg(jat, iat)) /= 0.0_wp)
         call get_damat_dir_3d(vec, alpha, dtrans, dGd, dSd)
         call get_damat_rec_3d(vec, ewald, dGr, dSr)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
            do ish = 1, nshell(iat)
               do jsh = 1, nshell(jat)
                  gam = hubbard(jsh, ish, jzp, izp)
                  stmp = 0.0_wp
                  if (need_weight_energy) call get_amat_wsc_3d(vec, gam, gexp, stmp)
                  call get_damat_wsc_3d(vec, gam, gexp, dGw, dSw)
                  dG = (dGd + dGr + dGw) * weight(img) + stmp*dwdr(:, img)
                  dS = (dSd + dSr + dSw) * weight(img) + stmp*dwdL(:, :, img)
                  call add_derivative(buffer, iat, jat, ii+ish, jj+jsh, qvec, dG, dS)
               end do
            end do
         end do
      end do

      if (.not.owns_pair(partition, iat, iat)) cycle
      vec = 0.0_wp
      dG = 0.0_wp
      call get_wignerseitz_weights(wsc, iat, iat, vec, weight, dwdr, dwdL)
      need_weight_energy = any(dwdL(:, :, :wsc%nimg(iat, iat)) /= 0.0_wp)
      call get_damat_dir_3d(vec, alpha, dtrans, dGd, dSd)
      call get_damat_rec_3d(vec, ewald, dGr, dSr)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         do ish = 1, nshell(iat)
            do jsh = 1, ish-1
               gam = hubbard(jsh, ish, izp, izp)
               stmp = 0.0_wp
               if (need_weight_energy) call get_amat_wsc_3d(vec, gam, gexp, stmp)
               call get_damat_wsc_3d(vec, gam, gexp, dGw, dSw)
               dS = (dSd + dSr + dSw) * weight(img) + stmp*dwdL(:, :, img)
               call add_derivative(buffer, iat, iat, ii+ish, ii+jsh, qvec, dG, dS)
            end do
            gam = hubbard(ish, ish, izp, izp)
            stmp = 0.0_wp
            if (need_weight_energy) call get_amat_wsc_3d(vec, gam, gexp, stmp)
            call get_damat_wsc_3d(vec, gam, gexp, dGw, dSw)
            dS = (dSd + dSr + dSw) * weight(img) + stmp*dwdL(:, :, img)
            call add_derivative(buffer, iat, iat, ii+ish, ii+ish, qvec, dG, dS)
         end do
      end do
   end do
   call reduce_derivatives(buffer, dadr, dadL, atrace, gradient, sigma)
   !$omp end parallel

end subroutine get_damat_3d

!> Derivatives of the periodic Klopman-Ohno generalized Ewald matrix
subroutine get_damat_ko_3d(mol, nshell, offset, hubbard, alpha, ewald, qvec, dadr, &
      & dadL, atrace, partition, gradient, sigma)
   type(structure_type), intent(in) :: mol
   type(ewald_cache), intent(in) :: ewald
   integer, intent(in) :: nshell(:), offset(:)
   real(wp), intent(in) :: hubbard(:, :, :, :), alpha, qvec(:)
   real(wp), intent(out), optional :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
   type(work_partition), intent(in), optional :: partition
   real(wp), intent(inout), optional :: gradient(:, :), sigma(:, :)

   integer :: iat, jat, izp, jzp, ii, jj, ish, jsh
   real(wp) :: vol, gam, vec(3), dG(3), dS(3, 3)
   real(wp) :: dG1(3), dS1(3, 3), dG3(3), dS3(3, 3), dGsr(3), dSsr(3, 3)
   type(derivative_buffer) :: buffer
   type(ko_distances) :: distances
   real(wp), allocatable :: dtrans(:, :), strans(:, :)

   if (present(dadr)) then
      atrace = 0.0_wp
      dadr = 0.0_wp
      dadL = 0.0_wp
   end if
   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, alpha, conv, dtrans)
   call get_lattice_points([.true.], mol%lattice, ko_cutoff, strans)

   !$omp parallel default(none) shared(atrace, dadr, dadL, gradient, sigma) &
   !$omp shared(mol, nshell, offset, hubbard, alpha, vol, qvec, dtrans, ewald, strans, partition) &
   !$omp private(iat, izp, ii, ish, jat, jzp, jj, jsh, gam, vec, dG, dS) &
   !$omp private(dG1, dS1, dG3, dS3, dGsr, dSsr, buffer, distances)
   call new_derivative_buffer(buffer, mol%nat, size(qvec), present(gradient))
   !$omp do schedule(static, 1)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = offset(iat)
      do jat = 1, iat-1
         if (.not.owns_pair(partition, iat, jat)) cycle
         jzp = mol%id(jat)
         jj = offset(jat)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         call get_ko_distances(vec, strans, distances)
         call get_damat_ko_ewald_3d(vec, alpha, vol, dtrans, ewald, distances, &
            & dG1, dS1, dG3, dS3)
         do ish = 1, nshell(iat)
            do jsh = 1, nshell(jat)
               gam = hubbard(jsh, ish, jzp, izp)
               call get_damat_ko_short_3d(gam, distances, dGsr, dSsr)
               dG = dG1 + dGsr - 0.5_wp*dG3/(gam*gam)
               dS = dS1 + dSsr - 0.5_wp*dS3/(gam*gam)
               call add_derivative(buffer, iat, jat, ii+ish, jj+jsh, qvec, dG, dS)
            end do
         end do
      end do

      if (.not.owns_pair(partition, iat, iat)) cycle
      vec = 0.0_wp
      dG = 0.0_wp
      call get_ko_distances(vec, strans, distances)
      call get_damat_ko_ewald_3d(vec, alpha, vol, dtrans, ewald, distances, &
         & dG1, dS1, dG3, dS3)
      do ish = 1, nshell(iat)
         do jsh = 1, ish-1
            gam = hubbard(jsh, ish, izp, izp)
            call get_damat_ko_short_3d(gam, distances, dGsr, dSsr)
            dS = dS1 + dSsr - 0.5_wp*dS3/(gam*gam)
            call add_derivative(buffer, iat, iat, ii+ish, ii+jsh, qvec, dG, dS)
         end do
         gam = hubbard(ish, ish, izp, izp)
         call get_damat_ko_short_3d(gam, distances, dGsr, dSsr)
         dS = dS1 + dSsr - 0.5_wp*dS3/(gam*gam)
         call add_derivative(buffer, iat, iat, ii+ish, ii+ish, qvec, dG, dS)
      end do
   end do
   call reduce_derivatives(buffer, dadr, dadL, atrace, gradient, sigma)
   if (allocated(distances%r2)) then
      deallocate(distances%vec, distances%r2, distances%invr, distances%invr3, distances%invr5)
   end if
   !$omp end parallel
end subroutine get_damat_ko_3d

!> Derivatives of the shell-independent generalized Ewald sums
subroutine get_damat_ko_ewald_3d(rij, alpha, vol, dtrans, ewald, distances, &
      & dg1, ds1, dg3, ds3)
   real(wp), intent(in) :: rij(3), alpha, vol
   real(wp), intent(in) :: dtrans(:, :)
   type(ko_distances), intent(in) :: distances
   type(ewald_cache), intent(in) :: ewald
   real(wp), intent(out) :: dg1(3), ds1(3, 3), dg3(3), ds3(3, 3)

   integer :: itr, b
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
      do b = 1, 3
         ds3(:, b) = ds3(:, b) + dtmp*vec*vec(b)
      end do
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
   real(wp), intent(in) :: gam
   type(ko_distances), intent(in) :: distances
   real(wp), intent(out) :: dg(3), ds(3, 3)
   integer :: itr, b
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
      do b = 1, 3
         ds(:, b) = ds(:, b) + dtmp*vec*vec(b)
      end do
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

   integer :: itr, b
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
      do b = 1, 3
         ds(:, b) = ds(:, b) + dtmp*vec*vec(b)
      end do
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

   integer :: b
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
   do b = 1, 3
      ds(:, b) = dtmp*rij*rij(b)
   end do

end subroutine get_damat_wsc_3d

!> Calculate reciprocal space contributions for a pair under 3D periodic boundary conditions
subroutine get_damat_rec_3d(rij, ewald, dg, ds)
   real(wp), intent(in) :: rij(3)
   type(ewald_cache), intent(in) :: ewald
   real(wp), intent(out) :: dg(3), ds(3, 3)
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
