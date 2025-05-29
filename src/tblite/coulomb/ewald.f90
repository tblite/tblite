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

!> @file tblite/coulomb/ewald.f90
!> Provides an utilities for implementing Ewald summations

!> Helper tools for dealing with Ewald summation related calculations
module tblite_coulomb_ewald
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   implicit none
   private

   public :: get_alpha, get_dir_cutoff, get_rec_cutoff
   public :: get_amat_ewald_3d, get_amat_ewald_mp_3d
   public :: get_energy_ewald_3d, get_energy_ewald_mp_3d, get_potential_ewald_3d

   real(wp), parameter :: twopi = 2 * pi
   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

   !> Evaluator for interaction term
   type, abstract :: term_type
   contains
      !> Interaction at specified distance
      procedure(get_value), deferred :: get_value
   end type term_type

   abstract interface
      !> Interaction at specified distance
      pure function get_value(self, dist, alpha, vol) result(val)
         import :: term_type, wp
         !> Instance of interaction
         class(term_type), intent(in) :: self
         !> Distance between the two atoms
         real(wp), intent(in) :: dist
         !> Parameter of the Ewald summation
         real(wp), intent(in) :: alpha
         !> Volume of the real space unit cell
         real(wp), intent(in) :: vol
         !> Value of the interaction
         real(wp) :: val
      end function get_value
   end interface

   !> Real space Coulombic interaction 1/R
   type, extends(term_type) :: dir_term
   contains
      procedure :: get_value => dir_value
   end type dir_term

   !> Reciprocal space Coulombic interaction 1/R in 3D
   type, extends(term_type) :: rec_3d_term
   contains
      procedure :: get_value => rec_3d_value
   end type rec_3d_term

   !> Real space Coulombic interaction 1/R^2
   type, extends(term_type) :: dir_mp_term
   contains
      procedure :: get_value => dir_mp_value
   end type dir_mp_term

   !> Reciprocal space Coulombic interaction 1/R^2 in 3D
   type, extends(term_type) :: rec_3d_mp_term
   contains
      procedure :: get_value => rec_3d_mp_value
   end type rec_3d_mp_term

contains


!> Convenience interface to determine Ewald splitting parameter
subroutine get_alpha(lattice, alpha, multipole)
   !> Lattice vectors
   real(wp), intent(in) :: lattice(:, :)
   !> Estimated Ewald splitting parameter
   real(wp), intent(out) :: alpha
   !> Multipole expansion is used
   logical, intent(in) :: multipole

   real(wp) :: vol, rec_lat(3, 3)
   class(term_type), allocatable :: dirv, recv

   vol = abs(matdet_3x3(lattice))
   rec_lat = twopi*transpose(matinv_3x3(lattice))
   if (multipole) then
      dirv = dir_mp_term()
      recv = rec_3d_mp_term()
   else
      dirv = dir_term()
      recv = rec_3d_term()
   end if

   call search_alpha(lattice, rec_lat, vol, eps, dirv, recv, alpha)
end subroutine get_alpha


!> Get optimal alpha-parameter for the Ewald summation by finding alpha, where
!> decline of real and reciprocal part of Ewald are equal.
subroutine search_alpha(lattice, rec_lat, volume, tolerance, dirv, recv, alpha)
   !> Lattice vectors
   real(wp), intent(in) :: lattice(:,:)
   !> Reciprocal vectors
   real(wp), intent(in) :: rec_lat(:,:)
   !> Volume of the unit cell
   real(wp), intent(in) :: volume
   !> Tolerance for difference in real and rec. part
   real(wp), intent(in) :: tolerance
   !> Real-space interaction term
   class(term_type), intent(in) :: dirv
   !> Reciprocal space interaction term
   class(term_type), intent(in) :: recv
   !> Optimal alpha
   real(wp), intent(out) :: alpha

   real(wp) :: alpl, alpr, rlen, dlen, diff
   real(wp), parameter :: alpha0 = sqrt(epsilon(0.0_wp))
   integer, parameter :: niter = 30
   integer :: ibs, stat

   rlen = sqrt(minval(sum(rec_lat(:,:)**2, dim=1)))
   dlen = sqrt(minval(sum(lattice(:,:)**2, dim=1)))

   stat = 0
   alpha = alpha0
   diff = rec_dir_diff(alpha, dirv, recv, rlen, dlen, volume)
   do while (diff < -tolerance .and. alpha <= huge(1.0_wp))
      alpha = 2.0_wp * alpha
      diff = rec_dir_diff(alpha, dirv, recv, rlen, dlen, volume)
   end do
   if (alpha > huge(1.0_wp)) then
      stat = 1
   elseif (alpha == alpha0) then
      stat = 2
   end if

   if (stat == 0) then
      alpl = 0.5_wp * alpha
      do while (diff < tolerance .and. alpha <= huge(1.0_wp))
         alpha = 2.0_wp * alpha
         diff = rec_dir_diff(alpha, dirv, recv, rlen, dlen, volume)
      end do
      if (alpha > huge(1.0_wp)) then
         stat = 3
      end if
   end if

   if (stat == 0) then
      alpr = alpha
      alpha = (alpl + alpr) * 0.5_wp
      ibs = 0
      diff = rec_dir_diff(alpha, dirv, recv, rlen, dlen, volume)
      do while (abs(diff) > tolerance .and. ibs <= niter)
         if (diff < 0) then
            alpl = alpha
         else
            alpr = alpha
         end if
         alpha = (alpl + alpr) * 0.5_wp
         diff = rec_dir_diff(alpha, dirv, recv, rlen, dlen, volume)
         ibs = ibs + 1
      end do
      if (ibs > niter) then
         stat = 4
      end if
   end if

   if (stat /= 0) then
      alpha = 0.25_wp
   end if

end subroutine search_alpha


!> Return cutoff for reciprocal contributions
function get_rec_cutoff(alpha, volume, conv) result(x)
   !> Parameter of Ewald summation
   real(wp), intent(in) :: alpha
   !> Volume of the unit cell
   real(wp), intent(in) :: volume
   !> Tolerance value
   real(wp), intent(in) :: conv
   !> Magnitude of reciprocal vector
   real(wp) :: x

   class(term_type), allocatable :: term

   term = rec_3d_term()
   x = search_cutoff(term, alpha, volume, conv)

end function get_rec_cutoff


!> Return cutoff for real-space contributions
function get_dir_cutoff(alpha, conv) result(x)
   !> Parameter of Ewald summation
   real(wp), intent(in) :: alpha
   !> Tolerance for Ewald summation
   real(wp), intent(in) :: conv
   !> Magnitude of real-space vector
   real(wp) :: x

   real(wp), parameter :: volume = 0.0_wp
   class(term_type), allocatable :: term

   term = dir_term()
   x = search_cutoff(term, alpha, volume, conv)

end function get_dir_cutoff


!> Search for cutoff value of interaction term
function search_cutoff(term, alpha, volume, conv) result(x)
   !> Interaction term
   class(term_type), intent(in) :: term
   !> Parameter of Ewald summation
   real(wp), intent(in) :: alpha
   !> Volume of the unit cell
   real(wp), intent(in) :: volume
   !> Tolerance value
   real(wp), intent(in) :: conv
   !> Magnitude of reciprocal vector
   real(wp) :: x

   real(wp), parameter :: init = sqrt(epsilon(0.0_wp))
   integer, parameter :: miter = 30
   real(wp) :: xl, xr, yl, yr, y
   integer :: iter

   x = init
   y = term%get_value(x, alpha, volume)
   do while (y > conv .and. x <= huge(1.0_wp))
      x = 2.0_wp * x
      y = term%get_value(x, alpha, volume)
   end do

   xl = 0.5_wp * x
   yl = term%get_value(xl, alpha, volume)
   xr = x
   yr = y

   do iter = 1, miter
      if (yl - yr <= conv) exit
      x = 0.5_wp * (xl + xr)
      y = term%get_value(x, alpha, volume)
      if (y >= conv) then
         xl = x
         yl = y
      else
         xr = x
         yr = y
      end if
   end do

end function search_cutoff


!> Returns the difference in the decrease of the real and reciprocal parts of the
!> Ewald sum. In order to make the real space part shorter than the reciprocal
!> space part, the values are taken at different distances for the real and the
!> reciprocal space parts.
pure function rec_dir_diff(alpha, dirv, recv, rlen, dlen, volume) result(diff)
   !> Parameter for the Ewald summation
   real(wp), intent(in) :: alpha
   !> Procedure pointer to real-space routine
   class(term_type), intent(in) :: dirv
   !> Procedure pointer to reciprocal routine
   class(term_type), intent(in) :: recv
   !> Length of the shortest reciprocal space vector in the sum
   real(wp), intent(in) :: rlen
   !> Length of the shortest real space vector in the sum
   real(wp), intent(in) :: dlen
   !> Volume of the real space unit cell
   real(wp), intent(in) :: volume
   !> Difference between changes in the two terms
   real(wp) :: diff

   diff = ((recv%get_value(4*rlen, alpha, volume) - recv%get_value(5*rlen, alpha, volume))) &
      & - (dirv%get_value(2*dlen, alpha, volume) - dirv%get_value(3*dlen, alpha, volume))

end function rec_dir_diff


!> Returns the max. value of a term in the reciprocal space part of the Ewald
!> summation for a given vector length.
pure function rec_3d_value(self, dist, alpha, vol) result(rval)
   !> Instance of interaction
   class(rec_3d_term), intent(in) :: self
   !> Length of the reciprocal space vector
   real(wp), intent(in) :: dist
   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha
   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol
   !> Reciprocal term
   real(wp) :: rval

   rval = 4.0_wp*pi*(exp(-0.25_wp*dist*dist/(alpha**2))/(vol*dist*dist))

end function rec_3d_value


!> Returns the max. value of a term in the reciprocal space part of the Ewald
!> summation for a given vector length.
pure function rec_3d_mp_value(self, dist, alpha, vol) result(rval)
   !> Instance of interaction
   class(rec_3d_mp_term), intent(in) :: self
   !> Length of the reciprocal space vector
   real(wp), intent(in) :: dist
   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha
   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol
   !> Reciprocal term
   real(wp) :: rval

   real(wp) :: g2

   g2 = dist*dist

   rval = 4.0_wp*pi*(exp(-0.25_wp*g2/(alpha**2))/vol)

end function rec_3d_mp_value


!> Direct space interaction at specified distance
pure function dir_value(self, dist, alpha, vol) result(val)
   !> Instance of interaction
   class(dir_term), intent(in) :: self
   !> Distance between the two atoms
   real(wp), intent(in) :: dist
   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha
   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol
   !> Value of the interaction
   real(wp) :: val

   val = erfc(alpha*dist)/dist
end function dir_value

!> Direct space interaction at specified distance
pure function dir_mp_value(self, dist, alpha, vol) result(val)
   !> Instance of interaction
   class(dir_mp_term), intent(in) :: self
   !> Distance between the two atoms
   real(wp), intent(in) :: dist
   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha
   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol
   !> Value of the interaction
   real(wp) :: val

   real(wp) :: arg

   arg = alpha * dist

   val = (2/sqrtpi*arg*exp(-arg*arg) + erfc(arg))/dist**3
end function dir_mp_value

!> Calculate the Coulomb matrix for a 3D periodic system using Ewald summation
subroutine get_amat_ewald_3d(mol, nshell, offset, alpha, vol, rtrans, amat)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells per atom
   integer, intent(in) :: nshell(:)
   !> Index offset for each atom
   integer, intent(in) :: offset(:)
   !> Convergence factor
   real(wp), intent(in) :: alpha
   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol
   !> Reciprocal lattice translations
   real(wp), intent(in) :: rtrans(:, :)
   !> Coulomb matrix
   real(wp), intent(inout) :: amat(:, :)

   integer :: iat, jat, izp, jzp, img, ii, jj, ish, jsh
   real(wp) :: vec(3), dtmp, rtmp, aval

   !$omp parallel do default(none) schedule(runtime) shared(amat) &
   !$omp shared(mol, nshell, offset, rtrans, alpha, vol) &
   !$omp private(iat, izp, jat, jzp, ii, jj, ish, jsh, vec, aval)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = offset(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         jj = offset(jat)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         call get_amat_rec_3d(vec, vol, alpha, rtrans, aval)
         aval = 2 * aval
         do ish = 1, nshell(iat)
            do jsh = 1, nshell(jat)
               !$omp atomic
               amat(jj+jsh, ii+ish) = amat(jj+jsh, ii+ish) + aval
               !$omp atomic
               amat(ii+ish, jj+jsh) = amat(ii+ish, jj+jsh) + aval
            end do
         end do
      end do

      vec = 0.0_wp
      call get_amat_rec_3d(vec, vol, alpha, rtrans, aval)
      aval = 2 * (aval - alpha / sqrtpi)
      do ish = 1, nshell(iat)
         do jsh = 1, ish-1
            !$omp atomic
            amat(ii+jsh, ii+ish) = amat(ii+jsh, ii+ish) + aval
            !$omp atomic
            amat(ii+ish, ii+jsh) = amat(ii+ish, ii+jsh) + aval
         end do
         !$omp atomic
         amat(ii+ish, ii+ish) = amat(ii+ish, ii+ish) + aval
      end do

   end do
end subroutine get_amat_ewald_3d

subroutine get_amat_ewald_mp_3d(mol, alpha, vol, rtrans, amat_sd, amat_dd, amat_sq)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Convergence parameter for Ewald sum
   real(wp), intent(in) :: alpha
   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol
   !> Reciprocal lattice translations
   real(wp), intent(in) :: rtrans(:, :)
   !> Interation matrix for charges and dipoles
   real(wp), intent(inout) :: amat_sd(:, :, :)
   !> Interation matrix for dipoles and dipoles
   real(wp), intent(inout) :: amat_dd(:, :, :, :)
   !> Interation matrix for charges and quadrupoles
   real(wp), intent(inout) :: amat_sq(:, :, :)

   integer :: iat, jat, k
   real(wp) :: vec(3), rr
   real(wp) :: r_sd(3), r_dd(3, 3), r_sq(6)

   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp shared(amat_sd, amat_dd, amat_sq, mol, vol, alpha, rtrans) &
   !$omp private(iat, jat, vec, r_sd, r_dd, r_sq)
   do iat = 1, mol%nat
      do jat = 1, mol%nat
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)

         call get_amat_sdq_rec_3d(vec, vol, alpha, rtrans, r_sd, r_dd, r_sq)

         amat_sd(:, jat, iat) = amat_sd(:, jat, iat) + r_sd
         amat_dd(:, jat, :, iat) = amat_dd(:, jat, :, iat) + r_dd
         amat_sq(:, jat, iat) = amat_sq(:, jat, iat) + r_sq
      end do
   end do

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(amat_sd, amat_dd, amat_sq, mol, vol, alpha) private(iat, rr, k)
   do iat = 1, mol%nat
      ! dipole-dipole selfenergy: -2/3·α³/sqrt(π) Σ(i) μ²(i)
      rr = -4.0_wp/3.0_wp * alpha**3 / sqrtpi
      do k = 1, 3
         amat_dd(k, iat, k, iat) = amat_dd(k, iat, k, iat) + rr
      end do

      ! charge-quadrupole selfenergy: 4/9·α³/sqrt(π) Σ(i) q(i)Tr(θi)
      ! (no actual contribution since quadrupoles are traceless)
      rr = 4.0_wp/9.0_wp * alpha**3 / sqrtpi
      amat_sq([1, 3, 6], iat, iat) = amat_sq([1, 3, 6], iat, iat) + rr
   end do
end subroutine get_amat_ewald_mp_3d

!> Calculate energies for each atom in the system using Ewald summation
subroutine get_energy_ewald_3d(mol, alpha, vol, rtrans, qat, energies)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha
   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol
   !> Reciprocal lattice translations
   real(wp), intent(in) :: rtrans(:, :)
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Energies for each atom
   real(wp), intent(inout) :: energies(:)

   integer :: iat, ii, ish, iG
   real(wp) :: fac, g2, expk, qcos, qsin, gv

   fac = 4 * pi / vol
   do iG = 1, size(rtrans, 2)
      g2 = dot_product(rtrans(:, iG), rtrans(:, iG))
      if (g2 < eps) cycle
      expk = fac * exp(-0.25_wp*g2/(alpha*alpha))/g2
      call get_fourier_transform(mol, rtrans(:, iG), qat, qcos, qsin)
      do iat = 1, mol%nat
         gv = dot_product(rtrans(:, iG), mol%xyz(:, iat))
         energies(iat) = energies(iat) + expk * (qcos * cos(gv) + qsin * sin(gv)) * qat(iat)
      end do
   end do
   ! Self energy correction
   do iat = 1, mol%nat
      energies(iat) = energies(iat) - alpha / (sqrt(pi)) * qat(iat)**2
   end do
end subroutine get_energy_ewald_3d

!> Calculate energies for each atom in the system using Ewald summation
subroutine get_energy_ewald_mp_3d(mol, alpha, vol, rtrans, qat, dpat, qpat, energies)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha
   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol
   !> Reciprocal lattice translations
   real(wp), intent(in) :: rtrans(:, :)
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Atomic dipole moments
   real(wp), intent(in) :: dpat(:, :)
   !> Atomic quadrupole moments
   real(wp), intent(in) :: qpat(:, :)
   !> Energies for each atom
   real(wp), intent(inout) :: energies(:)

   integer :: iat, ii, ish, iG
   real(wp) :: fac, g2, expk, qcos, qsin, dpcos, dpsin, qpcos, qpsin, gv
   real(wp) :: dpg, qpg2, cosgv, singv, dp2, qptr, gvec(3)

   fac = 4 * pi / vol
   do iG = 1, size(rtrans, 2)
      gvec = rtrans(:, iG)
      g2 = dot_product(gvec, gvec)
      if (g2 < eps) cycle
      expk = fac * exp(-0.25_wp*g2/(alpha*alpha))/g2
      call get_mp_fourier_transform(mol, gvec, qat, dpat, qpat, &
         & qcos, qsin, dpcos, dpsin, qpcos, qpsin)
      do iat = 1, mol%nat
         gv = dot_product(gvec, mol%xyz(:, iat))
         dpg = dot_product(dpat(:, iat), gvec)
         qpg2 = (qpat(1, iat) * gvec(1)**2 + &
            & qpat(3, iat) * gvec(2)**2 + qpat(6, iat) * gvec(3)**2 + &
            & 2 * (qpat(2, iat) * gvec(1) * gvec(2) + &
            & qpat(4, iat) * gvec(1) * gvec(3) + qpat(5, iat) * gvec(2) * gvec(3)))/3
         cosgv = cos(gv)
         singv = sin(gv)
         energies(iat) = energies(iat) + &
            & ((singv * dpcos - cosgv * dpsin) * qat(iat) &
            & + (qsin * cosgv - qcos * singv) * dpg &
            & + (dpcos * cosgv + dpsin * singv) * dpg &
            & + (cosgv * qpcos + singv * qpsin) * qat(iat) &
            & + (qcos * cosgv + qsin * singv) * qpg2) * expk
      end do
   end do
   ! Self energy correction
   do iat = 1, mol%nat
      dp2 = dot_product(dpat(:, iat), dpat(:, iat))
      qptr = qpat(1, iat) + qpat(3, iat) + qpat(6, iat)
      energies(iat) = energies(iat) &
         & - 2 * alpha**3 / (3 * sqrtpi) * dp2 &
         & + 4 * alpha**3 / (9 * sqrtpi) * qptr * qat(iat)
   end do
end subroutine get_energy_ewald_mp_3d

!> Calculate potential for each atom in the system using Ewald summation
subroutine get_potential_ewald_3d(mol, alpha, vol, rtrans, qat, vat)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha
   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol
   !> Reciprocal lattice translations
   real(wp), intent(in) :: rtrans(:, :)
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Energies for each atom
   real(wp), intent(inout) :: vat(:)

   integer :: iat, ii, ish, iG
   real(wp) :: fac, g2, expk, qcos, qsin, gv

   fac = 4 * pi / vol
   do iG = 1, size(rtrans, 2)
      g2 = dot_product(rtrans(:, iG), rtrans(:, iG))
      if (g2 < eps) cycle
      expk = fac * exp(-0.25_wp*g2/(alpha*alpha))/g2
      call get_fourier_transform(mol, rtrans(:, iG), qat, qcos, qsin)
      do iat = 1, mol%nat
         gv = dot_product(rtrans(:, iG), mol%xyz(:, iat))
         ! energies(iat) = energies(iat) + expk * (qcos * cos(gv) + qsin * sin(gv)) * qat(iat)
         vat(iat) = vat(iat) + expk * (qcos * cos(gv) + qsin * sin(gv))
         vat(iat) = vat(iat) + expk * 2 * (cos(gv) * cos(gv) + sin(gv) * sin(gv)) * qat(iat)
      end do
   end do
   ! Self energy correction
   do iat = 1, mol%nat
      vat(iat) = vat(iat) - 2 * alpha / (sqrt(pi)) * qat(iat)
   end do
end subroutine get_potential_ewald_3d

!> Calculate energies for each atom in the system using Ewald summation
subroutine get_potential_ewald_mp_3d(mol, alpha, vol, rtrans, qat, dpat, qpat, vat, vdpat, vqpat)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha
   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol
   !> Reciprocal lattice translations
   real(wp), intent(in) :: rtrans(:, :)
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Atomic dipole moments
   real(wp), intent(in) :: dpat(:, :)
   !> Atomic quadrupole moments
   real(wp), intent(in) :: qpat(:, :)
   !> Potential for each atom
   real(wp), intent(inout) :: vat(:)
   !> Dipole potential for each atom
   real(wp), intent(inout) :: vdpat(:, :)
   !> Quadrupole potential for each atom
   real(wp), intent(inout) :: vqpat(:, :)

   integer :: iat, ii, ish, iG
   real(wp) :: fac, g2, expk, qcos, qsin, dpcos, dpsin, qpcos, qpsin, gv
   real(wp) :: gg(6), cosgv, singv, dp2, qptr, gvec(3)

   fac = 4 * pi / vol
   do iG = 1, size(rtrans, 2)
      gvec = rtrans(:, iG)
      g2 = dot_product(gvec, gvec)
      if (g2 < eps) cycle
      expk = fac * exp(-0.25_wp*g2/(alpha*alpha))/g2
      call get_mp_fourier_transform(mol, gvec, qat, dpat, qpat, &
         & qcos, qsin, dpcos, dpsin, qpcos, qpsin)
      do iat = 1, mol%nat
         gv = dot_product(gvec, mol%xyz(:, iat))
         gg(:) = [gvec(1)**2, 2 * gvec(1) * gvec(2), gvec(2)**2, &
            & 2 * gvec(1) * gvec(3), 2 * gvec(2) * gvec(3), gvec(3)**2]/3
         cosgv = cos(gv)
         singv = sin(gv)
         vat(iat) = vat(iat) + &
            & ((singv * dpcos - cosgv * dpsin) &
            & + (cosgv * qpcos + singv * qpsin)) * expk
         vdpat(:, iat) = vdpat(:, iat) + &
            & ((qsin * cosgv - qcos * singv) &
            & + (dpcos * cosgv + dpsin * singv)) * expk * gvec
         vqpat(:, iat) = vqpat(:, iat) + &
            & + (qcos * cosgv + qsin * singv) * gg * expk
      end do
   end do
   ! Self energy correction
   do iat = 1, mol%nat
      dp2 = dot_product(dpat(:, iat), dpat(:, iat))
      qptr = qpat(1, iat) + qpat(3, iat) + qpat(6, iat)
      vat(iat) = vat(iat) &
         & + 4 * alpha**3 / (9 * sqrtpi) * qptr
      vdpat(:, iat) = vdpat(:, iat) &
         & - 4 * alpha**3 / (3 * sqrtpi) * dpat(:, iat)
      vqpat([1, 3, 6], iat) = vqpat([1, 3, 6], iat) &
         & + 4 * alpha**3 / (9 * sqrtpi) * qat(iat)
   end do
end subroutine get_potential_ewald_mp_3d

!> Calculate reciprocal space contributions for a pair under 3D periodic boundary conditions
subroutine get_amat_rec_3d(rij, vol, alp, trans, amat)
   !> Distance between pair
   real(wp), intent(in) :: rij(3)
   !> Volume of cell
   real(wp), intent(in) :: vol
   !> Convergence factor
   real(wp), intent(in) :: alp
   !> Translation vectors to consider
   real(wp), intent(in) :: trans(:, :)
   !> Interaction matrix element
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: fac, vec(3), g2, gv, expk, cosk

   amat = 0.0_wp
   fac = 4*pi/vol

   do itr = 1, size(trans, 2)
      vec(:) = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      gv = dot_product(rij, vec)
      expk = fac * exp(-0.25_wp*g2/(alp*alp))/g2
      cosk = cos(gv) * expk
      amat = amat + cosk
   end do

end subroutine get_amat_rec_3d

pure subroutine get_amat_sdq_rec_3d(rij, vol, alp, trans, amat_sd, amat_dd, amat_sq)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: vol
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: amat_sd(:)
   real(wp), intent(out) :: amat_dd(:, :)
   real(wp), intent(out) :: amat_sq(:)

   integer :: itr
   real(wp) :: fac, vec(3), g2, expk, sink, cosk, gv

   amat_sd = 0.0_wp
   amat_dd = 0.0_wp
   amat_sq = 0.0_wp
   fac = 8*pi/vol

   do itr = 1, size(trans, 2)
      vec(:) = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      gv = dot_product(rij, vec)
      expk = fac*exp(-0.25_wp*g2/(alp*alp))/g2
      sink = sin(gv)*expk
      cosk = cos(gv)*expk

      amat_sd(:) = amat_sd + vec*sink
      amat_dd(:, :) = amat_dd + spread(vec, 1, 3) * spread(vec, 2, 3) * cosk
      amat_sq(1) = amat_sq(1) +   vec(1)*vec(1)*cosk/3
      amat_sq(2) = amat_sq(2) + 2*vec(1)*vec(2)*cosk/3
      amat_sq(3) = amat_sq(3) +   vec(2)*vec(2)*cosk/3
      amat_sq(4) = amat_sq(4) + 2*vec(1)*vec(3)*cosk/3
      amat_sq(5) = amat_sq(5) + 2*vec(2)*vec(3)*cosk/3
      amat_sq(6) = amat_sq(6) +   vec(3)*vec(3)*cosk/3
   end do

end subroutine get_amat_sdq_rec_3d

!> Calculate the Fourier transform of atomic partial charges
subroutine get_fourier_transform(mol, gvec, qat, qcos, qsin)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reciprocal lattice translations
   real(wp), intent(in) :: gvec(:)
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Cosine part of the Fourier transform
   real(wp), intent(out) :: qcos
   !> Sine part of the Fourier transform
   real(wp), intent(out) :: qsin

   integer :: iat
   real(wp) :: gv

   qcos = 0.0_wp
   qsin = 0.0_wp
   do iat = 1, mol%nat
      gv = dot_product(gvec, mol%xyz(:, iat))
      qcos = qcos + qat(iat) * cos(gv)
      qsin = qsin + qat(iat) * sin(gv)
   end do
end subroutine get_fourier_transform

!> Calculate the Fourier transform of multipole moments
subroutine get_mp_fourier_transform(mol, gvec, qat, dpat, qpat, &
      & qcos, qsin, dpcos, dpsin, qpcos, qpsin)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reciprocal lattice translations
   real(wp), intent(in) :: gvec(:)
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Atomic dipole moments
   real(wp), intent(in) :: dpat(:, :)
   !> Atomic quadrupole moments
   real(wp), intent(in) :: qpat(:, :)
   !> Cosine part of the Fourier transform
   real(wp), intent(out) :: qcos
   !> Sine part of the Fourier transform
   real(wp), intent(out) :: qsin
   !> Cosine part of the Fourier transform for dipoles
   real(wp), intent(out) :: dpcos
   !> Sine part of the Fourier transform for dipoles
   real(wp), intent(out) :: dpsin
   !> Cosine part of the Fourier transform for quadrupoles
   real(wp), intent(out) :: qpcos
   !> Sine part of the Fourier transform for quadrupoles
   real(wp), intent(out) :: qpsin

   integer :: iat
   real(wp) :: gv, dpg, qpg2, cosgv, singv

   qcos = 0.0_wp
   qsin = 0.0_wp
   dpcos = 0.0_wp
   dpsin = 0.0_wp
   qpcos = 0.0_wp
   qpsin = 0.0_wp
   do iat = 1, mol%nat
      gv = dot_product(gvec, mol%xyz(:, iat))
      dpg = dot_product(dpat(:, iat), gvec)
      qpg2 = (qpat(1, iat) * gvec(1)**2 + &
         & qpat(3, iat) * gvec(2)**2 + qpat(6, iat) * gvec(3)**2 + &
         & 2 * (qpat(2, iat) * gvec(1) * gvec(2) + &
         & qpat(4, iat) * gvec(1) * gvec(3) + qpat(5, iat) * gvec(2) * gvec(3)))
      cosgv = cos(gv)
      singv = sin(gv)
      qcos = qcos + qat(iat) * cosgv
      qsin = qsin + qat(iat) * singv
      dpcos = dpcos + dpg * cosgv
      dpsin = dpsin + dpg * singv
      qpcos = qpcos + qpg2 * cosgv / 3
      qpsin = qpsin + qpg2 * singv / 3
   end do
end subroutine get_mp_fourier_transform

end module tblite_coulomb_ewald
