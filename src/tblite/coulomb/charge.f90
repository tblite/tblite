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

!> Isotropic second-order electrostatics using an effective Coulomb operator
module tblite_coulomb_charge
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   use mctc_io_constants, only : pi
   use tblite_blas, only : dot, gemv, symv, gemm
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_type, only : coulomb_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_wignerseitz, only : wignerseitz_cell
   implicit none
   private

   public :: effective_coulomb, new_effective_coulomb

   public :: average_interface, harmonic_average, arithmetic_average, geometric_average


   type, extends(coulomb_type) :: effective_coulomb
      integer, allocatable :: nshell(:)
      integer, allocatable :: offset(:)
      real(wp), allocatable :: hardness(:, :, :, :)
      real(wp) :: gexp
   contains
      procedure :: update
      procedure :: variable_info
      procedure :: get_energy
      procedure :: get_potential
      procedure :: get_gradient
      procedure :: get_coulomb_matrix
      procedure :: get_coulomb_derivs
   end type effective_coulomb


   abstract interface
      pure function average_interface(gi, gj) result(gij)
         import :: wp

         !> Hardness of shell i
         real(wp), intent(in) :: gi

         !> Hardness of shell j
         real(wp), intent(in) :: gj

         !> Averaged hardness
         real(wp) :: gij

      end function average_interface
   end interface

   real(wp), parameter :: twopi = 2 * pi
   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

contains

subroutine new_effective_coulomb(self, mol, gexp, hardness, average, nshell)
   !> Instance of the electrostatic container
   type(effective_coulomb), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   real(wp), intent(in) :: gexp
   procedure(average_interface) :: average
   real(wp), intent(in) :: hardness(:, :)
   integer, intent(in), optional :: nshell(:)

   integer :: mshell
   integer :: isp, jsp, ish, jsh, ind, iat

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
      allocate(self%hardness(mshell, mshell, mol%nid, mol%nid))
      do isp = 1, mol%nid
         do jsp = 1, mol%nid
            self%hardness(:, :, jsp, isp) = 0.0_wp
            do ish = 1, nshell(isp)
               do jsh = 1, nshell(jsp)
                  self%hardness(jsh, ish, jsp, isp) = &
                     & average(hardness(ish, isp), hardness(jsh, jsp))
               end do
            end do
         end do
      end do
   else
      allocate(self%hardness(1, 1, mol%nid, mol%nid))
      do isp = 1, mol%nid
         do jsp = 1, mol%nid
            self%hardness(1, 1, jsp, isp) = average(hardness(1, isp), hardness(1, jsp))
         end do
      end do
   end if

end subroutine new_effective_coulomb


!> Harmonic averaging functions for hardnesses in GFN1-xTB
pure function harmonic_average(gi, gj) result(gij)
   !> Hardness of shell i
   real(wp), intent(in) :: gi
   !> Hardness of shell j
   real(wp), intent(in) :: gj
   !> Averaged hardness
   real(wp) :: gij

   gij = 2.0_wp/(1.0_wp/gi+1.0_wp/gj)

end function harmonic_average


!> Arithmetic averaging functions for hardnesses in GFN2-xTB
pure function arithmetic_average(gi, gj) result(gij)
   !> Hardness of shell i
   real(wp), intent(in) :: gi
   !> Hardness of shell j
   real(wp), intent(in) :: gj
   !> Averaged hardness
   real(wp) :: gij

   gij = 0.5_wp*(gi+gj)

end function arithmetic_average


!> Geometric averaging functions for hardnesses
pure function geometric_average(gi, gj) result(gij)
   !> Hardness of shell i
   real(wp), intent(in) :: gi
   !> Hardness of shell j
   real(wp), intent(in) :: gj
   !> Averaged hardness
   real(wp) :: gij

   gij = sqrt(gi*gj)

end function geometric_average


subroutine update(self, mol, cache)
   !> Instance of the electrostatic container
   class(effective_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache

   call cache%update(mol)

   if (.not.allocated(cache%amat)) then
      allocate(cache%amat(sum(self%nshell), sum(self%nshell)))
   end if
   call self%get_coulomb_matrix(mol, cache, cache%amat)

   if (.not.allocated(cache%vvec)) then
      allocate(cache%vvec(sum(self%nshell)))
   end if

end subroutine update


subroutine get_energy(self, mol, cache, wfn, energy)
   !> Instance of the electrostatic container
   class(effective_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Electrostatic energy
   real(wp), intent(inout) :: energy

   call symv(cache%amat, wfn%qsh, cache%vvec, alpha=0.5_wp)
   energy = energy + dot(cache%vvec, wfn%qsh)

end subroutine get_energy


subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the electrostatic container
   class(effective_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   call symv(cache%amat, wfn%qsh, pot%vsh, beta=1.0_wp)

end subroutine get_potential


subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the electrostatic container
   class(effective_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular gradient of the repulsion energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the repulsion energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   integer :: ndim
   real(wp), allocatable :: dadr(:, :, :), dadL(:, :, :), datr(:, :)

   ndim = sum(self%nshell)
   allocate(dadr(3, mol%nat, ndim), dadL(3, 3, ndim), datr(3, ndim))

   call self%get_coulomb_derivs(mol, cache, wfn%qat, wfn%qsh, dadr, dadL, datr)

   call gemv(dadr, wfn%qsh, gradient, beta=1.0_wp)
   call gemv(dadL, wfn%qsh, sigma, beta=1.0_wp, alpha=0.5_wp)

end subroutine get_gradient


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
      call get_amat_3d(mol, self%nshell, self%offset, self%hardness, self%gexp, &
         & cache%wsc, cache%alpha, amat)
   else
      call get_amat_0d(mol, self%nshell, self%offset, self%hardness, self%gexp, amat)
   end if

end subroutine get_coulomb_matrix


subroutine get_dir_trans(lattice, trans)
   real(wp), intent(in) :: lattice(:, :)
   real(wp), allocatable, intent(out) :: trans(:, :)
   integer, parameter :: rep(3) = 2

   call get_lattice_points(lattice, rep, .true., trans)

end subroutine get_dir_trans

subroutine get_rec_trans(lattice, trans)
   real(wp), intent(in) :: lattice(:, :)
   real(wp), allocatable, intent(out) :: trans(:, :)
   integer, parameter :: rep(3) = 2
   real(wp) :: rec_lat(3, 3)

   rec_lat = twopi*transpose(matinv_3x3(lattice))
   call get_lattice_points(rec_lat, rep, .false., trans)

end subroutine get_rec_trans


subroutine get_amat_0d(mol, nshell, offset, hardness, gexp, amat)
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: nshell(:)
   integer, intent(in) :: offset(:)
   real(wp), intent(in) :: hardness(:, :, :, :)
   real(wp), intent(in) :: gexp
   real(wp), intent(inout) :: amat(:, :)

   integer :: iat, jat, izp, jzp, ii, jj, ish, jsh
   real(wp) :: vec(3), r1, r1g, gam, tmp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:amat) shared(mol, nshell, offset, hardness, gexp) &
   !$omp private(iat, izp, ii, ish, jat, jzp, jj, jsh, gam, vec, r1, r1g, tmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = offset(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         jj = offset(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         r1g = r1**gexp
         do ish = 1, nshell(iat)
            do jsh = 1, nshell(jat)
               gam = hardness(jsh, ish, jzp, izp)
               tmp = 1.0_wp/(r1g + gam**(-gexp))**(1.0_wp/gexp)
               amat(jj+jsh, ii+ish) = amat(jj+jsh, ii+ish) + tmp
               amat(ii+ish, jj+jsh) = amat(ii+ish, jj+jsh) + tmp
            end do
         end do
      end do
      do ish = 1, nshell(iat)
         do jsh = 1, ish-1
            gam = hardness(jsh, ish, izp, izp)
            amat(ii+jsh, ii+ish) = amat(ii+jsh, ii+ish) + gam
            amat(ii+ish, ii+jsh) = amat(ii+ish, ii+jsh) + gam
         end do
         amat(ii+ish, ii+ish) = amat(ii+ish, ii+ish) + hardness(ish, ish, izp, izp)
      end do
   end do

end subroutine get_amat_0d

subroutine get_amat_3d(mol, nshell, offset, hardness, gexp, wsc, alpha, amat)
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: nshell(:)
   integer, intent(in) :: offset(:)
   real(wp), intent(in) :: hardness(:, :, :, :)
   real(wp), intent(in) :: gexp
   type(wignerseitz_cell), intent(in) :: wsc
   real(wp), intent(in) :: alpha
   real(wp), intent(inout) :: amat(:, :)

   integer :: iat, jat, izp, jzp, img, ii, jj, ish, jsh
   real(wp) :: vec(3), gam, wsw, dtmp, rtmp, vol
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, dtrans)
   call get_rec_trans(mol%lattice, rtrans)

   !$omp parallel do default(none) schedule(runtime) reduction(+:amat) &
   !$omp shared(mol, nshell, offset, hardness, gexp, wsc, dtrans, rtrans, alpha, vol) &
   !$omp private(iat, izp, jat, jzp, ii, jj, ish, jsh, gam, wsw, vec, dtmp, rtmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = offset(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         jj = offset(jat)
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
            call get_amat_rec_3d(vec, vol, alpha, rtrans, rtmp)
            do ish = 1, nshell(iat)
               do jsh = 1, nshell(jat)
                  gam = hardness(jsh, ish, jzp, izp)
                  call get_amat_dir_3d(vec, gam, gexp, alpha, dtrans, dtmp)
                  amat(jj+jsh, ii+ish) = amat(jj+jsh, ii+ish) + (dtmp + rtmp) * wsw
                  amat(ii+ish, jj+jsh) = amat(ii+ish, jj+jsh) + (dtmp + rtmp) * wsw
               end do
            end do
         end do
      end do

      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         call get_amat_rec_3d(vec, vol, alpha, rtrans, rtmp)
         rtmp = rtmp - 2 * alpha / sqrtpi
         do ish = 1, nshell(iat)
            do jsh = 1, ish-1
               gam = hardness(jsh, ish, izp, izp)
               call get_amat_dir_3d(vec, gam, gexp, alpha, dtrans, dtmp)
               amat(ii+jsh, ii+ish) = amat(ii+jsh, ii+ish) + (dtmp + rtmp + gam) * wsw
               amat(ii+ish, ii+jsh) = amat(ii+ish, ii+jsh) + (dtmp + rtmp + gam) * wsw
            end do
            gam = hardness(ish, ish, izp, izp)
            call get_amat_dir_3d(vec, gam, gexp, alpha, dtrans, dtmp)
            amat(ii+ish, ii+ish) = amat(ii+ish, ii+ish) + (dtmp + rtmp + gam) * wsw
         end do
      end do

   end do

end subroutine get_amat_3d

subroutine get_amat_dir_3d(rij, gam, gexp, alp, trans, amat)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: gam
   real(wp), intent(in) :: gexp
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: vec(3), r1, tmp

   amat = 0.0_wp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      tmp = 1.0_wp/(r1**gexp + gam**(-gexp))**(1.0_wp/gexp) - erf(alp*r1)/r1
      amat = amat + tmp
   end do

end subroutine get_amat_dir_3d

subroutine get_amat_rec_3d(rij, vol, alp, trans, amat)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: vol
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: fac, vec(3), g2, tmp, gv, expk, cosk

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


subroutine get_coulomb_derivs(self, mol, cache, qat, qsh, dadr, dadL, atrace)

   !> Instance of the electrostatic container
   class(effective_coulomb), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache

   real(wp), intent(in) :: qat(:)
   real(wp), intent(in) :: qsh(:)
   real(wp), contiguous, intent(out) :: dadr(:, :, :)
   real(wp), contiguous, intent(out) :: dadL(:, :, :)
   real(wp), contiguous, intent(out) :: atrace(:, :)

   if (any(mol%periodic)) then
      call get_damat_3d(mol, self%nshell, self%offset, self%hardness, self%gexp, &
         & cache%wsc, cache%alpha, qsh, dadr, dadL, atrace)
   else
      call get_damat_0d(mol, self%nshell, self%offset, self%hardness, self%gexp, qsh, &
         & dadr, dadL, atrace)
   end if

end subroutine get_coulomb_derivs


subroutine get_damat_0d(mol, nshell, offset, hardness, gexp, qvec, dadr, dadL, atrace)
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: nshell(:)
   integer, intent(in) :: offset(:)
   real(wp), intent(in) :: hardness(:, :, :, :)
   real(wp), intent(in) :: gexp
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(out) :: dadr(:, :, :)
   real(wp), intent(out) :: dadL(:, :, :)
   real(wp), intent(out) :: atrace(:, :)

   integer :: iat, jat, izp, jzp, ii, jj, ish, jsh
   real(wp) :: vec(3), r1, gam, arg, dtmp, dG(3), dS(3, 3)

   atrace(:, :) = 0.0_wp
   dadr(:, :, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:atrace, dadr, dadL) shared(mol, qvec, hardness, nshell, offset, gexp) &
   !$omp private(iat, izp, ii, ish, jat, jzp, jj, jsh, gam, r1, vec, dG, dS, dtmp, arg)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = offset(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         jj = offset(jat)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         r1 = norm2(vec)
         do ish = 1, nshell(iat)
            do jsh = 1, nshell(jat)
               gam = hardness(jsh, ish, jzp, izp)
               dtmp = 1.0_wp / (r1**gexp + gam**(-gexp))
               dtmp = -r1**(gExp-2.0_wp) * dtmp * dtmp**(1.0_wp/gExp)
               dG = dtmp*vec
               dS = spread(dG, 1, 3) * spread(vec, 2, 3)
               atrace(:, ii+ish) = +dG*qvec(jj+jsh) + atrace(:, ii+ish)
               atrace(:, jj+jsh) = -dG*qvec(ii+ish) + atrace(:, jj+jsh)
               dadr(:, iat, jj+jsh) = +dG*qvec(ii+ish) + dadr(:, iat, jj+jsh)
               dadr(:, jat, ii+ish) = -dG*qvec(jj+jsh) + dadr(:, jat, ii+ish)
               dadL(:, :, jj+jsh) = +dS*qvec(ii+ish) + dadL(:, :, jj+jsh)
               dadL(:, :, ii+ish) = +dS*qvec(jj+jsh) + dadL(:, :, ii+ish)
            end do
         end do
      end do
   end do

end subroutine get_damat_0d

subroutine get_damat_3d(mol, nshell, offset, hardness, gexp, wsc, alpha, qvec, &
      & dadr, dadL, atrace)
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: nshell(:)
   integer, intent(in) :: offset(:)
   real(wp), intent(in) :: hardness(:, :, :, :)
   real(wp), intent(in) :: gexp
   type(wignerseitz_cell), intent(in) :: wsc
   real(wp), intent(in) :: alpha
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(out) :: dadr(:, :, :)
   real(wp), intent(out) :: dadL(:, :, :)
   real(wp), intent(out) :: atrace(:, :)

   integer :: iat, jat, izp, jzp, img, ii, jj, ish, jsh
   real(wp) :: vol, gam, wsw, vec(3), dG(3), dS(3, 3)
   real(wp) :: dGd(3), dSd(3, 3), dGr(3), dSr(3, 3)
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   atrace(:, :) = 0.0_wp
   dadr(:, :, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, dtrans)
   call get_rec_trans(mol%lattice, rtrans)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:atrace, dadr, dadL) &
   !$omp shared(mol, wsc, alpha, vol, dtrans, rtrans, qvec, hardness, nshell, offset, gexp) &
   !$omp private(iat, izp, jat, jzp, img, ii, jj, ish, jsh, gam, wsw, vec, dG, dS, &
   !$omp& dGr, dSr, dGd, dSd)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = offset(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         jj = offset(jat)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
            wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
            call get_damat_rec_3d(vec, vol, alpha, rtrans, dGr, dSr)
            do ish = 1, nshell(iat)
               do jsh = 1, nshell(jat)
                  gam = hardness(jsh, ish, jzp, izp)
                  call get_damat_dir_3d(vec, gam, gexp, alpha, dtrans, dGd, dSd)
                  dG = (dGd + dGr) * wsw
                  dS = (dSd + dSr) * wsw
                  atrace(:, ii+ish) = +dG*qvec(jj+jsh) + atrace(:, ii+ish)
                  atrace(:, jj+jsh) = -dG*qvec(ii+ish) + atrace(:, jj+jsh)
                  dadr(:, iat, jj+jsh) = +dG*qvec(ii+ish) + dadr(:, iat, jj+jsh)
                  dadr(:, jat, ii+ish) = -dG*qvec(jj+jsh) + dadr(:, jat, ii+ish)
                  dadL(:, :, jj+jsh) = +dS*qvec(ii+ish) + dadL(:, :, jj+jsh)
                  dadL(:, :, ii+ish) = +dS*qvec(jj+jsh) + dadL(:, :, ii+ish)
               end do
            end do
         end do
      end do

      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
         call get_damat_rec_3d(vec, vol, alpha, rtrans, dGr, dSr)
         do ish = 1, nshell(iat)
            do jsh = 1, ish-1
               gam = hardness(jsh, ish, izp, izp)
               call get_damat_dir_3d(vec, gam, gexp, alpha, dtrans, dGd, dSd)
               dS = (dSd + dSr) * wsw
               dadL(:, :, ii+jsh) = +dS*qvec(ii+ish) + dadL(:, :, ii+jsh)
               dadL(:, :, ii+ish) = +dS*qvec(ii+jsh) + dadL(:, :, ii+ish)
            end do
            gam = hardness(ish, ish, izp, izp)
            call get_damat_dir_3d(vec, gam, gexp, alpha, dtrans, dGd, dSd)
            dS = (dSd + dSr) * wsw
            dadL(:, :, ii+ish) = +dS*qvec(ii+ish) + dadL(:, :, ii+ish)
         end do
      end do
   end do

end subroutine get_damat_3d

subroutine get_damat_dir_3d(rij, gam, gexp, alp, trans, dg, ds)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: gam
   real(wp), intent(in) :: gexp
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: dg(3)
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: vec(3), r1, r2, gtmp, atmp, alp2

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   alp2 = alp*alp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      r2 = r1*r1
      gtmp = 1.0_wp / (r1**gexp + gam**(-gexp))
      gtmp = -r1**(gexp-2.0_wp) * gtmp * gtmp**(1.0_wp/gexp)
      atmp = -2*alp*exp(-r2*alp2)/(sqrtpi*r2) + erf(r1*alp)/(r2*r1)
      dg(:) = dg + (gtmp + atmp) * vec
      ds(:, :) = ds + (gtmp + atmp) * spread(vec, 1, 3) * spread(vec, 2, 3)
   end do

end subroutine get_damat_dir_3d

subroutine get_damat_rec_3d(rij, vol, alp, trans, dg, ds)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: vol
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: dg(3)
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: fac, vec(3), g2, gv, expk, sink, cosk, alp2
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp
   fac = 4*pi/vol
   alp2 = alp*alp

   do itr = 1, size(trans, 2)
      vec(:) = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      gv = dot_product(rij, vec)
      expk = fac * exp(-0.25_wp*g2/alp2)/g2
      cosk = cos(gv) * expk
      sink = sin(gv) * expk
      dg(:) = dg - sink * vec
      ds(:, :) = ds + cosk &
         & * ((2.0_wp/g2 + 0.5_wp/alp2) * spread(vec, 1, 3)*spread(vec, 2, 3) - unity)
   end do

end subroutine get_damat_rec_3d


pure function variable_info(self) result(info)
   use tblite_scf_info, only : scf_info, shell_resolved
   !> Instance of the electrostatic container
   class(effective_coulomb), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=shell_resolved)
end function variable_info


end module tblite_coulomb_charge
