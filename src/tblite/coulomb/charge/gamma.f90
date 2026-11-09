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

!> @file tblite/coulomb/charge/gamma.f90
!> Provides the DFTB Coulomb functional for isotropic electrostatic interactions.

!> Isotropic second-order electrostatics using the DFTB Coulomb functional.
module tblite_coulomb_charge_gamma
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_charge_type, only : coulomb_charge_type
   use tblite_coulomb_ewald, only : ewald_cache, get_amat_rec_3d, get_damat_rec_3d
   use tblite_partition, only : pair_list
   use tblite_wignerseitz, only : wignerseitz_cell, get_wignerseitz_weights
   implicit none
   private

   public :: new_gamma_coulomb


   !> DFTB Gamma functional second-order electrostatics
   type, public, extends(coulomb_charge_type) :: gamma_coulomb
      !> Chemical hardness for each shell and species
      real(wp), allocatable :: hubbard(:, :)
   contains
      !> Evaluate the Coulomb matrix blocks of the locally owned atom pairs
      procedure :: get_coulomb_matrix
      !> Evaluate the Coulomb matrix derivatives of one atom pair
      procedure :: get_pair_derivs
   end type gamma_coulomb


   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))
   real(wp), parameter :: conv = epsilon(0.0_wp)
   character(len=*), parameter :: label = "isotropic DFTB-γ electrostatics"

contains

!> Construct new effective electrostatic interaction container
subroutine new_gamma_coulomb(self, mol, hubbard, nshell)
   !> Instance of the electrostatic container
   type(gamma_coulomb), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Chemical hardness for all shells and species
   real(wp), intent(in) :: hubbard(:, :)
   !> Number of shells for each species
   integer, intent(in), optional :: nshell(:)

   integer :: mshell
   integer :: isp, ish, ind, iat

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

   if (present(nshell)) then
      allocate(self%hubbard(mshell, mol%nid))
      do isp = 1, mol%nid
         self%hubbard(:, isp) = 0.0_wp
         do ish = 1, nshell(isp)
            self%hubbard(ish, isp) = hubbard(ish, isp)
         end do
      end do
   else
      allocate(self%hubbard(1, mol%nid))
      do isp = 1, mol%nid
         self%hubbard(1, isp) = hubbard(1, isp)
      end do
   end if

end subroutine new_gamma_coulomb


!> Evaluate the Coulomb matrix blocks of the locally owned atom pairs
subroutine get_coulomb_matrix(self, mol, cache)
   !> Instance of the electrostatic container
   class(gamma_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container, local pair list and Coulomb blocks to overwrite
   type(coulomb_cache), intent(inout) :: cache

   cache%local_amat(:, :, :) = 0.0_wp
   if (any(mol%periodic)) then
      call cache%charge_ewald%update(mol%lattice, cache%alpha, conv)
      call get_amat_3d(mol, self%nshell, self%hubbard, &
         & cache%wsc, cache%alpha, cache%charge_ewald, cache%local_amat, cache%pairs)
   else
      call get_amat_0d(mol, self%nshell, self%hubbard, cache%local_amat, cache%pairs)
   end if
end subroutine get_coulomb_matrix


!> Evaluate Coulomb matrix for finite systems
subroutine get_amat_0d(mol, nshell, hubbard, amat, pairs)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each atom
   integer, intent(in) :: nshell(:)
   !> Hardness parameter for each shell
   real(wp), intent(in) :: hubbard(:, :)
   !> Coulomb blocks in local pair order
   real(wp), intent(inout) :: amat(:, :, :)
   !> Locally owned atom pairs in row order
   type(pair_list), intent(in) :: pairs

   integer :: iat, jat, izp, jzp, ipair, jpair, ish, jsh
   real(wp) :: vec(3), r1, gam, tmp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(amat, mol, nshell, hubbard, pairs) &
   !$omp private(iat, izp, ipair, ish, jat, jzp, jpair, jsh, gam, vec, r1, tmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do ipair = pairs%offset(iat-1)+1, pairs%offset(iat)
         jat = pairs%neighbour(ipair)
         if (jat > iat) cycle
         jpair = pairs%find(jat, iat)
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         tmp = 0.0_wp
         if (iat /= jat) tmp = 1.0_wp / r1
         do ish = 1, nshell(iat)
            do jsh = 1, merge(ish, nshell(jat), iat == jat)
               gam = tmp - exp_gamma(r1, hubbard(ish, izp), hubbard(jsh, jzp))
               if (iat /= jat .or. ish /= jsh) amat(jsh, ish, jpair) = amat(jsh, ish, jpair) + gam
               amat(ish, jsh, ipair) = amat(ish, jsh, ipair) + gam
            end do
         end do
      end do
   end do

end subroutine get_amat_0d

!> Evaluate the coulomb matrix for 3D systems
subroutine get_amat_3d(mol, nshell, hubbard, wsc, alpha, ewald, amat, pairs)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells per atom
   integer, intent(in) :: nshell(:)
   !> Hardness of the shells
   real(wp), intent(in) :: hubbard(:, :)
   !> Wigner-Seitz cell
   type(wignerseitz_cell), intent(in) :: wsc
   !> Convergence factor
   real(wp), intent(in) :: alpha
   !> Cached Ewald lattice data
   type(ewald_cache), intent(in) :: ewald
   !> Coulomb blocks in local pair order
   real(wp), intent(inout) :: amat(:, :, :)
   !> Locally owned atom pairs in row order
   type(pair_list), intent(in) :: pairs

   integer :: iat, jat, izp, jzp, img, ipair, jpair, ish, jsh
   real(wp) :: vec(3), ui, uj, dtmp, rtmp, aval
   real(wp) :: weight(size(wsc%tridx, 1))

   !$omp parallel do default(none) schedule(runtime) shared(amat) &
   !$omp shared(mol, nshell, hubbard, wsc, ewald, alpha, pairs) &
   !$omp private(iat, izp, jat, jzp, ipair, jpair, ish, jsh, ui, uj, weight, vec, dtmp, rtmp, aval)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do ipair = pairs%offset(iat-1)+1, pairs%offset(iat)
         jat = pairs%neighbour(ipair)
         if (jat > iat) cycle
         jpair = pairs%find(jat, iat)
         jzp = mol%id(jat)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         call get_wignerseitz_weights(wsc, jat, iat, vec, weight)
         call get_amat_rec_3d(vec, ewald, rtmp)
         if (iat == jat) rtmp = rtmp - 2 * alpha / sqrtpi
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
            do ish = 1, nshell(iat)
               do jsh = 1, merge(ish, nshell(jat), iat == jat)
                  ui = hubbard(ish, izp)
                  uj = hubbard(jsh, jzp)
                  call get_amat_dir_3d(vec, ui, uj, alpha, dtmp)
                  aval = dtmp + rtmp
                  if (iat == jat) aval = aval - exp_gamma(0.0_wp, ui, uj)
                  aval = aval * weight(img)
                  if (iat /= jat .or. ish /= jsh) amat(jsh, ish, jpair) = amat(jsh, ish, jpair) + aval
                  amat(ish, jsh, ipair) = amat(ish, jsh, ipair) + aval
               end do
            end do
         end do
      end do

   end do

end subroutine get_amat_3d

!> Calculate real space contributions for a pair under 3D periodic boundary conditions
subroutine get_amat_dir_3d(rij, ui, uj, alp, amat)
   !> Distance between pair
   real(wp), intent(in) :: rij(3)
   !> Chemical hardness
   real(wp), intent(in) :: ui, uj
   !> Convergence factor
   real(wp), intent(in) :: alp
   !> Interaction matrix element
   real(wp), intent(out) :: amat

   real(wp) :: r1

   amat = 0.0_wp

   r1 = norm2(rij)
   if (r1 < eps) return
   amat = 1.0_wp / r1 - exp_gamma(r1, ui, uj) &
      & - erf(alp*r1)/r1

end subroutine get_amat_dir_3d

!> Evaluate the Coulomb matrix derivatives of one atom pair. For iat == jat only the
!> strain derivatives of the shell pairs jsh <= ish are used, otherwise all shell pairs.
subroutine get_pair_derivs(self, mol, cache, iat, jat, dG, dS)
   !> Instance of the electrostatic container
   class(gamma_coulomb), intent(in) :: self
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
      call get_damat_pair_3d(mol, self%nshell, self%hubbard, cache%wsc, cache%alpha, &
         & cache%charge_ewald, iat, jat, dG, dS)
   else
      call get_damat_pair_0d(mol, self%nshell, self%hubbard, iat, jat, dG)
   end if
end subroutine get_pair_derivs


!> Coulomb matrix derivatives of one atom pair in a finite system
subroutine get_damat_pair_0d(mol, nshell, hubbard, iat, jat, dG)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each atom
   integer, intent(in) :: nshell(:)
   !> Chemical hardness for each shell and species
   real(wp), intent(in) :: hubbard(:, :)
   !> First atom
   integer, intent(in) :: iat
   !> Second atom
   integer, intent(in) :: jat
   !> Derivatives with respect to the position of the first atom
   real(wp), intent(out) :: dG(:, :, :)

   integer :: izp, jzp, ish, jsh
   real(wp) :: vec(3), r1, gam, dtmp

   izp = mol%id(iat)
   jzp = mol%id(jat)
   vec = mol%xyz(:, iat) - mol%xyz(:, jat)
   r1 = norm2(vec)
   dtmp = -1.0_wp / (r1*r1*r1)
   do ish = 1, nshell(iat)
      do jsh = 1, nshell(jat)
         gam = dtmp - dexp_gamma(r1, hubbard(ish, izp), hubbard(jsh, jzp))/r1
         dG(:, jsh, ish) = gam*vec
      end do
   end do
end subroutine get_damat_pair_0d

!> Coulomb matrix derivatives of one atom pair in a 3D system using Wigner-Seitz images
subroutine get_damat_pair_3d(mol, nshell, hubbard, wsc, alpha, ewald, iat, jat, dG, dS)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each atom
   integer, intent(in) :: nshell(:)
   !> Chemical hardness for each shell and species
   real(wp), intent(in) :: hubbard(:, :)
   !> Wigner-Seitz image information
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

   integer :: izp, jzp, img, ish, jsh
   logical :: need_weight_energy
   real(wp) :: ui, uj, dtmp, vec(3), dGd(3), dSd(3, 3), dGr(3), dSr(3, 3)
   real(wp) :: weight(size(wsc%tridx, 1))
   real(wp) :: dwdr(3, size(wsc%tridx, 1)), dwdL(3, 3, size(wsc%tridx, 1))

   izp = mol%id(iat)
   jzp = mol%id(jat)
   vec = mol%xyz(:, iat) - mol%xyz(:, jat)
   call get_wignerseitz_weights(wsc, jat, iat, vec, weight, dwdr, dwdL)
   need_weight_energy = any(dwdr(:, :wsc%nimg(jat, iat)) /= 0.0_wp) &
      & .or. any(dwdL(:, :, :wsc%nimg(jat, iat)) /= 0.0_wp)
   call get_damat_rec_3d(vec, ewald, dGr, dSr)
   dG = 0.0_wp
   dS = 0.0_wp
   do img = 1, wsc%nimg(jat, iat)
      vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
      do ish = 1, nshell(iat)
         ui = hubbard(ish, izp)
         do jsh = 1, merge(ish, nshell(jat), iat == jat)
            uj = hubbard(jsh, jzp)
            dtmp = 0.0_wp
            if (need_weight_energy) call get_amat_dir_3d(vec, ui, uj, alpha, dtmp)
            call get_damat_dir_3d(vec, ui, uj, alpha, dGd, dSd)
            dG(:, jsh, ish) = dG(:, jsh, ish) + (dGd + dGr)*weight(img) + dtmp*dwdr(:, img)
            dS(:, :, jsh, ish) = dS(:, :, jsh, ish) + (dSd + dSr)*weight(img) &
               & + dtmp*dwdL(:, :, img)
         end do
      end do
   end do
end subroutine get_damat_pair_3d

!> Calculate real space contributions for a pair under 3D periodic boundary conditions
subroutine get_damat_dir_3d(rij, ui, uj, alp, dg, ds)
   !> Distance between pair
   real(wp), intent(in) :: rij(3)
   !> Chemical hardness
   real(wp), intent(in) :: ui, uj
   !> Convergence factor
   real(wp), intent(in) :: alp
   !> Derivative with respect to cartesian displacements
   real(wp), intent(out) :: dg(3)
   !> Derivative with respect to strain deformations
   real(wp), intent(out) :: ds(3, 3)

   real(wp) :: r1, r2, gtmp, atmp, alp2

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   alp2 = alp*alp

   r1 = norm2(rij)
   if (r1 < eps) return
   r2 = r1*r1
   gtmp = -1.0_wp/(r2*r1) - dexp_gamma(r1, ui, uj)/r1
   atmp = -2*alp*exp(-r2*alp2)/(sqrtpi*r2) + erf(r1*alp)/(r2*r1)
   dg(:) = (gtmp + atmp) * rij
   ds(:, :) = (gtmp + atmp) * spread(rij, 1, 3) * spread(rij, 2, 3)

end subroutine get_damat_dir_3d


!> Determines the value of the short range contribution to gamma with the exponential form
pure function exp_gamma(r1, ui, uj)
   !> separation of sites a and b
   real(wp), intent(in) :: r1
   !> Hubbard U for site a
   real(wp), intent(in) :: ui
   !> Hubbard U for site b
   real(wp), intent(in) :: uj
   !> contribution
   real(wp) :: exp_gamma

   real(wp) :: taui, tauj, tauij

   ! 16/5 * U, see review papers / theses
   taui = 3.2_wp*ui
   tauj = 3.2_wp*uj

   if (r1 < eps) then
      ! on-site case with R~0
      if (abs(ui - uj) < eps) then
         ! same Hubbard U values, onsite
         exp_gamma = -0.5_wp*(ui + uj)
      else
         ! ui /= uj Hubbard U values - limiting case
         exp_gamma = &
            & -0.5_wp*((taui*tauj)/(taui+tauj) + (taui*tauj)**2/(taui+tauj)**3)
      end if
   else
      if (abs(ui - uj) < eps) then
         ! R > 0 and same Hubbard U values
         tauij = 0.5_wp*(taui + tauj)
         exp_gamma = &
            & exp(-tauij*r1) &
            & * (48.0_wp/r1 &
            &   + 33.0_wp*tauij &
            &   + 9.0_wp*r1*(tauij**2) &
            &   + (r1**2)*(tauij**3)) / 48.0_wp
      else
         exp_gamma = gamma_sf(r1,taui,tauj) + gamma_sf(r1,tauj,taui)
      end if
   end if

end function exp_gamma

!> Determines the value of a part of the short range contribution to the exponential gamma, when
!> ui /= uj and R > 0
pure function gamma_sf(r1, taui, tauj) result(sf)
   !> separation of sites a and b
   real(wp), intent(in) :: r1
   !> Charge fluctuation for site a
   real(wp), intent(in) :: taui
   !> Charge fluctuation U for site b
   real(wp), intent(in) :: tauj
   !> contribution
   real(wp) :: sf

   sf = exp(-taui * r1) &
      & * ((0.5_wp*tauj**4*taui/(taui**2-tauj**2)**2) &
      &   - (tauj**6-3.0_wp*tauj**4*taui**2)/(r1*(taui**2-tauj**2)**3))

end function gamma_sf

!> Determines the value of the derivative of the short range contribution to gamma with the
!> exponential form
pure function dexp_gamma(r1,ui,uj)
   !> separation of sites a and b
   real(wp), intent(in) :: r1
   !> Hubbard U for site a
   real(wp), intent(in) :: ui
   !> Hubbard U for site b
   real(wp), intent(in) :: uj
   !> returned contribution
   real(wp) :: dexp_gamma

   real(wp) :: taui, tauj, tauij

   dexp_gamma = 0.0_wp
   if (r1 < eps) return

   ! 16/5 * U, see review papers
   taui = 3.2_wp*ui
   tauj = 3.2_wp*uj

   if (abs(ui - uj) < eps) then
      ! R > 0 and same Hubbard U values
      tauij = 0.5_wp * (taui + tauj)
      dexp_gamma = &
         & -exp(-tauij*r1) &
         & * (tauij**4*r1**4 &
         &   + 7*tauij**3*r1**3 &
         &   + 24*tauij**2*r1**2 &
         &   + 48*tauij*r1 &
         &   + 48.0_wp) &
         & / (48.0_wp*r1**2)
   else
      dexp_gamma = dgamma_sf(r1,taui,tauj) + dgamma_sf(r1,tauj,taui)
   end if
end function dexp_gamma

!> Determines the derivative of the value of a part of the short range contribution to the
!> exponential gamma, when ui /= uj and R > 0
pure function dgamma_sf(r1, taui, tauj) result(dsf)
   !> separation of sites a and b
   real(wp), intent(in) :: r1
   !> Charge fluctuation for site a
   real(wp), intent(in) :: taui
   !> Charge fluctuation U for site b
   real(wp), intent(in) :: tauj
   !> contribution
   real(wp) :: dsf

   dsf = &
      & exp(- taui * r1) &
      & * (taui * tauj**4 * (tauj**2 - 3.0_wp*taui**2)/(r1*(taui**2-tauj**2)**3) &
      &   - tauj**4 * taui**2 / (2*(taui**2-tauj**2)**2) &
      &   + tauj**4*(tauj**2 - 3.0_wp*taui**2) / (r1**2 *(taui**2-tauj**2)**3))

end function dgamma_sf


end module tblite_coulomb_charge_gamma
