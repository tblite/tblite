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

!> Implementation of the effective core Hamiltonian used in the extended tight binding.
module tblite_xtb_h0
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_integral_multipole, only : multipole_cgto, multipole_grad_cgto, maxl, msao
   use tblite_scf_potential, only : potential_type
   use tblite_xtb_spec, only : tb_h0spec
   implicit none
   private

   public :: tb_hamiltonian, new_hamiltonian
   public :: get_selfenergy, get_hamiltonian, get_occupation, get_hamiltonian_gradient


   type :: tb_hamiltonian
      !> Atomic level information
      real(wp), allocatable :: selfenergy(:, :)
      !> Coordination number dependence of the atomic levels
      real(wp), allocatable :: kcn(:, :)
      !> Charge dependence of the atomic levels
      real(wp), allocatable :: kq1(:, :)
      !> Charge dependence of the atomic levels
      real(wp), allocatable :: kq2(:, :)
      !> Enhancement factor to scale the Hamiltonian elements
      real(wp), allocatable :: hscale(:, :, :, :)
      !> Polynomial coefficients for distance dependent enhancement factor
      real(wp), allocatable :: shpoly(:, :)
      !> Atomic radius for polynomial enhancement
      real(wp), allocatable :: rad(:)
      !> Reference occupation numbers
      real(wp), allocatable :: refocc(:, :)
   end type tb_hamiltonian


contains


!> Constructor for a new Hamiltonian object, consumes a Hamiltonian specification
subroutine new_hamiltonian(self, mol, bas, spec)
   type(tb_hamiltonian), intent(out) :: self
   type(structure_type), intent(in) :: mol
   type(basis_type), intent(in) :: bas
   class(tb_h0spec), intent(in) :: spec

   integer :: mshell

   mshell = maxval(bas%nsh_id)
   allocate(self%selfenergy(mshell, mol%nid), self%kcn(mshell, mol%nid), &
      & self%kq1(mshell, mol%nid), self%kq2(mshell, mol%nid))
   call spec%get_selfenergy(mol, bas, self%selfenergy)
   call spec%get_cnshift(mol, bas, self%kcn)
   call spec%get_q1shift(mol, bas, self%kq1)
   call spec%get_q2shift(mol, bas, self%kq2)

   allocate(self%hscale(mshell, mshell, mol%nid, mol%nid))
   call spec%get_hscale(mol, bas, self%hscale)

   allocate(self%shpoly(mshell, mol%nid), self%rad(mol%nid))
   call spec%get_rad(mol, bas, self%rad)
   call spec%get_shpoly(mol, bas, self%shpoly)

   allocate(self%refocc(mshell, mol%nid))
   call spec%get_reference_occ(mol, bas, self%refocc)
end subroutine new_hamiltonian


subroutine get_selfenergy(h0, id, ish_at, nshell, cn, qat, selfenergy, dsedcn, dsedq)
   type(tb_hamiltonian), intent(in) :: h0
   integer, intent(in) :: id(:)
   integer, intent(in) :: ish_at(:)
   integer, intent(in) :: nshell(:)
   real(wp), intent(in), optional :: cn(:)
   real(wp), intent(in), optional :: qat(:)
   real(wp), intent(out) :: selfenergy(:)
   real(wp), intent(out), optional :: dsedcn(:)
   real(wp), intent(out), optional :: dsedq(:)

   integer :: iat, izp, ish, ii

   selfenergy(:) = 0.0_wp
   if (present(dsedcn)) dsedcn(:) = 0.0_wp
   if (present(dsedq)) dsedq(:) = 0.0_wp
   do iat = 1, size(id)
      izp = id(iat)
      ii = ish_at(iat)
      do ish = 1, nshell(izp)
         selfenergy(ii+ish) = h0%selfenergy(ish, izp)
      end do
   end do
   if (present(cn)) then
      if (present(dsedcn)) then
         do iat = 1, size(id)
            izp = id(iat)
            ii = ish_at(iat)
            do ish = 1, nshell(izp)
               selfenergy(ii+ish) = selfenergy(ii+ish) - h0%kcn(ish, izp) * cn(iat)
               dsedcn(ii+ish) = -h0%kcn(ish, izp)
            end do
         end do
      else
         do iat = 1, size(id)
            izp = id(iat)
            ii = ish_at(iat)
            do ish = 1, nshell(izp)
               selfenergy(ii+ish) = selfenergy(ii+ish) - h0%kcn(ish, izp) * cn(iat)
            end do
         end do
      end if
   end if
   if (present(qat)) then
      if (present(dsedq)) then
         do iat = 1, size(id)
            izp = id(iat)
            ii = ish_at(iat)
            do ish = 1, nshell(izp)
               selfenergy(ii+ish) = selfenergy(ii+ish) &
                  & - h0%kq1(ish, izp)*qat(iat) - h0%kq2(ish, izp)*qat(iat)**2
               dsedq(ii+ish) = -h0%kq1(ish, izp) - h0%kq2(ish, izp)*2*qat(iat)
            end do
         end do
      else
         do iat = 1, size(id)
            izp = id(iat)
            ii = ish_at(iat)
            do ish = 1, nshell(izp)
               selfenergy(ii+ish) = selfenergy(ii+ish) &
                  & - h0%kq1(ish, izp)*qat(iat) - h0%kq2(ish, izp)*qat(iat)**2
            end do
         end do
      end if
   end if

end subroutine get_selfenergy


subroutine get_h0scale(mol, bas, en, enscale2, enscale4, zeta, wexp, valence, &
      & kpair, kshell, kdiff, hscale)
   type(structure_type), intent(in) :: mol
   type(basis_type), intent(in) :: bas
   real(wp), intent(in) :: en(:)
   real(wp), intent(in) :: enscale2(0:, 0:)
   real(wp), intent(in) :: enscale4(0:, 0:)
   real(wp), intent(in) :: zeta(:, :)
   real(wp), intent(in) :: wexp
   logical, intent(in) :: valence(:, :)
   real(wp), intent(in) :: kpair(:, :)
   real(wp), intent(in) :: kshell(0:, 0:)
   real(wp), intent(in) :: kdiff
   real(wp), intent(out) :: hscale(:, :, :, :)

   integer :: izp, jzp, ish, jsh, il, jl
   real(wp) :: zi, zj, zij, den, enp, km

   hscale(:, :, :, :) = 0.0_wp

   do izp = 1, mol%nid
      do jzp = 1, mol%nid
         den = (en(izp) - en(jzp))**2
         do ish = 1, bas%nsh_id(izp)
            il = bas%cgto(ish, izp)%ang
            do jsh = 1, bas%nsh_id(jzp)
               jl = bas%cgto(jsh, jzp)%ang
               zi = zeta(ish, izp)
               zj = zeta(jsh, jzp)
               zij = (2*sqrt(zi*zj)/(zi+zj))**wexp
               if (valence(ish, izp) .and. valence(jsh, jzp)) then
                  enp = 1.0_wp + enscale2(jl, il) * den * (1.0_wp + enscale4(jl, il) * den)
                  km = kpair(jzp, izp) * kshell(jl, il) * enp
               else if (valence(ish, izp)) then
                  km = 0.5_wp * (kshell(il, il) + kdiff)
               else if (valence(jsh, jzp)) then
                  km = 0.5_wp * (kshell(jl, jl) + kdiff)
               else
                  km = kdiff
               end if
               hscale(jsh, ish, jzp, izp) = zij * km
            end do
         end do
      end do
   end do

end subroutine get_h0scale


subroutine get_hamiltonian(mol, trans, cutoff, bas, h0, selfenergy, overlap, dpint, qpint, &
      & hamiltonian)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Hamiltonian interaction data
   type(tb_hamiltonian), intent(in) :: h0
   !> Diagonal elememts of the Hamiltonian
   real(wp), intent(in) :: selfenergy(:)
   !> Overlap integral matrix
   real(wp), intent(out) :: overlap(:, :)
   !> Dipole moment integral matrix
   real(wp), intent(out) :: dpint(:, :, :)
   !> Quadrupole moment integral matrix
   real(wp), intent(out) :: qpint(:, :, :)
   !> Effective Hamiltonian
   real(wp), intent(out) :: hamiltonian(:, :)

   integer :: iat, jat, izp, jzp, itr
   integer :: ish, jsh, is, js, ii, jj, iao, jao, nao
   real(wp) :: rr, r2, vec(3), cutoff2, hij, shpoly
   real(wp), allocatable :: stmp(:), dtmp(:, :), qtmp(:, :)

   overlap(:, :) = 0.0_wp
   dpint(:, :, :) = 0.0_wp
   qpint(:, :, :) = 0.0_wp
   hamiltonian(:, :) = 0.0_wp

   allocate(stmp(msao(bas%maxl)**2), dtmp(3, msao(bas%maxl)**2), qtmp(6, msao(bas%maxl)**2))
   cutoff2 = cutoff**2

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, bas, trans, cutoff2, overlap, dpint, qpint, hamiltonian, h0, selfenergy) &
   !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj, iao, jao, nao) &
   !$omp private(r2, vec, stmp, dtmp, qtmp, hij, shpoly, rr)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      do jat = 1, mol%nat
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         do itr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            if (r2 > cutoff2) cycle
            rr = sqrt(sqrt(r2) / (h0%rad(jzp) + h0%rad(izp)))
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is+ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js+jsh)
                  call multipole_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp, dtmp, qtmp)

                  shpoly = (1.0_wp + h0%shpoly(ish, izp)*rr) &
                     & * (1.0_wp + h0%shpoly(jsh, jzp)*rr)

                  hij = 0.5_wp * (selfenergy(is+ish) + selfenergy(js+jsh)) &
                     * merge(1.0_wp, h0%hscale(jsh, ish, jzp, izp), r2 < epsilon(cutoff2)) &
                     * shpoly


                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                           & + stmp(jao + nao*(iao-1))

                        dpint(:, jj+jao, ii+iao) = dpint(:, jj+jao, ii+iao) &
                           & + dtmp(:, jao + nao*(iao-1))

                        qpint(:, jj+jao, ii+iao) = qpint(:, jj+jao, ii+iao) &
                           & + qtmp(:, jao + nao*(iao-1))

                        hamiltonian(jj+jao, ii+iao) = hamiltonian(jj+jao, ii+iao) &
                           & + stmp(jao + nao*(iao-1)) * hij
                     end do
                  end do

               end do
            end do

         end do
      end do
   end do

end subroutine get_hamiltonian


subroutine get_hamiltonian_gradient(mol, trans, cutoff, bas, h0, selfenergy, dsedcn, &
      & pot, pmat, xmat, dEdcn, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Hamiltonian interaction data
   type(tb_hamiltonian), intent(in) :: h0
   !> Diagonal elememts of the Hamiltonian
   real(wp), intent(in) :: selfenergy(:)
   !> Derivative of the diagonal elements of the Hamiltonian w.r.t. the coordination number
   real(wp), intent(in) :: dsedcn(:)
   !> Density dependent potential shifts on the Hamiltonian
   type(potential_type), intent(in) :: pot
   !> Density matrix
   real(wp), intent(in) :: pmat(:, :)
   !> Energy weighted density matrix
   real(wp), intent(in) :: xmat(:, :)

   !> Derivative of the electronic energy w.r.t. the coordination number
   real(wp), intent(inout) :: dEdcn(:)
   !> Derivative of the electronic energy w.r.t. coordinate displacements
   real(wp), intent(inout) :: gradient(:, :)
   !> Derivative of the electronic energy w.r.t. strain deformations
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, izp, jzp, itr
   integer :: ish, jsh, is, js, ii, jj, iao, jao, nao, ij
   real(wp) :: rr, r2, vec(3), cutoff2, hij, shpoly, dshpoly, dG(3), hscale
   real(wp) :: sval, dcni, dcnj, dhdcni, dhdcnj, hpij, pij
   real(wp), allocatable :: stmp(:), dtmp(:, :), qtmp(:, :)
   real(wp), allocatable :: dstmp(:, :), ddtmpi(:, :, :), dqtmpi(:, :, :)
   real(wp), allocatable :: ddtmpj(:, :, :), dqtmpj(:, :, :)

   allocate(stmp(msao(bas%maxl)**2), dstmp(3, msao(bas%maxl)**2), &
      & dtmp(3, msao(bas%maxl)**2), ddtmpi(3, 3, msao(bas%maxl)**2), &
      & qtmp(6, msao(bas%maxl)**2), dqtmpi(3, 6, msao(bas%maxl)**2), &
      & ddtmpj(3, 3, msao(bas%maxl)**2), dqtmpj(3, 6, msao(bas%maxl)**2))
   cutoff2 = cutoff**2

   !$omp parallel do schedule(runtime) default(none) reduction(+:dEdcn, gradient, sigma) &
   !$omp shared(mol, bas, trans, cutoff2, h0, selfenergy, dsedcn, pot, pmat, xmat) &
   !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj, iao, jao, nao, ij, &
   !$omp& r2, vec, stmp, dtmp, qtmp, dstmp, ddtmpi, dqtmpi, ddtmpj, dqtmpj, hij, shpoly, &
   !$omp& dshpoly, dG, dcni, dcnj, dhdcni, dhdcnj, hpij, rr, sval, hscale, pij)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         do itr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            if (r2 > cutoff2 .or. r2 < epsilon(cutoff2)) cycle
            rr = sqrt(sqrt(r2) / (h0%rad(jzp) + h0%rad(izp)))
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is+ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js+jsh)
                  call multipole_grad_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp, dtmp, qtmp, dstmp, ddtmpj, dqtmpj, &
                     & ddtmpi, dqtmpi)

                  shpoly = (1.0_wp + h0%shpoly(ish, izp)*rr) &
                     & * (1.0_wp + h0%shpoly(jsh, jzp)*rr)
                  dshpoly = ((1.0_wp + h0%shpoly(ish, izp)*rr)*h0%shpoly(jsh, jzp)*rr &
                     & + (1.0_wp + h0%shpoly(jsh, jzp)*rr)*h0%shpoly(ish, izp)*rr) &
                     & * 0.5_wp / r2

                  hscale = h0%hscale(jsh, ish, jzp, izp)
                  hij = 0.5_wp * (selfenergy(is+ish) + selfenergy(js+jsh)) * hscale
                  dhdcni = dsedcn(is+ish) * shpoly * hscale
                  dhdcnj = dsedcn(js+jsh) * shpoly * hscale

                  dG(:) = 0.0_wp
                  dcni = 0.0_wp
                  dcnj = 0.0_wp
                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)
                        pij = pmat(jj+jao, ii+iao)
                        hpij = pij * hij * shpoly
                        sval = 2*hpij - 2*xmat(jj+jao, ii+iao) &
                           & - pij * (pot%vao(jj+jao) + pot%vao(ii+iao))

                        dG(:) = dG + sval * dstmp(:, ij) &
                           & + 2*hpij*stmp(ij) * dshpoly / shpoly * vec &
                           & - pij * matmul(ddtmpi(:, :, ij), pot%vdp(:, iat)) &
                           & - pij * matmul(ddtmpj(:, :, ij), pot%vdp(:, jat)) &
                           & - pij * matmul(dqtmpi(:, :, ij), pot%vqp(:, iat)) &
                           & - pij * matmul(dqtmpj(:, :, ij), pot%vqp(:, jat))

                        dcni = dcni + dhdcni * pmat(jj+jao, ii+iao) * stmp(ij)
                        dcnj = dcnj + dhdcnj * pmat(jj+jao, ii+iao) * stmp(ij)
                     end do
                  end do
                  dEdcn(iat) = dEdcn(iat) + dcni
                  dEdcn(jat) = dEdcn(jat) + dcnj
                  gradient(:, iat) = gradient(:, iat) + dG
                  gradient(:, jat) = gradient(:, jat) - dG
                  sigma(:, :) = sigma + 0.5_wp * (spread(vec, 1, 3) * spread(dG, 2, 3) &
                     & + spread(dG, 1, 3) * spread(vec, 2, 3))

               end do
            end do

         end do
      end do
   end do

   !$omp parallel do schedule(runtime) default(none) reduction(+:dEdcn) &
   !$omp shared(mol, bas, dsedcn, pmat) &
   !$omp private(iat, izp, jzp, is, ish, ii, iao, dcni, dhdcni)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      do ish = 1, bas%nsh_id(izp)
         ii = bas%iao_sh(is+ish)
         dhdcni = dsedcn(is+ish)
         dcni = 0.0_wp
         do iao = 1, msao(bas%cgto(ish, izp)%ang)
            dcni = dcni + dhdcni * pmat(ii+iao, ii+iao)
         end do
         dEdcn(iat) = dEdcn(iat) + dcni
      end do
   end do

end subroutine get_hamiltonian_gradient


subroutine get_occupation(mol, bas, h0, nocc, n0at, n0sh)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Hamiltonian interaction data
   type(tb_hamiltonian), intent(in) :: h0
   !> Occupation number
   real(wp), intent(out) :: nocc
   !> Reference occupation for each atom
   real(wp), intent(out) :: n0at(:)
   !> Reference occupation for each shell
   real(wp), intent(out) :: n0sh(:)

   integer :: iat, ish, izp, ii

   nocc = -mol%charge
   n0at(:) = 0.0_wp
   n0sh(:) = 0.0_wp
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = bas%ish_at(iat)
      do ish = 1, bas%nsh_id(izp)
         nocc = nocc + h0%refocc(ish, izp)
         n0at(iat) = n0at(iat) + h0%refocc(ish, izp)
         n0sh(ii+ish) = n0sh(ii+ish) + h0%refocc(ish, izp)
      end do
   end do

end subroutine get_occupation


end module tblite_xtb_h0
