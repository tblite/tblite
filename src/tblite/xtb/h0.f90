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

!> @file tblite/xtb/h0.f90
!> Provides the effective core Hamiltonian for xTB.

!> Implementation of the effective core Hamiltonian used in the extended tight binding.
module tblite_xtb_h0
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_adjlist, only : adjacency_list
   use tblite_basis_type, only : basis_type
   use tblite_integral_multipole, only : multipole_cgto, multipole_grad_cgto, maxl, msao
   use tblite_scf_potential, only : potential_type
   use tblite_xtb_spec, only : tb_h0spec
   implicit none
   private

   public ::  new_hamiltonian
   public :: get_selfenergy, get_hamiltonian, get_occupation, get_hamiltonian_gradient


   type, public :: tb_hamiltonian
      !> Atomic level information
      real(wp), allocatable :: selfenergy(:, :)
      !> Coordination number dependence of the atomic levels
      real(wp), allocatable :: kcn(:, :)
      !> Electronegativity scaled coordination number dependence of the atomic levels
      real(wp), allocatable :: kcn_en(:, :)
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
      !> Diatomic frame scaling of sigma bonding contribution
      real(wp), allocatable :: ksig(:, :)
      !> Diatomic frame scaling of pi bonding contribution
      real(wp), allocatable :: kpi(:, :)
      !> Diatomic frame scaling of delta bonding contribution
      real(wp), allocatable :: kdel(:, :)
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
   allocate(self%selfenergy(mshell, mol%nid), &
      & self%kcn(mshell, mol%nid), self%kcn_en(mshell, mol%nid), & 
      & self%kq1(mshell, mol%nid), self%kq2(mshell, mol%nid))
   call spec%get_selfenergy(mol, bas, self%selfenergy)
   call spec%get_cnshift(mol, bas, self%kcn)
   call spec%get_cnenshift(mol, bas, self%kcn_en)
   call spec%get_q1shift(mol, bas, self%kq1)
   call spec%get_q2shift(mol, bas, self%kq2)

   allocate(self%hscale(mshell, mshell, mol%nid, mol%nid))
   call spec%get_hscale(mol, bas, self%hscale)

   allocate(self%shpoly(mshell, mol%nid), self%rad(mol%nid))
   call spec%get_rad(mol, bas, self%rad)
   call spec%get_shpoly(mol, bas, self%shpoly)

   allocate(self%refocc(mshell, mol%nid))
   call spec%get_reference_occ(mol, bas, self%refocc)

   allocate(self%ksig(mol%nid, mol%nid), self%kpi(mol%nid, mol%nid), &
      & self%kdel(mol%nid, mol%nid))
   call spec%get_diat_scale(mol, bas, self%ksig, self%kpi, self%kdel)

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


subroutine get_hamiltonian(mol, trans, list, bas, h0, selfenergy, overlap, dpint, qpint, &
      & hamiltonian)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Neighbour list
   type(adjacency_list), intent(in) :: list
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

   integer :: iat, jat, izp, jzp, itr, k, img, inl
   integer :: ish, jsh, is, js, ii, jj, iao, jao, nao, ij
   real(wp) :: rr, r2, vec(3), hij, shpoly, dtmpj(3), qtmpj(6)
   real(wp), allocatable :: stmp(:), dtmpi(:, :), qtmpi(:, :)

   overlap(:, :) = 0.0_wp
   dpint(:, :, :) = 0.0_wp
   qpint(:, :, :) = 0.0_wp
   hamiltonian(:, :) = 0.0_wp

   allocate(stmp(msao(bas%maxl)**2), dtmpi(3, msao(bas%maxl)**2), qtmpi(6, msao(bas%maxl)**2))

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, bas, trans, list, overlap, dpint, qpint, hamiltonian, h0, selfenergy) &
   !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj, iao, jao, nao, ij, k) &
   !$omp private(r2, vec, stmp, dtmpi, qtmpi, dtmpj, qtmpj, hij, shpoly, rr, inl, img)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      inl = list%inl(iat)
      do img = 1, list%nnl(iat)
         jat = list%nlat(img+inl)
         itr = list%nltr(img+inl)
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         rr = sqrt(sqrt(r2) / (h0%rad(jzp) + h0%rad(izp)))
         do ish = 1, bas%nsh_id(izp)
            ii = bas%iao_sh(is+ish)
            do jsh = 1, bas%nsh_id(jzp)
               jj = bas%iao_sh(js+jsh)
               call multipole_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                  & r2, vec, bas%intcut, stmp, dtmpi, qtmpi)

               shpoly = (1.0_wp + h0%shpoly(ish, izp)*rr) &
                  * (1.0_wp + h0%shpoly(jsh, jzp)*rr)

               hij = 0.5_wp * (selfenergy(is+ish) + selfenergy(js+jsh)) &
                  * h0%hscale(jsh, ish, jzp, izp) * shpoly

               nao = msao(bas%cgto(jsh, jzp)%ang)
               do iao = 1, msao(bas%cgto(ish, izp)%ang)
                  do jao = 1, nao
                     ij = jao + nao*(iao-1)
                     call shift_operator(vec, stmp(ij), dtmpi(:, ij), qtmpi(:, ij), &
                        & dtmpj, qtmpj)
                     !$omp atomic
                     overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                        + stmp(ij)

                     do k = 1, 3
                        !$omp atomic
                        dpint(k, jj+jao, ii+iao) = dpint(k, jj+jao, ii+iao) &
                           + dtmpi(k, ij)
                     end do

                     do k = 1, 6
                        !$omp atomic
                        qpint(k, jj+jao, ii+iao) = qpint(k, jj+jao, ii+iao) &
                           + qtmpi(k, ij)
                     end do

                     !$omp atomic
                     hamiltonian(jj+jao, ii+iao) = hamiltonian(jj+jao, ii+iao) &
                        + stmp(ij) * hij

                     if (iat /= jat) then
                        !$omp atomic
                        overlap(ii+iao, jj+jao) = overlap(ii+iao, jj+jao) &
                           + stmp(ij)
                        do k = 1, 3
                           !$omp atomic
                           dpint(k, ii+iao, jj+jao) = dpint(k, ii+iao, jj+jao) &
                              + dtmpj(k)
                        end do

                        do k = 1, 6
                           !$omp atomic
                           qpint(k, ii+iao, jj+jao) = qpint(k, ii+iao, jj+jao) &
                              + qtmpj(k)
                        end do
                        !$omp atomic
                        hamiltonian(ii+iao, jj+jao) = hamiltonian(ii+iao, jj+jao) &
                           + stmp(ij) * hij
                     end if
                  end do
               end do

            end do
         end do

      end do
   end do

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, bas, overlap, dpint, qpint, hamiltonian, h0, selfenergy) &
   !$omp private(iat, izp, is, ish, jsh, ii, jj, iao, jao, nao, ij) &
   !$omp private(vec, stmp, dtmpi, qtmpi, hij)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      vec(:) = 0.0_wp
      do ish = 1, bas%nsh_id(izp)
         ii = bas%iao_sh(is+ish)
         do jsh = 1, bas%nsh_id(izp)
            jj = bas%iao_sh(is+jsh)
            call multipole_cgto(bas%cgto(jsh, izp), bas%cgto(ish, izp), &
               & 0.0_wp, vec, bas%intcut, stmp, dtmpi, qtmpi)

            ! shpoly is always 1.0, because rr is always 0.0
            hij = 0.5_wp * (selfenergy(is+ish) + selfenergy(is+jsh))

            nao = msao(bas%cgto(jsh, izp)%ang)
            do iao = 1, msao(bas%cgto(ish, izp)%ang)
               do jao = 1, nao
                  ij = jao + nao*(iao-1)
                  overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                     + stmp(ij)

                  dpint(:, jj+jao, ii+iao) = dpint(:, jj+jao, ii+iao) &
                     + dtmpi(:, ij)

                  qpint(:, jj+jao, ii+iao) = qpint(:, jj+jao, ii+iao) &
                     + qtmpi(:, ij)

                  hamiltonian(jj+jao, ii+iao) = hamiltonian(jj+jao, ii+iao) &
                     + stmp(ij) * hij
               end do
            end do

         end do
      end do
   end do

end subroutine get_hamiltonian


subroutine get_hamiltonian_gradient(mol, trans, list, bas, h0, selfenergy, dsedcn, &
      & pot, pmat, xmat, dEdcn, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Neighbour list
   type(adjacency_list), intent(in) :: list
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
   real(wp), intent(in) :: pmat(:, :, :)
   !> Energy weighted density matrix
   real(wp), intent(in) :: xmat(:, :, :)

   !> Derivative of the electronic energy w.r.t. the coordination number
   real(wp), intent(inout) :: dEdcn(:)
   !> Derivative of the electronic energy w.r.t. coordinate displacements
   real(wp), intent(inout) :: gradient(:, :)
   !> Derivative of the electronic energy w.r.t. strain deformations
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, izp, jzp, itr, img, inl, spin, nspin
   integer :: ish, jsh, is, js, ii, jj, iao, jao, nao, ij
   real(wp) :: rr, r2, vec(3), hij, dG(3), hscale, hs
   real(wp) :: shpolyi, shpolyj, shpoly, dshpoly, dsv(3)
   real(wp) :: sval, dcni, dcnj, dhdcni, dhdcnj, hpij, pij
   real(wp), allocatable :: stmp(:), dtmp(:, :), qtmp(:, :)
   real(wp), allocatable :: dstmp(:, :), ddtmpi(:, :, :), dqtmpi(:, :, :)
   real(wp), allocatable :: ddtmpj(:, :, :), dqtmpj(:, :, :)

   nspin = size(pmat, 3)

   allocate(stmp(msao(bas%maxl)**2), dstmp(3, msao(bas%maxl)**2), &
      & dtmp(3, msao(bas%maxl)**2), ddtmpi(3, 3, msao(bas%maxl)**2), &
      & qtmp(6, msao(bas%maxl)**2), dqtmpi(3, 6, msao(bas%maxl)**2), &
      & ddtmpj(3, 3, msao(bas%maxl)**2), dqtmpj(3, 6, msao(bas%maxl)**2))

   !$omp parallel do schedule(runtime) default(none) reduction(+:dEdcn, gradient, sigma) &
   !$omp shared(mol, bas, trans, h0, selfenergy, dsedcn, pot, pmat, xmat, list, nspin) &
   !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj, iao, jao, nao, ij, spin, &
   !$omp& r2, vec, stmp, dtmp, qtmp, dstmp, ddtmpi, dqtmpi, ddtmpj, dqtmpj, hij, &
   !$omp& dG, dcni, dcnj, dhdcni, dhdcnj, hpij, rr, sval, hscale, pij, inl, img, &
   !$omp& hs, shpolyi, shpolyj, shpoly, dshpoly, dsv)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      inl = list%inl(iat)
      do img = 1, list%nnl(iat)
         jat = list%nlat(img+inl)
         itr = list%nltr(img+inl)
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         if (iat == jat) cycle
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         rr = sqrt(sqrt(r2) / (h0%rad(jzp) + h0%rad(izp)))
         do ish = 1, bas%nsh_id(izp)
            ii = bas%iao_sh(is+ish)
            do jsh = 1, bas%nsh_id(jzp)
               jj = bas%iao_sh(js+jsh)
               call multipole_grad_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                  & r2, vec, bas%intcut, stmp, dtmp, qtmp, dstmp, ddtmpj, dqtmpj, &
                  & ddtmpi, dqtmpi)
               
               shpolyi = 1.0_wp + h0%shpoly(ish, izp)*rr
               shpolyj = 1.0_wp + h0%shpoly(jsh, jzp)*rr
               shpoly = shpolyi * shpolyj
               dshpoly = (shpolyi * h0%shpoly(jsh, jzp) + shpolyj * &
                  & h0%shpoly(ish, izp)) * 0.5_wp * rr / r2
               dsv(:) = dshpoly / shpoly * vec

               hscale = h0%hscale(jsh, ish, jzp, izp)
               hs = hscale * shpoly
               hij = 0.5_wp * (selfenergy(is+ish) + selfenergy(js+jsh)) * hs 

               dhdcni = dsedcn(is+ish) * hs
               dhdcnj = dsedcn(js+jsh) * hs

               dG(:) = 0.0_wp
               dcni = 0.0_wp
               dcnj = 0.0_wp
               nao = msao(bas%cgto(jsh, jzp)%ang)
               do iao = 1, msao(bas%cgto(ish, izp)%ang)
                  do jao = 1, nao
                     ij = jao + nao*(iao-1)
                     do spin = 1, nspin
                        pij = pmat(jj+jao, ii+iao, spin)
                        sval = - pij * (pot%vao(jj+jao, spin) + pot%vao(ii+iao, spin))

                        dG(:) = dG + sval * dstmp(:, ij)  
                     end do
                     pij = pmat(jj+jao, ii+iao, 1)
                     hpij = pij * hij 
                     sval = 2*hpij - 2*xmat(jj+jao, ii+iao, 1)

                     dG(:) = dG + sval * dstmp(:, ij) &
                        + 2*hpij*stmp(ij) * dsv &
                        - pij * matmul(ddtmpi(:, :, ij), pot%vdp(:, iat, 1)) &
                        - pij * matmul(ddtmpj(:, :, ij), pot%vdp(:, jat, 1)) &
                        - pij * matmul(dqtmpi(:, :, ij), pot%vqp(:, iat, 1)) &
                        - pij * matmul(dqtmpj(:, :, ij), pot%vqp(:, jat, 1))

                     dcni = dcni + dhdcni * pmat(jj+jao, ii+iao, 1) * stmp(ij)
                     dcnj = dcnj + dhdcnj * pmat(jj+jao, ii+iao, 1) * stmp(ij)
                  end do
               end do
               dEdcn(iat) = dEdcn(iat) + dcni
               dEdcn(jat) = dEdcn(jat) + dcnj
               gradient(:, iat) = gradient(:, iat) + dG
               gradient(:, jat) = gradient(:, jat) - dG
               sigma(:, :) = sigma + 0.5_wp * (spread(vec, 1, 3) * spread(dG, 2, 3) &
                  + spread(dG, 1, 3) * spread(vec, 2, 3))

            end do
         end do

      end do
   end do

   !$omp parallel do schedule(runtime) default(none) reduction(+:dEdcn) &
   !$omp shared(mol, bas, dsedcn, pmat, nspin) &
   !$omp private(iat, izp, jzp, is, ish, ii, iao, dcni, dhdcni, spin)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      do ish = 1, bas%nsh_id(izp)
         ii = bas%iao_sh(is+ish)
         dhdcni = dsedcn(is+ish)
         dcni = 0.0_wp
         do iao = 1, msao(bas%cgto(ish, izp)%ang)
            dcni = dcni + dhdcni * pmat(ii+iao, ii+iao, 1)
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


!> Shift multipole operator from Ket function (center i) to Bra function (center j),
!> the multipole operator on the Bra function can be assembled from the lower moments
!> on the Ket function and the displacement vector using horizontal shift rules.
pure subroutine shift_operator(vec, s, di, qi, dj, qj)
   !> Displacement vector of center i and j
   real(wp),intent(in) :: vec(:)
   !> Overlap integral between basis functions
   real(wp),intent(in) :: s
   !> Dipole integral with operator on Ket function (center i)
   real(wp),intent(in) :: di(:)
   !> Quadrupole integral with operator on Ket function (center i)
   real(wp),intent(in) :: qi(:)
   !> Dipole integral with operator on Bra function (center j)
   real(wp),intent(out) :: dj(:)
   !> Quadrupole integral with operator on Bra function (center j)
   real(wp),intent(out) :: qj(:)

   real(wp) :: tr

   ! Create dipole operator on Bra function from Ket function and shift contribution
   ! due to monopol displacement
   dj(1) = di(1) + vec(1)*s
   dj(2) = di(2) + vec(2)*s
   dj(3) = di(3) + vec(3)*s

   ! For the quadrupole operator on the Bra function we first construct the shift
   ! contribution from the dipole and monopol displacement, since we have to remove
   ! the trace contribution from the shift and the moment integral on the Ket function
   ! is already traceless
   qj(1) = 2*vec(1)*di(1) + vec(1)**2*s
   qj(3) = 2*vec(2)*di(2) + vec(2)**2*s
   qj(6) = 2*vec(3)*di(3) + vec(3)**2*s
   qj(2) = vec(1)*di(2) + vec(2)*di(1) + vec(1)*vec(2)*s
   qj(4) = vec(1)*di(3) + vec(3)*di(1) + vec(1)*vec(3)*s
   qj(5) = vec(2)*di(3) + vec(3)*di(2) + vec(2)*vec(3)*s
   ! Now collect the trace of the shift contribution
   tr = 0.5_wp * (qj(1) + qj(3) + qj(6))

   ! Finally, assemble the quadrupole operator on the Bra function from the operator
   ! on the Ket function and the traceless shift contribution
   qj(1) = qi(1) + 1.5_wp * qj(1) - tr
   qj(2) = qi(2) + 1.5_wp * qj(2)
   qj(3) = qi(3) + 1.5_wp * qj(3) - tr
   qj(4) = qi(4) + 1.5_wp * qj(4)
   qj(5) = qi(5) + 1.5_wp * qj(5)
   qj(6) = qi(6) + 1.5_wp * qj(6) - tr
end subroutine shift_operator


end module tblite_xtb_h0
