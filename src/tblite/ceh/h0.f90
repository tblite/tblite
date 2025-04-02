
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

!> @file tblite/ceh/h0.f90
!> Provides the Hamiltonian level functions for CEH.

module tblite_ceh_h0
   use mctc_env, only : wp
   use mctc_io, only: structure_type
   use tblite_basis_type, only:  basis_type
   use tblite_xtb_spec, only : tb_h0spec
   use tblite_xtb_h0, only : tb_hamiltonian
   use tblite_integral_dipole, only: get_dipole_integrals, dipole_cgto, &
   & maxl, msao, smap
   use tblite_adjlist, only : adjacency_list
   use tblite_integral_diat_trafo, only: diat_trafo

   implicit none
   private

   public :: get_scaled_selfenergy, get_hamiltonian, get_occupation

contains


   !> Scale the selfenergy parameters by the CN
   subroutine get_scaled_selfenergy(h0, id, ish_at, nshell, cn, cn_en, selfenergy, dsedcn, dsedcn_en)
      type(tb_hamiltonian), intent(in) :: h0
      integer, intent(in) :: id(:)
      integer, intent(in) :: ish_at(:)
      integer, intent(in) :: nshell(:)
      real(wp), intent(in), optional :: cn(:)
      real(wp), intent(in), optional :: cn_en(:)
      real(wp), intent(out) :: selfenergy(:)
      real(wp), intent(out), optional :: dsedcn(:)
      real(wp), intent(out), optional :: dsedcn_en(:)

      integer :: iat, izp, ish, ii

      selfenergy(:) = 0.0_wp
      if (present(dsedcn)) dsedcn(:) = 0.0_wp
      if (present(dsedcn_en)) dsedcn_en(:) = 0.0_wp
      do iat = 1, size(id)
         izp = id(iat)
         ii = ish_at(iat)
         do ish = 1, nshell(izp)
            selfenergy(ii+ish) = h0%selfenergy(ish, izp)
         end do
      end do
      if (.not.present(cn) .or. .not.present(cn_en)) return
      if (present(dsedcn) .and. present(dsedcn_en)) then
         do iat = 1, size(id)
            izp = id(iat)
            ii = ish_at(iat)
            do ish = 1, nshell(izp)
               selfenergy(ii+ish) = selfenergy(ii+ish) + h0%kcn(ish, izp) * cn(iat) + &
               &  h0%kcn_en(ish, izp) * cn_en(iat)
               dsedcn(ii+ish) = h0%kcn(ish, izp)
               dsedcn_en(ii+ish) = h0%kcn_en(ish, izp)
            end do
         end do
      else
         do iat = 1, size(id)
            izp = id(iat)
            ii = ish_at(iat)
            do ish = 1, nshell(izp)
               selfenergy(ii+ish) = selfenergy(ii+ish) + h0%kcn(ish, izp) * cn(iat) + &
               &  h0%kcn_en(ish, izp) * cn_en(iat)
            end do
         end do
      end if

   end subroutine get_scaled_selfenergy


   subroutine get_hamiltonian(mol, trans, list, bas, h0, selfenergy, &
   & overlap, overlap_diat, dpint, hamiltonian)
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
      !> Diatomic frame-scaled overlap integral matrix
      real(wp), intent(out) :: overlap_diat(:, :)
      !> Effective Hamiltonian
      real(wp), intent(out) :: hamiltonian(:, :)

      integer  :: itr, img, inl, ii, jj, is, js, jzp, izp, nao
      integer  :: iat, ish, jat, jsh, k, iao, jao, ij, iaosh, jaosh
      real(wp) :: hij, rr, r2, vec(3), dtmpj(3)
      real(wp), allocatable :: stmp(:), dtmpi(:, :), block_overlap(:,:)

      overlap(:, :) = 0.0_wp
      dpint(:, :, :) = 0.0_wp
      overlap_diat(:, :) = 0.0_wp
      hamiltonian(:, :) = 0.0_wp

      ! Allocate temporary matrices for overlap, diatomic frame scaled overlap and dipole moment integrals
      allocate(stmp(msao(bas%maxl)**2), dtmpi(3, msao(bas%maxl)**2), &
      & block_overlap(smap(bas%maxl+1),smap(bas%maxl+1)))

      !$omp parallel do schedule(runtime) default(none) &
      !$omp shared(mol, bas, trans, list, overlap, overlap_diat, dpint, hamiltonian, h0, selfenergy) &
      !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj, iao, jao, nao, ij, iaosh, jaosh, k) &
      !$omp private(r2, vec, stmp, block_overlap, dtmpi, dtmpj, hij, rr, inl, img)
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

            ! Get the overlap and dipole integrals for the current diatomic pair
            block_overlap = 0.0_wp
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is+ish)
               iaosh = smap(ish-1) ! Offset for the block overlap matrix
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js+jsh)
                  jaosh = smap(jsh-1) ! Offset for the block overlap matrix

                  call dipole_cgto(bas%cgto(jsh,jzp), bas%cgto(ish,izp), &
                  & r2, vec, bas%intcut, stmp, dtmpi)

                  ! Store the overlap and dipole matrix
                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)
                        ! Shift dipole operator from Ket function (center i)
                        ! to Bra function (center j) to save the redundant calculation
                        call shift_operator(vec, stmp(ij), dtmpi(:, ij), dtmpj)

                        !$omp atomic
                        overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                           + stmp(ij)

                        !$omp atomic
                        block_overlap(jaosh+jao, iaosh+iao) = block_overlap(jaosh+jao, iaosh+iao) &
                           + stmp(ij)

                        do k = 1, 3
                           !$omp atomic
                           dpint(k, jj+jao, ii+iao) = dpint(k, jj+jao, ii+iao) &
                              + dtmpi(k, ij)
                        end do

                        if (iat /= jat) then
                           !$omp atomic
                           overlap(ii+iao, jj+jao) = overlap(ii+iao, jj+jao) &
                              + stmp(ij)

                           do k = 1, 3
                              !$omp atomic
                              dpint(k, ii+iao, jj+jao) = dpint(k, ii+iao, jj+jao) &
                                 + dtmpj(k)
                           end do

                        end if
                     end do
                  end do
               end do
            end do
   
            ! Diatomic frame transformation and scaling of the overlap
            call diat_trafo(block_overlap, vec, h0%ksig(izp,jzp), h0%kpi(izp,jzp), h0%kdel(izp,jzp), &
               & bas%nsh_at(jat)-1, bas%nsh_at(iat)-1) 

            ! Setup the Hamiltonian and store the diatomic frame scaled overlap. 
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is+ish)
               iaosh = smap(ish-1) ! Offset for the block overlap matrix
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js+jsh)
                  jaosh = smap(jsh-1) ! Offset for the block overlap matrix
                  
                  hij = 0.5_wp * h0%hscale(jsh, ish, jzp, izp) * (selfenergy(is+ish) + selfenergy(js+jsh)) 

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)

                        !$omp atomic
                        overlap_diat(jj+jao, ii+iao) = overlap_diat(jj+jao, ii+iao) &
                           + block_overlap(jaosh+jao, iaosh+iao)

                        !$omp atomic
                        hamiltonian(jj+jao, ii+iao) = hamiltonian(jj+jao, ii+iao) &
                           + block_overlap(jaosh+jao, iaosh+iao) * hij

                        if (iat /= jat) then
                           !$omp atomic
                           overlap_diat(ii+iao, jj+jao) = overlap_diat(ii+iao, jj+jao) &
                              + block_overlap(jaosh+jao, iaosh+iao)

                           !$omp atomic
                           hamiltonian(ii+iao, jj+jao) = hamiltonian(ii+iao, jj+jao) &
                              + block_overlap(jaosh+jao, iaosh+iao) * hij
                        end if
                     end do
                  end do
               end do
            end do

         end do
      end do

      !$omp parallel do schedule(runtime) default(none) &
      !$omp shared(mol, bas, trans, overlap, overlap_diat, dpint, hamiltonian, h0, selfenergy) &
      !$omp private(iat, izp, itr, is, ish, jsh, ii, jj, iao, jao, nao, ij) &
      !$omp private(r2, vec, stmp, dtmpi, hij, rr)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         is = bas%ish_at(iat)
         vec(:) = 0.0_wp
         r2 = 0.0_wp
         rr = sqrt(sqrt(r2) / (h0%rad(izp) + h0%rad(izp)))
         do ish = 1, bas%nsh_id(izp)
            ii = bas%iao_sh(is+ish)
            do jsh = 1, bas%nsh_id(izp)
               jj = bas%iao_sh(is+jsh)
               call dipole_cgto(bas%cgto(jsh,izp), bas%cgto(ish,izp), &
               & r2, vec, bas%intcut, stmp, dtmpi)
               nao = msao(bas%cgto(jsh, izp)%ang)
               do iao = 1, msao(bas%cgto(ish, izp)%ang)
                  do jao = 1, nao
                     ij = jao + nao*(iao-1)
                     overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                        + stmp(ij)

                     overlap_diat(jj+jao, ii+iao) = overlap_diat(jj+jao, ii+iao) &
                        + stmp(ij)

                     dpint(:, jj+jao, ii+iao) = dpint(:, jj+jao, ii+iao) &
                        + dtmpi(:, ij)
                  end do
               end do
            end do
            ! diagonal term (AO(i) == AO(j)) as on-site off-diagonal the hamiltonian is zero
            do iao = 1, msao(bas%cgto(ish, izp)%ang)
               hamiltonian(ii+iao, ii+iao) = selfenergy(is + ish)
            enddo
         end do
      end do

   end subroutine get_hamiltonian


   subroutine get_occupation(mol, bas, hamiltonian, nocc, n0at, n0sh)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set information
      type(basis_type), intent(in) :: bas
      !> Hamiltonian interaction data
      type(tb_hamiltonian), intent(in) :: hamiltonian
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
            nocc = nocc + hamiltonian%refocc(ish, izp)
            n0at(iat) = n0at(iat) + hamiltonian%refocc(ish, izp)
            n0sh(ii+ish) = n0sh(ii+ish) + hamiltonian%refocc(ish, izp)
         end do
      end do
   end subroutine get_occupation


   !> Shift dipole operator from Ket function (center i) to Bra function (center j),
   !> the dipole operator on the Bra function can be assembled from the lower moments
   !> on the Ket function and the displacement vector using horizontal shift rules.
   pure subroutine shift_operator(vec, s, di, dj)
      !> Displacement vector of center i and j
      real(wp),intent(in) :: vec(:)
      !> Overlap integral between basis functions
      real(wp),intent(in) :: s
      !> Dipole integral with operator on Ket function (center i)
      real(wp),intent(in) :: di(:)
      !> Dipole integral with operator on Bra function (center j)
      real(wp),intent(out) :: dj(:)
      ! Create dipole operator on Bra function from Ket function and shift contribution
      ! due to monopol displacement
      dj(1) = di(1) + vec(1)*s
      dj(2) = di(2) + vec(2)*s
      dj(3) = di(3) + vec(3)*s
   end subroutine shift_operator

end module tblite_ceh_h0
