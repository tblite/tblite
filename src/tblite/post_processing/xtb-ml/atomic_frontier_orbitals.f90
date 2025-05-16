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

!> @file tblite/post-processing/xtb-ml/atomic_frontier_orbitals.f90

!> compute atomic response function and effective H-L gap
module tblite_xtbml_atomic_frontier
   use mctc_env, only : wp, sp
   implicit none
   private
   public :: atomic_frontier_orbitals

contains

!> compute atomic response function and atomic H-L gap and Fermi level
subroutine atomic_frontier_orbitals(focc, eps, aoat, C, S, response, egap, chempot, ehoao, eluao)
   !> get atom index of AO orbital
   integer, intent(in) :: aoat(:)
   !> occupation numbers
   real(wp), intent(in) :: focc(:)
   !> orbital energies
   real(wp), intent(in) :: eps(:)
   !> MO coefficients
   real(wp), intent(in) :: C(:,:)
   !> overlap matrix
   real(wp), intent(in) :: S(:,:)
   !> atomic response
   real(wp), intent(out) :: response(:)
   !> effective atomic H-L gap
   real(wp), intent(out) :: egap(:)
   !> effective Fermi level per atom
   real(wp), intent(out) :: chempot(:)
   !> highest occ. atomic orbital
   real(wp), intent(out) :: ehoao(:)
   !> lowest virt. atomic orbital
   real(wp), intent(out) :: eluao(:)

   integer :: i, j, k, jj, kk, m
   real(wp), allocatable :: po(:,:), pv(:,:)
   real(wp) :: occ, tmp, tmp2, ps, virt, tmp3, weight
   ! occupation cutoff for fractional occupation, below this threshold orbitals are not considered
   real(wp), parameter :: occ_cutoff = 0.0001_wp
   ! theshold below which we consider the sum or virtual/ocuupied orbitals as zero
   real(wp), parameter :: zero_cutoff = 1.0e-12_wp
   real(wp), parameter :: damp = 0.5_wp ! damping in response function (in eV)
   ! value used when either virtual or occupied space is nearly unoccupied
   real(wp), parameter :: near_infty = 1.0e100_wp
   real(wp), parameter :: epsilon = 1.0e-14_wp
   integer :: nat, nao

   ! Thread-private arrays for reduction
   real(wp), allocatable :: ptmp_local(:, :)
   real(wp), allocatable :: response_local(:), egap_local(:), chempot_local(:)

   nat = size(response)
   nao = size(focc)

   !allocate temp arrays
   allocate(po(nat, nao), source=0.0_wp)
   allocate(pv(nat, nao), source=0.0_wp)

   ! we collect occ & virt separately (and allow for fractional occupation)
   !$omp parallel default(none) &
   !$omp shared(aoat, nao, nat, S, C, focc, po, pv) &
   !$omp private(ps, occ, virt, jj, kk, j, i, k, ptmp_local)
   allocate(ptmp_local(nat, nao), source=0.0_wp)
   !$omp do schedule(runtime) 
   do i = 1, nao
      ! occupied part
      occ = focc(i)
      if (occ .gt. occ_cutoff) then
         do j = 1, nao
            jj = aoat(j)
            do k = 1, j - 1
               kk = aoat(k)
               ps = s(k, j) * C(j, i) * C(k, i) * occ
               ptmp_local(kk, i) = ptmp_local(kk, i) + ps
               ptmp_local(jj, i) = ptmp_local(jj, i) + ps
            end do
            ps = s(j, j) * C(j, i) * C(j, i) * occ
            ptmp_local(jj, i) = ptmp_local(jj, i) + ps
         end do
      end if
   end do
   !$omp end do
   !$omp critical (atomic_frontier_orbitals_)
   po(:, :) = po(:, :) + ptmp_local(:, :)
   !$omp end critical (atomic_frontier_orbitals_)

   ptmp_local(:, :) = 0.0_wp
   !$omp do schedule(runtime)
   do i = 1, nao
      ! virtual part
      virt = (1.0_wp - focc(i))
      if (virt .gt. occ_cutoff) then
         do j = 1, nao
            jj = aoat(j)
            do k = 1, j - 1
               kk = aoat(k)
               ps = s(k, j) * C(j, i) * C(k, i) * virt
               ptmp_local(kk, i) = ptmp_local(kk, i) + ps
               ptmp_local(jj, i) = ptmp_local(jj, i) + ps
            end do
            ps = s(j, j)*C(j, i)*C(j, i)*virt
            ptmp_local(jj, i) = ptmp_local(jj, i) + ps
         end do
      end if
   end do
   !$omp critical (atomic_frontier_orbitals_)
   pv(:, :) = pv(:, :) + ptmp_local(:, :)
   !$omp end critical (atomic_frontier_orbitals_)
   !$omp end parallel

   ! now accumulate the atomic H-L gaps and Fermi levels (we approximate it as 0.5*(eH + eL))
   ! we make use of an atomic response-type weightinhg:
   ! chempot=\sum_ai * wA_ai * 0.5 * (e_a + e_i)
   ! here wA_ai is the variable "response" computed as: wA_ai=  [p_i p_a / ( (e_a - e_i)**2 + damp**2) ]/ [\sum_ia p_j p_b / ((e_b - e_j)**2 + damp**2)]
   ! the "p"s are the MO densities, for gaps we accumulate the regularized inverse and later on invert it again

   response(:) = 0.0_wp
   egap(:) = 0.0_wp
   chempot(:) = 0.0_wp
   ! go through occ (including fractionally occupied)
   !$omp parallel default(none) & 
   !$omp shared(eps, nao, focc, po, pv, nat, response, egap, chempot) &
   !$omp private(occ, virt, tmp, tmp2, tmp3, weight, i, j, m, & 
   !$omp& response_local, egap_local, chempot_local)
   allocate(response_local, source=response)
   allocate(egap_local, source=egap)
   allocate(chempot_local, source=chempot)
   !$omp do schedule(runtime)
   do i = 1, nao
      occ = focc(i)
      if (occ .gt. occ_cutoff) then
         !virt part  (including fractionally occupied)
         do j = 1, nao
            virt = (1.0_wp - focc(j))
            if (virt .gt. occ_cutoff) then
               tmp = 1.0_wp / ((eps(j) - eps(i))**2 + damp**2) !divsion by zero not possible due to damping
               tmp2 = 0.5_wp * (eps(j) + eps(i))
               tmp3 = 1.0_wp / (eps(j) - eps(i) + damp) !divsion by zero not possible due to damping
               ! compute atomic response
               do m = 1, nat
                  weight = po(m, i) * pv(m, j) * tmp
                  response_local(m) = response_local(m) + weight
                  chempot_local(m) = chempot_local(m) + tmp2 * weight
                  egap_local(m) = egap_local(m) + weight * tmp3
               end do ! at
            end if ! if virt
         end do !  virt
      end if ! if occ
   end do ! occ
   !$omp end do
   !$omp critical (atomic_frontier_orbitals_)
   response(:) = response(:) + response_local(:)
   egap(:) = egap(:) + egap_local(:)
   chempot(:) = chempot(:) + chempot_local(:)
   !$omp end critical (atomic_frontier_orbitals_)
   deallocate(response_local, egap_local, chempot_local)
   !$omp end parallel

   ehoao(:) = 0.0_wp
   eluao(:) = 0.0_wp

   ! now get the atomic frontier orbitals
   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(eps, focc, nat, nao, po, pv, response, egap, chempot, ehoao, eluao) &
   !$omp private(m, j, occ, virt, tmp, tmp2, tmp3, weight)
   do m = 1, nat
      egap(m) = egap(m)/(response(m)+epsilon)
      egap(m) = 1.0_wp / (egap(m)+epsilon) - damp
      chempot(m) = chempot(m)/(response(m)+epsilon)
      ehoao(m) = (chempot(m) - 0.5_wp*egap(m))
      eluao(m) = (chempot(m) + 0.5_wp*egap(m))

      ! fix cases for missing HOAO/LUAO
      if (sum(po(m, :)) .lt. zero_cutoff) then ! there is no occupied orbital for this atom
         !virt part  (including fractionally occupied)
         egap(m) = 0.0_wp
         chempot(m) = 0.0_wp
         do j = 1, nao
            virt = (1.0_wp - focc(j))
            if (virt .gt. occ_cutoff) then
               tmp = 1.0_wp / ((eps(j))**2 + damp**2)
               tmp2 = eps(j)
               tmp3 = 1.0_wp / (eps(j) + near_infty + damp)
               ! compute atomic response
               weight = pv(m, j) * tmp
               chempot(m) = chempot(m) + tmp2 * weight
               egap(m) = egap(m) + weight * tmp3
            end if ! if virt
            egap(m) = 1.0_wp / (egap(m) + zero_cutoff) - damp
            eluao(m) = chempot(m)
            ehoao(m) = chempot(m) - egap(m)
         end do !  virt
      end if
      if (sum(pv(m, :)) .lt. zero_cutoff) then
         !occ part  (including fractionally occupied)
         egap(m) = 0.0_wp
         chempot(m) = 0.0_wp
         do j = 1, nao
            occ = focc(j)
            if (occ .gt. occ_cutoff) then
               tmp = 1.0_wp / ((eps(j))**2 + damp**2)
               tmp2 = eps(j)
               tmp3 = 1.0_wp / (near_infty - eps(j) + damp)
               ! compute atomic response
               weight = po(m, j) * tmp
               chempot(m) = chempot(m) + tmp2 * weight
               egap(m) = egap(m) + weight * tmp3
            end if ! if occ
            egap(m) = 1.0_wp/(egap(m) + zero_cutoff) - damp
            ehoao(m) = chempot(m)
            eluao(m) = chempot(m) + egap(m)
         end do !  occ
      end if
   end do
   !$omp end parallel do

end subroutine atomic_frontier_orbitals
end module
