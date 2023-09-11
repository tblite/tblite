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

!> @file tblite/xtb-ml/functions.f90
module tblite_xtbml_functions
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_results, only : results_type
   use tblite_blas
   implicit none
   private

   real(wp), allocatable :: inv_cn_a(:,:,:)
   real(wp), allocatable :: rcov(:)
   real(wp), parameter :: k1 = 16.0

contains


subroutine get_total_xtb_weights(nat, e_1e, e_rep, e_diff_atom, e_disp_3, e_ies_ixc, e_aes, e_axc, w_tot_xtb, e_tot)
   integer, INTENT(IN) :: nat
   real(wp), INTENT(IN) :: e_1e(nat), e_diff_atom(nat), e_disp_3(nat), e_ies_ixc(nat), e_aes(nat), e_axc(nat), e_rep(nat)
   real(wp), INTENT(OUT) :: w_tot_xtb(nat), e_tot
   real(wp) :: sum_energy, sum_atom(nat)

   sum_atom = e_1e + e_rep + e_diff_atom + e_disp_3 + e_ies_ixc + e_aes + e_axc
   sum_energy = sum(sum_atom(:))

   w_tot_xtb(:) = sum_atom/sum_energy
   e_tot = sum_energy

end subroutine

subroutine get_beta(nat, n_a, at, xyz, beta)
   implicit none
   intrinsic :: SUM
   integer, INTENT(IN) :: nat, at(nat), n_a
   real(wp), INTENT(IN) :: xyz(3, nat)
   real(wp) :: beta(nat, nat, n_a)

   real(wp) :: sigma_tot, sigma(nat, nat, n_a)
   real(wp) :: damp_func

   integer :: A, B, k
   do k = 1, n_a
      do A = 1, nat
         do B = 1, nat
            damp_func = inv_cn_a(A, B, k)
            sigma(A, B, k) = 1/damp_func
         end do
      end do
   end do

   do k = 1, n_a
      do A = 1, nat
         do B = 1, nat
            beta(A, B, k) = sigma(A, B, k)/sum(sigma(A, :, k))
         end do
      end do
   end do
end subroutine

subroutine get_chem_pot_ext(nat, n_a, beta, chempot, chempot_ext)
   implicit none
   integer, intent(in) :: nat, n_a
   real(wp), intent(in) :: beta(nat, nat, n_a)
   real(wp), intent(in) :: chempot(nat)
   real(wp) :: chempot_ext(nat, n_a)
   integer :: A, B, k
   do k = 1, n_a
      do A = 1, nat
         do B = 1, nat

            chempot_ext(A, k) = chempot_ext(A, k) + beta(A, B, k)*chempot(B)
         end do
      end do
   end do
end subroutine

subroutine get_e_gap_ext(nat, n_a, hl_gap, beta, e_gap, e_gap_ext)
   implicit none
   integer, intent(in) :: nat, n_a
   real(wp), INTENT(IN) :: hl_gap
   real(wp), intent(in) :: beta(nat, nat, n_a)
   real(wp), intent(in) :: e_gap(nat)
   real(wp) :: e_gap_ext(nat, n_a), e_gap_tot
   integer :: A, B, k

   e_gap_tot = 0.0_wp
   do k = 1, n_a
      do A = 1, nat
         e_gap_tot = e_gap_tot + e_gap(A)
         do B = 1, nat
            e_gap_ext(A, k) = e_gap_ext(A, k) + beta(A, B, k)*e_gap(B)
         end do
      end do
   end do
   ! correction factor for mol. properties
   !e_gap_ext = e_gap_ext * (hl_gap * nat/e_gap_tot)
end subroutine

subroutine get_ehoao_ext(nat, n_a, chempot_ext, e_gap_ext, ehoao_ext)
   implicit none
   integer, intent(in) :: nat, n_a
   real(wp), intent(in) :: chempot_ext(nat, n_a)
   real(wp), intent(in) :: e_gap_ext(nat, n_a)
   real(wp), intent(out) :: ehoao_ext(nat, n_a)
   integer :: A, k
   do k = 1, n_a
      do A = 1, nat
         ehoao_ext(A, k) = chempot_ext(A, k) - e_gap_ext(A, k)/2
      end do
   end do
end subroutine

subroutine get_eluao_ext(nat, n_a, chempot_ext, e_gap_ext, eluao_ext)
   implicit none
   integer, intent(in) :: nat, n_a
   real(wp), intent(in) :: chempot_ext(nat, n_a)
   real(wp), intent(in) :: e_gap_ext(nat, n_a)
   real(wp), intent(out) :: eluao_ext(nat, n_a)
   integer :: A, k
   do k = 1, n_a
      do A = 1, nat
         eluao_ext(A, k) = chempot_ext(A, k) + e_gap_ext(A, k)/2
      end do
   end do
end subroutine

end module tblite_xtbml_functions
