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
   real(wp), allocatable :: rcov(:), inv_cn_a(:, :, :)

   real(wp), parameter :: k1 = 16.0_wp

contains

subroutine populate_inv_cn_array(nat, at, xyz, a)
   integer, intent(in):: nat, at(nat)
   real(wp), intent(in) :: a(:), xyz(:, :)
   real(wp) :: result
   integer :: i, j, k
   n_a = size(a)
   if (allocated(inv_cn_a)) then
      deallocate (inv_cn_a)
   end if
   allocate (inv_cn_a(nat, nat, n_a), source=0.0_wp)
   !$omp parallel do default(none) collapse(2)&
   !$omp shared(a,nat, at,xyz,inv_cn_a,n_a)&
   !$omp private(result,i,j,k)
   do k = 1, n_a
      do i = 1, nat
         do j = 1, nat
            !if (i == j) cycle
            call inv_cn(nat, i, j, at, xyz, a(k), result)
            inv_cn_a(i, j, k) = result
         end do
      end do
   end do
   !$omp end parallel do

end subroutine populate_inv_cn_array

subroutine get_rcov(mol)
   use tblite_data_covrad, only : get_covalent_rad
   type(structure_type), intent(in) :: mol
   if (allocated(rcov)) then
      deallocate (rcov)
   end if
   allocate (rcov(mol%nid), source=0.0_wp)
   rcov(:) = get_covalent_rad(mol%num)
end subroutine

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

subroutine get_delta_cn(nat, n_a, cn, at, xyz, delta_cn)
   use tblite_timer, only : timer_type, format_time
   implicit none
   integer, INTENT(IN) :: nat, at(nat), n_a
   real(wp), INTENT(IN) :: cn(nat), xyz(3, nat)
   real(wp), INTENT(OUT) :: delta_cn(nat, n_a)
   real(wp):: delta_cn_tmp(nat, n_a), stime
   integer :: i, j, k
   type(timer_type) :: timer
   
   !$acc enter data create(delta_cn(:, :))
   !$acc kernels default(present)
   delta_cn = 0.0_wp
   !$acc end kernels
   !$acc enter data copyin( nat,n_a,cn(:), inv_cn_a(:,:,:))

   !$acc loop gang vector collapse(3)
   do k = 1, n_a
      do i = 1, nat
         do j = 1, nat
            if (i == j) cycle
            !$acc atomic
            delta_cn(i, k) = delta_cn(i, k) + cn(j)/inv_cn_a(i, j, k)

         end do
      end do
   end do

   !$acc end parallel
   !$acc exit data copyout(delta_cn(:, :))
end subroutine

subroutine compute_cn(nat, cn)
   real(wp) :: cn(nat)
   integer :: nat

   real(wp) :: mat(nat)

   mat = 1.0_wp
   call gemv(1.0_wp/inv_cn_a(:, :, 1), mat, cn)
   cn = cn - mat

end subroutine

subroutine inv_cn(nat, a, b, at, xyz, dampening_fact, result)
   implicit none
   integer, INTENT(IN) :: a, b, nat
   integer, INTENT(IN) :: at(:)
   real(wp), intent(in)  :: xyz(3, nat), dampening_fact
   real(wp), INTENT(OUT) :: result
   real(wp) :: rab(3), r, rco, r2

   result = 0.0_wp

   rab = xyz(1:3, a) - xyz(1:3, b)
   r2 = sum(rab**2)
   r = sqrt(r2)

   rco = dampening_fact*(rcov(at(a)) + rcov(at(b)))

   result = 1.0_wp/exp_count(k1, r, rco)

end subroutine

pure elemental function exp_count(k, r, r0) result(count)
   real(wp), intent(in) :: k
   real(wp), intent(in) :: r
   real(wp), intent(in) :: r0
   real(wp) :: count
   count = 1.0_wp/(1.0_wp + exp(-k*(r0/r - 1.0_wp)))
end function exp_count

subroutine sum_up_spin(q_2, q_1)
   real(wp), intent(in):: q_2(:, :)
   real(wp), intent(out) :: q_1(:)
   integer :: i
   do i = 1, size(q_2, 1)
      q_1(i) = sum(q_2(i, :))
   end do
end subroutine sum_up_spin

subroutine mulliken_shellwise(nao, nshell, ao2shell, p, s, charges_shell)
   implicit none
   integer, intent(in) :: nao, nshell, ao2shell(:)
   real(wp), intent(in) :: s(:, :)
   real(wp), intent(in) :: p(:, :)
   real(wp), intent(out) :: charges_shell(nshell)
   integer ::  mu, nu

   charges_shell = 0.0_wp
   !$omp parallel do default(none) collapse(2)&
   !$omp private(mu,nu) shared(charges_shell,nao,ao2shell,p,s)
   do mu = 1, nao
      do nu = 1, nao
         !$omp atomic
         charges_shell(ao2shell(mu)) = charges_shell(ao2shell(mu)) + p(mu, nu)*s(nu, mu)
      end do
   end do

end subroutine

subroutine get_delta_partial(nat, n_a, atom_partial, at, xyz, cn, delta_partial)
   implicit none
   integer, INTENT(IN) :: nat, at(nat), n_a
   real(wp), INTENT(IN) :: atom_partial(nat), xyz(3, nat), cn(nat)
   real(wp), INTENT(OUT) :: delta_partial(nat, n_a)
   real(wp) :: delta_partial_tmp(nat, n_a), f_log(nat, nat, n_a)
   integer :: a, b, k
   real(wp) :: result
   f_log = 0.0_wp

   delta_partial = 0.0_wp

   do k = 1, n_a
      do a = 1, nat
         do b = 1, nat
            delta_partial(a, k) = delta_partial(a, k) + atom_partial(b)/(inv_cn_a(a, b, k)*(cn(b) + 1))
            f_log(a, b, k) = 1.0_wp/(inv_cn_a(a, b, k)*(cn(b) + 1))
         end do
      end do
   end do
end subroutine

subroutine sum_up_mm(nat, nshell, aoat2, ash, dipm_shell, qm_shell, dipm_at, qm_at)
   implicit none
   integer, INTENT(IN) :: nat, nshell, aoat2(:), ash(:)
   real(wp), INTENT(IN) :: dipm_shell(:, :), qm_shell(:, :)
   real(wp), INTENT(OUT) :: dipm_at(3, nat), qm_at(6, nat)
   integer :: i

   dipm_at = 0.0_wp
   qm_at = 0.0_wp

   do i = 1, nshell
      dipm_at(:, ash(i)) = dipm_at(:, ash(i)) + dipm_shell(:, i)
      qm_at(:, ash(i)) = qm_at(:, ash(i)) + qm_shell(:, i)
   end do

end subroutine

subroutine get_delta_mm(nat, n_a, q, dipm, qp, at, xyz, cn, delta_dipm, delta_qp)
   implicit none
   integer, INTENT(IN) :: nat, at(nat), n_a
   real(wp), INTENT(IN) :: dipm(3, nat), xyz(3, nat), qp(6, nat), q(nat), cn(nat)
   real(wp), INTENT(OUT) :: delta_dipm(3, nat, n_a), delta_qp(6, nat, n_a)
   integer :: a, b, k
   real(wp) :: result, r_ab(3), qp_part(6, nat, n_a), sum_qm(6)

   delta_dipm = 0.0_wp
   delta_qp = 0.0_wp
   qp_part = 0.0_wp

   do k = 1, n_a
      do a = 1, nat
         do b = 1, nat

            result = inv_cn_a(a, b, k)

            r_ab = xyz(:, a) - xyz(:, b)

            delta_dipm(:, a, k) = delta_dipm(:, a, k) + (dipm(:, b) - r_ab(:)*q(b))/(result*(cn(b) + 1))
            !sorting of qp xx,xy,yy,xz,yz,zz

            delta_qp(1, a, k) = delta_qp(1, a, k) + &
               (1.5_wp*(-1*(r_ab(1)*dipm(1, b) + r_ab(1)*dipm(1, b)) + r_ab(1)*r_ab(1)*q(b)))/(result*(cn(b) + 1))

            delta_qp(2, a, k) = delta_qp(2, a, k) + &
               (1.5_wp*(-1*(r_ab(1)*dipm(2, b) + r_ab(2)*dipm(1, b)) + r_ab(1)*r_ab(2)*q(b)))/(result*(cn(b) + 1))

            delta_qp(3, a, k) = delta_qp(3, a, k) + &
               (1.5_wp*(-1*(r_ab(2)*dipm(2, b) + r_ab(2)*dipm(2, b)) + r_ab(2)*r_ab(2)*q(b)))/(result*(cn(b) + 1))

            delta_qp(4, a, k) = delta_qp(4, a, k) + &
               (1.5_wp*(-1*(r_ab(3)*dipm(1, b) + r_ab(1)*dipm(3, b)) + r_ab(1)*r_ab(3)*q(b)))/(result*(cn(b) + 1))

            delta_qp(5, a, k) = delta_qp(5, a, k) + &
               (1.5_wp*(-1*(r_ab(3)*dipm(2, b) + r_ab(2)*dipm(3, b)) + r_ab(2)*r_ab(3)*q(b)))/(result*(cn(b) + 1))

            delta_qp(6, a, k) = delta_qp(6, a, k) + &
               (1.5_wp*(-1*(r_ab(3)*dipm(3, b) + r_ab(3)*dipm(3, b)) + r_ab(3)*r_ab(3)*q(b)))/(result*(cn(b) + 1))

            qp_part(:, a, k) = qp_part(:, a, k) + &
               qp(:, b)/(result*(cn(b) + 1))
         end do
      end do
   end do

   call remove_trac_qp(nat, n_a, delta_qp, qp_part)

end subroutine

subroutine comp_norm_3(ndim, n_a, dipm, qm, dipm_norm, qm_norm)
   implicit none
   integer, INTENT(IN) :: ndim, n_a
   real(wp), INTENT(IN) :: dipm(3, ndim, n_a), qm(6, ndim, n_a)
   real(wp), INTENT(OUT) :: dipm_norm(ndim, n_a), qm_norm(ndim, n_a)
   real(wp) :: r(ndim, n_a), r2(ndim, n_a)
   INTEGER :: i

   do i = 1, ndim
      r2(i, :) = dipm(1, i, :)**2 + dipm(2, i, :)**2 + dipm(3, i, :)**2
   end do
   r = sqrt(r2)

   dipm_norm(:, :) = r(:, :)

   do i = 1, ndim
      r2(i, :) = qm(1, i, :)**2 + 2*qm(2, i, :)**2 + qm(3, i, :)**2 + 2*qm(4, i, :)**2 + 2*qm(5, i, :)**2 + qm(6, i, :)**2
   end do
   r = sqrt(r2)
   qm_norm(:, :) = r(:, :)

end subroutine

subroutine comp_norm(ndim, dipm, qm, dipm_norm, qm_norm)
   implicit none
   integer, INTENT(IN) :: ndim
   real(wp), INTENT(IN) :: dipm(3, ndim), qm(6, ndim)
   real(wp), INTENT(OUT) :: dipm_norm(ndim), qm_norm(ndim)
   real(wp) :: r(ndim), r2(ndim)
   INTEGER :: i

   do i = 1, ndim
      r2(i) = dipm(1, i)**2 + dipm(2, i)**2 + dipm(3, i)**2
   end do
   r = sqrt(r2)

   dipm_norm(:) = r(:)

   do i = 1, ndim
      r2(i) = qm(1, i)**2 + 2*qm(2, i)**2 + qm(3, i)**2 + 2*qm(4, i)**2 + 2*qm(5, i)**2 + qm(6, i)**2
   end do
   r = sqrt(r2)
   qm_norm(:) = r(:)

end subroutine

subroutine get_delta_mm_Z(nat, n_a, q, dipm, qp, at, xyz, cn, delta_dipm, delta_qp)! all effects due to the electrons are set to 0, only the distribution of positive charges is left
   implicit none
   integer, INTENT(IN) :: nat, at(nat), n_a
   real(wp), INTENT(IN) :: dipm(3, nat), xyz(3, nat), qp(6, nat), q(nat), cn(nat)
   real(wp), INTENT(OUT) :: delta_dipm(3, nat, n_a), delta_qp(6, nat, n_a)
   integer :: a, b, k
   real(wp) :: result, r_ab(3), qp_part(6, nat, n_a)

   delta_dipm = 0.0_wp
   delta_qp = 0.0_wp
   qp_part = 0.0_wp

   do k = 1, n_a
      do a = 1, nat
         do b = 1, nat
            result = inv_cn_a(a, b, k)
            r_ab = xyz(:, a) - xyz(:, b)
            delta_dipm(:, a, k) = delta_dipm(:, a, k) + (-r_ab(:)*q(b))/(result*(cn(b) + 1))
            !sorting of qp xx,xy,yy,xz,yz,zz

            delta_qp(1, a, k) = delta_qp(1, a, k) + &
               (1.5_wp*(r_ab(1)*r_ab(1)*q(b)))/(result*(cn(b) + 1))

            delta_qp(2, a, k) = delta_qp(2, a, k) + &
               (1.5_wp*(r_ab(1)*r_ab(2)*q(b)))/(result*(cn(b) + 1))

            delta_qp(3, a, k) = delta_qp(3, a, k) + &
               (1.5_wp*(r_ab(2)*r_ab(2)*q(b)))/(result*(cn(b) + 1))

            delta_qp(4, a, k) = delta_qp(4, a, k) + &
               (1.5_wp*(r_ab(1)*r_ab(3)*q(b)))/(result*(cn(b) + 1))

            delta_qp(5, a, k) = delta_qp(5, a, k) + &
               (1.5_wp*(r_ab(2)*r_ab(3)*q(b)))/(result*(cn(b) + 1))

            delta_qp(6, a, k) = delta_qp(6, a, k) + &
               (1.5_wp*(r_ab(3)*r_ab(3)*q(b)))/(result*(cn(b) + 1))
            !qp_part(:,a) = qp_part(:,a) + qp(:,b) / (result*(cn(b)+1))
         end do
         !delta_dipm(:,a) = dipm(:,a) + delta_dipm(:,a)
         !delta_qp(:,a) = qp(:,a) + delta_qp(:,a)
      end do
   end do

   call remove_trac_qp(nat, n_a, delta_qp, qp_part)

end subroutine

subroutine get_delta_mm_p(nat, n_a, q, dipm, qp, at, xyz, cn, delta_dipm, delta_qp) ! the sign of q was changed to respect the charge of the electrons
   implicit none
   integer, INTENT(IN) :: nat, at(nat), n_a
   real(wp), INTENT(IN) :: dipm(3, nat), xyz(3, nat), qp(6, nat), q(nat), cn(nat)
   real(wp), INTENT(OUT) :: delta_dipm(3, nat, n_a), delta_qp(6, nat, n_a)
   integer :: a, b, k
   real(wp) :: result, r_ab(3), qp_part(6, nat, n_a)

   delta_dipm = 0.0_wp
   delta_qp = 0.0_wp
   qp_part = 0.0_wp

   do k = 1, n_a
      do a = 1, nat
         do b = 1, nat
            result = inv_cn_a(a, b, k)
            r_ab = xyz(:, a) - xyz(:, b)

            delta_dipm(:, a, k) = delta_dipm(:, a, k) + &
               (dipm(:, b) + r_ab(:)*q(b))/(result*(cn(b) + 1))
            !sorting of qp xx,xy,yy,xz,yz,zz

            delta_qp(1, a, k) = delta_qp(1, a, k) + &
               (1.5_wp*(-1*(r_ab(1)*dipm(1, b) + r_ab(1)*dipm(1, b)) - r_ab(1)*r_ab(1)*q(b)))/(result*(cn(b) + 1))

            delta_qp(2, a, k) = delta_qp(2, a, k) + &
               (1.5_wp*(-1*(r_ab(1)*dipm(2, b) + r_ab(2)*dipm(1, b)) - r_ab(1)*r_ab(2)*q(b)))/(result*(cn(b) + 1))

            delta_qp(3, a, k) = delta_qp(3, a, k) + &
               (1.5_wp*(-1*(r_ab(2)*dipm(2, b) + r_ab(2)*dipm(2, b)) - r_ab(2)*r_ab(2)*q(b)))/(result*(cn(b) + 1))

            delta_qp(4, a, k) = delta_qp(4, a, k) + &
               (1.5_wp*(-1*(r_ab(3)*dipm(1, b) + r_ab(1)*dipm(3, b)) - r_ab(1)*r_ab(3)*q(b)))/(result*(cn(b) + 1))

            delta_qp(5, a, k) = delta_qp(5, a, k) + &
               (1.5_wp*(-1*(r_ab(3)*dipm(2, b) + r_ab(2)*dipm(3, b)) - r_ab(2)*r_ab(3)*q(b)))/(result*(cn(b) + 1))

            delta_qp(6, a, k) = delta_qp(6, a, k) + &
               (1.5_wp*(-1*(r_ab(3)*dipm(3, b) + r_ab(3)*dipm(3, b)) - r_ab(3)*r_ab(3)*q(b)))/(result*(cn(b) + 1))
            qp_part(:, a, k) = qp_part(:, a, k) + qp(:, b)/(result*(cn(b) + 1))
         end do
         !delta_dipm(:,a) = dipm(:,a) + delta_dipm(:,a)
         !delta_qp(:,a) = qp(:,a) + delta_qp(:,a)
      end do
   end do

   call remove_trac_qp(nat, n_a, delta_qp, qp_part)

end subroutine

subroutine remove_trac_qp(nat, n_a, qp_matrix, qp_part)
   implicit none
   integer, INTENT(IN) :: nat, n_a
   real(wp), INTENT(IN) :: qp_part(6, nat, n_a)
   real(wp), INTENT(INOUT) :: qp_matrix(6, nat, n_a)
   integer :: i
   real(wp) :: tii(n_a)

   do i = 1, nat
      tii = qp_matrix(1, i, :) + qp_matrix(3, i, :) + qp_matrix(6, i, :)
      tii = tii/3.0_wp
      !qp_matrix(1:6,i) = 1.50_wp*qp(1:6,i)
      qp_matrix(1, i, :) = qp_matrix(1, i, :) - tii
      qp_matrix(3, i, :) = qp_matrix(3, i, :) - tii
      qp_matrix(6, i, :) = qp_matrix(6, i, :) - tii
      qp_matrix(:, i, :) = qp_matrix(:, i, :) + qp_part(:, i, :)
   end do
end subroutine

subroutine sum_up_mulliken(nat, nshell, aoat2, ash, mull_shell, mull_at)
   implicit none
   integer, INTENT(IN) :: nat, nshell, aoat2(:), ash(:)
   real(wp), INTENT(IN) :: mull_shell(nshell)
   real(wp), INTENT(OUT) :: mull_at(nat)
   integer :: i

   mull_at = 0.0_wp

   do i = 1, nshell
      mull_at(ash(i)) = mull_at(ash(i)) + mull_shell(i)
   end do

end subroutine

subroutine mol_set_nuclear_charge(nat, at, id, z)
   implicit none
   integer, intent(in)::nat, at(nat), id(nat)
   real(wp), intent(out) :: z(nat)
   integer :: i
   do i = 1, nat
      z(i) = real(at(id(i)), wp) - real(ncore(at(id(i))))
      if (at(i) > 57 .and. at(i) < 72) z(i) = 3.0_wp
   end do
contains

pure elemental integer function ncore(at)
   integer, intent(in) :: at
   if (at .le. 2) then
      ncore = 0
   elseif (at .le. 10) then
      ncore = 2
   elseif (at .le. 18) then
      ncore = 10
   elseif (at .le. 29) then   !zn
      ncore = 18
   elseif (at .le. 36) then
      ncore = 28
   elseif (at .le. 47) then
      ncore = 36
   elseif (at .le. 54) then
      ncore = 46
   elseif (at .le. 71) then
      ncore = 54
   elseif (at .le. 79) then
      ncore = 68
   elseif (at .le. 86) then
      ncore = 78
   end if
end function ncore
end subroutine mol_set_nuclear_charge

subroutine pack_mult_xyz(mult_xyz, res, start_id, nat)
   real(wp), intent(in) :: mult_xyz(:, :)
   type(results_type), intent(inout) :: res
   integer, intent(in) :: start_id
   integer :: i, j, k
   do k = 1, nat
      j = 1
      do i = start_id, start_id + size(mult_xyz, dim=1) - 1
         res%ml_features(k, i) = mult_xyz(j, k)
         j = j + 1
      end do
   end do
end subroutine

subroutine pack_mult_xyz_shell(mult_xyz, res, start_id, nat, at2nsh)
   real(wp), intent(in) :: mult_xyz(:, :)
   type(results_type), intent(inout) :: res
   integer, intent(in) :: start_id, at2nsh(:)
   integer :: j, k, nsh, id_tmp
   nsh = 1

   do k = 1, nat
      id_tmp = start_id
      j = 1
      do i = id_tmp, (id_tmp + size(mult_xyz, dim=1) - 1)
         res%ml_features(k, i) = mult_xyz(j, nsh)
         j = j + 1
      end do
      nsh = nsh + 1

      if (at2nsh(k) == 2) then
         id_tmp = id_tmp + size(mult_xyz, dim=1)
         j = 1
         do i = id_tmp, id_tmp + size(mult_xyz, dim=1) - 1
            res%ml_features(k, i) = mult_xyz(j, nsh)
            j = j + 1
         end do
         nsh = nsh + 1
      elseif (at2nsh(k) == 3) then
         id_tmp = id_tmp + size(mult_xyz, dim=1)
         j = 1
         do i = id_tmp, id_tmp + size(mult_xyz, dim=1) - 1
            res%ml_features(k, i) = mult_xyz(j, nsh)
            j = j + 1
         end do
         nsh = nsh + 1
         id_tmp = id_tmp + size(mult_xyz, dim=1)
         j = 1
         do i = id_tmp, id_tmp + size(mult_xyz, dim=1) - 1
            res%ml_features(k, i) = mult_xyz(j, nsh)
            j = j + 1
         end do
         nsh = nsh + 1
      end if

   end do
end subroutine

subroutine pack_shellwise(shell_prop, res, start_id, at2nsh, nat)
   type(results_type), intent(inout) :: res
   real(wp), intent(in) :: shell_prop(:)
   integer, intent(in) ::  nat, at2nsh(:)
   integer, intent(in) :: start_id
   integer :: nsh
   nsh = 1
   do i = 1, nat
      res%ml_features(i, start_id) = shell_prop(nsh) !s shell always filled
      nsh = nsh + 1
      if (at2nsh(i) == 2) then
         res%ml_features(i, start_id + 1) = shell_prop(nsh)
         nsh = nsh + 1
      elseif (at2nsh(i) == 3) then
         res%ml_features(i, start_id + 1) = shell_prop(nsh)
         nsh = nsh + 1
         res%ml_features(i, start_id + 2) = shell_prop(nsh)
         nsh = nsh + 1
      end if
   end do
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
