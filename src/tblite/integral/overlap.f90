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

!> @file tblite/integral/overlap.f90
!> Provides evaluation of overlap integrals

!> Implementation of overlap integrals
module tblite_integral_overlap
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use tblite_basis_type, only : basis_type, cgto_type
   use tblite_integral_trafo, only : transform0, transform1, transform2
   use tblite_blas, only: gemm
   implicit none
   private

   public :: overlap_cgto, overlap_grad_cgto
   public :: get_overlap
   public :: maxl, msao

   interface get_overlap
      module procedure :: get_overlap_lat
      module procedure :: get_overlap_diatframe_lat
   end interface get_overlap

   interface eff_equality
      module procedure :: eff_equality_two_numbers
   end interface eff_equality

   integer, parameter :: maxl = 6
   integer, parameter :: maxl2 = maxl*2
   integer, parameter :: msao(0:maxl) = [1, 3, 5, 7, 9, 11, 13]
   integer, parameter :: mlao(0:maxl) = [1, 3, 6, 10, 15, 21, 28]
   integer, parameter :: lmap(0:maxl) = [0, 1, 4, 10, 20, 35, 56]
   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: sqrtpi3 = sqrtpi**3

   ! x (+1), y (-1), z (0) in [-1, 0, 1] sorting
   integer, parameter :: lx(3, 84) = reshape([&
      & 0, &
      & 0,0,1, &
      & 2,0,0,1,1,0, &
      & 3,0,0,2,2,1,0,1,0,1, &
      & 4,0,0,3,3,1,0,1,0,2,2,0,2,1,1, &
      & 5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1, &
      & 6,0,0,3,3,0,5,5,1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2, &
      & 0, &
      & 1,0,0, &
      & 0,2,0,1,0,1, &
      & 0,3,0,1,0,2,2,0,1,1, &
      & 0,4,0,1,0,3,3,0,1,2,0,2,1,2,1, &
      & 0,5,0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2, &
      & 0,6,0,3,0,3,1,0,0,1,5,5,2,0,0,2,4,4,2,1,3,1,3,2,1,4,1,2, &
      & 0, &
      & 0,1,0, &
      & 0,0,2,0,1,1, &
      & 0,0,3,0,1,0,1,2,2,1, &
      & 0,0,4,0,1,0,1,3,3,0,2,2,1,1,2, &
      & 0,0,5,0,2,0,3,2,3,0,1,0,1,4,4,3,1,1,1,2,2, &
      & 0,0,6,0,3,3,0,1,5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3,1,1,4,2], &
      & shape(lx), order=[2, 1])


contains


elemental function overlap_1d(moment, alpha) result(overlap)
   integer, intent(in) :: moment
   real(wp), intent(in) :: alpha
   real(wp) :: overlap
   real(wp), parameter :: dfactorial(0:7) = & ! see OEIS A001147
      & [1._wp,1._wp,3._wp,15._wp,105._wp,945._wp,10395._wp,135135._wp]

   if (modulo(moment, 2) == 0) then
      overlap = (0.5_wp/alpha)**(moment/2) * dfactorial(moment/2)
   else
      overlap = 0.0_wp
   end if
end function overlap_1d


pure subroutine horizontal_shift(ae, l, cfs)
   integer, intent(in) :: l
   real(wp), intent(in) :: ae
   real(wp), intent(inout) :: cfs(*)
   select case(l)
   case(0) ! s
      continue
   case(1) ! p
      cfs(1)=cfs(1)+ae*cfs(2)
   case(2) ! d
      cfs(1)=cfs(1)+ae*ae*cfs(3)
      cfs(2)=cfs(2)+ 2*ae*cfs(3)
   case(3) ! f
      cfs(1)=cfs(1)+ae*ae*ae*cfs(4)
      cfs(2)=cfs(2)+ 3*ae*ae*cfs(4)
      cfs(3)=cfs(3)+ 3*ae*cfs(4)
   case(4) ! g
      cfs(1)=cfs(1)+ae*ae*ae*ae*cfs(5)
      cfs(2)=cfs(2)+ 4*ae*ae*ae*cfs(5)
      cfs(3)=cfs(3)+ 6*ae*ae*cfs(5)
      cfs(4)=cfs(4)+ 4*ae*cfs(5)
   end select
end subroutine horizontal_shift

pure subroutine form_product(a, b, la, lb, d)
   integer, intent(in) :: la, lb
   real(wp), intent(in) :: a(*), b(*)
   real(wp), intent(inout) :: d(*)
   if(la.ge.4.or.lb.ge.4) goto 40
   if(la.ge.3.or.lb.ge.3) goto 30
   if(la.ge.2.or.lb.ge.2) goto 20
   ! <s|s> = <s>
   d(1)=a(1)*b(1)
   if(la.eq.0.and.lb.eq.0) return
   ! <s|p> = <s|*(|s>+|p>)
   !       = <s> + <p>
   d(2)=a(1)*b(2)+a(2)*b(1)
   if(la.eq.0.or.lb.eq.0) return
   ! <p|p> = (<s|+<p|)*(|s>+|p>)
   !       = <s> + <p> + <d>
   d(3)=a(2)*b(2)
   return
20 continue
   ! <s|d> = <s|*(|s>+|p>+|d>)
   !       = <s> + <p> + <d>
   d(1)=a(1)*b(1)
   d(2)=a(1)*b(2)+a(2)*b(1)
   d(3)=a(1)*b(3)+a(3)*b(1)
   if(la.eq.0.or.lb.eq.0) return
   ! <p|d> = (<s|+<p|)*(|s>+|p>+|d>)
   !       = <s> + <p> + <d> + <f>
   d(3)=d(3)+a(2)*b(2)
   d(4)=a(2)*b(3)+a(3)*b(2)
   if(la.le.1.or.lb.le.1) return
   ! <d|d> = (<s|+<p|+<d|)*(|s>+|p>+|d>)
   !       = <s> + <p> + <d> + <f> + <g>
   d(5)=a(3)*b(3)
   return
30 continue
   ! <s|f> = <s|*(|s>+|p>+|d>+|f>)
   !       = <s> + <p> + <d> + <f>
   d(1)=a(1)*b(1)
   d(2)=a(1)*b(2)+a(2)*b(1)
   d(3)=a(1)*b(3)+a(3)*b(1)
   d(4)=a(1)*b(4)+a(4)*b(1)
   if(la.eq.0.or.lb.eq.0) return
   ! <p|f> = (<s|+<p|)*(|s>+|p>+|d>+|f>)
   !       = <s> + <p> + <d> + <f> + <g>
   d(3)=d(3)+a(2)*b(2)
   d(4)=d(4)+a(2)*b(3)+a(3)*b(2)
   d(5)=a(2)*b(4)+a(4)*b(2)
   if(la.le.1.or.lb.le.1) return
   ! <d|f> = (<s|+<p|+<d|)*(|s>+|p>+|d>+|f>)
   !       = <s> + <p> + <d> + <f> + <g> + <h>
   d(5)=d(5)+a(3)*b(3)
   d(6)=a(3)*b(4)+a(4)*b(3)
   if(la.le.2.or.lb.le.2) return
   ! <f|f> = (<s|+<p|+<d|+<f|)*(|s>+|p>+|d>+|f>)
   !       = <s> + <p> + <d> + <f> + <g> + <h> + <i>
   d(7)=a(4)*b(4)
   return
40 continue
   ! <s|g> = <s|*(|s>+|p>+|d>+|f>+|g>)
   !       = <s> + <p> + <d> + <f> + <g>
   d(1)=a(1)*b(1)
   d(2)=a(1)*b(2)+a(2)*b(1)
   d(3)=a(1)*b(3)+a(3)*b(1)
   d(4)=a(1)*b(4)+a(4)*b(1)
   d(5)=a(1)*b(5)+a(5)*b(1)
   if(la.eq.0.or.lb.eq.0) return
   ! <p|g> = (<s|+<p|)*(|s>+|p>+|d>+|f>+|g>)
   !       = <s> + <p> + <d> + <f> + <g> + <h>
   d(3)=d(3)+a(2)*b(2)
   d(4)=d(4)+a(2)*b(3)+a(3)*b(2)
   d(5)=d(5)+a(2)*b(4)+a(4)*b(2)
   d(6)=a(2)*b(5)+a(5)*b(2)
   if(la.le.1.or.lb.le.1) return
   ! <d|g> = (<s|+<p|+<d|)*(|s>+|p>+|d>+|f>+|g>)
   !       = <s> + <p> + <d> + <f> + <g> + <h> + <i>
   d(5)=d(5)+a(3)*b(3)
   d(6)=d(5)+a(3)*b(4)+a(4)*b(3)
   d(7)=a(3)*b(5)+a(5)*b(3)
   if(la.le.2.or.lb.le.2) return
   ! <f|g> = (<s|+<p|+<d|+<f|)*(|s>+|p>+|d>+|f>+|g>)
   !       = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k>
   d(7)=d(7)+a(4)*b(4)
   d(8)=a(4)*b(5)+a(5)*b(4)
   if(la.le.3.or.lb.le.3) return
   ! <g|g> = (<s|+<p|+<d|+<f|+<g|)*(|s>+|p>+|d>+|f>+|g>)
   !       = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k> + <l>
   d(9)=a(5)*b(5)

end subroutine form_product


pure subroutine overlap_3d(rpj, rpi, aj, ai, lj, li, s1d, s3d)
   real(wp), intent(in) :: rpi(3)
   real(wp), intent(in) :: rpj(3)
   real(wp), intent(in) :: ai
   real(wp), intent(in) :: aj
   integer, intent(in) :: li(3)
   integer, intent(in) :: lj(3)
   real(wp), intent(in) :: s1d(0:)
   real(wp), intent(out) :: s3d

   integer :: k, l
   real(wp) :: vi(0:maxl), vj(0:maxl), vv(0:maxl2), v1d(3)

   v1d(:) = 0.0_wp

   do k = 1, 3
      vv(:) = 0.0_wp
      vi(:) = 0.0_wp
      vj(:) = 0.0_wp
      vi(li(k)) = 1.0_wp
      vj(lj(k)) = 1.0_wp

      call horizontal_shift(rpi(k), li(k), vi)
      call horizontal_shift(rpj(k), lj(k), vj)
      call form_product(vi, vj, li(k), lj(k), vv)
      do l = 0, li(k) + lj(k)
         v1d(k) = v1d(k) + s1d(l) * vv(l)
      end do
   end do

   s3d = v1d(1) * v1d(2) * v1d(3)

end subroutine overlap_3d


pure subroutine overlap_grad_3d(rpj, rpi, aj, ai, lj, li, s1d, s3d, ds3d)
   real(wp), intent(in) :: rpi(3)
   real(wp), intent(in) :: rpj(3)
   real(wp), intent(in) :: ai
   real(wp), intent(in) :: aj
   integer, intent(in) :: li(3)
   integer, intent(in) :: lj(3)
   real(wp), intent(in) :: s1d(0:)
   real(wp), intent(out) :: s3d
   real(wp), intent(out) :: ds3d(3)

   integer :: k, l
   real(wp) :: vi(0:maxl), vj(0:maxl), vv(0:maxl2), v1d(3)
   real(wp) :: gi(0:maxl), gg(0:maxl2), g1d(3)

   v1d(:) = 0.0_wp
   g1d(:) = 0.0_wp

   do k = 1, 3
      vv(:) = 0.0_wp
      gg(:) = 0.0_wp
      vi(:) = 0.0_wp
      vj(:) = 0.0_wp
      gi(:) = 0.0_wp

      vi(li(k)) = 1.0_wp
      vj(lj(k)) = 1.0_wp
      gi(li(k)+1) = 2*ai
      if (li(k) > 0) gi(li(k)-1) = -li(k)

      call horizontal_shift(rpi(k), li(k)-1, gi)
      call horizontal_shift(rpi(k), li(k)+1, gi)
      call horizontal_shift(rpi(k), li(k), vi)
      call horizontal_shift(rpj(k), lj(k), vj)
      call form_product(vi, vj, li(k), lj(k), vv)
      call form_product(gi, vj, li(k)+1, lj(k), gg)
      do l = 0, li(k) + lj(k) + 1
         v1d(k) = v1d(k) + s1d(l) * vv(l)
         g1d(k) = g1d(k) + s1d(l) * gg(l)
      end do
   end do

   s3d = v1d(1) * v1d(2) * v1d(3)
   ds3d(1) = g1d(1) * v1d(2) * v1d(3)
   ds3d(2) = v1d(1) * g1d(2) * v1d(3)
   ds3d(3) = v1d(1) * v1d(2) * g1d(3)

end subroutine overlap_grad_3d


pure subroutine overlap_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap)
   !> Description of contracted Gaussian function on center i
   type(cgto_type), intent(in) :: cgtoi
   !> Description of contracted Gaussian function on center j
   type(cgto_type), intent(in) :: cgtoj
   !> Square distance between center i and j
   real(wp), intent(in) :: r2
   !> Distance vector between center i and j, ri - rj
   real(wp), intent(in) :: vec(3)
   !> Maximum value of integral prefactor to consider
   real(wp), intent(in) :: intcut
   !> Overlap integrals for the given pair i  and j
   real(wp), intent(out) :: overlap(msao(cgtoj%ang), msao(cgtoi%ang))

   integer :: ip, jp, mli, mlj, l
   real(wp) :: eab, oab, est, s1d(0:maxl2), rpi(3), rpj(3), cc, val, pre
   real(wp) :: s3d(mlao(cgtoj%ang), mlao(cgtoi%ang))

   s3d(:, :) = 0.0_wp

   do ip = 1, cgtoi%nprim
      do jp = 1, cgtoj%nprim
         eab = cgtoi%alpha(ip) + cgtoj%alpha(jp)
         oab = 1.0_wp/eab
         est = cgtoi%alpha(ip) * cgtoj%alpha(jp) * r2 * oab
         if (est > intcut) cycle
         pre = exp(-est) * sqrtpi3*sqrt(oab)**3
         rpi = -vec * cgtoj%alpha(jp) * oab
         rpj = +vec * cgtoi%alpha(ip) * oab
         do l = 0, cgtoi%ang + cgtoj%ang
            s1d(l) = overlap_1d(l, eab)
         end do
         cc = cgtoi%coeff(ip) * cgtoj%coeff(jp) * pre
         do mli = 1, mlao(cgtoi%ang)
            do mlj = 1, mlao(cgtoj%ang)
               call overlap_3d(rpj, rpi, cgtoj%alpha(jp), cgtoi%alpha(ip), &
                  & lx(:, mlj+lmap(cgtoj%ang)), lx(:, mli+lmap(cgtoi%ang)), &
                  & s1d, val)
               s3d(mlj, mli) = s3d(mlj, mli) + cc*val
            end do
         end do
      end do
   end do

   call transform0(cgtoj%ang, cgtoi%ang, s3d, overlap)

end subroutine overlap_cgto


pure subroutine overlap_grad_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap, doverlap)
   !> Description of contracted Gaussian function on center i
   type(cgto_type), intent(in) :: cgtoi
   !> Description of contracted Gaussian function on center j
   type(cgto_type), intent(in) :: cgtoj
   !> Square distance between center i and j
   real(wp), intent(in) :: r2
   !> Distance vector between center i and j, ri - rj
   real(wp), intent(in) :: vec(3)
   !> Maximum value of integral prefactor to consider
   real(wp), intent(in) :: intcut
   !> Overlap integrals for the given pair i  and j
   real(wp), intent(out) :: overlap(msao(cgtoj%ang), msao(cgtoi%ang))
   !> Overlap integral gradient for the given pair i  and j
   real(wp), intent(out) :: doverlap(3, msao(cgtoj%ang), msao(cgtoi%ang))

   integer :: ip, jp, mli, mlj, l
   real(wp) :: eab, oab, est, s1d(0:maxl2), rpi(3), rpj(3), cc, val, grad(3), pre
   real(wp) :: s3d(mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: ds3d(3, mlao(cgtoj%ang), mlao(cgtoi%ang))

   s3d(:, :) = 0.0_wp
   ds3d(:, :, :) = 0.0_wp

   do ip = 1, cgtoi%nprim
      do jp = 1, cgtoj%nprim
         eab = cgtoi%alpha(ip) + cgtoj%alpha(jp)
         oab = 1.0_wp/eab
         est = cgtoi%alpha(ip) * cgtoj%alpha(jp) * r2 * oab
         if (est > intcut) cycle
         pre = exp(-est) * sqrtpi3*sqrt(oab)**3
         rpi = -vec * cgtoj%alpha(jp) * oab
         rpj = +vec * cgtoi%alpha(ip) * oab
         do l = 0, cgtoi%ang + cgtoj%ang + 1
            s1d(l) = overlap_1d(l, eab)
         end do
         cc = cgtoi%coeff(ip) * cgtoj%coeff(jp) * pre
         do mli = 1, mlao(cgtoi%ang)
            do mlj = 1, mlao(cgtoj%ang)
               call overlap_grad_3d(rpj, rpi, cgtoj%alpha(jp), cgtoi%alpha(ip), &
                  & lx(:, mlj+lmap(cgtoj%ang)), lx(:, mli+lmap(cgtoi%ang)), &
                  & s1d, val, grad)
               s3d(mlj, mli) = s3d(mlj, mli) + cc*val
               ds3d(:, mlj, mli) = ds3d(:, mlj, mli) + cc*grad
            end do
         end do
      end do
   end do

   call transform0(cgtoj%ang, cgtoi%ang, s3d, overlap)
   call transform1(cgtoj%ang, cgtoi%ang, ds3d, doverlap)

end subroutine overlap_grad_cgto


!> Evaluate overlap for a molecular structure
subroutine get_overlap_lat(mol, trans, cutoff, bas, overlap)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Overlap matrix
   real(wp), intent(out) :: overlap(:, :)

   integer :: iat, jat, izp, jzp, itr, is, js
   integer :: ish, jsh, ii, jj, iao, jao, nao
   real(wp) :: r2, vec(3), cutoff2
   real(wp), allocatable :: stmp(:)

   overlap(:, :) = 0.0_wp

   allocate(stmp(msao(bas%maxl)**2))
   cutoff2 = cutoff**2

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, bas, trans, cutoff2, overlap) private(r2, vec, stmp) &
   !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj, iao, jao, nao)
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
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is+ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js+jsh)
                  call overlap_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp)

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  !$omp simd collapse(2)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                           & + stmp(jao + nao*(iao-1))
                     end do
                  end do

               end do
            end do

         end do
      end do
   end do

end subroutine get_overlap_lat

!> Evaluate overlap for a molecular structure,
!> with scaled elements in the diatomic frame
   subroutine get_overlap_diatframe_lat(mol, trans, cutoff, bas, scal_fac, overlap, overlap_scaled)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Lattice points within a given realspace cutoff
      real(wp), intent(in) :: trans(:, :)
      !> Realspace cutoff
      real(wp), intent(in) :: cutoff
      !> Basis set information
      type(basis_type), intent(in) :: bas
      !> Scaling factors for the diatomic frame
      real(wp), intent(in) :: scal_fac(:,:)
      !> Overlap matrix
      real(wp), intent(out) :: overlap(:, :)
      !> Overlap matrix with scaled elements in the diatomic frame
      real(wp), intent(out) :: overlap_scaled(:, :)
      !> Transformation vector for the diatomic frame
      real(wp) :: vec_diat_trafo(3)

      integer :: iat, jat, izp, jzp, itr, is, js
      integer :: nao_ati, nao_atj, lbi, lbj, ubi, ubj
      integer :: ish, jsh, ii, jj, iao, jao, nao, maxl_ish_jsh
      real(wp) :: r2, vec(3), cutoff2, ksig, kpi, kdel
      real(wp), allocatable :: stmp(:)
      real(wp) :: trafomat(9,9), block_overlap(9,9), trans_block_s(9,9),tmp(9,9)

      if (size(scal_fac,1) /= 3) then
         write(*,*) 'Error: scal_fac must have 3 rows'
         stop
      end if

      allocate(stmp(msao(bas%maxl)**2))
      cutoff2 = cutoff**2

      !$omp parallel do schedule(runtime) default(none) &
      !$omp shared(mol, bas, trans, cutoff2, overlap, overlap_scaled) private(r2, vec, stmp) &
      !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj, iao, jao, nao) &
      !$omp private(nao_ati, nao_atj, lbi, lbj, ubi, ubj, maxl_ish_jsh, trafomat, block_overlap) &
      !$omp private(trans_block_s, ksig, kpi, kdel, vec_diat_trafo)
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
               if (iat /= jat) then
                  call relvec(vec, sqrt(r2), vec_diat_trafo)
               end if
               maxl_ish_jsh = 0
               nao_ati = 0
               do ish = 1, bas%nsh_id(izp)
                  ii = bas%iao_sh(is+ish)
                  nao_ati = nao_ati + bas%nao_sh(is+ish)
                  nao_atj = 0
                  do jsh = 1, bas%nsh_id(jzp)
                     jj = bas%iao_sh(js+jsh)
                     nao_atj = nao_atj + bas%nao_sh(js+jsh)
                     call overlap_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp)

                     maxl_ish_jsh = max(maxl_ish_jsh, bas%cgto(ish, izp)%ang, bas%cgto(jsh, jzp)%ang)
                     nao = msao(bas%cgto(jsh, jzp)%ang)
                     !$omp simd collapse(2)
                     do iao = 1, msao(bas%cgto(ish, izp)%ang)
                        do jao = 1, nao
                           overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                           & + stmp(jao + nao*(iao-1))
                        end do
                     end do

                  end do
               end do

               !> Transform 9x9 submatrix (in minimal basis case with s,p,d) to diatomic frame,
               !> scale the elements with the corresponding factor,
               !> transform them back and add them to the scaled overlap matrix

               !> 1. Define elements of overlap that should be transformed (or copied back plain)
               lbi = bas%iao_sh(is+1)+1
               ubi = bas%iao_sh(is+1)+nao_ati
               lbj = bas%iao_sh(js+1)+1
               ubj = bas%iao_sh(js+1)+nao_atj
               ! ### DEV ###
               !write(*,*) "lbi, ubi, lbj, ubj: ", lbi, ubi, lbj, ubj
               ! ### DEV ###
               if (iat /= jat) then
                  block_overlap = 0.0_wp
                  !> Fill a 9x9 submatrix (initialized with 0's) with the overlap matrix elements

                  !> 2. Fill the submatrix with the correct overlap matrix elements
                  block_overlap(1:nao_atj, 1:nao_ati) = overlap(lbj:ubj, lbi:ubi)

                  !> (optional) 1.1. Transpose the submatrix to correspond to MSINDO processing,
                  ! ### Transposing back and forth doesn't make a difference for the total overlap matrix
                  ! ### but it does for the comparison of submatrices with MSINDO
                  ! block_overlap = transpose(block_overlap)

                  ! ### DEV ###
                  ! write(*,*) "vec_diat_trafo:"
                  ! write(*,*) vec_diat_trafo
                  ! ### DEV ###

                  call harmtr(2, vec_diat_trafo, trafomat)
                  ! Doing it only for the maxl_ish_jsh elements (for s,p -> 4x4) would in
                  ! theory be more efficient, but is not implemented yet:
                  ! call harmtr(maxl_ish_jsh, vec_diat_trafo, trafomat) 

                  ! ### DEV ###
                  ! write(*,*) "trafomat:"
                  ! do i = 1, 9
                  !    write(*,'(9f10.6)') (trafomat(i,j), j=1,9)
                  ! end do
                  ! write(*,*) "block_overlap:"
                  ! do i = 1, 9
                  !    write(*,'(9f10.6)') (block_overlap(i,j), j=1,9)
                  ! end do
                  ! ### DEV ###

                  !> 3. Transform the submatrix
                  ! trans_block_s = matmul(matmul(transpose(trafomat), block_overlap),trafomat)
                  ! trans_block_s = O^T * S * O
                  call gemm(amat=trafomat,bmat=block_overlap,cmat=tmp,transa='T',transb='N')
                  call gemm(amat=tmp,bmat=trafomat,cmat=trans_block_s,transa='N',transb='N')

                  ! ### DEV ###
                  ! write(*,*) "transformed block overlap:"
                  ! do i = 1, 9
                  !    write(*,'(9f10.6)') (trans_block_s(i,j), j=1,9)
                  ! end do
                  ! ### DEV ###

                  !> 4. Scale elements in diatomic frame
                  !> 4.1. Determine scaling factors from atom parameters
                  ksig = 2.0_wp / (1.0_wp / scal_fac(1,mol%num(mol%id(iat))) &
                  & + 1.0_wp / scal_fac(1,mol%num(mol%id(jat))) )
                  kpi = 2.0_wp / (1.0_wp / scal_fac(2,mol%num(mol%id(iat))) &
                  & + 1.0_wp / scal_fac(2,mol%num(mol%id(jat))) )
                  kdel = 2.0_wp / (1.0_wp / scal_fac(3,mol%num(mol%id(iat))) &
                  & + 1.0_wp / scal_fac(3,mol%num(mol%id(jat))) )
                  !> 4.2. Scale elements with equivalent bonding situation in the
                  !>      diatomic frame.
                  trans_block_s(1,1) = trans_block_s(1,1)*ksig ! Sigma bond s   <-> s
                  trans_block_s(1,3) = trans_block_s(1,3)*ksig ! Sigma bond s   <-> pz
                  trans_block_s(3,1) = trans_block_s(3,1)*ksig ! Sigma bond pz  <-> s
                  trans_block_s(1,5) = trans_block_s(1,5)*ksig ! Sigma bond s   <-> dz2
                  trans_block_s(5,1) = trans_block_s(5,1)*ksig ! Sigma bond dz2 <-> s
                  trans_block_s(3,3) = trans_block_s(3,3)*ksig ! Sigma bond pz  <-> pz
                  trans_block_s(3,5) = trans_block_s(3,5)*ksig ! Sigma bond pz  <-> dz2
                  trans_block_s(5,3) = trans_block_s(5,3)*ksig ! Sigma bond dz2 <-> pz
                  trans_block_s(5,5) = trans_block_s(5,5)*ksig ! Sigma bond dz2 <-> dz2
                  trans_block_s(4,4) = trans_block_s(4,4)*kpi  ! Pi    bond px  <-> px
                  trans_block_s(4,6) = trans_block_s(4,6)*kpi  ! Pi    bond px  <-> dxz
                  trans_block_s(6,4) = trans_block_s(6,4)*kpi  ! Pi    bond dxz <-> px
                  trans_block_s(6,6) = trans_block_s(6,6)*kpi  ! Pi    bond dxz <-> dxz
                  trans_block_s(2,2) = trans_block_s(2,2)*kpi  ! Pi    bond py  <-> py
                  trans_block_s(2,7) = trans_block_s(2,7)*kpi  ! Pi    bond py  <-> dyz
                  trans_block_s(7,2) = trans_block_s(7,2)*kpi  ! Pi    bond dyz <-> py
                  trans_block_s(7,7) = trans_block_s(7,7)*kpi  ! Pi    bond dyz <-> dyz
                  trans_block_s(8,8) = trans_block_s(8,8)*kdel ! Delta bond dx2-y2 <-> dx2-y2
                  trans_block_s(9,9) = trans_block_s(9,9)*kdel ! Delta bond dxy <-> dxy

                  ! ### DEV ###
                  ! write(*,*) "transformed scaled block overlap:"
                  ! do i = 1, 9
                  !    write(*,'(9f10.6)') (trans_block_s(i,j), j=1,9)
                  ! end do
                  ! ### DEV ###

                  !> 5. Transform back to original frame
                  ! block_overlap = matmul(matmul(trafomat, trans_block_s),transpose(trafomat))
                  ! block_overlap = O * S * O^T 
                  call gemm(amat=trafomat,bmat=trans_block_s,cmat=tmp,transa='N',transb='N')
                  call gemm(amat=tmp,bmat=trafomat,cmat=block_overlap,transa='N',transb='T')

                  ! ### DEV ###
                  ! write(*,*) "Scaled block overlap in cartesian frame:"
                  ! do i = 1, 9
                  !    write(*,'(9f10.6)') (trans_block_s(i,j), j=1,9)
                  ! end do
                  ! ### DEV ###

                  !> 6. Fill the overlap_diat matrix with the back-transformed submatrix
                  ! (see 'optional (1.1)' above)
                  ! trans_block_s = transpose(trans_block_s)
                  overlap_scaled(lbj:ubj, lbi:ubi) = block_overlap(1:nao_atj, 1:nao_ati)
               else
                  !> 7. Fill the overlap_diat matrix with the plain overlap submatrix (no transformation)
                  overlap_scaled(lbj:ubj, lbi:ubi) = overlap(lbj:ubj, lbi:ubi)
               endif
            end do
         end do
      end do

   end subroutine get_overlap_diatframe_lat

   subroutine relvec(vec, rkl, veckl)

      !> Original vector between atoms A and B
      real(wp), intent(in)             :: vec(3)
      !> Distance between the two atoms
      real(wp), intent(in)             :: rkl
      !> Normalized vector from atom k to atom l
      real(wp), intent(out)            :: veckl(3)

      real(wp), parameter              :: eps = 4.0e-08_wp

      real(wp)                         :: sq

      veckl(1:3) = vec(1:3) / rkl
      if ( abs(1.0_wp-abs(veckl(1))) .lt. eps ) then
         veckl(1) = sign(1.0_wp,veckl(1))
         veckl(2) = 0.0_wp
         veckl(3) = 0.0_wp
      else if ( abs(1.0_wp-abs(veckl(2))) .lt. eps ) then
         veckl(1) = 0.0_wp
         veckl(2) = sign(1.0_wp,veckl(2))
         veckl(3) = 0.0_wp
      else if ( abs(1.0_wp-abs(veckl(3))) .lt. eps ) then
         veckl(1) = 0.0_wp
         veckl(2) = 0.0_wp
         veckl(3) = sign(1.0_wp,veckl(3))
      else if ( (abs(veckl(1)) .lt. eps) .and. .not. eff_equality(veckl(1),0.0_wp) ) then
         veckl(1) = 0.0_wp
         sq = sqrt( veckl(2)**2 + veckl(3)**2 )
         veckl(2) = veckl(2)/sq
         veckl(3) = veckl(3)/sq
      else if ( (abs(veckl(2)) .lt. eps) .and. .not. eff_equality(veckl(2),0.0_wp) ) then
         veckl(2) = 0.0_wp
         sq = sqrt( veckl(1)**2 + veckl(3)**2 )
         veckl(1) = veckl(1)/sq
         veckl(3) = veckl(3)/sq
      else if ( (abs(veckl(3)) .lt. eps) .and. .not. eff_equality(veckl(3),0.0_wp) ) then
         veckl(3) = 0.0_wp
         sq = sqrt(veckl(1)**2 + veckl(2)**2)
         veckl(1) = veckl(1)/sq
         veckl(2) = veckl(2)/sq
      endif

   end subroutine relvec

   logical pure function eff_equality_two_numbers(num1, num2)
      !> Numbers to compare
      real(wp), intent(in) :: num1, num2
      !> Logical deciding if numbers are (almost) equal or not
      eff_equality_two_numbers = (abs( num1 - num2 ) .le. 1.0e-12_wp)

   end function eff_equality_two_numbers

   subroutine harmtr(maxkl,veckl,trafomat)
      !> Maximum angular momentum
      integer, intent(in)  :: maxkl
      !> Normalized vector from atom k to atom l
      real(wp), intent(in) :: veckl(3)
      !> Transformation matrix
      real(wp), intent(out) :: trafomat(9,9)

      real(wp) :: cos2p, cos2t, cosp, cost, sin2p, sin2t, sinp, sint, sqrt3

      !     ------------------------------------------------------------------
      if (maxkl > 2) then
         write(*,*) "ERROR: f function or higher ang. mom. not implemented in harmtr"
         stop
      endif

      trafomat = 0.0_wp
      !     -----------------------------
      !     *** s functions (trafomat(1x1)) ***
      !     -----------------------------
      trafomat(1,1) = 1.0

      if ( maxkl == 0 ) return
      !     -----------------------------
      !     *** p functions (trafomat(4x4)) ***
      !     -----------------------------

      cost = veckl(3)
      if ( abs(cost) .eq. 1.0_wp ) then
         sint = 0.0_wp
         cosp = 1.0_wp
         sinp = 0.0_wp
      else if ( abs(cost) .eq. 0.0_wp ) then
         sint = 1.0_wp
         cosp = veckl(1)
         sinp = veckl(2)
      else
         sint = SQRT(1.0_wp-COST**2)
         cosp = veckl(1)/SINT
         sinp = veckl(2)/SINT
      endif

      !> Original MSINDO ordering (-> x, y, z)
      ! trafomat(2,2) = SINT*COSP
      ! trafomat(3,2) = SINT*SINP
      ! trafomat(4,2) = COST
      ! trafomat(2,3) = COST*COSP
      ! trafomat(3,3) = COST*SINP
      ! trafomat(4,3) = -SINT
      ! trafomat(2,4) = -SINP
      ! trafomat(3,4) = COSP
      ! trafomat(4,4) = 0.0_wp

      !> tblite ordering with adapted column ordering
      ! 1st index:
      ! MSINDO defintion of p function ordering is converted to
      ! tblite definition of p function ordering. E.g. for first entry:
      ! trafomat(2,:)_MSINDO -> trafomat(px,:) -> trafomat(4:)_tblite
      ! 2nd index:
      ! Final ordering of p functions (see below) corresponds to the
      ! tblite ordering of p functions. For the second index, the ordering
      ! 3, 4, 2 holds, corresponding (in MSINDO convention)
      ! to the tblite ordering of p functions, i.e.
      ! y, z, x
      trafomat(4,3) = SINT*COSP
      trafomat(2,3) = SINT*SINP
      trafomat(3,3) = COST
      trafomat(4,4) = COST*COSP
      trafomat(2,4) = COST*SINP
      trafomat(3,4) = -SINT
      trafomat(4,2) = -SINP
      trafomat(2,2) = COSP
      trafomat(3,2) = 0.0_wp

      if ( maxkl <= 1 ) return

!     -----------------------------
!     *** d functions (trafomat(9x9)) ***
!     -----------------------------

      COS2T = COST**2 - SINT**2
      SIN2T = 2.0_wp * SINT*COST
      COS2P = COSP**2 - SINP**2
      SIN2P = 2.0_wp * SINP*COSP
      SQRT3 = SQRT(3.0_wp)

      !> Original MSINDO ordering
      !> The MSINDO d SAO ordering apparently corresponds to the
      !> tblite ordering of d SAOs, which is why it doesn't have
      !> to be adapted
      trafomat(5,5) = (3.0_wp * COST**2 - 1.0_wp) * 0.5_wp
      trafomat(6,5) = SQRT3*SIN2T*COSP*0.5_wp
      trafomat(7,5) = SQRT3*SIN2T*SINP*0.5_wp
      trafomat(8,5) = SQRT3*SINT**2*COS2P*0.5_wp
      trafomat(9,5) = SQRT3*SINT**2*SIN2P*0.5_wp
      trafomat(5,6) = -SQRT3*SIN2T*0.5_wp
      trafomat(6,6) = COS2T*COSP
      trafomat(7,6) = COS2T*SINP
      trafomat(8,6) = SIN2T*COS2P*0.5_wp
      trafomat(9,6) = SIN2T*SIN2P*0.5_wp
      trafomat(5,7) = 0.0_wp
      trafomat(6,7) = -COST*SINP
      trafomat(7,7) = COST*COSP
      trafomat(8,7) = -SINT*SIN2P
      trafomat(9,7) = SINT*COS2P
      trafomat(5,8) = SQRT3*SINT**2 * 0.5_wp
      trafomat(6,8) = -SIN2T*COSP*0.5_wp
      trafomat(7,8) = -SIN2T*SINP*0.5_wp
      trafomat(8,8) = (1.0_wp + COST**2) * COS2P * 0.5_wp
      trafomat(9,8) = (1.0_wp + COST**2) * SIN2P * 0.5_wp
      trafomat(5,9) = 0.0_wp
      trafomat(6,9) = SINT*SINP
      trafomat(7,9) = -SINT*COSP
      trafomat(8,9) = -COST*SIN2P
      trafomat(9,9) = COST*COS2P

   end subroutine harmtr

end module tblite_integral_overlap
