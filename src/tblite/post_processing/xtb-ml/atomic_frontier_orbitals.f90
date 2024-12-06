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

!> @file tblite/xtb-ml/atomic_frontier_orbitals.f90

!> compute atomic response function and effective H-L gap
! nat:  # of atoms
! nao:  # of AOs/MOs (assuming spherical basis here)
! focca,b: (fractional) occupation numbers for (alpha, beta) part
! eps:  orbital energies
! aoat: ao-to-orbital indexer
! C: MO coefficient
! S: overlap matrix
! response: atomic response (measures level availability on atom)
! egap: effective atomic H-L gap
! chempot: efffective Fermi level per atom (as 0.5*(eh+el)
! ehoaoa,b: highest occ. atomic orbital for (alpha, beta) part
! eluaoa,b: lowest virt. atomic orbital for (alpha, beta) part
module tblite_xtbml_atomic_frontier
   use mctc_io_convert, only : autoev
   use tblite_timer, only : timer_type, format_time
   use tblite_context, only : context_type
   use mctc_env, only : wp, sp
   use tblite_output_format, only : format_string
   implicit none
   private
   public :: atomic_frontier_orbitals

contains
   subroutine atomic_frontier_orbitals(nat, nao, focca, foccb, eps, aoat, C, S, response, egap, &
      chempot, ehoaoa, eluaoa, ehoaob, eluaob, print, ctx)
   integer, intent(in)    :: nat, nao, aoat(nao)
   real(wp), intent(in)    :: focca(nao), foccb(nao), eps(nao)
   real(wp), intent(in)    :: C(nao, nao), S(nao, nao)
   type(context_type), intent(inout) :: ctx
   real(wp), intent(out)  :: ehoaoa(nat), eluaoa(nat), ehoaob(nat), eluaob(nat), egap(nat), response(nat), chempot(nat)
   integer  :: i, j, k, jj, kk, m
   real(wp):: po(nat, nao), pv(nat, nao)
   real(wp) :: occa, occb, tmp, tmp2, ps, virta, virtb, tmp3, weight, osite, vsite
   real(wp), parameter :: damp = 0.5_wp ! damping in response function (in eV)
   logical, intent(in) :: print
   type(timer_type) :: timer
   logical, parameter :: debug = .false.
   character(len=150)  :: tmp_str
   call timer%push("total")
   call timer%push("occupation a")

   po = 0.0_wp
   pv = 0.0_wp

   ! collect the atomic MO populations (alpha)
   ! we collect occ & virt separately (and allow for fractional occupation)
   !$omp parallel do default(private) schedule(runtime) reduction(+:po)&
   !$omp shared(aoat, nao, s, C, focca) &
   !$omp private(ps, occa, jj, kk, j, i, k)
   do i = 1, nao
      ! occ part
      occa = focca(i)
      if (occa .gt. 0.0001_wp) then
         do j = 1, nao
            jj = aoat(j)
            do k = 1, j - 1
               kk = aoat(k)
               ps = s(k, j)*C(j, i)*C(k, i)*occa
               po(kk, i) = po(kk, i) + ps
               po(jj, i) = po(jj, i) + ps
            end do
            ps = s(j, j)*C(j, i)*C(j, i)*occa
            po(jj, i) = po(jj, i) + ps
         end do
      end if
   end do
   
   !$omp parallel do default(private) schedule(runtime) reduction(+:pv)&
   !$omp shared(aoat, nao, s, C, focca) &
   !$omp private(ps, jj, kk, virta, j, i, k)
   do i = 1, nao
      virta = (1.0_wp - focca(i))
      if (virta .gt. 0.0001_wp) then
         do j = 1, nao
            jj = aoat(j)
            do k = 1, j - 1
               kk = aoat(k)
               ps = s(k, j)*C(j, i)*C(k, i)*virta
               pv(kk, i) = pv(kk, i) + ps
               pv(jj, i) = pv(jj, i) + ps
            end do
            ps = s(j, j)*C(j, i)*C(j, i)*virta
            pv(jj, i) = pv(jj, i) + ps
         end do
      end if
   end do

   call timer%pop()

   ! now accumulate the atomic H-L gaps and Fermi levels (we approximate it as 0.5*(eH + eL))
   ! we make use of an atomic response-type weightinhg:
   ! chempot=\sum_ai * wA_ai * 0.5 * (e_a + e_i)
   ! here wA_ai is the variable "response" computed as: wA_ai=  [p_i p_a / ( (e_a - e_i)**2 + damp**2) ]/ [\sum_ia p_j p_b / ((e_b - e_j)**2 + damp**2)]
   ! the "p"s are the MO densities, for gaps we accumulate the regularized inverse and later on invert it again
   call timer%push("prep a")

   response = 0.0_wp
   chempot = 0.0_wp
   egap = 0.0_wp
   ! go through occ (including fractionally occupied)
   !!$omp parallel do default(private) &
   !!$omp shared( eps,nao,focca,po,pv,response,chempot_tmp,egap_tmp,nat) &
   !!$omp private(occa,virta,tmp,tmp2,tmp3)&
   !!$omp private(weight,i,j,m)
   do i = 1, nao
      occa = focca(i)
      if (occa .gt. 0.0001_wp) then
         osite = sum(po(:, i))/maxval(po(:, i)) ! not used, but roughly gives the number of centers for the canonical MOs
         !virt part  (including fractionally occupied)
         do j = 1, nao
            virta = (1.0_wp - focca(j))
            if (virta .gt. 0.0001_wp) then
               vsite = sum(pv(:, j))/maxval(pv(:, j)) ! not used, but roughly gives the number of centers for the canonical MOs
               tmp = 1.0_wp/((eps(j) - eps(i))**2 + damp**2)
               tmp2 = 0.5_wp*(eps(j) + eps(i))
               tmp3 = 1.0_wp/(eps(j) - eps(i) + damp)
               ! compute atomic response
               do m = 1, nat
                  weight = po(m, i)*pv(m, j)*tmp
                  response(m) = response(m) + weight
                  chempot(m) = chempot(m) + tmp2*weight
                  egap(m) = egap(m) + weight*tmp3
               end do ! at
            end if ! if virt
         end do !  virt
      end if ! if occ
   end do ! occ
   !!$omp end parallel do
   call timer%pop()
   if (print) then
      call ctx%message("")
      call ctx%message( "  -------------------------")
      call ctx%message( "  atomic frontier MO (alpha) info ")
      call ctx%message( "  -------------------------")
      call ctx%message( "  atom   response (a.u.)   gap (eV)  chem.pot (eV)  HOAO (eV)    LUAO (eV)")
   end if
   call timer%push("AO A")

   ehoaoa = 0.0_wp
   eluaoa = 0.0_wp

   ! now get the atomic frontier orbitals
   tmp = 0.0_wp
   do m = 1, nat
      egap(m) = egap(m)/(response(m))
      egap(m) = 1.0_wp/egap(m) - damp
      chempot(m) = chempot(m)/(response(m))
      ehoaoa(m) = (chempot(m) - 0.5_wp*egap(m))
      eluaoa(m) = (chempot(m) + 0.5_wp*egap(m))

      ! fix cases for missing HOAO/LUAO
      if (sum(po(m, :)) .lt. 1.0d-12) then ! there is no occupied orbital for this atom
         !virt part  (including fractionally occupied)
         egap(m) = 0.0_wp
         chempot(m) = 0.0_wp
         do j = 1, nao
            virta = (1.0_wp - focca(j))
            if (virta .gt. 0.0001_wp) then
               tmp = 1.0_wp/((eps(j))**2 + damp**2)
               tmp2 = eps(j)
               tmp3 = 1.0_wp/(eps(j) + 1.0d100 + damp)
               ! compute atomic response
               weight = pv(m, j)*tmp
               chempot(m) = chempot(m) + tmp2*weight
               egap(m) = egap(m) + weight*tmp3
            end if ! if virt
            egap(m) = 1.0_wp/(egap(m) + 1.0d-12) - damp
            eluaoa(m) = chempot(m)
            ehoaoa(m) = chempot(m) - egap(m)
         end do !  virt
      end if
      if (sum(pv(m, :)) .lt. 1.0d-12) then
         !occ part  (including fractionally occupied)
         egap(m) = 0.0_wp
         chempot(m) = 0.0_wp
         do j = 1, nao
            occa = focca(j)
            if (occa .gt. 0.0001_wp) then
               tmp = 1.0_wp/((eps(j))**2 + damp**2)
               tmp2 = eps(j)
               tmp3 = 1.0_wp/(1.0d100 - eps(j) + damp)
               ! compute atomic response
               weight = po(m, j)*tmp
               chempot(m) = chempot(m) + tmp2*weight
               egap(m) = egap(m) + weight*tmp3
            end if ! if occ
            egap(m) = 1.0_wp/(egap(m) + 1.0d-12) - damp
            ehoaoa(m) = chempot(m)
            eluaoa(m) = chempot(m) + egap(m)
         end do !  occ
      end if
   end do

   if (print) then
      do m = 1, nat
         write (tmp_str, '(3x,i6,5x,f8.4,5x,f8.4,5x,f8.4,5x,f8.4,5x,f8.4)') m, &
            response(m)*(autoev**2), egap(m), chempot(m), ehoaoa(m), eluaoa(m)
         call ctx%message(trim(tmp_str))
      end do
      call ctx%message("")
   end if
   call timer%pop()
   call timer%push("occupation b")
   po = 0.0_wp
   pv = 0.0_wp
   
   !$omp parallel do default(private) schedule(runtime) reduction(+:po)&
   !$omp shared(aoat, nao, s, C, foccb) &
   !$omp private(ps, occb, jj, kk, j, i, k)
   do i = 1, nao
      ! occ part
      occb = foccb(i)
      if (occb .gt. 0.0001_wp) then
         do j = 1, nao
            jj = aoat(j)
            do k = 1, j - 1
               kk = aoat(k)
               ps = s(k, j)*C(j, i)*C(k, i)*occb
               po(kk, i) = po(kk, i) + ps
               po(jj, i) = po(jj, i) + ps
            end do
            ps = s(j, j)*C(j, i)*C(j, i)*occb
            po(jj, i) = po(jj, i) + ps
         end do
      end if
   end do
   
   !$omp parallel do default(private) schedule(runtime) reduction(+:pv)&
   !$omp shared(aoat, nao, s, C, foccb) &
   !$omp private(ps, jj, kk, virtb, j, i, k)
   do i = 1, nao
      virtb = (1.0_wp - foccb(i))
      if (virtb .gt. 0.0001_wp) then
         do j = 1, nao
            jj = aoat(j)
            do k = 1, j - 1
               kk = aoat(k)
               ps = s(k, j)*C(j, i)*C(k, i)*virtb
               pv(kk, i) = pv(kk, i) + ps
               pv(jj, i) = pv(jj, i) + ps
            end do
            ps = s(j, j)*C(j, i)*C(j, i)*virtb
            pv(jj, i) = pv(jj, i) + ps
         end do
      end if
   end do
   
   call timer%pop()
   call timer%push("prep b")

   response = 0.0_wp
   chempot = 0.0_wp
   egap = 0.0_wp
   !!$omp parallel do default(private) &
   !!$omp shared( eps,nao,foccb,po,pv) &
   !!$omp private(occb,virtb,tmp,tmp2,tmp3,weight,response,chempot_tmp,egap_tmp)
   ! go through occ (including fractionally occupied)
   do i = 1, nao
      occb = foccb(i)
      if (occb .gt. 0.0001_wp) then
         osite = sum(po(:, i))/maxval(po(:, i)) ! not used, but roughly gives the number of centers for the canonical MOs
         !virt part  (including fractionally occupied)
         do j = 1, nao
            virtb = (1.0_wp - foccb(j))
            if (virtb .gt. 0.0001_wp) then
               vsite = sum(pv(:, j))/maxval(pv(:, j)) ! not used, but roughly gives the number of centers for the canonical MOs
               tmp = 1.0_wp/((eps(j) - eps(i))**2 + damp**2)
               tmp2 = 0.5_wp*(eps(j) + eps(i))
               tmp3 = 1.0_wp/(eps(j) - eps(i) + damp)
               ! compute atomic response
               do m = 1, nat
                  weight = po(m, i)*pv(m, j)*tmp
                  response(m) = response(m) + weight
                  chempot(m) = chempot(m) + tmp2*weight
                  egap(m) = egap(m) + weight*tmp3
               end do ! at
            end if ! if virt
         end do !  virt
      end if ! if occ
   end do ! occ
   !!$omp end parallel do
   call timer%pop()
   call timer%push("AO B")
   if (print) then
      call ctx%message("")
      call ctx%message( "  -------------------------")
      call ctx%message( "  atomic frontier MO (beta) info ")
      call ctx%message( "  -------------------------")
      call ctx%message( "  atom   response (a.u.)   gap (eV)  chem.pot (eV)  HOAO (eV)    LUAO (eV)")
   end if

   ehoaob = 0.0_wp
   eluaob = 0.0_wp

   ! now get the atomic frontier orbitals
   tmp = 0.0_wp
   do m = 1, nat
      egap(m) = egap(m)/response(m)
      egap(m) = 1.0_wp/egap(m) - damp
      chempot(m) = chempot(m)/response(m)
      ehoaob(m) = (chempot(m) - 0.5_wp*egap(m))
      eluaob(m) = (chempot(m) + 0.5_wp*egap(m))

      ! fix cases for missing HOAO/LUAO
      if (sum(po(m, :)) .lt. 1.0d-12) then ! there is no occupied orbital for this atom
         !virt part  (including fractionally occupied)
         egap(m) = 0.0_wp
         chempot(m) = 0.0_wp
         do j = 1, nao
            virtb = (1.0_wp - foccb(j))
            if (virtb .gt. 0.0001_wp) then
               tmp = 1.0_wp/((eps(j))**2 + damp**2)
               tmp2 = eps(j)
               tmp3 = 1.0_wp/(eps(j) + 1.0d100 + damp)
               ! compute atomic response
               weight = pv(m, j)*tmp
               chempot(m) = chempot(m) + tmp2*weight
               egap(m) = egap(m) + weight*tmp3
            end if ! if virt
            egap(m) = 1.0_wp/(egap(m) + 1.0d-12) - damp
            eluaob(m) = chempot(m)
            ehoaob(m) = chempot(m) - egap(m)
         end do !  virt
      end if
      if (sum(pv(m, :)) .lt. 1.0d-12) then
         !occ part  (including fractionally occupied)
         egap(m) = 0.0_wp
         chempot(m) = 0.0_wp
         do j = 1, nao
            occb = foccb(j)
            if (occb .gt. 0.0001_wp) then
               tmp = 1.0_wp/((eps(j))**2 + damp**2)
               tmp2 = eps(j)
               tmp3 = 1.0_wp/(1.0d100 - eps(j) + damp)
               ! compute atomic response
               weight = po(m, j)*tmp
               chempot(m) = chempot(m) + tmp2*weight
               egap(m) = egap(m) + weight*tmp3
            end if ! if occ
            egap(m) = 1.0_wp/(egap(m) + 1.0d-12) - damp
            ehoaob(m) = chempot(m)
            eluaob(m) = chempot(m) + egap(m)
         end do !  occ
      end if
   end do

   call timer%pop()
   if (print) then
      do m = 1, nat
         write (tmp_str, '(3x,i6,5x,f8.4,5x,f8.4,5x,f8.4,5x,f8.4,5x,f8.4)') m, &
            response(m)*(autoev**2), egap(m), chempot(m), ehoaob(m), eluaob(m)
         call ctx%message(tmp_str)
      end do
      call ctx%message("")
   end if

   if (debug) then
      block
         
         integer :: it
         real(wp) :: ttime, stime
         character(len=*), parameter :: label(*) = [character(len=20):: &
            & "occupation a", "prep a", "AO A", "occupation b", "prep b", "AO B"]
         call ctx%message("")
         call ctx%message( "Frontier detail")
         call timer%pop()
         ttime = timer%get("total")
         call ctx%message(" total:"//repeat(" ", 16)//format_time(ttime))

         do it = 1, size(label)
            stime = timer%get(label(it))
            if (stime <= epsilon(0.0_wp)) cycle
            call ctx%message(" - "//label(it)//format_time(stime) &
               & //" ("//format_string(int(stime/ttime*100), '(i3)')//"%)")
         end do
         call ctx%message("")
      end block
   end if

   end subroutine atomic_frontier_orbitals
end module
