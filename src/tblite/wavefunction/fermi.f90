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

module tblite_wavefunction_fermi
   use mctc_env, only : wp
   implicit none
   private

   public :: get_alpha_beta_occupation, get_fermi_filling


contains

subroutine get_fermi_filling(nocc, nuhf, kt, emo, homoa, homob, focc, e_fermi, ts)
   real(wp), intent(in) :: nocc
   real(wp), intent(in) :: nuhf
   real(wp), intent(in) :: emo(:)
   real(wp), intent(in) :: kt
   integer, intent(out) :: homoa
   integer, intent(out) :: homob
   real(wp), intent(out) :: focc(:)
   real(wp), intent(out) :: e_fermi
   real(wp), intent(out) :: ts

   real(wp) :: nalp, nbet, etmp, stmp
   real(wp), allocatable :: occt(:)
   integer :: homo

   call get_alpha_beta_occupation(nocc, nuhf, nalp, nbet)

   allocate(occt(size(focc)))
   ts = 0.0_wp
   e_fermi = 0.0_wp
   focc(:) = 0.0_wp
   occt(:) = 0.0_wp
   homo = floor(nalp)
   occt(:homo) = 1.0_wp
   if (homo < size(focc)) occt(homo+1) = mod(nalp, 1.0_wp)
   homoa = merge(homo+1, homo, mod(nalp, 1.0_wp) > 0.5_wp)

   if (homoa > 0) then
      call get_fermi_filling_(homoa, kt, emo, occt, etmp)
      call get_electronic_entropy(occt, kt, stmp)
      focc(:) = focc + occt
      e_fermi = 0.5_wp * etmp
      ts = ts + stmp
   end if

   occt(:) = 0.0_wp
   homo = floor(nbet)
   occt(:homo) = 1.0_wp
   if (homo < size(focc)) occt(homo+1) = mod(nbet, 1.0_wp)
   homob = merge(homo+1, homo, mod(nbet, 1.0_wp) > 0.5_wp)

   if (homob > 0) then
      call get_fermi_filling_(homob, kt, emo, occt, etmp)
      call get_electronic_entropy(occt, kt, stmp)
      focc(:) = focc + occt
      e_fermi = e_fermi + 0.5_wp * etmp
      ts = ts + stmp
   end if

end subroutine get_fermi_filling

subroutine get_fermi_filling_(homo, kt, emo, occ, e_fermi)
   integer, intent(in) :: homo
   real(wp), intent(in) :: emo(:)
   real(wp), intent(in) :: kt
   real(wp), intent(out) :: occ(:)
   real(wp), intent(out) :: e_fermi

   real(wp) :: occt, total_number
   real(wp) :: total_dfermi, dfermifunct, fermifunct, change_fermi
   integer, parameter :: max_cycle = 200
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   integer :: ncycle, iao

   e_fermi = 0.5*(emo(max(homo, 1))+emo(min(homo+1, size(emo))))
   occt = homo

   do ncycle = 1, max_cycle
      total_number = 0.0
      total_dfermi = 0.0
      do iao = 1, size(emo)
         fermifunct = 0.0
         if((emo(iao)-e_fermi)/kt.lt.50) then
            fermifunct = 1.0/(exp((emo(iao)-e_fermi)/kt)+1.0)
            dfermifunct = exp((emo(iao)-e_fermi)/kt) / &
               & (kt*(exp((emo(iao)-e_fermi)/kt)+1.0)**2)
         else
            dfermifunct = 0.0
         end if
         occ(iao) = fermifunct
         total_number = total_number + fermifunct
         total_dfermi = total_dfermi + dfermifunct
      end do
      change_fermi = (occt-total_number)/total_dfermi
      e_fermi = e_fermi+change_fermi
      if (abs(occt-total_number) <= thr) exit
   end do

end subroutine get_fermi_filling_

subroutine get_electronic_entropy(occ, kt, s)
   real(wp), intent(in) :: occ(:)
   real(wp), intent(in) :: kt
   real(wp), intent(out) :: s

   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   integer :: iao

   s = 0.0_wp
   do iao = 1, size(occ)
      if (occ(iao) > thr .and. 1.0_wp - occ(iao) > thr) then
         s = s + occ(iao)*log(occ(iao)) + (1.0_wp - occ(iao))*log(1.0_wp - occ(iao))
      end if
   enddo
   s=s*kt

end subroutine get_electronic_entropy


!> Split an real occupation number into alpha and beta space.
!>
!> This routine does not perform any checks on the condition
!> ``mod(nocc, 2) == 0 .eqv. mod(nuhf, 2) == 0`` and will yield fractional
!> occupations in case those condtions are not fullfilled.
!> However, it will avoid creating negative occupation numbers.
subroutine get_alpha_beta_occupation(nocc, nuhf, nalp, nbet)
   real(wp), intent(in) :: nocc
   real(wp), intent(in) :: nuhf
   real(wp), intent(out) :: nalp
   real(wp), intent(out) :: nbet

   real(wp) :: ntmp, diff

   ! make sure we cannot get a negative occupation here
   diff = min(nuhf, nocc)
   ntmp = nocc - diff

   nalp = ntmp / 2 + diff
   nbet = ntmp / 2
end subroutine get_alpha_beta_occupation

end module tblite_wavefunction_fermi
