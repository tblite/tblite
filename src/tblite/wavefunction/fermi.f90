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

!> @file tblite/wavefunction/fermi.f90
!> Provides a Fermi based electronic filling function

!> Implementation of a Fermi distribution for filling the electronic levels
module tblite_wavefunction_fermi
   use mctc_env, only : wp
   implicit none
   private

   public :: get_fermi_filling


contains

subroutine get_fermi_filling(nel, kt, emo, homo, focc, e_fermi)
   !> Number of electrons
   real(wp), intent(in) :: nel
   !> Molecular orbital energy levels
   real(wp), intent(in) :: emo(:)
   !> Electronic temperature in atomic units
   real(wp), intent(in) :: kt
   !> Index of the highest occupied molecular orbital
   integer, intent(out) :: homo
   !> Occupation numbers
   real(wp), intent(out) :: focc(:)
   !> Fermi energy
   real(wp), intent(out) :: e_fermi

   real(wp) :: etmp

   e_fermi = 0.0_wp

   call get_aufbau_filling(nel, homo, focc)

   ! Optimize the Fermilevel if there is a finite temperature
   if (nel > 0.0_wp .and. kt > 0.0_wp) then
      call get_fermi_filling_(nel, homo, kt, emo, focc, e_fermi)
   else
      e_fermi = 0.0_wp
   end if

end subroutine get_fermi_filling

subroutine get_aufbau_filling(nel, homo, occ)
   !> Number of electrons
   real(wp), intent(in) :: nel
   !> Index of the highest occupied molecular orbital
   integer, intent(out) :: homo
   !> Occupation numbers
   real(wp), intent(out) :: occ(:)

   occ(:) = 0.0_wp
   homo = floor(nel)
   occ(:min(homo, size(occ))) = 1.0_wp
   if (homo < size(occ)) occ(homo+1) = mod(nel, 1.0_wp)
   homo = merge(homo+1, homo, mod(nel, 1.0_wp) > 0.5_wp)
end subroutine get_aufbau_filling

subroutine get_fermi_filling_(nel, homo, kt, emo, occ, e_fermi)
   !> Number of electrons
   real(wp), intent(in) :: nel
   !> Index of the highest occupied molecular orbital
   integer, intent(in) :: homo
   !> Molecular orbital energy levels
   real(wp), intent(in) :: emo(:)
   !> Electronic temperature in atomic units
   real(wp), intent(in) :: kt
   !> Occupation numbers
   real(wp), intent(out) :: occ(:)
   !> Fermi energy
   real(wp), intent(out) :: e_fermi

   real(wp) :: arg, total_number, total_dfermi, dfermifunct, fermifunct, change_fermi
   integer :: ncycle, iao

   integer, parameter :: max_cycle = 200
   real(wp), parameter :: thr = min(sqrt(epsilon(1.0_wp)), 1e+5_wp*epsilon(1.0_wp))
   real(wp), parameter :: sqrttiny = sqrt(tiny(1.0_wp))

   e_fermi = 0.5_wp * (emo(max(homo, 1)) + emo(min(homo + 1, size(emo))))
   do ncycle = 1, max_cycle
      total_number = 0.0_wp
      total_dfermi = 0.0_wp

      do iao = 1, size(emo)
         arg = (emo(iao) - e_fermi) / kt

         if (arg > 50.0_wp) then
            fermifunct = 0.0_wp
            dfermifunct = 0.0_wp
         else if (arg < -50.0_wp) then
            fermifunct = 1.0_wp
            dfermifunct = 0.0_wp
         else
            fermifunct = 1.0_wp / (exp(arg) + 1.0_wp)
            dfermifunct = fermifunct * (1.0_wp - fermifunct) / kt
         end if
         occ(iao) = fermifunct
         total_number = total_number + fermifunct
         total_dfermi = total_dfermi + dfermifunct
      end do
      if (abs(nel - total_number) <= thr) exit
      if (total_dfermi > sqrttiny) then
         change_fermi = (nel - total_number) / total_dfermi
         e_fermi = e_fermi + change_fermi
      else
         exit
      end if
   end do

end subroutine get_fermi_filling_

end module tblite_wavefunction_fermi
