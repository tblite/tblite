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

!> @dir tblite/output
!> Contains routines for terminal or string output generation

!> @file tblite/output/ascii.f90
!> Provides routines for terminal or text file output

!> Implementation of terminal or text file output
module tblite_output_ascii
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_convert, only : autoev
   use tblite_version, only : get_tblite_version
   implicit none
   private

   public :: ascii_levels, ascii_dipole_moments, ascii_quadrupole_moments
   public :: json_results, tagged_result, ascii_atomic_charges

contains


subroutine ascii_levels(unit, verbosity, emo, focc, range)
   integer, intent(in) :: unit
   integer, intent(in) :: verbosity
   real(wp), intent(in) :: emo(:, :)
   real(wp), intent(in) :: focc(:, :)
   integer, intent(in) :: range

   character(len=*), parameter :: hlfmt = '(a21, f21.7, 1x, "Eh", f18.4, 1x, "eV")'
   integer :: nao, maxorb, minorb, iorb, spin, homo
   real(wp) :: gap, nel

   do spin = 1, size(emo, 2)
   nel = sum(focc(:, spin))
   homo = floor(nel)
   homo = merge(homo+1, homo, mod(nel, 1.0_wp) > 0.5_wp)
   nao = size(emo, 1)
   minorb = max(homo - (range+1), 1)
   maxorb = min(homo +  range, nao)
   gap = emo(min(homo+1, nao), spin) - emo(max(homo, 1), spin)

   write(unit, '(66("-"))')
   write(unit, '(a7, a14, a21, a21)') "#", "Occupation", "Energy/Eh", "Energy/eV"
   write(unit, '(66("-"))')
   if (minorb > 1) then
      call write_line(1, focc(:, spin), emo(:, spin), homo)
      if (minorb > 2) &
         write(unit, '(a7, a14, a21, a21)') "...", "...", "...", "..."
   endif
   do iorb = minorb, maxorb
      call write_line(iorb, focc(:, spin), emo(:, spin), homo)
   enddo
   if (maxorb < nao) then
      if (maxorb < nao-1) then
         if (focc(maxorb, spin) > sqrt(epsilon(1.0_wp))) then
            write(unit, '(a7, a14, a21, a21)') "...", "...", "...", "..."
         else
            write(unit, '(a7, a14, a21, a21)') "...", "", "...", "..."
         endif
      endif
      call write_line(nao, focc(:, spin), emo(:, spin), homo)
   endif
   write(unit, '(66("-"))')
   write(unit, hlfmt) "HL-Gap", gap, gap*autoev
   write(unit, '(66("-"))')
   write(unit, '(a)')
   end do

contains
   subroutine write_line(iorb, focc, emo, ihomo)
      integer, intent(in) :: iorb
      integer, intent(in) :: ihomo
      real(wp), intent(in) :: focc(:)
      real(wp), intent(in) :: emo (:)
      character(len=*), parameter :: mofmt = '(i7, f14.4, f21.7, f21.4)'
      character(len=*), parameter :: vofmt = '(i7, 14x,  f21.7, f21.4)'
      if (focc(iorb) < 1.0e-7_wp) then
         write(unit, vofmt, advance='no') iorb,             emo(iorb), emo(iorb)*autoev
      else
         write(unit, mofmt, advance='no') iorb, focc(iorb), emo(iorb), emo(iorb)*autoev
      endif
      if (iorb == ihomo) then
         write(unit, '(1x, "(HOMO)")')
      elseif (iorb == ihomo+1) then
         write(unit, '(1x, "(LUMO)")')
      else
         write(unit, '(a)')
      endif
   end subroutine write_line
end subroutine ascii_levels

subroutine ascii_dipole_moments(unit, verbosity, mol, dpat, dpmom)
   integer, intent(in) :: unit
   integer, intent(in) :: verbosity
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: dpat(:, :)
   real(wp), intent(in) :: dpmom(:)

   integer :: iat, isp

   write(unit, '(a,":")') "Atomic dipole moments (in atomic units)"
   write(unit, '(57("-"))')
   write(unit, '(a6, 1x, a4, 5x, *(1x, a12))') &
      & "#", "Z", "x", "y", "z"
   write(unit, '(57("-"))')
   do iat = 1, mol%nat
      isp = mol%id(iat)
      write(unit, '(i6, 1x, i4, 1x, a4, *(1x, f12.6))') &
         & iat, mol%num(isp), mol%sym(isp), dpat(:, iat)
   end do
   write(unit, '(57("-"))')
   write(unit, '(1x, a15, *(1x, f12.6))') "total", dpmom
   write(unit, '(57("-"))')
   write(unit, '(a)')
end subroutine ascii_dipole_moments

subroutine ascii_atomic_charges(unit, verbosity, mol, qat)
   integer, intent(in) :: unit
   integer, intent(in) :: verbosity
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: qat(:)

   integer :: iat, isp

   write(unit, '(a,":")') "Atomic charges (in atomic units)"
   write(unit, '(57("-"))')
   write(unit, '(a6, 1x, a4, 7x, *(1x, a10))') &
      & "#", "Z", "q"
   write(unit, '(57("-"))')
   do iat = 1, mol%nat
      isp = mol%id(iat)
      write(unit, '(i6, 1x, i4, 1x, a4, *(1x, f12.6))') &
         & iat, mol%num(isp), mol%sym(isp), qat(iat)
   end do
   write(unit, '(57("-"))')
   write(unit, '(a)')
end subroutine ascii_atomic_charges

subroutine ascii_quadrupole_moments(unit, verbosity, mol, qpat, qpmom)
   integer, intent(in) :: unit
   integer, intent(in) :: verbosity
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: qpat(:, :)
   real(wp), intent(in) :: qpmom(:)

   integer :: iat, isp

   write(unit, '(a,":")') "Atomic quadrupole moments (in atomic units)"
   write(unit, '(83("-"))')
   write(unit, '(a6, 1x, a4, 5x, *(1x, a10))') &
      & "#", "Z", "xx", "xy", "yy", "xz", "yz", "zz"
   write(unit, '(83("-"))')
   do iat = 1, mol%nat
      isp = mol%id(iat)
      write(unit, '(i6, 1x, i4, 1x, a4, *(1x, f10.4))') &
         & iat, mol%num(isp), mol%sym(isp), qpat(:, iat)
   end do
   write(unit, '(83("-"))')
   write(unit, '(1x, a15, *(1x, f10.4))') "total", qpmom
   write(unit, '(83("-"))')
   write(unit, '(a)')
end subroutine ascii_quadrupole_moments


subroutine json_results(unit, indentation, energy, gradient, sigma, energies, charges)
   integer, intent(in) :: unit
   character(len=*), intent(in), optional :: indentation
   real(wp), intent(in), optional :: energy
   real(wp), intent(in), optional :: gradient(:, :)
   real(wp), intent(in), optional :: sigma(:, :)
   real(wp), intent(in), optional :: energies(:)
   real(wp), intent(in), optional :: charges(:)
   character(len=:), allocatable :: indent, version_string
   character(len=*), parameter :: jsonkey = "('""',a,'"":',1x)"
   real(wp), allocatable :: array(:)

   call get_tblite_version(string=version_string)

   if (present(indentation)) then
      indent = indentation
   end if

   write(unit, '("{")', advance='no')
   if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
   write(unit, jsonkey, advance='no') 'version'
   write(unit, '(1x,a)', advance='no') '"'//version_string//'"'
   if (present(energy)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'energy'
      write(unit, '(1x,es25.16)', advance='no') energy
   end if
   if (present(energies)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'energies'
      call write_json_array(unit, energies, indent)
   end if
   if (present(sigma)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'virial'
      array = reshape(sigma, [size(sigma)])
      call write_json_array(unit, array, indent)
   end if
   if (present(gradient)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'gradient'
      array = reshape(gradient, [size(gradient)])
      call write_json_array(unit, array, indent)
   end if
   if (present(charges)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'charges'
      array = reshape(charges, [size(charges)])
      call write_json_array(unit, array, indent)
   end if
   if (allocated(indent)) write(unit, '(/)', advance='no')
   write(unit, '("}")')

end subroutine json_results


subroutine write_json_array(unit, array, indent)
   integer, intent(in) :: unit
   real(wp), intent(in) :: array(:)
   character(len=:), allocatable, intent(in) :: indent
   integer :: i
   write(unit, '("[")', advance='no')
   do i = 1, size(array)
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 2)
      write(unit, '(es23.16)', advance='no') array(i)
      if (i /= size(array)) write(unit, '(",")', advance='no')
   end do
   if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
   write(unit, '("]")', advance='no')
end subroutine write_json_array


subroutine tagged_result(unit, energy, gradient, sigma, energies)
   integer, intent(in) :: unit
   real(wp), intent(in), optional :: energy
   real(wp), intent(in), optional :: gradient(:, :)
   real(wp), intent(in), optional :: sigma(:, :)
   real(wp), intent(in), optional :: energies(:)
   character(len=*), parameter :: tag_header = &
      & '(a,t20,":",a,":",i0,":",*(i0:,","))'

   if (present(energy)) then
      write(unit, tag_header) "energy", "real", 0
      write(unit, '(3es24.16)') energy
   end if
   if (present(energies)) then
      write(unit, tag_header) "energies", "real", 1, shape(energies)
      write(unit, '(3es24.16)') energies
   end if
   if (present(gradient)) then
      write(unit, tag_header) "gradient", "real", 2, shape(gradient)
      write(unit, '(3es24.16)') gradient
   end if
   if (present(sigma)) then
      write(unit, tag_header) "virial", "real", 2, shape(sigma)
      write(unit, '(3es24.16)') sigma
   end if

end subroutine tagged_result


end module tblite_output_ascii
