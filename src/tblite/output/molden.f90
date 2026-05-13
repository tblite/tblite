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

!> @file tblite/output/molden.f90
!> Provides routines for writing a Molden file with orbital and basis set information

!> This writer produces a Molden input with the following sections:
!>
!> ```
!> [Molden Format]
!> [Title]
!> [Cell] AU          non-standard extension for periodic systems
!> [Atoms] AU
!> [Core]             non-standard extension for core electron count
!> [Pseudo]           non-standard extension for valence electron count
!> [Nval]             non-standard extension for valence electron count
!> [GTO]
!> [6D]/[10F]/[15G]   optional explicit declaration of cartesian format
!> [MO]
!> ```
!> 
!> This format extends the official Molden specification documented at
!>
!>   https://www.theochem.ru.nl/molden/molden_format.html
!> 
!> The non-standard `[Cell]` section (consistent with Multiwfn and CP2K)
!> specifies periodic information lacking in the Molden Format. 
!> Contains cartesian lattice vectors assuming full three-dimensional periodicity:
!> 
!> ```
!> [Cell] AU
!>   ax ay az
!>   bx by bz
!>   cx cy cz
!> ```
!>
!> The number of valence electrons and the effective atomic charge is specified via
!> the non-standard `[Core]` section (consistent with Molden2AIM and PySCF), 
!> the non-standard `[Pseudo]` section (consistent with Molden and Orca),
!> the non-standard `[Nval]` section (consistent with CP2K). 
!>
!> The writer follows the official specification for the ordering of cartesian
!> atomic orbitals in the `[GTO]` section:
!>
!> ```
!> s    : s
!>
!> p    : px, py, pz
!>
!> 6D   : xx, yy, zz, xy, xz, yz
!>
!> 10F  : xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
!>
!> 15G  : xxxx, yyyy, zzzz, xxxy, xxxz, xyyy, yyyz, xzzz, yzzz,
!>        xxyy, xxzz, yyzz, xxyz, xyyz, xyzz
!> ```
!>
!> [6D]/[10F]/[15G] are specified to explicitly declare cartesian format. 

!> Implementation of Molden file output
module tblite_output_molden
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io_structure, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_version, only : get_tblite_version
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: molden_result

contains

   !> Writes the MOs, basis set, and geometry (including optional lattice parameters)
   !> to a Molden input file.
   subroutine molden_result(unit, mol, wfn, bas, error)
      !> Output file unit number for the Molden file
      integer, intent(in) :: unit
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Wavefunction structure data
      type(wavefunction_type), intent(in) :: wfn
      !> Basis set information
      class(basis_type), intent(in) :: bas
      !> Error handling
      type(error_type), intent(out), allocatable :: error

      integer :: ispin, ilat, iat, isp, ish, is, li, iprim, imo, ii
      integer :: iao_cart, jao, perm(15)
      character(len=1) :: aang
      character(len=:), allocatable :: version_string
      real(wp), allocatable :: coeff_cart(:, :, :)
      character(len=2) :: sym

      ! Convert MO coefficients from the spherical to cartesian AO basis
      allocate(coeff_cart(bas%nao_cart, bas%nao, wfn%nspin))
      do ispin = 1, wfn%nspin
         call bas%spherical_to_cartesian_trafo(mol, wfn%coeff(:, :, ispin), &
            & coeff_cart(:, :, ispin))
      end do

      ! Molden file header including version title
      write(unit,'(A)') '[Molden Format]'
      write(unit,'(A)') '[Title]'
      call get_tblite_version(string=version_string)
      write(unit,'(A)') 'tblite version '//trim(version_string)

      ! Non-standard periodic cell section
      if (any(mol%periodic)) then
         write(unit,'(A)') '[Cell] AU'
         do ilat = 1, 3
            write(unit,'(3(1x,ES24.16))') mol%lattice(:, ilat)
         end do
      end if

      ! Print symbols, atomic numbers, and coordinates
      write(unit,'(A)') '[Atoms] AU'
      do iat = 1, mol%nat
         isp = mol%id(iat)
         sym = trim(mol%sym(isp))
         ! Notation: symbol, atom, nuclear charge, coordinates
         write(unit,'(a2,1x,i6,1x,i6,3(1x,ES24.16))') sym, iat, mol%num(isp), &
            & mol%xyz(1, iat), mol%xyz(2, iat), mol%xyz(3, iat)
      end do

      ! Print number of core electrons, Molden2AIM / Multiwfn compatible
      write(unit,'(A)') '[Core]'
      do isp = 1, mol%nid
         sym = trim(mol%sym(isp))
         do iat = 1, mol%nat
            if (mol%id(iat) == isp) exit
         end do
         ! Notation: Element symbol : number of core electrons
         write(unit,'(a2,1x,a,1x,i0)') sym, ':', mol%num(isp) - nint(wfn%n0at(iat))
      end do

      ! Print the number of valence electrons, Molden / Orca compatible
      write(unit,'(A)') '[Pseudo]'
      do iat = 1, mol%nat
         isp = mol%id(iat)
         sym = trim(mol%sym(isp))
         ! Notation: Element symbol, atom number, number of valence electrons
         write(unit,'(a2,1x,i6,1x,i6)') sym, iat, nint(wfn%n0at(iat))
      end do

      write(unit,'(A)') '[Nval]'
      do isp = 1, mol%nid
         sym = trim(mol%sym(isp))
         do iat = 1, mol%nat
            if (mol%id(iat) == isp) exit
         end do
         ! Notation: Element symbol, number of valence electrons
         write(unit,'(a2,1x,i6)') sym, nint(wfn%n0at(iat))
      end do

      ! Print Gaussian-type orbital basis set data
      write(unit,'(A)') '[GTO]'
      do iat = 1, mol%nat
         write(unit,'(i6,a)') iat, ' 0'
         do ish = 1, bas%nsh_at(iat)
            associate(p_cgto => bas%cgto(ish, mol%id(iat)))
               select case(p_cgto%ang)
               case(0)
                  aang = 's'
               case(1)
                  aang = 'p'
               case(2)
                  aang = 'd'
               case(3)
                  aang = 'f'
               case(4)
                  aang = 'g'
               case default
                  call fatal_error(error, "Molden writer only supports angular momenta up to g")
                  return
               end select

               write(unit,'(a,1x,i6,1x,f8.2)') aang, p_cgto%nprim, 1.00_wp
               do iprim = 1, p_cgto%nprim
                  write(unit,'(2(1x,ES24.16))') &
                     & p_cgto%alpha(iprim), p_cgto%coeff(iprim)
               end do
            end associate
         end do
         write(unit,*)
      end do

      ! Explicitly declare Cartesian higher angular momentum functions
      if (bas%maxl >= 2) write(unit,'(A)') '[6D]'
      if (bas%maxl >= 3) write(unit,'(A)') '[10F]'
      if (bas%maxl >= 4) write(unit,'(A)') '[15G]'

      ! Print occupation numbers, orbital energies, and MO coefficients
      write(unit,'(A)') '[MO]'
      do ispin = 1, wfn%nspin

         do imo = 1, bas%nao
            ! Write symmetry (always `a` irreducible of C1)
            write(unit,'(A)', advance='no') 'Sym= '
            if (wfn%nspin == 2) then
               if (ispin == 1) then
                  write(unit,'(i5,a)') imo, 'a (alpha)'
               else
                  write(unit,'(i5,a)') imo, 'a (beta)'
               end if
            else
               write(unit,'(i5,a)') imo, 'a'
            end if

            ! Write the orbital energy
            write(unit,'(A)', advance='no') 'Ene= '
            write(unit,*) wfn%emo(imo, ispin)

            ! Write the spin channel
            write(unit,'(A)', advance='no') 'Spin= '
            if (wfn%nspin == 2) then
               if (ispin == 1) then
                  write(unit,'(A)') 'Alpha'
               else
                  write(unit,'(A)') 'Beta'
               end if
            else
               write(unit,'(A)') 'Alpha'
            end if

            ! Write the occupation numbers
            write(unit,'(A)', advance='no') 'Occup= '
            if (wfn%nspin == 1) then
               write(unit,'(F14.8)') wfn%focc(imo, 1) + wfn%focc(imo, 2)
            else
               write(unit,'(F14.8)') wfn%focc(imo, ispin)
            end if

            ! Write Cartesian coefficients for each shell in Molden order
            jao = 0
            do iat = 1, mol%nat
               isp = mol%id(iat)
               is = bas%ish_at(iat)
               do ish = 1, bas%nsh_at(iat)
                  li = bas%cgto(ish, isp)%ang
                  ii = bas%iao_cart_sh(is + ish)

                  ! Get index permutation for the current shell
                  ! to match Molden specific cartesian AO ordering
                  call get_molden_cart_perm(li, perm, error)
                  if (allocated(error)) return

                  do iao_cart = 1, bas%nao_cart_sh(is + ish)
                     jao = jao + 1
                     write(unit,'(i6,1x,ES24.16)') jao, &
                        & coeff_cart(ii + perm(iao_cart), imo, ispin)
                  end do
               end do
            end do

            if (jao /= bas%nao_cart) then
               call fatal_error(error, "Inconsistent AO count in Molden [MO] section")
               return
            end if
         end do
      end do

   end subroutine molden_result


   !> Return the index permutation to reach the Cartesian Molden ordering
   subroutine get_molden_cart_perm(l, perm, error)
      !> Angular momentum of the shell
      integer, intent(in)  :: l
      !> Index permutation relative to tblite internal ordering
      integer, intent(out) :: perm(15)
      !> Error handling
      type(error_type), intent(out), allocatable :: error

      perm = 0

      select case(l)
      case(0)
         ! s shell
         perm(1) = 1
      case(1)
         ! tblite Cartesian p order;
         ! py, pz, px
         ! Molden Cartesian p order:
         ! px, py, pz
         perm(1:3) = [3, 1, 2]
      case(2)
         ! tblite/Molden Cartesian d order:
         ! xx, yy, zz, xy, xz, yz
         perm(1:6) = [1, 2, 3, 4, 5, 6]
      case(3)
         ! tblite Cartesian f order:
         ! xxx, yyy, zzz, xxy, xxz, xyy, yyz, xzz, yzz, xyz
         ! Molden Cartesian f order:
         ! xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
         perm(1:10) = [1, 2, 3, 6, 4, 5, 8, 9, 7, 10]
      case(4)
         ! tblite/Molden Cartesian g order:
         ! xxxx, yyyy, zzzz, xxxy, xxxz, xyyy, yyyz, xzzz, yzzz,
         ! xxyy, xxzz, yyzz, xxyz, xyyz, xyzz
         perm(1:15) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
      case default
         call fatal_error(error, "Molden writer only supports angular momenta up to g")
         return
      end select

   end subroutine get_molden_cart_perm

end module tblite_output_molden