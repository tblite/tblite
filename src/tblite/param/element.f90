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

!> @file tblite/param/element.f90
!> Provides records for the element specific parameters

!> Definition of the element specific parameter records
module tblite_param_element
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io_symbols, only : to_number, symbol_length
   use tblite_data_paulingen, only : get_pauling_en
   use tblite_param_serde, only : serde_record
   use tblite_toml, only : toml_table, toml_array, get_value, set_value, add_array, len
   implicit none
   private

   public :: count


   !> The conversion factor from eV to Hartree is used for compatibility with older
   !> versions of xtb
   real(wp), parameter :: evtoau = 1.0_wp / 27.21138505_wp

   integer, parameter :: ngauss_default = 6
   real(wp), parameter :: kcn_default = 0.0_wp, lgam_default = 0.0_wp, &
      & mprad_default = 5.0_wp, xbond_default = 0.0_wp

   character(len=*), parameter :: k_shells = "shells", k_levels = "levels", &
      & k_shpoly = "shpoly", k_slater = "slater", k_refocc = "refocc", k_ngauss = "ngauss", &
      & k_gam = "gam", k_lgam = "lgam", k_gam3 = "gam3", k_kcn = "kcn", k_zeff = "zeff", &
      & k_arep = "arep", k_dkernel = "dkernel", k_qkernel = "qkernel", k_mprad = "mprad", &
      & k_mpvcn = "mpvcn", k_xbond = "xbond", k_en = "en"

   !> Representation of the element specific parameters
   type, public, extends(serde_record) :: element_record
      !> Element symbol of specie represented by this record
      character(len=symbol_length) :: sym = ''
      !> Atomic number of the specie represented by this record
      integer :: num = 0

      !> Effective nuclear charge used in repulsion
      real(wp) :: zeff = 0.0_wp
      !> Repulsion exponent
      real(wp) :: alpha = 0.0_wp

      !> Halogen bonding strength
      real(wp) :: xbond = 0.0_wp

      !> Number of valence and polarization shells
      integer :: nsh = 0
      !> Angular momentum for each shell
      integer, allocatable :: lsh(:)
      !> Principal quantum number for each shell
      integer, allocatable :: pqn(:)
      !> Number of primitive Gaussian functions used in the STO-NG expansion for each shell
      integer, allocatable :: ngauss(:)
      !> Atomic level energies for each shell
      real(wp), allocatable :: levels(:)
      !> Slater exponents of the STO-NG functions for each shell
      real(wp), allocatable :: slater(:)
      !> Reference occupation for each shell
      real(wp), allocatable :: refocc(:)
      !> CN dependent shift of the self energy for each shell
      real(wp), allocatable :: kcn(:)
      !> Polynomial enhancement for Hamiltonian elements
      real(wp), allocatable :: shpoly(:)

      !> Chemical hardness / Hubbard parameter
      real(wp) :: gam = 0.0_wp
      !> Relative chemical hardness for each shell
      real(wp), allocatable :: lgam(:)

      !> Atomic Hubbard derivative
      real(wp) :: gam3 = 0.0_wp

      !> Dipolar exchange-correlation kernel
      real(wp) :: dkernel = 0.0_wp
      !> Quadrupolar exchange-correlation kernel
      real(wp) :: qkernel = 0.0_wp
      !> Multipole radius
      real(wp) :: mprad = 0.0_wp
      !> Multipole valence CN
      real(wp) :: mpvcn = 0.0_wp

      !> Electronnegativity
      real(wp) :: en = 0.0_wp
   contains
      generic :: load => load_from_array
      generic :: dump => dump_to_array
      !> Read parametrization data from TOML data structure
      procedure :: load_from_toml
      !> Write parametrization data to TOML data structure
      procedure :: dump_to_toml
      !> Read parametrization data from parameter array
      procedure, private :: load_from_array
      !> Write parametrization data to parameter array
      procedure, private :: dump_to_array
   end type


   !> Masking for the element record
   type, public :: element_mask
      !> Element symbol of specie represented by this record
      character(len=symbol_length) :: sym = ''
      !> Atomic number of the specie represented by this record
      integer :: num = 0

      !> Effective nuclear charge used in repulsion
      logical :: zeff
      !> Repulsion exponent
      logical :: alpha

      !> Halogen bonding strength
      logical :: xbond

      !> Number of valence and polarization shells
      integer :: nsh
      !> Atomic level energies for each shell
      logical, allocatable :: levels(:)
      !> Slater exponents of the STO-NG functions for each shell
      logical, allocatable :: slater(:)
      !> CN dependent shift of the self energy for each shell
      logical, allocatable :: kcn(:)
      !> Polynomial enhancement for Hamiltonian elements
      logical, allocatable :: shpoly(:)

      !> Chemical hardness / Hubbard parameter
      logical :: gam
      !> Relative chemical hardness for each shell
      logical, allocatable :: lgam(:)

      !> Atomic Hubbard derivative
      logical :: gam3

      !> Dipolar exchange-correlation kernel
      logical :: dkernel
      !> Quadrupolar exchange-correlation kernel
      logical :: qkernel
   end type element_mask


   interface count
      module procedure :: count_mask
   end interface count


contains


!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   class(element_record), intent(inout) :: self
   type(toml_table), intent(inout) :: table
   type(error_type), allocatable, intent(out) :: error

   character(len=:), allocatable :: sym
   type(toml_array), pointer :: array
   integer :: stat

   call table%get_key(sym)
   self%sym = sym
   self%num = to_number(sym)

   call get_value(table, k_shells, array)
   if (.not.associated(array)) then
      call fatal_error(error, "No shells specified for "//trim(self%sym))
      return
   end if
   call get_shells(self, array, error)
   if (allocated(error)) return

   call get_value(table, k_levels, array)
   if (.not.associated(array)) then
      call fatal_error(error, "No atomic levels specified for "//trim(self%sym))
      return
   end if
   call get_levels(self, array, error)
   if (allocated(error)) return

   call get_value(table, k_slater, array)
   if (.not.associated(array)) then
      call fatal_error(error, "No slater exponents specified for "//trim(self%sym))
      return
   end if
   call get_slater(self, array, error)
   if (allocated(error)) return

   call get_value(table, k_refocc, array)
   if (.not.associated(array)) then
      call fatal_error(error, "No reference occupation specified for "//trim(self%sym))
      return
   end if
   call get_refocc(self, array, error)
   if (allocated(error)) return

   call get_value(table, k_ngauss, array, requested=.false.)
   if (associated(array)) then
      call get_ngauss(self, array, error)
      if (allocated(error)) return
   else
      allocate(self%ngauss(self%nsh), source=ngauss_default)
   end if

   call get_value(table, k_shpoly, array)
   if (.not.associated(array)) then
      call fatal_error(error, "No shell-polynomials specified for "//trim(self%sym))
      return
   end if
   call get_shpoly(self, array, error)
   if (allocated(error)) return

   call get_value(table, k_kcn, array, requested=.false.)
   if (associated(array)) then
      call get_kcn(self, array, error)
      if (allocated(error)) return
   else
      allocate(self%kcn(self%nsh), source=kcn_default)
   end if

   call get_value(table, k_gam, self%gam, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid atomic Hubbard parameter for "//trim(self%sym))
      return
   end if
   call get_value(table, k_lgam, array, requested=.false.)
   if (associated(array)) then
      call get_lgam(self, array, error)
      if (allocated(error)) return
   else
      allocate(self%lgam(self%nsh), source=lgam_default)
   end if

   call get_value(table, k_gam3, self%gam3, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid atomic Hubbard derivative for "//trim(self%sym))
      return
   end if

   call get_repulsion(self, table, error)
   if (allocated(error)) return

   call get_multipole(self, table, error)
   if (allocated(error)) return

   call get_value(table, k_xbond, self%xbond, xbond_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid halogen bonding strength for "//trim(self%sym))
      return
   end if

   call get_value(table, k_en, self%en, get_pauling_en(self%sym), stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid electronegativity for "//trim(self%sym))
      return
   end if

end subroutine load_from_toml

subroutine get_multipole(self, table, error)
   type(element_record), intent(inout) :: self
   type(toml_table), intent(inout) :: table
   type(error_type), allocatable, intent(out) :: error

   integer :: stat

   call get_value(table, k_dkernel, self%dkernel, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid dipole kernel for "//trim(self%sym))
      return
   end if
   call get_value(table, k_qkernel, self%qkernel, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid quadrupole kernel for "//trim(self%sym))
      return
   end if
   call get_value(table, k_mprad, self%mprad, mprad_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid multipole damping radius for "//trim(self%sym))
      return
   end if
   call get_value(table, k_mpvcn, self%mpvcn, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid multipole valence CN for "//trim(self%sym))
      return
   end if
end subroutine get_multipole

subroutine get_repulsion(self, table, error)
   type(element_record), intent(inout) :: self
   type(toml_table), intent(inout) :: table
   type(error_type), allocatable, intent(out) :: error

   integer :: stat

   call get_value(table, k_zeff, self%zeff, real(self%num, wp), stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid effective nuclear charge for "//trim(self%sym))
      return
   end if
   call get_value(table, k_arep, self%alpha, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid repulsion exponent for "//trim(self%sym))
      return
   end if
end subroutine get_repulsion

subroutine get_shells(self, array, error)
   type(element_record), intent(inout) :: self
   type(toml_array), intent(inout) :: array
   type(error_type), allocatable, intent(out) :: error

   integer :: i
   character(len=:), allocatable :: tmp

   self%nsh = len(array)
   if (self%nsh == 0) then
      call fatal_error(error, "No entries for "//trim(self%sym))
      return
   end if
   allocate(self%lsh(self%nsh))
   allocate(self%pqn(self%nsh))

   do i = 1, self%nsh
      call get_value(array, i, tmp)
      self%lsh(i) = get_lsh(tmp(2:2))
      self%pqn(i) = get_pqn(tmp(1:1))
   end do

   if (any(self%lsh < 0)) then
      call fatal_error(error, "Invalid angular momentum for "//trim(self%sym))
      return
   end if

   if (any(self%pqn < 1)) then
      call fatal_error(error, "Invalid principal quantum number for "//trim(self%sym))
      return
   end if

end subroutine get_shells

pure function get_pqn(str) result(pqn)
   character, intent(in) :: str
   integer :: pqn
   select case(str)
   case default
      pqn = 0
   case("1")
      pqn = 1
   case("2")
      pqn = 2
   case("3")
      pqn = 3
   case("4")
      pqn = 4
   case("5")
      pqn = 5
   case("6")
      pqn = 6
   end select
end function get_pqn

pure function get_lsh(str) result(lsh)
   character, intent(in) :: str
   integer :: lsh
   select case(str)
   case default
      lsh = -1
   case("s", "S")
      lsh = 0
   case("p", "P")
      lsh = 1
   case("d", "D")
      lsh = 2
   case("f", "F")
      lsh = 3
   case("g", "G")
      lsh = 4
   end select
end function get_lsh

subroutine get_levels(self, array, error)
   type(element_record), intent(inout) :: self
   type(toml_array), intent(inout) :: array
   type(error_type), allocatable, intent(out) :: error

   integer :: i
   real(wp) :: tmp

   if (len(array) < self%nsh) then
      call fatal_error(error, "Insufficient atomic level energies for "//trim(self%sym))
      return
   end if
   allocate(self%levels(self%nsh))

   do i = 1, self%nsh
      call get_value(array, i, tmp)
      self%levels(i) = tmp * evtoau
   end do
end subroutine get_levels

subroutine get_shpoly(self, array, error)
   type(element_record), intent(inout) :: self
   type(toml_array), intent(inout) :: array
   type(error_type), allocatable, intent(out) :: error

   integer :: i
   real(wp) :: tmp

   if (len(array) < self%nsh) then
      call fatal_error(error, "Insufficient shell-polynomials for "//trim(self%sym))
      return
   end if
   allocate(self%shpoly(self%nsh))

   do i = 1, self%nsh
      call get_value(array, i, tmp)
      self%shpoly(i) = tmp
   end do
end subroutine get_shpoly

subroutine get_slater(self, array, error)
   type(element_record), intent(inout) :: self
   type(toml_array), intent(inout) :: array
   type(error_type), allocatable, intent(out) :: error

   integer :: i
   real(wp) :: tmp

   if (len(array) < self%nsh) then
      call fatal_error(error, "Insufficient slater exponents for "//trim(self%sym))
      return
   end if
   allocate(self%slater(self%nsh))

   do i = 1, self%nsh
      call get_value(array, i, tmp)
      self%slater(i) = tmp
   end do

   if (any(self%slater < epsilon(0.0_wp))) then
      call fatal_error(error, "Invalid slater exponents for "//trim(self%sym))
   end if
end subroutine get_slater

subroutine get_refocc(self, array, error)
   type(element_record), intent(inout) :: self
   type(toml_array), intent(inout) :: array
   type(error_type), allocatable, intent(out) :: error

   integer :: i
   real(wp) :: tmp

   if (len(array) < self%nsh) then
      call fatal_error(error, "Insufficient reference occupations for "//trim(self%sym))
      return
   end if
   allocate(self%refocc(self%nsh))

   do i = 1, self%nsh
      call get_value(array, i, tmp)
      self%refocc(i) = tmp
   end do
end subroutine get_refocc

subroutine get_ngauss(self, array, error)
   type(element_record), intent(inout) :: self
   type(toml_array), intent(inout) :: array
   type(error_type), allocatable, intent(out) :: error

   integer :: i, tmp

   if (len(array) < self%nsh) then
      call fatal_error(error, "Insufficient primitive gaussians for "//trim(self%sym))
      return
   end if
   allocate(self%ngauss(self%nsh))

   do i = 1, self%nsh
      call get_value(array, i, tmp)
      self%ngauss(i) = tmp
   end do

   if (any(self%ngauss < 1) .or. any(self%ngauss > 6)) then
      call fatal_error(error, "Invalid number of primitive gaussian for "//trim(self%sym))
   end if
end subroutine get_ngauss

subroutine get_kcn(self, array, error)
   type(element_record), intent(inout) :: self
   type(toml_array), intent(inout) :: array
   type(error_type), allocatable, intent(out) :: error

   integer :: i
   real(wp) :: tmp

   if (len(array) < self%nsh) then
      call fatal_error(error, "Insufficient CN shift parameters for "//trim(self%sym))
      return
   end if
   allocate(self%kcn(self%nsh))

   do i = 1, self%nsh
      call get_value(array, i, tmp)
      self%kcn(i) = tmp * evtoau
   end do
end subroutine get_kcn

subroutine get_lgam(self, array, error)
   type(element_record), intent(inout) :: self
   type(toml_array), intent(inout) :: array
   type(error_type), allocatable, intent(out) :: error

   integer :: i
   real(wp) :: tmp

   if (len(array) < self%nsh) then
      call fatal_error(error, "Insufficient shell hardnesses for "//trim(self%sym))
      return
   end if
   allocate(self%lgam(self%nsh))

   do i = 1, self%nsh
      call get_value(array, i, tmp)
      self%lgam(i) = tmp
   end do
end subroutine get_lgam

!> Write parametrization data to TOML data structure
subroutine dump_to_toml(self, table, error)
   class(element_record), intent(in) :: self
   type(toml_table), intent(inout) :: table
   type(error_type), allocatable, intent(out) :: error

   integer :: i
   type(toml_array), pointer :: array
   character, parameter :: lsh(0:4) = ["s", "p", "d", "f", "g"]
   character, parameter :: pqn(1:6) = ["1", "2", "3", "4", "5", "6"]

   call add_array(table, k_shells, array)
   do i = 1, self%nsh
      call set_value(array, i, pqn(self%pqn(i))//lsh(self%lsh(i)))
   end do

   call add_array(table, k_levels, array)
   do i = 1, self%nsh
      call set_value(array, i, self%levels(i) / evtoau)
   end do

   call add_array(table, k_slater, array)
   do i = 1, self%nsh
      call set_value(array, i, self%slater(i))
   end do

   call add_array(table, k_ngauss, array)
   do i = 1, self%nsh
      call set_value(array, i, self%ngauss(i))
   end do

   call add_array(table, k_refocc, array)
   do i = 1, self%nsh
      call set_value(array, i, self%refocc(i))
   end do

   call add_array(table, k_shpoly, array)
   do i = 1, self%nsh
      call set_value(array, i, self%shpoly(i))
   end do

   call add_array(table, k_kcn, array)
   do i = 1, self%nsh
      call set_value(array, i, self%kcn(i) / evtoau)
   end do

   call set_value(table, k_gam, self%gam)
   call add_array(table, k_lgam, array)
   do i = 1, self%nsh
      call set_value(array, i, self%lgam(i))
   end do

   call set_value(table, k_gam3, self%gam3)

   call set_value(table, k_zeff, self%zeff)
   call set_value(table, k_arep, self%alpha)
   call set_value(table, k_xbond, self%xbond)
   call set_value(table, k_en, self%en)

   call set_value(table, k_dkernel, self%dkernel)
   call set_value(table, k_qkernel, self%qkernel)
   call set_value(table, k_mprad, self%mprad)
   call set_value(table, k_mpvcn, self%mpvcn)

end subroutine dump_to_toml

!> Read parametrization data from parameter array
subroutine load_from_array(self, array, offset, base, mask, error)
   class(element_record), intent(inout) :: self
   real(wp), intent(in) :: array(:)
   integer, intent(inout) :: offset
   type(element_record), intent(in) :: base
   type(element_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

   select type(self)
   type is (element_record)
      self = base
   end select

   call load_atom_par(self%zeff, mask%zeff, array, offset)
   call load_atom_par(self%alpha, mask%alpha, array, offset)

   call load_atom_par(self%xbond, mask%xbond, array, offset)

   call load_shell_par(self%levels, mask%levels, array, offset, scale=evtoau)
   call load_shell_par(self%slater, mask%slater, array, offset)
   call load_shell_par(self%kcn, mask%kcn, array, offset, scale=evtoau)
   call load_shell_par(self%shpoly, mask%shpoly, array, offset, scale=0.01_wp)

   call load_atom_par(self%gam, mask%gam, array, offset)
   call load_shell_par(self%lgam, mask%lgam, array, offset)

   call load_atom_par(self%gam3, mask%gam3, array, offset, scale=0.1_wp)

   call load_atom_par(self%dkernel, mask%dkernel, array, offset, scale=0.01_wp)
   call load_atom_par(self%qkernel, mask%qkernel, array, offset, scale=0.01_wp)
end subroutine load_from_array

!> Write parametrization data to parameter array
subroutine dump_to_array(self, array, offset, mask, error)
   class(element_record), intent(in) :: self
   real(wp), intent(inout) :: array(:)
   integer, intent(inout) :: offset
   type(element_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

   call dump_atom_par(self%zeff, mask%zeff, array, offset)
   call dump_atom_par(self%alpha, mask%alpha, array, offset)

   call dump_atom_par(self%xbond, mask%xbond, array, offset)

   call dump_shell_par(self%levels, mask%levels, array, offset, scale=evtoau)
   call dump_shell_par(self%slater, mask%slater, array, offset)
   call dump_shell_par(self%kcn, mask%kcn, array, offset, scale=evtoau)
   call dump_shell_par(self%shpoly, mask%shpoly, array, offset, scale=0.01_wp)

   call dump_atom_par(self%gam, mask%gam, array, offset)
   call dump_shell_par(self%lgam, mask%lgam, array, offset)

   call dump_atom_par(self%gam3, mask%gam3, array, offset, scale=0.1_wp)

   call dump_atom_par(self%dkernel, mask%dkernel, array, offset, scale=0.01_wp)
   call dump_atom_par(self%qkernel, mask%qkernel, array, offset, scale=0.01_wp)
end subroutine dump_to_array

pure subroutine load_atom_par(par, mask, array, ii, scale)
   real(wp), intent(inout) :: par
   logical, intent(in) :: mask
   real(wp), intent(in) :: array(:)
   integer, intent(inout) :: ii
   real(wp), intent(in), optional :: scale
   real(wp) :: scale_
   scale_ = 1.0_wp
   if (present(scale)) scale_ = scale
   if (mask) then
      ii = ii+1
      par = array(ii) * scale_
   end if
end subroutine load_atom_par

pure subroutine load_shell_par(par, mask, array, ii, scale)
   real(wp), intent(inout) :: par(:)
   logical, intent(in) :: mask(:)
   real(wp), intent(in) :: array(:)
   integer, intent(inout) :: ii
   real(wp), intent(in), optional :: scale
   real(wp) :: scale_
   integer :: ish
   scale_ = 1.0_wp
   if (present(scale)) scale_ = scale
   do ish = 1, size(par)
      if (mask(ish)) then
         ii = ii+1
         par(ish) = array(ii) * scale_
      end if
   end do
end subroutine load_shell_par

pure subroutine dump_atom_par(par, mask, array, ii, scale)
   real(wp), intent(in) :: par
   logical, intent(in) :: mask
   real(wp), intent(inout) :: array(:)
   integer, intent(inout) :: ii
   real(wp), intent(in), optional :: scale
   real(wp) :: scale_
   scale_ = 1.0_wp
   if (present(scale)) scale_ = scale
   if (mask) then
      ii = ii+1
      array(ii) = par / scale_
   end if
end subroutine dump_atom_par

pure subroutine dump_shell_par(par, mask, array, ii, scale)
   real(wp), intent(in) :: par(:)
   logical, intent(in) :: mask(:)
   real(wp), intent(inout) :: array(:)
   integer, intent(inout) :: ii
   real(wp), intent(in), optional :: scale
   real(wp) :: scale_
   integer :: ish
   scale_ = 1.0_wp
   if (present(scale)) scale_ = scale
   do ish = 1, size(par)
      if (mask(ish)) then
         ii = ii+1
         array(ii) = par(ish) / scale_
      end if
   end do
end subroutine dump_shell_par

elemental function count_mask(mask) result(ncount)
   type(element_mask), intent(in) :: mask
   integer :: ncount
   ncount = count([ &
      mask%zeff, &
      mask%alpha, &
      mask%xbond, &
      mask%levels, &
      mask%slater, &
      mask%kcn, &
      mask%shpoly, &
      mask%gam, &
      mask%lgam, &
      mask%gam3, &
      mask%dkernel, &
      mask%qkernel])
end function count_mask

end module tblite_param_element
