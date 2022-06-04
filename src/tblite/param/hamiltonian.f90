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

!> @file tblite/param/hamiltonian.f90
!> Provides parameter record for the Hamiltonian

!> Defines model for the Hamiltonian parameters
module tblite_param_hamiltonian
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io_symbols, only : symbol_length
   use tblite_param_serde, only : serde_record
   use tblite_toml, only : toml_table, get_value, set_value, add_table
   implicit none
   private

   public :: count

   character(len=*), parameter :: k_xtb = "xtb", k_ang(0:4) = ["s", "p", "d", "f", "g"], &
      & k_enscale = "enscale", k_shell = "shell", k_wexp = "wexp", k_cn = "cn", &
      & k_kpol = "kpol", k_kpair = "kpair"
   real(wp), parameter :: kpol_default = 2.0_wp


   !> Hamiltonian parametrization record
   type, public, extends(serde_record) :: hamiltonian_record
      character(len=symbol_length), allocatable :: sym(:)
      character(len=:), allocatable :: cn
      real(wp), allocatable :: kpair(:, :)
      real(wp) :: ksh(0:4, 0:4)
      real(wp) :: kpol
      real(wp) :: enscale
      real(wp) :: wexp
      integer :: lmax
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


   !> Masking for the Hamiltonian parameter record
   type, public :: hamiltonian_mask
   end type hamiltonian_mask


   interface count
      module procedure :: count_mask
   end interface count


contains


!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   !> Instance of the parametrization data
   class(hamiltonian_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child
   integer :: stat

   call get_value(table, k_xtb, child, requested=.false.)
   if (.not.associated(child)) then
      call fatal_error(error, "No entry for xTB Hamiltonian found")
      return
   end if

   call get_value(child, k_cn, self%cn, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for CN type in Hamiltonian")
      return
   end if

   call get_value(child, k_wexp, self%wexp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for Slater exponent weighting of Hamiltonian")
      return
   end if

   call get_value(child, k_kpol, self%kpol, kpol_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for scaling of non-valence Hamiltonian elements")
      return
   end if

   call get_ksh(self, child, error)
   if (allocated(error)) return

   call get_value(child, k_enscale, self%enscale, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for EN scaling parameter in Hamiltonian")
      return
   end if

   call get_kpair(self, child, error)
   if (allocated(error)) return

end subroutine load_from_toml


subroutine get_kpair(self, table, error)
   !> Instance of the parametrization data
   type(hamiltonian_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child
   character(len=:), allocatable :: key
   real(wp) :: tmp
   integer :: ir, jr

   call get_value(table, k_kpair, child)

   allocate(self%kpair(size(self%sym), size(self%sym)))
   self%kpair(:, :) = 1.0_wp

   do ir = 1, size(self%sym)
      do jr = 1, ir
         key = trim(self%sym(jr))//"-"//trim(self%sym(ir))
         if (.not.child%has_key(key)) key = trim(self%sym(ir))//"-"//trim(self%sym(jr))
         if (child%has_key(key)) then
            call get_value(child, key, tmp, 1.0_wp)
            self%kpair(jr, ir) = tmp
            self%kpair(ir, jr) = tmp
         end if
      end do
   end do

end subroutine get_kpair


subroutine get_ksh(self, table, error)
   !> Instance of the parametrization data
   type(hamiltonian_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child
   character(len=:), allocatable :: key
   real(wp), allocatable :: last, kav
   integer :: l, k, stat, lmax

   call get_value(table, k_shell, child, requested=.false.)
   if (associated(child)) then
      do l = 0, 4
         key = repeat(k_ang(l), 2)
         if (.not.child%has_key(key)) then
            if (allocated(last)) then
               self%ksh(l, l) = last
               cycle
            end if
            call fatal_error(error, "No entry for "//key//"-shell pair scaling parameter")
            exit
         end if
         call get_value(child, key, self%ksh(l, l), stat=stat)
         if (stat /= 0) then
            call fatal_error(error, "Cannot read "//key//"-shell pair scaling parameter")
            exit
         end if
         if (stat == 0) then
            last = self%ksh(l, l)
            lmax = l
         end if
      end do
      if (allocated(error)) return

      do l = 0, 4
         do k = 0, l-1
            key = k_ang(k)//k_ang(l)
            kav = 0.5_wp * (self%ksh(k, k) + self%ksh(l, l))
            call get_value(child, key, self%ksh(k, l), kav, stat=stat)
            if (stat /= 0) then
               call fatal_error(error, "Cannot read "//key//"-shell pair scaling parameter")
               exit
            end if
            self%ksh(l, k) = self%ksh(k, l)
         end do
      end do
   end if

   self%lmax = lmax

end subroutine get_ksh

!> Write parametrization data to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(hamiltonian_record), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child, sub
   character(len=:), allocatable :: key
   integer :: l, k, ir, jr

   call add_table(table, k_xtb, child)
   call set_value(child, k_wexp, self%wexp)
   call set_value(child, k_kpol, self%kpol)
   call set_value(child, k_enscale, self%enscale)
   if (allocated(self%cn)) call set_value(child, k_cn, self%cn)
   call add_table(child, k_shell, sub)
   do l = 0, self%lmax
      call set_value(sub, repeat(k_ang(l), 2), self%ksh(l, l))
   end do
   do l = 0, self%lmax
      do k = 0, l-1
         if (abs(self%ksh(k, l) - (self%ksh(l, l) + self%ksh(k, k))/2) > epsilon(1.0_wp)) then
            call set_value(sub, k_ang(k)//k_ang(l), self%ksh(k, l))
         end if
      end do
   end do

   if (any(abs(self%kpair - 1.0_wp) > epsilon(1.0_wp))) then
      call add_table(child, k_kpair, sub)
      do ir = 1, size(self%sym)
         do jr = 1, ir
            if (abs(self%kpair(jr, ir) - 1.0_wp) > epsilon(1.0_wp)) then
               key = trim(self%sym(ir))//"-"//trim(self%sym(jr))
               call set_value(sub, key, self%kpair(jr, ir))
            end if
         end do
      end do
   end if
end subroutine dump_to_toml


!> Read parametrization data from parameter array
subroutine load_from_array(self, array, offset, base, mask, error)
   class(hamiltonian_record), intent(inout) :: self
   real(wp), intent(in) :: array(:)
   integer, intent(inout) :: offset
   type(hamiltonian_record), intent(in) :: base
   type(hamiltonian_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

   select type(self)
   type is (hamiltonian_record)
      self = base
   end select

end subroutine load_from_array

!> Write parametrization data to parameter array
subroutine dump_to_array(self, array, offset, mask, error)
   class(hamiltonian_record), intent(in) :: self
   real(wp), intent(inout) :: array(:)
   integer, intent(inout) :: offset
   type(hamiltonian_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

end subroutine dump_to_array

elemental function count_mask(mask) result(ncount)
   type(hamiltonian_mask), intent(in) :: mask
   integer :: ncount
   ncount = 0
end function count_mask


end module tblite_param_hamiltonian
