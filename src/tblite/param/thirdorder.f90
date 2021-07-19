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

!> Definition of the isotropic third-order electrostatic contributions
module tblite_param_thirdorder
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_param_serde, only : serde_record
   use tblite_toml, only : toml_table, get_value, set_value, add_table
   implicit none
   private

   public :: thirdorder_record


   character(len=*), parameter :: k_shell = "shell", k_ang(0:4) = ["s", "p", "d", "f", "g"]

   !> Parametrization record for third-order electrostatic contributions
   type, extends(serde_record) :: thirdorder_record
      integer :: lmax
      logical :: shell
      real(wp) :: ksh(0:4)
   contains
      !> Read parametrization data from TOML data structure
      procedure :: load_from_toml
      !> Write parametrization data to TOML data structure
      procedure :: dump_to_toml
   end type


contains


!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   !> Instance of the parametrization data
   class(thirdorder_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child
   real(wp), allocatable :: last
   integer :: l, stat

   call get_value(table, k_shell, child, requested=.false.)
   self%shell = associated(child)
   if (self%shell) then
      do l = 0, 4
         if (.not.child%has_key(k_ang(l))) then
            if (allocated(last)) then
               self%ksh(l) = last
               cycle
            end if
            call fatal_error(error, "No entry for "//k_ang(l)//"-shell Hubbard derivative")
            exit
         end if
         call get_value(child, k_ang(l), self%ksh(l), stat=stat)
         if (stat /= 0) then
            call fatal_error(error, "Cannot read "//k_ang(l)//"-shell Hubbard derivative")
            exit
         end if
         if (stat == 0) then
            last = self%ksh(l)
            self%lmax = l
         end if
      end do
      if (allocated(error)) return
   end if
end subroutine load_from_toml


!> Write parametrization data to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(thirdorder_record), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child
   integer :: l

   if (self%shell) then
      call add_table(table, k_shell, child)
      do l = 0, self%lmax
         call set_value(child, k_ang(l), self%ksh(l))
      end do
   else
      call set_value(table, k_shell, .false.)
   end if
end subroutine dump_to_toml


end module tblite_param_thirdorder
