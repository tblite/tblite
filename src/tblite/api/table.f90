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

!> @file tblite/api/table.f90
!> Provides API exports for the #tblite_table handle.

!> API export for managing data tables
module tblite_api_table
   use, intrinsic :: iso_c_binding
   use mctc_env, only : fatal_error
   use tblite_api_error, only : vp_error
   use tblite_api_utils, only : c_f_character, f_c_character, strlen
   use tblite_api_version, only : namespace
   use tblite_toml, only : toml_table, toml_array, toml_key, toml_value, add_array, set_value, get_value, len, toml_error, toml_dump
   implicit none
   private

   public :: vp_table
   public :: vp_array
   public :: new_table_api, delete_table_api
   public :: table_set_double_api, table_set_int64_t_api, table_set_bool_api, &
      & table_set_char_api, table_add_table_api, dump_table_api
   public :: new_array_api, delete_array_api
   public :: array_push_back_double_api, array_push_back_int64_t_api, &
      & array_push_back_bool_api, array_push_back_char_api
   public :: array_size_api, array_get_type_api, array_get_double_api, &
      & array_get_int64_t_api, array_get_bool_api, array_get_char_api
   public :: table_set_array_api, table_get_type_api, table_get_bool_api, &
      & table_get_int64_t_api, table_get_double_api, table_get_char_api, &
      & table_get_table_api, table_get_array_api, table_get_n_keys_api, &
      & table_get_key_api

   !> Void pointer to manage general data tables
   type :: vp_table
      !> Actual payload
      type(toml_table), pointer :: ptr => null()
      !> Data is owned by the object
      logical :: owned
   end type vp_table

   !> Void pointer to manage general data arrays
   type :: vp_array
      !> Actual payload
      type(toml_array), pointer :: ptr => null()
      !> Data is owned by the object
      logical :: owned
   end type vp_array

   integer, parameter :: value_type_none = 0
   integer, parameter :: value_type_bool = 1
   integer, parameter :: value_type_int = 2
   integer, parameter :: value_type_double = 3
   integer, parameter :: value_type_char = 4
   integer, parameter :: value_type_array = 5
   integer, parameter :: value_type_table = 6

   logical, parameter :: debug = .false.

contains

!> Create data table reference object
function new_table_api(vtable) &
      & result(vval) &
      & bind(C, name=namespace//"new_table")
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   type(vp_table), pointer :: val
   type(c_ptr) :: vval
   type(toml_table), pointer :: dat

   if (debug) print '("[Info]", 1x, a)', "new_table"

   allocate(val)
   if (c_associated(vtable)) then
      call c_f_pointer(vtable, table)
      val%ptr => table%ptr
      val%owned = .false.
   else
      allocate(dat)
      dat = toml_table()
      val%ptr => dat
      val%owned = .true.
   end if
   vval = c_loc(val)
end function new_table_api

!> Delete data table object
subroutine delete_table_api(vtable) &
      & bind(C, name=namespace//"delete_table")
   type(c_ptr), intent(inout) :: vtable
   type(vp_table), pointer :: table

   if (debug) print '("[Info]", 1x, a)', "delete_table"

   if (c_associated(vtable)) then
      call c_f_pointer(vtable, table)

      if (table%owned) deallocate(table%ptr)
      deallocate(table)
      vtable = c_null_ptr
   end if

end subroutine delete_table_api

subroutine table_set_double_api(verror, vtable, ckey, val, n) &
      & bind(C, name=namespace//"table_set_double")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   character(kind=c_char), intent(in) :: ckey(*)
   character(len=:), allocatable :: key
   real(c_double), intent(in) :: val(*)
   integer(c_int), value :: n
   type(toml_array), pointer :: array
   integer :: i, stat

   if (debug) print '("[Info]", 1x, a)', "table_set_double"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if

   call c_f_pointer(vtable, table)
   call c_f_character(ckey, key)

   if (table%ptr%has_key(key)) call table%ptr%delete(key)

   if (n > 0) then
      call add_array(table%ptr, key, array)
      do i = 1, n
         call set_value(array, i, val(i), stat=stat)
         if (stat /= 0) exit
      end do
   else
      call set_value(table%ptr, key, val(1), stat=stat)
   end if

   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to push back double value(s) to data table")
   end if
end subroutine table_set_double_api

subroutine table_set_int64_t_api(verror, vtable, ckey, val, n) &
      & bind(C, name=namespace//"table_set_int64_t")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   character(kind=c_char), intent(in) :: ckey(*)
   character(len=:), allocatable :: key
   integer(c_int64_t), intent(in) :: val(*)
   integer(c_int), value :: n
   type(toml_array), pointer :: array
   integer :: i, stat

   if (debug) print '("[Info]", 1x, a)', "table_set_int64_t"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if

   call c_f_pointer(vtable, table)
   call c_f_character(ckey, key)

   if (table%ptr%has_key(key)) call table%ptr%delete(key)

   if (n > 0) then
      call add_array(table%ptr, key, array)
      do i = 1, n
         call set_value(array, i, val(i), stat=stat)
         if (stat /= 0) exit
      end do
   else
      call set_value(table%ptr, key, val(1), stat=stat)
   end if

   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to push back integer value(s) to data table")
   end if
end subroutine table_set_int64_t_api

subroutine table_set_bool_api(verror, vtable, ckey, val, n) &
      & bind(C, name=namespace//"table_set_bool")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   character(kind=c_char), intent(in) :: ckey(*)
   character(len=:), allocatable :: key
   logical(c_bool), intent(in) :: val(*)
   integer(c_int), value :: n
   type(toml_array), pointer :: array
   integer :: i, stat

   if (debug) print '("[Info]", 1x, a)', "table_set_bool"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if

   call c_f_pointer(vtable, table)
   call c_f_character(ckey, key)

   if (table%ptr%has_key(key)) call table%ptr%delete(key)

   if (n > 0) then
      call add_array(table%ptr, key, array)
      do i = 1, n
         call set_value(array, i, merge(.true., .false., val(i)), stat=stat)
         if (stat /= 0) exit
      end do
   else
      call set_value(table%ptr, key, merge(.true., .false., val(1)), stat=stat)
   end if

   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to push back boolean value(s) to data table")
   end if
end subroutine table_set_bool_api

subroutine table_set_char_api(verror, vtable, ckey, cval, n) &
      & bind(C, name=namespace//"table_set_char")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   character(kind=c_char), intent(in) :: ckey(*)
   character(len=:), allocatable :: key, val
   type(c_ptr), value :: cval
   character(kind=c_char), pointer :: carr(:, :)
   integer(c_int), value :: n
   type(toml_array), pointer :: array
   integer :: i, stat

   if (debug) print '("[Info]", 1x, a)', "table_set_char"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if

   call c_f_pointer(vtable, table)
   call c_f_character(ckey, key)
   call c_f_pointer(cval, carr, [strlen(cval)+1, max(n, 1)])

   if (table%ptr%has_key(key)) call table%ptr%delete(key)

   if (n > 0) then
      call add_array(table%ptr, key, array)
      do i = 1, n
         call c_f_character(carr(:, i), val)
         call set_value(array, i, val, stat=stat)
         if (stat /= 0) exit
      end do
   else
      call c_f_character(carr(:, 1), val)
      call set_value(table%ptr, key, val, stat=stat)
   end if

   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to push back character string to data table")
   end if
end subroutine table_set_char_api

function table_add_table_api(verror, vtable, ckey) &
      & result(vval) &
      & bind(C, name=namespace//"table_add_table")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   character(kind=c_char), intent(in) :: ckey(*)
   character(len=:), allocatable :: key
   type(c_ptr) :: vval
   type(vp_table), pointer :: val
   type(toml_table), pointer :: tmp
   integer :: stat

   if (debug) print '("[Info]", 1x, a)', "table_add_table"

   vval = c_null_ptr
   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if

   vval = new_table_api(vtable)
   call c_f_pointer(vtable, table)
   call c_f_pointer(vval, val)
   call c_f_character(ckey, key)

   call get_value(table%ptr, key, tmp, stat=stat)
   val%ptr => tmp

   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to push back subtable to data table")
      call delete_table_api(vval)
   end if
end function table_add_table_api

function new_array_api() &
      & result(vval) &
      & bind(C, name=namespace//"new_array")
   type(c_ptr) :: vval
   type(vp_array), pointer :: val
   type(toml_array), pointer :: dat

   if (debug) print '("[Info]", 1x, a)', "new_array"

   allocate(val)
   allocate(dat)
   dat = toml_array()
   val%ptr => dat
   val%owned = .true.
   vval = c_loc(val)
end function new_array_api

subroutine delete_array_api(varray) &
      & bind(C, name=namespace//"delete_array")
   type(c_ptr), intent(inout) :: varray
   type(vp_array), pointer :: array

   if (debug) print '("[Info]", 1x, a)', "delete_array"

   if (c_associated(varray)) then
      call c_f_pointer(varray, array)

      if (array%owned) deallocate(array%ptr)
      deallocate(array)
      varray = c_null_ptr
   end if
end subroutine delete_array_api

subroutine array_push_back_double_api(verror, varray, value) &
      & bind(C, name=namespace//"array_push_back_double")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: varray
   type(vp_array), pointer :: array
   real(c_double), value :: value
   integer :: stat

   if (debug) print '("[Info]", 1x, a)', "array_push_back_double"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(varray)) then
      call fatal_error(error%ptr, "Data array object is missing")
      return
   end if

   call c_f_pointer(varray, array)
   call set_value(array%ptr, len(array%ptr) + 1, value, stat=stat)

   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to push back double value to data array")
   end if
end subroutine array_push_back_double_api

subroutine array_push_back_int64_t_api(verror, varray, value) &
      & bind(C, name=namespace//"array_push_back_int64_t")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: varray
   type(vp_array), pointer :: array
   integer(c_int64_t), value :: value
   integer :: stat

   if (debug) print '("[Info]", 1x, a)', "array_push_back_int64_t"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(varray)) then
      call fatal_error(error%ptr, "Data array object is missing")
      return
   end if

   call c_f_pointer(varray, array)
   call set_value(array%ptr, len(array%ptr) + 1, value, stat=stat)

   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to push back integer value to data array")
   end if
end subroutine array_push_back_int64_t_api

subroutine array_push_back_bool_api(verror, varray, value) &
      & bind(C, name=namespace//"array_push_back_bool")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: varray
   type(vp_array), pointer :: array
   logical(c_bool), value :: value
   integer :: stat

   if (debug) print '("[Info]", 1x, a)', "array_push_back_bool"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(varray)) then
      call fatal_error(error%ptr, "Data array object is missing")
      return
   end if

   call c_f_pointer(varray, array)
   call set_value(array%ptr, len(array%ptr) + 1, merge(.true., .false., value), stat=stat)

   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to push back boolean value to data array")
   end if
end subroutine array_push_back_bool_api

subroutine array_push_back_char_api(verror, varray, cvalue) &
      & bind(C, name=namespace//"array_push_back_char")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: varray
   type(vp_array), pointer :: array
   character(kind=c_char), intent(in) :: cvalue(*)
   character(len=:), allocatable :: value
   integer :: stat

   if (debug) print '("[Info]", 1x, a)', "array_push_back_char"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(varray)) then
      call fatal_error(error%ptr, "Data array object is missing")
      return
   end if

   call c_f_pointer(varray, array)
   call c_f_character(cvalue, value)
   call set_value(array%ptr, len(array%ptr) + 1, value, stat=stat)

   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to push back character value to data array")
   end if
end subroutine array_push_back_char_api

function array_size_api(verror, varray) result(n) &
      & bind(C, name=namespace//"array_size")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: varray
   type(vp_array), pointer :: array
   integer(c_int) :: n

   if (debug) print '("[Info]", 1x, a)', "array_size"

   n = 0
   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(varray)) then
      call fatal_error(error%ptr, "Data array object is missing")
      return
   end if

   call c_f_pointer(varray, array)
   n = len(array%ptr)
end function array_size_api

function array_get_type_api(verror, varray, index) result(itype) &
      & bind(C, name=namespace//"array_get_type")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: varray
   type(vp_array), pointer :: array
   integer(c_int), value :: index
   integer(c_int) :: itype
   logical :: tmp_bool
   integer(c_int64_t) :: tmp_int
   real(c_double) :: tmp_double
   character(len=:), allocatable :: tmp_char
   integer :: stat

   if (debug) print '("[Info]", 1x, a)', "array_get_type"

   itype = value_type_none
   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(varray)) then
      call fatal_error(error%ptr, "Data array object is missing")
      return
   end if

   call c_f_pointer(varray, array)
   if (index < 1 .or. index > len(array%ptr)) return

   call get_value(array%ptr, index, tmp_bool, stat=stat)
   if (stat == 0) then
      itype = value_type_bool
      return
   end if

   call get_value(array%ptr, index, tmp_int, stat=stat)
   if (stat == 0) then
      itype = value_type_int
      return
   end if

   call get_value(array%ptr, index, tmp_double, stat=stat)
   if (stat == 0) then
      itype = value_type_double
      return
   end if

   call get_value(array%ptr, index, tmp_char, stat=stat)
   if (stat == 0) then
      itype = value_type_char
   end if
end function array_get_type_api

subroutine array_get_double_api(verror, varray, index, value) &
      & bind(C, name=namespace//"array_get_double")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: varray
   type(vp_array), pointer :: array
   integer(c_int), value :: index
   real(c_double) :: value
   integer :: stat

   if (debug) print '("[Info]", 1x, a)', "array_get_double"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(varray)) then
      call fatal_error(error%ptr, "Data array object is missing")
      return
   end if

   call c_f_pointer(varray, array)
   call get_value(array%ptr, index, value, stat=stat)
   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to read double value from data array")
   end if
end subroutine array_get_double_api

subroutine array_get_int64_t_api(verror, varray, index, value) &
      & bind(C, name=namespace//"array_get_int64_t")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: varray
   type(vp_array), pointer :: array
   integer(c_int), value :: index
   integer(c_int64_t) :: value
   integer :: stat

   if (debug) print '("[Info]", 1x, a)', "array_get_int64_t"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(varray)) then
      call fatal_error(error%ptr, "Data array object is missing")
      return
   end if

   call c_f_pointer(varray, array)
   call get_value(array%ptr, index, value, stat=stat)
   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to read integer value from data array")
   end if
end subroutine array_get_int64_t_api

subroutine array_get_bool_api(verror, varray, index, value) &
      & bind(C, name=namespace//"array_get_bool")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: varray
   type(vp_array), pointer :: array
   integer(c_int), value :: index
   logical(c_bool) :: value
   logical :: tmp
   integer :: stat

   if (debug) print '("[Info]", 1x, a)', "array_get_bool"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(varray)) then
      call fatal_error(error%ptr, "Data array object is missing")
      return
   end if

   call c_f_pointer(varray, array)
   call get_value(array%ptr, index, tmp, stat=stat)
   value = tmp
   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to read boolean value from data array")
   end if
end subroutine array_get_bool_api

subroutine array_get_char_api(verror, varray, index, cvalue, n) &
      & bind(C, name=namespace//"array_get_char")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: varray
   type(vp_array), pointer :: array
   integer(c_int), value :: index
   character(kind=c_char), intent(out) :: cvalue(*)
   integer(c_int), value :: n
   character(len=:), allocatable :: value
   integer :: stat

   if (debug) print '("[Info]", 1x, a)', "array_get_char"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(varray)) then
      call fatal_error(error%ptr, "Data array object is missing")
      return
   end if

   call c_f_pointer(varray, array)
   call get_value(array%ptr, index, value, stat=stat)
   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to read character value from data array")
   end if
   call f_c_character(value, cvalue, n)
end subroutine array_get_char_api

subroutine table_set_array_api(verror, vtable, ckey, varray) &
      & bind(C, name=namespace//"table_set_array")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   character(kind=c_char), intent(in) :: ckey(*)
   type(c_ptr), value :: varray
   type(vp_array), pointer :: array
   character(len=:), allocatable :: key
   type(toml_array), pointer :: values
   integer(c_int) :: n, i, itype, stat
   logical :: tmp_bool
   integer(c_int64_t) :: tmp_int
   real(c_double) :: tmp_double
   character(len=:), allocatable :: tmp_char

   if (debug) print '("[Info]", 1x, a)', "table_set_array"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if

   call c_f_pointer(vtable, table)
   call c_f_character(ckey, key)

   if (table%ptr%has_key(key)) call table%ptr%delete(key)

   call c_f_pointer(varray, array)
   call add_array(table%ptr, key, values)
   n = len(array%ptr)

   do i = 1, n
      call get_value(array%ptr, i, tmp_bool, stat=stat)
      if (stat == 0) then
         call set_value(values, i, tmp_bool, stat=stat)
         if (stat /= 0) exit
         cycle
      end if

      call get_value(array%ptr, i, tmp_int, stat=stat)
      if (stat == 0) then
         call set_value(values, i, tmp_int, stat=stat)
         if (stat /= 0) exit
         cycle
      end if

      call get_value(array%ptr, i, tmp_double, stat=stat)
      if (stat == 0) then
         call set_value(values, i, tmp_double, stat=stat)
         if (stat /= 0) exit
         cycle
      end if

      call get_value(array%ptr, i, tmp_char, stat=stat)
      if (stat == 0) then
         call set_value(values, i, tmp_char, stat=stat)
         if (stat /= 0) exit
      end if
   end do

   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to push back array value(s) to data table")
   end if
end subroutine table_set_array_api

function table_get_type_api(verror, vtable, ckey) result(itype) &
      & bind(C, name=namespace//"table_get_type")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   character(kind=c_char), intent(in) :: ckey(*)
   character(len=:), allocatable :: key
   integer(c_int) :: itype
   logical :: tmp_bool
   integer(c_int64_t) :: tmp_int
   real(c_double) :: tmp_double
   character(len=:), allocatable :: tmp_char
   type(toml_table), pointer :: tmp_table
   type(toml_array), pointer :: tmp_array
   integer :: stat

   if (debug) print '("[Info]", 1x, a)', "table_get_type"

   itype = value_type_none
   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if

   call c_f_pointer(vtable, table)
   call c_f_character(ckey, key)

   call get_value(table%ptr, key, tmp_table, stat=stat)
   if (stat == 0) then
      itype = value_type_table
      return
   end if

   call get_value(table%ptr, key, tmp_array, stat=stat)
   if (stat == 0) then
      itype = value_type_array
      return
   end if

   call get_value(table%ptr, key, tmp_bool, stat=stat)
   if (stat == 0) then
      itype = value_type_bool
      return
   end if

   call get_value(table%ptr, key, tmp_int, stat=stat)
   if (stat == 0) then
      itype = value_type_int
      return
   end if

   call get_value(table%ptr, key, tmp_double, stat=stat)
   if (stat == 0) then
      itype = value_type_double
      return
   end if

   call get_value(table%ptr, key, tmp_char, stat=stat)
   if (stat == 0) then
      itype = value_type_char
   end if
end function table_get_type_api

subroutine table_get_bool_api(verror, vtable, ckey, value) &
      & bind(C, name=namespace//"table_get_bool")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   character(kind=c_char), intent(in) :: ckey(*)
   character(len=:), allocatable :: key
   logical(c_bool) :: value
   logical :: tmp
   integer :: stat

   if (debug) print '("[Info]", 1x, a)', "table_get_bool"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if

   call c_f_pointer(vtable, table)
   call c_f_character(ckey, key)
   call get_value(table%ptr, key, tmp, stat=stat)
   value = tmp
   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to read boolean value from data table")
   end if
end subroutine table_get_bool_api

subroutine table_get_int64_t_api(verror, vtable, ckey, value) &
      & bind(C, name=namespace//"table_get_int64_t")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   character(kind=c_char), intent(in) :: ckey(*)
   character(len=:), allocatable :: key
   integer(c_int64_t) :: value
   integer :: stat

   if (debug) print '("[Info]", 1x, a)', "table_get_int64_t"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if

   call c_f_pointer(vtable, table)
   call c_f_character(ckey, key)
   call get_value(table%ptr, key, value, stat=stat)
   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to read integer value from data table")
   end if
end subroutine table_get_int64_t_api

subroutine table_get_double_api(verror, vtable, ckey, value) &
      & bind(C, name=namespace//"table_get_double")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   character(kind=c_char), intent(in) :: ckey(*)
   character(len=:), allocatable :: key
   real(c_double) :: value
   integer :: stat

   if (debug) print '("[Info]", 1x, a)', "table_get_double"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if

   call c_f_pointer(vtable, table)
   call c_f_character(ckey, key)
   call get_value(table%ptr, key, value, stat=stat)
   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to read double value from data table")
   end if
end subroutine table_get_double_api

subroutine table_get_char_api(verror, vtable, ckey, cvalue, n) &
      & bind(C, name=namespace//"table_get_char")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   character(kind=c_char), intent(in) :: ckey(*)
   character(len=:), allocatable :: key
   character(kind=c_char), intent(out) :: cvalue(*)
   integer(c_int), value :: n
   character(len=:), allocatable :: value
   integer :: stat

   if (debug) print '("[Info]", 1x, a)', "table_get_char"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if

   call c_f_pointer(vtable, table)
   call c_f_character(ckey, key)
   call get_value(table%ptr, key, value, stat=stat)
   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to read character value from data table")
   end if
   call f_c_character(value, cvalue, n)
end subroutine table_get_char_api

function table_get_table_api(verror, vtable, ckey) result(vval) &
      & bind(C, name=namespace//"table_get_table")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   character(kind=c_char), intent(in) :: ckey(*)
   character(len=:), allocatable :: key
   type(c_ptr) :: vval
   type(vp_table), pointer :: val
   type(toml_table), pointer :: tmp
   integer :: stat

   if (debug) print '("[Info]", 1x, a)', "table_get_table"

   vval = c_null_ptr
   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if

   call c_f_pointer(vtable, table)
   call c_f_character(ckey, key)
   call get_value(table%ptr, key, tmp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to read subtable from data table")
      return
   end if

   allocate(val)
   val%ptr => tmp
   val%owned = .false.
   vval = c_loc(val)
end function table_get_table_api

function table_get_array_api(verror, vtable, ckey) result(vval) &
      & bind(C, name=namespace//"table_get_array")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   character(kind=c_char), intent(in) :: ckey(*)
   character(len=:), allocatable :: key
   type(c_ptr) :: vval
   type(vp_array), pointer :: val
   type(toml_array), pointer :: tmp
   integer :: stat

   if (debug) print '("[Info]", 1x, a)', "table_get_array"

   vval = c_null_ptr
   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if

   call c_f_pointer(vtable, table)
   call c_f_character(ckey, key)
   call get_value(table%ptr, key, tmp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to read array from data table")
      return
   end if

   allocate(val)
   val%ptr => tmp
   val%owned = .false.
   vval = c_loc(val)
end function table_get_array_api

function table_get_n_keys_api(verror, vtable) result(n) &
      & bind(C, name=namespace//"table_get_n_keys")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   integer(c_int) :: n
   type(toml_key), allocatable :: list(:)

   if (debug) print '("[Info]", 1x, a)', "table_get_n_keys"

   n = 0
   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if

   call c_f_pointer(vtable, table)
   call table%ptr%get_keys(list)
   n = size(list)
end function table_get_n_keys_api

subroutine table_get_key_api(verror, vtable, index, ckey, n) &
      & bind(C, name=namespace//"table_get_key")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   integer(c_int), value :: index
   character(kind=c_char), intent(out) :: ckey(*)
   integer(c_int), value :: n
   type(toml_key), allocatable :: list(:)
   character(len=:), allocatable :: key

   if (debug) print '("[Info]", 1x, a)', "table_get_key"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if

   call c_f_pointer(vtable, table)
   call table%ptr%get_keys(list)
   if (index < 1 .or. index > size(list)) then
      call fatal_error(error%ptr, "Table key index out of range")
      return
   end if

   key = list(index)%key
   call f_c_character(key, ckey, n)
end subroutine table_get_key_api

subroutine dump_table_api(verror, vtable, cfilename) &
      & bind(C, name=namespace//"dump_table")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   character(kind=c_char), intent(in) :: cfilename(*)
   character(len=:), allocatable :: filename
   type(toml_error), allocatable :: ser_error
   integer :: unit

   if (debug) print '("[Info]", 1x, a)', "dump_table"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if
   call c_f_pointer(vtable, table)
   call c_f_character(cfilename, filename)

   open(newunit=unit, file=filename, status="replace", action="write")
   call toml_dump(table%ptr, unit, ser_error)
   close(unit)

   if (allocated(ser_error)) then
      call fatal_error(error%ptr, ser_error%message)
   end if
end subroutine dump_table_api


end module tblite_api_table
