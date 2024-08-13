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

!> @dir tblite/double_dictionary
!> Contains implemenation of a dictionary with strings as keys and arrays of doubles as entries.

!> @file tblite/double_dictionery.f90
!> Implements double dictionary type
module tblite_double_dictionary
   use mctc_env_accuracy, only : wp, i8
   use mctc_env, only : error_type, fatal_error
   use tblite_toml, only : toml_array, toml_table, toml_key, add_table, set_value, toml_error
   use tblite_toml, only : toml_dump, add_array,  get_value, toml_parse
   implicit none
   private

   public :: double_dictionary_type

   type :: double_record
      character(len=:), allocatable :: label
      real(wp), allocatable :: array1(:)
      real(wp), allocatable :: array2(:, :)
      real(wp), allocatable :: array3(:, :, :)
   contains
      generic :: assignment(=) => copy_record
      generic :: operator(==) => equal_record
      procedure :: copy_record
      procedure :: equal_record
   end type double_record

   type :: double_dictionary_type
      integer :: n = 0
      type(double_record), allocatable :: record(:)
   contains
      private
      generic, public :: initialize_entry => ini_label, ini_1d, ini_2d, ini_3d
      procedure :: ini_label
      procedure :: ini_1d
      procedure :: ini_2d
      procedure :: ini_3d
      generic, public :: add_entry =>  add_1d, add_2d, add_3d
      procedure :: add_1d
      procedure :: add_2d
      procedure :: add_3d
      !procedure :: update_entry ! label and array pair
      generic, public :: get_entry =>  get_1d_index, get_2d_index, get_3d_index, get_1d_label, get_2d_label, get_3d_label !check
      procedure :: get_1d_index
      procedure :: get_2d_index
      procedure :: get_3d_index
      procedure :: get_1d_label
      procedure :: get_2d_label
      procedure :: get_3d_label
      generic, public :: update => update_1d, update_2d, update_3d
      procedure :: update_1d
      procedure :: update_2d
      procedure :: update_3d
      procedure, public :: get_label ! return label label
      procedure :: push
      procedure, public :: get_n_entries
      generic, public :: concatenate => concatenate_overwrite
      procedure :: concatenate_overwrite
      generic, public :: assignment(=) => copy
      procedure :: copy
      generic, public :: operator(+) => combine_dict
      procedure :: combine_dict
      generic, public :: remove_entry => remove_entry_label, remove_entry_index
      procedure :: remove_entry_label
      procedure :: remove_entry_index
      generic, public :: get_index => return_label_index
      procedure :: return_label_index
      procedure :: dump_to_toml
      procedure :: dump_to_file
      procedure :: dump_to_unit
      generic, public :: dump => dump_to_file, dump_to_toml, dump_to_unit
      generic, public :: operator(==) => equal_dict
      procedure :: equal_dict
      generic, public :: load => load_from_file, load_from_unit, load_from_toml
      procedure :: load_from_toml
      procedure :: load_from_file
      procedure :: load_from_unit

   end type double_dictionary_type

contains

function equal_record(lhs, rhs) result(equal)
   class(double_record), intent(in) :: lhs, rhs
   integer :: i
   logical :: equal
   equal = .false.

   if (lhs%label /= rhs%label) then
      return
   end if

   if (allocated(lhs%array1) .and. allocated(rhs%array1)) then
      if (all(lhs%array1 == rhs%array1)) then
         equal = .true.
         return
      else
         return
      end if 
   endif
      
   if (allocated(lhs%array2) .and. allocated(rhs%array2)) then
      if (all(lhs%array2 == rhs%array2)) then
         equal = .true.
         return
      else
         return
      end if 
   endif

   if (allocated(lhs%array3) .and. allocated(rhs%array3)) then
      if (all(lhs%array3 == rhs%array3)) then
         equal = .true.
         return
      else
         return
      end if 
   endif

end function


function equal_dict(lhs, rhs) result(equal)
   class(double_dictionary_type), intent(in) :: lhs, rhs
   integer :: i
   logical :: equal
   equal = .false.

   if (lhs%get_n_entries() /= rhs%get_n_entries()) then
      return
   end if

   do i = 1, lhs%get_n_entries()
      if (.not.(lhs%record(i) == rhs%record(i))) return
   end do

   equal = .true.
end function

!> Read double dictionary data from file
subroutine load_from_file(self, file, error)
   !> Instance of the parametrization data
   class(double_dictionary_type), intent(inout) :: self
   !> File name
   character(len=*), intent(in) :: file
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: unit
   logical :: exist

   inquire(file=file, exist=exist)
   if (.not.exist) then
     call fatal_error(error, "Could not find toml file '"//file//"'")
     return
   end if

   open(file=file, newunit=unit)
   call self%load(unit, error)
   close(unit)
end subroutine load_from_file


!> Read double_dictionary data from file
subroutine load_from_unit(self, unit, error)
   !> Instance of the double dictionary data
   class(double_dictionary_type), intent(inout) :: self
   !> File name
   integer, intent(in) :: unit
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_error), allocatable :: parse_error
   type(toml_table), allocatable :: table

   call toml_parse(table, unit, parse_error)

   if (allocated(parse_error)) then
      allocate(error)
      call move_alloc(parse_error%message, error%message)
      return
   end if

   call self%load(table, error)
   if (allocated(error)) return

end subroutine load_from_unit

subroutine load_from_toml(self, table, error)
   use tblite_toml, only : len
   !iterate over entries and dump to toml
   class(double_dictionary_type) :: self
   !> toml table to add entries to
   type(toml_table), intent(inout) :: table
   type(toml_key), allocatable :: list_keys(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(toml_array), pointer :: array
   real(kind=wp), allocatable :: array1(:)

   integer :: i, stat

   call table%get_keys(list_keys)

   do i = 1, size(list_keys)
      call get_value(table, list_keys(i), array, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot read entry for array")
         return
      end if
      if (allocated(array1)) deallocate(array1)
      allocate(array1(len(array)))
      call get_value(array, array1)
      call self%add_entry(list_keys(i)%key, array1)
   end do
end subroutine

subroutine dump_to_toml(self, table, error)
   !iterate over entries and dump to toml
   class(double_dictionary_type) :: self
   !> toml table to add entries to
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(toml_array), pointer :: array
   real(kind=wp), allocatable :: array1(:), array2(:, :), array3(:, :, :)

   integer :: i, stat, j
   do i = 1, self%get_n_entries()
      call add_array(table, self%record(i)%label, array)

      if (allocated(array1)) deallocate(array1)
      call self%get_entry(i, array1)
      if (allocated(array1)) then
         do j = 1, size(array1)
            call set_value(array, j, array1(j), stat=stat)
         end do
         cycle
      end if

      call self%get_entry(i, array2)
      if (allocated(array2)) then
         array1 = reshape(array2, [size(array2, 1)*size(array2, 2)])
         deallocate(array2)
         do j = 1, size(array1)
            call set_value(array, j, array1(j), stat=stat)
         end do
         cycle
      end if

      call self%get_entry(i, array3)
      if (allocated(array3)) then
         array1 = reshape(array3, [size(array3, 1)*size(array3, 2)*size(array3, 3)])
         do j = 1, size(array1)
            call set_value(array, j, array1(j), stat=stat)
         end do
         deallocate(array3)
         cycle
      end if
   end do

   deallocate(array1)
end subroutine

subroutine dump_to_file(self, file, error)
   !> Instance of the parametrization data
   class(double_dictionary_type), intent(in) :: self
   !> File name
   character(len=*), intent(in) :: file
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: unit

   open(file=file, newunit=unit)
   call self%dump(unit, error)
   close(unit)
   if (allocated(error)) return

end subroutine dump_to_file


!> Write double dictionary data to file
subroutine dump_to_unit(self, unit, error)
   !> Instance of the parametrization data
   class(double_dictionary_type), intent(in) :: self
   !> Formatted unit
   integer, intent(in) :: unit
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(toml_error), allocatable :: ser_error

   table = toml_table()

   call self%dump(table, error)

   call toml_dump(table, unit, ser_error)
   if (allocated(ser_error)) then
      call fatal_error(error, ser_error%message)
   end if

end subroutine dump_to_unit

subroutine remove_entry_index(self, index)
   class(double_dictionary_type) :: self
   integer :: index, old_n, i, it
   type(double_dictionary_type) :: tmp

   if (index > self%n) return
   tmp = self
   old_n = self%n
   self%n = self%n - 1
   
   deallocate(self%record)
   allocate(self%record(self%n))
   it = 1
   do i = 1, old_n
      if (i == index) cycle
      self%record(it) = tmp%record(i)
      it = it + 1
   end do
 
end subroutine

subroutine remove_entry_label(self, label)
   class(double_dictionary_type) :: self
   character(len=*) :: label
   type(double_dictionary_type) :: tmp
   integer :: it 
   it = return_label_index(self, label)
   if (it /= 0) then
      call self%remove_entry_index(it) 
   else
      return
   end if
end subroutine 

function get_n_entries(self) result(n)
   class(double_dictionary_type) :: self
   integer :: n
   n = self%n
end function

subroutine copy(to, from)
   class(double_dictionary_type), intent(inout) :: to
   type(double_dictionary_type), intent(in) :: from
   integer :: n_entries, i
   to%n = from%n
   if (allocated(to%record)) deallocate(to%record)
   if (from%get_n_entries() > 0) then 
      allocate(to%record(size(from%record)))
      n_entries = from%get_n_entries()

      do i = 1, n_entries
         to%record(i) = from%record(i)
      end do
   end if

end subroutine

subroutine copy_record(to, from)
   class(double_record), intent(inout) :: to
   type(double_record), intent(in) :: from
   integer :: n_entries, it, i
   if (allocated(to%label)) deallocate(to%label)
   to%label = from%label
   if (allocated(from%array1)) to%array1 = from%array1
   if (allocated(from%array2)) to%array2 = from%array2
   if (allocated(from%array3)) to%array3 = from%array3
end subroutine

subroutine concatenate_overwrite(self, dict2)
   class(double_dictionary_type), intent(inout) :: self
   type(double_dictionary_type), intent(in) :: dict2
   type(double_dictionary_type) :: tmp_dict
   tmp_dict = self + dict2
   self = tmp_dict

end subroutine

function combine_dict(self, dict2) result(new_dict)
   class(double_dictionary_type), intent(in) :: self
   type(double_dictionary_type), intent(in) :: dict2
   type(double_dictionary_type) :: new_dict
   integer :: it, i, n_entries
   new_dict = self
   associate(dict => dict2)
      n_entries = dict%get_n_entries()
      do i = 1, n_entries
         call new_dict%push(dict%record(i)%label, it)
         new_dict%record(it) = dict%record(i)
      end do
   end associate

end function

subroutine update_1d(self, label, array)
   class(double_dictionary_type) :: self
   character(len=*) :: label
   real(wp) :: array(:)
   integer :: it

   it = return_label_index(self, label)
   if (it /= 0) then
      associate(record => self%record(it))
         if (allocated(record%array1)) deallocate(record%array1)
         if (allocated(record%array2)) deallocate(record%array2)
         if (allocated(record%array3)) deallocate(record%array3)
         record%array1 = array
      end associate
   else
      return
   end if
end subroutine

subroutine update_2d(self, label, array)
   class(double_dictionary_type) :: self
   character(len=*) :: label
   real(wp) :: array(:, :)
   integer :: it

   it = return_label_index(self, label)
   if (it /= 0) then
      associate(record => self%record(it))
         if (allocated(record%array1)) deallocate(record%array1)
         if (allocated(record%array2)) deallocate(record%array2)
         if (allocated(record%array3)) deallocate(record%array3)
         record%array2 = array
      end associate
   else
      return
   end if
end subroutine

subroutine update_3d(self, label, array)
   class(double_dictionary_type) :: self
   character(len=*) :: label
   real(wp) :: array(:, :, :)
   integer :: it

   it = return_label_index(self, label)
   if (it /= 0) then
      associate(record => self%record(it))
         if (allocated(record%array1)) deallocate(record%array1)
         if (allocated(record%array2)) deallocate(record%array2)
         if (allocated(record%array3)) deallocate(record%array3)
         record%array3 = array
      end associate
   else
      return
   end if
end subroutine

subroutine get_label(self, index, label)
   class(double_dictionary_type) :: self
   integer :: index
   character(len=:), allocatable :: label

   if (index > self%n) return
   if (allocated(label)) deallocate(label)

   label = self%record(index)%label
end subroutine

subroutine get_1d_label(self, label, array)
   class(double_dictionary_type) :: self
   character(len=*) :: label
   real(wp), allocatable :: array(:)
   integer :: it
   it = return_label_index(self, label)
   if (it == 0) return
   call self%get_entry(it, array)

end subroutine

subroutine get_2d_label(self, label, array)
   class(double_dictionary_type) :: self
   character(len=*) :: label
   real(wp), allocatable :: array(:,:)
   integer :: it
   it = return_label_index(self, label)
   if (it == 0) return
   call self%get_entry(it, array)

end subroutine

subroutine get_3d_label(self, label, array)
   class(double_dictionary_type) :: self
   character(len=*) :: label
   real(wp), allocatable :: array(:, :, :)
   integer :: it
   it = return_label_index(self, label)
   if (it == 0) return
   call self%get_entry(it, array)

end subroutine

subroutine get_1d_index(self, index, array)
   class(double_dictionary_type) :: self
   integer :: index
   real(wp), allocatable :: array(:)

   if (index > self%n) return
   if (index <= 0) return
   associate(rec => self%record(index))
      if (allocated(rec%array1)) then
         if (allocated(array)) deallocate(array)
         allocate(array(size(rec%array1, dim = 1)))
         array = rec%array1
      else
         return
      end if
   end associate

end subroutine

subroutine get_2d_index(self, index, array)
   class(double_dictionary_type) :: self
   integer :: index
   real(wp), allocatable :: array(:,:)

   if (index > self%n) return
   if (index <= 0) return
   associate(rec => self%record(index))
      if (allocated(rec%array2)) then
         if (allocated(array)) deallocate(array)
         allocate(array(size(rec%array2, dim = 1), size(rec%array2, dim = 2)))
         array = rec%array2
      else
         return
      end if
   end associate

end subroutine

subroutine get_3d_index(self, index, array)
   class(double_dictionary_type) :: self
   integer :: index
   real(wp), allocatable :: array(:, :, :)

   if (index > self%n) return
   if (index <= 0) return
   associate(rec => self%record(index))
      if (allocated(rec%array3)) then
         if (allocated(array)) deallocate(array)
         allocate(array(size(rec%array3, dim = 1), size(rec%array3, dim = 2), size(rec%array3, dim=3)))
         array = rec%array3
      else
         return
      end if
   end associate

end subroutine

subroutine push(self, label, it)
   class(double_dictionary_type), intent(inout) :: self
   character(len=*), intent(in) :: label

   integer, intent(out) :: it

   if (.not.allocated(self%record)) call resize(self%record)
   it = find(self%record(:self%n), label)

   if (it == 0) then
      if (self%n >= size(self%record)) then
         call resize(self%record)
      end if

      self%n = self%n + 1
      it = self%n
      self%record(it) = double_record(label=label)
   end if
end subroutine push

subroutine ini_label(self,label)
   class(double_dictionary_type) :: self
   character(len=*) :: label
   integer :: it

   call self%push(label, it)
end subroutine

subroutine ini_1d(self, label, ndim1)
   class(double_dictionary_type) :: self
   character(len=*) :: label
   integer :: ndim1
   integer :: it

   call self%push(label, it)

   associate(record => self%record(it))
      allocate(record%array1(ndim1), source = 0.0_wp)
   end associate
end subroutine

subroutine ini_2d(self, label, ndim1, ndim2)

   class(double_dictionary_type) :: self
   character(len=*) :: label
   integer :: ndim1, ndim2

   integer :: it

   call self%push(label, it)

   associate(record => self%record(it))
      allocate(record%array2(ndim1, ndim2), source = 0.0_wp)
   end associate
end subroutine

subroutine ini_3d(self, label, ndim1, ndim2, ndim3)
   class(double_dictionary_type) :: self
   character(len=*) :: label
   integer :: ndim1, ndim2, ndim3
   integer :: it

   call self%push(label, it)

   associate(record => self%record(it))
      allocate(record%array3(ndim1, ndim2, ndim3), source = 0.0_wp)
   end associate
end subroutine

subroutine add_1d(self, label, array)
   class(double_dictionary_type) :: self
   character(len=*) :: label
   real(wp) :: array(:)
   integer :: it

   call self%push(label, it)

   associate(record => self%record(it))
      record%array1 = array
   end associate
end subroutine

subroutine add_2d(self, label, array)
   class(double_dictionary_type) :: self
   character(len=*) :: label
   real(wp) :: array(:,:)
   integer :: it

   call self%push(label, it)

   associate(record => self%record(it))
      record%array2 = array
   end associate

end subroutine

subroutine add_3d(self, label, array)
   class(double_dictionary_type) :: self
   character(len=*) :: label
   real(wp) :: array(:, :, :)
   integer :: it

   call self%push(label, it)

   associate(record => self%record(it))
      record%array3 = array
   end associate
end subroutine


function return_label_index(self, label) result(it)
   class(double_dictionary_type), intent(in) :: self
   character(len=*), intent(in) :: label
   integer :: it
   it = 0
   if (self%n <= 0) return
   it = find(self%record(:self%n), label)
   if (it == 0) return

end function return_label_index


pure function find(record, label) result(pos)
   type(double_record), intent(in) :: record(:)
   character(len=*), intent(in) :: label
   integer :: pos

   integer :: i

   pos = 0

   do i = size(record), 1, -1
      if (allocated(record(i)%label)) then
         if (label == record(i)%label) then
            pos = i
            exit
         end if
      end if
   end do

end function find

   !> Reallocate list of double arrays
pure subroutine resize(var, n)
   !> Instance of the array to be resized
   type(double_record), allocatable, intent(inout) :: var(:)
   !> Dimension of the final array size
   integer, intent(in), optional :: n

   type(double_record), allocatable :: tmp(:)
   integer :: this_size, new_size
   integer, parameter :: initial_size = 20

   if (allocated(var)) then
      this_size = size(var, 1)
      call move_alloc(var, tmp)
   else
      this_size = initial_size
   end if

   if (present(n)) then
      new_size = n
   else
      new_size = this_size + this_size/2 + 1
   end if

   allocate(var(new_size))

   if (allocated(tmp)) then
      this_size = min(size(tmp, 1), size(var, 1))
      var(:this_size) = tmp(:this_size)
      deallocate(tmp)
   end if

end subroutine resize

end module tblite_double_dictionary
