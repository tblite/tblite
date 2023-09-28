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
      procedure :: copy_record
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
   end type double_dictionary_type


contains

subroutine remove_entry_index(self, index)
   class(double_dictionary_type) :: self
   integer :: index, old_n, i, it
   type(double_dictionary_type) :: tmp

   if (index > self%n) return
   old_n = self%n
   self%n = self%n - 1

   tmp = self
   deallocate(self%record)
   allocate(self%record(self%n))
   it = 1
   do i = 1, old_n
      if (i == index) cycle
      self%record(it) = tmp%record(it)
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
      self%record(it) = double_record(label)
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
   integer, parameter :: initial_size = 0

   if (allocated(var)) then
      this_size = size(var, 1)
      call move_alloc(var, tmp)
   else
      this_size = initial_size
   end if

   if (present(n)) then
      new_size = n
   else
      new_size = this_size + 1
   end if

   allocate(var(new_size))

   if (allocated(tmp)) then
      this_size = min(size(tmp, 1), size(var, 1))
      var(:this_size) = tmp(:this_size)
      deallocate(tmp)
   end if

end subroutine resize

end module tblite_double_dictionary
