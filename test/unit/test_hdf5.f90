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

module test_hdf5
   use mctc_env, only : dp, error_type, i4
   use mctc_env_testing, only : new_unittest, unittest_type, check
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_features, only : get_tblite_feature
   use tblite_io_hdf5, only : save_hdf5, load_hdf5, get_hdf5_descriptor, list_hdf5_datasets
   implicit none
   private

   public :: collect_hdf5

contains

subroutine collect_hdf5(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   if (get_tblite_feature("hdf5")) then
      testsuite = [ &
         new_unittest("write-hdf5-i4-r1", test_write_hdf5_i4_rank1), &
         new_unittest("write-hdf5-rdp-r1", test_write_hdf5_rdp_rank1), &
         new_unittest("write-hdf5-rdp-r2", test_write_hdf5_rdp_rank2), &
         new_unittest("write-hdf5-rdp-r3", test_write_hdf5_rdp_rank3), &
         new_unittest("write-hdf5-multi", test_write_hdf5_multi), &
         new_unittest("hdf5-descriptor", test_hdf5_descriptor), &
         new_unittest("double-dictionary-hdf5", test_double_dictionary_hdf5) &
      ]
   else
      testsuite = [new_unittest("hdf5-disabled", test_hdf5_disabled)]
   end if
end subroutine collect_hdf5

subroutine test_hdf5_disabled(error)
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=:), allocatable :: msg

   call save_hdf5(".hdf5-disabled.h5", "test", [1_i4], stat, msg)
   call check(error, stat /= 0, "HDF5 write unexpectedly succeeded without HDF5 support")
   if (allocated(error)) return
   call check(error, allocated(msg), "HDF5 disabled error message was not returned")
end subroutine test_hdf5_disabled

subroutine test_write_hdf5_i4_rank1(error)
   type(error_type), allocatable, intent(out) :: error

   integer :: stat, i
   character(len=*), parameter :: filename = ".write-i4-r1.h5", varname = "test"
   integer(i4), allocatable :: input(:), output(:)
   character(len=:), allocatable :: msg

   allocate(input(72))
   input = [(2*i, i=1, size(input))]
   call save_hdf5(filename, varname, input, stat, msg)
   call check(error, stat, "Writing of HDF5 file failed", msg)
   if (allocated(error)) return

   call load_hdf5(filename, varname, output, stat, msg)
   call delete_file(filename)
   call check(error, stat, "Reading of HDF5 file failed", msg)
   if (allocated(error)) return

   call check(error, size(output), size(input))
   if (allocated(error)) return
   call check(error, all((output - input) == 0), "Precision loss when rereading array")
end subroutine test_write_hdf5_i4_rank1

subroutine test_write_hdf5_rdp_rank1(error)
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=*), parameter :: filename = ".write-rdp-r1.h5", varname = "test"
   real(dp), allocatable :: input(:), output(:)
   character(len=:), allocatable :: msg

   allocate(input(51))
   call random_number(input)
   call save_hdf5(filename, varname, input, stat, msg)
   call check(error, stat, "Writing of HDF5 file failed", msg)
   if (allocated(error)) return

   call load_hdf5(filename, varname, output, stat, msg)
   call delete_file(filename)
   call check(error, stat, "Reading of HDF5 file failed", msg)
   if (allocated(error)) return

   call check(error, size(output), size(input))
   if (allocated(error)) return
   call check(error, all(abs(output - input) <= epsilon(1.0_dp)), &
      & "Precision loss when rereading array")
end subroutine test_write_hdf5_rdp_rank1

subroutine test_write_hdf5_rdp_rank2(error)
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=*), parameter :: filename = ".write-rdp-r2.h5", varname = "test"
   real(dp), allocatable :: input(:, :), output(:, :)
   character(len=:), allocatable :: msg

   allocate(input(5, 13))
   call random_number(input)
   call save_hdf5(filename, varname, input, stat, msg)
   call check(error, stat, "Writing of HDF5 file failed", msg)
   if (allocated(error)) return

   call load_hdf5(filename, varname, output, stat, msg)
   call delete_file(filename)
   call check(error, stat, "Reading of HDF5 file failed", msg)
   if (allocated(error)) return

   call check(error, all(shape(output) == shape(input)), "Shape changed when rereading array")
   if (allocated(error)) return
   call check(error, all(abs(output - input) <= epsilon(1.0_dp)), &
      & "Precision loss when rereading array")
end subroutine test_write_hdf5_rdp_rank2

subroutine test_write_hdf5_rdp_rank3(error)
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=*), parameter :: filename = ".write-rdp-r3.h5", varname = "test"
   real(dp), allocatable :: input(:, :, :), output(:, :, :)
   character(len=:), allocatable :: msg

   allocate(input(5, 7, 3))
   call random_number(input)
   call save_hdf5(filename, varname, input, stat, msg)
   call check(error, stat, "Writing of HDF5 file failed", msg)
   if (allocated(error)) return

   call load_hdf5(filename, varname, output, stat, msg)
   call delete_file(filename)
   call check(error, stat, "Reading of HDF5 file failed", msg)
   if (allocated(error)) return

   call check(error, all(shape(output) == shape(input)), "Shape changed when rereading array")
   if (allocated(error)) return
   call check(error, all(abs(output - input) <= epsilon(1.0_dp)), &
      & "Precision loss when rereading array")
end subroutine test_write_hdf5_rdp_rank3

subroutine test_write_hdf5_multi(error)
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=*), parameter :: filename = ".mwrite-rdp.h5"
   real(dp), allocatable :: input1(:), input2(:, :), output1(:), output2(:, :)
   character(len=:), allocatable :: msg

   allocate(input1(12), input2(3, 5))
   call random_number(input1)
   call random_number(input2)
   call save_hdf5(filename, "array1", input1, stat, msg)
   if (stat == 0) call save_hdf5(filename, "array2", input2, stat, msg)
   call check(error, stat, "Writing of HDF5 file failed", msg)
   if (allocated(error)) return

   call load_hdf5(filename, "array1", output1, stat, msg)
   if (stat == 0) call load_hdf5(filename, "array2", output2, stat, msg)
   call delete_file(filename)
   call check(error, stat, "Reading of HDF5 file failed", msg)
   if (allocated(error)) return

   call check(error, all(abs(output1 - input1) <= epsilon(1.0_dp)), &
      & "Precision loss when rereading rank-1 array")
   if (allocated(error)) return
   call check(error, all(abs(output2 - input2) <= epsilon(1.0_dp)), &
      & "Precision loss when rereading rank-2 array")
end subroutine test_write_hdf5_multi

subroutine test_hdf5_descriptor(error)
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=*), parameter :: filename = ".descriptor.h5"
   integer, allocatable :: vshape(:)
   real(dp), allocatable :: input(:, :, :)
   character(len=:), allocatable :: msg, paths(:)

   allocate(input(4, 6, 9), source=0.0_dp)
   call save_hdf5(filename, "tblite0_test", input, stat, msg)
   call check(error, stat, "Writing of HDF5 file failed", msg)
   if (allocated(error)) return

   call get_hdf5_descriptor(filename, "tblite0_test", vshape, stat, msg)
   call check(error, stat, "Reading HDF5 descriptor failed", msg)
   if (allocated(error)) return
   call check(error, all(vshape == shape(input)), "Descriptor shape does not match array shape")
   if (allocated(error)) return

   call list_hdf5_datasets(filename, paths, stat, msg)
   call delete_file(filename)
   call check(error, stat, "Listing HDF5 datasets failed", msg)
   if (allocated(error)) return
   call check(error, size(paths), 1)
   if (allocated(error)) return
   call check(error, trim(paths(1)) == "tblite0_test", "Dataset name not found in HDF5 file")
end subroutine test_hdf5_descriptor

subroutine test_double_dictionary_hdf5(error)
   type(error_type), allocatable, intent(out) :: error

   type(double_dictionary_type) :: dict1, dict2
   character(len=*), parameter :: filename = ".read-write-ddict.h5"

   call fill_test_dict(dict1)
   call dict1%dump(filename, error)
   if (allocated(error)) return

   call dict2%load(filename, error)
   call delete_file(filename)
   if (allocated(error)) return

   call check(error, dict1 == dict2)
end subroutine test_double_dictionary_hdf5

subroutine fill_test_dict(dict)
   type(double_dictionary_type), intent(inout) :: dict

   real(dp), allocatable :: array1(:), array2(:, :), array3(:, :, :)

   allocate(array1(4), source=1.0_dp)
   allocate(array2(4, 6), source=2.0_dp)
   allocate(array3(4, 6, 9), source=3.0_dp)

   call dict%add_entry("test1", array1)
   call dict%add_entry("test2", array2)
   call dict%add_entry("test3", array3)
end subroutine fill_test_dict

subroutine delete_file(filename)
   character(len=*), intent(in) :: filename

   integer :: io
   logical :: exist

   inquire(file=filename, exist=exist)
   if (exist) then
      open(newunit=io, file=filename)
      close(io, status="delete")
   end if
end subroutine delete_file

end module test_hdf5
