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

!> @file tblite/io/data/npz.f90
!> Provides NPZ data handling capabilities for tblite

!> Implementation of NPZ data handler
module tblite_io_data_npz
   use mctc_env, only: dp, i4, error_type, fatal_error
   use tblite_io_data_type, only: iodata_type, iodata_record
   use tblite_io_numpy, only : load_npz, save_npz
   use tblite_io_numpy_loadz, only : get_npz_descriptor
   use tblite_io_numpy_zip, only : zip_file, list_zip_file
   implicit none
   private

   type, extends(iodata_type), public :: npz_data
      !> NPZ file name handled by this backend
      character(len=:), allocatable :: filename
   contains
      !> List logical datasets stored in the NPZ archive.
      procedure :: list
      !> Get the shape of a logical dataset.
      procedure :: get_shape
      !> Load an integer rank-1 dataset.
      procedure :: load_i4_r1
      !> Load a real rank-1 dataset.
      procedure :: load_rdp_r1
      !> Load a real rank-2 dataset.
      procedure :: load_rdp_r2
      !> Load a real rank-3 dataset.
      procedure :: load_rdp_r3
      !> Save an integer rank-1 dataset.
      procedure :: save_i4_r1
      !> Save a real rank-1 dataset.
      procedure :: save_rdp_r1
      !> Save a real rank-2 dataset.
      procedure :: save_rdp_r2
      !> Save a real rank-3 dataset.
      procedure :: save_rdp_r3
   end type npz_data

contains

subroutine list(self, datasets, error)
   !> NPZ backend instance
   class(npz_data), intent(in) :: self
   !> Logical dataset records found in the archive
   type(iodata_record), allocatable, intent(out) :: datasets(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(zip_file) :: zip
   integer :: io, stat, irec
   character(len=:), allocatable :: msg

   open(newunit=io, file=self%filename, form="unformatted", access="stream", iostat=stat)
   call list_zip_file(io, self%filename, zip, stat, msg)
   close(io)
   if (stat /= 0) then
      if (.not.allocated(msg)) msg = "Failed to read zip file '"//self%filename//"'"
      call fatal_error(error, msg)
      return
   end if

   if (.not.allocated(zip%records)) then
      call fatal_error(error, "No records found in file '"//self%filename//"'")
      return
   end if

   allocate(datasets(size(zip%records)))
   do irec = 1, size(zip%records)
      associate(path => zip%records(irec)%path)
         if (len(path) > 4 .and. path(len(path)-3:) == ".npy") then
            datasets(irec)%name = path(:len(path)-4)
         else
            datasets(irec)%name = path
         end if
      end associate
   end do
end subroutine list

subroutine get_shape(self, dataset, shape, error)
   !> NPZ backend instance
   class(npz_data), intent(in) :: self
   !> Logical dataset name
   character(len=*), intent(in) :: dataset
   !> Dataset shape
   integer, allocatable, intent(out) :: shape(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=:), allocatable :: msg, vtype

   call get_npz_descriptor(self%filename, dataset, vtype, shape, stat, msg)
   if (stat /= 0) then
      call fatal_error(error, "Failed to read descriptor for dataset '"//dataset//"' from NPZ file '"//self%filename//"': "//msg)
      return
   end if
end subroutine get_shape

subroutine load_i4_r1(self, dataset, data, error)
   !> NPZ backend instance
   class(npz_data), intent(in) :: self
   !> Logical dataset name
   character(len=*), intent(in) :: dataset
   !> Data loaded from the dataset
   integer(i4), allocatable, intent(out) :: data(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=:), allocatable :: msg

   call load_npz(self%filename, dataset, data, stat, msg)
   if (stat /= 0) then
      call fatal_error(error, "Failed to load dataset '"//dataset//"' from npz file '"//self%filename//"': "//msg)
   end if
end subroutine load_i4_r1

subroutine load_rdp_r1(self, dataset, data, error)
   !> NPZ backend instance
   class(npz_data), intent(in) :: self
   !> Logical dataset name
   character(len=*), intent(in) :: dataset
   !> Data loaded from the dataset
   real(dp), allocatable, intent(out) :: data(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=:), allocatable :: msg

   call load_npz(self%filename, dataset, data, stat, msg)
   if (stat /= 0) then
      call fatal_error(error, "Failed to load dataset '"//dataset//"' from npz file '"//self%filename//"': "//msg)
   end if
end subroutine load_rdp_r1

subroutine load_rdp_r2(self, dataset, data, error)
   !> NPZ backend instance
   class(npz_data), intent(in) :: self
   !> Logical dataset name
   character(len=*), intent(in) :: dataset
   !> Data loaded from the dataset
   real(dp), allocatable, intent(out) :: data(:,:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=:), allocatable :: msg

   call load_npz(self%filename, dataset, data, stat, msg)
   if (stat /= 0) then
      call fatal_error(error, "Failed to load dataset '"//dataset//"' from npz file '"//self%filename//"': "//msg)
   end if
end subroutine load_rdp_r2

subroutine load_rdp_r3(self, dataset, data, error)
   !> NPZ backend instance
   class(npz_data), intent(in) :: self
   !> Logical dataset name
   character(len=*), intent(in) :: dataset
   !> Data loaded from the dataset
   real(dp), allocatable, intent(out) :: data(:,:,:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=:), allocatable :: msg

   call load_npz(self%filename, dataset, data, stat, msg)
   if (stat /= 0) then
      call fatal_error(error, "Failed to load dataset '"//dataset//"' from npz file '"//self%filename//"': "//msg)
   end if
end subroutine load_rdp_r3

subroutine save_i4_r1(self, dataset, data, error)
   !> NPZ backend instance
   class(npz_data), intent(in) :: self
   !> Logical dataset name
   character(len=*), intent(in) :: dataset
   !> Data to save
   integer(i4), allocatable, intent(in) :: data(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=:), allocatable :: msg

   call save_npz(self%filename, dataset, data, stat, msg)
   if (stat /= 0) then
      call fatal_error(error, "Failed to save dataset '"//dataset//"' to npz file '"//self%filename//"': "//msg)
   end if
end subroutine save_i4_r1

subroutine save_rdp_r1(self, dataset, data, error)
   !> NPZ backend instance
   class(npz_data), intent(in) :: self
   !> Logical dataset name
   character(len=*), intent(in) :: dataset
   !> Data to save
   real(dp), allocatable, intent(in) :: data(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=:), allocatable :: msg

   call save_npz(self%filename, dataset, data, stat, msg)
   if (stat /= 0) then
      call fatal_error(error, "Failed to save dataset '"//dataset//"' to npz file '"//self%filename//"': "//msg)
   end if
end subroutine save_rdp_r1

subroutine save_rdp_r2(self, dataset, data, error)
   !> NPZ backend instance
   class(npz_data), intent(in) :: self
   !> Logical dataset name
   character(len=*), intent(in) :: dataset
   !> Data to save
   real(dp), allocatable, intent(in) :: data(:,:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=:), allocatable :: msg

   call save_npz(self%filename, dataset, data, stat, msg)
   if (stat /= 0) then
      call fatal_error(error, "Failed to save dataset '"//dataset//"' to npz file '"//self%filename//"': "//msg)
   end if
end subroutine save_rdp_r2

subroutine save_rdp_r3(self, dataset, data, error)
   !> NPZ backend instance
   class(npz_data), intent(in) :: self
   !> Logical dataset name
   character(len=*), intent(in) :: dataset
   !> Data to save
   real(dp), allocatable, intent(in) :: data(:,:,:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=:), allocatable :: msg
   call save_npz(self%filename, dataset, data, stat, msg)
   if (stat /= 0) then
      call fatal_error(error, "Failed to save dataset '"//dataset//"' to npz file '"//self%filename//"': "//msg)
   end if
end subroutine save_rdp_r3
end module tblite_io_data_npz
