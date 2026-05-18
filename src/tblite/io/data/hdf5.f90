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

!> @file tblite/io/data/hdf5.f90
!> Provides HDF5 data handling capabilities for tblite

!> Implementation of HDF5 data handler
module tblite_io_data_hdf5
   use mctc_env, only: dp, i4, error_type, fatal_error
   use tblite_io_data_type, only: iodata_type, iodata_record
   use tblite_io_hdf5, only : get_hdf5_descriptor, list_hdf5_datasets, load_hdf5, save_hdf5
   implicit none
   private

   type, extends(iodata_type), public :: hdf5_data
      !> HDF5 file name handled by this backend
      character(len=:), allocatable :: filename
   contains
      !> List logical datasets stored in the HDF5 file.
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
   end type hdf5_data

contains

subroutine list(self, datasets, error)
   !> HDF5 backend instance
   class(hdf5_data), intent(in) :: self
   !> Logical dataset records found in the file
   type(iodata_record), allocatable, intent(out) :: datasets(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat, irec
   character(len=:), allocatable :: msg, paths(:)

   call list_hdf5_datasets(self%filename, paths, stat, msg)
   if (stat /= 0) then
      if (.not.allocated(msg)) msg = "Failed to read HDF5 file '"//self%filename//"'"
      call fatal_error(error, msg)
      return
   end if
   if (.not.allocated(paths)) then
      call fatal_error(error, "No records found in file '"//self%filename//"'")
      return
   end if

   allocate(datasets(size(paths)))
   do irec = 1, size(paths)
      block
         character(len=:), allocatable :: path
         path = trim(paths(irec))
         datasets(irec)%name = path
      end block
   end do
end subroutine list

subroutine get_shape(self, dataset, shape, error)
   !> HDF5 backend instance
   class(hdf5_data), intent(in) :: self
   !> Logical dataset name
   character(len=*), intent(in) :: dataset
   !> Dataset shape
   integer, allocatable, intent(out) :: shape(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=:), allocatable :: msg

   call get_hdf5_descriptor(self%filename, dataset, shape, stat, msg)
   if (stat /= 0) then
      call fatal_error(error, "Failed to read descriptor for dataset '"//dataset//"' from HDF5 file '"//self%filename//"': "//msg)
      return
   end if
end subroutine get_shape

subroutine load_i4_r1(self, dataset, data, error)
   !> HDF5 backend instance
   class(hdf5_data), intent(in) :: self
   !> Logical dataset name
   character(len=*), intent(in) :: dataset
   !> Data loaded from the dataset
   integer(i4), allocatable, intent(out) :: data(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=:), allocatable :: msg

   call load_hdf5(self%filename, dataset, data, stat, msg)
   if (stat /= 0) then
      call fatal_error(error, "Failed to load dataset '"//dataset//"' from HDF5 file '"//self%filename//"': "//msg)
   end if
end subroutine load_i4_r1

subroutine load_rdp_r1(self, dataset, data, error)
   !> HDF5 backend instance
   class(hdf5_data), intent(in) :: self
   !> Logical dataset name
   character(len=*), intent(in) :: dataset
   !> Data loaded from the dataset
   real(dp), allocatable, intent(out) :: data(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=:), allocatable :: msg

   call load_hdf5(self%filename, dataset, data, stat, msg)
   if (stat /= 0) then
      call fatal_error(error, "Failed to load dataset '"//dataset//"' from HDF5 file '"//self%filename//"': "//msg)
   end if
end subroutine load_rdp_r1

subroutine load_rdp_r2(self, dataset, data, error)
   !> HDF5 backend instance
   class(hdf5_data), intent(in) :: self
   !> Logical dataset name
   character(len=*), intent(in) :: dataset
   !> Data loaded from the dataset
   real(dp), allocatable, intent(out) :: data(:,:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=:), allocatable :: msg

   call load_hdf5(self%filename, dataset, data, stat, msg)
   if (stat /= 0) then
      call fatal_error(error, "Failed to load dataset '"//dataset//"' from HDF5 file '"//self%filename//"': "//msg)
   end if
end subroutine load_rdp_r2

subroutine load_rdp_r3(self, dataset, data, error)
   !> HDF5 backend instance
   class(hdf5_data), intent(in) :: self
   !> Logical dataset name
   character(len=*), intent(in) :: dataset
   !> Data loaded from the dataset
   real(dp), allocatable, intent(out) :: data(:,:,:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=:), allocatable :: msg

   call load_hdf5(self%filename, dataset, data, stat, msg)
   if (stat /= 0) then
      call fatal_error(error, "Failed to load dataset '"//dataset//"' from HDF5 file '"//self%filename//"': "//msg)
   end if
end subroutine load_rdp_r3

subroutine save_i4_r1(self, dataset, data, error)
   !> HDF5 backend instance
   class(hdf5_data), intent(in) :: self
   !> Logical dataset name
   character(len=*), intent(in) :: dataset
   !> Data to save
   integer(i4), allocatable, intent(in) :: data(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=:), allocatable :: msg

   call save_hdf5(self%filename, dataset, data, stat, msg)
   if (stat /= 0) then
      call fatal_error(error, "Failed to save dataset '"//dataset//"' to HDF5 file '"//self%filename//"': "//msg)
   end if
end subroutine save_i4_r1

subroutine save_rdp_r1(self, dataset, data, error)
   !> HDF5 backend instance
   class(hdf5_data), intent(in) :: self
   !> Logical dataset name
   character(len=*), intent(in) :: dataset
   !> Data to save
   real(dp), allocatable, intent(in) :: data(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=:), allocatable :: msg

   call save_hdf5(self%filename, dataset, data, stat, msg)
   if (stat /= 0) then
      call fatal_error(error, "Failed to save dataset '"//dataset//"' to HDF5 file '"//self%filename//"': "//msg)
   end if
end subroutine save_rdp_r1

subroutine save_rdp_r2(self, dataset, data, error)
   !> HDF5 backend instance
   class(hdf5_data), intent(in) :: self
   !> Logical dataset name
   character(len=*), intent(in) :: dataset
   !> Data to save
   real(dp), allocatable, intent(in) :: data(:,:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=:), allocatable :: msg

   call save_hdf5(self%filename, dataset, data, stat, msg)
   if (stat /= 0) then
      call fatal_error(error, "Failed to save dataset '"//dataset//"' to HDF5 file '"//self%filename//"': "//msg)
   end if
end subroutine save_rdp_r2

subroutine save_rdp_r3(self, dataset, data, error)
   !> HDF5 backend instance
   class(hdf5_data), intent(in) :: self
   !> Logical dataset name
   character(len=*), intent(in) :: dataset
   !> Data to save
   real(dp), allocatable, intent(in) :: data(:,:,:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=:), allocatable :: msg
   call save_hdf5(self%filename, dataset, data, stat, msg)
   if (stat /= 0) then
      call fatal_error(error, "Failed to save dataset '"//dataset//"' to HDF5 file '"//self%filename//"': "//msg)
   end if
end subroutine save_rdp_r3

end module tblite_io_data_hdf5