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

!> @file tblite/io/data/type.f90
!> Defines a common interface for structured array storage backends.

!> Abstract data handler API used by restart and post-processing output.
module tblite_io_data_type
   use mctc_env, only: dp, i4, error_type
   implicit none
   private

   !> Abstract interface for file-backed named array storage.
   !>
   !> Concrete implementations provide the same load, save, shape, and listing
   !> operations for different container formats such as NPZ and HDF5.
   type, abstract, public :: iodata_type
   contains
      !> List logical dataset names available in the file.
      procedure(list), deferred :: list
      !> Query the rank and extent of a named dataset.
      procedure(get_shape), deferred :: get_shape
      !> Load a named dataset.
      generic :: load => load_i4_r1, load_rdp_r1, load_rdp_r2, load_rdp_r3
      !> Load a named integer rank-1 dataset.
      procedure(load_i4_r1), deferred :: load_i4_r1
      !> Load a named real rank-1 dataset.
      procedure(load_rdp_r1), deferred :: load_rdp_r1
      !> Load a named real rank-2 dataset.
      procedure(load_rdp_r2), deferred :: load_rdp_r2
      !> Load a named real rank-3 dataset.
      procedure(load_rdp_r3), deferred :: load_rdp_r3
      !> Save a named dataset.
      generic :: save => save_i4_r1, save_rdp_r1, save_rdp_r2, save_rdp_r3
      !> Save a named integer rank-1 dataset.
      procedure(save_i4_r1), deferred :: save_i4_r1
      !> Save a named real rank-1 dataset.
      procedure(save_rdp_r1), deferred :: save_rdp_r1
      !> Save a named real rank-2 dataset.
      procedure(save_rdp_r2), deferred :: save_rdp_r2
      !> Save a named real rank-3 dataset.
      procedure(save_rdp_r3), deferred :: save_rdp_r3
   end type iodata_type

   !> Description of a dataset stored in an iodata container.
   type, public :: iodata_record
      !> Logical dataset name, without container-specific suffixes.
      character(len=:), allocatable :: name
   end type iodata_record

   abstract interface
      !> List logical datasets in a storage backend.
      subroutine list(self, datasets, error)
         import :: iodata_type, iodata_record, error_type
         !> Storage backend instance
         class(iodata_type), intent(in) :: self
         !> Dataset records found in the file
         type(iodata_record), allocatable, intent(out) :: datasets(:)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine list

      !> Get the shape of a named dataset.
      subroutine get_shape(self, dataset, shape, error)
         import :: iodata_type, error_type
         !> Storage backend instance
         class(iodata_type), intent(in) :: self
         !> Logical dataset name
         character(len=*), intent(in) :: dataset
         !> Dataset shape
         integer, allocatable, intent(out) :: shape(:)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_shape

      !> Load an integer rank-1 dataset.
      subroutine load_i4_r1(self, dataset, data, error)
         import :: iodata_type, error_type, i4
         !> Storage backend instance
         class(iodata_type), intent(in) :: self
         !> Logical dataset name
         character(len=*), intent(in) :: dataset
         !> Data loaded from the dataset
         integer(i4), allocatable, intent(out) :: data(:)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine load_i4_r1

      !> Load a real rank-1 dataset.
      subroutine load_rdp_r1(self, dataset, data, error)
         import :: iodata_type, error_type, dp
         !> Storage backend instance
         class(iodata_type), intent(in) :: self
         !> Logical dataset name
         character(len=*), intent(in) :: dataset
         !> Data loaded from the dataset
         real(dp), allocatable, intent(out) :: data(:)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine load_rdp_r1

      !> Load a real rank-2 dataset.
      subroutine load_rdp_r2(self, dataset, data, error)
         import :: iodata_type, error_type, dp
         !> Storage backend instance
         class(iodata_type), intent(in) :: self
         !> Logical dataset name
         character(len=*), intent(in) :: dataset
         !> Data loaded from the dataset
         real(dp), allocatable, intent(out) :: data(:,:)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine load_rdp_r2

      !> Load a real rank-3 dataset.
      subroutine load_rdp_r3(self, dataset, data, error)
         import :: iodata_type, error_type, dp
         !> Storage backend instance
         class(iodata_type), intent(in) :: self
         !> Logical dataset name
         character(len=*), intent(in) :: dataset
         !> Data loaded from the dataset
         real(dp), allocatable, intent(out) :: data(:,:,:)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine load_rdp_r3

      !> Save an integer rank-1 dataset.
      subroutine save_i4_r1(self, dataset, data, error)
         import :: iodata_type, error_type, i4
         !> Storage backend instance
         class(iodata_type), intent(in) :: self
         !> Logical dataset name
         character(len=*), intent(in) :: dataset
         !> Data to save
         integer(i4), allocatable, intent(in) :: data(:)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine save_i4_r1

      !> Save a real rank-1 dataset.
      subroutine save_rdp_r1(self, dataset, data, error)
         import :: iodata_type, error_type, dp
         !> Storage backend instance
         class(iodata_type), intent(in) :: self
         !> Logical dataset name
         character(len=*), intent(in) :: dataset
         !> Data to save
         real(dp), allocatable, intent(in) :: data(:)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine save_rdp_r1

      !> Save a real rank-2 dataset.
      subroutine save_rdp_r2(self, dataset, data, error)
         import :: iodata_type, error_type, dp
         !> Storage backend instance
         class(iodata_type), intent(in) :: self
         !> Logical dataset name
         character(len=*), intent(in) :: dataset
         !> Data to save
         real(dp), allocatable, intent(in) :: data(:,:)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine save_rdp_r2

      !> Save a real rank-3 dataset.
      subroutine save_rdp_r3(self, dataset, data, error)
         import :: iodata_type, error_type, dp
         !> Storage backend instance
         class(iodata_type), intent(in) :: self
         !> Logical dataset name
         character(len=*), intent(in) :: dataset
         !> Data to save
         real(dp), allocatable, intent(in) :: data(:,:,:)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine save_rdp_r3
   end interface
end module tblite_io_data_type