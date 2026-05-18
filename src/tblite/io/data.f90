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

!> @dir tblite/io/data/
!> Directory containing data handling modules for tblite

!> @file tblite/io/data.f90
!> Provides format-independent array storage handlers.

!> Factory module for selecting NPZ or HDF5 structured array handlers.
module tblite_io_data
   use mctc_env, only : error_type, fatal_error
   use tblite_io_data_type, only: iodata_type, iodata_record
   use tblite_io_data_hdf5, only: hdf5_data
   use tblite_io_data_npz, only: npz_data
   implicit none
   private

   public :: iodata_type, iodata_record
   public :: open_iodata_handler

contains

!> Open an array data handler appropriate for the file extension.
subroutine open_iodata_handler(iodata, filename, error)
   !> File name used to select the storage backend
   character(len=*), intent(in) :: filename
   !> Allocated storage backend instance
   class(iodata_type), allocatable, intent(out) :: iodata
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   logical :: is_hdf5

   is_hdf5 = is_hdf5_file(filename)

   if (is_hdf5) then
      block
         type(hdf5_data), allocatable :: hdf5_handler
         allocate(hdf5_handler)
         hdf5_handler%filename = filename
         call move_alloc(hdf5_handler, iodata)
      end block
   else
      block
         type(npz_data), allocatable :: npz_handler
         allocate(npz_handler)
         npz_handler%filename = filename
         call move_alloc(npz_handler, iodata)
      end block
   end if
end subroutine open_iodata_handler

pure function is_hdf5_file(filename) result(is_hdf5)
   !> File name to inspect
   character(len=*), intent(in) :: filename
   !> Whether the file name uses an HDF5 extension
   logical :: is_hdf5

   integer :: n

   n = len_trim(filename)
   is_hdf5 = (n >= 3 .and. filename(n-2:n) == ".h5") .or. &
      & (n >= 5 .and. filename(n-4:n) == ".hdf5")
end function is_hdf5_file

end module tblite_io_data