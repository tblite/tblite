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

#ifndef TBLITE_HAS_HDF5
#define TBLITE_HAS_HDF5 0
#endif

!> @file tblite/io/hdf5.f90
!> Provides HDF5 input and output routines

!> Implementation of HDF5 input and output routines
module tblite_io_hdf5
   use mctc_env, only : dp, i4
#if TBLITE_HAS_HDF5
   use hdf5, only : hid_t, hsize_t, h5open_f, h5close_f, h5fopen_f, h5fcreate_f, &
      & h5fclose_f, h5dopen_f, h5dcreate_f, h5dclose_f, h5dread_f, h5dwrite_f, &
      & h5dget_space_f, h5screate_simple_f, h5sclose_f, h5sget_simple_extent_ndims_f, &
      & h5sget_simple_extent_dims_f, h5lexists_f, h5ldelete_f, h5gn_members_f, &
      & h5gget_obj_info_idx_f, H5F_ACC_RDONLY_F, H5F_ACC_RDWR_F, H5F_ACC_TRUNC_F, &
      & H5P_DEFAULT_F, H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER, H5G_DATASET_F
#endif
   implicit none
   private

   public :: get_hdf5_descriptor, list_hdf5_datasets, load_hdf5, save_hdf5

   !> Interface for loading hdf5 files
   interface load_hdf5
      module procedure load_hdf5_i4_r1
      module procedure load_hdf5_rdp_r1
      module procedure load_hdf5_rdp_r2
      module procedure load_hdf5_rdp_r3
   end interface load_hdf5

   !> Interface for saving hdf5 files
   interface save_hdf5
      module procedure save_hdf5_i4_r1
      module procedure save_hdf5_rdp_r1
      module procedure save_hdf5_rdp_r2
      module procedure save_hdf5_rdp_r3
   end interface save_hdf5

   character(len=*), parameter :: nl = new_line('a')

contains

subroutine get_hdf5_descriptor(filename, varname, vshape, iostat, iomsg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> Name of the dataset to inspect
   character(len=*), intent(in) :: varname
   !> Shape of the HDF5 dataset
   integer, allocatable, intent(out) :: vshape(:)
   !> Status of the read operation
   integer, intent(out), optional :: iostat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out), optional :: iomsg

   integer :: stat
   character(len=:), allocatable :: msg
#if TBLITE_HAS_HDF5
   integer :: rank
   integer(hsize_t), allocatable :: dims(:), maxdims(:)
   integer(hid_t) :: file_id, dset_id, space_id

   file_id = -1_hid_t
   dset_id = -1_hid_t
   space_id = -1_hid_t
   call hdf5_open(stat, msg)
   if (stat == 0) call open_file_read(filename, file_id, stat, msg)
   if (stat == 0) call h5dopen_f(file_id, varname, dset_id, stat)
   if (stat /= 0 .and. .not.allocated(msg)) msg = "Dataset '"//filename//"::"//varname//"' not found"
   if (stat == 0) call h5dget_space_f(dset_id, space_id, stat)
   if (stat == 0) call h5sget_simple_extent_ndims_f(space_id, rank, stat)
   if (stat == 0) then
      allocate(dims(rank), maxdims(rank))
      call h5sget_simple_extent_dims_f(space_id, dims, maxdims, stat)
      if (stat >= 0) stat = 0
   end if
   if (stat == 0) vshape = int(dims)
   call close_space(space_id)
   call close_dataset(dset_id)
   call close_file(file_id)
   call hdf5_close()
#else
   stat = -1
   msg = "HDF5 support is not available in this build of tblite."
#endif
   call handle_iostat(stat, msg, filename//"::"//varname, iostat, iomsg, "read")
end subroutine get_hdf5_descriptor

subroutine list_hdf5_datasets(filename, paths, iostat, iomsg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> Names of root-level datasets in the file
   character(len=:), allocatable, intent(out) :: paths(:)
   !> Status of the read operation
   integer, intent(out), optional :: iostat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out), optional :: iomsg

   integer :: stat
   character(len=:), allocatable :: msg
#if TBLITE_HAS_HDF5
   integer :: irec, nrec, obj_type, nset
   integer(hid_t) :: file_id
   character(len=256), allocatable :: buffer(:)

   file_id = -1_hid_t
   call hdf5_open(stat, msg)
   if (stat == 0) call open_file_read(filename, file_id, stat, msg)
   if (stat == 0) call h5gn_members_f(file_id, "/", nrec, stat)
   if (stat == 0) then
      allocate(buffer(nrec))
      nset = 0
      do irec = 0, nrec - 1
         call h5gget_obj_info_idx_f(file_id, "/", irec, buffer(irec + 1), obj_type, stat)
         if (stat /= 0) exit
         if (obj_type == H5G_DATASET_F) nset = nset + 1
      end do
   end if
   if (stat == 0) then
      allocate(character(len=256) :: paths(nset))
      nset = 0
      do irec = 0, nrec - 1
         call h5gget_obj_info_idx_f(file_id, "/", irec, buffer(irec + 1), obj_type, stat)
         if (stat /= 0) exit
         if (obj_type == H5G_DATASET_F) then
            nset = nset + 1
            paths(nset) = trim(buffer(irec + 1))
         end if
      end do
   end if
   call close_file(file_id)
   call hdf5_close()
#else
   stat = -1
   msg = "HDF5 support is not available in this build of tblite."
#endif
   call handle_iostat(stat, msg, filename, iostat, iomsg, "read")
end subroutine list_hdf5_datasets

subroutine load_hdf5_i4_r1(filename, varname, array, iostat, iomsg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> Name of the dataset to load
   character(len=*), intent(in) :: varname
   !> Array to load the data into
   integer(i4), allocatable, intent(out) :: array(:)
   !> Status of the read operation
   integer, intent(out), optional :: iostat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out), optional :: iomsg

   integer :: stat
   character(len=:), allocatable :: msg
#if TBLITE_HAS_HDF5
   integer(hsize_t) :: dims(1)
   call read_shape(filename, varname, dims, stat, msg)
   if (stat == 0) allocate(array(dims(1)), stat=stat)
   if (stat == 0) call read_i4_r1(filename, varname, array, dims, stat, msg)
#else
   stat = -1
   msg = "HDF5 support is not available in this build of tblite."
#endif
   call handle_iostat(stat, msg, filename//"::"//varname, iostat, iomsg, "read")
end subroutine load_hdf5_i4_r1

subroutine load_hdf5_rdp_r1(filename, varname, array, iostat, iomsg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> Name of the dataset to load
   character(len=*), intent(in) :: varname
   !> Array to load the data into
   real(dp), allocatable, intent(out) :: array(:)
   !> Status of the read operation
   integer, intent(out), optional :: iostat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out), optional :: iomsg

   integer :: stat
   character(len=:), allocatable :: msg
#if TBLITE_HAS_HDF5
   integer(hsize_t) :: dims(1)
   call read_shape(filename, varname, dims, stat, msg)
   if (stat == 0) allocate(array(dims(1)), stat=stat)
   if (stat == 0) call read_rdp_r1(filename, varname, array, dims, stat, msg)
#else
   stat = -1
   msg = "HDF5 support is not available in this build of tblite."
#endif
   call handle_iostat(stat, msg, filename//"::"//varname, iostat, iomsg, "read")
end subroutine load_hdf5_rdp_r1

subroutine load_hdf5_rdp_r2(filename, varname, array, iostat, iomsg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> Name of the dataset to load
   character(len=*), intent(in) :: varname
   !> Array to load the data into
   real(dp), allocatable, intent(out) :: array(:, :)
   !> Status of the read operation
   integer, intent(out), optional :: iostat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out), optional :: iomsg

   integer :: stat
   character(len=:), allocatable :: msg
#if TBLITE_HAS_HDF5
   integer(hsize_t) :: dims(2)
   call read_shape(filename, varname, dims, stat, msg)
   if (stat == 0) allocate(array(dims(1), dims(2)), stat=stat)
   if (stat == 0) call read_rdp_r2(filename, varname, array, dims, stat, msg)
#else
   stat = -1
   msg = "HDF5 support is not available in this build of tblite."
#endif
   call handle_iostat(stat, msg, filename//"::"//varname, iostat, iomsg, "read")
end subroutine load_hdf5_rdp_r2

subroutine load_hdf5_rdp_r3(filename, varname, array, iostat, iomsg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> Name of the dataset to load
   character(len=*), intent(in) :: varname
   !> Array to load the data into
   real(dp), allocatable, intent(out) :: array(:, :, :)
   !> Status of the read operation
   integer, intent(out), optional :: iostat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out), optional :: iomsg

   integer :: stat
   character(len=:), allocatable :: msg
#if TBLITE_HAS_HDF5
   integer(hsize_t) :: dims(3)
   call read_shape(filename, varname, dims, stat, msg)
   if (stat == 0) allocate(array(dims(1), dims(2), dims(3)), stat=stat)
   if (stat == 0) call read_rdp_r3(filename, varname, array, dims, stat, msg)
#else
   stat = -1
   msg = "HDF5 support is not available in this build of tblite."
#endif
   call handle_iostat(stat, msg, filename//"::"//varname, iostat, iomsg, "read")
end subroutine load_hdf5_rdp_r3

subroutine save_hdf5_i4_r1(filename, varname, array, iostat, iomsg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> Name of the dataset to save
   character(len=*), intent(in) :: varname
   !> Array to save
   integer(i4), intent(in) :: array(:)
   !> Status of the write operation
   integer, intent(out), optional :: iostat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out), optional :: iomsg

   integer :: stat
   character(len=:), allocatable :: msg
#if TBLITE_HAS_HDF5
   integer(hsize_t) :: dims(1)
   dims = shape(array, kind=hsize_t)
   call write_i4_r1(filename, varname, array, dims, stat, msg)
#else
   stat = -1
   msg = "HDF5 support is not available in this build of tblite."
#endif
   call handle_iostat(stat, msg, filename//"::"//varname, iostat, iomsg, "write")
end subroutine save_hdf5_i4_r1

subroutine save_hdf5_rdp_r1(filename, varname, array, iostat, iomsg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> Name of the dataset to save
   character(len=*), intent(in) :: varname
   !> Array to save
   real(dp), intent(in) :: array(:)
   !> Status of the write operation
   integer, intent(out), optional :: iostat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out), optional :: iomsg

   integer :: stat
   character(len=:), allocatable :: msg
#if TBLITE_HAS_HDF5
   integer(hsize_t) :: dims(1)
   dims = shape(array, kind=hsize_t)
   call write_rdp_r1(filename, varname, array, dims, stat, msg)
#else
   stat = -1
   msg = "HDF5 support is not available in this build of tblite."
#endif
   call handle_iostat(stat, msg, filename//"::"//varname, iostat, iomsg, "write")
end subroutine save_hdf5_rdp_r1

subroutine save_hdf5_rdp_r2(filename, varname, array, iostat, iomsg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> Name of the dataset to save
   character(len=*), intent(in) :: varname
   !> Array to save
   real(dp), intent(in) :: array(:, :)
   !> Status of the write operation
   integer, intent(out), optional :: iostat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out), optional :: iomsg

   integer :: stat
   character(len=:), allocatable :: msg
#if TBLITE_HAS_HDF5
   integer(hsize_t) :: dims(2)
   dims = shape(array, kind=hsize_t)
   call write_rdp_r2(filename, varname, array, dims, stat, msg)
#else
   stat = -1
   msg = "HDF5 support is not available in this build of tblite."
#endif
   call handle_iostat(stat, msg, filename//"::"//varname, iostat, iomsg, "write")
end subroutine save_hdf5_rdp_r2

subroutine save_hdf5_rdp_r3(filename, varname, array, iostat, iomsg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> Name of the dataset to save
   character(len=*), intent(in) :: varname
   !> Array to save
   real(dp), intent(in) :: array(:, :, :)
   !> Status of the write operation
   integer, intent(out), optional :: iostat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out), optional :: iomsg

   integer :: stat
   character(len=:), allocatable :: msg
#if TBLITE_HAS_HDF5
   integer(hsize_t) :: dims(3)
   dims = shape(array, kind=hsize_t)
   call write_rdp_r3(filename, varname, array, dims, stat, msg)
#else
   stat = -1
   msg = "HDF5 support is not available in this build of tblite."
#endif
   call handle_iostat(stat, msg, filename//"::"//varname, iostat, iomsg, "write")
end subroutine save_hdf5_rdp_r3

#if TBLITE_HAS_HDF5
subroutine hdf5_open(stat, msg)
   !> Status of HDF5 library initialization
   integer, intent(out) :: stat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out) :: msg

   call h5open_f(stat)
   if (stat /= 0) msg = "Failed to initialize HDF5 library"
end subroutine hdf5_open

subroutine hdf5_close()
   integer :: stat

   call h5close_f(stat)
end subroutine hdf5_close

subroutine close_file(file_id)
   !> HDF5 file identifier to close
   integer(hid_t), intent(in) :: file_id

   integer :: stat

   if (file_id > 0) call h5fclose_f(file_id, stat)
end subroutine close_file

subroutine close_dataset(dset_id)
   !> HDF5 dataset identifier to close
   integer(hid_t), intent(in) :: dset_id

   integer :: stat

   if (dset_id > 0) call h5dclose_f(dset_id, stat)
end subroutine close_dataset

subroutine close_space(space_id)
   !> HDF5 dataspace identifier to close
   integer(hid_t), intent(in) :: space_id

   integer :: stat

   if (space_id > 0) call h5sclose_f(space_id, stat)
end subroutine close_space

subroutine open_file_read(filename, file_id, stat, msg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> HDF5 file identifier
   integer(hid_t), intent(out) :: file_id
   !> Status of the open operation
   integer, intent(out) :: stat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out) :: msg

   file_id = -1_hid_t
   call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, stat)
   if (stat /= 0) msg = "Failed to open HDF5 file '"//filename//"' for reading"
end subroutine open_file_read

subroutine open_file_write(filename, file_id, stat, msg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> HDF5 file identifier
   integer(hid_t), intent(out) :: file_id
   !> Status of the open or create operation
   integer, intent(out) :: stat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out) :: msg

   logical :: exist

   file_id = -1_hid_t
   inquire(file=filename, exist=exist)
   if (exist) then
      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, stat)
      if (stat /= 0) msg = "Failed to open HDF5 file '"//filename//"' for writing"
   else
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, stat)
      if (stat /= 0) msg = "Failed to create HDF5 file '"//filename//"'"
   end if
end subroutine open_file_write

subroutine create_dataset(file_id, varname, type_id, dims, dset_id, stat, msg)
   !> HDF5 file identifier
   integer(hid_t), intent(in) :: file_id
   !> Name of the dataset to create
   character(len=*), intent(in) :: varname
   !> HDF5 datatype identifier for the dataset
   integer(hid_t), intent(in) :: type_id
   !> Dataset dimensions
   integer(hsize_t), intent(in) :: dims(:)
   !> HDF5 dataset identifier
   integer(hid_t), intent(out) :: dset_id
   !> Status of the create operation
   integer, intent(out) :: stat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out) :: msg

   logical :: exists
   integer(hid_t) :: space_id

   dset_id = -1_hid_t
   space_id = -1_hid_t
   call h5lexists_f(file_id, varname, exists, stat)
   if (stat == 0 .and. exists) call h5ldelete_f(file_id, varname, stat)
   if (stat == 0) call h5screate_simple_f(size(dims), dims, space_id, stat)
   if (stat == 0) call h5dcreate_f(file_id, varname, type_id, space_id, dset_id, stat)
   if (stat /= 0) msg = "Failed to create HDF5 dataset '"//varname//"'"
   call close_space(space_id)
end subroutine create_dataset

subroutine read_shape(filename, varname, dims, stat, msg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> Name of the dataset to inspect
   character(len=*), intent(in) :: varname
   !> Dataset dimensions
   integer(hsize_t), intent(out) :: dims(:)
   !> Status of the read operation
   integer, intent(out) :: stat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out) :: msg

   integer :: rank
   integer(hsize_t) :: maxdims(size(dims))
   integer(hid_t) :: file_id, dset_id, space_id

   file_id = -1_hid_t
   dset_id = -1_hid_t
   space_id = -1_hid_t
   call hdf5_open(stat, msg)
   if (stat == 0) call open_file_read(filename, file_id, stat, msg)
   if (stat == 0) call h5dopen_f(file_id, varname, dset_id, stat)
   if (stat /= 0 .and. .not.allocated(msg)) msg = "Dataset '"//filename//"::"//varname//"' not found"
   if (stat == 0) call h5dget_space_f(dset_id, space_id, stat)
   if (stat == 0) call h5sget_simple_extent_ndims_f(space_id, rank, stat)
   if (stat == 0 .and. rank /= size(dims)) then
      stat = 502
      msg = "Unexpected rank for HDF5 dataset '"//filename//"::"//varname//"'"
   end if
   if (stat == 0) then
      call h5sget_simple_extent_dims_f(space_id, dims, maxdims, stat)
      if (stat >= 0) stat = 0
   end if
   call close_space(space_id)
   call close_dataset(dset_id)
   call close_file(file_id)
   call hdf5_close()
end subroutine read_shape

subroutine write_i4_r1(filename, varname, array, dims, stat, msg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> Name of the dataset to write
   character(len=*), intent(in) :: varname
   !> Array to write
   integer(i4), intent(in) :: array(:)
   !> Dataset dimensions
   integer(hsize_t), intent(in) :: dims(:)
   !> Status of the write operation
   integer, intent(out) :: stat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out) :: msg
   integer(hid_t) :: file_id, dset_id

   call hdf5_open(stat, msg)
   if (stat == 0) call open_file_write(filename, file_id, stat, msg)
   if (stat == 0) call create_dataset(file_id, varname, H5T_NATIVE_INTEGER, dims, dset_id, stat, msg)
   if (stat == 0) call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dims, stat)
   if (stat /= 0 .and. .not.allocated(msg)) msg = "Failed to write HDF5 dataset '"//filename//"::"//varname//"'"
   call close_dataset(dset_id)
   call close_file(file_id)
   call hdf5_close()
end subroutine write_i4_r1

subroutine write_rdp_r1(filename, varname, array, dims, stat, msg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> Name of the dataset to write
   character(len=*), intent(in) :: varname
   !> Array to write
   real(dp), intent(in) :: array(:)
   !> Dataset dimensions
   integer(hsize_t), intent(in) :: dims(:)
   !> Status of the write operation
   integer, intent(out) :: stat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out) :: msg
   integer(hid_t) :: file_id, dset_id

   call hdf5_open(stat, msg)
   if (stat == 0) call open_file_write(filename, file_id, stat, msg)
   if (stat == 0) call create_dataset(file_id, varname, H5T_NATIVE_DOUBLE, dims, dset_id, stat, msg)
   if (stat == 0) call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, stat)
   if (stat /= 0 .and. .not.allocated(msg)) msg = "Failed to write HDF5 dataset '"//filename//"::"//varname//"'"
   call close_dataset(dset_id)
   call close_file(file_id)
   call hdf5_close()
end subroutine write_rdp_r1

subroutine write_rdp_r2(filename, varname, array, dims, stat, msg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> Name of the dataset to write
   character(len=*), intent(in) :: varname
   !> Array to write
   real(dp), intent(in) :: array(:, :)
   !> Dataset dimensions
   integer(hsize_t), intent(in) :: dims(:)
   !> Status of the write operation
   integer, intent(out) :: stat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out) :: msg
   integer(hid_t) :: file_id, dset_id

   call hdf5_open(stat, msg)
   if (stat == 0) call open_file_write(filename, file_id, stat, msg)
   if (stat == 0) call create_dataset(file_id, varname, H5T_NATIVE_DOUBLE, dims, dset_id, stat, msg)
   if (stat == 0) call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, stat)
   if (stat /= 0 .and. .not.allocated(msg)) msg = "Failed to write HDF5 dataset '"//filename//"::"//varname//"'"
   call close_dataset(dset_id)
   call close_file(file_id)
   call hdf5_close()
end subroutine write_rdp_r2

subroutine write_rdp_r3(filename, varname, array, dims, stat, msg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> Name of the dataset to write
   character(len=*), intent(in) :: varname
   !> Array to write
   real(dp), intent(in) :: array(:, :, :)
   !> Dataset dimensions
   integer(hsize_t), intent(in) :: dims(:)
   !> Status of the write operation
   integer, intent(out) :: stat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out) :: msg
   integer(hid_t) :: file_id, dset_id

   call hdf5_open(stat, msg)
   if (stat == 0) call open_file_write(filename, file_id, stat, msg)
   if (stat == 0) call create_dataset(file_id, varname, H5T_NATIVE_DOUBLE, dims, dset_id, stat, msg)
   if (stat == 0) call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, stat)
   if (stat /= 0 .and. .not.allocated(msg)) msg = "Failed to write HDF5 dataset '"//filename//"::"//varname//"'"
   call close_dataset(dset_id)
   call close_file(file_id)
   call hdf5_close()
end subroutine write_rdp_r3

subroutine read_i4_r1(filename, varname, array, dims, stat, msg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> Name of the dataset to read
   character(len=*), intent(in) :: varname
   !> Array to read into
   integer(i4), intent(out) :: array(:)
   !> Dataset dimensions
   integer(hsize_t), intent(in) :: dims(:)
   !> Status of the read operation
   integer, intent(out) :: stat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out) :: msg
   integer(hid_t) :: file_id, dset_id

   call hdf5_open(stat, msg)
   if (stat == 0) call open_file_read(filename, file_id, stat, msg)
   if (stat == 0) call h5dopen_f(file_id, varname, dset_id, stat)
   if (stat == 0) call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dims, stat)
   if (stat /= 0 .and. .not.allocated(msg)) msg = "Failed to read HDF5 dataset '"//filename//"::"//varname//"'"
   call close_dataset(dset_id)
   call close_file(file_id)
   call hdf5_close()
end subroutine read_i4_r1

subroutine read_rdp_r1(filename, varname, array, dims, stat, msg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> Name of the dataset to read
   character(len=*), intent(in) :: varname
   !> Array to read into
   real(dp), intent(out) :: array(:)
   !> Dataset dimensions
   integer(hsize_t), intent(in) :: dims(:)
   !> Status of the read operation
   integer, intent(out) :: stat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out) :: msg
   integer(hid_t) :: file_id, dset_id

   call hdf5_open(stat, msg)
   if (stat == 0) call open_file_read(filename, file_id, stat, msg)
   if (stat == 0) call h5dopen_f(file_id, varname, dset_id, stat)
   if (stat == 0) call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, stat)
   if (stat /= 0 .and. .not.allocated(msg)) msg = "Failed to read HDF5 dataset '"//filename//"::"//varname//"'"
   call close_dataset(dset_id)
   call close_file(file_id)
   call hdf5_close()
end subroutine read_rdp_r1

subroutine read_rdp_r2(filename, varname, array, dims, stat, msg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> Name of the dataset to read
   character(len=*), intent(in) :: varname
   !> Array to read into
   real(dp), intent(out) :: array(:, :)
   !> Dataset dimensions
   integer(hsize_t), intent(in) :: dims(:)
   !> Status of the read operation
   integer, intent(out) :: stat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out) :: msg
   integer(hid_t) :: file_id, dset_id

   call hdf5_open(stat, msg)
   if (stat == 0) call open_file_read(filename, file_id, stat, msg)
   if (stat == 0) call h5dopen_f(file_id, varname, dset_id, stat)
   if (stat == 0) call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, stat)
   if (stat /= 0 .and. .not.allocated(msg)) msg = "Failed to read HDF5 dataset '"//filename//"::"//varname//"'"
   call close_dataset(dset_id)
   call close_file(file_id)
   call hdf5_close()
end subroutine read_rdp_r2

subroutine read_rdp_r3(filename, varname, array, dims, stat, msg)
   !> Filename of the HDF5 file
   character(len=*), intent(in) :: filename
   !> Name of the dataset to read
   character(len=*), intent(in) :: varname
   !> Array to read into
   real(dp), intent(out) :: array(:, :, :)
   !> Dataset dimensions
   integer(hsize_t), intent(in) :: dims(:)
   !> Status of the read operation
   integer, intent(out) :: stat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out) :: msg
   integer(hid_t) :: file_id, dset_id

   call hdf5_open(stat, msg)
   if (stat == 0) call open_file_read(filename, file_id, stat, msg)
   if (stat == 0) call h5dopen_f(file_id, varname, dset_id, stat)
   if (stat == 0) call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, stat)
   if (stat /= 0 .and. .not.allocated(msg)) msg = "Failed to read HDF5 dataset '"//filename//"::"//varname//"'"
   call close_dataset(dset_id)
   call close_file(file_id)
   call hdf5_close()
end subroutine read_rdp_r3
#endif

subroutine handle_iostat(stat, msg, filename, iostat, iomsg, action)
   !> Status of the HDF5 operation
   integer, intent(in) :: stat
   !> Error message from the HDF5 operation
   character(len=:), allocatable, intent(in) :: msg
   !> Filename and dataset used in diagnostics
   character(len=*), intent(in) :: filename
   !> Status returned to the caller
   integer, intent(out), optional :: iostat
   !> Error message returned to the caller
   character(len=:), allocatable, intent(out), optional :: iomsg
   !> Operation name used in diagnostics
   character(len=*), intent(in) :: action

   if (present(iostat)) then
      iostat = stat
   else if (stat /= 0) then
      if (allocated(msg)) then
         error stop "Failed to "//action//" array from file '"//filename//"'"//nl//msg
      else
         error stop "Failed to "//action//" array from file '"//filename//"'"
      end if
   end if
   if (present(iomsg) .and. allocated(msg) .and. stat /= 0) iomsg = msg
end subroutine handle_iostat

end module tblite_io_hdf5
