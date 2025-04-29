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

!> @file tblite/io/numpy/savez.f90
!> Provides npz output routines

!> Implementation of saving multidimensional arrays to npz files
module tblite_io_numpy_savez
   use iso_fortran_env, only : iostat_end
   use mctc_env, only : sp, dp, i2, i4, i1
   use tblite_io_numpy_constants, only : type_rdp, type_i4, &
      & zip_global_sig, zip_local_sig, zip_footer_sig, zip_min_version
   use tblite_io_numpy_crc32, only : crc32_hash
   use tblite_io_numpy_save, only : npy_header
   use tblite_output_format, only : format_string
   implicit none
   private

   public :: save_npz

   interface save_npz
      module procedure save_npz_i4_r1
      module procedure save_npz_rdp_r1
      module procedure save_npz_rdp_r2
      module procedure save_npz_rdp_r3
   end interface save_npz

contains

subroutine save_npz_i4_r1(filename, varname, array, iostat, iomsg)
   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: varname
   integer(i4), intent(in) :: array(:)
   integer, intent(out), optional :: iostat
   character(len=:), allocatable, intent(out), optional :: iomsg

   character(len=*), parameter :: vtype = type_i4
   integer, parameter :: vsize = 4

   integer(i2) :: nrecs
   integer(i4) :: global_header_offset, checksum
   character(len=:), allocatable :: file_header, global_header, local_header, footer, path
   integer(i4) :: nbytes

   logical :: exist
   integer :: io, stat
   character(len=:), allocatable :: msg

   path = varname // ".npy"
   nrecs = 0_i2
   global_header_offset = 0_i4
   global_header = ""

   inquire(file=filename, exist=exist)
   open(newunit=io, file=filename, form="unformatted", access="stream", iostat=stat)
   if (stat == 0 .and. exist) then
      call parse_zip_footer(io, filename, nrecs, global_header, global_header_offset, stat, msg)
   end if

   if (stat == 0) then
      file_header = npy_header(vtype, shape(array))

      checksum = crc32_hash(file_header)
      checksum = crc32_hash(array, checksum)

      nbytes = len(file_header) + size(array) * vsize

      local_header = get_local_header(path, checksum, nbytes)
      global_header = global_header // get_global_header(path, local_header, global_header_offset)
      footer = get_footer(nrecs + 1_i2, len(global_header), global_header_offset + len(local_header) + nbytes)

      write(io, pos=global_header_offset+1, iostat=stat) &
         & local_header, &
         & file_header, &
         & array, &
         & global_header, &
         & footer
   end if
   close(io)

   call handle_iostat(stat, filename, iostat, iomsg)
   if (allocated(msg) .and. present(iomsg)) call move_alloc(msg, iomsg)
end subroutine save_npz_i4_r1

subroutine save_npz_rdp_r1(filename, varname, array, iostat, iomsg)
   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: varname
   real(dp), intent(in) :: array(:)
   integer, intent(out), optional :: iostat
   character(len=:), allocatable, intent(out), optional :: iomsg

   character(len=*), parameter :: vtype = type_rdp
   integer, parameter :: vsize = 8

   integer(i2) :: nrecs
   integer(i4) :: global_header_offset, checksum
   character(len=:), allocatable :: file_header, global_header, local_header, footer, path
   integer(i4) :: nbytes

   logical :: exist
   integer :: io, stat
   character(len=:), allocatable :: msg

   path = varname // ".npy"
   nrecs = 0_i2
   global_header_offset = 0_i4
   global_header = ""

   inquire(file=filename, exist=exist)
   open(newunit=io, file=filename, form="unformatted", access="stream", iostat=stat)
   if (stat == 0 .and. exist) then
      call parse_zip_footer(io, filename, nrecs, global_header, global_header_offset, stat, msg)
   end if

   if (stat == 0) then
      file_header = npy_header(vtype, shape(array))

      checksum = crc32_hash(file_header)
      checksum = crc32_hash(array, checksum)

      nbytes = len(file_header) + size(array) * vsize

      local_header = get_local_header(path, checksum, nbytes)
      global_header = global_header // get_global_header(path, local_header, global_header_offset)
      footer = get_footer(nrecs + 1_i2, len(global_header), global_header_offset + len(local_header) + nbytes)

      write(io, pos=global_header_offset+1, iostat=stat) &
         & local_header, &
         & file_header, &
         & array, &
         & global_header, &
         & footer
   end if
   close(io)

   call handle_iostat(stat, filename, iostat, iomsg)
   if (allocated(msg) .and. present(iomsg)) call move_alloc(msg, iomsg)
end subroutine save_npz_rdp_r1

subroutine save_npz_rdp_r2(filename, varname, array, iostat, iomsg)
   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: varname
   real(dp), intent(in) :: array(:, :)
   integer, intent(out), optional :: iostat
   character(len=:), allocatable, intent(out), optional :: iomsg

   character(len=*), parameter :: vtype = type_rdp
   integer, parameter :: vsize = 8

   integer(i2) :: nrecs
   integer(i4) :: global_header_offset, checksum
   character(len=:), allocatable :: file_header, global_header, local_header, footer, path
   integer(i4) :: nbytes

   logical :: exist
   integer :: io, stat
   character(len=:), allocatable :: msg

   path = varname // ".npy"
   nrecs = 0_i2
   global_header_offset = 0_i4
   global_header = ""

   inquire(file=filename, exist=exist)
   open(newunit=io, file=filename, form="unformatted", access="stream", iostat=stat)
   if (stat == 0 .and. exist) then
      call parse_zip_footer(io, filename, nrecs, global_header, global_header_offset, stat, msg)
   end if

   if (stat == 0) then
      file_header = npy_header(vtype, shape(array))

      checksum = crc32_hash(file_header)
      checksum = crc32_hash(reshape(array, [size(array)]), checksum)

      nbytes = len(file_header) + size(array) * vsize

      local_header = get_local_header(path, checksum, nbytes)
      global_header = global_header // get_global_header(path, local_header, global_header_offset)
      footer = get_footer(nrecs + 1_i2, len(global_header), global_header_offset + len(local_header) + nbytes)

      write(io, pos=global_header_offset+1, iostat=stat) &
         & local_header, &
         & file_header, &
         & array, &
         & global_header, &
         & footer
   end if
   close(io)

   call handle_iostat(stat, filename, iostat, iomsg)
   if (allocated(msg) .and. present(iomsg)) call move_alloc(msg, iomsg)
end subroutine save_npz_rdp_r2

subroutine save_npz_rdp_r3(filename, varname, array, iostat, iomsg)
   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: varname
   real(dp), intent(in) :: array(:, :, :)
   integer, intent(out), optional :: iostat
   character(len=:), allocatable, intent(out), optional :: iomsg

   character(len=*), parameter :: vtype = type_rdp
   integer, parameter :: vsize = 8

   integer(i2) :: nrecs
   integer(i4) :: global_header_offset, checksum
   character(len=:), allocatable :: file_header, global_header, local_header, footer, path
   integer(i4) :: nbytes

   logical :: exist
   integer :: io, stat
   character(len=:), allocatable :: msg

   path = varname // ".npy"
   nrecs = 0_i2
   global_header_offset = 0_i4
   global_header = ""

   inquire(file=filename, exist=exist)
   open(newunit=io, file=filename, form="unformatted", access="stream", iostat=stat)
   if (stat == 0 .and. exist) then
      call parse_zip_footer(io, filename, nrecs, global_header, global_header_offset, stat, msg)
   end if

   if (stat == 0) then
      file_header = npy_header(vtype, shape(array))

      checksum = crc32_hash(file_header)
      checksum = crc32_hash(reshape(array, [size(array)]), checksum)

      nbytes = len(file_header) + size(array) * vsize

      local_header = get_local_header(path, checksum, nbytes)
      global_header = global_header // get_global_header(path, local_header, global_header_offset)
      footer = get_footer(nrecs + 1_i2, len(global_header), global_header_offset + len(local_header) + nbytes)

      write(io, pos=global_header_offset+1, iostat=stat) &
         & local_header, &
         & file_header, &
         & array, &
         & global_header, &
         & footer
   end if
   close(io)

   call handle_iostat(stat, filename, iostat, iomsg)
   if (allocated(msg) .and. present(iomsg)) call move_alloc(msg, iomsg)
end subroutine save_npz_rdp_r3


!> Create a local header for a zip file
pure function get_local_header(filename, checksum, nbytes) result(header)
   !> Name of the file to be saved
   character(len=*), intent(in) :: filename
   !> CRC32 checksum of the file
   integer(i4), intent(in) :: checksum
   !> Size of the file in bytes
   integer(i4), intent(in) :: nbytes

   character(len=2), parameter :: general_purpose_bit_flag = repeat(char(0), 2)
   character(len=2), parameter :: compression_method = repeat(char(0), 2)
   character(len=2), parameter :: file_last_mod_time = repeat(char(0), 2)
   character(len=2), parameter :: file_last_mod_date = repeat(char(0), 2)
   character(len=2), parameter :: extra_field_length = repeat(char(0), 2)

   character(len=30+len(filename)) :: header

   header(1:4) = transfer(zip_local_sig, header(1:4))
   header(5:6) = transfer(zip_min_version, header(5:6))
   header(7:8) = general_purpose_bit_flag
   header(9:10) = compression_method
   header(11:12) = file_last_mod_time
   header(13:14) = file_last_mod_date
   header(15:18) = transfer(checksum, header(15:18))
   header(19:22) = transfer(nbytes, header(19:22))
   header(23:26) = header(19:22)
   header(27:28) = transfer(len(filename), header(25:28))
   header(29:30) = extra_field_length
   header(31:30+len(filename)) = filename

end function get_local_header

pure function get_global_header(filename, local_header, global_header_offset) result(header)
   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: local_header
   integer(i4), intent(in) :: global_header_offset

   character(len=2), parameter :: file_comment_length = repeat(char(0), 2)
   character(len=2), parameter :: disk_no = repeat(char(0), 2)
   character(len=2), parameter :: internal_file_attributes = repeat(char(0), 2)
   character(len=4), parameter :: external_file_attributes = repeat(char(0), 4)

   character(len=46+len(filename)) :: header

   header(1:4) = transfer(zip_global_sig, header(1:4))
   header(5:6) = transfer(zip_min_version, header(5:6))
   header(7:32) = local_header(5:30)
   header(33:34) = file_comment_length
   header(35:36) = disk_no
   header(37:38) = internal_file_attributes
   header(39:42) = external_file_attributes
   header(43:46) = transfer(global_header_offset, header(43:46))
   header(47:46+len(filename)) = filename

end function get_global_header

pure function get_footer(nrecs, global_header_size, global_header_offset) result(footer)
   integer(i2), intent(in) :: nrecs
   integer(i4), intent(in) :: global_header_size
   integer(i4), intent(in) :: global_header_offset

   character(len=2), parameter :: disk_no = repeat(char(0), 2)
   character(len=2), parameter :: disk_start = repeat(char(0), 2)
   character(len=2), parameter :: comment_len = repeat(char(0), 2)

   character(len=22) :: footer

   footer(1:4) = transfer(zip_footer_sig, footer(1:4))
   footer(5:6) = disk_no
   footer(7:8) = disk_start
   footer(9:10) = transfer(nrecs, footer(9:10))
   footer(11:12) = footer(9:10)
   footer(13:16) = transfer(global_header_size, footer(13:16))
   footer(17:20) = transfer(global_header_offset, footer(17:20))
   footer(21:22) = comment_len

end function get_footer

subroutine parse_zip_footer(io, filename, nrecs, global_header, global_header_offset, stat, msg)
   integer, intent(in) :: io
   character(len=*), intent(in) :: filename
   integer(i2), intent(out) :: nrecs
   character(len=:), allocatable, intent(out) :: global_header
   integer(i4), intent(out) :: global_header_offset
   integer, intent(out) :: stat
   character(len=:), allocatable :: msg

   integer(i2) :: path_size, extra_field_size, comment_size
   integer(i2) :: disk_no, disk_start, nrecs_on_disk
   integer(i4) :: nbytes_compressed, global_header_size
   character(len=512) :: errmsg
   integer :: res, length, pos
   integer(i4) :: header_sig
   character(len=:), allocatable :: path

   stat = 0
   pos = 1
   read(io, pos=pos, iostat=stat, iomsg=errmsg) header_sig
   do while(stat == 0 .and. is_local_header(header_sig))

      if (stat == 0) &
         read(io, pos=pos+18, iostat=stat, iomsg=errmsg) nbytes_compressed
      if (stat == 0) &
         read(io, pos=pos+26, iostat=stat, iomsg=errmsg) path_size
      if (stat == 0) &
         read(io, pos=pos+28, iostat=stat, iomsg=errmsg) extra_field_size

      if (stat == 0) then
         if (allocated(path)) deallocate(path)
         allocate(character(len=path_size) :: path, stat=stat)
      end if
      if (stat == 0) &
         read(io, pos=pos+30, iostat=stat, iomsg=errmsg) path

      pos = pos + 30 + path_size + extra_field_size + nbytes_compressed
      read(io, pos=pos, iostat=stat, iomsg=errmsg) header_sig
   end do
   if (stat /= 0) then
      msg = "Failed to read local header block for '"//filename//"'"
      if (len_trim(errmsg) > 0) &
         msg = msg // " ("//trim(errmsg)//")"
      return
   end if

   if (.not.is_global_header(header_sig)) then
      stat = 400
      msg = "Expected global header signature for '"//filename//"' got "// &
         & format_string(header_sig, '(z0.8)')
      return
   end if

   ! global_header_offset = pos - 1
   do while(stat == 0 .and. is_global_header(header_sig))

      if (stat == 0) &
         read(io, pos=pos+28, iostat=stat, iomsg=errmsg) path_size
      if (stat == 0) &
         read(io, pos=pos+30, iostat=stat, iomsg=errmsg) extra_field_size
      if (stat == 0) &
         read(io, pos=pos+32, iostat=stat, iomsg=errmsg) comment_size

      if (stat == 0) then
         if (allocated(path)) deallocate(path)
         allocate(character(len=path_size) :: path, stat=stat)
      end if
      if (stat == 0) &
         read(io, pos=pos+46, iostat=stat, iomsg=errmsg) path

      pos = pos + 46 + path_size + extra_field_size + comment_size
      read(io, pos=pos, iostat=stat, iomsg=errmsg) header_sig
   end do
   if (stat /= 0) then
      msg = "Failed to read global header block for '"//filename//"'"
      if (len_trim(errmsg) > 0) &
         msg = msg // " ("//trim(errmsg)//")"
      return
   end if
   if (.not.is_footer_header(header_sig)) then
      stat = 401
      msg = "Expected footer signature for '"//filename//"' got "// &
         & format_string(header_sig, '(z0.8)')
      return
   end if
   ! global_header_size = pos - global_header_offset + 1

   read(io, pos=pos+4, iostat=stat, iomsg=errmsg) &
      & disk_no, disk_start, nrecs_on_disk, nrecs, &
      & global_header_size, global_header_offset, comment_size

   if (stat == 0) &
      allocate(character(len=global_header_size) :: global_header, stat=stat)
   if (stat == 0) &
      read(io, iostat=stat, pos=global_header_offset+1) global_header

   if (stat /= 0) then
      msg = "Failed to read footer for '"//filename//"'"
      if (len_trim(errmsg) > 0) &
         msg = msg // " ("//trim(errmsg)//")"
      return
   end if

   if (disk_no /= 0) then
      stat = 402
      msg = "Cannot read zip file with disk_no != 0"
   end if

   if (disk_start /= 0) then
      stat = 403
      msg = "Cannot read zip file with disk_start != 0"
   end if

   if (nrecs_on_disk /= nrecs) then
      stat = 404
      msg = "Cannot read zip file with nrecs_on_disk != nrecs"
   end if
end subroutine parse_zip_footer

pure function is_local_header(header_sig) result(is_local)
   integer(i4), intent(in) :: header_sig
   logical :: is_local

   is_local = header_sig == zip_local_sig
end function is_local_header

pure function is_global_header(header_sig) result(is_global)
   integer(i4), intent(in) :: header_sig
   logical :: is_global

   is_global = header_sig == zip_global_sig
end function is_global_header

pure function is_footer_header(header_sig) result(is_footer)
   integer(i4), intent(in) :: header_sig
   logical :: is_footer

   is_footer = header_sig == zip_footer_sig
end function is_footer_header

!> Handle iostat and iomsg of write operation
subroutine handle_iostat(stat, filename, iostat, iomsg)
   !> Error status of loading, zero on success
   integer, intent(in) :: stat
   !> Name of the npy file to load from
   character(len=*), intent(in) :: filename
   !> Error status of loading, zero on success
   integer, intent(out), optional :: iostat
   !> Associated error message in case of non-zero status code
   character(len=:), allocatable, intent(out), optional :: iomsg

   if (present(iostat)) then
      iostat = stat
   else if (stat /= 0) then
      error stop "Failed to write array to file '"//filename//"'"
   end if

   if (present(iomsg)) then
      if (stat /= 0) then
         iomsg = "Failed to write array to file '"//filename//"'"
      end if
   end if
end subroutine handle_iostat

end module tblite_io_numpy_savez
