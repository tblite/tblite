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

!> @file tblite/io/numpy/load.f90
!> Provides npy input routines

!> Implementation of cyclic redundancy check hashing function,
!> used to check the integrity of data in zip files.
module tblite_io_numpy_crc32
   use mctc_env, only : i4, dp
   implicit none
   private

   public :: crc32_hash

   !> Compute crc32 checksum
   interface crc32_hash
      module procedure crc32_hash_char_r0
      module procedure crc32_hash_char_r1
      module procedure crc32_hash_i4_r1
      module procedure crc32_hash_rdp_r1
   end interface crc32_hash

   integer(i4), parameter :: crc_table(0:255) = (/ &
      & int(z'00000000', kind=i4), int(z'77073096', kind=i4), int(z'ee0e612c', kind=i4), int(z'990951ba', kind=i4), &
      & int(z'076dc419', kind=i4), int(z'706af48f', kind=i4), int(z'e963a535', kind=i4), int(z'9e6495a3', kind=i4), &
      & int(z'0edb8832', kind=i4), int(z'79dcb8a4', kind=i4), int(z'e0d5e91e', kind=i4), int(z'97d2d988', kind=i4), &
      & int(z'09b64c2b', kind=i4), int(z'7eb17cbd', kind=i4), int(z'e7b82d07', kind=i4), int(z'90bf1d91', kind=i4), &
      & int(z'1db71064', kind=i4), int(z'6ab020f2', kind=i4), int(z'f3b97148', kind=i4), int(z'84be41de', kind=i4), &
      & int(z'1adad47d', kind=i4), int(z'6ddde4eb', kind=i4), int(z'f4d4b551', kind=i4), int(z'83d385c7', kind=i4), &
      & int(z'136c9856', kind=i4), int(z'646ba8c0', kind=i4), int(z'fd62f97a', kind=i4), int(z'8a65c9ec', kind=i4), &
      & int(z'14015c4f', kind=i4), int(z'63066cd9', kind=i4), int(z'fa0f3d63', kind=i4), int(z'8d080df5', kind=i4), &
      & int(z'3b6e20c8', kind=i4), int(z'4c69105e', kind=i4), int(z'd56041e4', kind=i4), int(z'a2677172', kind=i4), &
      & int(z'3c03e4d1', kind=i4), int(z'4b04d447', kind=i4), int(z'd20d85fd', kind=i4), int(z'a50ab56b', kind=i4), &
      & int(z'35b5a8fa', kind=i4), int(z'42b2986c', kind=i4), int(z'dbbbc9d6', kind=i4), int(z'acbcf940', kind=i4), &
      & int(z'32d86ce3', kind=i4), int(z'45df5c75', kind=i4), int(z'dcd60dcf', kind=i4), int(z'abd13d59', kind=i4), &
      & int(z'26d930ac', kind=i4), int(z'51de003a', kind=i4), int(z'c8d75180', kind=i4), int(z'bfd06116', kind=i4), &
      & int(z'21b4f4b5', kind=i4), int(z'56b3c423', kind=i4), int(z'cfba9599', kind=i4), int(z'b8bda50f', kind=i4), &
      & int(z'2802b89e', kind=i4), int(z'5f058808', kind=i4), int(z'c60cd9b2', kind=i4), int(z'b10be924', kind=i4), &
      & int(z'2f6f7c87', kind=i4), int(z'58684c11', kind=i4), int(z'c1611dab', kind=i4), int(z'b6662d3d', kind=i4), &
      & int(z'76dc4190', kind=i4), int(z'01db7106', kind=i4), int(z'98d220bc', kind=i4), int(z'efd5102a', kind=i4), &
      & int(z'71b18589', kind=i4), int(z'06b6b51f', kind=i4), int(z'9fbfe4a5', kind=i4), int(z'e8b8d433', kind=i4), &
      & int(z'7807c9a2', kind=i4), int(z'0f00f934', kind=i4), int(z'9609a88e', kind=i4), int(z'e10e9818', kind=i4), &
      & int(z'7f6a0dbb', kind=i4), int(z'086d3d2d', kind=i4), int(z'91646c97', kind=i4), int(z'e6635c01', kind=i4), &
      & int(z'6b6b51f4', kind=i4), int(z'1c6c6162', kind=i4), int(z'856530d8', kind=i4), int(z'f262004e', kind=i4), &
      & int(z'6c0695ed', kind=i4), int(z'1b01a57b', kind=i4), int(z'8208f4c1', kind=i4), int(z'f50fc457', kind=i4), &
      & int(z'65b0d9c6', kind=i4), int(z'12b7e950', kind=i4), int(z'8bbeb8ea', kind=i4), int(z'fcb9887c', kind=i4), &
      & int(z'62dd1ddf', kind=i4), int(z'15da2d49', kind=i4), int(z'8cd37cf3', kind=i4), int(z'fbd44c65', kind=i4), &
      & int(z'4db26158', kind=i4), int(z'3ab551ce', kind=i4), int(z'a3bc0074', kind=i4), int(z'd4bb30e2', kind=i4), &
      & int(z'4adfa541', kind=i4), int(z'3dd895d7', kind=i4), int(z'a4d1c46d', kind=i4), int(z'd3d6f4fb', kind=i4), &
      & int(z'4369e96a', kind=i4), int(z'346ed9fc', kind=i4), int(z'ad678846', kind=i4), int(z'da60b8d0', kind=i4), &
      & int(z'44042d73', kind=i4), int(z'33031de5', kind=i4), int(z'aa0a4c5f', kind=i4), int(z'dd0d7cc9', kind=i4), &
      & int(z'5005713c', kind=i4), int(z'270241aa', kind=i4), int(z'be0b1010', kind=i4), int(z'c90c2086', kind=i4), &
      & int(z'5768b525', kind=i4), int(z'206f85b3', kind=i4), int(z'b966d409', kind=i4), int(z'ce61e49f', kind=i4), &
      & int(z'5edef90e', kind=i4), int(z'29d9c998', kind=i4), int(z'b0d09822', kind=i4), int(z'c7d7a8b4', kind=i4), &
      & int(z'59b33d17', kind=i4), int(z'2eb40d81', kind=i4), int(z'b7bd5c3b', kind=i4), int(z'c0ba6cad', kind=i4), &
      & int(z'edb88320', kind=i4), int(z'9abfb3b6', kind=i4), int(z'03b6e20c', kind=i4), int(z'74b1d29a', kind=i4), &
      & int(z'ead54739', kind=i4), int(z'9dd277af', kind=i4), int(z'04db2615', kind=i4), int(z'73dc1683', kind=i4), &
      & int(z'e3630b12', kind=i4), int(z'94643b84', kind=i4), int(z'0d6d6a3e', kind=i4), int(z'7a6a5aa8', kind=i4), &
      & int(z'e40ecf0b', kind=i4), int(z'9309ff9d', kind=i4), int(z'0a00ae27', kind=i4), int(z'7d079eb1', kind=i4), &
      & int(z'f00f9344', kind=i4), int(z'8708a3d2', kind=i4), int(z'1e01f268', kind=i4), int(z'6906c2fe', kind=i4), &
      & int(z'f762575d', kind=i4), int(z'806567cb', kind=i4), int(z'196c3671', kind=i4), int(z'6e6b06e7', kind=i4), &
      & int(z'fed41b76', kind=i4), int(z'89d32be0', kind=i4), int(z'10da7a5a', kind=i4), int(z'67dd4acc', kind=i4), &
      & int(z'f9b9df6f', kind=i4), int(z'8ebeeff9', kind=i4), int(z'17b7be43', kind=i4), int(z'60b08ed5', kind=i4), &
      & int(z'd6d6a3e8', kind=i4), int(z'a1d1937e', kind=i4), int(z'38d8c2c4', kind=i4), int(z'4fdff252', kind=i4), &
      & int(z'd1bb67f1', kind=i4), int(z'a6bc5767', kind=i4), int(z'3fb506dd', kind=i4), int(z'48b2364b', kind=i4), &
      & int(z'd80d2bda', kind=i4), int(z'af0a1b4c', kind=i4), int(z'36034af6', kind=i4), int(z'41047a60', kind=i4), &
      & int(z'df60efc3', kind=i4), int(z'a867df55', kind=i4), int(z'316e8eef', kind=i4), int(z'4669be79', kind=i4), &
      & int(z'cb61b38c', kind=i4), int(z'bc66831a', kind=i4), int(z'256fd2a0', kind=i4), int(z'5268e236', kind=i4), &
      & int(z'cc0c7795', kind=i4), int(z'bb0b4703', kind=i4), int(z'220216b9', kind=i4), int(z'5505262f', kind=i4), &
      & int(z'c5ba3bbe', kind=i4), int(z'b2bd0b28', kind=i4), int(z'2bb45a92', kind=i4), int(z'5cb36a04', kind=i4), &
      & int(z'c2d7ffa7', kind=i4), int(z'b5d0cf31', kind=i4), int(z'2cd99e8b', kind=i4), int(z'5bdeae1d', kind=i4), &
      & int(z'9b64c2b0', kind=i4), int(z'ec63f226', kind=i4), int(z'756aa39c', kind=i4), int(z'026d930a', kind=i4), &
      & int(z'9c0906a9', kind=i4), int(z'eb0e363f', kind=i4), int(z'72076785', kind=i4), int(z'05005713', kind=i4), &
      & int(z'95bf4a82', kind=i4), int(z'e2b87a14', kind=i4), int(z'7bb12bae', kind=i4), int(z'0cb61b38', kind=i4), &
      & int(z'92d28e9b', kind=i4), int(z'e5d5be0d', kind=i4), int(z'7cdcefb7', kind=i4), int(z'0bdbdf21', kind=i4), &
      & int(z'86d3d2d4', kind=i4), int(z'f1d4e242', kind=i4), int(z'68ddb3f8', kind=i4), int(z'1fda836e', kind=i4), &
      & int(z'81be16cd', kind=i4), int(z'f6b9265b', kind=i4), int(z'6fb077e1', kind=i4), int(z'18b74777', kind=i4), &
      & int(z'88085ae6', kind=i4), int(z'ff0f6a70', kind=i4), int(z'66063bca', kind=i4), int(z'11010b5c', kind=i4), &
      & int(z'8f659eff', kind=i4), int(z'f862ae69', kind=i4), int(z'616bffd3', kind=i4), int(z'166ccf45', kind=i4), &
      & int(z'a00ae278', kind=i4), int(z'd70dd2ee', kind=i4), int(z'4e048354', kind=i4), int(z'3903b3c2', kind=i4), &
      & int(z'a7672661', kind=i4), int(z'd06016f7', kind=i4), int(z'4969474d', kind=i4), int(z'3e6e77db', kind=i4), &
      & int(z'aed16a4a', kind=i4), int(z'd9d65adc', kind=i4), int(z'40df0b66', kind=i4), int(z'37d83bf0', kind=i4), &
      & int(z'a9bcae53', kind=i4), int(z'debb9ec5', kind=i4), int(z'47b2cf7f', kind=i4), int(z'30b5ffe9', kind=i4), &
      & int(z'bdbdf21c', kind=i4), int(z'cabac28a', kind=i4), int(z'53b39330', kind=i4), int(z'24b4a3a6', kind=i4), &
      & int(z'bad03605', kind=i4), int(z'cdd70693', kind=i4), int(z'54de5729', kind=i4), int(z'23d967bf', kind=i4), &
      & int(z'b3667a2e', kind=i4), int(z'c4614ab8', kind=i4), int(z'5d681b02', kind=i4), int(z'2a6f2b94', kind=i4), &
      & int(z'b40bbe37', kind=i4), int(z'c30c8ea1', kind=i4), int(z'5a05df1b', kind=i4), int(z'2d02ef8d', kind=i4)  &
      /)

contains

!> Compute crc32 checksum for a character string
pure function crc32_hash_char_r0(val, crc_in) result(crc)
   !> Previous crc32 checksum to continue from
   integer(i4), intent(in), optional :: crc_in
   !> Value to hash
   character(len=*), intent(in) :: val
   !> Resulting crc32 checksum
   integer(i4) :: crc
   integer :: ii

   if (present(crc_in)) then
      crc = crc_in
   else
      crc = 0_i4
   end if
   crc = not(crc)
   do ii = 1, len(val)
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(val(ii:ii))), 255)))
   enddo
   crc = not(crc)
end function crc32_hash_char_r0

!> Compute crc32 checksum for a character array
pure function crc32_hash_char_r1(val, crc_in) result(crc)
   !> Previous crc32 checksum to continue from
   integer(i4), intent(in), optional :: crc_in
   !> Value to hash
   character(len=1), intent(in) :: val(:)
   !> Resulting crc32 checksum
   integer(i4) :: crc
   integer :: ii

   if (present(crc_in)) then
      crc = crc_in
   else
      crc = 0_i4
   end if
   crc = not(crc)
   do ii = 1, size(val)
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(val(ii))), 255)))
   enddo
   crc = not(crc)
end function crc32_hash_char_r1

!> Compute crc32 checksum for a 4-byte integer array
pure function crc32_hash_i4_r1(val, crc_in) result(crc)
   !> Previous crc32 checksum to continue from
   integer(i4), intent(in), optional :: crc_in
   !> Value to hash
   integer(i4), intent(in) :: val(:)
   !> Resulting crc32 checksum
   integer(i4) :: crc
   integer :: ii

   character(len=1) :: chunk(4)

   if (present(crc_in)) then
      crc = crc_in
   else
      crc = 0_i4
   end if
   crc = not(crc)
   do ii = 1, size(val)
      chunk = transfer(val(ii), chunk)
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(1))), 255)))
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(2))), 255)))
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(3))), 255)))
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(4))), 255)))
   enddo
   crc = not(crc)
end function crc32_hash_i4_r1

!> Compute crc32 checksum for a real array
pure function crc32_hash_rdp_r1(val, crc_in) result(crc)
   !> Previous crc32 checksum to continue from
   integer(i4), intent(in), optional :: crc_in
   !> Value to hash
   real(dp), intent(in) :: val(:)
   !> Resulting crc32 checksum
   integer(i4) :: crc
   integer :: ii

   character(len=1) :: chunk(8)

   if (present(crc_in)) then
      crc = crc_in
   else
      crc = 0_i4
   end if
   crc = not(crc)
   do ii = 1, size(val)
      chunk = transfer(val(ii), chunk)
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(1))), 255)))
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(2))), 255)))
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(3))), 255)))
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(4))), 255)))
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(5))), 255)))
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(6))), 255)))
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(7))), 255)))
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(chunk(8))), 255)))
   enddo
   crc = not(crc)
end function crc32_hash_rdp_r1

end module tblite_io_numpy_crc32
