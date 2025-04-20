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

   integer(i4), parameter :: crc_stub(256) = 0
   integer(i4), parameter :: crc_table(0:255) = transfer([ &
      & real(z'00000000'), real(z'77073096'), real(z'ee0e612c'), real(z'990951ba'), &
      & real(z'076dc419'), real(z'706af48f'), real(z'e963a535'), real(z'9e6495a3'), &
      & real(z'0edb8832'), real(z'79dcb8a4'), real(z'e0d5e91e'), real(z'97d2d988'), &
      & real(z'09b64c2b'), real(z'7eb17cbd'), real(z'e7b82d07'), real(z'90bf1d91'), &
      & real(z'1db71064'), real(z'6ab020f2'), real(z'f3b97148'), real(z'84be41de'), &
      & real(z'1adad47d'), real(z'6ddde4eb'), real(z'f4d4b551'), real(z'83d385c7'), &
      & real(z'136c9856'), real(z'646ba8c0'), real(z'fd62f97a'), real(z'8a65c9ec'), &
      & real(z'14015c4f'), real(z'63066cd9'), real(z'fa0f3d63'), real(z'8d080df5'), &
      & real(z'3b6e20c8'), real(z'4c69105e'), real(z'd56041e4'), real(z'a2677172'), &
      & real(z'3c03e4d1'), real(z'4b04d447'), real(z'd20d85fd'), real(z'a50ab56b'), &
      & real(z'35b5a8fa'), real(z'42b2986c'), real(z'dbbbc9d6'), real(z'acbcf940'), &
      & real(z'32d86ce3'), real(z'45df5c75'), real(z'dcd60dcf'), real(z'abd13d59'), &
      & real(z'26d930ac'), real(z'51de003a'), real(z'c8d75180'), real(z'bfd06116'), &
      & real(z'21b4f4b5'), real(z'56b3c423'), real(z'cfba9599'), real(z'b8bda50f'), &
      & real(z'2802b89e'), real(z'5f058808'), real(z'c60cd9b2'), real(z'b10be924'), &
      & real(z'2f6f7c87'), real(z'58684c11'), real(z'c1611dab'), real(z'b6662d3d'), &
      & real(z'76dc4190'), real(z'01db7106'), real(z'98d220bc'), real(z'efd5102a'), &
      & real(z'71b18589'), real(z'06b6b51f'), real(z'9fbfe4a5'), real(z'e8b8d433'), &
      & real(z'7807c9a2'), real(z'0f00f934'), real(z'9609a88e'), real(z'e10e9818'), &
      & real(z'7f6a0dbb'), real(z'086d3d2d'), real(z'91646c97'), real(z'e6635c01'), &
      & real(z'6b6b51f4'), real(z'1c6c6162'), real(z'856530d8'), real(z'f262004e'), &
      & real(z'6c0695ed'), real(z'1b01a57b'), real(z'8208f4c1'), real(z'f50fc457'), &
      & real(z'65b0d9c6'), real(z'12b7e950'), real(z'8bbeb8ea'), real(z'fcb9887c'), &
      & real(z'62dd1ddf'), real(z'15da2d49'), real(z'8cd37cf3'), real(z'fbd44c65'), &
      & real(z'4db26158'), real(z'3ab551ce'), real(z'a3bc0074'), real(z'd4bb30e2'), &
      & real(z'4adfa541'), real(z'3dd895d7'), real(z'a4d1c46d'), real(z'd3d6f4fb'), &
      & real(z'4369e96a'), real(z'346ed9fc'), real(z'ad678846'), real(z'da60b8d0'), &
      & real(z'44042d73'), real(z'33031de5'), real(z'aa0a4c5f'), real(z'dd0d7cc9'), &
      & real(z'5005713c'), real(z'270241aa'), real(z'be0b1010'), real(z'c90c2086'), &
      & real(z'5768b525'), real(z'206f85b3'), real(z'b966d409'), real(z'ce61e49f'), &
      & real(z'5edef90e'), real(z'29d9c998'), real(z'b0d09822'), real(z'c7d7a8b4'), &
      & real(z'59b33d17'), real(z'2eb40d81'), real(z'b7bd5c3b'), real(z'c0ba6cad'), &
      & real(z'edb88320'), real(z'9abfb3b6'), real(z'03b6e20c'), real(z'74b1d29a'), &
      & real(z'ead54739'), real(z'9dd277af'), real(z'04db2615'), real(z'73dc1683'), &
      & real(z'e3630b12'), real(z'94643b84'), real(z'0d6d6a3e'), real(z'7a6a5aa8'), &
      & real(z'e40ecf0b'), real(z'9309ff9d'), real(z'0a00ae27'), real(z'7d079eb1'), &
      & real(z'f00f9344'), real(z'8708a3d2'), real(z'1e01f268'), real(z'6906c2fe'), &
      & real(z'f762575d'), real(z'806567cb'), real(z'196c3671'), real(z'6e6b06e7'), &
      & real(z'fed41b76'), real(z'89d32be0'), real(z'10da7a5a'), real(z'67dd4acc'), &
      & real(z'f9b9df6f'), real(z'8ebeeff9'), real(z'17b7be43'), real(z'60b08ed5'), &
      & real(z'd6d6a3e8'), real(z'a1d1937e'), real(z'38d8c2c4'), real(z'4fdff252'), &
      & real(z'd1bb67f1'), real(z'a6bc5767'), real(z'3fb506dd'), real(z'48b2364b'), &
      & real(z'd80d2bda'), real(z'af0a1b4c'), real(z'36034af6'), real(z'41047a60'), &
      & real(z'df60efc3'), real(z'a867df55'), real(z'316e8eef'), real(z'4669be79'), &
      & real(z'cb61b38c'), real(z'bc66831a'), real(z'256fd2a0'), real(z'5268e236'), &
      & real(z'cc0c7795'), real(z'bb0b4703'), real(z'220216b9'), real(z'5505262f'), &
      & real(z'c5ba3bbe'), real(z'b2bd0b28'), real(z'2bb45a92'), real(z'5cb36a04'), &
      & real(z'c2d7ffa7'), real(z'b5d0cf31'), real(z'2cd99e8b'), real(z'5bdeae1d'), &
      & real(z'9b64c2b0'), real(z'ec63f226'), real(z'756aa39c'), real(z'026d930a'), &
      & real(z'9c0906a9'), real(z'eb0e363f'), real(z'72076785'), real(z'05005713'), &
      & real(z'95bf4a82'), real(z'e2b87a14'), real(z'7bb12bae'), real(z'0cb61b38'), &
      & real(z'92d28e9b'), real(z'e5d5be0d'), real(z'7cdcefb7'), real(z'0bdbdf21'), &
      & real(z'86d3d2d4'), real(z'f1d4e242'), real(z'68ddb3f8'), real(z'1fda836e'), &
      & real(z'81be16cd'), real(z'f6b9265b'), real(z'6fb077e1'), real(z'18b74777'), &
      & real(z'88085ae6'), real(z'ff0f6a70'), real(z'66063bca'), real(z'11010b5c'), &
      & real(z'8f659eff'), real(z'f862ae69'), real(z'616bffd3'), real(z'166ccf45'), &
      & real(z'a00ae278'), real(z'd70dd2ee'), real(z'4e048354'), real(z'3903b3c2'), &
      & real(z'a7672661'), real(z'd06016f7'), real(z'4969474d'), real(z'3e6e77db'), &
      & real(z'aed16a4a'), real(z'd9d65adc'), real(z'40df0b66'), real(z'37d83bf0'), &
      & real(z'a9bcae53'), real(z'debb9ec5'), real(z'47b2cf7f'), real(z'30b5ffe9'), &
      & real(z'bdbdf21c'), real(z'cabac28a'), real(z'53b39330'), real(z'24b4a3a6'), &
      & real(z'bad03605'), real(z'cdd70693'), real(z'54de5729'), real(z'23d967bf'), &
      & real(z'b3667a2e'), real(z'c4614ab8'), real(z'5d681b02'), real(z'2a6f2b94'), &
      & real(z'b40bbe37'), real(z'c30c8ea1'), real(z'5a05df1b'), real(z'2d02ef8d')],&
      & crc_stub)

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
