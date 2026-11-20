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

!> @file tblite/wavefunction/restart.f90
!> Read and write restart files

!> Implementation for reading and writing restart files
module tblite_wavefunction_restart
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_io_data, only : iodata_type, open_iodata_handler
   use tblite_io_hdf5, only : load_hdf5, save_hdf5
   use tblite_io_numpy, only : load_npz, save_npz
   use tblite_output_format, only : format_string
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: load_wavefunction, save_wavefunction

   character(len=*), parameter :: prefix = "tblite0_"
contains

subroutine load_wavefunction(wfn, filename, error, nat, nsh, nao, nspin)
   type(wavefunction_type), intent(out) :: wfn
   character(len=*), intent(in) :: filename
   type(error_type), allocatable, intent(out) :: error
   integer, intent(in), optional :: nat
   integer, intent(in), optional :: nsh
   integer, intent(in), optional :: nao
   integer, intent(in), optional :: nspin

   logical :: exist
   real(wp), allocatable :: kt(:)
   integer, allocatable :: nat_sizes(:), nsh_sizes(:), nao_sizes(:), nspin_sizes(:)
   class(iodata_type), allocatable :: iodata

   inquire(file=filename, exist=exist)
   if (.not.exist) then
      call fatal_error(error, "Could not find '"//filename//"'")
      return
   end if

   call open_iodata_handler(iodata, filename, error)
   if (allocated(error)) return

   call iodata%load(prefix//"kt", kt, error)
   if (allocated(error)) return

   call iodata%load(prefix//"nel", wfn%nel, error)
   if (allocated(error)) return

   call iodata%load(prefix//"n0at", wfn%n0at, error)
   if (allocated(error)) return

   call iodata%load(prefix//"n0sh", wfn%n0sh, error)
   if (allocated(error)) return

   call iodata%load(prefix//"density", wfn%density, error)
   if (allocated(error)) return

   call iodata%load(prefix//"coeff", wfn%coeff, error)
   if (allocated(error)) return

   call iodata%load(prefix//"emo", wfn%emo, error)
   if (allocated(error)) return

   call iodata%load(prefix//"focc", wfn%focc, error)
   if (allocated(error)) return

   call iodata%load(prefix//"qat", wfn%qat, error)
   if (allocated(error)) return

   call iodata%load(prefix//"qsh", wfn%qsh, error)
   if (allocated(error)) return

   call iodata%load(prefix//"dpat", wfn%dpat, error)
   if (allocated(error)) return

   call iodata%load(prefix//"qpat", wfn%qpat, error)
   if (allocated(error)) return

   nat_sizes = [size(wfn%n0at), size(wfn%qat, 1), size(wfn%dpat, 2), size(wfn%qpat, 2)]
   if (present(nat)) nat_sizes = [nat, nat_sizes]
   if (any(nat_sizes(1) /= nat_sizes)) then
      call fatal_error(error, "Dimension mismatch in '"//filename//&
         & "' for number of atoms")
      return
   end if

   nsh_sizes = [size(wfn%n0sh), size(wfn%qsh, 1)]
   if (present(nsh)) nsh_sizes = [nsh, nsh_sizes]
   if (any(nsh_sizes(1) /= nsh_sizes)) then
      call fatal_error(error, "Dimension mismatch in '"//filename//&
         & "' for number of shells")
      return
   end if

   nao_sizes = [size(wfn%density, 1), size(wfn%density, 2), &
      & size(wfn%coeff, 1), size(wfn%coeff, 2), size(wfn%emo, 1), size(wfn%focc, 1)]
   if (present(nao)) nao_sizes = [nao, nao_sizes]
   if (any(nao_sizes(1) /= nao_sizes)) then
      call fatal_error(error, "Dimension mismatch in '"//filename//&
         & "' for number of orbitals")
      return
   end if

   nspin_sizes = [size(wfn%density, 3), size(wfn%coeff, 3), size(wfn%emo, 2)]
   if (present(nspin)) nspin_sizes = [nspin, nspin_sizes]
   if (any(nspin_sizes(1) /= nspin_sizes) .or. (all(size(wfn%nel) /= [2, nspin_sizes(1)]))) then
      call fatal_error(error, "Dimension mismatch in '"//filename//&
         & "' for number of spin channels")
      return
   end if

   if (size(wfn%dpat, 1) /= 3) then
      call fatal_error(error, "Dimension mismatch in '"//filename//&
         & "' for atomic dipole moments "//&
         & format_string(size(wfn%dpat, 1), "(i0)")//" != 3")
      return
   end if
   if (size(wfn%qpat, 1) /= 6) then
      call fatal_error(error, "Dimension mismatch in '"//filename//&
         & "' for atomic quadrupole moments "//&
         & format_string(size(wfn%qpat, 1), "(i0)")//" != 6")
      return
   end if
   if (size(kt) /= 1) then
      call fatal_error(error, "Dimension mismatch in '"//filename//&
         & "' for electronic temperature")
      return
   end if

   wfn%kt = kt(1)
   wfn%nocc = sum(wfn%nel)
   wfn%nspin = size(wfn%density, 3)
   wfn%nuhf = wfn%nel(1) - wfn%nel(2)
end subroutine load_wavefunction

subroutine save_wavefunction(wfn, filename, error)
   type(wavefunction_type), intent(in) :: wfn
   character(len=*), intent(in) :: filename
   type(error_type), allocatable, intent(out) :: error

   integer :: io, stat
   logical :: exist
   real(wp), allocatable :: kt(:)
   character(len=:), allocatable :: msg
   class(iodata_type), allocatable :: iodata

   inquire(file=filename, exist=exist)
   if (exist) then
      open(file=filename, newunit=io)
      close(io, status="delete")
   end if

   call open_iodata_handler(iodata, filename, error)
   if (allocated(error)) return

   kt = [wfn%kt]
   call iodata%save(prefix//"kt", kt, error)
   if (allocated(error)) return

   call iodata%save(prefix//"nel", wfn%nel, error)
   if (allocated(error)) return

   call iodata%save(prefix//"n0at", wfn%n0at, error)
   if (allocated(error)) return

   call iodata%save(prefix//"n0sh", wfn%n0sh, error)
   if (allocated(error)) return

   call iodata%save(prefix//"density", wfn%density, error)
   if (allocated(error)) return
   call iodata%save(prefix//"coeff", wfn%coeff, error)
   if (allocated(error)) return
   call iodata%save(prefix//"emo", wfn%emo, error)
   if (allocated(error)) return
   call iodata%save(prefix//"focc", wfn%focc, error)
   if (allocated(error)) return

   call iodata%save(prefix//"qat", wfn%qat, error)
   if (allocated(error)) return
   call iodata%save(prefix//"qsh", wfn%qsh, error)
   if (allocated(error)) return
   call iodata%save(prefix//"dpat", wfn%dpat, error)
   if (allocated(error)) return
   call iodata%save(prefix//"qpat", wfn%qpat, error)
   if (allocated(error)) return

end subroutine save_wavefunction


end module tblite_wavefunction_restart
