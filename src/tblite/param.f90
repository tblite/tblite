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

!> @dir tblite/param
!> Definition of the individual components of the parametrization data.

!> @file tblite/param.f90
!> Reexport of the parametrization records and definition of the high-level representation
!> if the parametrization data format.

!> Defininition of the parametrization data format
module tblite_param
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io_symbols, only : to_number, symbol_length
   use tblite_param_charge, only : charge_record
   use tblite_param_dispersion, only : dispersion_record
   use tblite_param_element, only : element_record
   use tblite_param_halogen, only : halogen_record
   use tblite_param_hamiltonian, only : hamiltonian_record
   use tblite_param_mask, only : param_mask, count
   use tblite_param_multipole, only : multipole_record
   use tblite_param_repulsion, only : repulsion_record
   use tblite_param_serde, only : serde_record
   use tblite_param_thirdorder, only : thirdorder_record
   use tblite_param_post_processing, only :  post_processing_param_list, molecular_multipole_record
   use tblite_toml, only : toml_table, toml_key, get_value, set_value, add_table
   implicit none
   private

   public :: param_record, param_mask, count
   public :: charge_record, dispersion_record, element_record, halogen_record, &
      & hamiltonian_record, multipole_record, repulsion_record, thirdorder_record, &
      & post_processing_param_list, molecular_multipole_record


   character(len=*), parameter :: k_dispersion = "dispersion", k_repulsion = "repulsion", &
      & k_charge = "charge", k_thirdorder = "thirdorder", k_multipole = "multipole", &
      & k_halogen = "halogen", k_hamiltonian = "hamiltonian", k_element = "element", &
      & k_meta = "meta", k_version = "version", k_name = "name", k_reference = "reference", &
      & k_format = "format", k_post_proc = "post-processing"

   !> Current parameter format version
   integer, parameter :: current_format = 1


   !> Complete self-contained representation of a complete parametrization record.
   type, extends(serde_record) :: param_record
      !> Version of the represented method
      integer :: version = 0
      !> Name of the represented method
      character(len=:), allocatable :: name
      !> References relevant for the parametrization records
      character(len=:), allocatable :: reference
      !> Definition of the Hamiltonian, always required
      type(hamiltonian_record) :: hamiltonian
      !> Definition of the dispersion correction
      type(dispersion_record), allocatable :: dispersion
      !> Definition of the repulsion contribution
      type(repulsion_record), allocatable :: repulsion
      !> Definition of the isotropic second-order charge interactions
      type(charge_record), allocatable :: charge
      !> Definition of the anisotropic second-order multipolar interactions
      type(multipole_record), allocatable :: multipole
      !> Definition of the halogen bonding correction
      type(halogen_record), allocatable :: halogen
      !> Definition of the isotropic third-order charge interactions
      type(thirdorder_record), allocatable :: thirdorder
      !> Element specific parameter records
      type(element_record), allocatable :: record(:)
      !> Abstract post processing class 
      type(post_processing_param_list), allocatable :: post_proc
   contains
      generic :: load => load_from_array
      generic :: dump => dump_to_array
      !> Read parametrization data from TOML data structure
      procedure :: load_from_toml
      !> Write parametrization data to TOML data structure
      procedure :: dump_to_toml
      !> Read parametrization data from parameter array
      procedure, private :: load_from_array
      !> Write parametrization data to parameter array
      procedure, private :: dump_to_array
      !> Get element record
      procedure :: get
   end type param_record


contains


!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   !> Instance of the parametrization data
   class(param_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child
   integer :: this_format

   call get_value(table, k_meta, child)
   call get_value(child, k_format, this_format, current_format)
   if (this_format < 1) then
      call fatal_error(error, "Parametrer file format must be a positive version")
      return
   end if
   if (this_format > current_format) then
      call fatal_error(error, "Requested future parameter file format which is not yet supported")
      return
   end if
   call get_value(child, k_version, self%version, 0)
   call get_value(child, k_name, self%name)
   call get_value(child, k_reference, self%reference)

   call get_value(table, k_dispersion, child, requested=.false.)
   if (associated(child)) then
      allocate(self%dispersion)
      call self%dispersion%load(child, error)
      if (allocated(error)) return
   end if

   call get_value(table, k_repulsion, child, requested=.false.)
   if (associated(child)) then
      allocate(self%repulsion)
      call self%repulsion%load(child, error)
      if (allocated(error)) return
   end if

   call get_value(table, k_halogen, child, requested=.false.)
   if (associated(child)) then
      allocate(self%halogen)
      call self%halogen%load(child, error)
      if (allocated(error)) return
   end if

   call get_value(table, k_charge, child, requested=.false.)
   if (associated(child)) then
      allocate(self%charge)
      call self%charge%load(child, error)
      if (allocated(error)) return
   end if

   call get_value(table, k_thirdorder, child, requested=.false.)
   if (associated(child)) then
      allocate(self%thirdorder)
      call self%thirdorder%load(child, error)
      if (allocated(error)) return
   end if

   call get_value(table, k_multipole, child, requested=.false.)
   if (associated(child)) then
      allocate(self%multipole)
      call self%multipole%load(child, error)
      if (allocated(error)) return
   end if

   call get_value(table, k_post_proc, child, requested=.false.)
   if (associated(child)) then
      allocate(self%post_proc)
      call self%post_proc%load(child, error)
   end if 

   call get_value(table, k_element, child)
   call records_from_table(self%record, child, error)
   if (allocated(error)) return

   call get_value(table, k_hamiltonian, child, requested=.false.)
   if (associated(child)) then
      self%hamiltonian%sym = self%record%sym
      call self%hamiltonian%load(child, error)
      if (allocated(error)) return
   else
      call fatal_error(error, "Hamiltonian entry is missing")
      if (allocated(error)) return
   end if
end subroutine load_from_toml


!> Write parametrization data to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(param_record), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child

   call add_table(table, k_meta, child)
   if (allocated(self%name)) call set_value(child, k_name, self%name)
   if (allocated(self%reference)) call set_value(child, k_reference, self%reference)
   call set_value(child, k_version, self%version)
   call set_value(child, k_format, current_format)

   call add_table(table, k_hamiltonian, child)
   call self%hamiltonian%dump(child, error)
   if (allocated(error)) return

   if (allocated(self%dispersion)) then
      call add_table(table, k_dispersion, child)
      call self%dispersion%dump(child, error)
      if (allocated(error)) return
   end if

   if (allocated(self%repulsion)) then
      call add_table(table, k_repulsion, child)
      call self%repulsion%dump(child, error)
      if (allocated(error)) return
   end if

   if (allocated(self%halogen)) then
      call add_table(table, k_halogen, child)
      call self%halogen%dump(child, error)
      if (allocated(error)) return
   end if

   if (allocated(self%charge)) then
      call add_table(table, k_charge, child)
      call self%charge%dump(child, error)
      if (allocated(error)) return
   end if

   if (allocated(self%thirdorder)) then
      call add_table(table, k_thirdorder, child)
      call self%thirdorder%dump(child, error)
      if (allocated(error)) return
   end if

   if (allocated(self%multipole)) then
      call add_table(table, k_multipole, child)
      call self%multipole%dump(child, error)
      if (allocated(error)) return
   end if

   if (allocated(self%post_proc)) then
      call add_table(table, k_post_proc, child)
      call self%post_proc%dump(child, error)
      if (allocated(error)) return
   end if

   call add_table(table, k_element, child)
   call records_to_table(self%record, child, error)
   if (allocated(error)) return
end subroutine dump_to_toml


!> Deserialize records from a table by iterating over all entires
subroutine records_from_table(record, table, error)
   !> List of all element records
   type(element_record), allocatable, intent(out) :: record(:)
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ii
   type(toml_key), allocatable :: list(:)
   type(toml_table), pointer :: child

   call table%get_keys(list)
   allocate(record(size(list)))

   do ii = 1, size(list)
      call get_value(table, list(ii)%key, child)
      call record(ii)%load(child, error)
      if (allocated(error)) exit
   end do
   if (allocated(error)) return
end subroutine records_from_table


!> Serialize records to a table by iterating over all entries
subroutine records_to_table(record, table, error)
   !> List of all element records
   type(element_record), intent(in) :: record(:)
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ii
   type(toml_table), pointer :: child

   do ii = 1, size(record)
      call add_table(table, trim(record(ii)%sym), child)
      call record(ii)%dump(child, error)
      if (allocated(error)) exit
   end do
   if (allocated(error)) return
end subroutine records_to_table


!> Get the position of an element in the parameter records
pure subroutine get(self, sym, num, pos)
   !> Instance of the parametrization records
   class(param_record), target, intent(in) :: self
   !> Symbol of the element
   character(len=*), intent(in) :: sym
   !> Atomic number of the element
   integer, intent(in) :: num
   !> Position in the records
   integer, intent(out) :: pos

   integer :: ii

   pos = 0
   if (.not.allocated(self%record)) return

   do ii = 1, size(self%record)
      if (self%record(ii)%sym == sym) then
         pos = ii
         exit
      end if
   end do
   if (pos /= 0) return

   do ii = 1, size(self%record)
      if (self%record(ii)%num == num) then
         pos = ii
         exit
      end if
   end do
   if (pos /= 0) return
end subroutine get

!> Read parametrization data from parameter array
subroutine load_from_array(self, array, base, mask, error)
   class(param_record), intent(inout) :: self
   real(wp), intent(in) :: array(:)
   type(param_record), intent(in) :: base
   type(param_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

   integer :: ii, ir
   integer :: offset

   select type(self)
   type is (param_record)
      self = base
   end select

   offset = 0
   call self%hamiltonian%load(array, offset, base%hamiltonian, mask%hamiltonian, error)
   if (allocated(error)) return

   if (allocated(self%dispersion)) then
      call self%dispersion%load(array, offset, base%dispersion, mask%dispersion, error)
      if (allocated(error)) return
   end if

   if (allocated(self%repulsion)) then
      call self%repulsion%load(array, offset, base%repulsion, mask%repulsion, error)
      if (allocated(error)) return
   end if

   if (allocated(self%halogen)) then
      call self%halogen%load(array, offset, base%halogen, mask%halogen, error)
      if (allocated(error)) return
   end if

   if (allocated(self%charge)) then
      call self%charge%load(array, offset, base%charge, mask%charge, error)
      if (allocated(error)) return
   end if

   if (allocated(self%thirdorder)) then
      call self%thirdorder%load(array, offset, base%thirdorder, mask%thirdorder, error)
      if (allocated(error)) return
   end if

   if (allocated(self%multipole)) then
      call self%multipole%load(array, offset, base%multipole, mask%multipole, error)
      if (allocated(error)) return
   end if

   do ii = 1, size(mask%record)
      call self%get(mask%record(ii)%sym, mask%record(ii)%num, ir)
      call self%record(ir)%load(array, offset, base%record(ir), mask%record(ii), error)
      if (allocated(error)) exit
   end do
   if (allocated(error)) return
end subroutine load_from_array


!> Write parametrization data to parameter array
subroutine dump_to_array(self, array, mask, error)
   class(param_record), intent(in) :: self
   real(wp), intent(inout) :: array(:)
   type(param_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

   integer :: ii, ir
   integer :: offset

   offset = 0
   call self%hamiltonian%dump(array, offset, mask%hamiltonian, error)
   if (allocated(error)) return

   if (allocated(self%dispersion)) then
      call self%dispersion%dump(array, offset, mask%dispersion, error)
      if (allocated(error)) return
   end if

   if (allocated(self%repulsion)) then
      call self%repulsion%dump(array, offset, mask%repulsion, error)
      if (allocated(error)) return
   end if

   if (allocated(self%halogen)) then
      call self%halogen%dump(array, offset, mask%halogen, error)
      if (allocated(error)) return
   end if

   if (allocated(self%charge)) then
      call self%charge%dump(array, offset, mask%charge, error)
      if (allocated(error)) return
   end if

   if (allocated(self%thirdorder)) then
      call self%thirdorder%dump(array, offset, mask%thirdorder, error)
      if (allocated(error)) return
   end if

   if (allocated(self%multipole)) then
      call self%multipole%dump(array, offset, mask%multipole, error)
      if (allocated(error)) return
   end if

   do ii = 1, size(mask%record)
      call self%get(mask%record(ii)%sym, mask%record(ii)%num, ir)
      call self%record(ir)%dump(array, offset, mask%record(ii), error)
      if (allocated(error)) exit
   end do
   if (allocated(error)) return
end subroutine dump_to_array

end module tblite_param
