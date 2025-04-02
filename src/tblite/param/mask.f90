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

!> @file tblite/param/mask.f90
!> Provides a general mask for the complete parameters set

!> Collection of the parameter masking
module tblite_param_mask
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io_symbols, only : to_number, symbol_length
   use tblite_param_charge, only : charge_mask, count
   use tblite_param_dispersion, only : dispersion_mask, count
   use tblite_param_element, only : element_record, element_mask, count
   use tblite_param_halogen, only : halogen_mask, count
   use tblite_param_hamiltonian, only : hamiltonian_mask, count
   use tblite_param_multipole, only : multipole_mask, count
   use tblite_param_repulsion, only : repulsion_mask, count
   use tblite_param_serde, only : serde_record
   use tblite_param_thirdorder, only : thirdorder_mask, count
   use tblite_toml, only : toml_table, toml_array, toml_key, get_value, set_value, &
      & add_table, add_array
   implicit none
   private

   public :: count


   type :: allowed_records
      logical :: hamiltonian = .false.
      logical :: dispersion = .false.
      logical :: repulsion = .false.
      logical :: charge = .false.
      logical :: multipole = .false.
      logical :: halogen = .false.
      logical :: thirdorder = .false.
   end type allowed_records


   !> Definition of the complete parameter mask
   type, public, extends(serde_record) :: param_mask
      !> Definition of the Hamiltonian, always required
      type(hamiltonian_mask), allocatable :: hamiltonian
      !> Definition of the dispersion correction
      type(dispersion_mask), allocatable :: dispersion
      !> Definition of the repulsion contribution
      type(repulsion_mask), allocatable :: repulsion
      !> Definition of the isotropic second-order charge interactions
      type(charge_mask), allocatable :: charge
      !> Definition of the anisotropic second-order multipolar interactions
      type(multipole_mask), allocatable :: multipole
      !> Definition of the halogen bonding correction
      type(halogen_mask), allocatable :: halogen
      !> Definition of the isotropic third-order charge interactions
      type(thirdorder_mask), allocatable :: thirdorder
      !> Element specific parameter masks
      type(element_mask), allocatable :: record(:)
      !> Reference to base parametrization
      type(element_record), pointer :: ref(:) => null()
   contains
      !> Read parametrization mask from TOML data structure
      procedure :: load_from_toml
      !> Write parametrization mask to TOML data structure
      procedure :: dump_to_toml
   end type param_mask


   interface count
      module procedure :: count_mask
   end interface count


contains


!> Read parametrization mask from TOML data structure
subroutine load_from_toml(self, table, error)
   !> Instance of the parametrization data
   class(param_mask), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(allowed_records) :: allowed
   type(toml_table), pointer :: child

   call get_value(table, "hamiltonian", child, requested=.false.)
   allowed%hamiltonian = associated(child)
   if (associated(child)) then
      allocate(self%hamiltonian)
   end if

   call get_value(table, "dispersion", child, requested=.false.)
   allowed%dispersion = associated(child)
   if (associated(child)) then
      allocate(self%dispersion)
   end if

   call get_value(table, "repulsion", child, requested=.false.)
   allowed%repulsion = associated(child)
   if (associated(child)) then
      allocate(self%repulsion)
   end if

   call get_value(table, "halogen", child, requested=.false.)
   allowed%halogen = associated(child)
   if (associated(child)) then
      allocate(self%halogen)
   end if

   call get_value(table, "charge", child, requested=.false.)
   allowed%charge = associated(child)
   if (associated(child)) then
      allocate(self%charge)
   end if

   call get_value(table, "multipole", child, requested=.false.)
   allowed%multipole = associated(child)
   if (associated(child)) then
      allocate(self%multipole)
   end if

   call get_value(table, "thirdorder", child, requested=.false.)
   allowed%thirdorder = associated(child)
   if (associated(child)) then
      allocate(self%thirdorder)
   end if

   call get_value(table, "element", child)
   call records_from_table(self%record, child, self%ref, allowed, error)
   if (allocated(error)) return
end subroutine load_from_toml


!> Deserialize records from a table by iterating over all entires
subroutine records_from_table(record, table, ref, allowed, error)
   !> List of all element records
   type(element_mask), allocatable, intent(out) :: record(:)
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> List of all element records
   type(element_record), intent(in), optional :: ref(:)
   !> Allowed entries in record
   type(allowed_records), intent(in) :: allowed
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ii, nsh, ir, stat
   type(toml_key), allocatable :: list(:)
   type(toml_table), pointer :: child

   call table%get_keys(list)
   allocate(record(size(list)))

   do ii = 1, size(list)
      call get_value(table, list(ii)%key, child)
      if (present(ref)) then
         call get(ref, list(ii)%key, ir)
         if (ir == 0) cycle
         nsh = ref(ir)%nsh
      else
         call get_value(child, "nsh", nsh, stat=stat)
         if (stat /= 0) cycle
      end if
      call record_from_table(record(ii), child, nsh, allowed, error)
      if (allocated(error)) exit
   end do
   if (allocated(error)) return
end subroutine records_from_table


!> Get the position of an element in the parameter records
pure subroutine get(record, sym, pos)
   !> Instance of the parametrization records
   type(element_record), intent(in) :: record(:)
   !> Symbol of the element
   character(len=*), intent(in) :: sym
   !> Position in the records
   integer, intent(out) :: pos

   integer :: num
   integer :: ii

   num = to_number(sym)
   pos = 0

   do ii = 1, size(record)
      if (record(ii)%sym == sym) then
         pos = ii
         exit
      end if
   end do
   if (pos /= 0) return

   do ii = 1, size(record)
      if (record(ii)%num == num) then
         pos = ii
         exit
      end if
   end do
   if (pos /= 0) return
end subroutine get


!> Deserialize records from a table by iterating over all entires
subroutine record_from_table(record, table, nsh, allowed, error)
   !> List of all element records
   type(element_mask), intent(out) :: record
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Number of shells for this record
   integer, intent(in) :: nsh
   !> Allowed entries in record
   type(allowed_records), intent(in) :: allowed
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   record%sym = table%key
   record%num = to_number(table%key)

   allocate(record%levels(nsh), source=.false.)
   allocate(record%slater(nsh), source=.false.)
   allocate(record%kcn(nsh), source=.false.)
   allocate(record%shpoly(nsh), source=.false.)
   allocate(record%lgam(nsh), source=.false.)

   call read_shell_mask(table, "levels", record%levels, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_shell_mask(table, "slater", record%slater, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_shell_mask(table, "shpoly", record%shpoly, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_shell_mask(table, "kcn", record%kcn, allowed%hamiltonian, error)
   if (allocated(error)) return

   call read_atom_mask(table, "gam", record%gam, allowed%charge, error)
   if (allocated(error)) return
   call read_shell_mask(table, "lgam", record%lgam, allowed%charge, error)
   if (allocated(error)) return

   call read_atom_mask(table, "dkernel", record%dkernel, allowed%multipole, error)
   if (allocated(error)) return
   call read_atom_mask(table, "qkernel", record%qkernel, allowed%multipole, error)
   if (allocated(error)) return

   call read_atom_mask(table, "gam3", record%gam3, allowed%thirdorder, error)
   if (allocated(error)) return

   call read_atom_mask(table, "zeff", record%zeff, allowed%repulsion, error)
   if (allocated(error)) return
   call read_atom_mask(table, "arep", record%alpha, allowed%repulsion, error)
   if (allocated(error)) return

   call read_atom_mask(table, "xbond", record%xbond, allowed%halogen, error)
   if (allocated(error)) return
end subroutine record_from_table


subroutine read_atom_mask(table, key, mask, default, error)
   type(toml_table), intent(inout) :: table
   character(len=*), intent(in) :: key
   logical, intent(out) :: mask
   logical, intent(in) :: default
   type(error_type), allocatable, intent(out) :: error

   call get_value(table, key, mask, default)
end subroutine read_atom_mask


subroutine read_shell_mask(table, key, mask, default, error)
   type(toml_table), intent(inout) :: table
   character(len=*), intent(in) :: key
   logical, intent(out) :: mask(:)
   logical, intent(in) :: default
   type(error_type), allocatable, intent(out) :: error

   type(toml_array), pointer :: array
   integer :: ii

   call get_value(table, key, array, requested=.false.)
   if (associated(array)) then
      do ii = 1, size(mask)
         call get_value(array, ii, mask(ii))
      end do
   else
      call get_value(table, key, mask(1), default)
      mask(:) = mask(1)
   end if
end subroutine read_shell_mask


!> Write parametrization mask to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(param_mask), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child

   if (allocated(self%hamiltonian)) then
      call add_table(table, "hamiltonian", child)
   end if

   if (allocated(self%dispersion)) then
      call add_table(table, "dispersion", child)
   end if

   if (allocated(self%repulsion)) then
      call add_table(table, "repulsion", child)
   end if

   if (allocated(self%halogen)) then
      call add_table(table, "halogen", child)
   end if

   if (allocated(self%charge)) then
      call add_table(table, "charge", child)
   end if

   if (allocated(self%multipole)) then
      call add_table(table, "multipole", child)
   end if

   if (allocated(self%thirdorder)) then
      call add_table(table, "thirdorder", child)
   end if

   if (allocated(self%thirdorder)) then
      call add_table(table, "thirdorder", child)
   end if

   call add_table(table, "element", child)
   call records_to_table(self%record, child, error)
   if (allocated(error)) return

end subroutine dump_to_toml


!> Serialize records to a table by iterating over all entries
subroutine records_to_table(record, table, error)
   !> List of all element records
   type(element_mask), intent(in) :: record(:)
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ii
   type(toml_table), pointer :: child

   do ii = 1, size(record)
      call add_table(table, trim(record(ii)%sym), child)
      call record_to_table(record(ii), child, error)
      if (allocated(error)) exit
   end do
   if (allocated(error)) return
end subroutine records_to_table


!> Serialize records to a table by iterating over all entries
subroutine record_to_table(record, table, error)
   !> List of all element records
   type(element_mask), intent(in) :: record
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: i, nsh
   type(toml_array), pointer :: array

   nsh = size(record%levels)

   call add_array(table, "levels", array)
   do i = 1, nsh
      call set_value(array, i, record%levels(i))
   end do

   call add_array(table, "slater", array)
   do i = 1, nsh
      call set_value(array, i, record%slater(i))
   end do

   call add_array(table, "shpoly", array)
   do i = 1, nsh
      call set_value(array, i, record%shpoly(i))
   end do

   call add_array(table, "kcn", array)
   do i = 1, nsh
      call set_value(array, i, record%kcn(i))
   end do

   call set_value(table, "gam", record%gam)
   call add_array(table, "lgam", array)
   do i = 1, nsh
      call set_value(array, i, record%lgam(i))
   end do

   call set_value(table, "gam3", record%gam3)

   call set_value(table, "zeff", record%zeff)
   call set_value(table, "arep", record%alpha)
   call set_value(table, "xbond", record%xbond)

   call set_value(table, "dkernel", record%dkernel)
   call set_value(table, "qkernel", record%qkernel)
end subroutine record_to_table


elemental function count_mask(mask) result(ncount)
   type(param_mask), intent(in) :: mask
   integer :: ncount
   ncount = 0
   if (allocated(mask%hamiltonian)) ncount = ncount + count(mask%hamiltonian)
   if (allocated(mask%dispersion)) ncount = ncount + count(mask%dispersion)
   if (allocated(mask%repulsion)) ncount = ncount + count(mask%repulsion)
   if (allocated(mask%halogen)) ncount = ncount + count(mask%halogen)
   if (allocated(mask%charge)) ncount = ncount + count(mask%charge)
   if (allocated(mask%multipole)) ncount = ncount + count(mask%multipole)
   if (allocated(mask%thirdorder)) ncount = ncount + count(mask%thirdorder)
   if (allocated(mask%record)) ncount = ncount + sum(count(mask%record))
end function count_mask


end module tblite_param_mask
