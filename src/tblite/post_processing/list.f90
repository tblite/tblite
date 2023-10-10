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

!> @file tblite/post_processing/list.f90
!> Implements post processing conatiner list, collcting post processing methods.
module tblite_post_processing_list
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_post_processing_type, only : post_processing_type
   use tblite_param_post_processing, only : post_processing_param_list
   use tblite_toml, only : toml_error, toml_parse, toml_table, get_value
   use tblite_param, only : param_record
   use tblite_param_serde, only : serde_record
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_context, only : context_type
   use tblite_container, only : container_list
   use tblite_timer, only : timer_type, format_time
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_container_cache, only : container_cache
   use tblite_double_dictionary, only : double_dictionary_type
   use mctc_io, only : structure_type
   use tblite_integral_type, only : integral_type
   use tblite_results, only : results_type
   use tblite_post_processing_bond_orders, only : new_wbo, wiberg_bond_orders
   use tblite_post_processing_molecular_moments, only : new_molecular_moments, molecular_moments
   use tblite_param_molecular_moments, only : molecular_multipole_record
   implicit none
   private

   public :: new_post_processing, post_processing_type

   type :: post_processing_record
      class(post_processing_type), allocatable :: pproc
   end type

   type, public :: post_processing_list
      type(post_processing_record), allocatable :: list(:)
      integer :: n = 0
      type(double_dictionary_type) :: dict
   contains
      procedure :: compute
      procedure :: pack_res
      procedure :: print_csv
      procedure :: info
      procedure :: print_timer
      procedure :: push
   end type

   interface new_post_processing
      procedure :: new_post_processing_param
      procedure :: new_post_processing_cli
   end interface
   logical :: print_csv_bool =.false.
contains

subroutine print_timer(self, prlevel, ctx)
   !> Instance of the interaction container
   class(post_processing_list), intent(in) :: self
   integer :: prlevel
   type(context_type) :: ctx
   integer :: i
   do i = 1, self%n
      call self%list(i)%pproc%print_timer(prlevel, ctx)
   end do
end subroutine

subroutine pack_res(self, mol, res)
   class(post_processing_list),intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(results_type), intent(inout) :: res
   integer :: i, n
   real(wp), allocatable :: tmp_array(:)
   character(len=:), allocatable :: tmp_label

   !allocate(res%dict)
   res%dict = self%dict
end subroutine

subroutine print_csv(self, mol)
   class(post_processing_list),intent(in) :: self
   type(structure_type) :: mol
   integer :: n, i, out, j, nat
   character(len=:), allocatable :: tmp_label
   real(wp), allocatable :: tmp_array(:, :), array(:)
   integer, allocatable :: z_array(:)

   if (.not.(print_csv_bool)) return

   associate(dict => self%dict)
      n = dict%get_n_entries()
      out = 42

      open (file='post_processing_tblite.csv', newunit=out)

      allocate(z_array(mol%nat))
      do i=1, mol%nat
         z_array = mol%num(mol%id(i))
      end do
      write(out, '(a)', advance="no") trim("Atom")//','
      allocate(tmp_array(n, mol%nat))

      do i=1, n
         call dict%get_label(i, tmp_label)
         call dict%get_entry(i, array)
         tmp_array(i, :) = array
         if (i == n) then
            write(out, '(a)') trim(tmp_label)
            cycle
         end if
         write(out, '(a)', advance="no") trim(tmp_label)//','

      end do



      do j = 1, mol%nat
         write(out, '(i2,a)', advance="no") z_array(j), ','
         do i = 1, n
            if (i == n) then
               write(tmp_label, '(f14.8)') tmp_array(i, j)
               write(out, '(a)') trim(adjustl(tmp_label))
               cycle
            end if
            write(tmp_label, '(f14.8)') tmp_array(i, j)
            write(out, '(a)', advance='no') trim(adjustl(tmp_label))//","
         end do
      end do
   end associate

end subroutine

subroutine compute(self, mol, wfn, int, calc, c_list, ctx, prlevel)
   class(post_processing_list) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction strcuture data
   type(wavefunction_type), intent(in) :: wfn
   !> integral container
   type(integral_type) :: int
   !> calculator instance
   type(xtb_calculator), intent(in) :: calc
   !> Context container for writing to stdout
   type(context_type), intent(inout) :: ctx
   type(container_cache), intent(inout) :: c_list(:)
   integer :: prlevel
   integer :: i
   do i = 1, self%n
      call self%list(i)%pproc%compute(mol, wfn, int, calc, c_list, ctx, prlevel, self%dict)
   end do
end subroutine

pure function info(self, verbosity, indent) result(str)
   !> Instance of the interaction container
   class(post_processing_list), intent(in) :: self
   !> Verbosity level
   integer, intent(in) :: verbosity
   !> Indentation level
   character(len=*), intent(in) :: indent
   !> Information on the container
   character(len=:), allocatable :: str
   integer :: i
   character(len=*), parameter :: nl = new_line('a')

   str = "Post processing:"// nl

   do i = 1, self%n
      str = str // self%list(i)%pproc%info(verbosity, indent)
      str = str // nl
   end do

end function info

subroutine new_post_processing_param(self, param)
   !> Instance of the xTB evaluator
   class(post_processing_type), allocatable, intent(inout) :: self
   type(post_processing_param_list) :: param
   integer :: i
   do i = 1, size(param%list)
      select type(par => param%list(i)%record)
      type is (molecular_multipole_record)
         block
            type(molecular_moments), allocatable :: tmp
            allocate(tmp)
            call new_molecular_moments(tmp, par)
            call move_alloc(tmp, self)
         end block
      end select
   end do
end subroutine

subroutine new_post_processing_cli(self, config, error)
   class(post_processing_type), allocatable, intent(inout) :: self
   type(error_type), intent(inout) , allocatable:: error
   character(len=:), allocatable :: config
   type(post_processing_param_list), allocatable :: param
   allocate(param)
   select case(config)
   case("wbo")
      block
         type(wiberg_bond_orders), allocatable :: wbo_tmp
         allocate(wbo_tmp)
         call new_wbo(wbo_tmp)
         call move_alloc(wbo_tmp, self)
         return
      end block
   case("molmom")
      block
         type(molecular_moments), allocatable :: molmom_tmp
         allocate(molmom_tmp)
         call new_molecular_moments(molmom_tmp)
         call move_alloc(molmom_tmp, self)
      end block
   case default
      block
         type(toml_table), allocatable :: table
         integer :: io, stat
         type(toml_error), allocatable :: t_error
         type(param_record) :: record
         type(toml_table), pointer :: child

         open(file=config, newunit=io, status="old")

         call toml_parse(table, io, t_error)
         close(io)
         if (allocated(t_error)) then
            allocate(error)
            call move_alloc(t_error%message, error%message)
            return
         end if
         if (allocated(error)) then
            call fatal_error(error, "File name provided could not be parsed as a toml table")
         end if
         call get_value(table, "post-processing", child, requested=.false.)
         if (associated(child)) then
            call param%load(child, error)
         else
            call fatal_error(error, "Could not find post-processing key in toml file.")
         end if
         call new_post_processing(self, param)
      end block
   end select

   if (.not.allocated(self)) then 
      call fatal_error(error, "Outgoing post processing container is not allocated!")
   end if
end subroutine

subroutine push(self, record)
   class(post_processing_list), intent(inout) :: self
   class(post_processing_type), allocatable, intent(inout) :: record

   if (.not.allocated(self%list)) call resize(self%list)
   if (self%n >= size(self%list)) then
      call resize(self%list)
   end if
   write(*,*) size(self%list)
   self%n = self%n + 1
   write(*,*) self%n
   call move_alloc(record, self%list(self%n)%pproc)

end subroutine push

pure subroutine resize(list, n)
   !> Instance of the array to be resized
   type(post_processing_record), allocatable, intent(inout) :: list(:)
   !> Dimension of the final array size
   integer, intent(in), optional :: n

   type(post_processing_record), allocatable :: tmp(:)
   integer :: this_size, new_size, item
   integer, parameter :: initial_size = 20

   if (allocated(list)) then
      this_size = size(list, 1)
      call move_alloc(list, tmp)
   else
      this_size = initial_size
   end if

   if (present(n)) then
      new_size = n
   else
      new_size = this_size + this_size/2 + 1
   end if

   allocate(list(new_size))

   if (allocated(tmp)) then
      this_size = min(size(tmp, 1), size(list, 1))
      do item = 1, this_size
         call move_alloc(tmp(item)%pproc, list(item)%pproc)
      end do
      deallocate(tmp)
   end if

end subroutine resize
end module
