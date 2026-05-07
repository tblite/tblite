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
   use mctc_io, only : structure_type
   use tblite_context, only : context_type
   use tblite_container_list, only : cache_list
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_integral_type, only : integral_type
   use tblite_param_molecular_moments, only : molecular_multipole_record
   use tblite_param_post_processing, only : post_processing_param_list
   use tblite_param_serde, only : serde_record
   use tblite_param_xtbml_features, only : xtbml_features_record
   use tblite_post_processing_bond_orders, only : new_wbo, wiberg_bond_orders
   use tblite_post_processing_molecular_moments, only : new_molecular_moments, &
      & molecular_moments
   use tblite_post_processing_type, only : post_processing_type
   use tblite_post_processing_xtbml_type, only : xtbml_type, new_xtbml_features
   use tblite_results, only : results_type
   use tblite_timer, only : timer_type
   use tblite_toml, only : toml_error, toml_parse, toml_table, get_value
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: add_post_processing

   !> Wrapped post-processing type for creation of lists
   type :: post_processing_record
      !> Actual post-processing method
      class(post_processing_type), allocatable :: pproc
   end type post_processing_record

   !> List of post-processing methods
   type, public :: post_processing_list
      !> Raw list of post-processing methods
      type(post_processing_record), allocatable :: list(:)
      !> Number of post-processing methods in the list
      integer :: npp = 0
   contains
      !> Perform post-processing methods in the list
      procedure :: compute
      !> List the used post-processing methods
      procedure :: info
      !> Print timings
      procedure :: print_timer
      !> Add post-processing method to list
      procedure :: push
   end type post_processing_list

   interface add_post_processing
      procedure :: add_post_processing_param
      procedure :: add_post_processing_cli
   end interface add_post_processing

contains

subroutine print_timer(self, timer, prlevel, ctx)
   !> Instance of the interaction container
   class(post_processing_list), intent(in) :: self
   !> Timer instance
   type(timer_type), intent(in) :: timer
   !> Print level
   integer, intent(in) :: prlevel
   !> Context container for writing to stdout
   type(context_type), intent(inout) :: ctx

   integer :: ipp
   do ipp = 1, self%npp
      call self%list(ipp)%pproc%print_timer(timer, prlevel, ctx)
   end do
end subroutine print_timer

subroutine compute(self, mol, wfn, ints, calc, caches, ctx, timer, prlevel, &
   & results)
   !> Instance of the post-processing list
   class(post_processing_list), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction strcuture data
   type(wavefunction_type), intent(in) :: wfn
   !> Integral container
   type(integral_type), intent(in) :: ints
   !> Calculator instance
   type(xtb_calculator), intent(in) :: calc
   !> Context container for writing to stdout
   type(context_type), intent(inout) :: ctx
   !> Cache list for storing caches of various interactions
   type(cache_list), intent(inout) :: caches
   !> Timer instance
   type(timer_type), intent(inout) :: timer
   !> Print level
   integer, intent(in) :: prlevel
   !> Container for storing results
   type(results_type), intent(inout) :: results

   integer :: ipp

   do ipp = 1, self%npp
      call self%list(ipp)%pproc%compute(mol, wfn, ints, calc, caches, ctx, &
         & timer, prlevel, results%dict)
   end do

end subroutine compute

pure function info(self, verbosity, indent) result(str)
   !> Instance of the interaction container
   class(post_processing_list), intent(in) :: self
   !> Verbosity level
   integer, intent(in) :: verbosity
   !> Indentation level
   character(len=*), intent(in) :: indent
   !> Information on the container
   character(len=:), allocatable :: str

   integer :: ipp
   character(len=*), parameter :: nl = new_line('a')

   str = "Post processing:"// nl

   do ipp = 1, self%npp
      str = str // self%list(ipp)%pproc%info(verbosity, indent)
      str = str // nl
   end do

end function info

subroutine add_post_processing_param(self, mol, param, error)
   !> Instance of the post-processing list
   class(post_processing_list), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> List of post-processing parameterizations
   type(post_processing_param_list), intent(in) :: param
   !> Error handling
   type(error_type), intent(inout), allocatable :: error

   integer :: ipp
   do ipp = 1, param%get_n_records()
      select type(par => param%list(ipp)%record)
      type is (molecular_multipole_record)
         block
            type(molecular_moments), allocatable :: tmp
            class(post_processing_type), allocatable :: proc
            allocate(tmp)
            call new_molecular_moments(tmp, par)
            call move_alloc(tmp, proc)
            call self%push(proc)
         end block
      type is (xtbml_features_record)
         block
            type(xtbml_type), allocatable :: tmp
            class(post_processing_type), allocatable :: proc
            allocate(tmp)
            call new_xtbml_features(tmp, mol, par, error)
            if (allocated(error)) return
            call move_alloc(tmp, proc)
            call self%push(proc)
         end block
      end select
   end do

end subroutine add_post_processing_param

subroutine add_post_processing_cli(self, mol, config, error)
   !> Instance of the post-processing list
   class(post_processing_list), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Name of the post-processing method to be added
   character(len=:), intent(in), allocatable :: config
   !> Error handling
   type(error_type), intent(inout) , allocatable:: error

   type(post_processing_param_list) :: param
   class(post_processing_type), allocatable :: tmp_proc

   select case(config)
   case("bond-orders")
      block
         type(wiberg_bond_orders), allocatable :: wbo_tmp
         allocate(wbo_tmp)
         call new_wbo(wbo_tmp)
         call move_alloc(wbo_tmp, tmp_proc)
         call self%push(tmp_proc)
         return
      end block
   case("molmom")
      block
         type(molecular_multipole_record), allocatable :: molmom_tmp
         class(serde_record), allocatable :: tmp
         allocate(molmom_tmp)
         call molmom_tmp%populate_default_param()
         call move_alloc(molmom_tmp, tmp)
         call param%push(tmp)
      end block
   case("xtbml")
      block
         type(xtbml_features_record), allocatable :: ml_param
         class(serde_record), allocatable :: cont
         allocate(ml_param)
         call ml_param%populate_default_param(.false.)
         call move_alloc(ml_param, cont)
         call param%push(cont)
      end block
   case("xtbml-xyz","xtbml_xyz")
      block
         type(xtbml_features_record), allocatable :: ml_param
         class(serde_record), allocatable :: cont 
         allocate(ml_param)
         call ml_param%populate_default_param(.true.)
         call move_alloc(ml_param, cont)
         call param%push(cont)
      end block
   case default
      block
         type(toml_table), allocatable :: table
         integer :: io
         type(toml_error), allocatable :: t_error
         type(toml_table), pointer :: child

         open(file=config, newunit=io, status="old")
         
         call toml_parse(table, io, t_error)
         close(io)
         if (allocated(t_error)) then
            call fatal_error(error, "Could not parse TOML file due to "//t_error%message)
            return
         end if
         call get_value(table, "post-processing", child, requested=.false.)
         if (associated(child)) then
            call param%load(child, error)
         else
            call fatal_error(error, "Could not find post-processing key in toml file.")
         end if
         
      end block   
   end select

   call add_post_processing(self, mol, param, error)

end subroutine add_post_processing_cli

subroutine push(self, record)
   !> Instance of the post-processing list
   class(post_processing_list), intent(inout) :: self
   !> Instance of the post-processing method to be added
   class(post_processing_type), allocatable, intent(inout) :: record
   
   if (.not.allocated(self%list)) call resize(self%list)
   if (is_duplicate(self, record)) return
   if (self%npp >= size(self%list)) then
      call resize(self%list)
   end if
   
   self%npp = self%npp + 1
   call move_alloc(record, self%list(self%npp)%pproc)

end subroutine push

function is_duplicate(self, record) result(duplicate)
   !> Instance of the post-processing list
   class(post_processing_list), intent(inout) :: self
   !> Instance of the post-processing method to be checked
   class(post_processing_type), allocatable, intent(inout) :: record
   logical :: duplicate 
   integer :: ipp
   duplicate = .false.
   do ipp = 1, self%npp
      if (record%label == self%list(ipp)%pproc%label) duplicate = .true.
   end do
end function is_duplicate

subroutine resize(list, n)
   !> Instance of the array to be resized
   type(post_processing_record), allocatable, intent(inout) :: list(:)
   !> Dimension of the final array size
   integer, intent(in), optional :: n

   type(post_processing_record), allocatable :: tmp(:)
   integer :: this_size, new_size, item
   integer, parameter :: initial_size = 1

   if (allocated(list)) then
      this_size = size(list, 1)
      call move_alloc(list, tmp)
   else
      this_size = initial_size
   end if

   if (present(n)) then
      new_size = n
   else
      new_size = this_size + 1
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

end module tblite_post_processing_list
