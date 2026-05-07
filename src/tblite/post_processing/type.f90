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

!> @file tblite/post_processing/type.f90
!> Implements post processing container abstract class, and the collection of computing caches.
module tblite_post_processing_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_container_cache, only : container_cache
   use tblite_container_list, only : cache_list
   use tblite_context, only : context_type
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_integral_type, only : integral_type
   use tblite_results, only : results_type
   use tblite_timer, only : timer_type, format_time
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: post_processing_type, collect_containers_caches

   !> Abstract base class for post-processing methods
   type, abstract :: post_processing_type
      !> Post-processing label
      character(len=:), allocatable :: label
   contains
      !> Perform post-processing method
      procedure(compute), deferred :: compute
      !> Information on the post-processing method
      procedure :: info
      !> Print timings
      procedure(print_timer), deferred :: print_timer
   end type  post_processing_type

   abstract interface
      subroutine compute(self, mol, wfn, ints, calc, caches, ctx, timer, prlevel, dict)
         import :: post_processing_type, structure_type, wavefunction_type, &
            & integral_type, xtb_calculator, context_type, timer_type, &
            & cache_list, double_dictionary_type
         !> Instance of a post-processing method
         class(post_processing_type),intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Wavefunction strcuture data
         type(wavefunction_type), intent(in) :: wfn
         !> Integral container
         type(integral_type), intent(in) :: ints
         !> Calculator instance
         type(xtb_calculator), intent(in) :: calc
         !> Cache list for storing caches of various interactions
         type(cache_list), intent(inout) :: caches
         !> Context container for writing to stdout
         type(context_type), intent(inout) :: ctx
         !> Timer instance
         type(timer_type), intent(inout) :: timer
         !> Print level
         integer, intent(in) :: prlevel
         !> Dictionary for storing results
         type(double_dictionary_type), intent(inout) :: dict
      end subroutine compute

      subroutine print_timer(self, timer, prlevel, ctx)
         import :: post_processing_type, timer_type, context_type
         !> Instance of a post-processing method
         class(post_processing_type), intent(in) :: self
         !> Timer instance
         type(timer_type), intent(in) :: timer
         !> Print level
         integer, intent(in) :: prlevel
         !> Context container for writing to stdout
         type(context_type), intent(inout) :: ctx
      end subroutine print_timer

   end interface

contains


pure function info(self, verbosity, indent) result(str)
   !> Instance of a post-processing method
   class(post_processing_type), intent(in) :: self
   !> Verbosity level
   integer, intent(in) :: verbosity
   !> Indentation level
   character(len=*), intent(in) :: indent
   !> Information on the container
   character(len=:), allocatable :: str

   if (allocated(self%label)) then
      str = indent // self%label
   else
      str = "Unknown"
   end if

end function info

subroutine collect_containers_caches(rcache, ccache, hcache, dcache, icache, calc, caches)
   !> Container cache for the repulsion interaction
   type(container_cache), allocatable, intent(inout) :: rcache
   !> Container cache for the Coulomb interaction
   type(container_cache), allocatable, intent(inout) :: ccache
   !> Container cache for the halogen correction
   type(container_cache), allocatable, intent(inout) :: hcache
   !> Container cache for the dispersion interaction
   type(container_cache), allocatable, intent(inout) :: dcache
   !> Container cache for general interactions
   type(container_cache), allocatable, intent(inout) :: icache
   !> Instance of the calculator
   type(xtb_calculator), intent(in) :: calc
   !> List of caches for use in the post-processing
   type(cache_list), intent(inout) :: caches

   allocate(caches%list(5))
   if (allocated(calc%repulsion)) call move_alloc(rcache%raw, caches%list(1)%raw)
   if (allocated(calc%coulomb)) call move_alloc(ccache%raw, caches%list(2)%raw)
   if (allocated(calc%halogen)) call move_alloc(hcache%raw, caches%list(3)%raw)
   if (allocated(calc%dispersion)) call move_alloc(dcache%raw, caches%list(4)%raw)
   if (allocated(calc%interactions)) call move_alloc(icache%raw, caches%list(5)%raw)

end subroutine collect_containers_caches

end module tblite_post_processing_type
