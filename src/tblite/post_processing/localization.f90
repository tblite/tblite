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

!> @file tblite/post_processing/localization.f90
!> Implements localized molecular orbitals as a post processing method.
module tblite_post_processing_localization
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_container_list, only : cache_list
   use tblite_context, only : context_type
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_integral_type, only : integral_type
   use tblite_post_processing_type, only : post_processing_type
   use tblite_timer, only : timer_type, format_time
   use tblite_wavefunction_localization, only : localization_type, &
      & localization_method, fosterboys_localization_type, &
      & new_fosterboys_localization, get_orbital_centers
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: new_orbital_localization, orbital_localization

   !> Localized molecular orbitals as a post-processing method
   type, extends(post_processing_type) :: orbital_localization
      !> Localization object used to construct the localized orbitals
      class(localization_type), allocatable :: localizer
   contains
      !> Calculate the localized molecular orbitals
      procedure :: compute
      !> Print timings
      procedure :: print_timer
   end type orbital_localization

   character(len=*), parameter :: label = "Localized molecular orbitals"

contains


!> Create a new orbital localization post-processing method.
subroutine new_orbital_localization(self, error, method_id)
   !> Instance of the orbital localization post-processing
   type(orbital_localization), intent(out) :: self
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Integer identifier of the localization method to use
   integer, intent(in) :: method_id

   self%label = label

   select case(method_id)
   case(localization_method%fosterboys)
      block
         type(fosterboys_localization_type), allocatable :: tmp
         allocate(tmp)
         call new_fosterboys_localization(tmp)
         call move_alloc(tmp, self%localizer)
      end block
   case default
      call fatal_error(error, "Unhandled orbital localization method identifier")
   end select

end subroutine new_orbital_localization


subroutine compute(self, mol, wfn, ints, calc, caches, accuracy, ctx, timer, &
   & prlevel, dict)
   !> Instance of the orbital localization post-processing
   class(orbital_localization), intent(in) :: self
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
   !> Accuracy for computation
   real(wp), intent(in) :: accuracy
   !> Context container for writing to stdout
   type(context_type), intent(inout) :: ctx
   !> Timer instance
   type(timer_type), intent(inout) :: timer
   !> Print level
   integer, intent(in) :: prlevel
   !> Dictionary for storing results
   type(double_dictionary_type), intent(inout) :: dict

   real(wp), allocatable :: coeff_local(:, :, :), centers(:, :, :)
   type(error_type), allocatable :: error

   call timer%push("localization")

   allocate(coeff_local(calc%bas%nao, calc%bas%nao, wfn%nspin))
   call self%localizer%localize(mol, calc%bas, ints%overlap, ints%dipole, wfn%coeff, &
      & wfn%emo, wfn%nel(:wfn%nspin), accuracy, coeff_local, error)
   if (allocated(error)) then
      call ctx%set_error(error)
      call timer%pop()
      return
   end if

   call dict%add_entry("localized-orbitals", coeff_local)

   ! Setup localized orbital centers
   allocate(centers(3, calc%bas%nao, wfn%nspin))
   call get_orbital_centers(mol, calc%bas, ints%overlap, ints%dipole, coeff_local, &
      & wfn%nel(:wfn%nspin), centers)
   call dict%add_entry("localized-centers", centers)

   call timer%pop()

end subroutine compute

subroutine print_timer(self, timer, prlevel, ctx)
   !> Instance of the orbital localization post-processing
   class(orbital_localization), intent(in) :: self
   !> Timer instance
   type(timer_type), intent(in) :: timer
   !> Print level
   integer, intent(in) :: prlevel
   !> Context container for writing to stdout
   type(context_type), intent(inout) :: ctx

   real(wp) :: ttime

   if (prlevel > 2) then
      call ctx%message(label//" timing details:")
      ttime = timer%get("localization")
      call ctx%message(" total:"//repeat(" ", 16)//format_time(ttime))
      call ctx%message("")
   end if

end subroutine print_timer

end module tblite_post_processing_localization
