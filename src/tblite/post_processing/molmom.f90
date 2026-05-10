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

!> @file tblite/post_processing/molmom.f90
!> Implements the calculation of molecular moments as post processing methods.
module tblite_post_processing_molmom
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_container_list, only : cache_list
   use tblite_context, only : context_type
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_integral_type, only : integral_type
   use tblite_output_format, only : format_string
   use tblite_param_post_processing_molmom, only : molmom_record
   use tblite_post_processing_type, only : post_processing_type
   use tblite_results, only : results_type
   use tblite_timer, only : timer_type, format_time
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_wavefunction_mulliken, only : get_molecular_dipole_moment, &
      & get_molecular_quadrupole_moment
   implicit none
   private

   public :: molecular_moments, new_molecular_moments

   !> Molecular moments as post-processing method  
   type, extends(post_processing_type) :: molecular_moments
      !> Perform dipole moment calculation
      logical :: comp_dipm = .true.
      !> Perform quadrupole moment calculation
      logical :: comp_qm = .true.
   contains
      !> Calculate molecular moments
      procedure :: compute
      !> Print timings
      procedure :: print_timer
   end type molecular_moments

   character(len=27), parameter :: label = "Molecular Multipole Moments"

contains

subroutine new_molecular_moments(self, param)
   !> Instance of the molecular moments post-processing
   type(molecular_moments), intent(out) :: self
   !> Molecular multipole parameterization
   type(molmom_record), intent(in) :: param
   
   self%label = label

   self%comp_dipm = param%moldipm
   self%comp_qm = param%molqp

end subroutine new_molecular_moments

subroutine compute(self, mol, wfn, ints, calc, caches, ctx, timer, prlevel, dict)
   !> Instance of the molecular moments post-processing
   class(molecular_moments),intent(in) :: self
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
   real(wp) :: dipm(3), qp(6)

   call timer%push("molmom")
   if (self%comp_dipm) then
      call timer%push("dipole")
      call get_molecular_dipole_moment(mol, wfn%qat(:, 1), wfn%dpat(:, :, 1), dipm)
      call dict%add_entry("molecular-dipole", dipm)
      call timer%pop()
   end if
   if (self%comp_qm) then
      call timer%push("quadrupole")
      call get_molecular_quadrupole_moment(mol, wfn%qat(:, 1), wfn%dpat(:, :, 1), &
         & wfn%qpat(:, :, 1), qp)
      call dict%add_entry("molecular-quadrupole", qp)
      call timer%pop()
   end if
   call timer%pop()

end subroutine compute

subroutine print_timer(self, timer, prlevel, ctx)
   !> Instance of the molecular moments post-processing
   class(molecular_moments), intent(in) :: self
   !> Timer instance
   type(timer_type), intent(in) :: timer
   !> Print level
   integer, intent(in) :: prlevel
   !> Context container for writing to stdout
   type(context_type), intent(inout) :: ctx

   real(wp) :: ttime, stime
   integer :: it
   character(len=*), parameter :: labels(*) = [character(len=20):: &
      "dipole", "quadrupole" ]

   if (prlevel > 2) then
      call ctx%message(label//" timing details:")
      ttime = timer%get("molmom")
      call ctx%message(" total:"//repeat(" ", 16)//format_time(ttime))
      do it = 1, size(labels)
         stime = timer%get(labels(it))
         if (stime <= epsilon(0.0_wp)) cycle
         call ctx%message(" - "//labels(it)//format_time(stime) &
            & //" ("//format_string(int(stime/ttime*100), '(i3)')//"%)")
      end do
      call ctx%message("")
   end if

end subroutine print_timer

end module tblite_post_processing_molmom
