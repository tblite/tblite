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

!> @dir tblite/post-processing/xtb-ml
!> Explicit implementation of different xTB-ML features

!> @file tblite/post-processing/xtb-ml.f90
!> xTB-ML features as a post-processing method
module tblite_post_processing_xtbml
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use tblite_container_list, only : cache_list
   use tblite_context, only : context_type
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_integral_type, only : integral_type
   use tblite_param_post_processing_xtbml, only : xtbml_record
   use tblite_post_processing_type, only : post_processing_type
   use tblite_output_format, only : format_string
   use tblite_timer, only : timer_type, format_time
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_post_processing_xtbml_cache, only : xtbml_cache
   use tblite_post_processing_xtbml_convolution, only : xtbml_convolution, &
      & new_xtbml_convolution
   use tblite_post_processing_xtbml_density, only : xtbml_density_features, &
      & new_xtbml_density_features
   use tblite_post_processing_xtbml_energy, only : xtbml_energy_features, &
      & new_xtbml_energy_features
   use tblite_post_processing_xtbml_geometry, only : xtbml_geometry_features, &
      & new_xtbml_geometry_features
   use tblite_post_processing_xtbml_orbital, only : xtbml_orbital_features, &
      & new_xtbml_orbital_features
   implicit none
   private

   public :: xtbml_type, new_xtbml_features

   !> xTB-ML features as post-processing method
   type, extends(post_processing_type) :: xtbml_type
      !> Geometry-based xTB-ML features
      type(xtbml_geometry_features), allocatable :: geom
      !> Density-based xTB-ML features
      type(xtbml_density_features), allocatable :: dens
      !> Orbital energy-based xTB-ML features
      type(xtbml_orbital_features), allocatable :: orb
      !> Energy-based xTB-ML features
      type(xtbml_energy_features), allocatable :: energy
      !> Convolution of xTB-ML features
      type(xtbml_convolution), allocatable :: conv
   contains
      !> Calculate xTB-ML features
      procedure :: compute
      !> Information on the post-processing method
      procedure :: info
      !> Print timings
      procedure :: print_timer
   end type xtbml_type

   character(len=*), parameter :: label = "  xtbml features:"

contains

subroutine new_xtbml_features(self, mol, param, error)
   !> Instance of the xTB-ML features post-processing
   type(xtbml_type), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parameterization for the xTB-ML features
   type(xtbml_record), intent(in) :: param
   !> Error handling
   type(error_type), intent(inout) , allocatable:: error

   real(wp), allocatable :: rcov_scale(:)

   self%label = label

   if (param%xtbml_geometry) then
      allocate(self%geom)
      call new_xtbml_geometry_features(self%geom, mol, error)
      if (allocated(error)) return
   end if

   if (param%xtbml_density) then
      allocate(self%dens)
      call new_xtbml_density_features(self%dens, param%xtbml_tensor)
   end if

   if (param%xtbml_orbital_energy) then
      allocate(self%orb)
      call new_xtbml_orbital_features(self%orb)
   end if

   if (param%xtbml_energy) then
      allocate(self%energy)
      call new_xtbml_energy_features(self%energy)
   end if

   if (param%xtbml_convolution) then
      allocate(self%conv)
      if (allocated(param%xtbml_a)) then
         rcov_scale = param%xtbml_a
      else 
         rcov_scale = (/1.0_wp/)
      end if

      call new_xtbml_convolution(self%conv, mol, rcov_scale, error)
   end if

end subroutine new_xtbml_features

subroutine compute(self, mol, wfn, ints, calc, caches, ctx, timer, &
   & prlevel, dict)
   !> Instance of the xTB-ML features post-processing
   class(xtbml_type),intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction strcuture data
   type(wavefunction_type), intent(in) :: wfn
   !> Integral container
   type(integral_type), intent(in) :: ints
   !> Single-point calculator conatiner
   type(xtb_calculator), intent(in) :: calc
   !> Context container for writing to stdout
   type(context_type), intent(inout) :: ctx
   !> Cache list for storing caches of various interactions
   type(cache_list), intent(inout) :: caches
   !> Timer instance
   type(timer_type), intent(inout) :: timer
   !> Print level
   integer, intent(in) :: prlevel
   !> Dictionary to store the features
   type(double_dictionary_type), intent(inout) :: dict

   type(xtbml_cache) :: mlcache
   type(error_type), allocatable :: error
   integer :: n_features

   call timer%push("xtbML")

   n_features = 0

   if (allocated(self%geom)) then
      call timer%push("geometry")
      call self%geom%compute_features(mol, wfn, ints, calc, caches, &
         & mlcache, dict, n_features, error)
      if (allocated(error)) then
         call ctx%set_error(error)
         return
      end if
      call timer%pop()
   end if

   if (allocated(self%dens)) then
      call timer%push("density")
      call self%dens%compute_features(mol, wfn, ints, calc, caches, &
         & mlcache, dict, n_features, error)
      if (allocated(error)) then
         call ctx%set_error(error)
         return
      end if
      call timer%pop()
   end if

   if (allocated(self%orb)) then
      call timer%push("orbital energy")
      call self%orb%compute_features(mol, wfn, ints, calc, caches, &
         & mlcache, dict, n_features, error)
      if (allocated(error)) then
         call ctx%set_error(error)
         return
      end if
      call timer%pop()
   end if

   if (allocated(self%energy)) then
      call timer%push("energy")
      call self%energy%compute_features(mol, wfn, ints, calc, caches, &
         & mlcache, dict, n_features, error)
      if (allocated(error)) then
         call ctx%set_error(error)
         return
      end if
      call timer%pop()
   end if

   if (allocated(self%conv)) then

      call timer%push("convolution")
      call self%conv%compute_kernel(mol, mlcache)
      call timer%pop()

      if (allocated(self%geom)) then
         call timer%push("geometry convolution")
         call self%geom%compute_extended(mol, wfn, ints, calc, caches, &
            & mlcache, self%conv, dict, n_features, error)
         if (allocated(error)) then
            call ctx%set_error(error)
            return
         end if
         call timer%pop()
      end if

      if (allocated(self%dens)) then
         call timer%push("density convolution")
         call self%dens%compute_extended(mol, wfn, ints, calc, caches, &
            & mlcache, self%conv, dict, n_features, error)
         if (allocated(error)) then
            call ctx%set_error(error)
            return
         end if
         call timer%pop()
      end if

      if (allocated(self%orb)) then
         call timer%push("orbital energy convolution")
         call self%orb%compute_extended(mol, wfn, ints, calc, caches, &
            & mlcache, self%conv, dict, n_features, error)
         if (allocated(error)) then
            call ctx%set_error(error)
            return
         end if
         call timer%pop()
      end if
   end if
   
   ! Store the number of xTB-ML features also in the dictionary
   call dict%add_entry("n_features", [real(n_features, wp)])

   call timer%pop()

end subroutine compute

pure function info(self, verbosity, indent) result(str)
   !> Instance of the xTB-ML features post-processing
   class(xtbml_type), intent(in) :: self
   !> Verbosity level
   integer, intent(in) :: verbosity
   !> Indentation level
   character(len=*), intent(in) :: indent
   !> Information on the container
   character(len=:), allocatable :: str

   character(len=*), parameter :: nl = new_line('a')
   character(len=:), allocatable :: category_indent

   if (allocated(self%label)) then
      str = self%label
   else
      str = "Unknown"
   end if

   category_indent = indent // " * "

   if (allocated(self%geom)) then
      str = str // nl // self%geom%info(verbosity, category_indent)
   end if

   if (allocated(self%dens)) then
      str = str // nl // self%dens%info(verbosity, category_indent)
   end if

   if (allocated(self%orb)) then
      str = str // nl // self%orb%info(verbosity, category_indent)
   end if

   if (allocated(self%energy)) then
      str = str // nl // self%energy%info(verbosity, category_indent)
   end if

   if (allocated(self%conv)) then
      str = str // nl // self%conv%info(verbosity, category_indent)
   end if

end function info

subroutine print_timer(self, timer, prlevel, ctx)
   !> Instance of the xTB-ML features post-processing
   class(xtbml_type), intent(in) :: self
   !> Timer instance
   type(timer_type), intent(in) :: timer
   !> Print level
   integer, intent(in) :: prlevel
   !> Context container for writing to stdout
   type(context_type), intent(inout) :: ctx

   real(wp) :: ttime, stime
   integer :: it
   character(len=*), parameter :: labels(*) = [character(len=20):: &
   & "geometry", "density", "orbital energy", "energy", &
   & "geometry convolution", "density convolution", &
   & "orbital energy convolution"]

   if (prlevel > 2) then
      call ctx%message("ML features timing details:")
      ttime = timer%get("xtbML")
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

end module tblite_post_processing_xtbml
