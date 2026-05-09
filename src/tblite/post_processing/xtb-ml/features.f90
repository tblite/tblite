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

!> @file tblite/post-processing/xtb-ml/features.f90
!> Abstract base class for xtb based features
module tblite_post_processing_xtbml_features
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use tblite_container_list, only : cache_list
   use tblite_context , only : context_type
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_integral_type, only : integral_type
   use tblite_post_processing_xtbml_cache, only : xtbml_cache
   use tblite_post_processing_xtbml_convolution, only : xtbml_convolution
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   !> Abstract base class for xTB-ML features
   type, public, abstract :: xtbml_feature_type
      !> Name of the feature
      character(len=:), allocatable :: label
   contains
      !> Compute xtb-ML features
      procedure(compute_features), deferred :: compute_features
      !> Compute xtb-ML features with convolution
      procedure(compute_extended), deferred :: compute_extended
      !> Write the xtb-ML feature information
      procedure :: info
   end type xtbml_feature_type

   character(len=*), parameter :: label = "General feature class"

   abstract interface
      subroutine compute_features(self, mol, wfn, ints, calc, caches, &
         & mlcache, dict, n_features, error)
         import :: wp, wavefunction_type, structure_type, integral_type, &
            & xtb_calculator, cache_list, context_type, xtbml_feature_type, &
            & xtbml_cache, double_dictionary_type, error_type
         !> Instance of a xTB-ML feature
         class(xtbml_feature_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Wavefunction strcuture data
         type(wavefunction_type), intent(in) :: wfn
         !> Integral container
         type(integral_type), intent(in) :: ints
         !> Single-point calculator
         type(xtb_calculator), intent(in) :: calc
         !> Cache list for storing caches of various interactions
         type(cache_list), intent(inout) :: caches
         !> Cache for xTB-ML features
         type(xtbml_cache), intent(inout) :: mlcache
         !> Dictionary for storing results
         type(double_dictionary_type), intent(inout) :: dict
         !> Number of features
         integer, intent(inout) :: n_features
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_features

      subroutine compute_extended(self, mol, wfn, ints, calc, caches, &
         & mlcache, convolution, dict, n_features, error)
         import :: wp, wavefunction_type, structure_type, integral_type, &
            & xtb_calculator, cache_list, context_type, xtbml_feature_type, &
            & xtbml_cache, xtbml_convolution, double_dictionary_type, error_type
         !> Instance of a xTB-ML feature
         class(xtbml_feature_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Wavefunction strcuture data
         type(wavefunction_type), intent(in) :: wfn
         !> Integral container
         type(integral_type), intent(in) :: ints
         !> Single-point calculator
         type(xtb_calculator), intent(in) :: calc
         !> Cache list for storing caches of various interactions
         type(cache_list), intent(inout) :: caches
         !> Cache for xTB-ML features
         type(xtbml_cache), intent(inout) :: mlcache
         !> Convolution container
         type(xtbml_convolution), intent(in) :: convolution
         !> Dictionary for storing results
         type(double_dictionary_type), intent(inout) :: dict
         !> Number of features
         integer, intent(inout) :: n_features
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_extended

   end interface

contains

pure function info(self, verbosity, indent) result(str)
   !> Instance of a xTB-ML feature
   class(xtbml_feature_type), intent(in) :: self
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

end module tblite_post_processing_xtbml_features
