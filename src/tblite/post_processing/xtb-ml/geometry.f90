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

!> @file tblite/post-processing/xtb-ml/geometry.f90
!> Geometry based xTB-ML features
module tblite_post_processing_xtbml_geometry
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use mctc_ncoord, only : ncoord_type, new_ncoord, cn_count
   use tblite_container_list, only : cache_list
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_integral_type, only : integral_type
   use tblite_output_format, only : format_string
   use tblite_post_processing_xtbml_cache, only : xtbml_cache
   use tblite_post_processing_xtbml_convolution, only : xtbml_convolution
   use tblite_post_processing_xtbml_features, only : xtbml_feature_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: new_xtbml_geometry_features, xtbml_geometry_features

   !> Geometry-based xTB-ML features
   type, extends(xtbml_feature_type) :: xtbml_geometry_features
      !> Coordination number
      class(ncoord_type), allocatable :: ncoord
   contains
      !> Compute xtb-ML geometry-based features
      procedure :: compute_features
      !> Compute extended xtb-ML geometry-based features with convolution
      procedure :: compute_extended
   end type xtbml_geometry_features

   character(len=23), parameter :: label = "geometry-based features"

contains

!> Factory for new geometry-based xTB-ML feature container
subroutine new_xtbml_geometry_features(self, mol, error)
   !> Instance of the xTB-ML geometry features
   class(xtbml_geometry_features), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), intent(inout) , allocatable:: error

   self%label = label

   call new_ncoord(self%ncoord, mol, cn_count%exp, error)

end subroutine new_xtbml_geometry_features


subroutine compute_features(self, mol, wfn, ints, calc, caches, mlcache, &
   & dict, n_features)
   !> Instance of the xTB-ML geometry features
   class(xtbml_geometry_features), intent(in) :: self
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

   ! Count number of features
   n_features = n_features + 1

   ! Compute the coordination number
   if (.not. allocated(mlcache%cn)) then
      allocate(mlcache%cn(mol%nat))
   end if
   call self%ncoord%get_cn(mol, mlcache%cn)

   call dict%add_entry("CN_A", mlcache%cn)

end subroutine compute_features


subroutine compute_extended(self, mol, wfn, ints, calc, caches, mlcache, &
   & convolution, dict, n_features)
   !> Instance of the xTB-ML geometry features
   class(xtbml_geometry_features), intent(in) :: self
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

   character(len=:), allocatable :: tmp_label
   real(wp), allocatable ::  ext_cn(:, :)
   integer :: isc
   
   allocate(ext_cn(mol%nat, convolution%nscale))
   call convolve_cn(mol, mlcache%cn, convolution%nscale, &
      & mlcache%conv_kernel, ext_cn)
   
   do isc = 1, convolution%nscale
      tmp_label = trim("ext_CN_A"//'_'//adjustl(format_string( &
         & convolution%rcov_scale(isc), '(f12.2)')))
      if (tmp_label .eq. "ext_CN_A_1.00") tmp_label = "ext_CN"
      call dict%add_entry(tmp_label, ext_cn(:, isc))
      ! Count number of features
      n_features = n_features + 1
   end do

end subroutine compute_extended


! Convolution of the scalar coordination numbers for different length scales
subroutine convolve_cn(mol, cn, nscale, kernel, ext_cn)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Coordinated number
   real(wp), intent(in) :: cn(:)
   !> Number of convolution length scales
   integer, intent(in) :: nscale
   !> Convolution kernel
   real(wp), intent(in) :: kernel(:, :, :)
   !> Delta CN
   real(wp), intent(out) :: ext_cn(:, :)

   integer :: iat, jat, ksc
   real(wp) :: result

   ext_cn = 0.0_wp
   !$omp parallel do default(none) collapse(2)&
   !$omp shared(mol, nscale, cn, kernel, ext_cn)&
   !$omp private(iat, jat , ksc, result)
   do ksc = 1, nscale
      do iat = 1, mol%nat
         result = 0.0_wp
         do jat = 1, mol%nat
            if (iat == jat) cycle
            result = result + cn(jat) / kernel(iat, jat, ksc)
         end do
         ext_cn(iat, ksc) = ext_cn(iat, ksc) + result
      end do
   end do
   !$omp end parallel do

end subroutine convolve_cn

end module tblite_post_processing_xtbml_geometry
