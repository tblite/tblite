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

!> @file tblite/post-processing/xtb-ml/convolution.f90
!> Coordination number based convolution
module tblite_post_processing_xtbml_convolution
   use mctc_data_covrad, only : get_covalent_rad
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use mctc_ncoord, only : ncoord_type, new_ncoord, cn_count
   use tblite_post_processing_xtbml_cache, only : xtbml_cache
   implicit none
   private

   public :: xtbml_convolution, new_xtbml_convolution

   !> Wrapped coordination number type for the creation of lists
   type :: ncoord_container
      !> Actual coordination number
      class(ncoord_type), allocatable :: raw
   end type ncoord_container

   !> Convolution xTB-ML features
   type :: xtbml_convolution
      !> Covalent radii
      real(wp), allocatable :: rcov(:)
      !> Linear scaling factor for the covalent radii to generate
      !> coordination numbers of different length scales
      real(wp), allocatable :: rcov_scale(:)
      !> Number of different length scales for the convolution
      integer :: nscale
      !> Label for the convolution container
      character(len=:), allocatable :: label
      !> List of coordination number objects for different length scales
      type(ncoord_container), allocatable :: ncoord(:)
   contains
      !> Update the xTB-ML feature cache with the convolution kernel
      procedure :: compute_kernel
      !> Populate the convolution kernel
      procedure, private :: populate_kernel
      !> Write the xtb-ML convolution information
      procedure :: info
   end type xtbml_convolution

   character(len=*), parameter :: label = "CN-based convolution"

   ! Coordination number exponent
   real(wp), parameter :: kcn = 16.0_wp

contains

subroutine new_xtbml_convolution(self, mol, rcov_scale, error)
   !> Instance of the xTB-ML convolution features
   class(xtbml_convolution), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Scaling factors for the covalent radii to generate
   !> coordination numbers of different length scales
   real(wp), intent(in) :: rcov_scale(:)
   !> Error handling
   type(error_type), intent(inout) , allocatable:: error

   integer :: isc

   self%label = label

   ! Base covalent radii for convolution
   allocate (self%rcov(mol%nid))
   self%rcov(:) = get_covalent_rad(mol%num)

   ! Covalent radii scaling factors for different length scales
   self%rcov_scale = rcov_scale
   self%nscale = size(self%rcov_scale)

   ! Add DFT-D3-like coordination numbers with scaled radii
   allocate(self%ncoord(self%nscale))
   do isc = 1, self%nscale
      call new_ncoord(self%ncoord(isc)%raw, mol, cn_count%exp, error, kcn=kcn, &
         & rcov=self%rcov*self%rcov_scale(isc))
      if (allocated(error)) return
   end do

end subroutine new_xtbml_convolution


!> Update xTB-ML feature cache with the convolution kernel
subroutine compute_kernel(self, mol, mlcache)
   !> Convolution container
   class(xtbml_convolution), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cache for xTB-ML features
   type(xtbml_cache), intent(inout) :: mlcache

   integer :: isc

   if (.not. allocated(mlcache%conv_kernel)) then
      allocate(mlcache%conv_kernel(mol%nat, mol%nat, self%nscale), source=0.0_wp)
   end if
   call self%populate_kernel(mol, mlcache%conv_kernel)

   if (.not. allocated(mlcache%conv_cn)) then
      allocate(mlcache%conv_cn(mol%nat, self%nscale), source=0.0_wp)
   end if
   do isc = 1, self%nscale
      call self%ncoord(isc)%raw%get_cn(mol, mlcache%conv_cn(:, isc))
   end do

end subroutine compute_kernel


!> Calculate the actual convolution kernel
subroutine populate_kernel(self, mol, kernel)
   !> Convolution container
   class(xtbml_convolution), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Convolution kernel
   real(wp), intent(out) :: kernel(:, :, :)

   integer :: iat, jat, izp, jzp, ksc
   real(wp) :: rij(3), r, count

   !$omp parallel do default(none) collapse(3) &
   !$omp shared(self, mol, kernel) &
   !$omp private(iat, jat, izp, jzp, ksc, rij, r, count)
   do ksc = 1, self%nscale
      do iat = 1, mol%nat
         do jat = 1, mol%nat
            if (iat /= jat) then
               izp = mol%id(iat)
               jzp = mol%id(jat)
               rij = mol%xyz(:, iat) - mol%xyz(:, jat)
               r = sqrt(sum(rij**2))

               ! Coordination number counting function
               count = self%ncoord(ksc)%raw%ncoord_count(izp, jzp, r)
               ! Inverse of the counting function as the kernel value
               kernel(iat, jat, ksc) = 1.0_wp / count
            else
               kernel(iat, jat, ksc) = 1.0_wp
            end if
         end do
      end do
   end do

end subroutine populate_kernel


!> Information on the container
pure function info(self, verbosity, indent) result(str)
   !> Instance of the interaction container
   class(xtbml_convolution), intent(in) :: self
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

end module tblite_post_processing_xtbml_convolution
