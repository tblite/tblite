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

!> @file tblite/post-processing/xtb-ml/orbital_energy.f90
!> Orbital energy based xtbml features
module tblite_post_processing_xtbml_orbital
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_convert, only : autoev
   use tblite_container_list, only : cache_list
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_integral_type, only : integral_type
   use tblite_post_processing_xtbml_atomic_frontier, only : atomic_frontier_orbitals
   use tblite_post_processing_xtbml_cache, only : xtbml_cache
   use tblite_post_processing_xtbml_convolution, only : xtbml_convolution
   use tblite_post_processing_xtbml_features, only : xtbml_feature_type
   use tblite_output_format, only : format_string
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: new_xtbml_orbital_features, xtbml_orbital_features

   !> Orbital energy-based xTB-ML features
   type, extends(xtbml_feature_type) :: xtbml_orbital_features

   contains
      !> Compute xTB-ML density-based features
      procedure :: compute_features
      !> Compute extended xTB-ML density-based features with convolution
      procedure :: compute_extended
   end type xtbml_orbital_features

   character(len=*), parameter :: label = "orbital energy-based features"

contains


subroutine new_xtbml_orbital_features(self)
   !> Instance of the xTB-ML orbital energy features
   class(xtbml_orbital_features) :: self

   self%label = label

end subroutine new_xtbml_orbital_features


subroutine compute_features(self, mol, wfn, ints, calc, caches, mlcache, &
   & dict, n_features)
   !> Instance of the xTB-ML orbital energy features
   class(xtbml_orbital_features), intent(in) :: self
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

   integer :: nspin, spin
   character(len=6), allocatable :: spin_label(:)
   real(wp), allocatable :: response(:), ehoao(:), eluao(:)

   ! Two spin channels for unrestricted or restricted open-shell calculations
   if (wfn%nuhf > 0 .or. wfn%nspin > 1) then
      nspin = 2
   else
      nspin = 1
   end if

   ! Select the spin label
   if (nspin > 1) then
      spin_label = [character(len=6) :: "_alpha", "_beta"]
   else
      spin_label = [""]
   end if

   allocate(response(mol%nat), ehoao(mol%nat), eluao(mol%nat))
   if (.not. allocated(mlcache%egap)) allocate(mlcache%egap(mol%nat, nspin))
   if (.not. allocated(mlcache%chempot)) allocate(mlcache%chempot(mol%nat, nspin))

   do spin = 1, nspin
      if (wfn%nspin > 1) then
         ! Unrestricted case
         call atomic_frontier_orbitals(wfn%focc(:, spin), wfn%emo(:, spin)*autoev, &
            calc%bas%ao2at, wfn%coeff(:, :, spin), ints%overlap, &
            response, mlcache%egap(:, spin), mlcache%chempot(:, spin), ehoao, eluao)
      else
         ! Restricted (or restricted open-shell) case
         call atomic_frontier_orbitals(wfn%focc(:, spin), wfn%emo(:, 1)*autoev, &
            calc%bas%ao2at, wfn%coeff(:, :, 1), ints%overlap, &
            response, mlcache%egap(:, spin), mlcache%chempot(:, spin), ehoao, eluao)
      endif

      ! Store features in the dictionary
      call dict%add_entry(trim("response"//spin_label(spin)), response)
      call dict%add_entry(trim("gap"//spin_label(spin)), mlcache%egap(:, spin))
      call dict%add_entry(trim("chem_pot"//spin_label(spin)), mlcache%chempot(:, spin))
      call dict%add_entry(trim("HOAO"//spin_label(spin)), ehoao)
      call dict%add_entry(trim("LUAO"//spin_label(spin)), eluao)

      ! Count number of features
      n_features = n_features + 5
   end do

end subroutine compute_features

subroutine compute_extended(self, mol, wfn, ints, calc, caches, mlcache, &
   & convolution, dict, n_features)
   !> Instance of the xTB-ML orbital energy features
   class(xtbml_orbital_features), intent(in) :: self
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

   integer :: isc, nspin, spin
   character(len=:), allocatable :: a_label, conv_label, spin_label(:)
   real(wp), allocatable :: ext_chempot(:, :), ext_egap(:, :), tmp_ext_egap(:, :), tmp_ext_chempot(:, :)
   real(wp), allocatable :: ext_ehoao(:, :), ext_eluao(:, :)

   nspin = size(mlcache%chempot, dim=2)

   ! Select the spin label
   if (nspin > 1) then
      spin_label = [character(len=6) :: "_alpha", "_beta"]
   else
      spin_label = [""]
   end if

   allocate(ext_chempot(mol%nat, convolution%nscale))
   allocate(ext_egap(mol%nat, convolution%nscale))
   allocate(ext_ehoao(mol%nat, convolution%nscale))
   allocate(ext_eluao(mol%nat, convolution%nscale))

   allocate(tmp_ext_egap(mol%nat, convolution%nscale))
   allocate(tmp_ext_chempot(mol%nat, convolution%nscale))

   do spin = 1, nspin

      ! Compute convolution of the effective Fermi level per atom
      call convolve_scalar(mol, mlcache%chempot(:, spin), convolution%nscale, &
         & mlcache%conv_kernel, mlcache%conv_cn, ext_chempot)

      ! Compute convolution of the effective atomic HOMO-LUMO gap
      call convolve_scalar(mol, mlcache%egap(:, spin), convolution%nscale, &
         & mlcache%conv_kernel, mlcache%conv_cn, ext_egap)

      ! Compute the highest occupied atomic orbital based on the extended features
      call get_ehoao_ext(ext_chempot, ext_egap, ext_ehoao)

      ! Compute the lowest unoccupied atomic orbital based on the extended features
      call get_eluao_ext(ext_chempot, ext_egap, ext_eluao)

      do isc = 1, convolution%nscale
         a_label = "_"//trim(adjustl(format_string(convolution%rcov_scale(isc), '(f12.2)')))
         if (a_label .eq. "_1.00") a_label = ""
         conv_label = trim(spin_label(spin))//a_label

         ! Add convoluted features to the dictionary
         call dict%add_entry(trim("ext_gap"//conv_label), ext_egap(:, isc))
         call dict%add_entry(trim("ext_chem_pot"//conv_label), ext_chempot(:, isc))
         call dict%add_entry(trim("ext_HOAO"//conv_label), ext_ehoao(:, isc))
         call dict%add_entry(trim("ext_LUAO"//conv_label), ext_eluao(:, isc))

         ! Count number of features
         n_features = n_features + 4
      end do
   end do

end subroutine compute_extended


!> Convolution for scalar atom-resolved properties
!> with rescaling based on the target coordination number
subroutine convolve_scalar(mol, atom_prop, nscale, kernel, cn, ext_prop)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Atom-resolved property data
   real(wp), intent(in) :: atom_prop(:)
   !> Number of convolution length scales
   integer, intent(in) :: nscale
   !> Convolution kernel
   real(wp), intent(in) :: kernel(:, :, :)
   !> Convolution coordination number
   real(wp), intent(in) :: cn(:, :)
   !> Convolved atom-resolved property data
   real(wp), intent(out) :: ext_prop(:, :)

   integer :: iat, jat, ksc
   real(wp) :: tmp

   ext_prop(:, :) = 0.0_wp

   !$omp parallel do default(none) collapse(2) schedule(runtime) &
   !$omp shared(mol, atom_prop, nscale, kernel, cn, ext_prop) &
   !$omp private(iat, jat, ksc, tmp)
   do ksc = 1, nscale
      do iat = 1, mol%nat
         tmp = 0.0_wp
         do jat = 1, mol%nat
            ! Convolution with rescaling based on target coordination number
            tmp = tmp + atom_prop(jat) / kernel(iat, jat, ksc)
         end do
         ext_prop(iat, ksc) = tmp / (cn(iat, ksc) + 1.0_wp)
      end do
   end do

end subroutine convolve_scalar


!> Compute the extended highest occupied atomic orbital
subroutine get_ehoao_ext(ext_chempot, ext_egap, ext_ehoao)
   !> Extended chemical potential
   real(wp), intent(in) :: ext_chempot(:, :)
   !> Extended gap 
   real(wp), intent(in) :: ext_egap(:, :)
   !> Extended highest occupied atomic orbital
   real(wp), intent(out) :: ext_ehoao(:, :)

   integer :: iat, ksc

   do ksc = 1, size(ext_chempot, 2)
      do iat = 1, size(ext_chempot, 1)
         ext_ehoao(iat, ksc) = ext_chempot(iat, ksc) - ext_egap(iat, ksc) / 2.0_wp
      end do
   end do

end subroutine get_ehoao_ext


!> Compute the extended lowest unoccupied atomic orbital
subroutine get_eluao_ext(ext_chempot, ext_egap, ext_eluao)
   !> Extended chemical potential
   real(wp), intent(in) :: ext_chempot(:, :)
   !> Extended gap 
   real(wp), intent(in) :: ext_egap(:, :)
   !> Extended lowest unoccupied atomic orbital
   real(wp), intent(out) :: ext_eluao(:, :)

   integer :: iat, ksc

   do ksc = 1, size(ext_chempot, 2)
      do iat = 1, size(ext_chempot, 1)
         ext_eluao(iat, ksc) = ext_chempot(iat, ksc) + ext_egap(iat, ksc) / 2.0_wp
      end do
   end do

end subroutine get_eluao_ext

end module tblite_post_processing_xtbml_orbital
