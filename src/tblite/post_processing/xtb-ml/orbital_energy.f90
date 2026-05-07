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
   use mctc_io_convert, only : evtoau, autoev
   use tblite_basis_type, only : basis_type
   use tblite_blas, only : symm, gemm, dot
   use tblite_container_list, only : cache_list
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_integral_type, only : integral_type
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

   !> Occupation cutoff for fractional occupations
   real(wp), parameter :: occ_cutoff = 1.0e-4_wp

   !> Threshold for missing occupied/virtual atomic population
   real(wp), parameter :: population_cutoff = 1.0e-12_wp

   !> Damping for the resonse function in Hartree
   real(wp), parameter :: damp = 0.5_wp * evtoau

   !> Effective infinite frontier energy
   real(wp), parameter :: near_infty = 1.0e100_wp

   !> Regularizer for divisions by response-like quantities
   real(wp), parameter :: regularizer = 1.0e-14_wp

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
      ! Compute atomic response, effective HOMO-LUMO gap, chemical potential, HOAO and LUAO
      call atomic_frontier_orbitals(spin, mol, calc%bas, wfn, ints%overlap, &
         & response, mlcache%egap(:, spin), mlcache%chempot(:, spin), ehoao, eluao)

      ! Store features in the dictionary
      call dict%add_entry(trim("response"//spin_label(spin)), response * evtoau**2)
      call dict%add_entry(trim("gap"//spin_label(spin)), mlcache%egap(:, spin) * autoev)
      call dict%add_entry(trim("chem_pot"//spin_label(spin)), mlcache%chempot(:, spin) * autoev)
      call dict%add_entry(trim("HOAO"//spin_label(spin)), ehoao * autoev)
      call dict%add_entry(trim("LUAO"//spin_label(spin)), eluao * autoev)

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
   real(wp), allocatable :: ext_chempot(:, :), ext_egap(:, :)
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
         call dict%add_entry(trim("ext_gap"//conv_label), ext_egap(:, isc) * autoev)
         call dict%add_entry(trim("ext_chem_pot"//conv_label), ext_chempot(:, isc) * autoev)
         call dict%add_entry(trim("ext_HOAO"//conv_label), ext_ehoao(:, isc) * autoev)
         call dict%add_entry(trim("ext_LUAO"//conv_label), ext_eluao(:, isc) * autoev)

         ! Count number of features
         n_features = n_features + 4
      end do
   end do

end subroutine compute_extended


!> Compute atomic response, effective H-L gap, chemical potential, HOAO and LUAO
subroutine atomic_frontier_orbitals(spin, mol, bas, wfn, overlap, &
   & response, egap, chempot, ehoao, eluao)
   !> Spin channel
   integer, intent(in) :: spin
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Overlap integral matrix
   real(wp), intent(in) :: overlap(:,:)
   !> Atomic response
   real(wp), intent(out) :: response(:)
   !> Effective atomic H-L gap
   real(wp), intent(out) :: egap(:)
   !> Effective atomic chemical potential
   real(wp), intent(out) :: chempot(:)
   !> Highest occupied atomic orbital
   real(wp), intent(out) :: ehoao(:)
   !> Lowest unoccupied atomic orbital
   real(wp), intent(out) :: eluao(:)

   integer :: ispin, nocc, nvirt
   integer, allocatable :: iocc(:), ivirt(:)
   real(wp), allocatable :: po(:, :), pv(:, :)
   real(wp), allocatable :: po_sum(:), pv_sum(:)

   ! Select spin channel for the coefficients
   if (wfn%nspin > 1) then
      ispin = spin
   else
      ispin = 1
   end if

   ! Select occupied and virtual orbitals based on the occupation numbers
   call get_occ_virt_orbitals(wfn%focc(:, spin), iocc, ivirt)
   nocc = size(iocc)
   nvirt = size(ivirt)

   ! Project occupied and virtual MO populations onto atoms
   allocate(po(nocc, mol%nat), pv(nvirt, mol%nat))
   allocate(po_sum(mol%nat), pv_sum(mol%nat))
   call atomic_mo_projection(bas, iocc, ivirt, wfn%focc(:, spin), &
      & wfn%coeff(:, :, ispin), overlap, po, pv, po_sum, pv_sum)

   ! Accumulate response, chemical potential, inverse-gap, and HOAO/LUAO
   call accumulate_properties(nocc, nvirt, iocc, ivirt, po, pv, po_sum, pv_sum, &
      & wfn%emo(:, ispin), response, egap, chempot, ehoao, eluao)

end subroutine atomic_frontier_orbitals


!> Collect occupied and virtual orbital indices
subroutine get_occ_virt_orbitals(focc, iocc, ivirt)

   !> Occupation numbers
   real(wp), intent(in) :: focc(:)
   !> Active occupied orbitals
   integer, allocatable, intent(out) :: iocc(:)
   !> Active virtual orbitals
   integer, allocatable, intent(out) :: ivirt(:)

   integer :: imo, io, iv

   allocate(iocc(count(focc > occ_cutoff)))
   allocate(ivirt(count(1.0_wp - focc > occ_cutoff)))

   io = 0
   iv = 0
   do imo = 1, size(focc)
      if (focc(imo) > occ_cutoff) then
         io = io + 1
         iocc(io) = imo
      end if

      if (1.0_wp - focc(imo) > occ_cutoff) then
         iv = iv + 1
         ivirt(iv) = imo
      end if
   end do

end subroutine get_occ_virt_orbitals


!> Project occupied and virtual MO populations onto atoms
subroutine atomic_mo_projection(bas, iocc, ivirt, focc, coeff, overlap, po, pv, &
   & po_sum, pv_sum)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Indices of occupied orbitals
   integer, intent(in) :: iocc(:)
   !> Indices of virtual orbitals
   integer, intent(in) :: ivirt(:)
   !> MO occupation numbers
   real(wp), intent(in) :: focc(:)
   !> MO coefficients
   real(wp), contiguous, intent(in) :: coeff(:, :)
   !> Overlap matrix
   real(wp), contiguous, intent(in) :: overlap(:, :)
   !> Occupied atom-resolved MO population
   real(wp), intent(out) :: po(:, :)
   !> Virtual atom-resolved MO population
   real(wp), intent(out) :: pv(:, :)
   !> Total occupied population per atom
   real(wp), intent(out) :: po_sum(:)
   !> Total virtual population per atom
   real(wp), intent(out) :: pv_sum(:)

   integer :: io, iv, imo, iao, iat
   real(wp) :: ps
   real(wp), allocatable :: scoeff(:, :)

   po(:, :) = 0.0_wp
   pv(:, :) = 0.0_wp

   ! Precompute overlap and coefficient contraction for all molecular orbitals
   allocate(scoeff(bas%nao, bas%nao))
   call symm(overlap, coeff, scoeff)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(bas, iocc, focc, coeff, scoeff, po) &
   !$omp private(io, imo, iao, iat, ps)
   do io = 1, size(iocc)
      imo = iocc(io)
      do iao = 1, bas%nao
         iat = bas%ao2at(iao)
         ps = coeff(iao, imo) * scoeff(iao, imo)
         po(io, iat) = po(io, iat) + focc(imo) * ps
      end do
   end do

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(bas, ivirt, focc, coeff, scoeff, pv) &
   !$omp private(iv, imo, iao, iat, ps)
   do iv = 1, size(ivirt)
      imo = ivirt(iv)
      do iao = 1, bas%nao
         iat = bas%ao2at(iao)
         ps = coeff(iao, imo) * scoeff(iao, imo)
         pv(iv, iat) = pv(iv, iat) + (1.0_wp - focc(imo)) * ps
      end do
   end do

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(po, pv, po_sum, pv_sum) private(iat)
   do iat = 1, size(po_sum)
      po_sum(iat) = sum(po(:, iat))
      pv_sum(iat) = sum(pv(:, iat))
   end do

end subroutine atomic_mo_projection


!> Compute response, effective H-L gap, chemical potential, HOAO and LUAO
subroutine accumulate_properties(nocc, nvirt, iocc, ivirt, po, pv, &
      & po_sum, pv_sum, eps, response, egap, chempot, ehoao, eluao)
   !> Number of occupied orbitals
   integer, intent(in) :: nocc
   !> Number of virtual orbitals
   integer, intent(in) :: nvirt
   !> Indices of occupied orbitals
   integer, intent(in) :: iocc(:)
   !> Indices of virtual orbitals
   integer, intent(in) :: ivirt(:)
   !> Occupied atom-resolved MO population
   real(wp), intent(in) :: po(:, :)
   !> Virtual atom-resolved MO population
   real(wp), intent(in) :: pv(:, :)
   !> Total occupied population per atom
   real(wp), intent(in) :: po_sum(:)
   !> Total virtual population per atom
   real(wp), intent(in) :: pv_sum(:)
   !> MO energies
   real(wp), intent(in) :: eps(:)
   !> Atomic response
   real(wp), intent(out) :: response(:)
   !> Effective atomic HOMO-LUMO gap
   real(wp), intent(out) :: egap(:)
   !> Effective atomic chemical potential
   real(wp), intent(out) :: chempot(:)
   !> Highest occupied atomic orbital
   real(wp), intent(out) :: ehoao(:)
   !> Lowest unoccupied atomic orbital
   real(wp), intent(out) :: eluao(:)

   integer :: io, iv, iocc_mo, ivirt_mo, iat, nat
   real(wp) :: delta_eps, sum_eps, tmp
   real(wp), allocatable :: response_kernel(:, :), response_weight(:, :) 
   real(wp), allocatable :: chempot_kernel(:, :), chempot_weight(:, :)
   real(wp), allocatable :: gap_kernel(:, :), gap_weight(:, :)

   response(:) = 0.0_wp
   egap(:) = 0.0_wp
   chempot(:) = 0.0_wp
   ehoao(:) = 0.0_wp
   eluao(:) = 0.0_wp

   nat = size(response)

   ! Generate intermediates if both occupied and virtual populations are present
   if (nocc > 0 .and. nvirt > 0) then
      allocate(response_kernel(nocc, nvirt), response_weight(nocc, nat))
      allocate(chempot_kernel(nocc, nvirt), chempot_weight(nocc, nat))
      allocate(gap_kernel(nocc, nvirt), gap_weight(nocc, nat))

      !$omp parallel do default(none) collapse(2) schedule(runtime) &
      !$omp shared(nocc, nvirt, eps, iocc, ivirt, response_kernel) &
      !$omp shared(chempot_kernel, gap_kernel) &
      !$omp private(io, iv, ivirt_mo, iocc_mo, delta_eps, sum_eps, tmp)
      do iv = 1, nvirt
         do io = 1, nocc
            ivirt_mo = ivirt(iv)
            iocc_mo = iocc(io)

            delta_eps = eps(ivirt_mo) - eps(iocc_mo)
            sum_eps = eps(ivirt_mo) + eps(iocc_mo)

            tmp = 1.0_wp / (delta_eps**2 + damp**2)

            ! Construct kernels for response, chemical potential and gap
            response_kernel(io, iv) = tmp
            chempot_kernel(io, iv) = 0.5_wp * sum_eps * tmp
            gap_kernel(io, iv) = tmp / (delta_eps + damp)
         end do
      end do

      ! Contract virtual populations with the kernels
      call gemm(response_kernel, pv, response_weight)
      call gemm(chempot_kernel, pv, chempot_weight)
      call gemm(gap_kernel, pv, gap_weight)
   end if

   ! Contract with occupied populations and finalize atom-wise features.
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(nat, nocc, nvirt, iocc, ivirt, eps, po, pv, po_sum, pv_sum) &
   !$omp shared(response_weight, chempot_weight, gap_weight) &
   !$omp shared(response, chempot, egap, ehoao, eluao) &
   !$omp private(iat, io, iv, iocc_mo, ivirt_mo, tmp)
   do iat = 1, nat

      if (po_sum(iat) < population_cutoff) then
         ! No occupied population on this atom
         egap(iat) = 0.0_wp
         chempot(iat) = 0.0_wp

         ! Calculate intermediates only based on virtual orbitals
         do iv = 1, nvirt
            ivirt_mo = ivirt(iv)
            tmp = pv(iv, iat) / (eps(ivirt_mo)**2 + damp**2)

            chempot(iat) = chempot(iat) + eps(ivirt_mo) * tmp
            egap(iat) = egap(iat) + tmp / (eps(ivirt_mo) + near_infty + damp)
         end do

         ! Normalize gap with response and invert after regularization
         egap(iat) = 1.0_wp / (egap(iat) + regularizer) - damp

         ! Finite LUAO set to the chemical potential
         eluao(iat) = chempot(iat)
         ehoao(iat) = chempot(iat) - egap(iat)

      else if (pv_sum(iat) < population_cutoff) then
         ! No virtual population on this atom
         egap(iat) = 0.0_wp
         chempot(iat) = 0.0_wp

         ! Calculate intermediates only based on occupied orbitals
         do io = 1, nocc
            iocc_mo = iocc(io)
            tmp = po(io, iat) / (eps(iocc_mo)**2 + damp**2)

            chempot(iat) = chempot(iat) + eps(iocc_mo) * tmp
            egap(iat) = egap(iat) + tmp / (near_infty - eps(iocc_mo) + damp)
         end do

         ! Normalize gap with response and invert after regularization
         egap(iat) = 1.0_wp / (egap(iat) + regularizer) - damp

         ! Finite HOAO set to the chemical potential
         ehoao(iat) = chempot(iat)
         eluao(iat) = chempot(iat) + egap(iat)

      else
         ! Contract kernel weight with occupied populations
         response(iat) = dot(po(:, iat), response_weight(:, iat))
         chempot(iat) = dot(po(:, iat), chempot_weight(:, iat))
         egap(iat) = dot(po(:, iat), gap_weight(:, iat))
      
         ! Normalize gap with response and invert after regularization
         egap(iat) = egap(iat) / (response(iat) + regularizer)
         egap(iat) = 1.0_wp / (egap(iat) + regularizer) - damp

         ! Normalize chemical potential with response
         chempot(iat) = chempot(iat) / (response(iat) + regularizer)

         ! Compute HOAO and LUAO
         ehoao(iat) = chempot(iat) - 0.5_wp * egap(iat)
         eluao(iat) = chempot(iat) + 0.5_wp * egap(iat)
      end if
   end do

end subroutine accumulate_properties


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
