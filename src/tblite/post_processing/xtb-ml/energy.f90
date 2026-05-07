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

!> @file tblite/post-processing/xtb-ml/energy.f90
!> Energy-based xTB-ML features
module tblite_post_processing_xtbml_energy
   use, intrinsic :: iso_fortran_env, only : error_unit
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_container, only : container_type, container_cache
   use tblite_container_list, only : cache_list
   use tblite_disp_d3, only : d3_dispersion, new_d3_dispersion
   use tblite_disp_d4, only : d4_dispersion, new_d4_dispersion
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_integral_type, only : integral_type
   use tblite_post_processing_xtbml_cache, only : xtbml_cache
   use tblite_post_processing_xtbml_convolution, only : xtbml_convolution
   use tblite_post_processing_xtbml_features, only : xtbml_feature_type
   use tblite_scf_iterator, only : get_electronic_energy, reduce
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_coulomb, only : tb_coulomb
   implicit none
   private

   public :: xtbml_energy_features, new_xtbml_energy_features

   !> Energy-based xTB-ML features
   type, extends(xtbml_feature_type) :: xtbml_energy_features

   contains
      !> Compute xtb-ML energy-based features
      procedure :: compute_features
      !> Compute extended xtb-ML energy-based features with convolution
      procedure :: compute_extended
   end type xtbml_energy_features

   character(len=*), parameter :: label = "energy-based features"

contains


!> Factory for new energy-based xTB-ML feature container
subroutine new_xtbml_energy_features(self)
   !> Instance of the xTB-ML energy features
   class(xtbml_energy_features), intent(out) :: self

   self%label = label

end subroutine new_xtbml_energy_features


!> Compute energy-based features
!> all contributions to the molecular energy are computed here
!> and stored in the feature container
subroutine compute_features(self, mol, wfn, ints, calc, caches, mlcache, &
   & dict, n_features)
   !> Instance of the xTB-ML energy features
   class(xtbml_energy_features), intent(in) :: self
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

   type(container_cache), allocatable :: cache
   type(d3_dispersion), allocatable :: d3
   type(d4_dispersion), allocatable :: d4
   type(error_type), allocatable :: error
   class(container_type), allocatable :: cont
   real(wp), allocatable :: tmp_energy(:), e_ao(:), e_disp_tot(:), e_disp_ATM(:), tot_energy(:)

   allocate(e_ao(calc%bas%nao), source=0.0_wp)
   allocate(tmp_energy(mol%nat), tot_energy(mol%nat),source=0.0_wp)
   call get_electronic_energy(ints%hamiltonian, wfn%density, e_ao)
   call reduce(tmp_energy, e_ao, calc%bas%ao2at)
   call dict%add_entry("E_eht", tmp_energy)
   n_features = n_features + 1
   tot_energy = tmp_energy
   tmp_energy = 0.0_wp

   if (allocated(calc%repulsion)) then
      associate(cont => calc%repulsion)
         if (.not. allocated(cache)) allocate(cache)
         call move_alloc(caches%list(1)%raw, cache%raw)
         cache = caches%list(1)
         call cont%update(mol, cache)
         call cont%get_engrad(mol, cache, tmp_energy)
         call dict%add_entry("E_rep", tmp_energy)
         n_features = n_features + 1
      end associate
   end if

   if (allocated(cache%raw)) deallocate(cache%raw)
   tot_energy = tot_energy + tmp_energy
   tmp_energy = 0.0_wp
   if (allocated(calc%coulomb)) then
      call move_alloc(caches%list(2)%raw, cache%raw)
      associate(cont => calc%coulomb)
         call cont%update(mol, cache)
         if (allocated(cont%es2)) then
            call cont%es2%update(mol, cache)
            call cont%es2%get_energy(mol, cache, wfn, tmp_energy)
         end if
         if (allocated(cont%es3)) then
            call cont%es3%update(mol, cache)
            call cont%es3%get_energy(mol, cache, wfn, tmp_energy)
         end if
         tot_energy = tot_energy + tmp_energy
         call dict%add_entry("E_ies_ixc", tmp_energy)
         n_features = n_features + 1
         if (allocated(cont%aes2)) then
            tmp_energy = 0.0_wp
            call cont%aes2%get_energy_axc(mol, wfn, tmp_energy)
            call dict%add_entry("E_axc", tmp_energy)
            n_features = n_features + 1
            tot_energy = tot_energy + tmp_energy
            tmp_energy = 0.0_wp
            call cont%aes2%get_energy_aes(mol, cache, wfn, tmp_energy)
            call dict%add_entry("E_aes", tmp_energy)
            n_features = n_features + 1
            tot_energy = tot_energy + tmp_energy
         end if
      end associate
   end if

   if (allocated(cache%raw)) deallocate(cache%raw)
   tmp_energy = 0.0_wp
   if (allocated(calc%halogen)) then
      call move_alloc(caches%list(3)%raw, cache%raw)
      associate(cont => calc%halogen)
         call cont%update(mol, cache)
         call cont%get_engrad(mol, cache, tmp_energy)
         call dict%add_entry("E_hx", tmp_energy)
         n_features = n_features + 1
      end associate
   end if

   if (allocated(cache%raw)) deallocate(cache%raw)
   tot_energy = tot_energy + tmp_energy
   tmp_energy = 0.0_wp
   if (allocated(calc%dispersion)) then
      call move_alloc(caches%list(4)%raw, cache%raw)
      associate(cont => calc%dispersion)
         select type(cont)
         type is (d3_dispersion)
            allocate(e_disp_tot(mol%nat), e_disp_ATM(mol%nat), source=0.0_wp)
            call cont%update(mol, cache)
            call cont%get_engrad(mol, cache, e_disp_tot)

            allocate(d3)
            call new_d3_dispersion(d3, mol, s6=0.0_wp, s8=0.0_wp, a1=cont%param%a1, &
               & a2=cont%param%a2, s9=cont%param%s9, error=error, &
               & disp2_width=cont%cutoff%width2, disp3_width=cont%cutoff%width3)
            if(allocated(error)) then
               write(error_unit, '("[Error]:", 1x, a)') error%message
               error stop
            end if
            call d3%update(mol, cache)
            call d3%get_engrad(mol, cache, e_disp_ATM)
            call dict%add_entry("E_disp2", e_disp_tot-e_disp_ATM)
            call dict%add_entry("E_disp3", e_disp_ATM)
            n_features = n_features + 2
            tmp_energy = e_disp_tot
         type is (d4_dispersion)
            allocate(e_disp_tot(mol%nat), e_disp_ATM(mol%nat), source=0.0_wp)
            call cont%update(mol, cache)
            call cont%get_engrad(mol, cache, e_disp_tot)
            call cont%get_energy(mol, cache, wfn, e_disp_tot)

            allocate(d4)
            call new_d4_dispersion(d4, mol, s6=0.0_wp, s8=0.0_wp, a1=cont%param%a1, &
               & a2=cont%param%a2, s9=cont%param%s9, error=error, &
               & disp2_width=cont%cutoff%width2, disp3_width=cont%cutoff%width3)
            if(allocated(error)) then
               write(error_unit, '("[Error]:", 1x, a)') error%message
               error stop
            end if
            call d4%update(mol, cache)
            call d4%get_engrad(mol, cache, e_disp_ATM)
            call dict%add_entry("E_disp2", e_disp_tot-e_disp_ATM)
            call dict%add_entry("E_disp3", e_disp_ATM)
            n_features = n_features + 2
            tmp_energy = e_disp_tot
         end select
      end associate
   end if

   if (allocated(cache%raw)) deallocate(cache%raw)
   tot_energy = tot_energy + tmp_energy
   tmp_energy = 0.0_wp
   if (allocated(calc%interactions)) then
      call move_alloc(caches%list(5)%raw, cache%raw)
      associate(cont => calc%interactions)
         call cont%get_energy(mol, cache, wfn, tmp_energy)
         call cont%update(mol, cache)
         call cont%get_engrad(mol, cache, tmp_energy)
         call dict%add_entry(cont%info(0, ""), tmp_energy)
      end associate
   end if

   if (allocated(cache%raw)) deallocate(cache%raw)
   deallocate(cache)
   tot_energy = tot_energy + tmp_energy
   call dict%add_entry("E_tot", tot_energy)
   call dict%add_entry("w_tot", tot_energy/sum(tot_energy))
   n_features = n_features + 2

end subroutine compute_features

!> Compute extended energy-based features, empty for now
subroutine compute_extended(self, mol, wfn, ints, calc, caches, mlcache, &
   & convolution, dict, n_features)
   !> Instance of the xTB-ML energy features
   class(xtbml_energy_features), intent(in) :: self
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
end subroutine compute_extended

end module tblite_post_processing_xtbml_energy
