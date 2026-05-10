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

!> @file tblite/post-processing/xtb-ml/density.f90
!> Density-based xTB-ML features using Mulliken partitioning
module tblite_post_processing_xtbml_density
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_container_list, only : cache_list
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_integral_type, only : integral_type
   use tblite_post_processing_xtbml_cache, only : xtbml_cache
   use tblite_post_processing_xtbml_convolution, only : xtbml_convolution
   use tblite_post_processing_xtbml_features, only : xtbml_feature_type
   use tblite_output_format, only : format_string
   use tblite_wavefunction_mulliken, only : get_mulliken_shell_multipoles
   use tblite_wavefunction_spin, only : magnet_to_updown
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: new_xtbml_density_features, xtbml_density_features

   !> Density-based xTB-ML features
   type, extends(xtbml_feature_type) :: xtbml_density_features
      !> Flag to return cartesian components
      logical :: return_xyz
   contains
      !> Compute xTB-ML density-based features
      procedure :: compute_features
      !> Compute extended xTB-ML density-based features with convolution
      procedure :: compute_extended
   end type xtbml_density_features

   character(len=*), parameter :: label = "density-based features"

   ! Labels for the shell angular momenta
   character(len=1), parameter :: shell_labels(0:4) = &
      [character(len=1) :: "s", "p", "d", "f", "g"]

   ! Labels for the dipole moment components
   character(len=1), parameter :: dp_labels(3) = &
      [character(len=1) :: "x", "y", "z"]

   ! Labels for the quadrupole moment components
   character(len=2), parameter :: qp_labels(6) = &
      [character(len=2) :: "xx", "xy", "yy", "xz", "yz", "zz"]

   character(len=0), parameter :: empty = ""

contains


subroutine new_xtbml_density_features(self, return_xyz)
   !> Instance of the xTB-ML density features
   class(xtbml_density_features), intent(out) :: self
   !> Flag to return cartesian components
   logical, intent(in) :: return_xyz

   self%label = label
   self%return_xyz = return_xyz

end subroutine new_xtbml_density_features


subroutine compute_features(self, mol, wfn, ints, calc, caches, mlcache, &
   & dict, n_features, error)
   !> Instance of the xTB-ML density features
   class(xtbml_density_features), intent(in) :: self
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

   integer :: spin, ic
   real(wp), allocatable :: psh(:, :), qat(:, :), dpsh(:, :, :), qpsh(:, :, :)
   real(wp), allocatable :: dpsh_norm(:, :), qpsh_norm(:, :)
   real(wp), allocatable :: dpat_norm(:, :), qpat_norm(:, :)
   character(len=:), allocatable :: spin_label(:)

   ! Select the spin label
   if (wfn%nspin > 1) then
      spin_label = [character(len=6) ::"_alpha", "_beta"]
   else
      spin_label = [""]
   end if

   ! Restore Mulliken alpha/beta shell populations from Mulliken shell charges
   if (.not. allocated(mlcache%psh)) then
      allocate(mlcache%psh(calc%bas%nsh, wfn%nspin), source=-wfn%qsh)
   else
      mlcache%psh = -wfn%qsh
   end if
   mlcache%psh(:, 1) = wfn%n0sh + mlcache%psh(:, 1)
   call magnet_to_updown(mlcache%psh)

   ! Convert atomic charges and magnetization to alpha/beta representation
   if (.not. allocated(mlcache%qat)) then
      allocate(mlcache%qat(mol%nat, wfn%nspin), source=wfn%qat)
   else
      mlcache%qat = wfn%qat
   end if
   call magnet_to_updown(mlcache%qat)

   ! Compute shell multipoles missing in the wavefunction
   if (.not. allocated(mlcache%dpsh)) then
      allocate(mlcache%dpsh(3, calc%bas%nsh, wfn%nspin), source=0.0_wp)
   end if
   if (.not. allocated(mlcache%qpsh)) then
      allocate(mlcache%qpsh(6, calc%bas%nsh, wfn%nspin), source=0.0_wp)
   end if
   call get_mulliken_shell_multipoles(calc%bas, ints%dipole, wfn%density, &
      & mlcache%dpsh)
   call get_mulliken_shell_multipoles(calc%bas, ints%quadrupole, wfn%density, &
      & mlcache%qpsh)

   ! Calculate shell multipole norms
   if (.not. allocated(mlcache%dpsh_norm)) then
      allocate(mlcache%dpsh_norm(calc%bas%nsh, wfn%nspin), source=0.0_wp)
   end if
   if (.not. allocated(mlcache%qpsh_norm)) then
      allocate(mlcache%qpsh_norm(calc%bas%nsh, wfn%nspin), source=0.0_wp)
   end if
   call multipole_norm(mlcache%dpsh, mlcache%qpsh, mlcache%dpsh_norm, &
      & mlcache%qpsh_norm)

   ! Calculate atomic multipole norms
   if (.not. allocated(mlcache%dpat_norm)) then
      allocate(mlcache%dpat_norm(mol%nat, wfn%nspin), source=0.0_wp)
   end if
   if (.not. allocated(mlcache%qpat_norm)) then
      allocate(mlcache%qpat_norm(mol%nat, wfn%nspin), source=0.0_wp)
   end if
   call multipole_norm(wfn%dpat, wfn%qpat, mlcache%dpat_norm, mlcache%qpat_norm)

   ! Store features in dictionary
   do spin = 1, wfn%nspin

      ! Store shell Mulliken populations
      call add_shell_feature(mol, calc%bas, mlcache%psh(:, spin), "p", empty, &
         & spin_label(spin), dict, n_features)

      ! Store shell dipole moment norms
      call add_shell_feature(mol, calc%bas, mlcache%dpsh_norm(:, spin), "dipm", &
         & empty, spin_label(spin), dict, n_features)
      
      ! Optionally store shell dipole moment components
      if (self%return_xyz) then
         do ic = 1, size(dp_labels)
            call add_shell_feature(mol, calc%bas, mlcache%dpsh(ic, :, spin), &
               & "dipm", "_"//trim(dp_labels(ic)), spin_label(spin), dict, &
               & n_features)
         end do
      end if

      ! Store shell quadrupole moment norms
      call add_shell_feature(mol, calc%bas, mlcache%qpsh_norm(:, spin), "qm", &
         & empty, spin_label(spin), dict, n_features)

      ! Optionally store shell quadrupole moment components
      if (self%return_xyz) then
         do ic = 1, size(qp_labels)
            call add_shell_feature(mol, calc%bas, mlcache%qpsh(ic, :, spin), &
               & "qm", "_"//trim(qp_labels(ic)), spin_label(spin), dict, &
               & n_features)
         end do
      end if

      ! Store atomic Mulliken charges
      call add_atom_feature(mol, mlcache%qat(:, spin), "q_A", empty, &
         & spin_label(spin), dict, n_features)

      ! Store atomic dipole moment norms
      call add_atom_feature(mol, mlcache%dpat_norm(:, spin), "dipm_A", empty, &
         & spin_label(spin), dict, n_features)

      ! Optionally store atomic dipole moment components
      if (self%return_xyz) then
         do ic = 1, size(dp_labels)
            call add_atom_feature(mol, wfn%dpat(ic, :, spin), "dipm_A", &
               & "_"//trim(dp_labels(ic)), spin_label(spin), dict, n_features)
         end do
      end if

      ! Store atomic quadrupole moment norms
      call add_atom_feature(mol, mlcache%qpat_norm(:, spin), "qm_A", empty, &
         & spin_label(spin), dict, n_features)

      ! Optionally store atomic quadrupole moment components
      if (self%return_xyz) then
         do ic = 1, size(qp_labels)
            call add_atom_feature(mol, wfn%qpat(ic, :, spin), "qm_A", &
               & "_"//trim(qp_labels(ic)), spin_label(spin), dict, n_features)
         end do
      end if

   end do

end subroutine compute_features


subroutine compute_extended(self, mol, wfn, ints, calc, caches, mlcache, &
   & convolution, dict, n_features, error)
   !> Instance of the xTB-ML density features
   class(xtbml_density_features), intent(in) :: self
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

   integer :: spin, ic, isc
   real(wp), allocatable :: ext_qat(:, :), ext_dpat(:, :, :), ext_qpat(:, :, :)
   real(wp), allocatable :: ext_dpat_norm(:, :), ext_qpat_norm(:, :)
   real(wp), allocatable :: ext_dpat_z(:, :, :), ext_qpat_z(:, :, :)
   real(wp), allocatable :: ext_dpat_z_norm(:, :), ext_qpat_z_norm(:, :)
   real(wp), allocatable :: ext_dpat_p(:, :, :), ext_qpat_p(:, :, :)
   real(wp), allocatable :: ext_dpat_p_norm(:, :), ext_qpat_p_norm(:, :)
   real(wp), allocatable :: pat(:, :)
   character(len=:), allocatable :: a_label, conv_label, spin_label(:)

   ! Select the spin label
   if (wfn%nspin > 1) then
      spin_label = [character(len=6) :: "_alpha", "_beta"]
   else
      spin_label = [""]
   end if

   allocate(ext_qat(mol%nat, convolution%nscale))
   allocate(ext_dpat(3, mol%nat, convolution%nscale))
   allocate(ext_qpat(6, mol%nat, convolution%nscale))
   allocate(ext_dpat_z(3, mol%nat, convolution%nscale))
   allocate(ext_qpat_z(6, mol%nat, convolution%nscale))
   allocate(ext_dpat_p(3, mol%nat, convolution%nscale))
   allocate(ext_qpat_p(6, mol%nat, convolution%nscale))
   allocate(ext_dpat_norm(mol%nat, convolution%nscale))
   allocate(ext_qpat_norm(mol%nat, convolution%nscale))
   allocate(ext_dpat_z_norm(mol%nat, convolution%nscale))
   allocate(ext_qpat_z_norm(mol%nat, convolution%nscale))
   allocate(ext_dpat_p_norm(mol%nat, convolution%nscale))
   allocate(ext_qpat_p_norm(mol%nat, convolution%nscale))

   ! Atomic Mulliken populations
   allocate(pat(mol%nat, wfn%nspin), source=-wfn%qat)
   pat(:, 1) = wfn%n0at + pat(:, 1)
   call magnet_to_updown(pat)

   ! Calculate convolutions and store features in dictionary
   do spin = 1, wfn%nspin

      ! Compute convolution of atomic Mulliken charges
      call convolve_scalar(mol, mlcache%qat(:, spin), convolution%nscale, &
         & mlcache%conv_kernel, mlcache%conv_cn, ext_qat)

      ! Compute convolution of atomic multipole moments and Mulliken charges
      call convolve_multipole(mol, mlcache%qat(:, spin), convolution%nscale, &
         & mlcache%conv_kernel, mlcache%conv_cn, ext_dpat, ext_qpat, &
         & wfn%dpat(:, :, spin), wfn%qpat(:, :, spin))

      ! Compute norms of the extended atomic multipole moment convolution
      call multipole_norm(ext_dpat, ext_qpat, ext_dpat_norm, ext_qpat_norm)

      ! Compute convolution of atomic multipole moments and electronic populations
      call convolve_multipole(mol, -pat(:, spin), convolution%nscale, &
         & mlcache%conv_kernel, mlcache%conv_cn, ext_dpat_p, ext_qpat_p, &
         & wfn%dpat(:, :, spin), wfn%qpat(:, :, spin))

      ! Compute norms of the extended atomic multipole moment convolution
      call multipole_norm(ext_dpat_p, ext_qpat_p, ext_dpat_p_norm, ext_qpat_p_norm)

      ! Compute convolution of atomic multipole moments due to screened nuclear charges
      call convolve_multipole(mol, wfn%n0at, convolution%nscale, &
         & mlcache%conv_kernel, mlcache%conv_cn, ext_dpat_z, ext_qpat_z)

      ! Compute norms of the extended atomic multipole moment convolution
      call multipole_norm(ext_dpat_z, ext_qpat_z, ext_dpat_z_norm, ext_qpat_z_norm)

      do isc = 1, convolution%nscale

         ! Suffix for combined spin and convolution
         a_label = "_"//trim(adjustl(format_string(convolution%rcov_scale(isc), '(f12.2)')))
         if (a_label == "_1.00") a_label = empty
         conv_label = trim(spin_label(spin))//a_label

         ! Store convolved atomic charges
         call add_atom_feature(mol, ext_qat(:, isc), "ext_q_A", empty, &
            & conv_label, dict, n_features)

         ! Store convolved dipole moment norms
         call add_atom_feature(mol, ext_dpat_norm(:, isc), "ext_dipm_A", &
            & empty, conv_label, dict, n_features)

         ! Optionally store convolved dipole moment components
         if (self%return_xyz) then
            do ic = 1, size(dp_labels)
               call add_atom_feature(mol, ext_dpat(ic, :, isc), "ext_dipm_A", &
                  & "_"//trim(dp_labels(ic)), conv_label, dict, n_features)
            end do
         end if

         ! Store convolved quadrupole moment norms
         call add_atom_feature(mol, ext_qpat_norm(:, isc), "ext_qm_A", &
            & empty, conv_label, dict, n_features)

         ! Optionally store convolved quadrupole moment components
         if (self%return_xyz) then
            do ic = 1, size(qp_labels)
               call add_atom_feature(mol, ext_qpat(ic, :, isc), "ext_qm_A", &
                  & "_"//trim(qp_labels(ic)), conv_label, dict, n_features)
            end do
         end if

         ! Store convolved dipole moment (population) norms
         call add_atom_feature(mol, ext_dpat_p_norm(:, isc), "ext_dipm_e", &
            & empty, conv_label, dict, n_features)

         ! Optionally store convolved dipole moment (population) components
         if (self%return_xyz) then
            do ic = 1, size(dp_labels)
               call add_atom_feature(mol, ext_dpat_p(ic, :, isc), "ext_dipm_e", &
                  & "_"//trim(dp_labels(ic)), conv_label, dict, n_features)
            end do
         end if

         ! Store convolved quadrupole moment (population) norms
         call add_atom_feature(mol, ext_qpat_p_norm(:, isc), "ext_qm_e", &
            & empty, conv_label, dict, n_features)

         ! Optionally store convolved quadrupole moment (population) components
         if (self%return_xyz) then
            do ic = 1, size(qp_labels)
               call add_atom_feature(mol, ext_qpat_p(ic, :, isc), "ext_qm_e", &
                  & "_"//trim(qp_labels(ic)), conv_label, dict, n_features)
            end do
         end if

         ! Store convolved dipole moment (nuclear charge) norms
         call add_atom_feature(mol, ext_dpat_z_norm(:, isc), "ext_dipm_Z", &
            & empty, conv_label, dict, n_features)

         ! Optionally store convolved dipole moment (nuclear charge) components
         if (self%return_xyz) then
            do ic = 1, size(dp_labels)
               call add_atom_feature(mol, ext_dpat_z(ic, :, isc), "ext_dipm_Z", &
                  & "_"//trim(dp_labels(ic)), conv_label, dict, n_features)
            end do
         end if

         ! Store convolved quadrupole moment (nuclear charge) norms
         call add_atom_feature(mol, ext_qpat_z_norm(:, isc), "ext_qm_Z", &
            & empty, conv_label, dict, n_features)

         ! Optionally store convolved quadrupole moment (nuclear charge) components
         if (self%return_xyz) then
            do ic = 1, size(qp_labels)
               call add_atom_feature(mol, ext_qpat_z(ic, :, isc), "ext_qm_Z", &
                  & "_"//trim(qp_labels(ic)), conv_label, dict, n_features)
            end do
         end if

      end do
   end do

end subroutine compute_extended


!> Compute the norm for dipole and quadrupole moments
subroutine multipole_norm(dp, qp, dp_norm, qp_norm)
   !> Dipole moment
   real(wp), intent(in) :: dp(:, :, :)
   !> Quadrupole moment
   real(wp), intent(in) :: qp(:, :, :)
   !> Frobenius norm of the dipole moment tensor
   real(wp), intent(out) :: dp_norm(:, :)
   !> Frobenius norm of the quadrupole moment tensor
   real(wp), intent(out) :: qp_norm(:, :)

   dp_norm(:, :) = sqrt(dp(1, :, :)**2 + dp(2, :, :)**2 + dp(3, :, :)**2)

   qp_norm(:, :) = sqrt( &
      & qp(1, :, :)**2 + 2.0_wp * qp(2, :, :)**2 + &
      & qp(3, :, :)**2 + 2.0_wp * qp(4, :, :)**2 + &
      & 2.0_wp * qp(5, :, :)**2 + qp(6, :, :)**2)

end subroutine multipole_norm


!> Add shell-resolved properties to dictionary in angular-momentum categories 
subroutine add_shell_feature(mol, bas, shell_prop, prefix, postfix, &
   & spin_label, dict, n_features)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Shell-resolved property data
   real(wp), intent(in) :: shell_prop(:)
   !> Feature name prefix
   character(len=*), intent(in) :: prefix
   !> Feature name postfix
   character(len=*), intent(in) :: postfix
   !> Spin label suffix
   character(len=*), intent(in) :: spin_label
   !> Dictionary for storing results
   type(double_dictionary_type), intent(inout) :: dict
   !> Number of features
   integer, intent(inout) :: n_features

   real(wp), allocatable :: angmom_prop(:, :)
   integer :: il

   ! Count number of features
   n_features = n_features + bas%maxl + 1

   ! Reorganize into angular-momentum categories
   allocate(angmom_prop(mol%nat, bas%maxl+1))
   call reorganize_shell_to_angmom(mol, bas, shell_prop, angmom_prop)

   do il = 0, bas%maxl
      call dict%add_entry( &
         trim(prefix)//"_"//trim(shell_labels(il))//trim(postfix)// &
            & trim(spin_label), angmom_prop(:, il+1))
   end do

end subroutine add_shell_feature


!> Add atom-resolved properties to dictionary
subroutine add_atom_feature(mol, atom_prop, prefix, postfix, spin_label, &
   & dict, n_features)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Atom-resolved property data
   real(wp), intent(in) :: atom_prop(:)
   !> Feature name prefix
   character(len=*), intent(in) :: prefix
   !> Feature name postfix
   character(len=*), intent(in) :: postfix
   !> Spin label suffix
   character(len=*), intent(in) :: spin_label
   !> Dictionary for storing results
   type(double_dictionary_type), intent(inout) :: dict
   !> Number of features
   integer, intent(inout) :: n_features

   ! Count number of features
   n_features = n_features + 1

   call dict%add_entry(trim(prefix)//trim(postfix)// &
      & trim(spin_label), atom_prop)

end subroutine add_atom_feature


!> Reorganize shell-resolved properties into angular-momentum categories
subroutine reorganize_shell_to_angmom(mol, bas, shell_prop, angmom_prop)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Shell-resolved property data
   real(wp), intent(in) :: shell_prop(:)
   !> Atomwise angular-momentum-resolved property data
   real(wp), intent(out) :: angmom_prop(:, :)

   integer :: iat, izp, is, ish, li

   angmom_prop(:, :) = 0.0_wp
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      do ish = 1, bas%nsh_at(iat)
         li = bas%cgto(ish, izp)%ang

         angmom_prop(iat, li+1) = shell_prop(is+ish)
      end do
   end do

end subroutine reorganize_shell_to_angmom


!> Convolution for scalar atom-resolved properties
!> with rescaling based on the source coordination number
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
            ! Convolution with rescaling based on source coordination number
            tmp = tmp + atom_prop(jat) / (kernel(iat, jat, ksc) &
               & * (cn(jat, ksc) + 1.0_wp))
         end do
         ext_prop(iat, ksc) = tmp
      end do
   end do

end subroutine convolve_scalar


!> Convolution of multipole moments induced by atom-resolved charges
!> with an optional additional multipole moment contribution
subroutine convolve_multipole(mol, qat, nscale, kernel, cn, &
   & ext_dp, ext_qp, dp, qp)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Number of convolution length scales
   integer, intent(in) :: nscale
   !> Convolution kernel
   real(wp), intent(in) :: kernel(:, :, :)
   !> Convolution coordination number
   real(wp), intent(in) :: cn(:, :)
   !> Convolved dipole moments
   real(wp), intent(out) :: ext_dp(:, :, :)
   !> Convolved quadrupole moments
   real(wp), intent(out) :: ext_qp(:, :, :)
   !> Optional atom-resolved dipole moments
   real(wp), intent(in), optional :: dp(:, :)
   !> Optional atom-resolved quadrupole moments
   real(wp), intent(in), optional :: qp(:, :)

   integer :: iat, jat, ksc
   real(wp) :: weight, trace
   real(wp) :: rij(3), dj(3)
   real(wp) :: tmp_dp(3), tmp_qp_trans(6), tmp_qp(6)

   ext_dp(:, :, :) = 0.0_wp
   ext_qp(:, :, :) = 0.0_wp

   !$omp parallel do default(none) collapse(2) schedule(runtime) &
   !$omp shared(mol, qat, nscale, kernel, cn, ext_dp, ext_qp, dp, qp) &
   !$omp private(iat, jat, ksc, weight, trace, rij, dj) &
   !$omp private(tmp_dp, tmp_qp_trans, tmp_qp)
   do ksc = 1, nscale
      do iat = 1, mol%nat

         tmp_dp(:) = 0.0_wp
         tmp_qp_trans(:) = 0.0_wp
         tmp_qp(:) = 0.0_wp

         do jat = 1, mol%nat

            ! Convolution weight with rescaling based on source coordination number
            weight = 1.0_wp / (kernel(iat, jat, ksc) * (cn(jat, ksc) + 1.0_wp))

            rij(:) = mol%xyz(:, iat) - mol%xyz(:, jat)

            ! Check for dipole moments
            if (present(dp)) then
               dj(:) = dp(:, jat)
            else
               dj(:) = 0.0_wp
            end if

            ! Dipole moment and translated dipole moment due to charges
            tmp_dp(:) = tmp_dp + weight * (dj(:) - rij(:) * qat(jat))

            ! Translated quadrupole moment due to charges and dipole moments
            ! Cartesian quadrupole ordering: xx, xy, yy, xz, yz, zz
            tmp_qp_trans(1) = tmp_qp_trans(1) + 1.5_wp * weight &
               & * (-2.0_wp * rij(1) * dj(1) + rij(1) * rij(1) * qat(jat))
            tmp_qp_trans(2) = tmp_qp_trans(2) + 1.5_wp * weight &
               & * (-(rij(1) * dj(2) + rij(2) * dj(1)) + rij(1) * rij(2) * qat(jat))
            tmp_qp_trans(3) = tmp_qp_trans(3) + 1.5_wp * weight &
               & * (-2.0_wp * rij(2) * dj(2) + rij(2) * rij(2) * qat(jat))
            tmp_qp_trans(4) = tmp_qp_trans(4) + 1.5_wp * weight &
               & * (-(rij(3) * dj(1) + rij(1) * dj(3)) + rij(1) * rij(3) * qat(jat))
            tmp_qp_trans(5) = tmp_qp_trans(5) + 1.5_wp * weight &
               & * (-(rij(3) * dj(2) + rij(2) * dj(3)) + rij(2) * rij(3) * qat(jat))
            tmp_qp_trans(6) = tmp_qp_trans(6) + 1.5_wp * weight &
               & * (-2.0_wp * rij(3) * dj(3) + rij(3) * rij(3) * qat(jat))

            if (present(qp)) then
               ! Quadrupole moment contribution
               tmp_qp(:) = tmp_qp(:) + qp(:, jat) * weight
            end if
         end do

         ! Remove trace from quadrupole moment _transuced by point charges and dipoles
         trace = (tmp_qp_trans(1) + tmp_qp_trans(3) + tmp_qp_trans(6)) / 3.0_wp
         tmp_qp_trans(1) = tmp_qp_trans(1) - trace
         tmp_qp_trans(3) = tmp_qp_trans(3) - trace
         tmp_qp_trans(6) = tmp_qp_trans(6) - trace

         ext_dp(:, iat, ksc) = tmp_dp(:)
         ext_qp(:, iat, ksc) = tmp_qp_trans(:) + tmp_qp(:)

      end do
   end do

end subroutine convolve_multipole

end module tblite_post_processing_xtbml_density
