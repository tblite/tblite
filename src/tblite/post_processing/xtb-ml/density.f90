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
!> Desnity based fetaures using Mulliken partitioning
module tblite_xtbml_density_based
   use mctc_env, only : wp
   use tblite_xtbml_feature_type, only : xtbml_feature_type
   use tblite_wavefunction_type, only : wavefunction_type
   use mctc_io, only : structure_type
   use tblite_integral_type, only : integral_type
   use tblite_container, only : container_cache
   use tblite_context , only : context_type
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_xtbml_convolution, only : xtbml_convolution_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_ncoord_exp, only : new_exp_ncoord, exp_ncoord_type
   use tblite_wavefunction_mulliken, only : get_mulliken_shell_multipoles
   use tblite_output_format, only : format_string
   implicit none
   private
   character(len=*), parameter :: label = "density-based features"
   real(wp), parameter :: one = 1.0_wp

   type, public, extends(xtbml_feature_type) :: xtbml_density_features_type

      real(wp),allocatable ::  mulliken_shell(:)
      real(wp),allocatable ::  dipm_shell(:)
      real(wp),allocatable ::  qm_shell(:)
      real(wp),allocatable ::  partial_charge_atom(:)
      real(wp),allocatable ::  delta_partial_charge(:, :)
      real(wp),allocatable ::  dipm_atom(:)
      real(wp),allocatable ::  delta_dipm(:, :)
      real(wp),allocatable ::  qm_atom(:)
      real(wp),allocatable ::  delta_qm(:, :)
      real(wp),allocatable ::  delta_dipm_e(:, :)
      real(wp),allocatable ::  delta_qm_e(:, :)
      real(wp),allocatable ::  delta_dipm_Z(:, :)
      real(wp),allocatable ::  delta_qm_Z(:, :)
      !> shell dipm xyz
      real(wp),allocatable ::  dipm_shell_xyz(:, :)
      !> dipm xyz
      real(wp),allocatable ::  dipm_atom_xyz(:, :)
      !> delta dipm xyz
      real(wp),allocatable ::  delta_dipm_xyz(:, :, :)
      !> shell qm xyz
      real(wp),allocatable ::  qm_shell_xyz(:, :)
      !> qm xyz
      real(wp),allocatable ::  qm_atom_xyz(:, :)
      !> delta qm xyz
      real(wp),allocatable ::  delta_qm_xyz(:, :, :)
      !> delta dipm only electron effect
      real(wp),allocatable ::  delta_dipm_e_xyz(:, :, :)
      !> delta qm only electron effect
      real(wp),allocatable ::  delta_qm_e_xyz(:, :, :)
      !> delta dipm only nuclear effect
      real(wp),allocatable ::  delta_dipm_Z_xyz(:, :, :)
      !> delta qm only nuclear effect
      real(wp),allocatable ::  delta_qm_Z_xyz(:, :, :)

      logical :: return_xyz

   contains
      procedure :: compute_features
      procedure :: compute_extended
      procedure, private :: allocate
      procedure, private :: allocate_extended
      procedure :: setup
   end type

contains

subroutine setup(self)
   class(xtbml_density_features_type) :: self
   self%label = label
   if (allocated(self%dict)) deallocate(self%dict)
   allocate(self%dict)
   if (allocated(self%dict_ext)) deallocate(self%dict_ext)
   allocate(self%dict_ext)
end subroutine

subroutine compute_features(self, mol, wfn, integrals, calc, cache_list)
   !> Instance of the feature container
   class(xtbml_density_features_type), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction strcuture data
   type(wavefunction_type), intent(in) :: wfn
   !> Integral container
   type(integral_type), intent(in) :: integrals
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Cache list for storing caches of various interactions
   type(container_cache), intent(inout) :: cache_list(:)
   character(len=20), allocatable :: tmp_labels(:)
   real(wp), allocatable :: tmp_s_array(:), tmp_p_array(:), tmp_d_array(:), tmp_array(:, :)
   real(wp) :: z(mol%nat), dipm_shell_tmp(3, calc%bas%nsh, 1), qm_shell_tmp(6, calc%bas%nsh, 1)
   integer :: i, mu

   ! allocate features
   call self%allocate(calc%bas%nsh, mol%nat)
   ! compute shellwise mulliken charges
   call mulliken_shellwise(calc%bas%ao2sh, wfn%density(:, :, 1), integrals%overlap, self%mulliken_shell)
   z = 0.0_wp
   ! get screened nuclear charge
   call mol_set_nuclear_charge(mol%num, mol%id, z)
   ! sum up shelll charges to atomic charges
   call sum_up_mulliken(calc%bas%sh2at, self%mulliken_shell, self%partial_charge_atom)
   ! get partial charge by subtracting the core charge
   self%partial_charge_atom = -self%partial_charge_atom+z

   !dipole moments shellwise and then atomwise
   call get_mulliken_shell_multipoles(calc%bas, integrals%dipole, wfn%density, dipm_shell_tmp)
   self%dipm_shell_xyz=dipm_shell_tmp(:, :, 1)
   !quadrupole moments shellwise and then atomwise
   call get_mulliken_shell_multipoles(calc%bas, integrals%quadrupole, wfn%density, &
   & qm_shell_tmp)
   self%qm_shell_xyz=qm_shell_tmp(:, :, 1)
   ! get atomic dipole moments from wavefunctin object
   self%dipm_atom_xyz=wfn%dpat(:, :, 1)
   ! get atomic quadrupole moments from wavefunctin object
   self%qm_atom_xyz=wfn%qpat(:, :, 1)

   call comp_norm(calc%bas%nsh, dipm_shell_tmp(:, :, 1), qm_shell_tmp(:, :, 1), self%dipm_shell, self%qm_shell)
   call comp_norm(mol%nat, wfn%dpat(:, :, 1), wfn%qpat(:, :, 1), self%dipm_atom, self%qm_atom)
   associate(dict => self%dict)
      call resolve_shellwise(self%mulliken_shell, tmp_s_array, tmp_p_array, tmp_d_array, calc%bas%nsh_at)
      call dict%add_entry("p_s", tmp_s_array)
      call dict%add_entry("p_p", tmp_p_array)
      call dict%add_entry("p_d", tmp_d_array)
      call resolve_shellwise(self%dipm_shell, tmp_s_array, tmp_p_array, tmp_d_array, calc%bas%nsh_at)
      call dict%add_entry("dipm_s", tmp_s_array)
      call dict%add_entry("dipm_p", tmp_p_array)
      call dict%add_entry("dipm_d", tmp_d_array)
      if (self%return_xyz) then
         call resolve_xyz_shell(self%dipm_shell_xyz, tmp_array, calc%bas%nsh_at)
         tmp_labels = [ character(len=20) :: &
         &"dipm_s_x", "dipm_s_y", "dipm_s_z",&
         &"dipm_p_x", "dipm_p_y", "dipm_p_z",&
         &"dipm_d_x", "dipm_d_y", "dipm_d_z"]

         do i = 1, size(tmp_labels)
            call dict%add_entry(trim(tmp_labels(i)), tmp_array(:, i))
         end do
      end if
      call resolve_shellwise(self%qm_shell, tmp_s_array, tmp_p_array, tmp_d_array, calc%bas%nsh_at)
      call dict%add_entry("qm_s", tmp_s_array)
      call dict%add_entry("qm_p", tmp_p_array)
      call dict%add_entry("qm_d", tmp_d_array)
      if (self%return_xyz) then
         call resolve_xyz_shell(self%qm_shell_xyz, tmp_array, calc%bas%nsh_at)
         tmp_labels = [ character(len=20) :: "qm_s_xx", &
         & "qm_s_xy", "qm_s_yy", "qm_s_xz", "qm_s_yz", "qm_s_zz",&
         &"qm_p_xx", "qm_p_xy", "qm_p_yy", "qm_p_xz", "qm_p_yz", "qm_p_zz",&
         &"qm_d_xx", "qm_d_xy", "qm_d_yy", "qm_d_xz", "qm_d_yz", "qm_d_zz"]

         do i = 1, size(tmp_labels)
            call dict%add_entry(trim(tmp_labels(i)), tmp_array(:, i))
         end do
      end if

      call dict%add_entry("q_A", self%partial_charge_atom)
      call dict%add_entry("dipm_A", self%dipm_atom)
      if (self%return_xyz) then
         tmp_labels = [ character(len=20) :: &
         &"dipm_A_x", "dipm_A_y", "dipm_A_z"]
         do i = 1, size(tmp_labels)
            call dict%add_entry(trim(tmp_labels(i)), self%dipm_atom_xyz(i, :))
         end do
      end if

      call dict%add_entry("qm_A", self%qm_atom)
      if (self%return_xyz) then
         tmp_labels = [ character(len=20) :: &
         &"qm_A_xx", "qm_A_xy", "qm_A_yy", "qm_A_xz", "qm_A_yz", "qm_A_zz"]

         do i = 1, size(tmp_labels)
            call dict%add_entry(trim(tmp_labels(i)), self%qm_atom_xyz(i, :))
         end do
      end if
   end associate
end subroutine

!> get the screened nuclear charge for each atom
! which means the nuclear charge minus the core charge
subroutine mol_set_nuclear_charge(at, id, z)
   !> atomic number, as function of the index
   integer, intent(in) :: at(:)
   !> list to retrieve index for each atom
   integer, intent(in) :: id(:)
   !> screened nuclear charge
   real(wp), intent(out) :: z(:)
   integer :: i
   do i = 1, size(z)
      z(i) = real(at(id(i)), wp) - real(ncore(at(id(i))), wp)
      ! if the atom is a lanthanide
      if (at(id(i)) > 57 .and. at(id(i)) < 72) z(i) = 3.0_wp
   end do
contains
   !> get core nuclear charge
   pure elemental integer function ncore(at)
      integer, intent(in) :: at
      if (at .le. 2) then
         ncore = 0
      elseif (at .le. 10) then
         ncore = 2
      elseif (at .le. 18) then
         ncore = 10
      elseif (at .le. 29) then   !zn
         ncore = 18
      elseif (at .le. 36) then
         ncore = 28
      elseif (at .le. 47) then
         ncore = 36
      elseif (at .le. 54) then
         ncore = 46
      elseif (at .le. 71) then
         ncore = 54
      elseif (at .le. 79) then
         ncore = 68
      elseif (at .le. 86) then
         ncore = 78
      end if
   end function ncore
end subroutine mol_set_nuclear_charge

!> seperate a vector of properties for all shells into sperate arrays for each shell type
subroutine resolve_shellwise(shell_prop, array_s, array_p, array_d, at2nsh)
   !> array containing the shell properties
   real(wp), intent(in) :: shell_prop(:)
   !> array for s shell properties
   real(wp), allocatable, intent(out) :: array_s(:)
   !> array for p shell properties
   real(wp), allocatable, intent(out) :: array_p(:)
   !> array for d shell properties
   real(wp), allocatable, intent(out) :: array_d(:)
   !> atom 2 number of shell mapping
   integer, intent(in) :: at2nsh(:)
   integer :: nsh, i, nat
   nsh = 1
   nat = size(at2nsh)
   if (allocated(array_s)) deallocate(array_s)
   if (allocated(array_p)) deallocate(array_p)
   if (allocated(array_d)) deallocate(array_d)
   allocate(array_s(nat), array_p(nat), array_d(nat), source= 0.0_wp)
   !in this loop we are filling shellwise properties 
   do i = 1, nat
      array_s(i) = shell_prop(nsh) !s shell always filled
      nsh = nsh + 1
      ! if s and p are filled
      if (at2nsh(i) == 2) then
         array_p(i) = shell_prop(nsh)
         nsh = nsh + 1
      ! if s, p and d are filled
      elseif (at2nsh(i) == 3) then
         array_p(i) = shell_prop(nsh)
         nsh = nsh + 1
         array_d(i) = shell_prop(nsh)
         nsh = nsh + 1
      end if
   end do
end subroutine

!> seperate a vector of properties for all shells into sperate arrays for each shell type, for tensorial properties
subroutine resolve_xyz_shell(mult_xyz, array, at2nsh)
   !> array containing the tensorial shell properties
   real(wp), intent(in) :: mult_xyz(:, :)
   !> array for xyz shell properties
   real(wp), allocatable :: array(:, :)
   !> atom 2 number of shell mapping
   integer, intent(in) :: at2nsh(:)
   integer :: j, k, nsh, id_tmp, nat, i
   nsh = 1
   nat = size(at2nsh)
   ! we reuse array so if we have dipm we need 3 columns, if we have qm we need 6 columns
   if (allocated(array)) deallocate(array)
   allocate(array(nat, 3*size(mult_xyz, dim=1)), source = 0.0_wp)

   do k = 1, nat
      ! id_tmp is used to keep track of the current index in the array
      id_tmp = 1
      j = 1
      do i = id_tmp, (id_tmp + size(mult_xyz, dim=1) - 1)
         array(k, i) = mult_xyz(j, nsh)
         j = j + 1
      end do
      nsh = nsh + 1
      if (at2nsh(k) == 2) then
         id_tmp = id_tmp + size(mult_xyz, dim=1)
         j = 1
         do i = id_tmp, id_tmp + size(mult_xyz, dim=1) - 1
            array(k, i) = mult_xyz(j, nsh)
            j = j + 1
         end do
         nsh = nsh + 1
      elseif (at2nsh(k) == 3) then
         id_tmp = id_tmp + size(mult_xyz, dim=1)
         j = 1
         do i = id_tmp, id_tmp + size(mult_xyz, dim=1) - 1
            array(k, i) = mult_xyz(j, nsh)
            j = j + 1
         end do
         nsh = nsh + 1
         id_tmp = id_tmp + size(mult_xyz, dim=1)
         j = 1
         do i = id_tmp, id_tmp + size(mult_xyz, dim=1) - 1
            array(k, i) = mult_xyz(j, nsh)
            j = j + 1
         end do
         nsh = nsh + 1
      end if
   end do
end subroutine

subroutine compute_extended(self, mol, wfn, integrals, calc, cache_list, convolution) 
   !> Instance of feature container
   class(xtbml_density_features_type), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction strcuture data
   type(wavefunction_type), intent(in) :: wfn
   !> Integral container
   type(integral_type), intent(in) :: integrals
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Cache list for storing caches of various interactions
   type(container_cache), intent(inout) :: cache_list(:)
   !> Convolution container
   type(xtbml_convolution_type), intent(in) :: convolution
   real(wp), allocatable :: mull_charge_atomic(:)
   real(wp), allocatable :: z(:)
   integer :: i, j, n 
   character(len=20), allocatable :: tmp_labels(:)
   character(len=:), allocatable :: tmp_label, a_label
   n = convolution%n_a
   allocate(z(mol%nat), source = 0.0_wp)
   allocate(mull_charge_atomic(mol%nat), source = 0.0_wp)

   ! allocate extended features
   call self%allocate_extended(calc%bas%nsh, mol%nat, n)

   ! compute extended partial charges
   call get_delta_partial_charge(self%partial_charge_atom, mol%xyz, &
      convolution, self%delta_partial_charge)
   ! compute extended multipole moments
   call get_delta_mm(self%partial_charge_atom, wfn%dpat(:, :, 1), wfn%qpat(:, :, 1), mol%id, mol%xyz, &
      convolution, self%delta_dipm_xyz, self%delta_qm_xyz)
   ! copmute the norm of the extended multipole moments
   call comp_norm_3(self%delta_dipm_xyz, self%delta_qm_xyz, self%delta_dipm, self%delta_qm)
   ! get screened nuclear charge
   call mol_set_nuclear_charge(mol%num, mol%id, z)
   ! compute the extended multipole moments only considering the nuclear charge
   call get_delta_mm_Z(z, wfn%dpat(:, :, 1), wfn%qpat(:, :, 1), mol%id, mol%xyz, convolution, self%delta_dipm_Z_xyz, &
      self%delta_qm_Z_xyz)
   ! sum up shelll charges to atomic charges
   call sum_up_mulliken(calc%bas%sh2at, self%mulliken_shell, mull_charge_atomic)
   ! compute the extended multipole moments only considering the electronic effects
   call get_delta_mm_p(mull_charge_atomic, wfn%dpat(:, :, 1), wfn%qpat(:, :, 1), mol%id, mol%xyz, convolution, &
      self%delta_dipm_e_xyz, self%delta_qm_e_xyz)
   !compute norms
   call comp_norm_3(self%delta_dipm_e_xyz, self%delta_qm_e_xyz, self%delta_dipm_e, self%delta_qm_e)
   call comp_norm_3(self%delta_dipm_Z_xyz, self%delta_qm_Z_xyz, self%delta_dipm_Z, self%delta_qm_Z)
   associate( dict => self%dict_ext)
      do j = 1, n
         a_label = "_"//trim(adjustl(format_string(convolution%a(j), '(f12.2)')))
         if (a_label .eq. "_1.00") a_label = ""
         tmp_label = "delta_q_A"//a_label
         call dict%add_entry(tmp_label, self%delta_partial_charge(:, j))

         tmp_label = "delta_dipm_A"//a_label
         call dict%add_entry(tmp_label, self%delta_dipm(:, j))
         if (self%return_xyz) then
            tmp_labels = [ character(len=20) :: &
            &"delta_dipm_A_x", "delta_dipm_A_y", "delta_dipm_A_z"]
            do i = 1, size(tmp_labels)
               tmp_label = trim(tmp_labels(i))//a_label
               call dict%add_entry(tmp_label, self%delta_dipm_xyz(i, :, j))
            end do
         end if

         tmp_label = "delta_qm_A"//a_label
         call dict%add_entry(tmp_label, self%delta_qm(:, j))
         if (self%return_xyz) then
            tmp_labels = [ character(len=20) :: &
            &"delta_qm_A_xx", "delta_qm_A_xy", "delta_qm_A_yy", "delta_qm_A_xz", "delta_qm_A_yz", "delta_qm_A_zz"]

            do i = 1, size(tmp_labels)
               tmp_label = trim(tmp_labels(i))//a_label
               call dict%add_entry(tmp_label, self%delta_qm_xyz(i, :, j))
            end do
         end if

         tmp_label = "delta_dipm_e"//a_label
         call dict%add_entry(tmp_label, self%delta_dipm_e(:, j))
         if (self%return_xyz) then
            tmp_labels = [ character(len=20) :: &
            &"delta_dipm_e_x", "delta_dipm_e_y", "delta_dipm_e_z"]
            do i = 1, size(tmp_labels)
               tmp_label = trim(tmp_labels(i))//a_label
               call dict%add_entry(tmp_label, self%delta_dipm_e_xyz(i, :, j))
            end do
         end if

         tmp_label = "delta_qm_e"//a_label
         call dict%add_entry(tmp_label, self%delta_qm_e(:, j))
         if (self%return_xyz) then
            tmp_labels = [ character(len=20) :: &
            &"delta_qm_e_xx", "delta_qm_e_xy", "delta_qm_e_yy", "delta_qm_e_xz", "delta_qm_e_yz", "delta_qm_e_zz"]

            do i = 1, size(tmp_labels)
               tmp_label = trim(tmp_labels(i))//a_label
               call dict%add_entry(tmp_label, self%delta_qm_e_xyz(i, :, j))
            end do
         end if
         tmp_label = "delta_dipm_Z"//a_label
         call dict%add_entry(tmp_label, self%delta_dipm_Z(:, j))
         if (self%return_xyz) then
            tmp_labels = [ character(len=20) :: &
            &"delta_dipm_Z_x", "delta_dipm_Z_y", "delta_dipm_Z_z"]
            do i = 1, size(tmp_labels)
               tmp_label = trim(tmp_labels(i))//a_label
               call dict%add_entry(tmp_label, self%delta_dipm_Z_xyz(i, :, j))
            end do
         end if
         tmp_label = "delta_qm_Z"//a_label
         call dict%add_entry(tmp_label, self%delta_qm_Z(:, j))
         if (self%return_xyz) then
            tmp_labels = [ character(len=20) :: &
            &"delta_qm_Z_xx", "delta_qm_Z_xy", "delta_qm_Z_yy", "delta_qm_Z_xz", "delta_qm_Z_yz", "delta_qm_Z_zz"]

            do i = 1, size(tmp_labels)
               tmp_label = trim(tmp_labels(i))//a_label
               call dict%add_entry(tmp_label, self%delta_qm_Z_xyz(i, :, j))
            end do
         end if
      end do
   end associate
end subroutine

!> Allocate features
subroutine allocate(self, nsh, nat)
   !> Instance of feature container
   class(xtbml_density_features_type), intent(inout) :: self
   !> Number of shells
   integer, intent(in) :: nsh
   !> Number of atoms
   integer, intent(in) :: nat

   allocate(self%mulliken_shell(nsh), source=0.0_wp)
   allocate(self%dipm_shell(nsh), source=0.0_wp)
   allocate(self%qm_shell(nsh), source=0.0_wp)
   allocate(self%partial_charge_atom(nat), source=0.0_wp)
   allocate(self%dipm_atom(nat), source=0.0_wp)
   allocate(self%qm_atom(nat), source=0.0_wp)
   !
   allocate(self%dipm_shell_xyz(3, nsh), source=0.0_wp)
   allocate(self%dipm_atom_xyz(3, nat), source=0.0_wp)
   allocate(self%qm_shell_xyz(6, nsh), source=0.0_wp)
   allocate(self%qm_atom_xyz(6, nat), source=0.0_wp)
end subroutine
   
!> Allocate extended features
subroutine allocate_extended(self, nsh, nat, n_a)
   !> Instance of feature container
   class(xtbml_density_features_type), intent(inout) :: self
   !> Number of shells
   integer, intent(in) :: nsh
   !> Number of atoms
   integer, intent(in) :: nat
   !> Number of convolution values
   integer, intent(in) :: n_a

   allocate(self%delta_partial_charge(nat, n_a), source=0.0_wp)
   allocate(self%delta_dipm(nat, n_a), source=0.0_wp)
   allocate(self%delta_qm(nat, n_a), source=0.0_wp)
   allocate(self%delta_dipm_e(nat, n_a), source=0.0_wp)
   allocate(self%delta_qm_e(nat, n_a), source=0.0_wp)
   allocate(self%delta_dipm_Z(nat, n_a), source=0.0_wp)
   allocate(self%delta_qm_Z(nat, n_a), source=0.0_wp)
   !
   allocate(self%delta_dipm_xyz(3, nat, n_a), source=0.0_wp)
   allocate(self%delta_qm_xyz(6, nat, n_a), source=0.0_wp)
   !
   allocate(self%delta_dipm_e_xyz(3, nat, n_a), source=0.0_wp)
   allocate(self%delta_qm_e_xyz(6, nat, n_a), source=0.0_wp)
   allocate(self%delta_dipm_Z_xyz(3, nat, n_a), source=0.0_wp)
   allocate(self%delta_qm_Z_xyz(6, nat, n_a), source=0.0_wp)

end subroutine

!> Compute shellwise mulliken charges
subroutine mulliken_shellwise(ao2shell, p, s, charges_shell)
   !> AO to shell mapping
   integer, intent(in) :: ao2shell(:)
   !> Overlap matrix
   real(wp), intent(in) :: s(:, :)
   !> Density matrix
   real(wp), intent(in) :: p(:, :)
   !> Shellwise mulliken charges
   real(wp), intent(out) :: charges_shell(:)
   integer ::  mu, nu, nao
   nao = size(p, 1)

   charges_shell = 0.0_wp
   !$omp parallel do default(none) collapse(2)&
   !$omp private(mu, nu) shared(charges_shell, nao, ao2shell, p, s)
   do mu = 1, nao
      do nu = 1, nao
         !$omp atomic
         charges_shell(ao2shell(mu)) = charges_shell(ao2shell(mu)) + p(mu, nu)*s(nu, mu)
      end do
   end do

end subroutine

!> Copmute the extended atomic partial charges
subroutine get_delta_partial_charge(atom_partial, xyz, conv, delta_partial)
   !> Convolution container
   type(xtbml_convolution_type), intent(in) :: conv
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Atomic partial charges
   real(wp), intent(in) :: atom_partial(:)
   !> Delta partial charges
   real(wp), intent(out) :: delta_partial(:, :)
   integer :: a, b, k
   integer :: nat, n_a

   nat = size(atom_partial)
   n_a = size(conv%a)
   delta_partial = 0.0_wp

   !$omp parallel do default(none) collapse(2)&
   !$omp private(a, b, k) shared(delta_partial, nat, n_a, conv, atom_partial)
   do k = 1, n_a
      do a = 1, nat
         do b = 1, nat
            !$omp atomic
            delta_partial(a, k) = delta_partial(a, k) + atom_partial(b) / (conv%kernel(a, b, k) * (conv%cn(b, k)+1))
         enddo
      enddo
   enddo

end subroutine

!summing up the mulliken charges to get atomic charges
subroutine sum_up_mulliken(ash, mull_shell, mull_at)
   !> shell to atom mapping
   integer, intent(in) :: ash(:)
   !> shellwise mulliken charges
   real(wp), intent(in) :: mull_shell(:)
   !> atomic charges
   real(wp), intent(out) :: mull_at(:)
   integer :: i

   mull_at = 0.0_wp

   do i = 1, size(ash)
      mull_at(ash(i)) = mull_at(ash(i)) + mull_shell(i)
   end do

end subroutine

!> Compute the norm for extended dipole and quadrupole moments
subroutine get_delta_mm(q, dipm, qp, at, xyz, conv, delta_dipm, delta_qp)
   !> Convolution container
   type(xtbml_convolution_type), intent(in) :: conv
   !> atom list
   integer, intent(in) :: at(:)
   !> dipole moment
   real(wp), intent(in) :: dipm(:, :)
   !> quadrupole moment
   real(wp), intent(in) :: qp(:, :)
   !> xyz coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> atomic charges
   real(wp), intent(in) :: q(:)
   !> delta dipole moment
   real(wp), intent(out) :: delta_dipm(:, :, :)
   !> delta quadrupole moment
   real(wp), intent(out) :: delta_qp(:, :, :)
   integer :: a, b, k, n_a, nat
   real(wp) :: result 
   real(wp), allocatable :: r_ab(:), qp_part(:, :, :)
   n_a = size(conv%a)
   nat = size(at)

   allocate(r_ab(3), source = 0.0_wp)
   allocate(qp_part(6, nat, n_a), source = 0.0_wp)

   delta_dipm = 0.0_wp
   delta_qp = 0.0_wp
   qp_part = 0.0_wp


   !$omp parallel do default(none) collapse(2)&
   !$omp private(a, b, k, result, r_ab) shared(delta_dipm, delta_qp, nat, n_a, conv, xyz, q, qp_part, dipm, qp)
   do k = 1, n_a
      do a = 1, nat
         do b = 1, nat

            result = conv%kernel(a, b, k)
            r_ab = xyz(:, a) - xyz(:, b)

            delta_dipm(:, a, k) = delta_dipm(:, a, k) + (dipm(:, b) - r_ab(:) * q(b)) / (result * (conv%cn(b, k) + one))
            !sorting of qp xx,xy,yy,xz,yz,zz
            !$omp atomic
            delta_qp(1, a, k) = delta_qp(1, a, k) + &
               (1.5_wp * (-one * (r_ab(1)*dipm(1, b) + r_ab(1) * dipm(1, b)) + r_ab(1)*r_ab(1) * q(b))) / &
               (result * (conv%cn(b, k) + one))
            !$omp atomic
            delta_qp(2, a, k) = delta_qp(2, a, k) + &
               (1.5_wp * (-one * (r_ab(1)*dipm(2, b) + r_ab(2) * dipm(1, b)) + r_ab(1) * r_ab(2) * q(b))) / &
               (result * (conv%cn(b, k) + one))
            !$omp atomic
            delta_qp(3, a, k) = delta_qp(3, a, k) + &
               (1.5_wp * (-one * (r_ab(2) * dipm(2, b) + r_ab(2) * dipm(2, b)) + r_ab(2) * r_ab(2) * q(b))) / &
               (result*(conv%cn(b, k) + one))
            !$omp atomic
            delta_qp(4, a, k) = delta_qp(4, a, k) + &
               (1.5_wp * (-one * (r_ab(3) * dipm(1, b) + r_ab(1) * dipm(3, b)) + r_ab(1) * r_ab(3) * q(b))) / &
               (result * (conv%cn(b, k) + one))
            !$omp atomic
            delta_qp(5, a, k) = delta_qp(5, a, k) + &
               (1.5_wp * (-one * (r_ab(3) * dipm(2, b) + r_ab(2) * dipm(3, b)) + r_ab(2) * r_ab(3) * q(b))) / &
               (result * (conv%cn(b, k) + one))
            !$omp atomic
            delta_qp(6, a, k) = delta_qp(6, a, k) + &
               (1.5_wp * (-one * (r_ab(3) * dipm(3, b) + r_ab(3) * dipm(3, b)) + r_ab(3) * r_ab(3) * q(b))) / &
               (result * (conv%cn(b, k) + one))

            qp_part(:, a, k) = qp_part(:, a, k) + &
               qp(:, b) / (result * (conv%cn(b, k) + one))
         end do
      end do
   end do

   call remove_trac_qp(delta_qp, qp_part)

end subroutine

!> Compute the norm for extended dipole and quadrupole moments
subroutine comp_norm_3(dipm, qm, dipm_norm, qm_norm)
   !> Dipole moment
   real(wp), intent(in) :: dipm(:, :, :)
   !> Quadrupole moment
   real(wp), intent(in) :: qm(:, :, :)
   !> Dipole moment norm
   real(wp), intent(out) :: dipm_norm(:, :)
   !> Quadrupole moment norm
   real(wp), intent(out) :: qm_norm(:, :)
   real(wp), allocatable :: r(:, :), r2(:, :)
   integer :: i, ndim
   ndim = size(dipm, 2)

   allocate(r(ndim, size(dipm, 3)), source = 0.0_wp)
   allocate(r2(ndim, size(dipm, 3)), source = 0.0_wp)

   do i = 1, ndim
      r2(i, :) = dipm(1, i, :)**2 + dipm(2, i, :)**2 + dipm(3, i, :)**2
   end do
   r = sqrt(r2)

   dipm_norm(:, :) = r(:, :)

   do i = 1, ndim
      r2(i, :) = qm(1, i, :)**2 + 2.0_wp * qm(2, i, :)**2 + qm(3, i, :)**2 + &
                  2.0_wp * qm(4, i, :)**2 + 2.0_wp * qm(5, i, :)**2 + qm(6, i, :)**2
   end do
   r = sqrt(r2)
   qm_norm(:, :) = r(:, :)

end subroutine

!> Compute the norm for dipole and quadrupole moments
subroutine comp_norm(ndim, dipm, qm, dipm_norm, qm_norm)
   !> Dipole moment
   real(wp), intent(in) :: dipm(:, :)
   !> Quadrupole moment
   real(wp), intent(in) :: qm(:, :)
   !> Dipole moment norm
   real(wp), intent(out) :: dipm_norm(:)
   !> Quadrupole moment norm
   real(wp), intent(out) :: qm_norm(:)

   real(wp), allocatable :: r(:), r2(:)
   integer :: i, ndim
   ndim = size(dipm, 2)

   allocate(r(ndim), source = 0.0_wp)
   allocate(r2(ndim), source = 0.0_wp)
   do i = 1, ndim
      r2(i) = dipm(1, i)**2 + dipm(2, i)**2 + dipm(3, i)**2
   end do
   r = sqrt(r2)

   dipm_norm(:) = r(:)

   do i = 1, ndim
      r2(i) = qm(1, i)**2 + 2.0_wp * qm(2, i)**2 + qm(3, i)**2 + &
               2.0_wp * qm(4, i)**2 + 2.0_wp * qm(5, i)**2 + qm(6, i)**2
   end do
   r = sqrt(r2)
   qm_norm(:) = r(:)

end subroutine

!> Compute extended quadrupole moment only due to nuclear effects
! all effects due to the electrons are set to 0, only the distribution of positive charges is left
subroutine get_delta_mm_Z(q, dipm, qp, at, xyz, conv, delta_dipm, delta_qp)
   !> Convolution container
   type(xtbml_convolution_type), intent(in) :: conv
   !> atom list
   integer, intent(in) :: at(:)
   !> dipole moment
   real(wp), intent(in) :: dipm(:, :)
   !> quadrupole moment
   real(wp), intent(in) :: qp(:, :)
   !> xyz coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> atomic charges
   real(wp), intent(in) :: q(:)
   !> delta dipole moment
   real(wp), intent(out) :: delta_dipm(:, :, :)
   !> delta quadrupole moment
   real(wp), intent(out) :: delta_qp(:, :, :)
   integer :: a, b, k, n_a, nat
   real(wp) :: result 
   real(wp), allocatable :: r_ab(:), qp_part(:, :, :)
   n_a = size(conv%a)
   nat = size(at)

   allocate(r_ab(3), source = 0.0_wp)
   allocate(qp_part(6, nat, n_a), source = 0.0_wp)

   delta_dipm = 0.0_wp
   delta_qp = 0.0_wp
   qp_part = 0.0_wp

   !$omp parallel do default(none) collapse(2)&
   !$omp private(a, b, k, result, r_ab) shared(delta_dipm, delta_qp, nat, n_a, conv, xyz, q)
   do k = 1, n_a
      do a = 1, nat
         do b = 1, nat
            result = conv%kernel(a, b, k)
            r_ab = xyz(:, a) - xyz(:, b)
            delta_dipm(:, a, k) = delta_dipm(:, a, k) + (-r_ab(:)*q(b))/(result*(conv%cn(b, k) + one))
            !sorting of qp xx,xy,yy,xz,yz,zz
            !$omp atomic
            delta_qp(1, a, k) = delta_qp(1, a, k) + &
               (1.5_wp * (r_ab(1) * r_ab(1) * q(b))) / (result * (conv%cn(b, k) + one))
            !$omp atomic
            delta_qp(2, a, k) = delta_qp(2, a, k) + &
               (1.5_wp * (r_ab(1) * r_ab(2) * q(b))) / (result * (conv%cn(b, k) + one))
            !$omp atomic
            delta_qp(3, a, k) = delta_qp(3, a, k) + &
               (1.5_wp * (r_ab(2) * r_ab(2) * q(b))) / (result * (conv%cn(b, k) + one))
            !$omp atomic
            delta_qp(4, a, k) = delta_qp(4, a, k) + &
               (1.5_wp * (r_ab(1) * r_ab(3) * q(b))) / (result * (conv%cn(b, k) + one))
            !$omp atomic
            delta_qp(5, a, k) = delta_qp(5, a, k) + &
               (1.5_wp * (r_ab(2) * r_ab(3) * q(b))) / (result * (conv%cn(b, k) + one))
            !$omp atomic
            delta_qp(6, a, k) = delta_qp(6, a, k) + &
               (1.5_wp * (r_ab(3) * r_ab(3) * q(b))) / (result * (conv%cn(b, k) + one))
         end do
      end do
   end do

   call remove_trac_qp(delta_qp, qp_part)

end subroutine

!> Compute extended quadrupole moment only due to electronic effects
subroutine get_delta_mm_p(q, dipm, qp, at, xyz, conv, delta_dipm, delta_qp) ! the sign of q was changed to respect the charge of the electrons
   !> Convolution container
   type(xtbml_convolution_type), intent(in) :: conv
   !> atom list
   integer, intent(in) :: at(:)
   !> dipole moment
   real(wp), intent(in) :: dipm(:, :)
   !> quadrupole moment
   real(wp), intent(in) :: qp(:, :)
   !> xyz coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> atomic charges
   real(wp), intent(in) :: q(:)
   !> delta dipole moment
   real(wp), intent(out) :: delta_dipm(:, :, :)
   !> delta quadrupole moment
   real(wp), intent(out) :: delta_qp(:, :, :)
   integer :: a, b, k, n_a, nat
   real(wp) :: result 
   real(wp), allocatable :: r_ab(:), qp_part(:, :, :)
   n_a = size(conv%a)
   nat = size(at)

   allocate(r_ab(3), source = 0.0_wp)
   allocate(qp_part(6, nat, n_a), source = 0.0_wp)

   delta_dipm = 0.0_wp
   delta_qp = 0.0_wp
   qp_part = 0.0_wp
   !$omp parallel do default(none) collapse(2)&
   !$omp private(a, b, k, result, r_ab) shared(delta_dipm, delta_qp, nat, n_a, conv, xyz, q, dipm, qp, qp_part)
   do k = 1, n_a
      do a = 1, nat
         do b = 1, nat
            result = conv%kernel(a, b, k)
            r_ab = xyz(:, a) - xyz(:, b)
            delta_dipm(:, a, k) = delta_dipm(:, a, k) + &
               (dipm(:, b) + r_ab(:) * q(b)) / (result * (conv%cn(b, k) + one))
            !sorting of qp xx,xy,yy,xz,yz,zz
            delta_qp(1, a, k) = delta_qp(1, a, k) + &
               (1.5_wp * (-one * (r_ab(1)*dipm(1, b) + r_ab(1) * dipm(1, b)) - r_ab(1) * r_ab(1) * q(b))) / &
               (result * (conv%cn(b, k) + one))
            !$omp atomic
            delta_qp(2, a, k) = delta_qp(2, a, k) + &
               (1.5_wp * (-one * (r_ab(1) * dipm(2, b) + r_ab(2)*dipm(1, b)) - r_ab(1) * r_ab(2) * q(b))) / &
               (result * (conv%cn(b, k) + one))
            !$omp atomic
            delta_qp(3, a, k) = delta_qp(3, a, k) + &
               (1.5_wp * (-one * (r_ab(2) * dipm(2, b) + r_ab(2) * dipm(2, b)) - r_ab(2) * r_ab(2) * q(b))) / &
               (result * (conv%cn(b, k) + one))
            !$omp atomic
            delta_qp(4, a, k) = delta_qp(4, a, k) + &
               (1.5_wp * (-one * (r_ab(3) * dipm(1, b) + r_ab(1) * dipm(3, b)) - r_ab(1) * r_ab(3) * q(b))) / &
               (result * (conv%cn(b, k) + one))
            !$omp atomic
            delta_qp(5, a, k) = delta_qp(5, a, k) + &
               (1.5_wp * (-one * (r_ab(3) * dipm(2, b) + r_ab(2) * dipm(3, b)) - r_ab(2) * r_ab(3) * q(b))) / &
               (result * (conv%cn(b, k) + one))
            !$omp atomic
            delta_qp(6, a, k) = delta_qp(6, a, k) + &
               (1.5_wp * (-one * (r_ab(3) * dipm(3, b) + r_ab(3) * dipm(3, b)) - r_ab(3) * r_ab(3) * q(b))) / &
               (result*(conv%cn(b, k) + one))
            qp_part(:, a, k) = qp_part(:, a, k) + qp(:, b) / (result * (conv%cn(b, k) + one))
         end do
      end do
   end do

   call remove_trac_qp(delta_qp, qp_part)

end subroutine

!> Remove the trace of the quadrupole moment
subroutine remove_trac_qp(qp_matrix, qp_part)
   !> Quadrupole moment on site
   real(wp), intent(in) :: qp_part(:, :, :)
   !> Quadrupole moment matrix
   real(wp), intent(inout) :: qp_matrix(:, :, :)
   integer :: i
   real(wp), allocatable :: tii(:)

   allocate(tii(size(qp_matrix, 3)))

   do i = 1, size(qp_matrix, 2)
      tii = qp_matrix(1, i, :) + qp_matrix(3, i, :) + qp_matrix(6, i, :)
      tii = tii / 3.0_wp
      qp_matrix(1, i, :) = qp_matrix(1, i, :) - tii
      qp_matrix(3, i, :) = qp_matrix(3, i, :) - tii
      qp_matrix(6, i, :) = qp_matrix(6, i, :) - tii
      qp_matrix(:, i, :) = qp_matrix(:, i, :) + qp_part(:, i, :)
   end do
end subroutine

end module tblite_xtbml_density_based
