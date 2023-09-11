module tblite_xtbml_density_based
   use mctc_env, only : wp
   use tblite_xtbml_feature_type, only : xtbml_feature_type
   use mctc_env, only : wp
   use tblite_wavefunction, only : wavefunction_type
   use mctc_io, only : structure_type
   use tblite_integral_type, only : integral_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_container, only : container_cache
   use tblite_context , only : context_type
   use tblite_double_dictionary, only : double_dictionary_type
   implicit none
   private

   integer, parameter :: features = 12
   integer, parameter :: ext_features = 40



   type, public, extends(xtbml_feature_type) :: xtbml_geometry_features_type
      character(len=*) :: label = "density-based features"
      type(double_dictionary_type) :: dict, dict_ext
      real(wp),allocatable ::  mulliken_shell(:)
      real(wp),allocatable ::  dipm_shell(:)
      real(wp),allocatable ::  qm_shell(:)
      real(wp),allocatable ::  partial_charge_atom(:)
      real(wp),allocatable ::  delta_partial_charge(:,:)
      real(wp),allocatable ::  dipm_atom(:)
      real(wp),allocatable ::  delta_dipm(:,:)
      real(wp),allocatable ::  qm_atom(:)
      real(wp),allocatable ::  delta_qm(:,:)
      real(wp),allocatable ::  delta_dipm_e(:,:)
      real(wp),allocatable ::  delta_qm_e(:,:)
      real(wp),allocatable ::  delta_dipm_Z(:,:)
      real(wp),allocatable ::  delta_qm_Z(:,:)
      !
      !> shell dipm xyz
      real(wp),allocatable ::  dipm_shell_xyz(:,:)
      !> dipm xyz
      real(wp),allocatable ::  dipm_atom_xyz(:,:)
      !> delta dipm xyz
      real(wp),allocatable ::  delta_dipm_xyz(:,:,:)
      !> shell qm xyz
      real(wp),allocatable ::  qm_shell_xyz(:,:)
      !> qm xyz
      real(wp),allocatable ::  qm_atom_xyz(:,:)
      !> delta qm xyz
      real(wp),allocatable ::  delta_qm_xyz(:,:,:)
      !
      !> delta dipm only electron effect
      real(wp),allocatable ::  delta_dipm_e_xyz(:,:,:)
      !> delta qm only electron effect
      real(wp),allocatable ::  delta_qm_e_xyz(:,:,:)
      !> delta dipm only nuclear effect
      real(wp),allocatable ::  delta_dipm_Z_xyz(:,:,:)
      !> delta qm only nuclear effect
      real(wp),allocatable ::  delta_qm_Z_xyz(:,:,:)

      integer :: n_features
      type(enum_density_features) :: labels = enum_geometry_features()
      type(enum_ext_density_features) :: ext_labels = enum_geometry_features()
      logical :: return_xyz

   contains
      procedure :: compute_features
      procedure :: compute_extended
      procedure :: get_n_features
      procedure :: get_feature_labels
      procedure, private :: allocate
      procedure, private :: allocate_extended
   end type

contains

subroutine compute_features(self, mol, wfn, integrals, calc, cache, prlevel, ctx)
   use tblite_ncoord_exp, only : new_exp_ncoord, exp_ncoord_type
   class(xtbml_feature_type), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction strcuture data
   type(wavefunction_type), intent(in) :: wfn
   !> Integral container
   type(integral_type) :: integrals
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Container
   type(container_cache), intent(inout) :: cache
   !> Context type
   type(context_type),intent(inout) :: ctx
   !> Print Level
   integer, intent(in) :: prlevel
   character(len=20), allocatable :: tmp_labels(:)
   real(wp), allocatable :: tmp_s_array(:), tmp_p_array(:), tmp_d_array(:), tmp_array(:, :)

    !shellwise mulliken charges
   call mulliken_shellwise(calc%bas%nao, calc%bas%nsh, calc%bas%ao2sh, wfn%density(:, :, wfn%nspin), &
      integrals%overlap, self%mulliken_shell)
   

   call mol_set_nuclear_charge(mol%nat, mol%num, mol%id, z)

   do mu=1, calc%bas%nsh
      self%partial_charge_atom(calc%bas%sh2at(mu))=self%partial_charge_atom(calc%bas%sh2at(mu))+self%mulliken_shell(mu)
   enddo
   self%partial_charge_atom=-self%partial_charge_atom+z
   
   !multipole moments shellwise und then atomwise
   call get_mulliken_shell_multipoles(calc%bas, integrals%dipole, wfn%density, &
      & dipm_shell_tmp)
   self%dipm_shell_xyz=dipm_shell_tmp(:, :, 1)
   
  
   call get_mulliken_shell_multipoles(calc%bas, integrals%quadrupole, wfn%density, &
      & qm_shell_tmp)
   self%qm_shell_xyz=qm_shell_tmp(:, :, 1)
   
   self%dipm_atom_xyz=wfn%dpat(:, :, 1)

   self%qm_atom_xyz=wfn%qpat(:, :, 1)

   call comp_norm(calc%bas%nsh, dipm_shell_tmp, qm_shell_tmp, self%dipm_shell, self%qm_shell)
   call comp_norm(mol%nat, wfn%dpat, wfn%qpat, self%dipm_atom, self%qm_atom)
   call resolve_shellwise(self%mulliken_shell, tmp_s_array, tmp_p_array, tmp_d_array, calc%bas%nsh_at, mol%nat)
   call self%dict%add_entry("p_s", tmp_s_array)
   call self%dict%add_entry("p_p", tmp_p_array)
   call self%dict%add_entry("p_d", tmp_d_array)
   call resolve_shellwise(self%dipm_shell, tmp_s_array, tmp_p_array, tmp_d_array, calc%bas%nsh_at, mol%nat)
   call self%dict%add_entry("dipm_s", tmp_s_array)
   call self%dict%add_entry("dipm_p", tmp_p_array)
   call self%dict%add_entry("dipm_d", tmp_d_array)
   
   call resolve_xyz_shell(self%dipm_shell_xyz, tmp_array, calc%bas%nsh_at, mol%nat)
   tmp_labels = [ character(len=20) ::
   &"dipm_s_x", "dipm_s_y", "dipm_s_z",&
   &"dipm_p_x", "dipm_p_y", "dipm_p_z",&
   &"dipm_d_x", "dipm_d_y", "dipm_d_z",&
   ]
   do i = 1, size(tmp_labels)
      call self%dict%add_entry(trim(tmp_labels(i)), tmp_array(i))
   end do
   call resolve_shellwise(self%qm_shell, tmp_s_array, tmp_p_array, tmp_d_array, calc%bas%nsh_at, mol%nat)
   call self%dict%add_entry("qm_s", tmp_s_array)
   call self%dict%add_entry("qm_p", tmp_p_array)
   call self%dict%add_entry("qm_d", tmp_d_array)
   call resolve_xyz_shell(self%qm_shell_xyz, tmp_array, calc%bas%nsh_at, mol%nat)
   tmp_labels = [ character(len=20) ::
   &"qm_s_xx", "qm_s_xy", "qm_s_yy", "qm_s_xz", "qm_s_yz", "qm_s_zz",&
   &"qm_p_xx", "qm_p_xy", "qm_p_yy", "qm_p_xz", "qm_p_yz", "qm_p_zz",&
   &"qm_d_xx", "qm_d_xy", "qm_d_yy", "qm_d_xz", "qm_d_yz", "qm_d_zz",& 
   ]
   do i = 1, size(tmp_labels)
      call self%dict%add_entry(trim(tmp_labels(i)), tmp_array(i))
   end do
   call self%dict%add_entry("q_A", self%partial_charge_atom)
   call self%dict%add_entry("dipm_A", self%dipm_atom)
   tmp_labels = [ character(len=20) ::
   &"dipm_A_x", "dipm_A_y", "dipm_A_z",&
   ]
   do i = 1, size(tmp_labels)
      call self%dict%add_entry(trim(tmp_labels(i)), self%dimp_atom_xyz(i, :))
   end do
   call self%dict%add_entry("qm_A", self%dipm_atom)
   tmp_labels = [ character(len=20) ::
   &"qm_A_xx", "qm_A_xy", "qm_A_yy", "qm_A_xz", "qm_A_yz", "qm_A_zz",&
   ]
   do i = 1, size(tmp_labels)
      call self%dict%add_entry(trim(tmp_labels(i)), self%qm_atom_xyz(i, :))
   end do

   self%n_features = self%n_features + features
end subroutine

subroutine mol_set_nuclear_charge(nat, at, id, z)
   implicit none
   integer, intent(in)::nat, at(nat), id(nat)
   real(wp), intent(out) :: z(nat)
   integer :: i
   do i = 1, nat
      z(i) = real(at(id(i)), wp) - real(ncore(at(id(i))))
      if (at(i) > 57 .and. at(i) < 72) z(i) = 3.0_wp
   end do
contains

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

subroutine resolve_shellwise(shell_prop, array_s, array_p, array_d, at2nsh, nat)
   real(wp), allocatable, intent(out) :: array_s(nat), array_p(nat), array_d(nat) 
   real(wp), intent(in) :: shell_prop(:)
   integer, intent(in) ::  nat, at2nsh(:)
   integer :: nsh, i
   nsh = 1
   if (allocated(array_s)) deallocate(array_s)
   if (allocated(array_p)) deallocate(array_p)
   if (allocated(array_d)) deallocate(array_d)
   allocate(array_s(nat), array_p(nat), array_d(nat))
   do i = 1, nat
      array_s(i) = shell_prop(nsh) !s shell always filled
      nsh = nsh + 1
      if (at2nsh(i) == 2) then
         array_p(i) = shell_prop(nsh)
         nsh = nsh + 1
      elseif (at2nsh(i) == 3) then
         array_p(i) = shell_prop(nsh)
         nsh = nsh + 1
         array_d(i) = shell_prop(nsh)
         nsh = nsh + 1
      end if
   end do
end subroutine

subroutine resolve_xyz_shell(mult_xyz, array, at2nsh, nat)
   real(wp), intent(in) :: mult_xyz(:, :)
   real(wp), allocatable :: array(:, :)
   integer, intent(in) :: at2nsh(:)
   integer :: j, k, nsh, id_tmp, nat, i
   nsh = 1

   if (allocated(array)) deallocate(array)
   allocate(array(nat, size(mult_xyz, dim=1)), source = 0.0_wp)

   do k = 1, nat
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

subroutine compute_extended(self, mol, wfn, integrals, calc, cache, prlevel, ctx)
   class(xtbml_feature_type), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction strcuture data
   type(wavefunction_type), intent(in) :: wfn
   !> Integral container
   type(integral_type) :: integrals
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Container
   type(container_cache), intent(inout) :: cache
   !> Context type
   type(context_type),intent(inout) :: ctx
   !> Print Level
   integer, intent(in) :: prlevel
   
   
   
   self%n_features = self%n_features + ext_features
end subroutine

subroutine allocate(self, nsh, nat)
    class(xtbml_feature_type), intent(inout) :: self
    integer :: nsh, nat
    
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

subroutine allocate_extended(self, nsh, nat, n_a)
    class(xtbml_feature_type), intent(inout) :: self
    integer :: nsh, nat, n_a
    
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

subroutine mulliken_shellwise(nao, nshell, ao2shell, p, s, charges_shell)
   implicit none
   integer, intent(in) :: nao, nshell, ao2shell(:)
   real(wp), intent(in) :: s(:, :)
   real(wp), intent(in) :: p(:, :)
   real(wp), intent(out) :: charges_shell(nshell)
   integer ::  mu, nu

   charges_shell = 0.0_wp
   !$omp parallel do default(none) collapse(2)&
   !$omp private(mu,nu) shared(charges_shell,nao,ao2shell,p,s)
   do mu = 1, nao
      do nu = 1, nao
         !$omp atomic
         charges_shell(ao2shell(mu)) = charges_shell(ao2shell(mu)) + p(mu, nu)*s(nu, mu)
      end do
   end do

end subroutine

subroutine get_delta_partial(nat, n_a, atom_partial, at, xyz, cn, delta_partial)
   implicit none
   integer, INTENT(IN) :: nat, at(nat), n_a
   real(wp), INTENT(IN) :: atom_partial(nat), xyz(3, nat), cn(nat)
   real(wp), INTENT(OUT) :: delta_partial(nat, n_a)
   real(wp) :: delta_partial_tmp(nat, n_a), f_log(nat, nat, n_a)
   integer :: a, b, k
   real(wp) :: result
   f_log = 0.0_wp

   delta_partial = 0.0_wp

       
   do k =1,n_a
      do a = 1, nat
         do b = 1, nat
                    
            delta_partial(a,k) = delta_partial(a,k) + atom_partial(b) / (inv_cn_a(a,b,k)*(cn(b)+1))
               f_log(a,b,k) = 1.0_wp / (inv_cn_a(a,b,k)*(cn(b)+1))
         enddo
      enddo
   enddo
        
end subroutine

subroutine sum_up_mm(nat, nshell, aoat2, ash, dipm_shell, qm_shell, dipm_at, qm_at)
   implicit none
   integer, INTENT(IN) :: nat, nshell, aoat2(:), ash(:)
   real(wp), INTENT(IN) :: dipm_shell(:, :), qm_shell(:, :)
   real(wp), INTENT(OUT) :: dipm_at(3, nat), qm_at(6, nat)
   integer :: i

   dipm_at = 0.0_wp
   qm_at = 0.0_wp

   do i = 1, nshell
      dipm_at(:, ash(i)) = dipm_at(:, ash(i)) + dipm_shell(:, i)
      qm_at(:, ash(i)) = qm_at(:, ash(i)) + qm_shell(:, i)
   end do

end subroutine

subroutine sum_up_mulliken(nat, nshell, aoat2, ash, mull_shell, mull_at)
   implicit none
   integer, INTENT(IN) :: nat, nshell, aoat2(:), ash(:)
   real(wp), INTENT(IN) :: mull_shell(nshell)
   real(wp), INTENT(OUT) :: mull_at(nat)
   integer :: i

   mull_at = 0.0_wp

   do i = 1, nshell
      mull_at(ash(i)) = mull_at(ash(i)) + mull_shell(i)
   end do

end subroutine

subroutine get_delta_mm(nat, n_a, q, dipm, qp, at, xyz, cn, delta_dipm, delta_qp)
   implicit none
   integer, INTENT(IN) :: nat, at(nat), n_a
   real(wp), INTENT(IN) :: dipm(3, nat), xyz(3, nat), qp(6, nat), q(nat), cn(nat)
   real(wp), INTENT(OUT) :: delta_dipm(3, nat, n_a), delta_qp(6, nat, n_a)
   integer :: a, b, k
   real(wp) :: result, r_ab(3), qp_part(6, nat, n_a), sum_qm(6)

   delta_dipm = 0.0_wp
   delta_qp = 0.0_wp
   qp_part = 0.0_wp

   do k = 1, n_a
      do a = 1, nat
         do b = 1, nat

            result = inv_cn_a(a, b, k)

            r_ab = xyz(:, a) - xyz(:, b)

            delta_dipm(:, a, k) = delta_dipm(:, a, k) + (dipm(:, b) - r_ab(:)*q(b))/(result*(cn(b) + 1))
            !sorting of qp xx,xy,yy,xz,yz,zz

            delta_qp(1, a, k) = delta_qp(1, a, k) + &
               (1.5_wp*(-1*(r_ab(1)*dipm(1, b) + r_ab(1)*dipm(1, b)) + r_ab(1)*r_ab(1)*q(b)))/(result*(cn(b) + 1))

            delta_qp(2, a, k) = delta_qp(2, a, k) + &
               (1.5_wp*(-1*(r_ab(1)*dipm(2, b) + r_ab(2)*dipm(1, b)) + r_ab(1)*r_ab(2)*q(b)))/(result*(cn(b) + 1))

            delta_qp(3, a, k) = delta_qp(3, a, k) + &
               (1.5_wp*(-1*(r_ab(2)*dipm(2, b) + r_ab(2)*dipm(2, b)) + r_ab(2)*r_ab(2)*q(b)))/(result*(cn(b) + 1))

            delta_qp(4, a, k) = delta_qp(4, a, k) + &
               (1.5_wp*(-1*(r_ab(3)*dipm(1, b) + r_ab(1)*dipm(3, b)) + r_ab(1)*r_ab(3)*q(b)))/(result*(cn(b) + 1))

            delta_qp(5, a, k) = delta_qp(5, a, k) + &
               (1.5_wp*(-1*(r_ab(3)*dipm(2, b) + r_ab(2)*dipm(3, b)) + r_ab(2)*r_ab(3)*q(b)))/(result*(cn(b) + 1))

            delta_qp(6, a, k) = delta_qp(6, a, k) + &
               (1.5_wp*(-1*(r_ab(3)*dipm(3, b) + r_ab(3)*dipm(3, b)) + r_ab(3)*r_ab(3)*q(b)))/(result*(cn(b) + 1))

            qp_part(:, a, k) = qp_part(:, a, k) + &
               qp(:, b)/(result*(cn(b) + 1))
         end do
      end do
   end do

   call remove_trac_qp(nat, n_a, delta_qp, qp_part)

end subroutine

subroutine comp_norm_3(ndim, n_a, dipm, qm, dipm_norm, qm_norm)
   implicit none
   integer, INTENT(IN) :: ndim, n_a
   real(wp), INTENT(IN) :: dipm(3, ndim, n_a), qm(6, ndim, n_a)
   real(wp), INTENT(OUT) :: dipm_norm(ndim, n_a), qm_norm(ndim, n_a)
   real(wp) :: r(ndim, n_a), r2(ndim, n_a)
   INTEGER :: i

   do i = 1, ndim
      r2(i, :) = dipm(1, i, :)**2 + dipm(2, i, :)**2 + dipm(3, i, :)**2
   end do
   r = sqrt(r2)

   dipm_norm(:, :) = r(:, :)

   do i = 1, ndim
      r2(i, :) = qm(1, i, :)**2 + 2*qm(2, i, :)**2 + qm(3, i, :)**2 + 2*qm(4, i, :)**2 + 2*qm(5, i, :)**2 + qm(6, i, :)**2
   end do
   r = sqrt(r2)
   qm_norm(:, :) = r(:, :)

end subroutine

subroutine comp_norm(ndim, dipm, qm, dipm_norm, qm_norm)
   implicit none
   integer, INTENT(IN) :: ndim
   real(wp), INTENT(IN) :: dipm(3, ndim), qm(6, ndim)
   real(wp), INTENT(OUT) :: dipm_norm(ndim), qm_norm(ndim)
   real(wp) :: r(ndim), r2(ndim)
   INTEGER :: i

   do i = 1, ndim
      r2(i) = dipm(1, i)**2 + dipm(2, i)**2 + dipm(3, i)**2
   end do
   r = sqrt(r2)

   dipm_norm(:) = r(:)

   do i = 1, ndim
      r2(i) = qm(1, i)**2 + 2*qm(2, i)**2 + qm(3, i)**2 + 2*qm(4, i)**2 + 2*qm(5, i)**2 + qm(6, i)**2
   end do
   r = sqrt(r2)
   qm_norm(:) = r(:)

end subroutine

subroutine get_delta_mm_Z(nat, n_a, q, dipm, qp, at, xyz, cn, delta_dipm, delta_qp)! all effects due to the electrons are set to 0, only the distribution of positive charges is left
   implicit none
   integer, INTENT(IN) :: nat, at(nat), n_a
   real(wp), INTENT(IN) :: dipm(3, nat), xyz(3, nat), qp(6, nat), q(nat), cn(nat)
   real(wp), INTENT(OUT) :: delta_dipm(3, nat, n_a), delta_qp(6, nat, n_a)
   integer :: a, b, k
   real(wp) :: result, r_ab(3), qp_part(6, nat, n_a)

   delta_dipm = 0.0_wp
   delta_qp = 0.0_wp
   qp_part = 0.0_wp

   do k = 1, n_a
      do a = 1, nat
         do b = 1, nat
            result = inv_cn_a(a, b, k)
            r_ab = xyz(:, a) - xyz(:, b)
            delta_dipm(:, a, k) = delta_dipm(:, a, k) + (-r_ab(:)*q(b))/(result*(cn(b) + 1))
            !sorting of qp xx,xy,yy,xz,yz,zz

            delta_qp(1, a, k) = delta_qp(1, a, k) + &
               (1.5_wp*(r_ab(1)*r_ab(1)*q(b)))/(result*(cn(b) + 1))

            delta_qp(2, a, k) = delta_qp(2, a, k) + &
               (1.5_wp*(r_ab(1)*r_ab(2)*q(b)))/(result*(cn(b) + 1))

            delta_qp(3, a, k) = delta_qp(3, a, k) + &
               (1.5_wp*(r_ab(2)*r_ab(2)*q(b)))/(result*(cn(b) + 1))

            delta_qp(4, a, k) = delta_qp(4, a, k) + &
               (1.5_wp*(r_ab(1)*r_ab(3)*q(b)))/(result*(cn(b) + 1))

            delta_qp(5, a, k) = delta_qp(5, a, k) + &
               (1.5_wp*(r_ab(2)*r_ab(3)*q(b)))/(result*(cn(b) + 1))

            delta_qp(6, a, k) = delta_qp(6, a, k) + &
               (1.5_wp*(r_ab(3)*r_ab(3)*q(b)))/(result*(cn(b) + 1))
            !qp_part(:,a) = qp_part(:,a) + qp(:,b) / (result*(cn(b)+1))
         end do
         !delta_dipm(:,a) = dipm(:,a) + delta_dipm(:,a)
         !delta_qp(:,a) = qp(:,a) + delta_qp(:,a)
      end do
   end do

   call remove_trac_qp(nat, n_a, delta_qp, qp_part)

end subroutine

subroutine get_delta_mm_p(nat, n_a, q, dipm, qp, at, xyz, cn, delta_dipm, delta_qp) ! the sign of q was changed to respect the charge of the electrons
   implicit none
   integer, INTENT(IN) :: nat, at(nat), n_a
   real(wp), INTENT(IN) :: dipm(3, nat), xyz(3, nat), qp(6, nat), q(nat), cn(nat)
   real(wp), INTENT(OUT) :: delta_dipm(3, nat, n_a), delta_qp(6, nat, n_a)
   integer :: a, b, k
   real(wp) :: result, r_ab(3), qp_part(6, nat, n_a)

   delta_dipm = 0.0_wp
   delta_qp = 0.0_wp
   qp_part = 0.0_wp

   do k = 1, n_a
      do a = 1, nat
         do b = 1, nat
            result = inv_cn_a(a, b, k)
            r_ab = xyz(:, a) - xyz(:, b)

            delta_dipm(:, a, k) = delta_dipm(:, a, k) + &
               (dipm(:, b) + r_ab(:)*q(b))/(result*(cn(b) + 1))
            !sorting of qp xx,xy,yy,xz,yz,zz

            delta_qp(1, a, k) = delta_qp(1, a, k) + &
               (1.5_wp*(-1*(r_ab(1)*dipm(1, b) + r_ab(1)*dipm(1, b)) - r_ab(1)*r_ab(1)*q(b)))/(result*(cn(b) + 1))

            delta_qp(2, a, k) = delta_qp(2, a, k) + &
               (1.5_wp*(-1*(r_ab(1)*dipm(2, b) + r_ab(2)*dipm(1, b)) - r_ab(1)*r_ab(2)*q(b)))/(result*(cn(b) + 1))

            delta_qp(3, a, k) = delta_qp(3, a, k) + &
               (1.5_wp*(-1*(r_ab(2)*dipm(2, b) + r_ab(2)*dipm(2, b)) - r_ab(2)*r_ab(2)*q(b)))/(result*(cn(b) + 1))

            delta_qp(4, a, k) = delta_qp(4, a, k) + &
               (1.5_wp*(-1*(r_ab(3)*dipm(1, b) + r_ab(1)*dipm(3, b)) - r_ab(1)*r_ab(3)*q(b)))/(result*(cn(b) + 1))

            delta_qp(5, a, k) = delta_qp(5, a, k) + &
               (1.5_wp*(-1*(r_ab(3)*dipm(2, b) + r_ab(2)*dipm(3, b)) - r_ab(2)*r_ab(3)*q(b)))/(result*(cn(b) + 1))

            delta_qp(6, a, k) = delta_qp(6, a, k) + &
               (1.5_wp*(-1*(r_ab(3)*dipm(3, b) + r_ab(3)*dipm(3, b)) - r_ab(3)*r_ab(3)*q(b)))/(result*(cn(b) + 1))
            qp_part(:, a, k) = qp_part(:, a, k) + qp(:, b)/(result*(cn(b) + 1))
         end do
         !delta_dipm(:,a) = dipm(:,a) + delta_dipm(:,a)
         !delta_qp(:,a) = qp(:,a) + delta_qp(:,a)
      end do
   end do

   call remove_trac_qp(nat, n_a, delta_qp, qp_part)

end subroutine

subroutine remove_trac_qp(nat, n_a, qp_matrix, qp_part)
   implicit none
   integer, INTENT(IN) :: nat, n_a
   real(wp), INTENT(IN) :: qp_part(6, nat, n_a)
   real(wp), INTENT(INOUT) :: qp_matrix(6, nat, n_a)
   integer :: i
   real(wp) :: tii(n_a)

   do i = 1, nat
      tii = qp_matrix(1, i, :) + qp_matrix(3, i, :) + qp_matrix(6, i, :)
      tii = tii/3.0_wp
      !qp_matrix(1:6,i) = 1.50_wp*qp(1:6,i)
      qp_matrix(1, i, :) = qp_matrix(1, i, :) - tii
      qp_matrix(3, i, :) = qp_matrix(3, i, :) - tii
      qp_matrix(6, i, :) = qp_matrix(6, i, :) - tii
      qp_matrix(:, i, :) = qp_matrix(:, i, :) + qp_part(:, i, :)
   end do
end subroutine

end module tblite_xtbml_density_based
