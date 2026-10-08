! This file is part of tblite.
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

!> @file tblite/coulomb/multipole.f90
!> Provides an implemenation of a multipole based second-order electrostatic

!> Anisotropic second-order electrostatics using a damped multipole expansion
module tblite_coulomb_multipole
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_ncoord, only : new_ncoord, ncoord_type, cn_count
   use tblite_blas, only : dot, gemv, symv, gemm
   use tblite_container_cache, only : container_cache
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_ewald, only : ewald_cache
   use tblite_coulomb_type, only : coulomb_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_partition, only : pair_list, work_partition, owns_pair, owns_index
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_wignerseitz, only : wignerseitz_cell, get_wignerseitz_weights
   implicit none
   private

   public :: new_damped_multipole


   !> Container to handle multipole electrostatics
   type, public, extends(coulomb_type) :: damped_multipole
      !> Damping function for inverse quadratic contributions
      real(wp) :: kdmp3 = 0.0_wp
      !> Damping function for inverse cubic contributions
      real(wp) :: kdmp5 = 0.0_wp
      !> Kernel for on-site dipole exchange-correlation
      real(wp), allocatable :: dkernel(:)
      !> Kernel for on-site quadrupolar exchange-correlation
      real(wp), allocatable :: qkernel(:)

      !> Shift for the generation of the multipolar damping radii
      real(wp) :: shift = 0.0_wp
      !> Exponent for the generation of the multipolar damping radii
      real(wp) :: kexp = 0.0_wp
      !> Maximum radius for the multipolar damping radii
      real(wp) :: rmax = 0.0_wp
      !> Base radii for the multipolar damping radii
      real(wp), allocatable :: rad(:)
      !> Valence coordination number
      real(wp), allocatable :: valence_cn(:)

      !> Coordination number container for multipolar damping radii
      class(ncoord_type), allocatable :: ncoord
   contains
      !> Update cache from container
      procedure :: update
      !> Return dependency on density
      procedure :: variable_info
      !> Get anisotropic electrostatic energy
      procedure :: get_energy
      !> Get anisotropic electrostatic potential
      procedure :: get_potential
      !> Get derivatives of anisotropic electrostatics
      procedure :: get_gradient
      ! These additional functions are necessary as the atomic contributions to AES are not
      ! the same in this implemntation and in the GFN2 paper,
      ! for details see https://github.com/tblite/tblite/pull/224/files#r1970341792
      !> Get only AXC part of the anisotropic electrostatics
      procedure :: get_energy_axc
      !> Get AES energy of the anisotropic electrostatics
      procedure :: get_energy_aes
   end type damped_multipole

   real(wp), parameter :: unity(3, 3) = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))
   real(wp), parameter :: conv = 100*eps
   character(len=*), parameter :: label = "anisotropic electrostatics"

contains


!> Create a new anisotropic electrostatics container
subroutine new_damped_multipole(self, mol, kdmp3, kdmp5, dkernel, qkernel, &
      & shift, kexp, rmax, rad, vcn, error)
   !> Instance of the multipole container
   type(damped_multipole), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Damping function for inverse quadratic contributions
   real(wp), intent(in) :: kdmp3
   !> Damping function for inverse cubic contributions
   real(wp), intent(in) :: kdmp5
   !> Kernel for on-site dipole exchange-correlation
   real(wp), intent(in) :: dkernel(:)
   !> Kernel for on-site quadrupolar exchange-correlation
   real(wp), intent(in) :: qkernel(:)
   !> Shift for the generation of the multipolar damping radii
   real(wp), intent(in) :: shift
   !> Exponent for the generation of the multipolar damping radii
   real(wp), intent(in) :: kexp
   !> Maximum radius for the multipolar damping radii
   real(wp), intent(in) :: rmax
   !> Base radii for the multipolar damping radii
   real(wp), intent(in) :: rad(:)
   !> Valence coordination number
   real(wp), intent(in) :: vcn(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   self%label = label
   self%kdmp3 = kdmp3
   self%kdmp5 = kdmp5
   self%dkernel = dkernel
   self%qkernel = qkernel

   self%shift = shift
   self%kexp = kexp
   self%rmax = rmax
   self%rad = rad
   self%valence_cn = vcn

   call new_ncoord(self%ncoord, mol, cn_count%dexp, error)
end subroutine new_damped_multipole


!> Update cache from container
subroutine update(self, mol, cache)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   type(coulomb_cache), pointer :: ptr
   integer :: nlocal

   call taint(cache, ptr)
   call ptr%update(mol)
   call ptr%pairs%update(self%partition, mol%nat)

   if (.not.allocated(ptr%mrad)) then
      allocate(ptr%mrad(mol%nat))
   end if
   if (.not.allocated(ptr%dmrdcn)) then
      allocate(ptr%dmrdcn(mol%nat))
   end if

   if (.not.allocated(ptr%cn)) then
      allocate(ptr%cn(mol%nat))
   end if
   if (allocated(self%ncoord)) then
      call self%ncoord%get_cn(mol, ptr%cn)
   else
      ptr%cn(:) = self%valence_cn(mol%id)
   end if
   ptr%cn_derivs_valid = .false.

   call get_mrad(mol, self%shift, self%kexp, self%rmax, self%rad, self%valence_cn, &
      & ptr%cn, ptr%mrad, ptr%dmrdcn)

   nlocal = size(ptr%pairs%neighbour)
   if (allocated(ptr%local_sd)) then
      if (size(ptr%local_sd, 3) /= nlocal) deallocate(ptr%local_sd, ptr%local_dd, ptr%local_sq)
   end if
   if (.not.allocated(ptr%local_sd)) then
      allocate(ptr%local_sd(3, 2, nlocal), ptr%local_dd(3, 3, nlocal), ptr%local_sq(6, 2, nlocal))
   end if
   call get_multipole_matrix(self, mol, ptr)
end subroutine update


!> Get anisotropic electrostatic energy
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Electrostatic energy
   real(wp), intent(inout) :: energies(:)
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   real(wp), allocatable :: vd(:, :), vq(:, :)
   type(coulomb_cache), pointer :: ptr

   call view(cache, ptr)

   allocate(vd(3, mol%nat), vq(6, mol%nat))

   vd = 0.0_wp
   vq = 0.0_wp
   call contract_local(ptr, wfn, 0.5_wp, vd, vq)

   energies(:) = energies + sum(wfn%dpat(:, :, 1) * vd, 1) + sum(wfn%qpat(:, :, 1) * vq, 1)

   call get_kernel_energy(mol, self%dkernel, wfn%dpat(:, :, 1), energies, self%partition)
   call get_kernel_energy(mol, self%qkernel, wfn%qpat(:, :, 1), energies, self%partition)
end subroutine get_energy

!> Get anisotropic electrostatic energy
subroutine get_energy_aes(self, mol, cache, wfn, energies)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Electrostatic energy
   real(wp), intent(inout) :: energies(:)
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   real(wp), allocatable :: vs(:), vd(:, :), vq(:, :)
   type(coulomb_cache), pointer :: ptr

   call view(cache, ptr)
   allocate(vs(mol%nat), vd(3, mol%nat), vq(6, mol%nat), source=0.0_wp)
   call contract_local(ptr, wfn, 1.0_wp, vd, vq, vs)
   energies = energies + 0.5_wp*(wfn%qat(:, 1)*vs + &
      & sum(wfn%dpat(:, :, 1)*vd, 1) + sum(wfn%qpat(:, :, 1)*vq, 1))
end subroutine get_energy_aes

!> Get multipolar anisotropic exchange-correlation kernel
subroutine get_kernel_energy(mol, kernel, mpat, energies, partition)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multipole kernel
   real(wp), intent(in) :: kernel(:)
   !> Atomic multipole momemnt
   real(wp), intent(in) :: mpat(:, :)
   !> Electrostatic energy
   real(wp), intent(inout) :: energies(:)
   !> Share of the on-site terms evaluated here, absent selects the complete work
   type(work_partition), intent(in), optional :: partition

   integer :: iat, izp
   real(wp) :: mpt(size(mpat, 1)), mpscale(size(mpat, 1))

   mpscale(:) = 1
   if (size(mpat, 1) == 6) mpscale([2, 4, 5]) = 2

   do iat = 1, mol%nat
      if (.not.owns_index(partition, iat)) cycle
      izp = mol%id(iat)
      mpt(:) = mpat(:, iat) * mpscale
      energies(iat) = energies(iat) + kernel(izp) * dot_product(mpt, mpat(:, iat))
   end do
end subroutine get_kernel_energy


!> Get anisotropic electrostatic potential
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   type(coulomb_cache), pointer :: ptr

   call view(cache, ptr)

   call contract_local(ptr, wfn, 1.0_wp, pot%vdp(:, :, 1), pot%vqp(:, :, 1), pot%vat(:, 1))

   call get_kernel_potential(mol, self%dkernel, wfn%dpat(:, :, 1), pot%vdp(:, :, 1), &
      & self%partition)
   call get_kernel_potential(mol, self%qkernel, wfn%qpat(:, :, 1), pot%vqp(:, :, 1), &
      & self%partition)
end subroutine get_potential


!> Contract local multipole blocks, sharing the pair traversal of all moments.
subroutine contract_local(ptr, wfn, ddscale, vd, vq, vs)
   !> Local interaction blocks and their row-wise pair list
   type(coulomb_cache), intent(in) :: ptr
   !> Atomic charges and multipole moments
   type(wavefunction_type), intent(in) :: wfn
   !> Dipole-dipole weight: one half for the multipole-only energy contraction, one otherwise
   real(wp), intent(in) :: ddscale
   !> Dipolar potential to accumulate into (3, nat)
   real(wp), intent(inout) :: vd(:, :)
   !> Quadrupolar potential to accumulate into (6, nat)
   real(wp), intent(inout) :: vq(:, :)
   !> Charge potential to accumulate into, omitted for the multipole-only energy contraction
   real(wp), intent(inout), optional :: vs(:)
   integer :: iat, jat, ipair
   real(wp) :: dval(3), qval(6), sval

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(ptr, wfn, ddscale, vd, vq, vs) &
   !$omp private(iat, jat, ipair, dval, qval, sval)
   do iat = 1, size(vd, 2)
      dval = 0.0_wp
      qval = 0.0_wp
      sval = 0.0_wp
      do ipair = ptr%pairs%offset(iat-1)+1, ptr%pairs%offset(iat)
         jat = ptr%pairs%neighbour(ipair)
         dval = dval + ptr%local_sd(:, 1, ipair)*wfn%qat(jat, 1) &
            & + ddscale*(ptr%local_dd(:, 1, ipair)*wfn%dpat(1, jat, 1) &
            & + ptr%local_dd(:, 2, ipair)*wfn%dpat(2, jat, 1) &
            & + ptr%local_dd(:, 3, ipair)*wfn%dpat(3, jat, 1))
         qval = qval + ptr%local_sq(:, 1, ipair)*wfn%qat(jat, 1)
         if (present(vs)) then
            sval = sval + dot_product(ptr%local_sd(:, 2, ipair), wfn%dpat(:, jat, 1)) &
               & + dot_product(ptr%local_sq(:, 2, ipair), wfn%qpat(:, jat, 1))
         end if
      end do
      vd(:, iat) = vd(:, iat) + dval
      vq(:, iat) = vq(:, iat) + qval
      if (present(vs)) vs(iat) = vs(iat) + sval
   end do
end subroutine contract_local


!> Get multipolar anisotropic potential contribution
subroutine get_kernel_potential(mol, kernel, mpat, vm, partition)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multipole kernel
   real(wp), intent(in) :: kernel(:)
   !> Atomic multipole momemnt
   real(wp), intent(in) :: mpat(:, :)
   !> Potential shoft on atomic multipole moment
   real(wp), intent(inout) :: vm(:, :)
   !> Share of the on-site terms evaluated here, absent selects the complete work
   type(work_partition), intent(in), optional :: partition

   integer :: iat, izp
   real(wp) :: mpscale(size(mpat, 1))

   mpscale(:) = 1
   if (size(mpat, 1) == 6) mpscale([2, 4, 5]) = 2

   do iat = 1, mol%nat
      if (.not.owns_index(partition, iat)) cycle
      izp = mol%id(iat)
      vm(:, iat) = vm(:, iat) + 2*kernel(izp) * mpat(:, iat) * mpscale
   end do
end subroutine get_kernel_potential


!> Get derivatives of anisotropic electrostatics
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular gradient of the repulsion energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the repulsion energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   ! allow(C061): see https://github.com/PlasmaFAIR/fortitude/issues/695
   real(wp), allocatable :: dEdr(:)
   type(coulomb_cache), pointer :: ptr

   call view(cache, ptr)

   if (.not.ptr%cn_derivs_valid) then
      if (.not.allocated(ptr%dcndr)) allocate(ptr%dcndr(3, mol%nat, mol%nat))
      if (.not.allocated(ptr%dcndL)) allocate(ptr%dcndL(3, 3, mol%nat))
      if (allocated(self%ncoord)) then
         call self%ncoord%get_cn(mol, ptr%cn, ptr%dcndr, ptr%dcndL)
      else
         ptr%dcndr = 0.0_wp
         ptr%dcndL = 0.0_wp
      end if
      ptr%cn_derivs_valid = .true.
   end if

   allocate(dEdr(mol%nat))
   dEdr = 0.0_wp

   call get_multipole_gradient(self, mol, ptr, &
      & wfn%qat(:, 1), wfn%dpat(:, :, 1), wfn%qpat(:, :, 1), &
      & dEdr, gradient, sigma)

   dEdr(:) = dEdr * ptr%dmrdcn

   call gemv(ptr%dcndr, dEdr, gradient, beta=1.0_wp)
   call gemv(ptr%dcndL, dEdr, sigma, beta=1.0_wp)
end subroutine get_gradient


!> Calculate multipole damping radii
subroutine get_mrad(mol, shift, kexp, rmax, rad, valence_cn, cn, mrad, dmrdcn)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Shift for the generation of the multipolar damping radii
   real(wp), intent(in) :: shift
   !> Exponent for the generation of the multipolar damping radii
   real(wp), intent(in) :: kexp
   !> Maximum radius for the multipolar damping radii
   real(wp), intent(in) :: rmax
   !> Base radii for the multipolar damping radii
   real(wp), intent(in) :: rad(:)
   !> Valence coordination number
   real(wp), intent(in) :: valence_cn(:)
   !> Coordination numbers for all atoms
   real(wp), intent(in) :: cn(:)
   !> Multipole damping radii for all atoms
   real(wp), intent(out) :: mrad(:)
   !> Derivative of multipole damping radii with repect to the coordination numbers
   real(wp), intent(out) :: dmrdcn(:)

   integer :: iat, izp
   real(wp) :: arg, t1, t2

   do iat = 1, mol%nat
      izp = mol%id(iat)
      arg = cn(iat) - valence_cn(izp) - shift
      t1 = exp(-kexp*arg)
      t2 = (rmax - rad(izp)) / (1.0_wp + t1)
      mrad(iat) = rad(izp) + t2
      dmrdcn(iat) = -t2 * kexp * t1 / (1.0_wp + t1)
   end do
end subroutine get_mrad


!> Get real lattice vectors
subroutine get_dir_trans(lattice, trans)
   !> Lattice parameters
   real(wp), intent(in) :: lattice(:, :)
   !> Translation vectors
   real(wp), allocatable, intent(out) :: trans(:, :)

   call get_lattice_points([.true.], lattice, 100.0_wp, trans)

end subroutine get_dir_trans

!> Get interaction matrix for all multipole moments up to inverse cubic order
subroutine get_multipole_matrix(self, mol, cache)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Geometry data, local pair list and interaction blocks to overwrite
   type(coulomb_cache), intent(inout) :: cache

   cache%local_sd = 0.0_wp
   cache%local_dd = 0.0_wp
   cache%local_sq = 0.0_wp
   if (any(mol%periodic)) then
      call cache%multipole_ewald%update(mol%lattice, cache%alpha_multipole, conv)
      call get_multipole_matrix_3d(mol, cache%mrad, self%kdmp3, self%kdmp5, &
         & cache%wsc, cache%alpha_multipole, cache%multipole_ewald, cache%pairs, &
         & cache%local_sd, cache%local_dd, cache%local_sq)
   else
      call get_multipole_matrix_0d(mol, cache%mrad, self%kdmp3, self%kdmp5, &
         & cache%pairs, cache%local_sd, cache%local_dd, cache%local_sq)
   end if
end subroutine get_multipole_matrix

!> Calculate the multipole interaction matrix for finite systems
!> Overwrite off-site blocks; the caller initializes on-site blocks to zero.
subroutine get_multipole_matrix_0d(mol, rad, kdmp3, kdmp5, pairs, amat_sd, amat_dd, amat_sq)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multipole damping radii for all atoms
   real(wp), intent(in) :: rad(:)
   !> Damping function for inverse quadratic contributions
   real(wp), intent(in) :: kdmp3
   !> Damping function for inverse cubic contributions
   real(wp), intent(in) :: kdmp5
   !> Row-wise local pairs, including both orientations of each owned pair
   type(pair_list), intent(in) :: pairs
   !> Charge-dipole blocks (3, 2, npair); multipole on row atom (1) or neighbour (2)
   real(wp), intent(inout) :: amat_sd(:, :, :)
   !> Dipole-dipole blocks (3, 3, npair), with row-atom components first
   real(wp), intent(inout) :: amat_dd(:, :, :)
   !> Charge-quadrupole blocks (6, 2, npair), with the same directions as amat_sd
   real(wp), intent(inout) :: amat_sq(:, :, :)

   integer :: iat, jat, ipair, jpair
   real(wp) :: r1, vec(3), g1, g3, g5, fdmp3, fdmp5, tc(6), rr

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, rad, kdmp3, kdmp5, pairs, amat_sd, amat_dd, amat_sq) &
   !$omp private(iat, jat, r1, vec, g1, g3, g5, fdmp3, fdmp5, tc, rr, ipair, jpair)
   do iat = 1, mol%nat
      do jpair = pairs%offset(iat-1)+1, pairs%offset(iat)
         jat = pairs%neighbour(jpair)
         if (iat == jat) cycle
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         r1 = norm2(vec)
         g1 = 1.0_wp / r1
         g3 = g1 * g1 * g1
         g5 = g3 * g1 * g1

         rr = 0.5_wp * (rad(jat) + rad(iat)) * g1
         fdmp3 = 1.0_wp / (1.0_wp + 6.0_wp * rr**kdmp3)
         fdmp5 = 1.0_wp / (1.0_wp + 6.0_wp * rr**kdmp5)

         ipair = pairs%find(jat, iat)
         amat_sd(:, 1, ipair) = vec*g3*fdmp3
         amat_sd(:, 2, jpair) = amat_sd(:, 1, ipair)
         amat_dd(:, :, ipair) = unity*g3*fdmp5 &
            & - spread(vec, 1, 3)*spread(vec, 2, 3)*3*g5*fdmp5
         tc(2) = 2*vec(1)*vec(2)*g5*fdmp5
         tc(4) = 2*vec(1)*vec(3)*g5*fdmp5
         tc(5) = 2*vec(2)*vec(3)*g5*fdmp5
         tc(1) = vec(1)*vec(1)*g5*fdmp5
         tc(3) = vec(2)*vec(2)*g5*fdmp5
         tc(6) = vec(3)*vec(3)*g5*fdmp5
         amat_sq(:, 1, ipair) = tc
         amat_sq(:, 2, jpair) = tc
      end do
   end do
end subroutine get_multipole_matrix_0d

!> Evaluate multipole interaction matrix under 3D periodic boundary conditions
!> Accumulate Ewald and self-interaction terms into the local blocks.
subroutine get_multipole_matrix_3d(mol, rad, kdmp3, kdmp5, wsc, alpha, ewald, pairs, amat_sd, amat_dd, amat_sq)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multipole damping radii for all atoms
   real(wp), intent(in) :: rad(:)
   !> Damping function for inverse quadratic contributions
   real(wp), intent(in) :: kdmp3
   !> Damping function for inverse cubic contributions
   real(wp), intent(in) :: kdmp5
   !> Wigner-Seitz cell images
   type(wignerseitz_cell), intent(in) :: wsc
   !> Convergence parameter for Ewald sum
   real(wp), intent(in) :: alpha
   !> Cached reciprocal-space vectors and Ewald coefficients
   type(ewald_cache), intent(in) :: ewald
   !> Row-wise local pairs, including both orientations of each owned pair
   type(pair_list), intent(in) :: pairs
   !> Charge-dipole blocks (3, 2, npair); multipole on row atom (1) or neighbour (2)
   real(wp), intent(inout) :: amat_sd(:, :, :)
   !> Dipole-dipole blocks (3, 3, npair), with row-atom components first
   real(wp), intent(inout) :: amat_dd(:, :, :)
   !> Charge-quadrupole blocks (6, 2, npair), with the same directions as amat_sd
   real(wp), intent(inout) :: amat_sq(:, :, :)

   integer :: iat, jat, img, k, ipair, jpair
   real(wp) :: vec(3), rij(3), rr
   real(wp) :: d_sd(3), d_dd(3, 3), d_sq(6), r_sd(3), r_dd(3, 3), r_sq(6)
   real(wp) :: weight(size(wsc%tridx, 1))
   real(wp), allocatable :: dtrans(:, :)

   call get_dir_trans(mol%lattice, dtrans)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, wsc, rad, alpha, ewald, dtrans, kdmp3, kdmp5, pairs, amat_sd, amat_dd, amat_sq) &
   !$omp private(iat, jat, img, vec, rij, rr, weight, d_sd, d_dd, d_sq, r_sd, r_dd, r_sq, ipair, jpair)
   do iat = 1, mol%nat
      do jpair = pairs%offset(iat-1)+1, pairs%offset(iat)
         jat = pairs%neighbour(jpair)
         ipair = pairs%find(jat, iat)
         rij = mol%xyz(:, iat) - mol%xyz(:, jat)
         call get_wignerseitz_weights(wsc, jat, iat, rij, weight)
         do img = 1, wsc%nimg(jat, iat)
            vec = rij - wsc%trans(:, wsc%tridx(img, jat, iat))

            rr = 0.5_wp * (rad(jat) + rad(iat))
            call get_amat_sdq_rec_3d(vec, ewald, r_sd, r_dd, r_sq)
            call get_amat_sdq_dir_3d(vec, rr, kdmp3, kdmp5, alpha, dtrans, d_sd, d_dd, d_sq)

            amat_sd(:, 1, ipair) = amat_sd(:, 1, ipair) + weight(img)*(d_sd+r_sd)
            amat_dd(:, :, ipair) = amat_dd(:, :, ipair) + weight(img)*(r_dd+d_dd)
            amat_sq(:, 1, ipair) = amat_sq(:, 1, ipair) + weight(img)*(r_sq+d_sq)
            amat_sd(:, 2, jpair) = amat_sd(:, 2, jpair) + weight(img)*(d_sd+r_sd)
            amat_sq(:, 2, jpair) = amat_sq(:, 2, jpair) + weight(img)*(r_sq+d_sq)
         end do
      end do
   end do

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, alpha, pairs, amat_dd, amat_sq) private(iat, rr, k, ipair)
   do iat = 1, mol%nat
      ipair = pairs%find(iat, iat)
      if (ipair == 0) cycle
      ! dipole-dipole selfenergy: -2/3·α³/sqrt(π) Σ(i) μ²(i)
      rr = -2.0_wp/3.0_wp * alpha**3 / sqrtpi
      do k = 1, 3
         amat_dd(k, k, ipair) = amat_dd(k, k, ipair) + 2*rr
      end do

      ! charge-quadrupole selfenergy: 4/9·α³/sqrt(π) Σ(i) q(i)Tr(θi)
      ! (no actual contribution since quadrupoles are traceless)
      rr = 4.0_wp/9.0_wp * alpha**3 / sqrtpi
      amat_sq([1, 3, 6], :, ipair) = amat_sq([1, 3, 6], :, ipair) + rr
   end do
end subroutine get_multipole_matrix_3d

pure subroutine get_amat_sdq_rec_3d(rij, ewald, amat_sd, amat_dd, amat_sq)
   real(wp), intent(in) :: rij(3)
   type(ewald_cache), intent(in) :: ewald
   real(wp), intent(out) :: amat_sd(:)
   real(wp), intent(out) :: amat_dd(:, :)
   real(wp), intent(out) :: amat_sq(:)

   integer :: itr
   real(wp) :: vec(3), sink, cosk, gv, k_q(6)

   amat_sd = 0.0_wp
   amat_dd = 0.0_wp
   amat_sq = 0.0_wp
   do itr = 1, size(ewald%weight)
      vec = ewald%vec(:, itr)
      gv = dot_product(rij, vec)
      sink = sin(gv)*ewald%weight(itr)
      cosk = cos(gv)*ewald%weight(itr)

      ! packed quadratic basis
      k_q(1) =          vec(1) * vec(1)
      k_q(2) = 2.0_wp * vec(1) * vec(2)
      k_q(3) =          vec(2) * vec(2)
      k_q(4) = 2.0_wp * vec(1) * vec(3)
      k_q(5) = 2.0_wp * vec(2) * vec(3)
      k_q(6) =          vec(3) * vec(3)

      amat_sd(:) = amat_sd + vec * sink
      amat_dd(:, 1) = amat_dd(:, 1) + vec * vec(1) * cosk
      amat_dd(:, 2) = amat_dd(:, 2) + vec * vec(2) * cosk
      amat_dd(:, 3) = amat_dd(:, 3) + vec * vec(3) * cosk
      amat_sq(:) = amat_sq - cosk / 3.0_wp * k_q
   end do

end subroutine get_amat_sdq_rec_3d

pure subroutine get_amat_sdq_dir_3d(rij, rr, kdmp3, kdmp5, alp, trans, &
      & amat_sd, amat_dd, amat_sq)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: rr
   real(wp), intent(in) :: kdmp3
   real(wp), intent(in) :: kdmp5
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: amat_sd(:)
   real(wp), intent(out) :: amat_dd(:, :)
   real(wp), intent(out) :: amat_sq(:)

   integer :: itr
   real(wp) :: vec(3), r1, tmp, fdmp3, fdmp5, g1, g3, g5, arg, arg2, alp2, e1, e2, erft, expt

   amat_sd = 0.0_wp
   amat_dd = 0.0_wp
   amat_sq = 0.0_wp
   alp2 = alp*alp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      g1 = 1.0_wp/r1
      g3 = g1 * g1 * g1
      g5 = g3 * g1 * g1
      fdmp3 = 1.0_wp / (1.0_wp + 6.0_wp * (rr/r1)**kdmp3)
      fdmp5 = 1.0_wp / (1.0_wp + 6.0_wp * (rr/r1)**kdmp5)

      arg = r1*alp
      arg2 = arg*arg
      expt = exp(-arg2)/sqrtpi
      erft = -erf(arg)*g1
      e1 = g1*g1 * (erft + 2*expt*alp)
      e2 = g1*g1 * (e1 + 4*expt*alp2*alp/3)

      tmp = fdmp3 * g3 + e1
      amat_sd = amat_sd + vec * tmp
      tmp = fdmp5 * g3 + e1
      amat_dd(1, 1) = amat_dd(1, 1) + tmp
      amat_dd(2, 2) = amat_dd(2, 2) + tmp
      amat_dd(3, 3) = amat_dd(3, 3) + tmp
      tmp = fdmp5 * g5 + e2
      amat_dd(:, 1) = amat_dd(:, 1) - vec * vec(1) * (3 * tmp)
      amat_dd(:, 2) = amat_dd(:, 2) - vec * vec(2) * (3 * tmp)
      amat_dd(:, 3) = amat_dd(:, 3) - vec * vec(3) * (3 * tmp)
      amat_sq(1) = amat_sq(1) +   vec(1)*vec(1)*(g5*fdmp5 + e2) - (fdmp5*g3 + e1)/3.0_wp
      amat_sq(2) = amat_sq(2) + 2*vec(1)*vec(2)*(g5*fdmp5 + e2)
      amat_sq(3) = amat_sq(3) +   vec(2)*vec(2)*(g5*fdmp5 + e2) - (fdmp5*g3 + e1)/3.0_wp
      amat_sq(4) = amat_sq(4) + 2*vec(1)*vec(3)*(g5*fdmp5 + e2)
      amat_sq(5) = amat_sq(5) + 2*vec(2)*vec(3)*(g5*fdmp5 + e2)
      amat_sq(6) = amat_sq(6) +   vec(3)*vec(3)*(g5*fdmp5 + e2) - (fdmp5*g3 + e1)/3.0_wp
   end do

end subroutine get_amat_sdq_dir_3d


!> Calculate derivatives of multipole interactions
subroutine get_multipole_gradient(self, mol, cache, qat, dpat, qpat, dEdr, gradient, sigma)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   !> Atomic partial charges
   real(wp), contiguous, intent(in) :: qat(:)
   !> Atomic dipole moments
   real(wp), contiguous, intent(in) :: dpat(:, :)
   !> Atomic quadrupole moments
   real(wp), contiguous, intent(in) :: qpat(:, :)
   !> Derivative of the energy w.r.t. the critical radii
   real(wp), contiguous, intent(inout) :: dEdr(:)
   !> Derivative of the energy w.r.t. atomic displacements
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Derivative of the energy w.r.t. strain deformations
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   if (any(mol%periodic)) then
      call cache%multipole_ewald%update(mol%lattice, cache%alpha_multipole, conv)
      call get_multipole_gradient_3d(mol, cache%mrad, self%kdmp3, self%kdmp5, &
         & qat, dpat, qpat, cache%wsc, cache%alpha_multipole, cache%multipole_ewald, dEdr, gradient, sigma, &
         & self%partition)
   else
      call get_multipole_gradient_0d(mol, cache%mrad, self%kdmp3, self%kdmp5, &
         & qat, dpat, qpat, dEdr, gradient, sigma, self%partition)
   end if
end subroutine get_multipole_gradient

!> Evaluate multipole derivatives for finite systems
subroutine get_multipole_gradient_0d(mol, rad, kdmp3, kdmp5, qat, dpat, qpat, &
      & dEdr, gradient, sigma, partition)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multipole damping radii for all atoms
   real(wp), intent(in) :: rad(:)
   !> Damping function for inverse quadratic contributions
   real(wp), intent(in) :: kdmp3
   !> Damping function for inverse cubic contributions
   real(wp), intent(in) :: kdmp5
   !> Atomic partial charges
   real(wp), contiguous, intent(in) :: qat(:)
   !> Atomic dipole moments
   real(wp), contiguous, intent(in) :: dpat(:, :)
   !> Atomic quadrupole moments
   real(wp), contiguous, intent(in) :: qpat(:, :)
   !> Derivative of the energy w.r.t. the critical radii
   real(wp), contiguous, intent(inout) :: dEdr(:)
   !> Derivative of the energy w.r.t. atomic displacements
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Derivative of the energy w.r.t. strain deformations
   real(wp), contiguous, intent(inout) :: sigma(:, :)
   !> Share of the atom pairs evaluated here, absent selects the complete work
   type(work_partition), intent(in), optional :: partition

   integer :: iat, jat
   real(wp) :: r1, r2, vec(3), rr, fdmp3, fdmp5, g1, g3, g5, g7, dG(3), dS(3, 3)
   real(wp) :: ddmp3, ddmp5, fddr, eq, edd, dpidpj, dpiv, dpjv, dpiqj, qidpj

   !$omp parallel do default(none) schedule(static, 1) reduction(+:dEdr, gradient, sigma) &
   !$omp shared(mol, kdmp3, kdmp5, rad, qat, dpat, qpat, partition) &
   !$omp private(iat, jat, r1, r2, vec, rr, fdmp3, fdmp5, g1, g3, g5, g7, dG, dS, &
   !$omp& ddmp3, ddmp5, fddr, eq, edd, dpidpj, dpiv, dpjv, dpiqj, qidpj)
   do iat = 1, mol%nat
      do jat = 1, iat - 1
         if (.not.owns_pair(partition, iat, jat)) cycle
         rr = 0.5_wp*(rad(iat)+rad(jat))
         vec(:) = mol%xyz(:, jat)-mol%xyz(:, iat)
         r1 = norm2(vec)
         r2 = r1 * r1
         g1 = 1.0_wp/r1
         g3 = g1 * g1 * g1
         g5 = g3 * g1 * g1
         g7 = g5 * g1 * g1

         fdmp3 = 1.0_wp/(1.0_wp+6.0_wp*(rr*g1)**kdmp3)
         ddmp3 = -3*g5*fdmp3 - kdmp3*fdmp3*(fdmp3-1.0_wp)*g5
         fdmp5 = 1.0_wp/(1.0_wp+6.0_wp*(rr*g1)**kdmp5)
         ddmp5 = -5*fdmp5 - kdmp5*(fdmp5*fdmp5-fdmp5)

         dpiqj = dot_product(vec, dpat(:, iat))*qat(jat)
         qidpj = dot_product(vec, dpat(:, jat))*qat(iat)
         fddr = 3.0_wp*(dpiqj - qidpj)*kdmp3*fdmp3*g3*(fdmp3/rr)*(rr*g1)**kdmp3
         dg(:) = - ddmp3*vec * (dpiqj - qidpj) &
            & + fdmp3*g3*(qat(iat)*dpat(:, jat) - qat(jat)*dpat(:, iat))
         ds(:, :) = - 0.5_wp * (spread(vec, 1, 3) * spread(dG, 2, 3) &
            & + spread(dG, 1, 3) * spread(vec, 2, 3))

         dEdr(iat) = dEdr(iat) + fddr
         dEdr(jat) = dEdr(jat) + fddr
         gradient(:, iat) = gradient(:, iat) + dG
         gradient(:, jat) = gradient(:, jat) - dG
         sigma(:, :) = sigma + dS

         dpidpj = dot_product(dpat(:, jat), dpat(:, iat))
         dpiv = dot_product(dpat(:, iat), vec)
         dpjv = dot_product(dpat(:, jat), vec)
         edd = dpidpj*r2 - 3*dpjv*dpiv
         fddr = 3.0_wp*edd*kdmp5*fdmp5*g5*(fdmp5/rr)*(rr*g1)**kdmp5
         dg(:) = - 2.0_wp*fdmp5*g5*dpidpj*vec &
            & + 3.0_wp*fdmp5*g5*(dpiv*dpat(:, jat) + dpjv*dpat(:, iat)) &
            & - edd*ddmp5*g7*vec
         ds(:, :) = - 0.5_wp * (spread(vec, 1, 3) * spread(dg, 2, 3) &
            & + spread(dg, 1, 3) * spread(vec, 2, 3))

         dEdr(iat) = dEdr(iat) + fddr
         dEdr(jat) = dEdr(jat) + fddr
         gradient(:, iat) = gradient(:, iat) + dg
         gradient(:, jat) = gradient(:, jat) - dg
         sigma(:, :) = sigma + ds

         eq = &
            & + 2*(qat(jat)*qpat(2,iat) + qpat(2,jat)*qat(iat))*vec(1)*vec(2) &
            & + 2*(qat(jat)*qpat(4,iat) + qpat(4,jat)*qat(iat))*vec(1)*vec(3) &
            & + 2*(qat(jat)*qpat(5,iat) + qpat(5,jat)*qat(iat))*vec(2)*vec(3) &
            & + (qat(jat)*qpat(1,iat) + qpat(1,jat)*qat(iat))*vec(1)*vec(1) &
            & + (qat(jat)*qpat(3,iat) + qpat(3,jat)*qat(iat))*vec(2)*vec(2) &
            & + (qat(jat)*qpat(6,iat) + qpat(6,jat)*qat(iat))*vec(3)*vec(3)

         fddr = eq * 3.0_wp*kdmp5*fdmp5*g5*fdmp5/rr*(rr*g1)**kdmp5
         dg(:) = - eq*ddmp5*g7*vec &
            & - 2.0_wp*fdmp5*g5*qat(iat) * &
            &[vec(1)*qpat(1,jat) + vec(2)*qpat(2,jat) + vec(3)*qpat(4,jat), &
            & vec(1)*qpat(2,jat) + vec(2)*qpat(3,jat) + vec(3)*qpat(5,jat), &
            & vec(1)*qpat(4,jat) + vec(2)*qpat(5,jat) + vec(3)*qpat(6,jat)] &
            & - 2.0_wp*fdmp5*g5*qat(jat) * &
            &[vec(1)*qpat(1,iat) + vec(2)*qpat(2,iat) + vec(3)*qpat(4,iat), &
            & vec(1)*qpat(2,iat) + vec(2)*qpat(3,iat) + vec(3)*qpat(5,iat), &
            & vec(1)*qpat(4,iat) + vec(2)*qpat(5,iat) + vec(3)*qpat(6,iat)]
         ds(:, :) = - 0.5_wp * (spread(vec, 1, 3) * spread(dG, 2, 3) &
            & + spread(dG, 1, 3) * spread(vec, 2, 3))

         dEdr(iat) = dEdr(iat) + fddr
         dEdr(jat) = dEdr(jat) + fddr
         gradient(:, iat) = gradient(:, iat) + dg
         gradient(:, jat) = gradient(:, jat) - dg
         sigma(:, :) = sigma + ds
      end do
   end do
end subroutine get_multipole_gradient_0d

!> Evaluate multipole derivatives under 3D periodic boundary conditions
subroutine get_multipole_gradient_3d(mol, rad, kdmp3, kdmp5, qat, dpat, qpat, wsc, alpha, ewald, &
      & dEdr, gradient, sigma, partition)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multipole damping radii for all atoms
   real(wp), intent(in) :: rad(:)
   !> Damping function for inverse quadratic contributions
   real(wp), intent(in) :: kdmp3
   !> Damping function for inverse cubic contributions
   real(wp), intent(in) :: kdmp5
   !> Atomic partial charges
   real(wp), contiguous, intent(in) :: qat(:)
   !> Atomic dipole moments
   real(wp), contiguous, intent(in) :: dpat(:, :)
   !> Atomic quadrupole moments
   real(wp), contiguous, intent(in) :: qpat(:, :)
   !> Wigner-Seitz cell images
   type(wignerseitz_cell), intent(in) :: wsc
   !> Convergence parameter for Ewald sum
   real(wp), intent(in) :: alpha
   type(ewald_cache), intent(in) :: ewald
   !> Derivative of the energy w.r.t. the critical radii
   real(wp), contiguous, intent(inout) :: dEdr(:)
   !> Derivative of the energy w.r.t. atomic displacements
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Derivative of the energy w.r.t. strain deformations
   real(wp), contiguous, intent(inout) :: sigma(:, :)
   !> Share of the atom pairs evaluated here, absent selects the complete work
   type(work_partition), intent(in), optional :: partition

   integer :: iat, jat, img
   real(wp) :: vec(3), rij(3), dG(3), dGr(3), dGd(3), dS(3, 3), dSr(3, 3), dSd(3, 3)
   real(wp) :: rr, dE, dEd, eimg, erec, edir
   real(wp) :: weight(size(wsc%tridx, 1)), dwdr(3, size(wsc%tridx, 1))
   real(wp) :: dwdL(3, 3, size(wsc%tridx, 1))
   real(wp), allocatable :: dtrans(:, :)

   call get_dir_trans(mol%lattice, dtrans)

   !$omp parallel do default(none) schedule(static, 1) reduction(+:dEdr, gradient, sigma) &
   !$omp shared(mol, wsc, kdmp3, kdmp5, alpha, ewald, dtrans, rad, qat, dpat, qpat) &
   !$omp shared(partition) &
   !$omp private(iat, jat, dE, dG, dS, img, vec, rij, rr, dEd, dGd, dGr, dSd, dSr) &
   !$omp private(eimg, erec, edir, weight, dwdr, dwdL)
   do iat = 1, mol%nat
      do jat = 1, iat - 1
         if (.not.owns_pair(partition, iat, jat)) cycle
         dE = 0.0_wp
         dG(:) = 0.0_wp
         dS(:, :) = 0.0_wp
         rij = mol%xyz(:, iat) - mol%xyz(:, jat)
         call get_wignerseitz_weights(wsc, jat, iat, rij, weight, dwdr, dwdL)
         do img = 1, wsc%nimg(jat, iat)
            vec(:) = -rij + wsc%trans(:, wsc%tridx(img, jat, iat))
            rr = 0.5_wp * (rad(jat) + rad(iat))

            call get_damat_sdq_rec_3d(vec, qat(iat), qat(jat), dpat(:, iat), dpat(:, jat), &
               & qpat(:, iat), qpat(:, jat), ewald, dGr, dSr, erec)
            call get_damat_sdq_dir_3d(vec, qat(iat), qat(jat), dpat(:, iat), dpat(:, jat), &
               & qpat(:, iat), qpat(:, jat), rr, kdmp3, kdmp5, alpha, dtrans, dEd, dGd, dSd, edir)
            eimg = erec + edir
            dE = dE + dEd * weight(img)
            dG = dG + (dGd + dGr) * weight(img) + eimg*dwdr(:, img)
            dS = dS + (dSd + dSr) * weight(img) + eimg*dwdL(:, :, img)
         end do
         dEdr(iat) = dEdr(iat) + dE
         dEdr(jat) = dEdr(jat) + dE
         gradient(:, iat) = gradient(:, iat) + dG
         gradient(:, jat) = gradient(:, jat) - dG
         sigma = sigma + dS
      end do
   end do

   !$omp parallel do default(none) schedule(static, 1) reduction(+:dEdr, sigma) &
   !$omp shared(mol, wsc, kdmp3, kdmp5, alpha, ewald, dtrans, rad, qat, dpat, qpat) &
   !$omp shared(partition) &
   !$omp private(iat, jat, dE, dG, dS, img, vec, rij, rr, dEd, dGd, dGr, dSd, dSr) &
   !$omp private(eimg, erec, edir, weight, dwdr, dwdL)
   do iat = 1, mol%nat
      if (.not.owns_pair(partition, iat, iat)) cycle
      dE = 0.0_wp
      dS(:, :) = 0.0_wp
      rij(:) = 0.0_wp
      call get_wignerseitz_weights(wsc, iat, iat, rij, weight, dwdr, dwdL)
      do img = 1, wsc%nimg(iat, iat)
         vec(:) = wsc%trans(:, wsc%tridx(img, iat, iat))
         rr = rad(iat)

         call get_damat_sdq_rec_3d(vec, qat(iat), qat(iat), dpat(:, iat), dpat(:, iat), &
            & qpat(:, iat), qpat(:, iat), ewald, dGr, dSr, erec)
         call get_damat_sdq_dir_3d(vec, qat(iat), qat(iat), dpat(:, iat), dpat(:, iat), &
            & qpat(:, iat), qpat(:, iat), rr, kdmp3, kdmp5, alpha, dtrans, dEd, dGd, dSd, edir)
         eimg = erec + edir
         dE = dE + dEd * weight(img)
         dS = dS + (dSd + dSr) * weight(img) + eimg*dwdL(:, :, img)
      end do
      dEdr(iat) = dEdr(iat) + dE
      sigma = sigma + 0.5_wp * dS
   end do
end subroutine get_multipole_gradient_3d


pure subroutine get_damat_sdq_rec_3d(rij, qi, qj, mi, mj, ti, tj, ewald, dg, ds, energy)
   real(wp), intent(in) :: rij(3), qi, qj, mi(3), mj(3), ti(6), tj(6)
   type(ewald_cache), intent(in) :: ewald
   real(wp), intent(out) :: dg(3), ds(3, 3), energy

   integer :: itr, b
   real(wp) :: vec(3), phase, sink, cosk, dip(3), quad(6), qv(3), cross(3)
   real(wp) :: dpiv, dpjv, dv, qvv, value, weight

   dg = 0.0_wp
   ds = 0.0_wp
   energy = 0.0_wp
   dip = qj*mi - qi*mj
   quad = qi*tj + qj*ti
   do itr = 1, size(ewald%weight)
      vec = ewald%vec(:, itr)
      weight = ewald%weight(itr)
      phase = dot_product(rij, vec)
      sink = sin(phase)
      cosk = cos(phase)
      dv = dot_product(dip, vec)
      dpiv = dot_product(mi, vec)
      dpjv = dot_product(mj, vec)
      qv = [quad(1)*vec(1) + quad(2)*vec(2) + quad(4)*vec(3), &
         & quad(2)*vec(1) + quad(3)*vec(2) + quad(5)*vec(3), &
         & quad(4)*vec(1) + quad(5)*vec(2) + quad(6)*vec(3)]
      qvv = dot_product(vec, qv)/3.0_wp
      value = sink*dv + cosk*(dpiv*dpjv - qvv)
      energy = energy + weight*value
      dg = dg + weight*vec*(-cosk*dv + sink*(dpiv*dpjv - qvv))
      cross = weight*(sink*dip + cosk*(mi*dpjv + mj*dpiv - 2.0_wp/3.0_wp*qv))
      do b = 1, 3
         ds(:, b) = ds(:, b) + value*ewald%strain(:, b, itr) &
            & - 0.5_wp*(vec*cross(b) + cross*vec(b))
      end do
   end do
end subroutine get_damat_sdq_rec_3d

pure subroutine get_damat_sdq_dir_3d(rij, qi, qj, mi, mj, ti, tj, rr, kdmp3, kdmp5, &
      & alp, trans, de, dg, ds, energy)
   real(wp), intent(in) :: rij(3), qi, qj, mi(3), mj(3), ti(6), tj(6)
   real(wp), intent(in) :: rr, kdmp3, kdmp5, alp, trans(:, :)
   real(wp), intent(out) :: de, dg(3), ds(3, 3), energy

   integer :: itr, b
   real(wp) :: vec(3), r1, r2, g1, g3, g5, g7, fdmp3, fdmp5, ddmp3, ddmp5
   real(wp) :: alp2, arg, erft, expt, e1, e2, e3, damp3, damp5
   real(wp) :: dip(3), quad(6), trace, qv(3), dv, dpidpj, dpiv, dpjv, edd, eq
   real(wp) :: g_sd(3), g_dd(3), g_sq(3), force(3)

   de = 0.0_wp
   dg = 0.0_wp
   ds = 0.0_wp
   energy = 0.0_wp
   alp2 = alp*alp
   dip = qj*mi - qi*mj
   quad = qi*tj + qj*ti
   trace = quad(1) + quad(3) + quad(6)
   dpidpj = dot_product(mi, mj)

   do itr = 1, size(trans, 2)
      vec = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      r2 = r1*r1
      g1 = 1.0_wp/r1
      g3 = g1*g1*g1
      g5 = g3*g1*g1
      g7 = g5*g1*g1
      arg = r1*alp
      erft = -erf(arg)*g1
      expt = exp(-arg*arg)/sqrtpi
      e1 = g1*g1*(erft + expt*(2*alp2)/alp)
      e2 = g1*g1*(e1 + expt*(2*alp2)**2/(3*alp))
      e3 = g1*g1*(e2 + expt*(2*alp2)**3/(15*alp))

      damp3 = (rr*g1)**kdmp3
      damp5 = (rr*g1)**kdmp5
      fdmp3 = 1.0_wp/(1.0_wp + 6.0_wp*damp3)
      ddmp3 = -3*fdmp3 - kdmp3*fdmp3*(fdmp3-1.0_wp)
      fdmp5 = 1.0_wp/(1.0_wp + 6.0_wp*damp5)
      ddmp5 = -5*fdmp5 - kdmp5*(fdmp5*fdmp5-fdmp5)

      dv = dot_product(vec, dip)
      dpiv = dot_product(mi, vec)
      dpjv = dot_product(mj, vec)
      edd = dpidpj*r2 - 3*dpjv*dpiv
      qv = [quad(1)*vec(1) + quad(2)*vec(2) + quad(4)*vec(3), &
         & quad(2)*vec(1) + quad(3)*vec(2) + quad(5)*vec(3), &
         & quad(4)*vec(1) + quad(5)*vec(2) + quad(6)*vec(3)]
      eq = dot_product(vec, qv)

      ! Analytic contraction of the second and third derivatives of 1/r.
      g_sd = (-ddmp3*g5 + 3*e2)*dv*vec - (fdmp3*g3 + e1)*dip
      g_dd = (-2*fdmp5*g5*dpidpj - edd*ddmp5*g7 - 15*e3*dpiv*dpjv)*vec &
         & + 3*(fdmp5*g5 + e2)*(dpiv*mj + dpjv*mi) + 3*e2*dpidpj*vec
      g_sq = (-eq*ddmp5*g7 + 5*e3*eq - e2*trace)*vec &
         & - 2*(fdmp5*g5 + e2)*qv
      force = g_sd + g_dd + g_sq
      dg = dg + force
      do b = 1, 3
         ds(:, b) = ds(:, b) - 0.5_wp*(vec*force(b) + force*vec(b))
      end do
      de = de + 3.0_wp*dv*kdmp3*fdmp3*g3*(fdmp3/rr)*damp3 &
         & + 3.0_wp*edd*kdmp5*fdmp5*g5*(fdmp5/rr)*damp5 &
         & + eq*3.0_wp*kdmp5*fdmp5*g5*fdmp5/rr*damp5
      energy = energy + dv*(fdmp3*g3 + e1) &
         & + (dpidpj - trace/3.0_wp)*(fdmp5*g3 + e1) &
         & + (eq - 3*dpiv*dpjv)*(fdmp5*g5 + e2)
   end do
end subroutine get_damat_sdq_dir_3d


!> Get information about density dependent quantities used in the energy
pure function variable_info(self) result(info)
   use tblite_scf_info, only : scf_info, atom_resolved
   !> Instance of the electrostatic container
   class(damped_multipole), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(dipole=atom_resolved, quadrupole=atom_resolved)
end function variable_info


!> Inspect container cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(coulomb_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(coulomb_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

!> Return reference to container cache after resolving its type
subroutine view(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(coulomb_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(coulomb_cache)
      ptr => target
   end select
end subroutine view

subroutine get_energy_axc(self, mol, wfn, energies)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Electrostatic energy
   real(wp), intent(inout) :: energies(:)

   call get_kernel_energy(mol, self%dkernel, wfn%dpat(:, :, 1), energies)
   call get_kernel_energy(mol, self%qkernel, wfn%qpat(:, :, 1), energies)

end subroutine get_energy_axc

end module tblite_coulomb_multipole
