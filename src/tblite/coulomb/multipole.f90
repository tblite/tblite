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
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   use mctc_io_constants, only : pi
   use tblite_blas, only : dot, gemv, symv, gemm
   use tblite_container_cache, only : container_cache
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_ewald, only : get_dir_cutoff, get_rec_cutoff
   use tblite_coulomb_type, only : coulomb_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_ncoord_gfn, only : gfn_ncoord_type, new_gfn_ncoord
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_wignerseitz, only : wignerseitz_cell
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
      type(gfn_ncoord_type), allocatable :: ncoord
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
   real(wp), parameter :: twopi = 2 * pi
   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))
   real(wp), parameter :: conv = 100*eps
   character(len=*), parameter :: label = "anisotropic electrostatics"

contains


!> Create a new anisotropic electrostatics container
subroutine new_damped_multipole(self, mol, kdmp3, kdmp5, dkernel, qkernel, &
      & shift, kexp, rmax, rad, vcn)
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

   allocate(self%ncoord)
   call new_gfn_ncoord(self%ncoord, mol)
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

   call taint(cache, ptr)

   if (.not.allocated(ptr%mrad)) then
      allocate(ptr%mrad(mol%nat))
   end if
   if (.not.allocated(ptr%dmrdcn)) then
      allocate(ptr%dmrdcn(mol%nat))
   end if

   if (.not.allocated(ptr%amat_sd)) then
      allocate(ptr%amat_sd(3, mol%nat, mol%nat))
   end if
   if (.not.allocated(ptr%amat_dd)) then
      allocate(ptr%amat_dd(3, mol%nat, 3, mol%nat))
   end if
   if (.not.allocated(ptr%amat_sq)) then
      allocate(ptr%amat_sq(6, mol%nat, mol%nat))
   end if

   if (.not.allocated(ptr%cn)) then
      allocate(ptr%cn(mol%nat))
   end if
   if (.not.allocated(ptr%dcndr)) then
      allocate(ptr%dcndr(3, mol%nat, mol%nat))
   end if
   if (.not.allocated(ptr%dcndL)) then
      allocate(ptr%dcndL(3, 3, mol%nat))
   end if

   if (allocated(self%ncoord)) then
      call self%ncoord%get_cn(mol, ptr%cn, ptr%dcndr, ptr%dcndL)
   else
      ptr%cn(:) = self%valence_cn(mol%id)
      ptr%dcndr(:, :, :) = 0.0_wp
      ptr%dcndL(:, :, :) = 0.0_wp
   end if

   call get_mrad(mol, self%shift, self%kexp, self%rmax, self%rad, self%valence_cn, &
      & ptr%cn, ptr%mrad, ptr%dmrdcn)

   call get_multipole_matrix(self, mol, ptr, ptr%amat_sd, ptr%amat_dd, ptr%amat_sq)
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

   real(wp), allocatable :: vs(:), vd(:, :), vq(:, :)
   type(coulomb_cache), pointer :: ptr

   call view(cache, ptr)

   allocate(vs(mol%nat), vd(3, mol%nat), vq(6, mol%nat))

   call gemv(ptr%amat_sd, wfn%qat(:, 1), vd)
   call gemv(ptr%amat_dd, wfn%dpat(:, :, 1), vd, beta=1.0_wp, alpha=0.5_wp)
   call gemv(ptr%amat_sq, wfn%qat(:, 1), vq)

   energies(:) = energies + sum(wfn%dpat(:, :, 1) * vd, 1) + sum(wfn%qpat(:, :, 1) * vq, 1)

   call get_kernel_energy(mol, self%dkernel, wfn%dpat(:, :, 1), energies)
   call get_kernel_energy(mol, self%qkernel, wfn%qpat(:, :, 1), energies)
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

   real(wp), allocatable :: vs(:), vd(:, :), vq(:, :),mur(:)
   real(wp), allocatable :: e01(:),e11(:),e02(:)
   real(wp), allocatable :: t1(:)
   type(coulomb_cache), pointer :: ptr
   integer :: i,j

   call view(cache, ptr)

   allocate(vs(mol%nat), vd(3, mol%nat), vq(6, mol%nat))

   allocate(mur(mol%nat), source=0.0_wp)
   allocate(e01(mol%nat),e11(mol%nat),e02(mol%nat),source=0.0_wp)
   allocate(t1(mol%nat),source=0.0_wp)

   call gemv(ptr%amat_sd, wfn%qat(:, 1), vd)
   do i = 1,3
      call gemv(ptr%amat_sd(i,:,:),wfn%dpat(i,:,1),mur,beta=1.0_wp, alpha=1.0_wp,trans="T")
   end do

   e01 = (mur)*wfn%qat(:,1) + sum(wfn%dpat(:, :, 1) * vd, 1)

   vd = 0.0_wp
   call gemv(ptr%amat_dd, wfn%dpat(:, :, 1), vd, beta=1.0_wp, alpha=0.5_wp)
   e11 = sum(wfn%dpat(:, :, 1) * vd, 1)

   call gemv(ptr%amat_sq, wfn%qat(:, 1), vq)

   do i = 1,6
      call gemv(ptr%amat_sq(i,:,:),wfn%qpat(i,:,1),t1,beta=1.0_wp, alpha=1.0_wp,trans="T")
   end do
   e02 = t1*wfn%qat(:,1) + sum(wfn%qpat(:, :, 1) * vq, 1)

   energies(:) = energies + 0.5_wp * e01 + e11 + 0.5_wp * e02

end subroutine get_energy_aes

!> Get multipolar anisotropic exchange-correlation kernel
subroutine get_kernel_energy(mol, kernel, mpat, energies)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multipole kernel
   real(wp), intent(in) :: kernel(:)
   !> Atomic multipole momemnt
   real(wp), intent(in) :: mpat(:, :)
   !> Electrostatic energy
   real(wp), intent(inout) :: energies(:)

   integer :: iat, izp
   real(wp) :: mpt(size(mpat, 1)), mpscale(size(mpat, 1))

   mpscale(:) = 1
   if (size(mpat, 1) == 6) mpscale([2, 4, 5]) = 2

   do iat = 1, mol%nat
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

   call gemv(ptr%amat_sd, wfn%qat(:, 1), pot%vdp(:, :, 1), beta=1.0_wp)
   call gemv(ptr%amat_sd, wfn%dpat(:, :, 1), pot%vat(:, 1), beta=1.0_wp, trans="T")

   call gemv(ptr%amat_dd, wfn%dpat(:, :, 1), pot%vdp(:, :, 1), beta=1.0_wp)

   call gemv(ptr%amat_sq, wfn%qat(:, 1), pot%vqp(:, :, 1), beta=1.0_wp)
   call gemv(ptr%amat_sq, wfn%qpat(:, :, 1), pot%vat(:, 1), beta=1.0_wp, trans="T")

   call get_kernel_potential(mol, self%dkernel, wfn%dpat(:, :, 1), pot%vdp(:, :, 1))
   call get_kernel_potential(mol, self%qkernel, wfn%qpat(:, :, 1), pot%vqp(:, :, 1))
end subroutine get_potential


!> Get multipolar anisotropic potential contribution
subroutine get_kernel_potential(mol, kernel, mpat, vm)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multipole kernel
   real(wp), intent(in) :: kernel(:)
   !> Atomic multipole momemnt
   real(wp), intent(in) :: mpat(:, :)
   !> Potential shoft on atomic multipole moment
   real(wp), intent(inout) :: vm(:, :)

   integer :: iat, izp
   real(wp) :: mpscale(size(mpat, 1))

   mpscale(:) = 1
   if (size(mpat, 1) == 6) mpscale([2, 4, 5]) = 2

   do iat = 1, mol%nat
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

   real(wp), allocatable :: dEdr(:)
   type(coulomb_cache), pointer :: ptr

   call view(cache, ptr)

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
subroutine get_dir_trans(lattice, alpha, conv, trans)
   !> Lattice parameters
   real(wp), intent(in) :: lattice(:, :)
   !> Parameter for Ewald summation
   real(wp), intent(in) :: alpha
   !> Tolerance for Ewald summation
   real(wp), intent(in) :: conv
   !> Translation vectors
   real(wp), allocatable, intent(out) :: trans(:, :)

   call get_lattice_points([.true.], lattice, get_dir_cutoff(alpha, conv), trans)

end subroutine get_dir_trans

!> Get reciprocal lattice translations
subroutine get_rec_trans(lattice, alpha, volume, conv, trans)
   !> Lattice parameters
   real(wp), intent(in) :: lattice(:, :)
   !> Parameter for Ewald summation
   real(wp), intent(in) :: alpha
   !> Cell volume
   real(wp), intent(in) :: volume
   !> Tolerance for Ewald summation
   real(wp), intent(in) :: conv
   !> Translation vectors
   real(wp), allocatable, intent(out) :: trans(:, :)

   real(wp) :: rec_lat(3, 3)

   rec_lat = twopi*transpose(matinv_3x3(lattice))
   call get_lattice_points([.true.], rec_lat, get_rec_cutoff(alpha, volume, conv), trans)
   trans = trans(:, 2:)

end subroutine get_rec_trans


!> Get interaction matrix for all multipole moments up to inverse cubic order
subroutine get_multipole_matrix(self, mol, cache, amat_sd, amat_dd, amat_sq)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   !> Interation matrix for charges and dipoles
   real(wp), intent(inout) :: amat_sd(:, :, :)
   !> Interation matrix for dipoles and dipoles
   real(wp), intent(inout) :: amat_dd(:, :, :, :)
   !> Interation matrix for charges and quadrupoles
   real(wp), intent(inout) :: amat_sq(:, :, :)

   amat_sd(:, :, :) = 0.0_wp
   amat_dd(:, :, :, :) = 0.0_wp
   amat_sq(:, :, :) = 0.0_wp
   if (any(mol%periodic)) then
      call get_multipole_matrix_3d(mol, cache%mrad, self%kdmp3, self%kdmp5, &
         & cache%wsc, cache%alpha, amat_sd, amat_dd, amat_sq)
   else
      call get_multipole_matrix_0d(mol, cache%mrad, self%kdmp3, self%kdmp5, &
         & amat_sd, amat_dd, amat_sq)
   end if
end subroutine get_multipole_matrix

!> Calculate the multipole interaction matrix for finite systems
subroutine get_multipole_matrix_0d(mol, rad, kdmp3, kdmp5, amat_sd, amat_dd, amat_sq)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multipole damping radii for all atoms
   real(wp), intent(in) :: rad(:)
   !> Damping function for inverse quadratic contributions
   real(wp), intent(in) :: kdmp3
   !> Damping function for inverse cubic contributions
   real(wp), intent(in) :: kdmp5
   !> Interaction matrix for charges and dipoles
   real(wp), intent(inout) :: amat_sd(:, :, :)
   !> Interaction matrix for dipoles and dipoles
   real(wp), intent(inout) :: amat_dd(:, :, :, :)
   !> Interaction matrix for charges and quadrupoles
   real(wp), intent(inout) :: amat_sq(:, :, :)

   integer :: iat, jat
   real(wp) :: r1, vec(3), g1, g3, g5, fdmp3, fdmp5, tc(6), rr

   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp shared(amat_sd, amat_dd, amat_sq, mol, rad, kdmp3, kdmp5) &
   !$omp private(r1, vec, g1, g3, g5, fdmp3, fdmp5, tc, rr)
   do iat = 1, mol%nat
      do jat = 1, mol%nat
         if (iat == jat) cycle
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         r1 = norm2(vec)
         g1 = 1.0_wp / r1
         g3 = g1 * g1 * g1
         g5 = g3 * g1 * g1

         rr = 0.5_wp * (rad(jat) + rad(iat)) * g1
         fdmp3 = 1.0_wp / (1.0_wp + 6.0_wp * rr**kdmp3)
         fdmp5 = 1.0_wp / (1.0_wp + 6.0_wp * rr**kdmp5)

         amat_sd(:, jat, iat) = amat_sd(:, jat, iat) + vec * g3 * fdmp3
         amat_dd(:, jat, :, iat) = amat_dd(:, jat, :, iat) &
            & + unity * g3*fdmp5 - spread(vec, 1, 3) * spread(vec, 2, 3) * 3*g5*fdmp5
         tc(2) = 2*vec(1)*vec(2)*g5*fdmp5
         tc(4) = 2*vec(1)*vec(3)*g5*fdmp5
         tc(5) = 2*vec(2)*vec(3)*g5*fdmp5
         tc(1) = vec(1)*vec(1)*g5*fdmp5
         tc(3) = vec(2)*vec(2)*g5*fdmp5
         tc(6) = vec(3)*vec(3)*g5*fdmp5
         amat_sq(:, jat, iat) = amat_sq(:, jat, iat) + tc
      end do
   end do
end subroutine get_multipole_matrix_0d

!> Evaluate multipole interaction matrix under 3D periodic boundary conditions
subroutine get_multipole_matrix_3d(mol, rad, kdmp3, kdmp5, wsc, alpha, &
      & amat_sd, amat_dd, amat_sq)
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
   !> Interation matrix for charges and dipoles
   real(wp), intent(inout) :: amat_sd(:, :, :)
   !> Interation matrix for dipoles and dipoles
   real(wp), intent(inout) :: amat_dd(:, :, :, :)
   !> Interation matrix for charges and quadrupoles
   real(wp), intent(inout) :: amat_sq(:, :, :)

   integer :: iat, jat, img, k
   real(wp) :: vec(3), rr, wsw, vol
   real(wp) :: d_sd(3), d_dd(3, 3), d_sq(6), r_sd(3), r_dd(3, 3), r_sq(6)
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, alpha, conv, dtrans)
   call get_rec_trans(mol%lattice, alpha, vol, conv, rtrans)

   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp shared(amat_sd, amat_dd, amat_sq) &
   !$omp shared(mol, wsc, rad, vol, alpha, rtrans, dtrans, kdmp3, kdmp5) &
   !$omp private(iat, jat, img, vec, rr, wsw, d_sd, d_dd, d_sq, r_sd, r_dd, r_sq)
   do iat = 1, mol%nat
      do jat = 1, mol%nat
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))

            rr = 0.5_wp * (rad(jat) + rad(iat))
            call get_amat_sdq_rec_3d(vec, vol, alpha, rtrans, r_sd, r_dd, r_sq)
            call get_amat_sdq_dir_3d(vec, rr, kdmp3, kdmp5, alpha, dtrans, d_sd, d_dd, d_sq)

            amat_sd(:, jat, iat) = amat_sd(:, jat, iat) + wsw * (d_sd + r_sd)
            amat_dd(:, jat, :, iat) = amat_dd(:, jat, :, iat) + wsw * (r_dd + d_dd)
            amat_sq(:, jat, iat) = amat_sq(:, jat, iat) + wsw * (r_sq + d_sq)
         end do
      end do
   end do

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(amat_sd, amat_dd, amat_sq, mol, vol, alpha) private(iat, rr, k)
   do iat = 1, mol%nat
      ! dipole-dipole selfenergy: -2/3·α³/sqrt(π) Σ(i) μ²(i)
      rr = -2.0_wp/3.0_wp * alpha**3 / sqrtpi
      do k = 1, 3
         amat_dd(k, iat, k, iat) = amat_dd(k, iat, k, iat) + 2*rr
      end do

      ! charge-quadrupole selfenergy: 4/9·α³/sqrt(π) Σ(i) q(i)Tr(θi)
      ! (no actual contribution since quadrupoles are traceless)
      rr = 4.0_wp/9.0_wp * alpha**3 / sqrtpi
      amat_sq([1, 3, 6], iat, iat) = amat_sq([1, 3, 6], iat, iat) + rr
   end do
end subroutine get_multipole_matrix_3d

pure subroutine get_amat_sdq_rec_3d(rij, vol, alp, trans, amat_sd, amat_dd, amat_sq)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: vol
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: amat_sd(:)
   real(wp), intent(out) :: amat_dd(:, :)
   real(wp), intent(out) :: amat_sq(:)

   integer :: itr
   real(wp) :: fac, vec(3), g2, expk, sink, cosk, gv

   amat_sd = 0.0_wp
   amat_dd = 0.0_wp
   amat_sq = 0.0_wp
   fac = 4*pi/vol

   amat_dd(1, 1) = fac/6.0_wp
   amat_dd(2, 2) = fac/6.0_wp
   amat_dd(3, 3) = fac/6.0_wp
   amat_sq([1, 3, 6]) = -fac/9.0_wp

   do itr = 1, size(trans, 2)
      vec(:) = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      gv = dot_product(rij, vec)
      expk = fac*exp(-0.25_wp*g2/(alp*alp))/g2
      sink = sin(gv)*expk
      cosk = cos(gv)*expk

      amat_sd(:) = amat_sd + 2*vec*sink
      amat_dd(:, :) = amat_dd + spread(vec, 1, 3) * spread(vec, 2, 3) * cosk
      amat_sq(1) = amat_sq(1) +   vec(1)*vec(1)*cosk
      amat_sq(2) = amat_sq(2) + 2*vec(1)*vec(2)*cosk
      amat_sq(3) = amat_sq(3) +   vec(2)*vec(2)*cosk
      amat_sq(4) = amat_sq(4) + 2*vec(1)*vec(3)*cosk
      amat_sq(5) = amat_sq(5) + 2*vec(2)*vec(3)*cosk
      amat_sq(6) = amat_sq(6) +   vec(3)*vec(3)*cosk
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
      amat_dd(:, :) = amat_dd(:, :) + unity * (fdmp5*g3 + e1) &
         & - spread(vec, 1, 3) * spread(vec, 2, 3) * (3 * (g5*fdmp5 + e2))
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
      call get_multipole_gradient_3d(mol, cache%mrad, self%kdmp3, self%kdmp5, &
         & qat, dpat, qpat, cache%wsc, cache%alpha, dEdr, gradient, sigma)
   else
      call get_multipole_gradient_0d(mol, cache%mrad, self%kdmp3, self%kdmp5, &
         & qat, dpat, qpat, dEdr, gradient, sigma)
   end if
end subroutine get_multipole_gradient

!> Evaluate multipole derivatives for finite systems
subroutine get_multipole_gradient_0d(mol, rad, kdmp3, kdmp5, qat, dpat, qpat, &
      & dEdr, gradient, sigma)
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

   integer :: iat, jat
   real(wp) :: r1, r2, vec(3), rr, fdmp3, fdmp5, g1, g3, g5, g7, dG(3), dS(3, 3)
   real(wp) :: ddmp3, ddmp5, fddr, eq, edd, dpidpj, dpiv, dpjv, dpiqj, qidpj

   !$omp parallel do default(none) schedule(runtime) reduction(+:dEdr, gradient, sigma) &
   !$omp shared(mol, kdmp3, kdmp5, rad, qat, dpat, qpat) &
   !$omp private(iat, jat, r1, r2, vec, rr, fdmp3, fdmp5, g1, g3, g5, g7, dG, dS, &
   !$omp& ddmp3, ddmp5, fddr, eq, edd, dpidpj, dpiv, dpjv, dpiqj, qidpj)
   do iat = 1, mol%nat
      do jat = 1, iat - 1
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
subroutine get_multipole_gradient_3d(mol, rad, kdmp3, kdmp5, qat, dpat, qpat, wsc, alpha, &
      & dEdr, gradient, sigma)
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
   !> Derivative of the energy w.r.t. the critical radii
   real(wp), contiguous, intent(inout) :: dEdr(:)
   !> Derivative of the energy w.r.t. atomic displacements
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Derivative of the energy w.r.t. strain deformations
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   integer :: iat, jat, img
   real(wp) :: vec(3), dG(3), dGr(3), dGd(3), dS(3, 3), dSr(3, 3), dSd(3, 3)
   real(wp) :: wsw, vol, rr, dE, dEd
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, alpha, conv, dtrans)
   call get_rec_trans(mol%lattice, alpha, vol, conv, rtrans)

   !$omp parallel do default(none) schedule(runtime) reduction(+:dEdr, gradient, sigma) &
   !$omp shared(mol, wsc, kdmp3, kdmp5, vol, alpha, rtrans, dtrans, rad, qat, dpat, qpat) &
   !$omp private(iat, jat, dE, dG, dS, wsw, img, vec, rr, dEd, dGd, dGr, dSd, dSr)
   do iat = 1, mol%nat
      do jat = 1, iat - 1
         dE = 0.0_wp
         dG(:) = 0.0_wp
         dS(:, :) = 0.0_wp
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec(:) = mol%xyz(:, jat)-mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, iat))
            rr = 0.5_wp * (rad(jat) + rad(iat))

            call get_damat_sdq_rec_3d(vec, qat(iat), qat(jat), dpat(:, iat), dpat(:, jat), &
               & qpat(:, iat), qpat(:, jat), vol, alpha, rtrans, dGr, dSr)
            call get_damat_sdq_dir_3d(vec, qat(iat), qat(jat), dpat(:, iat), dpat(:, jat), &
               & qpat(:, iat), qpat(:, jat), rr, kdmp3, kdmp5, alpha, dtrans, dEd, dGd, dSd)
            dE = dE + dEd * wsw
            dG = dG + (dGd + dGr) * wsw
            dS = dS + (dSd + dSr) * wsw
         end do
         dEdr(iat) = dEdr(iat) + dE
         dEdr(jat) = dEdr(jat) + dE
         gradient(:, iat) = gradient(:, iat) + dG
         gradient(:, jat) = gradient(:, jat) - dG
         sigma = sigma + dS
      end do
   end do

   !$omp parallel do default(none) schedule(runtime) reduction(+:dEdr, sigma) &
   !$omp shared(mol, wsc, kdmp3, kdmp5, vol, alpha, rtrans, dtrans, rad, qat, dpat, qpat) &
   !$omp private(iat, jat, dE, dG, dS, wsw, img, vec, rr, dEd, dGd, dGr, dSd, dSr)
   do iat = 1, mol%nat
      dE = 0.0_wp
      dS(:, :) = 0.0_wp
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec(:) = wsc%trans(:, wsc%tridx(img, iat, iat))
         rr = rad(iat)

         call get_damat_sdq_rec_3d(vec, qat(iat), qat(iat), dpat(:, iat), dpat(:, iat), &
            & qpat(:, iat), qpat(:, iat), vol, alpha, rtrans, dGr, dSr)
         call get_damat_sdq_dir_3d(vec, qat(iat), qat(iat), dpat(:, iat), dpat(:, iat), &
            & qpat(:, iat), qpat(:, iat), rr, kdmp3, kdmp5, alpha, dtrans, dEd, dGd, dSd)
         dE = dE + dEd * wsw
         dS = dS + (dSd + dSr) * wsw
      end do
      dEdr(iat) = dEdr(iat) + dE
      sigma = sigma + 0.5_wp * dS
   end do
end subroutine get_multipole_gradient_3d

pure subroutine get_damat_sdq_rec_3d(rij, qi, qj, mi, mj, ti, tj, vol, alp, trans, dg, ds)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: qi, qj, mi(3), mj(3), ti(6), tj(6)
   real(wp), intent(in) :: vol
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: dg(3)
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: fac, vec(3), g2, gv, expk, alp2, sink, cosk, dpiqj, qidpj
   real(wp) :: qpiqj, qiqpj, dpiv, dpjv

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp
   fac = 4*pi/vol
   alp2 = alp*alp

   do itr = 1, size(trans, 2)
      vec(:) = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      gv = dot_product(rij, vec)
      expk = fac*exp(-0.25_wp*g2/alp2)/g2
      cosk = cos(gv)*expk
      sink = sin(gv)*expk

      dpiqj = dot_product(vec, mi)*qj
      qidpj = dot_product(vec, mj)*qi

      dg(:) = dg - 2*vec*cosk * (dpiqj - qidpj)
      ds(:, :) = ds + 2 * sink * (dpiqj - qidpj) &
         & * ((2.0_wp/g2 + 0.5_wp/alp2) * spread(vec, 1, 3)*spread(vec, 2, 3) - unity)

      dpiv = dot_product(mi, vec)
      dpjv = dot_product(mj, vec)

      dg(:) = dg + vec*sink*dpiv*dpjv
      ds(:, :) = ds + cosk * dpiv*dpjv &
         & * ((2.0_wp/g2 + 0.5_wp/alp2) * spread(vec, 1, 3)*spread(vec, 2, 3) - unity)

      qiqpj = qi*(tj(1)*vec(1)*vec(1) + tj(3)*vec(2)*vec(2) + tj(6)*vec(3)*vec(3) &
         & + 2*tj(2)*vec(1)*vec(2) + 2*tj(4)*vec(1)*vec(3) + 2*tj(5)*vec(2)*vec(3))
      qpiqj = qj*(ti(1)*vec(1)*vec(1) + ti(3)*vec(2)*vec(2) + ti(6)*vec(3)*vec(3) &
         & + 2*ti(2)*vec(1)*vec(2) + 2*ti(4)*vec(1)*vec(3) + 2*ti(5)*vec(2)*vec(3))

      dg(:) = dg + vec * sink * (qiqpj + qpiqj)
      ds(:, :) = ds + cosk * (qiqpj + qpiqj) &
         & * ((2.0_wp/g2 + 0.5_wp/alp2) * spread(vec, 1, 3)*spread(vec, 2, 3) - unity)
   end do

end subroutine get_damat_sdq_rec_3d

pure subroutine get_damat_sdq_dir_3d(rij, qi, qj, mi, mj, ti, tj, rr, kdmp3, kdmp5, &
      & alp, trans, de, dg, ds)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: qi, qj, mi(3), mj(3), ti(6), tj(6)
   real(wp), intent(in) :: rr
   real(wp), intent(in) :: kdmp3
   real(wp), intent(in) :: kdmp5
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: de
   real(wp), intent(out) :: dg(3)
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr, k
   real(wp) :: vec(3), r1, r2, g1, g3, g5, g7, fdmp3, fdmp5, ddmp3, ddmp5, dpiqj, qidpj
   real(wp) :: alp2, arg, arg2, erft, expt, e1, e2, e3, tabc(3, 3, 3), tab(3, 3)
   real(wp) :: dpidpj, dpiv, dpjv, edd, eq, g_sd(3), g_dd(3), g_sq(3)

   de = 0.0_wp
   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   alp2 = alp*alp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      r2 = r1 * r1
      g1 = 1.0_wp/r1
      g3 = g1 * g1 * g1
      g5 = g3 * g1 * g1
      g7 = g5 * g1 * g1

      arg = r1*alp
      arg2 = arg*arg
      erft = -erf(arg)*g1
      expt = exp(-arg2)/sqrtpi
      e1 = g1*g1 * (erft + expt*(2*alp2)/alp)
      e2 = g1*g1 * (e1 + expt*(2*alp2)**2/(3*alp))
      e3 = g1*g1 * (e2 + expt*(2*alp2)**3/(15*alp))

      fdmp3 = 1.0_wp/(1.0_wp+6.0_wp*(rr*g1)**kdmp3)
      ddmp3 = -3*fdmp3 - kdmp3*fdmp3*(fdmp3-1.0_wp)
      fdmp5 = 1.0_wp/(1.0_wp+6.0_wp*(rr*g1)**kdmp5)
      ddmp5 = -5*fdmp5 - kdmp5*(fdmp5*fdmp5-fdmp5)

      dpiqj = dot_product(vec, mi)*qj
      qidpj = dot_product(vec, mj)*qi
      dpidpj = dot_product(mj, mi)
      dpiv = dot_product(mi, vec)
      dpjv = dot_product(mj, vec)
      edd = dpidpj*r2 - 3*dpjv*dpiv

      ! Charge - dipole
      g_sd(:) = - (ddmp3*g5)*vec * (dpiqj - qidpj) + fdmp3*g3*(qi*mj - qj*mi)

      block
         integer :: a, b
         do b = 1, 3
            do a = 1, 3
               tab(a, b) = 3*vec(a)*vec(b)*e2
            end do
            tab(b, b) = tab(b, b) - e1
         end do
      end block

      do k = 1, 3
         g_sd(k) = g_sd(k) &
            & + qj * (tab(1, k) * mi(1) &
            &       + tab(2, k) * mi(2) &
            &       + tab(3, k) * mi(3))&
            & - qi * (tab(1, k) * mj(1) &
            &       + tab(2, k) * mj(2) &
            &       + tab(3, k) * mj(3))
      end do

      de = de + 3.0_wp*(dpiqj - qidpj)*kdmp3*fdmp3*g3*(fdmp3/rr)*(rr*g1)**kdmp3
      dg(:) = dg + g_sd
      ds(:, :) = - 0.5_wp * (spread(vec, 1, 3) * spread(g_sd, 2, 3) &
         & + spread(g_sd, 1, 3) * spread(vec, 2, 3))

      ! Dipole - dipole
      g_dd(:) = - 2.0_wp*fdmp5*g5*dpidpj*vec &
         & + 3.0_wp*fdmp5*g5*(dpiv*mj + dpjv*mi) &
         & - edd*ddmp5*g7*vec

      block
         integer :: a, b, c
         do c = 1, 3
            do b = 1, 3
               do a = 1, 3
                  tabc(a, b, c) = - 15*vec(a)*vec(b)*vec(c)*e3
               end do
            end do
            do a = 1, 3
               tabc(a, a, c) = tabc(a, a, c) + 3*e2*vec(c)
               tabc(c, a, c) = tabc(c, a, c) + 3*e2*vec(a)
               tabc(a, c, c) = tabc(a, c, c) + 3*e2*vec(a)
            end do
         end do
      end block

      do k = 1, 3
         g_dd(k) = g_dd(k) &
            & + mi(1)*(tabc(1, 1, k) * mj(1) &
            &        + tabc(2, 1, k) * mj(2) &
            &        + tabc(3, 1, k) * mj(3))&
            & + mi(2)*(tabc(1, 2, k) * mj(1) &
            &        + tabc(2, 2, k) * mj(2) &
            &        + tabc(3, 2, k) * mj(3))&
            & + mi(3)*(tabc(1, 3, k) * mj(1) &
            &        + tabc(2, 3, k) * mj(2) &
            &        + tabc(3, 3, k) * mj(3))
      end do

      de = de + 3.0_wp*edd*kdmp5*fdmp5*g5*(fdmp5/rr)*(rr*g1)**kdmp5
      dg(:) = dg + g_dd
      ds(:, :) = - 0.5_wp * (spread(vec, 1, 3) * spread(g_dd, 2, 3) &
         & + spread(g_dd, 1, 3) * spread(vec, 2, 3))

      ! Charge - quadrupole
      eq = &
         & + 2*(qj*ti(2) + tj(2)*qi)*vec(1)*vec(2) &
         & + 2*(qj*ti(4) + tj(4)*qi)*vec(1)*vec(3) &
         & + 2*(qj*ti(5) + tj(5)*qi)*vec(2)*vec(3) &
         & + (qj*ti(1) + tj(1)*qi)*vec(1)*vec(1) &
         & + (qj*ti(3) + tj(3)*qi)*vec(2)*vec(2) &
         & + (qj*ti(6) + tj(6)*qi)*vec(3)*vec(3)

      g_sq(:) = - eq*ddmp5*g7*vec &
         & - 2.0_wp*fdmp5*g5*qi * &
         &[vec(1)*tj(1) + vec(2)*tj(2) + vec(3)*tj(4), &
         & vec(1)*tj(2) + vec(2)*tj(3) + vec(3)*tj(5), &
         & vec(1)*tj(4) + vec(2)*tj(5) + vec(3)*tj(6)] &
         & - 2.0_wp*fdmp5*g5*qj * &
         &[vec(1)*ti(1) + vec(2)*ti(2) + vec(3)*ti(4), &
         & vec(1)*ti(2) + vec(2)*ti(3) + vec(3)*ti(5), &
         & vec(1)*ti(4) + vec(2)*ti(5) + vec(3)*ti(6)]

      do k = 1, 3
         g_sq(k) = g_sq(k) + (&
            & - qi * (tabc(1, 1, k) * tj(1) &
            &     + 2*tabc(2, 1, k) * tj(2) &
            &     + 2*tabc(3, 1, k) * tj(4) &
            &     +   tabc(2, 2, k) * tj(3) &
            &     + 2*tabc(3, 2, k) * tj(5) &
            &     +   tabc(3, 3, k) * tj(6))&
            & - qj * (tabc(1, 1, k) * ti(1) &
            &     + 2*tabc(2, 1, k) * ti(2) &
            &     + 2*tabc(3, 1, k) * ti(4) &
            &     +   tabc(2, 2, k) * ti(3) &
            &     + 2*tabc(3, 2, k) * ti(5) &
            &     +   tabc(3, 3, k) * ti(6)))/3.0_wp
      end do

      de = de + eq * 3.0_wp*kdmp5*fdmp5*g5*fdmp5/rr*(rr*g1)**kdmp5
      dg(:) = dg + g_sq
      ds(:, :) = ds - 0.5_wp * (spread(vec, 1, 3) * spread(g_sq, 2, 3) &
         & + spread(g_sq, 1, 3) * spread(vec, 2, 3))
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
