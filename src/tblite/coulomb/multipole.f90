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

!> Anisotropic second-order electrostatics using a damped multipole expansion
module tblite_coulomb_multipole
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   use mctc_io_constants, only : pi
   use tblite_blas, only : dot, gemv, symv, gemm
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_type, only : coulomb_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_ncoord_gfn, only : gfn_ncoord_type, new_gfn_ncoord
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_wignerseitz, only : wignerseitz_cell
   implicit none
   private

   public :: damped_multipole, new_damped_multipole


   type, extends(coulomb_type) :: damped_multipole
      real(wp) :: kdmp3 = 0.0_wp
      real(wp) :: kdmp5 = 0.0_wp
      real(wp), allocatable :: dkernel(:)
      real(wp), allocatable :: qkernel(:)

      real(wp) :: shift = 0.0_wp
      real(wp) :: kexp = 0.0_wp
      real(wp) :: rmax = 0.0_wp
      real(wp), allocatable :: rad(:)
      real(wp), allocatable :: valence_cn(:)

      type(gfn_ncoord_type), allocatable :: ncoord
   contains
      procedure :: update
      procedure :: variable_info
      procedure :: get_energy
      procedure :: get_potential
      procedure :: get_gradient
   end type damped_multipole

   real(wp), parameter :: unity(3, 3) = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   real(wp), parameter :: twopi = 2 * pi
   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

contains


subroutine new_damped_multipole(self, mol, kdmp3, kdmp5, dkernel, qkernel, &
      & shift, kexp, rmax, rad, vcn)
   !> Instance of the multipole container
   type(damped_multipole), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: kdmp3
   real(wp), intent(in) :: kdmp5
   real(wp), intent(in) :: dkernel(:)
   real(wp), intent(in) :: qkernel(:)
   real(wp), intent(in) :: shift
   real(wp), intent(in) :: kexp
   real(wp), intent(in) :: rmax
   real(wp), intent(in) :: rad(:)
   real(wp), intent(in) :: vcn(:)

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


subroutine update(self, mol, cache)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache

   if (.not.allocated(cache%mrad)) then
      allocate(cache%mrad(mol%nat))
   end if
   if (.not.allocated(cache%dmrdcn)) then
      allocate(cache%dmrdcn(mol%nat))
   end if

   if (.not.allocated(cache%amat_sd)) then
      allocate(cache%amat_sd(3, mol%nat, mol%nat))
   end if
   if (.not.allocated(cache%amat_dd)) then
      allocate(cache%amat_dd(3, mol%nat, 3, mol%nat))
   end if
   if (.not.allocated(cache%amat_sq)) then
      allocate(cache%amat_sq(6, mol%nat, mol%nat))
   end if

   if (.not.allocated(cache%cn)) then
      allocate(cache%cn(mol%nat))
   end if
   if (.not.allocated(cache%dcndr)) then
      allocate(cache%dcndr(3, mol%nat, mol%nat))
   end if
   if (.not.allocated(cache%dcndL)) then
      allocate(cache%dcndL(3, 3, mol%nat))
   end if

   if (allocated(self%ncoord)) then
      call self%ncoord%get_cn(mol, cache%cn, cache%dcndr, cache%dcndL)
   else
      cache%cn(:) = self%valence_cn(mol%id)
      cache%dcndr(:, :, :) = 0.0_wp
      cache%dcndL(:, :, :) = 0.0_wp
   end if

   call get_mrad(mol, self%shift, self%kexp, self%rmax, self%rad, self%valence_cn, &
      & cache%cn, cache%mrad, cache%dmrdcn)

   cache%amat_sd(:, :, :) = 0.0_wp
   call get_charge_dipole_matrix(self, mol, cache, cache%amat_sd)
   cache%amat_dd(:, :, :, :) = 0.0_wp
   call get_dipole_dipole_matrix(self, mol, cache, cache%amat_dd)
   cache%amat_sq(:, :, :) = 0.0_wp
   call get_charge_quadrupole_matrix(self, mol, cache, cache%amat_sq)
end subroutine update

subroutine get_energy(self, mol, cache, wfn, energy)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Electrostatic energy
   real(wp), intent(inout) :: energy
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache

   real(wp), allocatable :: vs(:), vd(:, :), vq(:, :)
   real(wp) :: ees, exc

   allocate(vs(mol%nat), vd(3, mol%nat), vq(6, mol%nat))

   call gemv(cache%amat_sd, wfn%qat, vd)
   call gemv(cache%amat_dd, wfn%dpat, vd, beta=1.0_wp, alpha=0.5_wp)
   call gemv(cache%amat_sq, wfn%qat, vq)

   ees = dot(wfn%dpat, vd) + dot(wfn%qpat, vq)
   exc = 0.0_wp

   call get_kernel_energy(mol, self%dkernel, wfn%dpat, exc)
   call get_kernel_energy(mol, self%qkernel, wfn%qpat, exc)

   energy = energy + ees + exc
end subroutine get_energy

subroutine get_kernel_energy(mol, kernel, mpat, energy)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multipole kernel
   real(wp), intent(in) :: kernel(:)
   !> Atomic multipole momemnt
   real(wp), intent(in) :: mpat(:, :)
   !> Electrostatic energy
   real(wp), intent(inout) :: energy

   integer :: iat, izp
   real(wp) :: mpt(size(mpat, 1)), mpscale(size(mpat, 1))

   mpscale(:) = 1
   if (size(mpat, 1) == 6) mpscale([2, 4, 5]) = 2

   do iat = 1, mol%nat
      izp = mol%id(iat)
      mpt(:) = mpat(:, iat) * mpscale
      energy = energy + kernel(izp) * dot_product(mpt, mpat(:, iat))
   end do
end subroutine get_kernel_energy

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
   type(coulomb_cache), intent(inout) :: cache

   call gemv(cache%amat_sd, wfn%qat, pot%vdp, beta=1.0_wp)
   call gemv(cache%amat_sd, wfn%dpat, pot%vat, beta=1.0_wp, trans="T")

   call gemv(cache%amat_dd, wfn%dpat, pot%vdp, beta=1.0_wp)

   call gemv(cache%amat_sq, wfn%qat, pot%vqp, beta=1.0_wp)
   call gemv(cache%amat_sq, wfn%qpat, pot%vat, beta=1.0_wp, trans="T")

   call get_kernel_potential(mol, self%dkernel, wfn%dpat, pot%vdp)
   call get_kernel_potential(mol, self%qkernel, wfn%qpat, pot%vqp)
end subroutine get_potential

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

subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular gradient of the repulsion energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the repulsion energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   real(wp), allocatable :: dEdr(:)

   allocate(dEdr(mol%nat))
   dEdr = 0.0_wp

   call get_charge_dipole_gradient(self, mol, cache, wfn%qat, wfn%dpat, dEdr, &
      & gradient, sigma)
   call get_dipole_dipole_gradient(self, mol, cache, wfn%dpat, dEdr, gradient, sigma)
   call get_charge_quadrupole_gradient(self, mol, cache, wfn%qat, wfn%qpat, dEdr, &
      & gradient, sigma)

   dEdr(:) = dEdr * cache%dmrdcn

   call gemv(cache%dcndr, dEdr, gradient, beta=1.0_wp)
   call gemv(cache%dcndL, dEdr, sigma, beta=1.0_wp)
end subroutine get_gradient


subroutine get_mrad(mol, shift, kexp, rmax, rad, valence_cn, cn, mrad, dmrdcn)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: shift
   real(wp), intent(in) :: kexp
   real(wp), intent(in) :: rmax
   real(wp), intent(in) :: rad(:)
   real(wp), intent(in) :: valence_cn(:)
   real(wp), intent(in) :: cn(:)
   real(wp), intent(out) :: mrad(:)
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


subroutine get_dir_trans(lattice, trans)
   real(wp), intent(in) :: lattice(:, :)
   real(wp), allocatable, intent(out) :: trans(:, :)
   integer, parameter :: rep(3) = 2

   call get_lattice_points(lattice, rep, .true., trans)

end subroutine get_dir_trans

subroutine get_rec_trans(lattice, trans)
   real(wp), intent(in) :: lattice(:, :)
   real(wp), allocatable, intent(out) :: trans(:, :)
   integer, parameter :: rep(3) = 2
   real(wp) :: rec_lat(3, 3)

   rec_lat = twopi*transpose(matinv_3x3(lattice))
   call get_lattice_points(rec_lat, rep, .false., trans)

end subroutine get_rec_trans


subroutine get_charge_dipole_matrix(self, mol, cache, amat)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   real(wp), intent(inout) :: amat(:, :, :)

   if (any(mol%periodic)) then
      call get_charge_dipole_matrix_3d(mol, cache%mrad, self%kdmp3, &
         & cache%wsc, cache%alpha, amat)
   else
      call get_charge_dipole_matrix_0d(mol, cache%mrad, self%kdmp3, amat)
   end if
end subroutine get_charge_dipole_matrix

subroutine get_charge_dipole_matrix_0d(mol, rad, kdmp, amat)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: rad(:)
   real(wp), intent(in) :: kdmp
   real(wp), intent(inout) :: amat(:, :, :)

   integer :: iat, jat
   real(wp) :: r1, r2, vec(3), rr, fdmp, g3

   do iat = 1, mol%nat
      do jat = 1, mol%nat
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         r1 = norm2(vec)
         if (r1 < epsilon(1.0_wp)) cycle
         r2 = r1*r1

         rr = 0.5_wp * (rad(jat) + rad(iat)) / r1
         fdmp = 1.0_wp / (1.0_wp + 6.0_wp * rr**kdmp)

         g3 = fdmp / (r1 * r2)
         amat(:, jat, iat) = amat(:, jat, iat) + vec * g3
      end do
   end do

end subroutine get_charge_dipole_matrix_0d

subroutine get_charge_dipole_matrix_3d(mol, rad, kdmp, wsc, alpha, amat)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: rad(:)
   real(wp), intent(in) :: kdmp
   type(wignerseitz_cell), intent(in) :: wsc
   real(wp), intent(in) :: alpha
   real(wp), intent(inout) :: amat(:, :, :)

   integer :: iat, jat, img
   real(wp) :: vec(3), rr, wsw, dtmp(3), rtmp(3), vol
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, dtrans)
   call get_rec_trans(mol%lattice, rtrans)

   do iat = 1, mol%nat
      do jat = 1, mol%nat
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))

            rr = 0.5_wp * (rad(jat) + rad(iat))
            call get_amat_sd_rec_3d(vec, vol, alpha, rtrans, rtmp)
            call get_amat_sd_dir_3d(vec, rr, kdmp, alpha, dtrans, dtmp)

            amat(:, jat, iat) = amat(:, jat, iat) + wsw * (dtmp + rtmp)
         end do
      end do
   end do

end subroutine get_charge_dipole_matrix_3d

pure subroutine get_amat_sd_dir_3d(rij, rr, kdmp, alp, trans, amat)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: rr
   real(wp), intent(in) :: kdmp
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: amat(:)

   integer :: itr
   real(wp) :: vec(3), r1, tmp, fdmp, g1, g3, arg, arg2, e1

   amat = 0.0_wp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      g1 = 1.0_wp/r1
      g3 = g1 * g1 * g1
      fdmp = 1.0_wp / (1.0_wp + 6.0_wp * (rr/r1)**kdmp)

      arg = r1*alp
      arg2 = arg*arg
      e1 = 2/sqrtpi*arg*exp(-arg2) - erf(arg)

      tmp = fdmp * g3 + e1 * g3
      amat = amat + vec * tmp
   end do

end subroutine get_amat_sd_dir_3d

pure subroutine get_amat_sd_rec_3d(rij, vol, alp, trans, amat)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: vol
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: amat(:)

   integer :: itr
   real(wp) :: fac, vec(3), g2, expk, sink

   amat = 0.0_wp
   fac = 4*pi/vol

   do itr = 1, size(trans, 2)
      vec(:) = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      expk = fac*exp(-0.25_wp*g2/(alp*alp))/g2
      sink = sin(dot_product(rij, vec))*expk

      amat(:) = amat + 2*vec*sink
   end do

end subroutine get_amat_sd_rec_3d


subroutine get_dipole_dipole_matrix(self, mol, cache, amat)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   real(wp), intent(inout) :: amat(:, :, :, :)

   if (any(mol%periodic)) then
      call get_dipole_dipole_matrix_3d(mol, cache%mrad, self%kdmp5, &
         & cache%wsc, cache%alpha, amat)
   else
      call get_dipole_dipole_matrix_0d(mol, cache%mrad, self%kdmp5, amat)
   end if
end subroutine get_dipole_dipole_matrix

subroutine get_dipole_dipole_matrix_0d(mol, rad, kdmp, amat)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: rad(:)
   real(wp), intent(in) :: kdmp
   real(wp), intent(out) :: amat(:, :, :, :)

   integer :: iat, jat
   real(wp) :: r1, r2, vec(3), rr, fdmp, g3, g5

   do iat = 1, mol%nat
      do jat = 1, mol%nat
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         r1 = norm2(vec)
         if (r1 < epsilon(1.0_wp)) cycle
         r2 = r1*r1

         rr = 0.5_wp * (rad(jat) + rad(iat)) / r1
         fdmp = 1.0_wp / (1.0_wp + 6.0_wp * rr**kdmp)

         g3 = fdmp / (r1 * r2)
         g5 = 3 * g3 / r2
         amat(:, jat, :, iat) = amat(:, jat, :, iat) &
            & + unity * g3 - spread(vec, 1, 3) * spread(vec, 2, 3) * g5
      end do
   end do
end subroutine get_dipole_dipole_matrix_0d

subroutine get_dipole_dipole_matrix_3d(mol, rad, kdmp, wsc, alpha, amat)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: rad(:)
   real(wp), intent(in) :: kdmp
   type(wignerseitz_cell), intent(in) :: wsc
   real(wp), intent(in) :: alpha
   real(wp), intent(out) :: amat(:, :, :, :)

   integer :: iat, jat, img, k
   real(wp) :: vec(3), rr, wsw, dtmp(3, 3), rtmp(3, 3), vol
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, dtrans)
   call get_rec_trans(mol%lattice, rtrans)

   do iat = 1, mol%nat
      do jat = 1, mol%nat
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))

            rr = 0.5_wp * (rad(jat) + rad(iat))
            call get_amat_dd_rec_3d(vec, vol, alpha, rtrans, rtmp)
            call get_amat_dd_dir_3d(vec, rr, kdmp, alpha, dtrans, dtmp)

            amat(:, jat, :, iat) = amat(:, jat, :, iat) + wsw * (rtmp + dtmp)
         end do
      end do
      ! dipole-dipole selfenergy: 2/3·α³/sqrt(π) Σ(i) μ²(i)
      rr = 2.0_wp/3.0_wp * alpha**3 / sqrtpi
      do k = 1, 3
         amat(k, iat, k, iat) = amat(k, iat, k, iat) + 2*rr
      end do
   end do
end subroutine get_dipole_dipole_matrix_3d

pure subroutine get_amat_dd_dir_3d(rij, rr, kdmp, alp, trans, amat)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: rr
   real(wp), intent(in) :: kdmp
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: amat(:, :)

   integer :: itr
   real(wp) :: vec(3), r1, tmp, fdmp, g1, g3, g5, arg, arg2, e1

   amat = 0.0_wp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      g1 = 1.0_wp/r1
      g3 = g1 * g1 * g1
      g5 = g3 * g1 * g1
      fdmp = 1.0_wp / (1.0_wp + 6.0_wp * (rr/r1)**kdmp)

      arg = r1*alp
      arg2 = arg*arg
      e1 = 2/sqrtpi*arg*exp(-arg2) - erf(arg)

      amat(:, :) = amat(:, :) + unity * (fdmp + e1) * g3 &
         & - spread(vec, 1, 3) * spread(vec, 2, 3) * 3*g5 * (fdmp + e1)
   end do

end subroutine get_amat_dd_dir_3d

pure subroutine get_amat_dd_rec_3d(rij, vol, alp, trans, amat)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: vol
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: amat(:, :)

   integer :: itr
   real(wp) :: fac, vec(3), g2, expk, cosk, gv

   amat = 0.0_wp
   fac = 4*pi/vol

   do itr = 1, size(trans, 2)
      vec(:) = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      gv = dot_product(rij, vec)
      expk = fac*exp(-0.25_wp*g2/(alp*alp))/g2
      cosk = cos(gv)*expk

      amat(:, :) = amat + spread(vec, 1, 3) * spread(vec, 2, 3) * cosk
   end do

end subroutine get_amat_dd_rec_3d


subroutine get_charge_quadrupole_matrix(self, mol, cache, amat)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   real(wp), intent(inout) :: amat(:, :, :)

   if (any(mol%periodic)) then
      call get_charge_quadrupole_matrix_3d(mol, cache%mrad, self%kdmp5, &
         & cache%wsc, cache%alpha, amat)
   else
      call get_charge_quadrupole_matrix_0d(mol, cache%mrad, self%kdmp5, amat)
   end if
end subroutine get_charge_quadrupole_matrix

subroutine get_charge_quadrupole_matrix_0d(mol, rad, kdmp, amat)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: rad(:)
   real(wp), intent(in) :: kdmp
   real(wp), intent(inout) :: amat(:, :, :)

   integer :: iat, jat
   real(wp) :: r1, r2, vec(3), g1, g3, g5, fdmp, arg, tc(6), rr

   do iat = 1, mol%nat
      do jat = 1, mol%nat
         if (iat == jat) cycle
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         g1 = 1.0_wp / r1
         g3 = g1 * g1 * g1
         g5 = g3 * g1 * g1

         rr = 0.5_wp * (rad(jat) + rad(iat)) * g1
         fdmp = 1.0_wp / (1.0_wp + 6.0_wp * rr**kdmp)
         tc(2) = 2*vec(1)*vec(2)*g5*fdmp
         tc(4) = 2*vec(1)*vec(3)*g5*fdmp
         tc(5) = 2*vec(2)*vec(3)*g5*fdmp
         tc(1) = vec(1)*vec(1)*g5*fdmp
         tc(3) = vec(2)*vec(2)*g5*fdmp
         tc(6) = vec(3)*vec(3)*g5*fdmp
         amat(:, jat, iat) = amat(:, jat, iat) + tc
      end do
   end do
end subroutine get_charge_quadrupole_matrix_0d

subroutine get_charge_quadrupole_matrix_3d(mol, rad, kdmp, wsc, alpha, amat)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: rad(:)
   real(wp), intent(in) :: kdmp
   type(wignerseitz_cell), intent(in) :: wsc
   real(wp), intent(in) :: alpha
   real(wp), intent(inout) :: amat(:, :, :)

   integer :: iat, jat, img
   real(wp) :: vec(3), rr, wsw, dtmp(6), rtmp(6), vol
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, dtrans)
   call get_rec_trans(mol%lattice, rtrans)

   do iat = 1, mol%nat
      do jat = 1, mol%nat
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))

            rr = 0.5_wp * (rad(jat) + rad(iat))
            call get_amat_sq_rec_3d(vec, vol, alpha, rtrans, rtmp)
            call get_amat_sq_dir_3d(vec, rr, kdmp, alpha, dtrans, dtmp)

            amat(:, jat, iat) = amat(:, jat, iat) + wsw * (rtmp + dtmp)
         end do
      end do
      ! charge-quadrupole selfenergy: 2/3·α³/sqrt(π) Σ(i) q(i)Tr(θi)
      rr = 4.0_wp/9.0_wp * alpha**3 / sqrtpi
      amat([1, 3, 6], iat, iat) = amat([1, 3, 6], iat, iat) + rr
   end do
end subroutine get_charge_quadrupole_matrix_3d

pure subroutine get_amat_sq_dir_3d(rij, rr, kdmp, alp, trans, amat)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: rr
   real(wp), intent(in) :: kdmp
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: amat(:)

   integer :: itr
   real(wp) :: vec(3), r1, tmp, fdmp, g1, g3, g5, arg, arg2, e1

   amat = 0.0_wp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      g1 = 1.0_wp/r1
      g3 = g1 * g1 * g1
      g5 = g3 * g1 * g1
      fdmp = 1.0_wp / (1.0_wp + 6.0_wp * (rr/r1)**kdmp)

      arg = r1*alp
      arg2 = arg*arg
      e1 = 2/sqrtpi*arg*exp(-arg2) - erf(arg)

      amat(1) = amat(1) +   vec(1)*vec(1)*g5*(fdmp + e1)
      amat(2) = amat(2) + 2*vec(1)*vec(2)*g5*(fdmp + e1)
      amat(3) = amat(3) +   vec(2)*vec(2)*g5*(fdmp + e1)
      amat(4) = amat(4) + 2*vec(1)*vec(3)*g5*(fdmp + e1)
      amat(5) = amat(5) + 2*vec(2)*vec(3)*g5*(fdmp + e1)
      amat(6) = amat(6) +   vec(3)*vec(3)*g5*(fdmp + e1)
   end do

end subroutine get_amat_sq_dir_3d

pure subroutine get_amat_sq_rec_3d(rij, vol, alp, trans, amat)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: vol
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: amat(:)

   integer :: itr
   real(wp) :: fac, vec(3), g2, expk, cosk, gtmp(3), gv

   amat = 0.0_wp
   fac = 4*pi/vol

   do itr = 1, size(trans, 2)
      vec(:) = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      gv = dot_product(rij, vec)
      expk = fac*exp(-0.25_wp*g2/(alp*alp))/g2
      cosk = cos(gv)*expk

      amat(1) = amat(1) +   vec(1)*vec(1)*cosk
      amat(2) = amat(2) + 2*vec(1)*vec(2)*cosk
      amat(3) = amat(3) +   vec(2)*vec(2)*cosk
      amat(4) = amat(4) + 2*vec(1)*vec(3)*cosk
      amat(5) = amat(5) + 2*vec(2)*vec(3)*cosk
      amat(6) = amat(6) +   vec(3)*vec(3)*cosk
   end do

end subroutine get_amat_sq_rec_3d


subroutine get_charge_dipole_gradient(self, mol, cache, qat, dpat, dEdr, gradient, sigma)
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
   !> Derivative of the energy w.r.t. the critical radii
   real(wp), contiguous, intent(inout) :: dEdr(:)
   !> Derivative of the energy w.r.t. atomic displacements
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Derivative of the energy w.r.t. strain deformations
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   if (any(mol%periodic)) then
      call get_charge_dipole_gradient_3d(mol, cache%mrad, self%kdmp3, qat, dpat, &
         & cache%wsc, cache%alpha, dEdr, gradient, sigma)
   else
      call get_charge_dipole_gradient_0d(mol, cache%mrad, self%kdmp3, qat, dpat, &
         & dEdr, gradient, sigma)
   end if
end subroutine get_charge_dipole_gradient

subroutine get_charge_dipole_gradient_0d(mol, rad, kdmp, qat, dpat, dEdr, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: rad(:)
   real(wp), intent(in) :: kdmp
   !> Atomic partial charges
   real(wp), contiguous, intent(in) :: qat(:)
   !> Atomic dipole moments
   real(wp), contiguous, intent(in) :: dpat(:, :)
   !> Derivative of the energy w.r.t. the critical radii
   real(wp), contiguous, intent(inout) :: dEdr(:)
   !> Derivative of the energy w.r.t. atomic displacements
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Derivative of the energy w.r.t. strain deformations
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   integer :: iat, jat
   real(wp) :: r1, vec(3), rr, fdmp, g1, g3, g5, dG(3), dS(3, 3)
   real(wp) :: dpiqj, qidpj, ddmp, fddr

   do iat = 1, mol%nat
      do jat = 1, iat - 1
         rr = 0.5_wp*(rad(iat)+rad(jat))
         vec(:) = mol%xyz(:, jat)-mol%xyz(:, iat)
         r1 = norm2(vec)
         g1 = 1.0_wp/r1
         g3 = g1 * g1 * g1
         g5 = g3 * g1 * g1
         fdmp = 1.0_wp/(1.0_wp+6.0_wp*(rr*g1)**kdmp)
         ddmp = -3*g5*fdmp - kdmp*fdmp*(fdmp-1.0_wp)*g5
         dpiqj = dot_product(vec, dpat(:, iat))*qat(jat)
         qidpj = dot_product(vec, dpat(:, jat))*qat(iat)
         fddr = 3.0_wp*(dpiqj - qidpj)*kdmp*fdmp*g3*(fdmp/rr)*(rr*g1)**kdmp
         dg(:) = - ddmp*vec * (dpiqj - qidpj) &
            & + fdmp*g3*(qat(iat)*dpat(:, jat) - qat(jat)*dpat(:, iat))
         ds(:, :) = - 0.5_wp * (spread(vec, 1, 3) * spread(dG, 2, 3) &
            & + spread(dG, 1, 3) * spread(vec, 2, 3))

         dEdr(iat) = dEdr(iat) + fddr
         dEdr(jat) = dEdr(jat) + fddr
         gradient(:, iat) = gradient(:, iat) + dG
         gradient(:, jat) = gradient(:, jat) - dG
         sigma(:, :) = sigma + dS
      end do
   end do

end subroutine get_charge_dipole_gradient_0d

subroutine get_charge_dipole_gradient_3d(mol, rad, kdmp, qat, dpat, wsc, alpha, &
      & dEdr, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: rad(:)
   real(wp), intent(in) :: kdmp
   !> Atomic partial charges
   real(wp), contiguous, intent(in) :: qat(:)
   !> Atomic dipole moments
   real(wp), contiguous, intent(in) :: dpat(:, :)
   type(wignerseitz_cell), intent(in) :: wsc
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
   call get_dir_trans(mol%lattice, dtrans)
   call get_rec_trans(mol%lattice, rtrans)

   do iat = 1, mol%nat
      do jat = 1, iat - 1
         dE = 0.0_wp
         dG(:) = 0.0_wp
         dS(:, :) = 0.0_wp
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec(:) = mol%xyz(:, jat)-mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, iat))

            rr = 0.5_wp * (rad(jat) + rad(iat))
            call get_damat_sd_rec_3d(vec, qat(iat), qat(jat), dpat(:, iat), dpat(:, jat), &
               & vol, alpha, rtrans, dGr, dSr)
            call get_damat_sd_dir_3d(vec, qat(iat), qat(jat), dpat(:, iat), dpat(:, jat), &
               & rr, kdmp, alpha, dtrans, dEd, dGd, dSd)
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

      dE = 0.0_wp
      dS(:, :) = 0.0_wp
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec(:) = wsc%trans(:, wsc%tridx(img, iat, iat))

         rr = rad(iat)
         call get_damat_sd_rec_3d(vec, qat(iat), qat(iat), dpat(:, iat), dpat(:, iat), &
            & vol, alpha, rtrans, dGr, dSr)
         call get_damat_sd_dir_3d(vec, qat(iat), qat(iat), dpat(:, iat), dpat(:, iat), &
            & rr, kdmp, alpha, dtrans, dEd, dGd, dSd)
         dE = dE + dEd * wsw
         dS = dS + (dSd + dSr) * wsw
      end do
      dEdr(iat) = dEdr(iat) + dE
      sigma = sigma + 0.5_wp * dS
   end do

end subroutine get_charge_dipole_gradient_3d

pure subroutine get_damat_sd_dir_3d(rij, qi, qj, mi, mj, rr, kdmp, alp, trans, de, dg, ds)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: qi, qj, mi(3), mj(3)
   real(wp), intent(in) :: rr
   real(wp), intent(in) :: kdmp
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: de
   real(wp), intent(out) :: dg(3)
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: vec(3), r1, g1, g3, g5, fdmp, ddmp, dpiqj, qidpj, gtmp(3), atmp, alp2
   real(wp) :: arg, arg2, erft, expt, e1, de1

   de = 0.0_wp
   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   alp2 = alp*alp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      g1 = 1.0_wp/r1
      g3 = g1 * g1 * g1
      g5 = g3 * g1 * g1

      arg = r1*alp
      arg2 = arg*arg
      erft = erf(arg)
      expt = exp(-arg2)
      e1 = 2/sqrtpi*arg*expt - erft
      de1 = 3*erft - expt*(4/sqrtpi*arg*arg2 + 6/sqrtpi*arg)

      fdmp = 1.0_wp/(1.0_wp+6.0_wp*(rr*g1)**kdmp)
      ddmp = -3*fdmp - kdmp*fdmp*(fdmp-1.0_wp)
      dpiqj = dot_product(vec, mi)*qj
      qidpj = dot_product(vec, mj)*qi
      gtmp(:) = - (ddmp+de1)*g5*vec * (dpiqj - qidpj) &
         & + (fdmp+e1)*g3*(qi*mj - qj*mi)

      de = de + 3.0_wp*(dpiqj - qidpj)*kdmp*fdmp*g3*(fdmp/rr)*(rr*g1)**kdmp
      dg(:) = dg + gtmp
      ds(:, :) = - 0.5_wp * (spread(vec, 1, 3) * spread(gtmp, 2, 3) &
         & + spread(gtmp, 1, 3) * spread(vec, 2, 3))
   end do

end subroutine get_damat_sd_dir_3d

pure subroutine get_damat_sd_rec_3d(rij, qi, qj, mi, mj, vol, alp, trans, dg, ds)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: qi, qj, mi(3), mj(3)
   real(wp), intent(in) :: vol
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: dg(3)
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: fac, vec(3), g2, gv, expk, alp2, sink, cosk, dpiqj, qidpj

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
   end do

end subroutine get_damat_sd_rec_3d


subroutine get_dipole_dipole_gradient(self, mol, cache, dpat, dEdr, gradient, sigma)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   !> Atomic dipole moments
   real(wp), contiguous, intent(in) :: dpat(:, :)
   !> Derivative of the energy w.r.t. the critical radii
   real(wp), contiguous, intent(inout) :: dEdr(:)
   !> Derivative of the energy w.r.t. atomic displacements
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Derivative of the energy w.r.t. strain deformations
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   if (any(mol%periodic)) then
      call get_dipole_dipole_gradient_3d(mol, cache%mrad, self%kdmp5, dpat, &
         & cache%wsc, cache%alpha, dEdr, gradient, sigma)
   else
      call get_dipole_dipole_gradient_0d(mol, cache%mrad, self%kdmp5, dpat, &
         & dEdr, gradient, sigma)
   end if
end subroutine get_dipole_dipole_gradient

subroutine get_dipole_dipole_gradient_0d(mol, rad, kdmp, dpat, dEdr, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: rad(:)
   real(wp), intent(in) :: kdmp
   !> Atomic dipole moments
   real(wp), contiguous, intent(in) :: dpat(:, :)
   !> Derivative of the energy w.r.t. the critical radii
   real(wp), contiguous, intent(inout) :: dEdr(:)
   !> Derivative of the energy w.r.t. atomic displacements
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Derivative of the energy w.r.t. strain deformations
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   integer :: iat, jat
   real(wp) :: r1, r2, vec(3), rr, fdmp, g1, g3, g5, g7, dg(3), ds(3, 3)
   real(wp) :: ddmp, edd, dpidpj, dpiv, dpjv, fddr

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
         fdmp = 1.0_wp/(1.0_wp+6.0_wp*(rr*g1)**kdmp)
         ddmp = -5*fdmp - kdmp*(fdmp*fdmp-fdmp)
         dpidpj = dot_product(dpat(:, jat), dpat(:, iat))
         dpiv = dot_product(dpat(:, iat), vec)
         dpjv = dot_product(dpat(:, jat), vec)
         edd = dpidpj*r2 - 3*dpjv*dpiv
         fddr = 3.0_wp*edd*kdmp*fdmp*g5*(fdmp/rr)*(rr*g1)**kdmp
         dg(:) = - 2.0_wp*fdmp*g5*dpidpj*vec &
            & + 3.0_wp*fdmp*g5*(dpiv*dpat(:, jat) + dpjv*dpat(:, iat)) &
            & - edd*ddmp*g7*vec
         ds(:, :) = - 0.5_wp * (spread(vec, 1, 3) * spread(dg, 2, 3) &
            & + spread(dg, 1, 3) * spread(vec, 2, 3))

         dEdr(iat) = dEdr(iat) + fddr
         dEdr(jat) = dEdr(jat) + fddr
         gradient(:, iat) = gradient(:, iat) + dg
         gradient(:, jat) = gradient(:, jat) - dg
         sigma(:, :) = sigma + ds
      end do
   end do
end subroutine get_dipole_dipole_gradient_0d

subroutine get_dipole_dipole_gradient_3d(mol, rad, kdmp, dpat, wsc, alpha, dEdr, &
      & gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: rad(:)
   real(wp), intent(in) :: kdmp
   !> Atomic dipole moments
   real(wp), contiguous, intent(in) :: dpat(:, :)
   type(wignerseitz_cell), intent(in) :: wsc
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
   call get_dir_trans(mol%lattice, dtrans)
   call get_rec_trans(mol%lattice, rtrans)

   do iat = 1, mol%nat
      do jat = 1, iat - 1
         dE = 0.0_wp
         dG(:) = 0.0_wp
         dS(:, :) = 0.0_wp
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, iat))

            rr = 0.5_wp * (rad(jat) + rad(iat))
            call get_damat_dd_rec_3d(vec, dpat(:, iat), dpat(:, jat), &
               & vol, alpha, rtrans, dGr, dSr)
            call get_damat_dd_dir_3d(vec, dpat(:, iat), dpat(:, jat), &
               & rr, kdmp, alpha, dtrans, dEd, dGd, dSd)
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

      dE = 0.0_wp
      dS(:, :) = 0.0_wp
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))

         rr = rad(iat)
         call get_damat_dd_rec_3d(vec, dpat(:, iat), dpat(:, iat), &
            & vol, alpha, rtrans, dGr, dSr)
         call get_damat_dd_dir_3d(vec, dpat(:, iat), dpat(:, iat), &
            & rr, kdmp, alpha, dtrans, dEd, dGd, dSd)
         dE = dE + dEd * wsw
         dS = dS + (dSd + dSr) * wsw
      end do
      dEdr(iat) = dEdr(iat) + dE
      sigma = sigma + 0.5_wp * dS
   end do
end subroutine get_dipole_dipole_gradient_3d

pure subroutine get_damat_dd_dir_3d(rij, mi, mj, rr, kdmp, alp, trans, de, dg, ds)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: mi(3)
   real(wp), intent(in) :: mj(3)
   real(wp), intent(in) :: rr
   real(wp), intent(in) :: kdmp
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: de
   real(wp), intent(out) :: dg(3)
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: vec(3), r1, r2, g1, g3, g5, g7, fdmp, ddmp, dpidpj, dpiv, dpjv, edd
   real(wp) :: gtmp(3), atmp, alp2, arg, arg2, erft, expt, e1, de1

   de = 0.0_wp
   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   alp2 = alp*alp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      arg = r1*alp
      arg2 = arg*arg
      erft = erf(arg)
      expt = exp(-arg2)
      e1 = 2/sqrtpi*arg*expt - erft
      de1 = 5*erft - expt*(4/sqrtpi*arg*arg2 + 10/sqrtpi*arg)

      r2 = r1 * r1
      g1 = 1.0_wp/r1
      g3 = g1 * g1 * g1
      g5 = g3 * g1 * g1
      g7 = g5 * g1 * g1
      fdmp = 1.0_wp/(1.0_wp+6.0_wp*(rr*g1)**kdmp)
      ddmp = -5*fdmp - kdmp*(fdmp*fdmp-fdmp)
      dpidpj = dot_product(mj, mi)
      dpiv = dot_product(mi, vec)
      dpjv = dot_product(mj, vec)
      edd = dpidpj*r2 - 3*dpjv*dpiv
      gtmp(:) = - 2.0_wp*(fdmp+e1)*g5*dpidpj*vec &
         & + 3.0_wp*(fdmp+e1)*g5*(dpiv*mj + dpjv*mi) &
         & - edd*(ddmp+de1)*g7*vec

      de = de + 3.0_wp*edd*kdmp*fdmp*g5*(fdmp/rr)*(rr*g1)**kdmp
      dg(:) = dg + gtmp
      ds(:, :) = - 0.5_wp * (spread(vec, 1, 3) * spread(gtmp, 2, 3) &
         & + spread(gtmp, 1, 3) * spread(vec, 2, 3))
   end do

end subroutine get_damat_dd_dir_3d

pure subroutine get_damat_dd_rec_3d(rij, mi, mj, vol, alp, trans, dg, ds)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: mi(3)
   real(wp), intent(in) :: mj(3)
   real(wp), intent(in) :: vol
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: dg(3)
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: fac, vec(3), g2, gv, expk, alp2, mikmjk, cosk, sink, gtmp(3)
   real(wp) :: dpiv, dpjv

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

      dpiv = dot_product(mi, vec)
      dpjv = dot_product(mj, vec)

      dg(:) = dg + vec*sink*dpiv*dpjv
      ds(:, :) = ds + cosk * dpiv*dpjv &
         & * ((2.0_wp/g2 + 0.5_wp/alp2) * spread(vec, 1, 3)*spread(vec, 2, 3) - unity)
   end do

end subroutine get_damat_dd_rec_3d


subroutine get_charge_quadrupole_gradient(self, mol, cache, qat, qpat, dEdr, gradient, sigma)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   !> Atomic partial charges
   real(wp), contiguous, intent(in) :: qat(:)
   !> Atomic quadrupole moments
   real(wp), contiguous, intent(in) :: qpat(:, :)
   !> Derivative of the energy w.r.t. the critical radii
   real(wp), contiguous, intent(inout) :: dEdr(:)
   !> Derivative of the energy w.r.t. atomic displacements
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Derivative of the energy w.r.t. strain deformations
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   if (any(mol%periodic)) then
      call get_charge_quadrupole_gradient_3d(mol, cache%mrad, self%kdmp5, qat, qpat, &
         & cache%wsc, cache%alpha, dEdr, gradient, sigma)
   else
      call get_charge_quadrupole_gradient_0d(mol, cache%mrad, self%kdmp5, qat, qpat, &
         & dEdr, gradient, sigma)
   end if
end subroutine get_charge_quadrupole_gradient

subroutine get_charge_quadrupole_gradient_0d(mol, rad, kdmp, qat, qpat, dEdr, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: rad(:)
   real(wp), intent(in) :: kdmp
   !> Atomic partial charges
   real(wp), contiguous, intent(in) :: qat(:)
   !> Atomic quadrupole moments
   real(wp), contiguous, intent(in) :: qpat(:, :)
   !> Derivative of the energy w.r.t. the critical radii
   real(wp), contiguous, intent(inout) :: dEdr(:)
   !> Derivative of the energy w.r.t. atomic displacements
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Derivative of the energy w.r.t. strain deformations
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   integer :: iat, jat
   real(wp) :: r1, r2, vec(3), rr, fdmp, g1, g3, g5, g7, dG(3), dS(3, 3), ddmp, fddr, eq

   do iat = 1, mol%nat
      do jat = 1, iat - 1
         vec(:) = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         r2 = r1 * r1
         rr = 0.5_wp*(rad(iat)+rad(jat))
         g1 = 1.0_wp / r1
         g3 = g1 * g1 * g1
         g5 = g3 * g1 * g1
         g7 = g5 * g1 * g1
         fdmp = 1.0_wp/(1.0_wp + 6.0_wp*(rr*g1)**kdmp)
         ddmp = -5*fdmp - kdmp*(fdmp*fdmp - fdmp)
         eq = &
            & + 2*(qat(jat)*qpat(2,iat) + qpat(2,jat)*qat(iat))*vec(1)*vec(2) &
            & + 2*(qat(jat)*qpat(4,iat) + qpat(4,jat)*qat(iat))*vec(1)*vec(3) &
            & + 2*(qat(jat)*qpat(5,iat) + qpat(5,jat)*qat(iat))*vec(2)*vec(3) &
            & + (qat(jat)*qpat(1,iat) + qpat(1,jat)*qat(iat))*vec(1)*vec(1) &
            & + (qat(jat)*qpat(3,iat) + qpat(3,jat)*qat(iat))*vec(2)*vec(2) &
            & + (qat(jat)*qpat(6,iat) + qpat(6,jat)*qat(iat))*vec(3)*vec(3)

         fddr = eq * 3.0_wp*kdmp*fdmp*g5*fdmp/rr*(rr*g1)**kdmp
         dg(:) = - eq*ddmp*g7*vec &
            & - 2.0_wp*fdmp*g5*qat(iat) * &
            &[vec(1)*qpat(1,jat) + vec(2)*qpat(2,jat) + vec(3)*qpat(4,jat), &
            & vec(1)*qpat(2,jat) + vec(2)*qpat(3,jat) + vec(3)*qpat(5,jat), &
            & vec(1)*qpat(4,jat) + vec(2)*qpat(5,jat) + vec(3)*qpat(6,jat)] &
            & - 2.0_wp*fdmp*g5*qat(jat) * &
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
end subroutine get_charge_quadrupole_gradient_0d

subroutine get_charge_quadrupole_gradient_3d(mol, rad, kdmp, qat, qpat, wsc, alpha, &
      & dEdr, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: rad(:)
   real(wp), intent(in) :: kdmp
   !> Atomic partial charges
   real(wp), contiguous, intent(in) :: qat(:)
   !> Atomic quadrupole moments
   real(wp), contiguous, intent(in) :: qpat(:, :)
   type(wignerseitz_cell), intent(in) :: wsc
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
   call get_dir_trans(mol%lattice, dtrans)
   call get_rec_trans(mol%lattice, rtrans)

   do iat = 1, mol%nat
      do jat = 1, iat - 1
         dE = 0.0_wp
         dG(:) = 0.0_wp
         dS(:, :) = 0.0_wp
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, iat))

            rr = 0.5_wp * (rad(jat) + rad(iat))
            call get_damat_sq_rec_3d(vec, qat(iat), qat(jat), qpat(:, iat), qpat(:, jat), &
               & vol, alpha, rtrans, dGr, dSr)
            call get_damat_sq_dir_3d(vec, qat(iat), qat(jat), qpat(:, iat), qpat(:, jat), &
               & rr, kdmp, alpha, dtrans, dEd, dGd, dSd)
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

      dE = 0.0_wp
      dS(:, :) = 0.0_wp
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))

         rr = rad(iat)
         call get_damat_sq_rec_3d(vec, qat(iat), qat(iat), qpat(:, iat), qpat(:, iat), &
            & vol, alpha, rtrans, dGr, dSr)
         call get_damat_sq_dir_3d(vec, qat(iat), qat(iat), qpat(:, iat), qpat(:, iat), &
            & rr, kdmp, alpha, dtrans, dEd, dGd, dSd)
         dE = dE + dEd * wsw
         dS = dS + (dSd + dSr) * wsw
      end do
      dEdr(iat) = dEdr(iat) + dE
      sigma = sigma + 0.5_wp * dS
   end do
end subroutine get_charge_quadrupole_gradient_3d

pure subroutine get_damat_sq_dir_3d(rij, qi, qj, ti, tj, rr, kdmp, alp, trans, de, dg, ds)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: qi, qj, ti(6), tj(6)
   real(wp), intent(in) :: rr
   real(wp), intent(in) :: kdmp
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: de
   real(wp), intent(out) :: dg(3)
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: vec(3), r1, r2, g1, g3, g5, g7, fdmp, ddmp, eq
   real(wp) :: gtmp(3), atmp, alp2, arg, arg2, erft, expt, e1, de1

   de = 0.0_wp
   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   alp2 = alp*alp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      arg = r1*alp
      arg2 = arg*arg
      erft = erf(arg)
      expt = exp(-arg2)
      e1 = 2/sqrtpi*arg*expt - erft
      de1 = 5*erft - expt*(4/sqrtpi*arg*arg2 + 10/sqrtpi*arg)

      r2 = r1 * r1
      g1 = 1.0_wp / r1
      g3 = g1 * g1 * g1
      g5 = g3 * g1 * g1
      g7 = g5 * g1 * g1
      fdmp = 1.0_wp/(1.0_wp + 6.0_wp*(rr*g1)**kdmp)
      ddmp = -5*fdmp - kdmp*(fdmp*fdmp - fdmp)
      eq = &
         & + 2*(qj*ti(2) + tj(2)*qi)*vec(1)*vec(2) &
         & + 2*(qj*ti(4) + tj(4)*qi)*vec(1)*vec(3) &
         & + 2*(qj*ti(5) + tj(5)*qi)*vec(2)*vec(3) &
         & + (qj*ti(1) + tj(1)*qi)*vec(1)*vec(1) &
         & + (qj*ti(3) + tj(3)*qi)*vec(2)*vec(2) &
         & + (qj*ti(6) + tj(6)*qi)*vec(3)*vec(3)
      gtmp(:) = - eq*(ddmp+de1)*g7*vec & 
         & - 2.0_wp*(fdmp+e1)*g5*qi * &
         &[vec(1)*tj(1) + vec(2)*tj(2) + vec(3)*tj(4), &
         & vec(1)*tj(2) + vec(2)*tj(3) + vec(3)*tj(5), &
         & vec(1)*tj(4) + vec(2)*tj(5) + vec(3)*tj(6)] &
         & - 2.0_wp*(fdmp+e1)*g5*qj * &
         &[vec(1)*ti(1) + vec(2)*ti(2) + vec(3)*ti(4), &
         & vec(1)*ti(2) + vec(2)*ti(3) + vec(3)*ti(5), &
         & vec(1)*ti(4) + vec(2)*ti(5) + vec(3)*ti(6)]

      de = de + eq * 3.0_wp*kdmp*fdmp*g5*fdmp/rr*(rr*g1)**kdmp
      dg(:) = dg + gtmp
      ds(:, :) = ds - 0.5_wp * (spread(vec, 1, 3) * spread(gtmp, 2, 3) &
         & + spread(gtmp, 1, 3) * spread(vec, 2, 3))
   end do

end subroutine get_damat_sq_dir_3d

pure subroutine get_damat_sq_rec_3d(rij, qi, qj, ti, tj, vol, alp, trans, dg, ds)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: qi, qj, ti(6), tj(6)
   real(wp), intent(in) :: vol
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: dg(3)
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: fac, vec(3), g2, gv, expk, alp2, mikmjk, cosk, sink, gtmp(3)
   real(wp) :: qpiqj, qiqpj

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

      qiqpj = qi*(tj(1)*vec(1)*vec(1) + tj(3)*vec(2)*vec(2) + tj(6)*vec(3)*vec(3) &
         & + 2*tj(2)*vec(1)*vec(2) + 2*tj(4)*vec(1)*vec(3) + 2*tj(5)*vec(2)*vec(3))
      qpiqj = qj*(ti(1)*vec(1)*vec(1) + ti(3)*vec(2)*vec(2) + ti(6)*vec(3)*vec(3) &
         & + 2*ti(2)*vec(1)*vec(2) + 2*ti(4)*vec(1)*vec(3) + 2*ti(5)*vec(2)*vec(3))

      dg(:) = dg + vec * sink * (qiqpj + qpiqj)
      ds(:, :) = ds + cosk * (qiqpj + qpiqj) &
         & * ((2.0_wp/g2 + 0.5_wp/alp2) * spread(vec, 1, 3)*spread(vec, 2, 3) - unity)
   end do

end subroutine get_damat_sq_rec_3d


pure function variable_info(self) result(info)
   use tblite_scf_info, only : scf_info, atom_resolved
   !> Instance of the electrostatic container
   class(damped_multipole), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(dipole=atom_resolved, quadrupole=atom_resolved)
end function variable_info

end module tblite_coulomb_multipole
