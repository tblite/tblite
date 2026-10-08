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

!> @file tblite/disp/d4.f90
!> Provides a proxy for the [DFT-D4 dispersion correction](https://dftd4.readthedocs.io)

!> Generally applicable charge-dependent London-dispersion correction, DFT-D4.
module tblite_disp_d4
   use dftd4, only : dispersion_model, d4_model, d4s_model, &
      & rational_damping_param, realspace_cutoff, &
      & new_d4_model, new_d4s_model
   use dftd4_cutoff, only : smooth_cutoff
   use dftd4_model, only : d4_qmod
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_ncoord, only : new_ncoord, ncoord_type, cn_count
   use tblite_blas, only : gemv
   use tblite_container_cache, only : container_cache
   use tblite_cutoff, only : get_lattice_points
   use tblite_disp_cache, only : dispersion_cache
   use tblite_disp_type, only : dispersion_type
   use tblite_partition, only : pair_list
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: new_d4_dispersion, new_d4s_dispersion


   !> Container for self-consistent D4 dispersion interactions
   type, public, extends(dispersion_type) :: d4_dispersion
      !> Instance of the actual D4 dispersion model
      class(dispersion_model), allocatable :: model
      !> Rational damping parameters
      type(rational_damping_param) :: param
      !> Selected real space cutoffs for this instance
      type(realspace_cutoff) :: cutoff
      !> Coordination number instance
      class(ncoord_type), allocatable :: ncoord
   contains
      !> Update dispersion cache
      procedure :: update
      !> Get information about density dependent quantities used in the energy
      procedure :: variable_info
      !> Evaluate non-selfconsistent part of the dispersion correction
      procedure :: get_engrad
      !> Evaluate selfconsistent energy of the dispersion correction
      procedure :: get_energy
      !> Evaluate charge dependent potential shift from the dispersion correction
      procedure :: get_potential
      !> Evaluate gradient contributions from the selfconsistent dispersion correction
      procedure :: get_gradient
   end type d4_dispersion

   character(len=*), parameter :: label_d4 = "self-consistent DFT-D4 dispersion"
   character(len=*), parameter :: label_d4s = "self-consistent DFT-D4S dispersion"
   real(wp), parameter :: default_disp2_width = 0.0_wp
   real(wp), parameter :: default_disp3_width = 0.0_wp


contains


!> Create a new instance of a self-consistent D4 dispersion correction
subroutine new_d4_dispersion(self, mol, s6, s8, a1, a2, s9, error, disp2_width, disp3_width)
   !> Instance of the dispersion correction
   type(d4_dispersion), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Damping parameters
   real(wp), intent(in) :: s6, s8, a1, a2, s9
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Width of smooth two-body interaction cutoff
   real(wp), intent(in), optional :: disp2_width
   !> Width of smooth three-body interaction cutoff
   real(wp), intent(in), optional :: disp3_width

   type(d4_model), allocatable :: tmp
   real(wp) :: width2, width3

   width2 = default_disp2_width
   width3 = default_disp3_width
   if (present(disp2_width)) width2 = disp2_width
   if (present(disp3_width)) width3 = disp3_width

   self%label = label_d4

   ! Create a new instance of the D4 model
   allocate(tmp)
   call new_d4_model(error, tmp, mol, qmod=d4_qmod%gfn2)
   if(allocated(error)) return
   call move_alloc(tmp, self%model)

   self%param = rational_damping_param(s6=s6, s8=s8, s9=s9, a1=a1, a2=a2)
   self%cutoff = realspace_cutoff(disp3=25.0_wp, disp2=50.0_wp, &
      & width2=width2, width3=width3)

   call new_ncoord(self%ncoord, mol, cn_count%dftd4, error, &
      & cutoff=self%cutoff%cn, rcov=self%model%rcov, en=self%model%en)
end subroutine new_d4_dispersion


!> Create a new instance of a self-consistent D4S dispersion correction
subroutine new_d4s_dispersion(self, mol, s6, s8, a1, a2, s9, error, disp2_width, disp3_width)
   !> Instance of the dispersion correction
   type(d4_dispersion), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Damping parameters
   real(wp), intent(in) :: s6, s8, a1, a2, s9
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Width of smooth two-body interaction cutoff
   real(wp), intent(in), optional :: disp2_width
   !> Width of smooth three-body interaction cutoff
   real(wp), intent(in), optional :: disp3_width

   type(d4s_model), allocatable :: tmp
   real(wp) :: width2, width3

   width2 = default_disp2_width
   width3 = default_disp3_width
   if (present(disp2_width)) width2 = disp2_width
   if (present(disp3_width)) width3 = disp3_width

   self%label = label_d4s

   ! Create a new instance of the D4S model
   allocate(tmp)
   call new_d4s_model(error, tmp, mol, qmod=d4_qmod%gfn2)
   if(allocated(error)) return
   call move_alloc(tmp, self%model)

   self%param = rational_damping_param(s6=s6, s8=s8, s9=s9, a1=a1, a2=a2)
   self%cutoff = realspace_cutoff(disp3=25.0_wp, disp2=50.0_wp, &
      & width2=width2, width3=width3)

   call new_ncoord(self%ncoord, mol, cn_count%dftd4, error, &
      & cutoff=self%cutoff%cn, rcov=self%model%rcov, en=self%model%en)
end subroutine new_d4s_dispersion


!> Update dispersion cache
subroutine update(self, mol, cache)
   !> Instance of the dispersion correction
   class(d4_dispersion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(container_cache), intent(inout) :: cache

   real(wp), allocatable :: lattr(:, :)
   type(dispersion_cache), pointer :: ptr
   integer :: mref, nlocal

   call taint(cache, ptr)
   mref = maxval(self%model%ref)

   if (.not.allocated(ptr%cn)) allocate(ptr%cn(mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%cn, lattr)
   call self%ncoord%get_coordination_number(mol, lattr, ptr%cn)
   ptr%cn_derivs_valid = .false.
   call ptr%pairs%update(self%partition, mol%nat)

   if (.not.allocated(ptr%gwvec)) allocate(ptr%gwvec(mref, mol%nat, self%model%ncoup))
   if (.not.allocated(ptr%dgwdcn)) allocate(ptr%dgwdcn(mref, mol%nat, self%model%ncoup))
   if (.not.allocated(ptr%dgwdq)) allocate(ptr%dgwdq(mref, mol%nat, self%model%ncoup))

   ! Keep both symmetric blocks of each owned pair. Energy and potential then
   ! follow the same partition without reducing the matrix itself.
   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%disp2, lattr)
   nlocal = size(ptr%pairs%neighbour)
   if (allocated(ptr%dispmat)) then
      if (any(shape(ptr%dispmat) /= [mref, mref, nlocal])) deallocate(ptr%dispmat)
   end if
   if (.not.allocated(ptr%dispmat)) allocate(ptr%dispmat(mref, mref, nlocal))
   call get_dispersion_matrix(mol, self%model, self%param, lattr, self%cutoff%disp2, &
      & self%cutoff%width2, self%model%r4r2, ptr%pairs, ptr%dispmat)
end subroutine update


!> Form the CN Jacobian only when derivatives are requested, once per geometry.
!> Keep the ncoord post-processing chain rule in the existing library routine.
subroutine ensure_cn_derivs(self, mol, ptr)
   class(d4_dispersion), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(dispersion_cache), intent(inout) :: ptr
   real(wp), allocatable :: lattr(:, :)

   if (ptr%cn_derivs_valid) return
   if (.not.allocated(ptr%dcndr)) allocate(ptr%dcndr(3, mol%nat, mol%nat))
   if (.not.allocated(ptr%dcndL)) allocate(ptr%dcndL(3, 3, mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%cn, lattr)
   call self%ncoord%get_coordination_number(mol, lattr, ptr%cn, ptr%dcndr, ptr%dcndL)
   ptr%cn_derivs_valid = .true.
end subroutine ensure_cn_derivs


!> Evaluate non-selfconsistent part of the dispersion correction
subroutine get_engrad(self, mol, cache, energies, gradient, sigma)
   !> Instance of the dispersion correction
   class(d4_dispersion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(container_cache), intent(inout) :: cache
   !> Dispersion energy
   real(wp), intent(inout) :: energies(:)
   !> Dispersion gradient
   real(wp), contiguous, intent(inout), optional :: gradient(:, :)
   !> Dispersion virial
   real(wp), contiguous, intent(inout), optional :: sigma(:, :)

   type(dispersion_cache), pointer :: ptr

   logical :: grad
   integer :: mref
   real(wp), allocatable :: qat(:)
   real(wp), allocatable :: gwvec(:, :, :), gwdcn(:, :, :), gwdq(:, :, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
   real(wp), allocatable :: dEdcn(:), dEdq(:)
   real(wp), allocatable :: lattr(:, :)

   if (abs(self%param%s9) < epsilon(1.0_wp)) return
   call view(cache, ptr)

   mref = maxval(self%model%ref)
   grad = present(gradient).and.present(sigma)
   if (grad) call ensure_cn_derivs(self, mol, ptr)

   allocate(gwvec(mref, mol%nat, self%model%ncoup), qat(mol%nat), c6(mol%nat, mol%nat))
   if (grad) then
      allocate(gwdcn(mref, mol%nat, self%model%ncoup), gwdq(mref, mol%nat, self%model%ncoup), &
         & dc6dcn(mol%nat, mol%nat), dc6dq(mol%nat, mol%nat))
   end if
   qat(:) = 0.0_wp
   if (grad) then
      allocate(dEdcn(mol%nat), dEdq(mol%nat))
      dEdcn(:) = 0.0_wp
      dEdq(:) = 0.0_wp
   end if

   call self%model%weight_references(mol, ptr%cn, qat, gwvec, gwdcn, gwdq)
   call self%model%get_atomic_c6(mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)

   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%disp3, lattr)
   call self%param%get_dispersion3(mol, lattr, self%cutoff%disp3, self%cutoff%width3, &
      & self%model%r4r2, c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient, sigma, &
      & partition=self%partition%get_d4())
   if (grad) then
      call gemv(ptr%dcndr, dEdcn, gradient, beta=1.0_wp)
      call gemv(ptr%dcndL, dEdcn, sigma, beta=1.0_wp)
   end if

end subroutine get_engrad


!> Evaluate selfconsistent energy of the dispersion correction
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the dispersion correction
   class(d4_dispersion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Dispersion energy
   real(wp), intent(inout) :: energies(:)

   type(dispersion_cache), pointer :: ptr

   call view(cache, ptr)

   call self%model%weight_references(mol, ptr%cn, wfn%qat(:, 1), ptr%gwvec)

   call contract_dispersion(self, mol, ptr%pairs, ptr%dispmat, ptr%gwvec, ptr%gwvec, &
      & energies, 0.5_wp)

end subroutine get_energy


!> Evaluate charge dependent potential shift from the dispersion correction
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the dispersion correction
   class(d4_dispersion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   type(dispersion_cache), pointer :: ptr

   call view(cache, ptr)

   call self%model%weight_references(mol, ptr%cn, wfn%qat(:, 1), ptr%gwvec, ptr%dgwdcn, ptr%dgwdq)

   call contract_dispersion(self, mol, ptr%pairs, ptr%dispmat, ptr%dgwdq, ptr%gwvec, &
      & pot%vat(:, 1), 1.0_wp)

end subroutine get_potential


!> Contract only locally owned D4/D4S blocks for energy or charge potential.
subroutine contract_dispersion(self, mol, pairs, dispmat, left, right, values, scale)
   class(d4_dispersion), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(pair_list), intent(in) :: pairs
   real(wp), intent(in) :: dispmat(:, :, :), left(:, :, :), right(:, :, :)
   real(wp), intent(inout) :: values(:)
   real(wp), intent(in) :: scale

   integer :: iat, jat, izp, jzp, iref, jref, ipair, icoup, jcoup
   real(wp) :: val, tmp

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(self, mol, pairs, dispmat, left, right, values, scale) &
   !$omp private(iat, jat, izp, jzp, iref, jref, val, tmp, ipair, icoup, jcoup)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      icoup = min(iat, self%model%ncoup)
      val = 0.0_wp
      do ipair = pairs%offset(iat-1)+1, pairs%offset(iat)
         jat = pairs%neighbour(ipair)
         jcoup = min(jat, self%model%ncoup)
         jzp = mol%id(jat)
         do jref = 1, self%model%ref(jzp)
            tmp = 0.0_wp
            do iref = 1, self%model%ref(izp)
               tmp = tmp + dispmat(iref, jref, ipair) * left(iref, iat, jcoup)
            end do
            val = val + scale * tmp * right(jref, jat, icoup)
         end do
      end do
      values(iat) = values(iat) + val
   end do
end subroutine contract_dispersion


!> Evaluate gradient contributions from the selfconsistent dispersion correction
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the dispersion correction
   class(d4_dispersion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Dispersion gradient
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Dispersion virial
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   integer :: mref
   real(wp), allocatable :: gwvec(:, :, :), gwdcn(:, :, :), gwdq(:, :, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
   real(wp), allocatable :: dEdcn(:), dEdq(:), energies(:)
   real(wp), allocatable :: lattr(:, :)
   type(dispersion_cache), pointer :: ptr

   call view(cache, ptr)
   mref = maxval(self%model%ref)
   call ensure_cn_derivs(self, mol, ptr)

   allocate(gwvec(mref, mol%nat, self%model%ncoup), gwdcn(mref, mol%nat, self%model%ncoup), &
      &  gwdq(mref, mol%nat, self%model%ncoup))
   call self%model%weight_references(mol, ptr%cn, wfn%qat(:, 1), gwvec, gwdcn, gwdq)

   allocate(c6(mol%nat, mol%nat), dc6dcn(mol%nat, mol%nat), dc6dq(mol%nat, mol%nat))
   call self%model%get_atomic_c6(mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)

   allocate(energies(mol%nat), dEdcn(mol%nat), dEdq(mol%nat))
   energies(:) = 0.0_wp
   dEdcn(:) = 0.0_wp
   dEdq(:) = 0.0_wp
   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%disp2, lattr)
   call self%param%get_dispersion2(mol, lattr, self%cutoff%disp2, self%cutoff%width2, &
      & self%model%r4r2, c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient, sigma, &
      & partition=self%partition%get_d4())
   call gemv(ptr%dcndr, dEdcn, gradient, beta=1.0_wp)
   call gemv(ptr%dcndL, dEdcn, sigma, beta=1.0_wp)
end subroutine get_gradient


subroutine get_dispersion_matrix(mol, disp, param, trans, cutoff, width, r4r2, pairs, dispmat)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Damping parameters
   type(rational_damping_param), intent(in) :: param
   !> Instance of the dispersion model
   class(dispersion_model), intent(in) :: disp
   !> Lattice points
   real(wp), intent(in) :: trans(:, :)
   !> Real space cutoff
   real(wp), intent(in) :: cutoff
   !> Width of smooth cutoff
   real(wp), intent(in) :: width
   !> Expectation values for r4 over r2 operator
   real(wp), intent(in) :: r4r2(:)
   !> Dispersion matrix
   real(wp), intent(out) :: dispmat(:, :, :)
   type(pair_list), intent(in) :: pairs
   integer :: iat, jat, izp, jzp, jtr, iref, jref, ipair, jpair
   real(wp) :: vec(3), r2, r, cutoff2, r0ij, rrij, t6, t8
   real(wp) :: edisp, dE, sw, dswdr

   dispmat = 0.0_wp
   cutoff2 = cutoff**2

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, param, disp, trans, cutoff, width, cutoff2, r4r2, pairs, dispmat) &
   !$omp private(iat, jat, izp, jzp, jtr, vec, r2, r0ij, rrij, &
   !$omp& t6, t8, edisp, dE, r, sw, dswdr, iref, jref, ipair, jpair)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do ipair = pairs%offset(iat-1)+1, pairs%offset(iat)
         jat = pairs%neighbour(ipair)
         if (jat > iat) cycle
         jzp = mol%id(jat)
         rrij = 3*r4r2(izp)*r4r2(jzp)
         r0ij = param%a1 * sqrt(rrij) + param%a2
         dE = 0.0_wp
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle
            r = sqrt(r2)
            call smooth_cutoff(r, cutoff, width, sw, dswdr)
            if (sw <= 0.0_wp) cycle

            t6 = 1.0_wp/(r2**3 + r0ij**6)
            t8 = 1.0_wp/(r2**4 + r0ij**8)

            edisp = sw * (param%s6*t6 + param%s8*rrij*t8)

            dE = dE - edisp
         end do

         jpair = pairs%find(jat, iat)
         do iref = 1, disp%ref(izp)
            do jref = 1, disp%ref(jzp)
               dispmat(iref, jref, ipair) = dE * disp%c6(iref, jref, izp, jzp)
               dispmat(jref, iref, jpair) = dE * disp%c6(jref, iref, jzp, izp)
            end do
         end do
      end do
   end do
end subroutine get_dispersion_matrix


!> Get information about density dependent quantities used in the energy
pure function variable_info(self) result(info)
   use tblite_scf_info, only : scf_info, atom_resolved
   !> Instance of the electrostatic container
   class(d4_dispersion), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=atom_resolved)
end function variable_info


!> Inspect container cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(dispersion_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(dispersion_cache), allocatable :: tmp
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
   type(dispersion_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(dispersion_cache)
      ptr => target
   end select
end subroutine view

end module tblite_disp_d4
