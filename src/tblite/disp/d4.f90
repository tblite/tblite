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
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use dftd4, only : d4_model, damping_param, rational_damping_param, realspace_cutoff, &
      & get_coordination_number, new_d4_model
   use dftd4_model, only : d4_ref
   use dftd4_charge, only : get_eeq_charges => get_charges
   use tblite_blas, only : dot, gemv
   use tblite_container_cache, only : container_cache
   use tblite_disp_cache, only : dispersion_cache
   use tblite_disp_type, only : dispersion_type
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_cutoff, only : get_lattice_points
   implicit none
   private

   public :: new_d4_dispersion, get_eeq_charges


   !> Container for self-consistent D4 dispersion interactions
   type, public, extends(dispersion_type) :: d4_dispersion
      !> Instance of the actual D4 dispersion model
      type(d4_model) :: model
      !> Rational damping parameters
      type(rational_damping_param) :: param
      !> Selected real space cutoffs for this instance
      type(realspace_cutoff) :: cutoff
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

   character(len=*), parameter :: label = "self-consistent DFT-D4 dispersion"


contains


!> Create a new instance of a self-consistent D4 dispersion correction
subroutine new_d4_dispersion(self, mol, s6, s8, a1, a2, s9)
   !> Instance of the dispersion correction
   type(d4_dispersion), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: s6, s8, a1, a2, s9

   self%label = label
   call new_d4_model(self%model, mol, ref=d4_ref%gfn2)
   self%param = rational_damping_param(s6=s6, s8=s8, s9=s9, a1=a1, a2=a2)
   self%cutoff = realspace_cutoff(disp3=25.0_wp, disp2=50.0_wp)
end subroutine new_d4_dispersion


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
   integer :: mref

   call taint(cache, ptr)
   mref = maxval(self%model%ref)

   if (.not.allocated(ptr%cn)) allocate(ptr%cn(mol%nat))
   if (.not.allocated(ptr%dcndr)) allocate(ptr%dcndr(3, mol%nat, mol%nat))
   if (.not.allocated(ptr%dcndL)) allocate(ptr%dcndL(3, 3, mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, self%cutoff%cn, self%model%rcov, self%model%en, &
      & ptr%cn, ptr%dcndr, ptr%dcndL)

   if (.not.allocated(ptr%dispmat)) allocate(ptr%dispmat(mref, mol%nat, mref, mol%nat))
   if (.not.allocated(ptr%vvec)) allocate(ptr%vvec(mref, mol%nat))
   if (.not.allocated(ptr%gwvec)) allocate(ptr%gwvec(mref, mol%nat))
   if (.not.allocated(ptr%dgwdq)) allocate(ptr%dgwdq(mref, mol%nat))

   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%disp2, lattr)
   call get_dispersion_matrix(mol, self%model, self%param, lattr, self%cutoff%disp2, &
      & self%model%r4r2, ptr%dispmat)
end subroutine update


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

   call view(cache, ptr)

   call get_dispersion_nonsc(mol, self%model, self%param, self%cutoff, ptr, &
      & energies, gradient, sigma)
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

   call gemv(ptr%dispmat, ptr%gwvec, ptr%vvec, alpha=0.5_wp)
   ptr%vvec(:, :) = ptr%vvec * ptr%gwvec
   energies(:) = energies + sum(ptr%vvec, 1)
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

   call self%model%weight_references(mol, ptr%cn, wfn%qat(:, 1), ptr%gwvec, ptr%vvec, &
      & ptr%dgwdq)

   call gemv(ptr%dispmat, ptr%gwvec, ptr%vvec)
   ptr%vvec(:, :) = ptr%vvec * ptr%dgwdq
   pot%vat(:, 1) = pot%vat(:, 1) + sum(ptr%vvec, 1)
end subroutine get_potential


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
   real(wp), allocatable :: gwvec(:, :), gwdcn(:, :), gwdq(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
   real(wp), allocatable :: dEdcn(:), dEdq(:), energies(:)
   real(wp), allocatable :: lattr(:, :)
   type(dispersion_cache), pointer :: ptr

   call view(cache, ptr)
   mref = maxval(self%model%ref)

   allocate(gwvec(mref, mol%nat), gwdcn(mref, mol%nat), gwdq(mref, mol%nat))
   call self%model%weight_references(mol, ptr%cn, wfn%qat(:, 1), gwvec, gwdcn, gwdq)

   allocate(c6(mol%nat, mol%nat), dc6dcn(mol%nat, mol%nat), dc6dq(mol%nat, mol%nat))
   call self%model%get_atomic_c6(mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)

   allocate(energies(mol%nat), dEdcn(mol%nat), dEdq(mol%nat))
   energies(:) = 0.0_wp
   dEdcn(:) = 0.0_wp
   dEdq(:) = 0.0_wp
   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%disp2, lattr)
   call self%param%get_dispersion2(mol, lattr, self%cutoff%disp2, self%model%r4r2, &
      & c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient, sigma)
   call gemv(ptr%dcndr, dEdcn, gradient, beta=1.0_wp)
   call gemv(ptr%dcndL, dEdcn, sigma, beta=1.0_wp)
end subroutine get_gradient


subroutine get_dispersion_matrix(mol, disp, param, trans, cutoff, r4r2, dispmat)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Damping parameters
   type(rational_damping_param), intent(in) :: param
   !> Instance of the dispersion model
   type(d4_model), intent(in) :: disp
   !> Lattice points
   real(wp), intent(in) :: trans(:, :)
   !> Real space cutoff
   real(wp), intent(in) :: cutoff
   !> Expectation values for r4 over r2 operator
   real(wp), intent(in) :: r4r2(:)
   !> Dispersion matrix
   real(wp), intent(out) :: dispmat(:, :, :, :)

   integer :: iat, jat, izp, jzp, jtr, iref, jref
   real(wp) :: vec(3), r2, cutoff2, r0ij, rrij, t6, t8, edisp, dE

   dispmat(:, :, :, :) = 0.0_wp
   cutoff2 = cutoff**2

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, param, disp, trans, cutoff2, r4r2, dispmat) &
   !$omp private(iat, jat, izp, jzp, jtr, vec, r2, r0ij, rrij, &
   !$omp& t6, t8, edisp, dE)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         rrij = 3*r4r2(izp)*r4r2(jzp)
         r0ij = param%a1 * sqrt(rrij) + param%a2
         dE = 0.0_wp
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle

            t6 = 1.0_wp/(r2**3 + r0ij**6)
            t8 = 1.0_wp/(r2**4 + r0ij**8)

            edisp = param%s6*t6 + param%s8*rrij*t8

            dE = dE - edisp
         end do

         do iref = 1, disp%ref(izp)
            do jref = 1, disp%ref(jzp)
               dispmat(iref, iat, jref, jat) = dE * disp%c6(iref, jref, izp, jzp)
               dispmat(jref, jat, iref, iat) = dE * disp%c6(jref, iref, jzp, izp)
            end do
         end do
      end do
   end do
end subroutine get_dispersion_matrix


!> Wrapper to handle the evaluation of dispersion energy and derivatives
subroutine get_dispersion_nonsc(mol, disp, param, cutoff, cache, energies, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Dispersion model
   type(d4_model), intent(in) :: disp
   !> Damping parameters
   class(damping_param), intent(in) :: param
   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff
   !> Cached data between different dispersion runs
   type(dispersion_cache), intent(inout) :: cache
   !> Dispersion energy
   real(wp), intent(inout) :: energies(:)
   !> Dispersion gradient
   real(wp), intent(inout), contiguous, optional :: gradient(:, :)
   !> Dispersion virial
   real(wp), intent(inout), contiguous, optional :: sigma(:, :)

   logical :: grad
   integer :: mref
   real(wp), allocatable :: qat(:)
   real(wp), allocatable :: gwvec(:, :), gwdcn(:, :), gwdq(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
   real(wp), allocatable :: dEdcn(:), dEdq(:)
   real(wp), allocatable :: lattr(:, :)

   mref = maxval(disp%ref)
   grad = present(gradient).and.present(sigma)

   allocate(gwvec(mref, mol%nat), qat(mol%nat), c6(mol%nat, mol%nat))
   if (grad) then
      allocate(gwdcn(mref, mol%nat), gwdq(mref, mol%nat), &
         & dc6dcn(mol%nat, mol%nat), dc6dq(mol%nat, mol%nat))
   end if
   qat(:) = 0.0_wp
   if (grad) then
      allocate(dEdcn(mol%nat), dEdq(mol%nat))
      dEdcn(:) = 0.0_wp
      dEdq(:) = 0.0_wp
   end if

   call disp%weight_references(mol, cache%cn, qat, gwvec, gwdcn, gwdq)
   call disp%get_atomic_c6(mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
   call param%get_dispersion3(mol, lattr, cutoff%disp3, disp%r4r2, &
      & c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient, sigma)
   if (grad) then
      call gemv(cache%dcndr, dEdcn, gradient, beta=1.0_wp)
      call gemv(cache%dcndL, dEdcn, sigma, beta=1.0_wp)
   end if

end subroutine get_dispersion_nonsc


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
