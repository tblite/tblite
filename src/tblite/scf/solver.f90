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

!> @file tblite/scf/solver.f90
!> Provides a base class for defining electronic solvers

!> Declaration of the abstract base class for electronic solvers
module tblite_scf_solver
   use mctc_env, only : sp, dp, error_type, wp 
   use tblite_blas, only : gemm
   use tblite_integral_type, only : integral_type
   use tblite_wavefunction_fermi, only : get_fermi_filling
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   !> Abstract base class for electronic solvers
   type, public, abstract :: solver_type
   contains
      !> Solve the eigenvalue problem for the given Hamiltonian and overlap matrix
      generic :: solve => solve_sp, solve_dp
      procedure(solve_sp), deferred :: solve_sp
      procedure(solve_dp), deferred :: solve_dp
      !> Delete the inner components of a solver instance
      procedure :: delete
      !> Build the denisty matrix 
      procedure :: get_density
      !> Build density matrix, based on coefficients and occupation numbers
      procedure :: get_density_matrix
      !> Build energy weighted density matrix
      procedure :: get_energy_w_density_matrix
      !> Reset internal data of the solver
      procedure :: reset
   end type solver_type

   abstract interface
      subroutine solve_sp(self, hmat, smat, eval, error)
         import :: solver_type, error_type, sp
         class(solver_type), intent(inout) :: self
         real(sp), contiguous, intent(inout) :: hmat(:, :)
         real(sp), contiguous, intent(in) :: smat(:, :)
         real(sp), contiguous, intent(inout) :: eval(:)
         type(error_type), allocatable, intent(out) :: error
      end subroutine solve_sp
      subroutine solve_dp(self, hmat, smat, eval, error)
         import :: solver_type, error_type, dp
         class(solver_type), intent(inout) :: self
         real(dp), contiguous, intent(inout) :: hmat(:, :)
         real(dp), contiguous, intent(in) :: smat(:, :)
         real(dp), contiguous, intent(inout) :: eval(:)
         type(error_type), allocatable, intent(out) :: error
      end subroutine solve_dp
   end interface

contains

!> Delete the solver instance
subroutine delete(self)
   !> Instance of the solver
   class(solver_type), intent(inout) :: self
end subroutine

!> Reset the solver instance, in case of a geometry change for example
subroutine reset(self)
   class(solver_type) :: self
end subroutine

subroutine get_density(self, wfn, ints, ts, error)
   !> Solver for the general eigenvalue problem
   class(solver_type), intent(inout) :: self
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Integral container
   type(integral_type), intent(in) :: ints
   !> Electronic entropy
   real(wp), intent(out) :: ts
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp) :: e_fermi, stmp(2)
   real(wp), allocatable :: focc(:)
   integer :: spin

   select case(wfn%nspin)
   case default
      call self%solve(wfn%coeff(:, :, 1), ints%overlap, wfn%emo(:, 1), error)
      if (allocated(error)) return

      allocate(focc(size(wfn%focc, 1)))
      wfn%focc(:, :) = 0.0_wp
      do spin = 1, 2
         call get_fermi_filling(wfn%nel(spin), wfn%kt, wfn%emo(:, 1), &
            & wfn%homo(spin), focc, e_fermi)
         call get_electronic_entropy(focc, wfn%kt, stmp(spin))
         wfn%focc(:, 1) = wfn%focc(:, 1) + focc
      end do
      ts = sum(stmp)

      call self%get_density_matrix(wfn%focc(:, 1), wfn%coeff(:, :, 1), wfn%density(:, :, 1))
   case(2)
      wfn%coeff = 2*wfn%coeff
      do spin = 1, 2
         call self%solve(wfn%coeff(:, :, spin), ints%overlap, wfn%emo(:, spin), error)
         if (allocated(error)) return

         call get_fermi_filling(wfn%nel(spin), wfn%kt, wfn%emo(:, spin), &
            & wfn%homo(spin), wfn%focc(:, spin), e_fermi)
         call get_electronic_entropy(wfn%focc(:, spin), wfn%kt, stmp(spin))
         call self%get_density_matrix(wfn%focc(:, spin), wfn%coeff(:, :, spin), &
            & wfn%density(:, :, spin))
      end do
      ts = sum(stmp)
   end select
end subroutine get_density

subroutine get_electronic_entropy(occ, kt, s)
   real(wp), intent(in) :: occ(:)
   real(wp), intent(in) :: kt
   real(wp), intent(out) :: s

   s = sum(log(occ ** occ * (1 - occ) ** (1 - occ))) * kt
end subroutine get_electronic_entropy

subroutine get_density_matrix(self, focc, coeff, pmat)
   class(solver_type) :: self
   real(wp), intent(in) :: focc(:)
   real(wp), contiguous, intent(in) :: coeff(:, :)
   real(wp), contiguous, intent(out) :: pmat(:, :)

   real(wp), allocatable :: scratch(:, :)
   integer :: iao, jao
   allocate(scratch(size(pmat, 1), size(pmat, 2)))
   !$omp parallel do collapse(2) default(none) schedule(runtime) &
   !$omp shared(scratch, coeff, focc, pmat) private(iao, jao)
   do iao = 1, size(pmat, 1)
      do jao = 1, size(pmat, 2)
         scratch(jao, iao) = coeff(jao, iao) * focc(iao)
      end do
   end do
   call gemm(scratch, coeff, pmat, transb='t')
end subroutine get_density_matrix

subroutine get_energy_w_density_matrix(self, wfn, wdensity)
   class(solver_type) :: self
   type(wavefunction_type), intent(inout) :: wfn
   real(wp) :: wdensity(:,:,:)
   real(wp), allocatable :: tmp(:)
   integer :: spin

   do spin = 1, wfn%nspin
      tmp = wfn%focc(:, spin) * wfn%emo(:, spin)
      call self%get_density_matrix(tmp, wfn%coeff(:, :, spin), wdensity(:, :, spin))
   end do
end subroutine

end module tblite_scf_solver
