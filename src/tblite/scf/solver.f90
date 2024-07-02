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
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   !> Abstract base class for electronic solvers
   type, public, abstract :: solver_type
   contains
      generic :: solve => solve_sp, solve_dp
      procedure(solve_sp), deferred :: solve_sp
      procedure(solve_dp), deferred :: solve_dp
      procedure :: delete
      procedure :: get_density_matrix
      procedure :: get_energy_w_density_matrix
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

subroutine delete(self)
   class(solver_type) :: self
end subroutine

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
