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

!> @file tblite/scf/diag.f90
!> Provides a base class for defining diagonalization based electronic solvers

!> Declaration of the abstract base class for electronic solvers based on diagonalization
module tblite_scf_diag
   use mctc_env, only : sp, dp, wp, error_type
   use tblite_blas, only : syrk
   use tblite_scf_solver, only : solver_type
   use tblite_wavefunction_fermi, only : get_fermi_filling
   implicit none
   private

   !> Abstract base class for electronic solvers
   type, public, abstract, extends(solver_type) :: diag_solver_type
      !> Occupation-weighted coefficients reused for density construction
      real(wp), allocatable, private :: scratch(:, :)
   contains
      generic :: solve => solve_sp, solve_dp
      procedure(solve_sp), deferred :: solve_sp
      procedure(solve_dp), deferred :: solve_dp
      procedure :: get_density
      procedure :: get_wdensity
      procedure :: delete => delete_diag_solver
      !> Build a density matrix in the storage layout supported by the solver
      procedure :: get_density_matrix
   end type diag_solver_type

   abstract interface
      subroutine solve_sp(self, hmat, smat, eval, error)
         import :: diag_solver_type, error_type, sp
         class(diag_solver_type), intent(inout) :: self
         real(sp), contiguous, intent(inout) :: hmat(:, :)
         real(sp), contiguous, intent(in) :: smat(:, :)
         real(sp), contiguous, intent(inout) :: eval(:)
         type(error_type), allocatable, intent(out) :: error
      end subroutine solve_sp
      subroutine solve_dp(self, hmat, smat, eval, error)
         import :: diag_solver_type, error_type, dp
         class(diag_solver_type), intent(inout) :: self
         real(dp), contiguous, intent(inout) :: hmat(:, :)
         real(dp), contiguous, intent(in) :: smat(:, :)
         real(dp), contiguous, intent(inout) :: eval(:)
         type(error_type), allocatable, intent(out) :: error
      end subroutine solve_dp
   end interface


contains

subroutine get_density(self, hmat, smat, eval, focc, density, error)
   !> Solver for the general eigenvalue problem
   class(diag_solver_type), intent(inout) :: self
   !> Overlap matrix
   real(wp), contiguous, intent(in) :: smat(:, :)
   !> Hamiltonian matrix, contains eigenvectors on output
   real(wp), contiguous, intent(inout) :: hmat(:, :, :)
   !> Eigenvalues
   real(wp), contiguous, intent(inout) :: eval(:, :)
   !> Occupation numbers
   real(wp), contiguous, intent(inout) :: focc(:, :)
   !> Density matrix
   real(wp), contiguous, intent(inout) :: density(:, :, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp) :: e_fermi
   integer :: nspin, spin, homo

   nspin = min(size(eval, 2), size(hmat, 3), size(density, 3))

   select case(nspin)
   case default
      call self%solve(hmat(:, :, 1), smat, eval(:, 1), error)
      if (allocated(error)) return

      focc(:, :) = 0.0_wp
      do spin = 1, 2
         call get_fermi_filling(self%nel(spin), self%kt, eval(:, 1), &
            & homo, focc(:, spin), e_fermi)
      end do

      focc(:, 1) = focc(:, 1) + focc(:, 2)
      call self%get_density_matrix(focc(:, 1), hmat(:, :, 1), density(:, :, 1), error)
      focc(:, 1) = focc(:, 1) - focc(:, 2)
   case(2)
      hmat(:, :, :) = 2*hmat
      do spin = 1, 2
         call self%solve(hmat(:, :, spin), smat, eval(:, spin), error)
         if (allocated(error)) return

         call get_fermi_filling(self%nel(spin), self%kt, eval(:, spin), &
            & homo, focc(:, spin), e_fermi)
         call self%get_density_matrix(focc(:, spin), hmat(:, :, spin), density(:, :, spin), error)
         if (allocated(error)) return
      end do
   end select
end subroutine get_density

subroutine get_wdensity(self, hmat, smat, eval, focc, density, error)
   !> Solver for the general eigenvalue problem
   class(diag_solver_type), intent(inout) :: self
   !> Overlap matrix
   real(wp), contiguous, intent(in) :: smat(:, :)
   !> Hamiltonian matrix containing eigenvectors from SCF
   real(wp), contiguous, intent(inout) :: hmat(:, :, :)
   !> Eigenvalues
   real(wp), contiguous, intent(inout) :: eval(:, :)
   !> Occupation numbers
   real(wp), contiguous, intent(inout) :: focc(:, :)
   !> Density matrix
   real(wp), contiguous, intent(inout) :: density(:, :, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: spin, nspin
   real(wp), allocatable :: tmp(:)

   nspin = min(size(eval, 2), size(hmat, 3), size(density, 3))

   if (nspin == 1 .and. size(focc, 2) == 2) then
      focc(:, 1) = focc(:, 1) + focc(:, 2)
   end if

   do spin = 1, nspin
      tmp = focc(:, spin) * eval(:, spin)
      call self%get_density_matrix(tmp, hmat(:, :, spin), density(:, :, spin), error)
      if (allocated(error)) exit
   end do

   if (nspin == 1 .and. size(focc, 2) == 2) then
      focc(:, 1) = focc(:, 1) - focc(:, 2)
   end if
end subroutine get_wdensity

!> Get the density matrix from the coefficients and occupation numbers
subroutine get_density_matrix(self, focc, coeff, pmat, error)
   !> Solver and reusable density workspace
   class(diag_solver_type), intent(inout) :: self
   !> Occupation numbers
   real(wp), intent(in) :: focc(:)
   !> Coefficients of the wavefunction
   real(wp), contiguous, intent(in) :: coeff(:, :)
   !> Density matrix to be computed
   real(wp), contiguous, intent(inout) :: pmat(:, :)
   !> Error handling for solver-specific density construction
   type(error_type), allocatable, intent(inout) :: error

   integer :: nao, nmo, iao, jao, imo, k, nk, isign
   integer :: idx(min(size(coeff, 2), size(focc)))
   real(wp) :: beta, scale, sgn

   if (allocated(self%scratch)) then
      if (any(shape(self%scratch) /= shape(coeff))) deallocate(self%scratch)
   end if
   if (.not.allocated(self%scratch)) allocate(self%scratch(size(coeff, 1), size(coeff, 2)))

   nao = size(coeff, 1)
   nmo = min(size(coeff, 2), size(focc))
   if (all(focc(:nmo) == 0.0_wp)) then
      pmat(:, :) = 0.0_wp
      return
   end if

   ! Upper-triangle SYRK over occupied orbitals (following g-xTB)
   beta = 0.0_wp
   do isign = 1, -1, -2
      sgn = real(isign, wp)
      nk = 0
      do imo = 1, nmo
         if (sgn * focc(imo) > 0.0_wp) then
            nk = nk + 1
            idx(nk) = imo
         end if
      end do
      if (nk == 0) cycle
      !$omp parallel do default(none) schedule(static) &
      !$omp shared(self, coeff, focc, idx, nao, nk, sgn) private(k, imo, iao, scale)
      do k = 1, nk
         imo = idx(k)
         scale = sqrt(sgn * focc(imo))
         do iao = 1, nao
            self%scratch(iao, k) = coeff(iao, imo) * scale
         end do
      end do
      call syrk(self%scratch, pmat, uplo="u", trans="n", alpha=sgn, beta=beta, n=nao, k=nk)
      beta = 1.0_wp
   end do

   !$omp parallel do default(none) schedule(static, 1) &
   !$omp shared(pmat, nao) private(iao, jao)
   do jao = 1, nao
      do iao = jao + 1, nao
         pmat(iao, jao) = pmat(jao, iao)
      end do
   end do
end subroutine get_density_matrix

!> Delete the solver instance
subroutine delete_diag_solver(self)
   !> Solver for the general eigenvalue problem
   class(diag_solver_type), intent(inout) :: self

   if (allocated(self%scratch)) deallocate(self%scratch)
end subroutine delete_diag_solver

end module tblite_scf_diag
