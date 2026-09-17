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

!> @file tblite/wavefunction/localization/type.f90
!> Provides an abstract base class for orbital localization methods

!> Declaration of the abstract orbital localization method
module tblite_wavefunction_localization_type
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_blas, only : gemm, dot
   use tblite_lapack, only : geqp3
   use tblite_wavefunction_localization_jacobi, only : jacobi_type
   implicit none
   private

   public :: localization_type
   public :: get_common_origin_position, get_orbital_centers
   public :: localization_method, get_localization_id

   !> Abstract base class for wavefunction orbital localization methods
   type, abstract :: localization_type
      !> Optimizer used to solve the localization problem
      type(jacobi_type) :: optimizer
   contains
      !> Construct a pivoted projected-AO guess for the occupied-space rotation
      procedure :: guess
      !> Localize the occupied molecular orbitals of a wavefunction
      procedure :: localize
      !> Prepare the method-specific localization problem for one spin channel
      procedure(prepare), deferred :: prepare
      !> Postprocess the optimized transformation regarding order
      procedure(postprocess), deferred :: postprocess
   end type localization_type

   abstract interface
      !> Prepare the localization objective to be maximized by the optimizer.
      subroutine prepare(self, mol, bas, overlap, dipole, coeff_occ, opmat, error)
         import :: localization_type, structure_type, basis_type, wp, error_type
         !> Instance of the localization method
         class(localization_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Basis set information
         type(basis_type), intent(in) :: bas
         !> Overlap integrals
         real(wp), intent(in) :: overlap(:, :)
         !> Dipole integrals with moment operator centered on last index
         real(wp), intent(in) :: dipole(:, :, :)
         !> Guess occupied orbital coefficients
         real(wp), intent(in) :: coeff_occ(:, :)
         !> Occupied-space operator matrices whose squared diagonal is maximized
         real(wp), allocatable, intent(out) :: opmat(:, :, :)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine prepare

      !> Sort the localized orbitals according to their localized (rotated) orbital energy
      subroutine postprocess(self, emo_occ, trafo)
         import :: localization_type, wp
         !> Instance of the localization method
         class(localization_type), intent(in) :: self
         !> Canonical occupied orbital energies
         real(wp), intent(in) :: emo_occ(:)
         !> Occupied-space orthogonal transformation
         real(wp), intent(inout) :: trafo(:, :)
      end subroutine postprocess
   end interface

   !> Possible orbital localization methods
   type :: enum_localization_method
      !> Foster-Boys localization
      integer :: fosterboys = 1
   end type enum_localization_method

   !> Actual enumerator for orbital localization methods
   type(enum_localization_method), parameter :: localization_method = enum_localization_method()

contains


!> Localize occupied block of canonical orbital coefficients for each spin channel
subroutine localize(self, mol, bas, overlap, dipole, coeff, emo, nel, accuracy, &
   & coeff_local, converged, error)
   !> Instance of the localization method
   class(localization_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Overlap integral
   real(wp), intent(in) :: overlap(:, :)
   !> Dipole integrals with moment operator centered on last index
   real(wp), intent(in) :: dipole(:, :, :)
   !> Canonical orbital coefficients
   real(wp), intent(in) :: coeff(:, :, :)
   !> Canonical orbital energies
   real(wp), intent(in) :: emo(:, :)
   !> Number of electrons in each spin channel
   real(wp), intent(in) :: nel(:)
   !> Accuracy setting of the calculation
   real(wp), intent(in) :: accuracy
   !> Localized orbital coefficients
   real(wp), intent(out) :: coeff_local(:, :, :)
   !> Whether the localization converged within the allowed number of iterations
   logical, intent(out) :: converged
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), allocatable :: opmat(:, :, :), trafo(:, :), coeff_guess(:, :)
   logical :: spin_converged
   integer :: nao, nspin, spin, nocc

   nao = size(coeff, 1)
   nspin = size(coeff, 3)

   converged = .true.
   coeff_local(:, :, :) = coeff

   ! Localize each spin channel independently in its own occupied subspace
   do spin = 1, nspin
      nocc = nint(nel(spin))
      if (nocc < 2) cycle

      ! Seed the occupied-space rotation with an atom-centered guess so the
      ! optimizer only has to refine an already localized starting point
      allocate(trafo(nocc, nocc))
      call self%guess(coeff(:, :nocc, spin), overlap, trafo, error)
      if (allocated(error)) then
         converged = .false.
         return
      end if

      ! Transform the canonical occupied orbitals into the guess localized orbitals
      allocate(coeff_guess(nao, nocc))
      call gemm(coeff(:, :nocc, spin), trafo, coeff_guess)

      ! Build localization criterion operator matrix relative to the guess
      call self%prepare(mol, bas, overlap, dipole, coeff_guess, opmat, error)
      if (allocated(error)) then
         converged = .false.
         return
      end if

      ! Refine the existing guess transformation via Jacobi sweeps
      call self%optimizer%optimize(opmat, accuracy, trafo, spin_converged)
      if (allocated(error)) then
         converged = .false.
         return
      end if

      ! Postprocess the transformation regarding the order of the localized orbitals
      call self%postprocess(emo(:nocc, spin), trafo)

      ! Apply localization transformation to occupied orbital coefficients
      call gemm(coeff(:, :nocc, spin), trafo, coeff_local(:, :nocc, spin))
      converged = converged .and. spin_converged
      deallocate(trafo, coeff_guess, opmat)
   end do

end subroutine localize


!> Construct a pivoted projected-AO guess for the occupied-space rotation,
!> by QR factorization of the MOs projected onto the atomic orbital basis.
subroutine guess(self, coeff_occ, overlap, guess_trafo, error)
   !> Instance of the localization method
   class(localization_type), intent(in) :: self
   !> Canonical occupied orbital coefficients
   real(wp), intent(in) :: coeff_occ(:, :)
   !> Overlap integrals
   real(wp), intent(in) :: overlap(:, :)
   !> Orthogonal occupied-space guess transformation
   real(wp), intent(out) :: guess_trafo(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: nao, nocc, info, i
   integer, allocatable :: jpvt(:)
   real(wp), allocatable :: projected(:, :), rmat(:, :)
   real(wp) :: rmax, rmin

   nao = size(coeff_occ, 1)
   nocc = size(coeff_occ, 2)

   allocate(projected(nocc, nao), jpvt(nao), rmat(nocc, nao))

   ! Project the AOs into the occupied MO basis
   call gemm(coeff_occ, overlap, projected, transa="T")

   ! Perform column-pivoted QR factorization of the projected AOs
   ! for the ordered nocc most linearly independent AOs in the occupied MO basis.
   call geqp3(projected, jpvt, guess_trafo, rmat, info)
   if (info /= 0) then
      call fatal_error(error, "Pivoted QR factorization of projected AOs failed")
      return
   end if

   ! Fallback to identity if the rank of the projected AO matrix is too low
   rmax = maxval([(abs(rmat(i, i)), i = 1, nocc)])
   rmin = minval([(abs(rmat(i, i)), i = 1, nocc)])
   if (rmax <= 0.0_wp .or. rmin <= epsilon(1.0_wp) * rmax) then
      guess_trafo(:, :) = 0.0_wp
      do i = 1, nocc
         guess_trafo(i, i) = 1.0_wp
      end do
      return
   end if

   ! Remove the column-sign ambiguity by making r non-negative
   do i = 1, nocc
      if (rmat(i, i) < 0.0_wp) guess_trafo(:, i) = -guess_trafo(:, i)
   end do

end subroutine guess


!> Construct the Cartesian position operators in a common origin
subroutine get_common_origin_position(mol, bas, overlap, dipole, position, error)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Overlap integrals
   real(wp), intent(in) :: overlap(:, :)
   !> Dipole integrals with moment operator centered on last index
   real(wp), intent(in) :: dipole(:, :, :)
   !> Common-origin Cartesian position operators
   real(wp), allocatable, intent(out) :: position(:, :, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: nao, iao, jao, k

   nao = size(overlap, 1)
   allocate(position(nao, nao, 3))

   ! Shift each AO-centered dipole integral to the common origin
   do k = 1, 3
      do jao = 1, nao
         do iao = 1, nao
            position(iao, jao, k) = dipole(k, iao, jao) &
               & + mol%xyz(k, bas%ao2at(jao)) * overlap(iao, jao)
         end do
      end do
      ! Symmetrize common-origin position operator
      position(:, :, k) = 0.5_wp * (position(:, :, k) + transpose(position(:, :, k)))
   end do

end subroutine get_common_origin_position


!> Compute charge center position of each occupied orbital of a given coefficient set
subroutine get_orbital_centers(mol, bas, overlap, dipole, coeff, nel, centers, error)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Overlap integrals
   real(wp), intent(in) :: overlap(:, :)
   !> Dipole integrals with moment operator centered on last index
   real(wp), intent(in) :: dipole(:, :, :)
   !> Orbital coefficients
   real(wp), intent(in) :: coeff(:, :, :)
   !> Number of electrons in each spin channel
   real(wp), intent(in) :: nel(:)
   !> Orbital centers with virtual columns left at zero
   real(wp), intent(out) :: centers(:, :, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), allocatable :: position(:, :, :), proj(:, :)
   integer :: nao, nspin, spin, nocc, i, k

   nao = size(coeff, 1)
   nspin = size(coeff, 3)

   centers(:, :, :) = 0.0_wp

   call get_common_origin_position(mol, bas, overlap, dipole, position, error)
   if (allocated(error)) return

   allocate(proj(nao, nao))
   do spin = 1, nspin
      nocc = nint(nel(spin))
      if (nocc < 1) cycle
      ! The orbital center is the diagonal of the occupied-space position operator
      do k = 1, 3
         call gemm(position(:, :, k), coeff(:, :nocc, spin), proj(:, :nocc))
         do i = 1, nocc
            centers(k, i, spin) = dot(coeff(:, i, spin), proj(:, i))
         end do
      end do
   end do

end subroutine get_orbital_centers


!> Translate a localization method name to its integer identifier
subroutine get_localization_id(name, id, error)
   !> Name of the localization method
   character(len=*), intent(in) :: name
   !> Integer identifier of the localization method
   integer, intent(out) :: id
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   id = -1

   select case(trim(name))
   case("foster-boys")
      id = localization_method%fosterboys
   case default
      call fatal_error(error, "Unknown orbital localization method '"//trim(name)//"'")
   end select

end subroutine get_localization_id

end module tblite_wavefunction_localization_type
