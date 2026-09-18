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

!> @file tblite/wavefunction/localization/jacobi.f90
!> Implements the Jacobi pair-rotation orbital-rotation optimizer

!> Jacobi pair-rotation optimizer repeatedly rotates pairs of orbitals to maximize
!> the sum of squared diagonal elements of the stacked operators
module tblite_wavefunction_localization_jacobi
   use mctc_env, only : wp, error_type, fatal_error
   implicit none
   private

   public :: jacobi_type, new_jacobi_optimizer

   !> Jacobi pair-rotation optimizer
   type :: jacobi_type
      !> Maximum number of Jacobi sweeps over all orbital pairs
      integer :: max_sweeps
   contains
      !> Optimize the orbital rotation by repeated pairwise Jacobi rotations
      procedure :: optimize
   end type jacobi_type

   !> Default convergence accuracy for optimizer
   real(wp), parameter :: default_accuracy = 1.0e-6_wp

   !> Initial screening threshold for skipping small pair rotations
   real(wp), parameter :: initial_screen = 0.25_wp

   !> Maximum number of Jacobi sweeps over all orbital pairs
   integer, parameter :: default_max_sweeps = 20000

   !> Minimum total pairs below which the parallel region runs single-threaded
   integer, parameter :: omp_min_work = 256

contains


!> Create a new Jacobi pair-rotation optimizer
subroutine new_jacobi_optimizer(self, max_sweeps)
   !> Instance of the Jacobi optimizer
   type(jacobi_type), intent(out) :: self
   !> Maximum number of Jacobi sweeps over all orbital pairs
   integer, intent(in), optional :: max_sweeps

   if (present(max_sweeps)) then
      self%max_sweeps = max_sweeps
   else
      self%max_sweeps = default_max_sweeps
   end if

end subroutine new_jacobi_optimizer


!> Build a round-robin schedule of disjoint occupied-orbital pairs
subroutine build_round_robin_schedule(nocc, pairs)
   !> Number of occupied orbitals
   integer, intent(in) :: nocc
   !> Orbital indices in every pair, zero denotes a dummy pair
   integer, allocatable, intent(out) :: pairs(:, :, :)

   integer :: nall, iround, ipair, tmp
   integer, allocatable :: players(:)

   ! No orbital pair exists for fewer than two orbitals.
   if (nocc < 2) then
      allocate(pairs(2, 0, 0))
      return
   end if

   ! Ensure an even number of participants for the round-robin construction.
   nall = nocc + mod(nocc, 2)

   ! There are nall-1 rounds containing nall/2 disjoint pairs each.
   allocate(pairs(2, nall/2, nall-1), players(nall))

   ! Initialize the real orbital indices.
   players(:nocc) = [(ipair, ipair = 1, nocc)]

   ! Add the dummy player for odd nocc.
   if (nall > nocc) players(nall) = 0

   do iround = 1, nall-1
      do ipair = 1, nall/2
         ! Pair a player with a corresponding player counted from the end.
         if (players(ipair) == 0 .or. players(nall-ipair+1) == 0) then
            ! A pair involving the dummy player performs no rotation.
            pairs(:, ipair, iround) = 0
         else
            ! Store the larger orbital index first
            pairs(1, ipair, iround) = max(players(ipair), players(nall-ipair+1))
            pairs(2, ipair, iround) = min(players(ipair), players(nall-ipair+1))
         end if
      end do
      ! Keep the first player fixed and rotate remaining players by one position.
      tmp = players(nall)
      players(3:nall) = players(2:nall-1)
      players(2) = tmp
   end do
end subroutine build_round_robin_schedule


!> Jacobi pair-rotation optimizer maximizes sum of squared diagonal elements
!> of the stacked operators by repeated pairwise rotations
subroutine optimize(self, opmat, accuracy, trafo, error)
   !> Instance of the Jacobi optimizer
   class(jacobi_type), intent(in) :: self
   !> Pair operator matrices to maximize
   real(wp), intent(inout) :: opmat(:, :, :)
   !> Accuracy setting of the calculation
   real(wp), intent(in) :: accuracy
   !> Accumulated orthogonal transformation the pair rotations
   real(wp), intent(inout) :: trafo(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: nocc, nop, npairs, nactive, isw, iround, ipair, imo, jmo, lmo, kop, iap
   integer, allocatable :: pairs(:, :, :), active_pairs(:, :)
   real(wp) :: delta, off, sin4_coeff, cos4_coeff, theta, s, c, rmax, tan4screen
   real(wp) :: tmpi, tmpj, opnorm, threshold, screen, maxangle, ntri
   real(wp), allocatable :: sine(:), cosine(:), active_sine(:), active_cosine(:)
   logical, allocatable :: rotate(:)
   logical :: use_omp, converged

   nocc = size(opmat, 1)
   nop = size(opmat, 3)

   if (nocc < 2) return

   ! Scale the base convergence threshold by calculation accuracy
   threshold = default_accuracy * accuracy

   ! Calculate the summed operator norm
   opnorm = 0.0_wp
   do kop = 1, nop
      do imo = 1, nocc
         opnorm = opnorm + opmat(imo, imo, kop)**2
         do jmo = 1, imo - 1
            opnorm = opnorm + opmat(imo, jmo, kop)**2
         end do
      end do
   end do
   if (opnorm <= 0.0_wp) return
   ! Operator-scale cutoff below which a pair rotation is skipped as insignificant
   ntri = 0.5_wp * real(nocc, wp) * real(nocc + 1, wp)
   opnorm = sqrt(opnorm / (ntri * nop)) * threshold

   ! Build the round-robin schedule of disjoint occupied-orbital pairs
   call build_round_robin_schedule(nocc, pairs)
   npairs = size(pairs, 2)

   allocate(sine(npairs), cosine(npairs), active_sine(npairs), active_cosine(npairs))
   allocate(rotate(npairs))
   allocate(active_pairs(2, npairs))

   screen = initial_screen
   converged = .false.

   ! Fall back to single thread for small problems
   use_omp = npairs * nop >= omp_min_work

   !$omp parallel if(use_omp) default(none) &
   !$omp& shared(self, nocc, nop, opmat, trafo, pairs, npairs, sine, cosine) &
   !$omp& shared(rotate, maxangle, screen, threshold, opnorm, converged, rmax) &
   !$omp& shared(tan4screen, active_pairs, active_sine, active_cosine, nactive) &
   !$omp& private(isw, iround, ipair, imo, jmo, lmo, kop, sin4_coeff, cos4_coeff) &
   !$omp& private(delta, off, theta, s, c, tmpi, tmpj, iap)
   do isw = 1, self%max_sweeps
      !$omp single
      maxangle = 0.0_wp
      rmax = 0.0_wp
      ! Screening threshhold avoiding atan2 calls for small angles, since:
      tan4screen = tan(4.0_wp * screen)
      !$omp end single

      ! Loop over individual rounds of disjoint pairs each visited once per sweep
      do iround = 1, size(pairs, 3)
         ! Loop over disjoint pairs, compute rotation angles, and rotate columns
         !$omp do schedule(static) reduction(max:maxangle, rmax)
         do ipair = 1, npairs
            sine(ipair) = 0.0_wp
            cosine(ipair) = 1.0_wp
            rotate(ipair) = .false.
            imo = pairs(1, ipair, iround)
            jmo = pairs(2, ipair, iround)
            ! Skip dummy pairs
            if (imo == 0) cycle

            ! Accumulate the sin(4*theta) and cos(4*theta) coefficients
            ! of the pair-localization objective
            sin4_coeff = 0.0_wp
            cos4_coeff = 0.0_wp
            do kop = 1, nop
               delta = opmat(jmo, jmo, kop) - opmat(imo, imo, kop)
               off = opmat(imo, jmo, kop)
               sin4_coeff = sin4_coeff + delta * off
               cos4_coeff = cos4_coeff + 0.25_wp * delta**2 - off**2
            end do
            if (abs(sin4_coeff) + abs(cos4_coeff) < opnorm) cycle

            ! Skip small angles and track the largest skipped |a|/b
            if (cos4_coeff > 0.0_wp) then
               if (abs(sin4_coeff) < cos4_coeff * tan4screen) then
                  rmax = max(rmax, abs(sin4_coeff) / cos4_coeff)
                  cycle
               end if
            end if

            ! Fourth-angle rotation that diagonalizes this pair criterion
            theta = 0.25_wp * atan2(sin4_coeff, cos4_coeff)

            ! Track the largest angle among all non-negligible pairs
            ! for next update of the screening threshold
            maxangle = max(maxangle, abs(theta))
            if (abs(theta) < screen) cycle

            ! Calculate and store sine and cosine for this pair
            sine(ipair) = sin(theta)
            cosine(ipair) = cos(theta)
            rotate(ipair) = .true.

            ! Add the jacobi rotation to the accumulated transformation matrix
            s = sine(ipair)
            c = cosine(ipair)
            do lmo = 1, nocc
               tmpj = trafo(lmo, jmo)
               trafo(lmo, jmo) = c * trafo(lmo, jmo) + s * trafo(lmo, imo)
               trafo(lmo, imo) = c * trafo(lmo, imo) - s * tmpj
            end do

            ! Right multiplication owns disjoint columns, including intersections.
            do kop = 1, nop
               do lmo = 1, nocc
                  tmpi = opmat(lmo, imo, kop)
                  tmpj = opmat(lmo, jmo, kop)
                  opmat(lmo, imo, kop) = c * tmpi - s * tmpj
                  opmat(lmo, jmo, kop) = s * tmpi + c * tmpj
               end do
            end do
         end do
         !$omp end do

         ! Compact the accepted rotations once per round for the left multiplication
         !$omp single
         nactive = 0
         do ipair = 1, npairs
            if (.not. rotate(ipair)) cycle
            nactive = nactive + 1
            active_pairs(1, nactive) = pairs(1, ipair, iround)
            active_pairs(2, nactive) = pairs(2, ipair, iround)
            active_sine(nactive) = sine(ipair)
            active_cosine(nactive) = cosine(ipair)
         end do
         !$omp end single

         ! Complete the transformation by applying the active rotations from the left
         !$omp do collapse(2) schedule(static)
         do kop = 1, nop
            do lmo = 1, nocc
               do iap = 1, nactive
                  imo = active_pairs(1, iap)
                  jmo = active_pairs(2, iap)
                  s = active_sine(iap)
                  c = active_cosine(iap)
                  tmpi = opmat(imo, lmo, kop)
                  tmpj = opmat(jmo, lmo, kop)
                  opmat(imo, lmo, kop) = c * tmpi - s * tmpj
                  opmat(jmo, lmo, kop) = s * tmpi + c * tmpj
               end do
            end do
         end do
         !$omp end do
      end do

      !$omp single
      ! After sweep recover the overall larges angle among applied and skipped pairs
      maxangle = max(maxangle, 0.25_wp * atan(rmax))
      ! Check for convergence
      converged = maxangle <= threshold
      ! Update screening threshold for next sweep limited by convergence threshold
      screen = min(maxangle**2, 0.25_wp * screen)
      screen = max(threshold, screen)
      !$omp end single
      if (converged) exit
   end do
   !$omp end parallel

   if (.not. converged) then
      call fatal_error(error, "Jacobi sweeps orbital localization did not converge")
   end if

end subroutine optimize

end module tblite_wavefunction_localization_jacobi
