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

!> @file tblite/post_processing/wbo.f90
!> Implements the calculation of Wiberg-Mayer bond orders as post processing method.
module tblite_post_processing_wbo
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_container_list, only : cache_list
   use tblite_context, only : context_type
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_integral_type, only : integral_type
   use tblite_mpi_utils, only : mpi_allreduce_sum, mpi_multiply_columns, &
      & mpi_transpose_columns, mpi_density_columns
   use tblite_post_processing_type, only : post_processing_type
   use tblite_timer, only : timer_type, format_time
   use tblite_wavefunction_mulliken, only : get_mayer_bond_orders, contract_mayer_columns
   use tblite_wavefunction_spin, only : updown_to_magnet
   use tblite_wavefunction_type, only : wavefunction_type, get_density_matrix
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: new_wiberg_bond_orders, wiberg_bond_orders

   !> Wiberg-Mayer bond orders as post-processing method
   type, extends(post_processing_type) :: wiberg_bond_orders
   contains
      !> Calculate Wiberg-Mayer bond orders
      procedure :: compute
      !> Print timings
      procedure :: print_timer
   end type wiberg_bond_orders

   character(len=24), parameter :: label = "Mayer-Wiberg bond orders"

contains

subroutine new_wiberg_bond_orders(self)
   !> Instance of the Wiberg-Mayer bond order post-processing
   type(wiberg_bond_orders), intent(out) :: self

   self%label = label
   self%local_matrices = .true.

end subroutine new_wiberg_bond_orders

subroutine compute(self, mol, wfn, ints, calc, caches, accuracy, ctx, timer, &
   & prlevel, dict)
   !> Instance of the Wiberg-Mayer bond order post-processing
   class(wiberg_bond_orders),intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction strcuture data
   type(wavefunction_type), intent(in) :: wfn
   !> Integral container
   type(integral_type), intent(in) :: ints
   !> Calculator instance
   type(xtb_calculator), intent(in) :: calc
   !> Cache list for storing caches of various interactions
   type(cache_list), intent(inout) :: caches
   !> Accuracy for computation
   real(wp), intent(in) :: accuracy
   !> Context container for writing to stdout
   type(context_type), intent(inout) :: ctx
   !> Timer instance
   type(timer_type), intent(inout) :: timer
   !> Print level
   integer, intent(in) :: prlevel
   !> Dictionary for storing results
   type(double_dictionary_type), intent(inout) :: dict

   real(wp), allocatable :: wbo(:, :, :), pmat(:, :, :), ps(:, :), pst(:, :)
   type(error_type), allocatable :: error
   integer :: spin, nspin, i
   integer, allocatable :: columns(:)
   logical :: restricted_open

   call timer%push("wbo")
   restricted_open = wfn%nspin == 1 .and. wfn%nel(1) /= wfn%nel(2)
   nspin = wfn%nspin
   if (restricted_open) nspin = 2
   allocate(wbo(mol%nat, mol%nat, nspin))
   if (restricted_open) allocate(pmat(calc%bas%nao, size(ints%overlap, 2), 1))
   if (ints%local) then
      columns = [(i, i=ints%columns(1), ints%columns(2))]
      allocate(ps(calc%bas%nao, size(columns)), pst(calc%bas%nao, size(columns)))
   end if

   do spin = 1, nspin
      if (restricted_open) then
         if (ints%local) then
            call mpi_density_columns(error, wfn%focc(:, spin), wfn%coeff(:, :, 1), pmat(:, :, 1), ctx%comm)
         else
            call get_density_matrix(wfn%focc(:, spin), wfn%coeff(:, :, 1), pmat(:, :, 1))
         end if
      end if
      if (ints%local) then
         if (restricted_open) then
            call mpi_multiply_columns(error, pmat(:, :, 1), ints%overlap, ps, ctx%comm)
         else
            call mpi_multiply_columns(error, wfn%density(:, :, spin), ints%overlap, ps, ctx%comm)
         end if
         call mpi_transpose_columns(error, ps, pst, ctx%comm)
         if (allocated(error)) exit
         call contract_mayer_columns(calc%bas, ps, pst, wbo(:, :, spin), columns, transposed=.true.)
      else
         if (restricted_open) then
            call get_mayer_bond_orders(mol, calc%bas, ints%overlap, pmat, wbo(:, :, spin:spin), ctx%partition)
         else
            call get_mayer_bond_orders(mol, calc%bas, ints%overlap, &
               & wfn%density(:, :, spin:spin), wbo(:, :, spin:spin), ctx%partition)
         end if
      end if
   end do
   if (nspin == 2) wbo = 2*wbo
   call updown_to_magnet(wbo)

   call mpi_allreduce_sum(error, wbo, ctx%comm)
   call timer%pop()
   if (allocated(error)) then
      call ctx%set_error(error)
      return
   end if
   call dict%add_entry("bond-orders", wbo)

end subroutine compute

subroutine print_timer(self, timer, prlevel, ctx)
   !> Instance of the Wiberg-Mayer bond order post-processing
   class(wiberg_bond_orders), intent(in) :: self
   !> Timer instance
   type(timer_type), intent(in) :: timer
   !> Print level
   integer, intent(in) :: prlevel
   !> Context container for writing to stdout
   type(context_type), intent(inout) :: ctx

   real(wp) :: ttime

   if (prlevel > 2) then
      call ctx%message(label//" timing details:")
      ttime = timer%get("wbo")
      call ctx%message(" total:"//repeat(" ", 16)//format_time(ttime))
      call ctx%message("")
   end if

end subroutine print_timer

end module tblite_post_processing_wbo
