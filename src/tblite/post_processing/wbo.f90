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
module tblite_post_processing_bond_orders
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_container, only : container_cache
   use tblite_context, only : context_type
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_integral_type, only : integral_type
   use tblite_output_format, only : format_string
   use tblite_post_processing_type, only : post_processing_type
   use tblite_results, only : results_type
   use tblite_timer, only : timer_type, format_time
   use tblite_wavefunction_type, only : wavefunction_type, get_density_matrix
   use tblite_wavefunction_mulliken, only : get_mayer_bond_orders, get_mayer_bond_orders_uhf
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: new_wbo, wiberg_bond_orders

   type, extends(post_processing_type) :: wiberg_bond_orders
   contains
      procedure :: compute
      procedure :: print_timer
   end type wiberg_bond_orders

   character(len=24), parameter :: label = "Mayer-Wiberg bond orders"
   type(timer_type) :: timer
contains

subroutine new_wbo(new_wbo_type)
   type(wiberg_bond_orders), intent(inout) :: new_wbo_type
   new_wbo_type%label = label
end subroutine new_wbo

subroutine compute(self, mol, wfn, integrals, calc, cache_list, ctx, prlevel, dict)
   class(wiberg_bond_orders),intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction strcuture data
   type(wavefunction_type), intent(in) :: wfn
   !> integral container
   type(integral_type), intent(in) :: integrals
   !> calculator instance
   type(xtb_calculator), intent(in) :: calc
   !> Cache list for storing caches of various interactions
   type(container_cache), intent(inout) :: cache_list(:)
   !> Context container for writing to stdout
   type(context_type), intent(inout) :: ctx
   !> Print level
   integer, intent(in) :: prlevel
   !> Dictionary for storing results
   type(double_dictionary_type), intent(inout) :: dict
   real(kind=wp), allocatable :: wbo(:, :, :), wbo_2d(:, :)
   real(kind=wp), allocatable :: focc_(:, :)
   integer :: nspin, i, j
   real(kind=wp) :: nel_
   real(wp), allocatable :: pmat(:, :, :)

   call timer%push("total")
   nspin = size(wfn%density, dim=3)

   if ((nspin == 1) .and. (wfn%nel(1) /= wfn%nel(2))) then
      allocate(wbo(mol%nat, mol%nat, 2), source=0.0_wp)
      allocate(focc_(calc%bas%nao, 2), source=0.0_wp)
      allocate(pmat(calc%bas%nao, calc%bas%nao, 2))
      do j = 1,2
         nel_ = wfn%nel(j)
         do i = 1, size(wfn%focc)
            if (nel_ > 1.0_wp) then
               focc_(i,j) = 1.0_wp
               nel_ = nel_ - 1.0_wp
            else
               focc_(i,j) = nel_
               exit
            end if
         end do
      end do
      do j = 1, 2
         call get_density_matrix(focc_(:, j), wfn%coeff(:, :, nspin), pmat(:, :, j))
      end do
      call get_mayer_bond_orders_uhf(calc%bas, integrals%overlap, pmat, wbo)
      wbo_2d = 2*wbo(:, :, 1)
      call dict%add_entry("bond-orders", wbo_2d)
   else
      allocate(wbo(mol%nat, mol%nat, nspin), source=0.0_wp)
      call get_mayer_bond_orders(calc%bas, integrals%overlap, wfn%density, wbo)
      if (size(wbo, dim = 3) == 1) then 
         wbo_2d = wbo(:, : , 1)
         call dict%add_entry("bond-orders", wbo_2d)
      else
         call dict%add_entry("bond-orders", wbo(:, :, :))
      end if
   end if

   
   call timer%pop()
end subroutine compute

subroutine print_timer(self, prlevel, ctx)
   !> Instance of the interaction container
   class(wiberg_bond_orders), intent(in) :: self
   integer :: prlevel
   type(context_type) :: ctx
   real(wp) :: ttime, stime
   integer :: it
   character(len=*), parameter :: labels(*) = [character(len=20):: &
      & ]
   if (prlevel > 2) then
      call ctx%message(label//" timing details:")
      ttime = timer%get("total")
      call ctx%message(" total:"//repeat(" ", 16)//format_time(ttime))
      do it = 1, size(labels)
         stime = timer%get(labels(it))
         if (stime <= epsilon(0.0_wp)) cycle
         call ctx%message(" - "//labels(it)//format_time(stime) &
            & //" ("//format_string(int(stime/ttime*100), '(i3)')//"%)")
      end do
      call ctx%message("")
   end if
end subroutine print_timer

end module tblite_post_processing_bond_orders
