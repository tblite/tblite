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

!> Distributing a single point calculation over MPI ranks must reproduce the
!> serial result on every rank.
program test_mpi_singlepoint
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use mpi, only : MPI_COMM_WORLD, MPI_Abort, MPI_Finalize, MPI_Init
   use mstore, only : get_structure
   use tblite_context, only : context_type
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none

   real(wp), parameter :: thr = 1.0e-8_wp

   type(structure_type) :: mol
   type(error_type), allocatable :: error
   real(wp) :: serial_energy, mpi_energy
   real(wp) :: serial_sigma(3, 3), mpi_sigma(3, 3)
   real(wp), allocatable :: serial_gradient(:, :), mpi_gradient(:, :)
   integer :: stat

   call MPI_Init(stat)

   call get_structure(mol, "MB16-43", "01")
   allocate(serial_gradient(3, mol%nat), mpi_gradient(3, mol%nat))

   call run(mol, .false., serial_energy, serial_gradient, serial_sigma, error)
   call check_error(error)
   call run(mol, .true., mpi_energy, mpi_gradient, mpi_sigma, error)
   call check_error(error)

   call assert(abs(mpi_energy - serial_energy) < thr, "energy")
   call assert(all(abs(mpi_gradient - serial_gradient) < thr), "gradient")
   call assert(all(abs(mpi_sigma - serial_sigma) < thr), "virial")

   call MPI_Finalize(stat)

contains

   subroutine run(mol, distributed, energy, gradient, sigma, error)
      type(structure_type), intent(in) :: mol
      logical, intent(in) :: distributed
      real(wp), intent(out) :: energy, gradient(:, :), sigma(:, :)
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn

      ctx%verbosity = 0
      call new_gfn2_calculator(calc, mol, error)
      if (allocated(error)) return

      if (distributed) then
         call ctx%set_mpi(error)
         if (allocated(error)) return
         call calc%set_partition(ctx%partition)
      end if

      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, 300.0_wp)
      call xtb_singlepoint(ctx, mol, calc, wfn, 1.0_wp, energy, gradient, sigma)
      if (ctx%failed()) call ctx%get_error(error)
   end subroutine run

   subroutine check_error(error)
      type(error_type), allocatable, intent(in) :: error
      if (allocated(error)) then
         write(*, '(2a)') "[Fatal] ", error%message
         call MPI_Abort(MPI_COMM_WORLD, 1, stat)
      end if
   end subroutine check_error

   subroutine assert(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not.condition) then
         write(*, '(3a)') "[Fatal] Distributed ", label, " does not match the serial result"
         call MPI_Abort(MPI_COMM_WORLD, 1, stat)
      end if
   end subroutine assert

end program test_mpi_singlepoint
