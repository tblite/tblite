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

!> Implementation of the driver entry points for the singlepoint runner
module tblite_driver_run
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, input_unit
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type, read_structure, filetype
   use mctc_io_constants, only : codata
   use mctc_io_convert, only : aatoau, ctoau
   use tblite_cli, only : run_config
   use tblite_container, only : container_type
   use tblite_context, only : context_type
   use tblite_external_field, only : electric_field
   use tblite_output_ascii
   use tblite_param, only : param_record
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator, new_xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator, export_gfn2_param
   use tblite_xtb_gfn1, only : new_gfn1_calculator, export_gfn1_param
   use tblite_xtb_ipea1, only : new_ipea1_calculator, export_ipea1_param
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: main

   interface main
      module procedure :: run_main
   end interface

   real(wp), parameter :: kt = 3.166808578545117e-06_wp
   real(wp), parameter :: jtoau = 1.0_wp / (codata%me*codata%c**2*codata%alpha**2)
   !> Convert V/Å = J/(C·Å) to atomic units
   real(wp), parameter :: vatoau = jtoau / (ctoau * aatoau)

contains


subroutine run_main(config, error)
   type(run_config), intent(in) :: config
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   character(len=:), allocatable :: method
   integer :: spin, charge, stat, unit
   logical :: exist
   real(wp) :: energy
   real(wp), allocatable :: gradient(:, :), sigma(:, :)
   type(param_record) :: param
   type(context_type) :: ctx
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn

   if (config%input == "-") then
      if (allocated(config%input_format)) then
         call read_structure(mol, input_unit, config%input_format, error)
      else
         call read_structure(mol, input_unit, filetype%xyz, error)
      end if
   else
      call read_structure(mol, config%input, error, config%input_format)
   end if
   if (allocated(error)) return

   if (allocated(config%charge)) then
      mol%charge = config%charge
   else
      inquire(file='.CHRG', exist=exist)
      if (exist) then
         open(file='.CHRG', newunit=unit)
         read(unit, *, iostat=stat) charge
         if (stat == 0) then
            mol%charge = charge
            if (config%verbosity > 0) write(output_unit, '(a)') &
               "[Info] Molecular charge read from .CHRG"
         else
            if (config%verbosity > 0) write(output_unit, '(a)') &
               "[Warn] Could not read molecular charge read from .CHRG"
         end if
         close(unit)
      end if
   end if

   if (allocated(config%spin)) then
      mol%uhf = config%spin
   else
      inquire(file='.UHF', exist=exist)
      if (exist) then
         open(file='.UHF', newunit=unit)
         read(unit, *, iostat=stat) spin
         if (stat == 0) then
            mol%uhf = spin
            if (config%verbosity > 0) write(output_unit, '(a)') &
               "[Info] Molecular spin read from .UHF"
         else
            if (config%verbosity > 0) write(output_unit, '(a)') &
               "[Warn] Could not read molecular spin read from .UHF"
         end if
         close(unit)
      end if
   end if

   if (config%grad) then
      allocate(gradient(3, mol%nat), sigma(3, 3))
   end if

   if (allocated(config%param)) then
      call param%load(config%param, error)
      if (.not. allocated(error)) then
         call new_xtb_calculator(calc, mol, param, error)
      end if
   else
      method = "gfn2"
      if (allocated(config%method)) method = config%method
      select case(method)
      case default
         call fatal_error(error, "Unknown method '"//method//"' requested")
      case("gfn2")
         call new_gfn2_calculator(calc, mol)
      case("gfn1")
         call new_gfn1_calculator(calc, mol)
      case("ipea1")
         call new_ipea1_calculator(calc, mol)
      end select
   end if
   if (allocated(error)) return

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, config%etemp * kt)

   if (allocated(config%efield)) then
      block
         class(container_type), allocatable :: cont
         cont = electric_field(config%efield*vatoau)
         call calc%push_back(cont)
      end block
   end if

   call xtb_singlepoint(ctx, mol, calc, wfn, config%accuracy, energy, gradient, sigma, &
      & config%verbosity)
   if (ctx%failed()) then
      write(error_unit, '("[Fatal]", 1x, a)') "Singlepoint calculation failed"
      do while(ctx%failed())
         call ctx%get_error(error)
         write(error_unit, '("->", 1x, a)') error%message
      end do
      error stop
   end if

   if (config%verbosity > 2) then
      call ascii_levels(ctx%unit, config%verbosity, wfn%homoa, wfn%emo, wfn%focc, 7)
   end if

   if (allocated(config%grad_output)) then
      open(file=config%grad_output, newunit=unit)
      call tagged_result(unit, energy, gradient, sigma)
      close(unit)
      if (config%verbosity > 0) then
         write(output_unit, '(a)') &
            & "[Info] Tight-binding results written to '"//config%grad_output//"'"
      end if
   end if

   if (config%json) then
      open(file=config%json_output, newunit=unit)
      call json_results(unit, "  ", energy=energy, gradient=gradient, sigma=sigma)
      close(unit)
      if (config%verbosity > 0) then
         write(output_unit, '(a)') &
            & "[Info] JSON dump of results written to '"//config%json_output//"'"
      end if
   end if
end subroutine run_main


end module tblite_driver_run
