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

!> Implementation of the driver entry points for the different commands
module tblite_driver
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, input_unit
   use mctc_env, only : error_type, fatal_error, get_argument, wp
   use mctc_io, only : structure_type, read_structure, write_structure, &
      & filetype, get_filetype, to_symbol
   use tblite_cli, only : run_config, param_config
   use tblite_context_type, only : context_type
   use tblite_output_ascii
   use tblite_param, only : param_record
   use tblite_toml, only : toml_table, merge_table
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator, new_xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator, export_gfn2_param
   use tblite_xtb_gfn1, only : new_gfn1_calculator, export_gfn1_param
   use tblite_xtb_ipea1, only : new_ipea1_calculator, export_ipea1_param
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   use tblite_version, only : get_tblite_version
   implicit none
   private

   public :: main

   interface main
      module procedure :: run_main
      module procedure :: param_main
   end interface

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
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

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
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 300.0_wp * 3.166808578545117e-06_wp)

   call xtb_singlepoint(ctx, mol, calc, wfn, config%accuracy, energy, gradient, sigma, 2)
   if (ctx%failed()) then
      write(error_unit, '("[Fatal]", 1x, a)') "Singlepoint calculation failed"
      do while(ctx%failed())
         call ctx%get_error(error)
         write(error_unit, '("->", 1x, a)') error%message
      end do
      error stop
   end if

   call ascii_levels(ctx%unit, config%verbosity, wfn%homoa, wfn%emo, wfn%focc, 7)

   if (config%grad) then
      open(file=config%grad_output, newunit=unit)
      call tagged_result(unit, energy, gradient, sigma)
      close(unit)
      if (config%verbosity > 0) then
         write(output_unit, '(a)') &
            & "[Info] Dispersion results written to '"//config%grad_output//"'"
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


subroutine param_main(config, error)
   type(param_config), intent(in) :: config
   type(error_type), allocatable, intent(out) :: error

   type(param_record), allocatable :: param, base
   type(toml_table) :: table1, table2

   if (allocated(config%method)) then
      allocate(base)
      select case(config%method)
      case default
         call fatal_error(error, "Unknown method '"//config%method//"' requested")
      case("gfn2")
         call export_gfn2_param(base)
      case("gfn1")
         call export_gfn1_param(base)
      case("ipea1")
         call export_ipea1_param(base)
      end select
   end if
   if (allocated(error)) return

   if (allocated(config%input)) then
      allocate(param)
      call param%load(config%input, error)

      if (.not.allocated(error) .and. allocated(base)) then
         table1 = toml_table()
         call param%dump(table1, error)
         if (.not. allocated(error)) then
            table2 = toml_table()
            call base%dump(table2, error)
         end if
         if (.not. allocated(error)) then
            call merge_table(table1, table2)
            call table2%destroy
            deallocate(param, base)
            allocate(param)
            call param%load(table1, error)
            call table1%destroy
         end if
      end if
   else
      call move_alloc(base, param)
   end if
   if (allocated(error)) return

   if (allocated(config%output)) then
      call param%dump(config%output, error)
   end if
   if (allocated(error)) return
end subroutine param_main


end module tblite_driver
