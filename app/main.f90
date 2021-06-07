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

program main_driver
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, input_unit
   use mctc_env, only : error_type, fatal_error, get_argument, wp
   use mctc_io, only : structure_type, read_structure, write_structure, &
      & filetype, get_filetype, to_symbol
   use tblite_context_type, only : context_type
   use tblite_output_ascii
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_ipea1, only : new_ipea1_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   use tblite_version, only : get_tblite_version
   implicit none
   character(len=*), parameter :: prog_name = "tblite"

   type :: driver_config
      character(len=:), allocatable :: input
      integer, allocatable :: input_format
      logical :: grad = .false.
      character(len=:), allocatable :: grad_output
      integer :: verbosity = 2
      integer, allocatable :: charge
      integer, allocatable :: spin
      character(len=:), allocatable :: method
      logical :: json = .false.
      character(len=:), allocatable :: json_output
      real(wp) :: accuracy = 1.0_wp
   end type driver_config

   type(driver_config) :: config
   type(structure_type) :: mol
   type(error_type), allocatable :: error
   integer :: spin, charge, stat, unit
   logical :: exist
   real(wp) :: energy
   real(wp), allocatable :: gradient(:, :), sigma(:, :)
   type(context_type) :: ctx
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn

   call get_arguments(config, error)
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   if (config%input == "-") then
      if (.not.allocated(config%input_format)) config%input_format = filetype%xyz
      call read_structure(mol, input_unit, config%input_format, error)
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

   if (.not.allocated(config%method)) config%method = "gfn2"
   select case(config%method)
   case default
      write(error_unit, '(a)') "Unknown method '"//config%method//"' requested"
      error stop
   case("gfn2")
      call new_gfn2_calculator(calc, mol)
   case("gfn1")
      call new_gfn1_calculator(calc, mol)
   case("ipea1")
      call new_ipea1_calculator(calc, mol)
   end select

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

contains


subroutine help(unit)
   integer, intent(in) :: unit

   write(unit, '(a, *(1x, a))') &
      "Usage: "//prog_name//" [options] <input>"

   write(unit, '(a)') &
      "", &
      ""

   write(unit, '(2x, a, t25, a)') &
      "-c, --charge <real>", "Set charge to molecule, overwrites .CHRG file", &
      "    --spin <int>", "Set number of unpaired electrons, overwrites .UHF file", &
      "    --method <name>", "Parametrization of the xTB Hamiltonian to use", &
      "", "Available methods: gfn1, gfn2 (Default: gfn2)", &
      "    --grad [file]", "Evaluate molecular gradient and virial", &
      "", "Results are stored in file (default: tblite.txt)", &
      "    --json [file]", "Dump results as JSON output (default: tblite.json)", &
      "-i, --input <format>", "Hint for the format of the input file", &
      "    --version", "Print program version and exit", &
      "    --help", "Show this help message"

   write(unit, '(a)')

end subroutine help


subroutine version(unit)
   integer, intent(in) :: unit
   character(len=:), allocatable :: version_string

   call get_tblite_version(string=version_string)
   write(unit, '(a, *(1x, a))') &
      & prog_name, "version", version_string

end subroutine version


subroutine get_argument_as_real(iarg, val, error)
   !> Index of command line argument, range [0:command_argument_count()]
   integer, intent(in) :: iarg
   !> Real value
   real(wp), intent(out) :: val
   !> Error handling
   type(error_type), allocatable :: error

   integer :: stat
   character(len=:), allocatable :: arg

   call get_argument(iarg, arg)
   if (.not.allocated(arg)) then
      call fatal_error(error, "Cannot read real value, argument missing")
      return
   end if
   read(arg, *, iostat=stat) val
   if (stat /= 0) then
      call fatal_error(error, "Cannot read real value from '"//arg//"'")
      return
   end if

end subroutine get_argument_as_real


subroutine get_argument_as_int(iarg, val, error)
   !> Index of command line argument, range [0:command_argument_count()]
   integer, intent(in) :: iarg
   !> Real value
   integer, intent(out) :: val
   !> Error handling
   type(error_type), allocatable :: error

   integer :: stat
   character(len=:), allocatable :: arg

   call get_argument(iarg, arg)
   if (.not.allocated(arg)) then
      call fatal_error(error, "Cannot read integer value, argument missing")
      return
   end if
   read(arg, *, iostat=stat) val
   if (stat /= 0) then
      call fatal_error(error, "Cannot read integer value from '"//arg//"'")
      return
   end if

end subroutine get_argument_as_int


subroutine get_arguments(config, error)

   type(driver_config), intent(out) :: config
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg
   logical :: getopts
   character(len=:), allocatable :: arg

   iarg = 0
   getopts = .true.
   narg = command_argument_count()
   do while(iarg < narg)
      iarg = iarg + 1
      call get_argument(iarg, arg)
      if (.not.getopts) then
         if (.not.allocated(config%input)) then
            call move_alloc(arg, config%input)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      end if
      select case(arg)
      case("--")
         getopts = .false.
      case("--help")
         call help(output_unit)
         stop
      case("--version")
         call version(output_unit)
         stop
      case default
         if (.not.allocated(config%input)) then
            call move_alloc(arg, config%input)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      case("-i", "--input")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for input format")
            exit
         end if
         config%input_format = get_filetype("."//arg)
      case("-c", "--charge")
         iarg = iarg + 1
         allocate(config%charge)
         call get_argument_as_int(iarg, config%charge, error)
         if (allocated(error)) exit
      case("--spin")
         iarg = iarg + 1
         allocate(config%spin)
         call get_argument_as_int(iarg, config%spin, error)
         if (allocated(error)) exit
      case("--method")
         iarg = iarg + 1
         call get_argument(iarg, config%method)
         if (.not.allocated(config%method)) then
            call fatal_error(error, "Missing argument for method")
            exit
         end if
      case("--acc")
         iarg = iarg + 1
         call get_argument_as_real(iarg, config%accuracy, error)
         if (allocated(error)) exit
      case("--grad")
         config%grad = .true.
         config%grad_output = "tblite.txt"
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (allocated(arg)) then
            if (arg(1:1) == "-") then
               iarg = iarg - 1
               cycle
            end if
            call move_alloc(arg, config%grad_output)
         end if
      case("--json")
         config%json = .true.
         config%json_output = "tblite.json"
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (allocated(arg)) then
            if (arg(1:1) == "-") then
               iarg = iarg - 1
               cycle
            end if
            call move_alloc(arg, config%json_output)
         end if
      end select
   end do

   if (.not.(allocated(config%input))) then
      if (.not.allocated(error)) then
         call help(output_unit)
         error stop
      end if
   end if

end subroutine get_arguments


end program main_driver
