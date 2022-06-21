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

!> Definition of the command-line interface of the tblite framework
module tblite_cli
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, input_unit
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : filetype, get_filetype, to_symbol
   use tblite_argument, only : argument_list, len
   use tblite_cli_help, only : prog_name, help_text, help_text_run, help_text_param, &
      & help_text_fit, help_text_tagdiff
   use tblite_features, only : get_tblite_feature
   use tblite_lapack_solver, only : lapack_algorithm
   use tblite_solvation, only : solvation_input, cpcm_input, alpb_input, &
      & solvent_data, get_solvent_data
   use tblite_version, only : get_tblite_version
   implicit none
   private

   public :: get_arguments, driver_config, run_config, param_config, fit_config, &
      & tagdiff_config


   type, abstract :: driver_config
      logical :: color = .false.
   end type driver_config

   !> Configuration for evaluating tight binding model on input structure
   type, extends(driver_config) :: run_config
      !> Geometry input file
      character(len=:), allocatable :: input
      !> Format for reading the input file
      integer, allocatable :: input_format
      !> Evaluate gradient
      logical :: grad = .false.
      !> File for output of gradient information
      character(len=:), allocatable :: grad_output
      !> Verbosity of calculation
      integer :: verbosity = 2
      !> Total charge of the system
      integer, allocatable :: charge
      !> Number of unpaired electrons
      integer, allocatable :: spin
      !> Parametrization of the xTB Hamiltonian to use
      character(len=:), allocatable :: method
      !> Parametrization file to use for calculation
      character(len=:), allocatable :: param
      !> Guess for SCF calculation
      character(len=:), allocatable :: guess
      !> Create JSON dump
      logical :: json = .false.
      !> File for output of JSON dump
      character(len=:), allocatable :: json_output
      !> Input for solvation model
      type(solvation_input), allocatable :: solvation
      !> Numerical accuracy for self-consistent iterations
      real(wp) :: accuracy = 1.0_wp
      !> Maximum number of iterations for SCF
      integer, allocatable :: max_iter
      !> Electronic temperature
      real(wp) :: etemp = 300.0_wp
      !> Electric field
      real(wp), allocatable :: efield(:)
      !> Spin polarization
      logical :: spin_polarized = .false.
      !> Algorithm for electronic solver
      integer :: solver = lapack_algorithm%gvd
   end type run_config

   type, extends(driver_config) :: param_config
      integer :: verbosity = 2
      character(len=:), allocatable :: input
      character(len=:), allocatable :: output
      character(len=:), allocatable :: method
   end type param_config

   type, extends(driver_config) :: fit_config
      integer :: verbosity = 2
      character(len=:), allocatable :: prog
      character(len=:), allocatable :: param
      character(len=:), allocatable :: input
      character(len=:), allocatable :: copy_input
      logical :: dry_run = .false.
   end type fit_config

   type, extends(driver_config) :: tagdiff_config
      character(len=:), allocatable :: actual
      character(len=:), allocatable :: reference
      logical :: fit = .false.
   end type tagdiff_config

contains

subroutine get_arguments(config, error)

   !> Command line options
   class(driver_config), allocatable, intent(out) :: config

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(argument_list) :: list
   integer :: iarg, narg
   logical :: getopts
   character(len=:), allocatable :: arg
   logical, allocatable :: color

   list = argument_list()

   iarg = 0
   getopts = .true.
   narg = len(list)
   do while(iarg < narg)
      iarg = iarg + 1
      call list%get(iarg, arg)
      select case(arg)
      case("--help")
         write(output_unit, '(a)') help_text
         stop
      case("--version")
         call version(output_unit)
         stop
      case default
         iarg = iarg - 1
         allocate(run_config :: config)
         exit
      case("--color")
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for color options")
            exit
         end if
         select case(arg)
         case default
            call fatal_error(error, "Unknown option for color support '"//arg//"'")
         case("always")
            color = .true.
         case("never")
            color = .false.
         case("auto")
            color = get_tblite_feature("color")
         end select
      case("fit")
         allocate(fit_config :: config)
         exit
      case("param")
         allocate(param_config :: config)
         exit
      case("run")
         allocate(run_config :: config)
         exit
      case("tagdiff")
         allocate(tagdiff_config :: config)
         exit
      end select
   end do
   if (allocated(error)) return

   if (.not.allocated(config)) then
      write(output_unit, '(a)') help_text
      error stop
   end if

   select type(config)
   type is (fit_config)
      call get_fit_arguments(config, list, iarg, error)
   type is (param_config)
      call get_param_arguments(config, list, iarg, error)
   type is (run_config)
      call get_run_arguments(config, list, iarg, error)
   type is (tagdiff_config)
      call get_tagdiff_arguments(config, list, iarg, error)
   end select
   if (allocated(error)) return

   if (allocated(color)) config%color = color

end subroutine get_arguments


subroutine get_run_arguments(config, list, start, error)

   !> Command line options
   type(run_config), intent(out) :: config

   !> List of command line arguments
   type(argument_list), intent(in) :: list

   !> First command line argument
   integer, intent(in) :: start

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg
   logical :: getopts
   character(len=:), allocatable :: arg
   logical :: alpb
   type(solvent_data) :: solvent

   iarg = start
   getopts = .true.
   narg = len(list)
   do while(iarg < narg)
      iarg = iarg + 1
      call list%get(iarg, arg)
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
         write(output_unit, '(a)') help_text_run
         stop

      case("--version")
         call version(output_unit)
         stop

      case("-v", "-vv", "--verbose")
         config%verbosity = config%verbosity + 1
         if (arg == "-vv") config%verbosity = config%verbosity + 1

      case("-s", "-ss", "--silent")
         config%verbosity = config%verbosity - 1
         if (arg == "-ss") config%verbosity = config%verbosity - 1

      case default
         if (.not.allocated(config%input)) then
            call move_alloc(arg, config%input)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit

      case("-i", "--input")
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for input format")
            exit
         end if
         config%input_format = get_filetype("."//arg)

      case("-c", "--charge")
         iarg = iarg + 1
         allocate(config%charge)
         call list%get(iarg, arg)
         call get_argument_as_int(arg, config%charge, error)
         if (allocated(error)) exit

      case("--spin")
         iarg = iarg + 1
         allocate(config%spin)
         call list%get(iarg, arg)
         call get_argument_as_int(arg, config%spin, error)
         if (allocated(error)) exit

      case("--spin-polarized")
         config%spin_polarized = .true.

      case("--method")
         if (allocated(config%param)) then
            call fatal_error(error, "Cannot specify method if parameter file is provided")
            exit
         end if
         iarg = iarg + 1
         call list%get(iarg, config%method)
         if (.not.allocated(config%method)) then
            call fatal_error(error, "Missing argument for method")
            exit
         end if

      case("--guess")
         iarg = iarg + 1
         call list%get(iarg, config%guess)

      case("--cpcm")
         if (allocated(config%solvation)) then
            call fatal_error(error, "Cannot use CPCM if ALPB/GBSA is enabled")
            exit
         end if
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for CPCM")
            exit
         end if

         solvent = get_solvent_data(arg)
         if (solvent%eps <= 0.0_wp) then
            call get_argument_as_real(arg, solvent%eps, error)
         end if
         if (allocated(error)) exit
         allocate(config%solvation)
         config%solvation%cpcm = cpcm_input(solvent%eps)

      case("--alpb", "--gbsa")
         if (allocated(config%solvation)) then
            call fatal_error(error, "Cannot use ALPB/GBSA if CPCM is enabled")
            exit
         end if
         alpb = arg == "--alpb"
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for ALPB/GBSA")
            exit
         end if

         solvent = get_solvent_data(arg)
         if (solvent%eps <= 0.0_wp) then
            call get_argument_as_real(arg, solvent%eps, error)
         end if
         if (allocated(error)) exit
         allocate(config%solvation)
         config%solvation%alpb = alpb_input(solvent%eps, alpb=alpb)

      case("--param")
         if (allocated(config%param)) then
            call fatal_error(error, "Cannot specify parameter file if method is provided")
            exit
         end if
         iarg = iarg + 1
          call list%get(iarg, config%param)
         if (.not.allocated(config%param)) then
            call fatal_error(error, "Missing argument for param")
            exit
         end if

      case("--acc")
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%accuracy, error)
         if (allocated(error)) exit

      case("--iterations")
         iarg = iarg + 1
         call list%get(iarg, arg)
         allocate(config%max_iter)
         call get_argument_as_int(arg, config%max_iter, error)
         if (allocated(error)) exit

      case("--solver")
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for solver")
            exit
         end if

         select case(arg)
         case default
            call fatal_error(error, "Unknown electronic solver '"//arg//"' specified")
            exit
         case("gvd")
            config%solver = lapack_algorithm%gvd
         case("gvr")
            config%solver = lapack_algorithm%gvr
         end select

      case("--etemp")
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%etemp, error)
         if (allocated(error)) exit

      case("--efield")
         iarg = iarg + 1
         call list%get(iarg, arg)
         allocate(config%efield(3))
         call get_argument_as_realv(arg, config%efield, error)
         if (allocated(error)) exit

      case("--grad")
         config%grad = .true.
         iarg = iarg + 1
         call list%get(iarg, arg)
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
         call list%get(iarg, arg)
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
         write(output_unit, '(a)') help_text_run
         error stop
      end if
   end if

   if (.not.allocated(config%guess)) config%guess = "sad"

   if (config%grad .and. .not.allocated(config%grad_output) .and. .not.config%json) then
      config%grad_output = "tblite.txt"
   end if
end subroutine get_run_arguments


subroutine get_param_arguments(config, list, start, error)

   !> Command line options
   type(param_config), intent(out) :: config

   !> List of command line arguments
   type(argument_list), intent(in) :: list

   !> First command line argument
   integer, intent(in) :: start

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg
   logical :: getopts
   character(len=:), allocatable :: arg

   iarg = start
   getopts = .true.
   narg = len(list)
   do while(iarg < narg)
      iarg = iarg + 1
      call list%get(iarg, arg)
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
         write(output_unit, '(a)') help_text_param
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
      case("--method")
         iarg = iarg + 1
         call list%get(iarg, config%method)
         if (.not.allocated(config%method)) then
            call fatal_error(error, "Missing argument for method")
            exit
         end if
      case("--output")
         iarg = iarg + 1
         call list%get(iarg, config%output)
         if (.not.allocated(config%output)) then
            call fatal_error(error, "Missing argument for output")
            exit
         end if
      end select
   end do

   if (.not.any([allocated(config%input), allocated(config%method)])) then
      if (.not.allocated(error)) then
         write(output_unit, '(a)') help_text_param
         error stop
      end if
   end if
end subroutine get_param_arguments


subroutine get_fit_arguments(config, list, start, error)

   !> Command line options
   type(fit_config), intent(out) :: config

   !> List of command line arguments
   type(argument_list), intent(in) :: list

   !> First command line argument
   integer, intent(in) :: start

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg
   logical :: getopts
   character(len=:), allocatable :: arg

   config%prog = list%prog

   iarg = start
   getopts = .true.
   narg = len(list)
   do while(iarg < narg)
      iarg = iarg + 1
      call list%get(iarg, arg)
      if (.not.getopts) then
         if (.not.allocated(config%param)) then
            call move_alloc(arg, config%param)
            cycle
         end if
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
         write(output_unit, '(a)') help_text_fit
         stop
      case("--version")
         call version(output_unit)
         stop
      case("--dry-run")
         config%dry_run = .true.
      case("-v", "-vv", "--verbose")
         config%verbosity = config%verbosity + 1
         if (arg == "-vv") config%verbosity = config%verbosity + 1
      case("-s", "-ss", "--silent")
         config%verbosity = config%verbosity - 1
         if (arg == "-ss") config%verbosity = config%verbosity - 1
      case("--copy")
         iarg = iarg + 1
         call list%get(iarg, config%copy_input)
         if (.not.allocated(config%copy_input)) then
            call fatal_error(error, "Missing argument for copy")
            exit
         end if
      case default
         if (.not.allocated(config%param)) then
            call move_alloc(arg, config%param)
            cycle
         end if
         if (.not.allocated(config%input)) then
            call move_alloc(arg, config%input)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      end select
   end do

   if (.not.allocated(config%input)) then
      if (.not.allocated(error)) then
         write(output_unit, '(a)') help_text_fit
         error stop
      end if
   end if
end subroutine get_fit_arguments


subroutine get_tagdiff_arguments(config, list, start, error)

   !> Command line options
   type(tagdiff_config), intent(out) :: config

   !> List of command line arguments
   type(argument_list), intent(in) :: list

   !> First command line argument
   integer, intent(in) :: start

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg
   logical :: getopts
   character(len=:), allocatable :: arg

   iarg = start
   getopts = .true.
   narg = len(list)
   do while(iarg < narg)
      iarg = iarg + 1
      call list%get(iarg, arg)
      if (.not.getopts) then
         if (.not.allocated(config%actual)) then
            call move_alloc(arg, config%actual)
            cycle
         end if
         if (.not.allocated(config%reference)) then
            call move_alloc(arg, config%reference)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      end if
      select case(arg)
      case("--")
         getopts = .false.
      case("--help")
         write(output_unit, '(a)') help_text_tagdiff
         stop
      case("--version")
         call version(output_unit)
         stop
      case("--fit")
         config%fit = .true.
      case default
         if (.not.allocated(config%actual)) then
            call move_alloc(arg, config%actual)
            cycle
         end if
         if (.not.allocated(config%reference)) then
            call move_alloc(arg, config%reference)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      end select
   end do

   if (.not.(allocated(config%actual).and.allocated(config%reference))) then
      if (.not.allocated(error)) then
         write(output_unit, '(a)') help_text_tagdiff
         error stop
      end if
   end if
end subroutine get_tagdiff_arguments


subroutine get_argument_as_real(arg, val, error)
   !> Index of command line argument, range [0:command_argument_count()]
   character(len=:), intent(in), allocatable :: arg
   !> Real value
   real(wp), intent(out) :: val
   !> Error handling
   type(error_type), allocatable :: error

   integer :: stat

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


subroutine get_argument_as_realv(arg, val, error)
   !> Index of command line argument, range [0:command_argument_count()]
   character(len=:), intent(in), allocatable :: arg
   !> Real value
   real(wp), intent(out) :: val(:)
   !> Error handling
   type(error_type), allocatable :: error

   integer :: stat
   integer :: i
   character(len=*), parameter :: sep = ","
   character(len=:), allocatable :: targ

   if (.not.allocated(arg)) then
      call fatal_error(error, "Cannot read real value, argument missing")
      return
   end if
   allocate(character(len=len(arg)) :: targ)
   do i = 1, len(arg)
      if (arg(i:i) == sep) then
         targ(i:i) = " "
      else
         targ(i:i) = arg(i:i)
      end if
   end do
   read(targ, *, iostat=stat) val
   if (stat /= 0) then
      call fatal_error(error, "Cannot read real value from '"//arg//"'")
      return
   end if

end subroutine get_argument_as_realv


subroutine get_argument_as_int(arg, val, error)
   !> Index of command line argument, range [0:command_argument_count()]
   character(len=:), intent(in), allocatable :: arg
   !> Real value
   integer, intent(out) :: val
   !> Error handling
   type(error_type), allocatable :: error

   integer :: stat

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


subroutine version(unit)
   integer, intent(in) :: unit
   character(len=:), allocatable :: version_string

   call get_tblite_version(string=version_string)
   write(unit, '(a, *(1x, a))') &
      & prog_name, "version", version_string

end subroutine version


end module tblite_cli
