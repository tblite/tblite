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

!> Module to collect all help printout for the tblite program
module tblite_cli_help
   implicit none
   private

   public :: prog_name, help_text, help_text_run, help_text_param

   character(len=*), parameter :: nl = new_line('a')

   !> Name of the main executable
   character(len=*), parameter :: prog_name = "tblite"

   !> General help text for version and help commands
   character(len=*), parameter :: help_text_general = &
      "      --version           Print program version and exit"//nl//&
      "      --help              Show this help message"

   !> Help text for main driver
   character(len=*), parameter :: help_text = &
      "Usage: "//prog_name//" <run|param> [options]"//nl//&
      ""//nl//&
      "Commands"//nl//&
      ""//nl//&
      "  run       Evaluate tight binding module on the provided input structure."//nl//&
      "            If no command is specified run is selected by default."//nl//&
      ""//nl//&
      "  param     Inspect and manipulate tight binding parametrization data."//nl//&
      ""//nl//&
      "Options"//nl//&
      ""//nl//&
      help_text_general//nl//&
      ""

   !> Help text for param command
   character(len=*), parameter :: help_text_param = &
      "Usage: "//prog_name//" param [options] <input>"//nl//&
      ""//nl//&
      "Import, export and manipulate tight-binding parameter files."//nl//&
      ""//nl//&
      "Options"//nl//&
      ""//nl//&
      "      --method <name>     Base parametrization of the xTB Hamiltonian to use"//nl//&
      "                          Available methods: gfn1, gfn2, ipea1"//nl//&
      "      --output <file>     Output file for writing parametrization"//nl//&
      help_text_general//nl//&
      ""

   !> Help text for run command
   character(len=*), parameter :: help_text_run = &
      "Usage: "//prog_name//" [run] [options] <input>"//nl//&
      ""//nl//&
      "Evaluates the tight binding model on the provided input structure."//nl//&
      "Periodic calculations are performed automatically for periodic input formats."//nl//&
      ""//nl//&
      "Supported geometry input formats are:"//nl//&
      ""//nl//&
      "- Xmol/xyz files (xyz, log)"//nl//&
      "- Turbomole's coord, riper's periodic coord (tmol, coord)"//nl//&
      "- DFTB+ genFormat geometry inputs as cluster, supercell or fractional (gen)"//nl//&
      "- VASP's POSCAR/CONTCAR input files (vasp, poscar, contcar)"//nl//&
      "- Protein Database files, only single files (pdb)"//nl//&
      "- Connection table files, molfile (mol) and structure data format (sdf)"//nl//&
      "- Gaussian's external program input (ein)"//nl//&
      ""//nl//&
      "Options"//nl//&
      ""//nl//&
      "  -c, --charge <real>     Set charge to molecule, overwrites .CHRG file"//nl//&
      "      --spin <int>        Set number of unpaired electrons, overwrites .UHF file"//nl//&
      "      --method <name>     Parametrization of the xTB Hamiltonian to use"//nl//&
      "                          Available methods: gfn1, gfn2, ipea1 (Default: gfn2)"//nl//&
      "      --param <file>      Parametrization file to use for calculation"//nl//&
      "      --grad [file]       Evaluate molecular gradient and virial"//nl//&
      "                          Results are stored in file (default: tblite.txt)"//nl//&
      "      --json [file]       Dump results as JSON output (default: tblite.json)"//nl//&
      "  -i, --input <format>    Hint for the format of the input file"//nl//&
      help_text_general//nl//&
      ""

end module tblite_cli_help
