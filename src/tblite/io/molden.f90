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

!> @file tblite/io/molden.f90
!> Provides Molden wavefunction input and output support.

!> The Molden writer produces a Molden input with the following sections:
!>
!> ```
!> [Molden Format]
!> [Title]
!> [Cell] AU          non-standard extension for periodic systems
!> [Atoms] AU
!> [Core]             non-standard extension for core electron count
!> [Pseudo]           non-standard extension for valence electron count
!> [Nval]             non-standard extension for valence electron count
!> [GTO]
!> [6D]/[10F]/[15G]   optional explicit declaration of cartesian format
!> [MO]
!> ```
!>
!> This format extends the official Molden specification documented at
!>
!>   https://www.theochem.ru.nl/molden/molden_format.html
!>
!> The non-standard `[Cell]` section (consistent with Multiwfn and CP2K)
!> specifies periodic information lacking in the Molden Format, and
!> contains cartesian lattice vectors assuming three-dimensional periodicity:
!>
!> ```
!> [Cell] AU
!>   ax ay az
!>   bx by bz
!>   cx cy cz
!> ```
!>
!> The Molden writer and reader follow the official specification for the ordering
!> of cartesian/spherical atomic orbitals in the `[GTO]` section:
!>
!> ```
!> s    : s
!>
!> p    : px, py, pz
!>
!> 6D   : xx, yy, zz, xy, xz, yz
!> 5D   : 0, +1, -1, +2, -2
!>
!> 10F  : xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
!> 7F   : 0, +1, -1, +2, -2, +3, -3
!>
!> 15G  : xxxx, yyyy, zzzz, xxxy, xxxz, xyyy, yyyz, xzzz, yzzz,
!>        xxyy, xxzz, yyzz, xxyz, xyyz, xyzz
!> 9G   : 0, +1, -1, +2, -2, +3, -3, +4, -4
!> ```
!>
!> The writer always emits Cartesian AO coefficients in `[MO]` and declares
!> Cartesian higher angular momentum shells with `[6D]`, `[10F]`, and `[15G]`.
!> The reader accepts both Cartesian and spherical Molden MO coefficient ordering.
!>
!> Primitive coefficients written in the `[GTO]` section contain normalization
!> of the spherical harmonic Gaussian primitives, as well as normalization
!> of the contracted basis functions (the used expansion by Stewart of the STO-nG
!> functions inherently contains contracted normalization). The reader assumes
!> that the coefficients contain both primitive and contracted normalization.
!>
!> The number of valence electrons and the effective atomic charge is specified via
!> the non-standard `[Core]` section (consistent with Molden2AIM and PySCF),
!> the non-standard `[Pseudo]` section (consistent with Molden, CP2K, and Orca),
!> the non-standard `[Nval]` section (consistent with Multiwfn).
!>
!> In the `[MO]` section, restricted closed-shell wavefunctions are written as
!> spatial orbitals, while restricted open-shell and unrestricted wavefunctions
!> are written as separate alpha and beta spin orbitals.
!>
!> The Molden reader constructs tblite internal `structure_type`, `basis_type`, and
!> `wavefunction_type` objects from the corresponding sections of the Molden file.

!> Implementation of Molden file I/O
module tblite_io_molden
   use, intrinsic :: iso_fortran_env, only : iostat_end
   use mctc_env, only : error_type, fatal_error, wp, i8
   use mctc_io, only : structure_type, new
   use mctc_io_convert, only : aatoau
   use tblite_basis_type, only : basis_type, cgto_type, new_basis, new_cgto
   use tblite_version, only : get_tblite_version
   use tblite_wavefunction_spin, only : magnet_to_updown
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction, &
      & get_density_matrix
   implicit none
   private

   public :: load_molden, save_molden

   !> Maximum contraction length allowed in Molden basis sets
   integer, parameter :: maxg = 12

   !> Line buffer for intermediate storage of section content during parsing
   type :: line_buffer_type
      !> Actual line content
      character(len=512), allocatable :: lines(:)
      !> Current capacity of the buffer
      integer :: capacity = 32
      !> Number of lines currently stored
      integer :: count = 0
   contains
      !> Add line to the buffer with automatic resizing
      procedure :: push_back
   end type line_buffer_type

contains

!> Add line to the buffer with automatic resizing
subroutine push_back(self, line)
   !> Line buffer instance
   class(line_buffer_type), intent(inout) :: self
   !> Line to add to the buffer
   character(len=*), intent(in) :: line

   character(len=512), allocatable :: tmp(:)
   integer :: new_capacity

   ! Allocate lines upon first insertion
   if (.not. allocated(self%lines)) then
      allocate(self%lines(self%capacity))
      self%count = 0
   end if

   ! Double buffer size if needed
   if (self%count >= self%capacity) then
      new_capacity = max(1, 2 * self%capacity)

      allocate(tmp(new_capacity))
      tmp(1:self%count) = self%lines(1:self%count)

      call move_alloc(tmp, self%lines)
      self%capacity = new_capacity
   end if

   self%count = self%count + 1
   self%lines(self%count) = line
end subroutine push_back


!> Read Molden file and construct structure, basis, and wavefunction objects.
subroutine load_molden(filename, mol, bas, wfn, error)
   !> Input Molden file name
   character(len=*), intent(in) :: filename
   !> Molecular structure data loaded from Molden file
   type(structure_type), intent(out) :: mol
   !> Basis set information loaded from Molden file
   type(basis_type), intent(out) :: bas
   !> Wavefunction data loaded from Molden file
   type(wavefunction_type), intent(out) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=512) :: line, pending_line, errmsg
   character(len=16) :: spec, atoms_spec, cell_spec
   integer :: unit, stat
   integer(i8) :: mo_position
   type(line_buffer_type) :: cell, atoms, core, pseudo, nval, gto
   logical :: pending, has_mo, use_6d, use_10f, use_15g

   open(newunit=unit, file=filename, status="old", action="read", &
      & access="stream", form="formatted", iostat=stat, iomsg=errmsg)
   if (stat /= 0) then
      call fatal_error(error, "Failed to open Molden input '"//filename//"': "//trim(errmsg))
      return
   end if

   pending = .false.
   has_mo = .false.
   use_6d = .true.
   use_10f = .true.
   use_15g = .true.
   atoms_spec = ""
   cell_spec = ""
   do
      if (pending) then
         line = pending_line
         pending = .false.
      else
         ! Read the next line keeping track of the file position
         read(unit, "(A)", iostat=stat) line
         ! Check for EOF and I/O errors
         if (stat == iostat_end) exit
         if (stat /= 0) then
            call fatal_error(error, "I/O error reading: "//trim(line))
            exit
         end if
      end if

      ! Check for section header
      if (len_trim(section_id(line)) == 0) cycle

      ! Check for unit specificiations in the header line
      spec = section_spec(line)

      ! Read individual sections into line buffers excluding the section header
      select case(trim(section_id(line)))
      case("cell")
         cell_spec = spec
         call collect_section_lines(unit, cell, pending, pending_line, error)
      case("atoms")
         atoms_spec = spec
         call collect_section_lines(unit, atoms, pending, pending_line, error)
      case("gto")
         call collect_section_lines(unit, gto, pending, pending_line, error)
      case("core")
         call collect_section_lines(unit, core, pending, pending_line, error)
      case("pseudo")
         call collect_section_lines(unit, pseudo, pending, pending_line, error)
      case("nval")
         call collect_section_lines(unit, nval, pending, pending_line, error)
      case("mo")
         ! Record file position for on-the-fly parsing the MO section
         ! after the structure and basis set have been build
         has_mo = .true.
         inquire(unit=unit, pos=mo_position)
         call skip_section(unit, pending, pending_line, error)
      case("5d")
         use_6d = .false.
      case("6d")
         use_6d = .true.
      case("7f")
         use_10f = .false.
      case("10f")
         use_10f = .true.
      case("9g")
         use_15g = .false.
      case("15g")
         use_15g = .true.
      case default
         ! Skipping sections like [Molden Format], [Title], [Freq], and [Geom]
         call skip_section(unit, pending, pending_line, error)
      end select
      if (allocated(error)) exit
   end do
   ! Check for missing required [MO] section
   if (.not. has_mo) then
      call fatal_error(error, "Molden input is missing required [MO] section")
   end if
   if (allocated(error)) then
      close(unit, iostat=stat, iomsg=errmsg)
      if (stat /= 0 .and. .not.allocated(error)) then
         call fatal_error(error, "Failed to close Molden input '"//filename//"': "//trim(errmsg))
      end if
      return
   end if

   ! Construct structure and basis set with the buffered section content
   call build_structure(cell, cell_spec, atoms, atoms_spec, mol, error)
   if (.not.allocated(error)) call build_basis(gto, mol, bas, error)

   ! Read the wavefunction information now that the dimensions are known
   if (.not.allocated(error)) call build_wavefunction(unit, mo_position, &
      & core, pseudo, nval, use_6d, use_10f, use_15g, mol, bas, wfn, error)

   close(unit, iostat=stat, iomsg=errmsg)
   if (stat /= 0 .and. .not.allocated(error)) then
      call fatal_error(error, "Failed to close Molden input '"//filename//"': "//trim(errmsg))
      return
   end if
end subroutine load_molden

subroutine collect_section_lines(unit, buffer, pending, pending_line, error)
   !> Open Molden file unit
   integer, intent(in) :: unit
   !> Line buffer to collect the section content
   type(line_buffer_type), intent(inout) :: buffer
   !> Flag to indicate that the last read line was a pending section header
   logical, intent(inout) :: pending
   !> Storage for a pending section header line
   character(len=512), intent(inout) :: pending_line
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=512) :: line
   integer :: stat

   do
      read(unit, "(A)", iostat=stat) line
      ! Check for EOF and I/O errors
      if (stat == iostat_end) exit
      if (stat /= 0) then
         call fatal_error(error, "I/O error reading: "//trim(line))
         return
      end if
      ! Check for empty lines and skip them
      if (len_trim(line) == 0) cycle
      ! Check if current line is a section header
      if (len_trim(section_id(line)) > 0) then
         pending_line = line
         pending = .true.
         exit
      end if
      ! Store line in buffer
      call buffer%push_back(line)
   end do
end subroutine collect_section_lines

subroutine skip_section(unit, pending, pending_line, error)
   !> Open Molden file unit
   integer, intent(in) :: unit
   !> Flag to indicate that the last read line was a pending section header
   logical, intent(inout) :: pending
   !> Storage for a pending section header line
   character(len=512), intent(inout) :: pending_line
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=512) :: line
   integer :: stat

   do
      read(unit, "(A)", iostat=stat) line
      ! Check for EOF and I/O errors
      if (stat == iostat_end) exit
      if (stat /= 0) then
         call fatal_error(error, "I/O error reading: "//trim(line))
         return
      end if
      ! Check if current line is a section header
      if (len_trim(section_id(line)) > 0) then
         pending_line = line
         pending = .true.
         exit
      end if
   end do
end subroutine skip_section


!> Write tblite singlepoint data to a Molden file.
subroutine save_molden(filename, mol, bas, wfn, error)
   !> Output Molden file name
   character(len=*), intent(in) :: filename
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis metadata used by the calculation
   class(basis_type), intent(in) :: bas
   !> Converged wavefunction
   type(wavefunction_type), intent(in) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=512) :: errmsg
   integer :: unit, stat
   character(len=:), allocatable :: version_string

   open(newunit=unit, file=filename, status="replace", action="write", &
      & iostat=stat, iomsg=errmsg)
   if (stat /= 0) then
      call fatal_error(error, "Failed to open Molden output '"//filename//"': "//trim(errmsg))
      return
   end if

   write(unit,"(A)") "[Molden Format]"
   write(unit,"(A)") "[Title]"
   call get_tblite_version(string=version_string)
   write(unit,"(A)") "tblite version "//trim(version_string)

   call write_cell(unit, mol, error)
   if (.not.allocated(error)) call write_atoms(unit, mol, error)
   if (.not.allocated(error)) call write_core(unit, mol, wfn, error)
   if (.not.allocated(error)) call write_pseudo(unit, mol, wfn, error)
   if (.not.allocated(error)) call write_nval(unit, mol, wfn, error)
   if (.not.allocated(error)) call write_gto(unit, mol, bas, error)
   if (.not.allocated(error)) then
      if (bas%maxl >= 2) write(unit,"(A)") "[6D]"
      if (bas%maxl >= 3) write(unit,"(A)") "[10F]"
      if (bas%maxl >= 4) write(unit,"(A)") "[15G]"
   end if
   if (.not.allocated(error)) call write_mo(unit, mol, bas, wfn, error)

   close(unit, iostat=stat, iomsg=errmsg)
   if (stat /= 0 .and. .not.allocated(error)) then
      call fatal_error(error, "Failed to close Molden output '"//filename//"': "//trim(errmsg))
   end if
end subroutine save_molden


!> Build molecular structure from buffered [Cell] and [Atoms] sections.
subroutine build_structure(cell, cell_spec, atoms, atoms_spec, mol, error)
   !> Buffered [Cell] section
   type(line_buffer_type), intent(in) :: cell
   !> Unit specifier from [Cell] header
   character(len=*), intent(in) :: cell_spec
   !> Buffered [Atoms] section
   type(line_buffer_type), intent(in) :: atoms
   !> Unit specifier from [Atoms] header
   character(len=*), intent(in) :: atoms_spec
   !> Molecular structure data loaded from Molden file
   type(structure_type), intent(out) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: nat
   integer, allocatable :: atomic_number(:)
   real(wp), allocatable :: nuclear_coord(:, :), lattice(:, :)
   logical, allocatable :: periodic(:)

   call read_cell(cell, cell_spec, periodic, lattice, error)
   if (.not. allocated(error)) then
      call read_atoms(atoms, atoms_spec, nat, atomic_number, nuclear_coord, error)
   end if
   if (allocated(error)) return

   ! Set up new molecular structure object with dummy charge and unpaired electrons
   call new(mol, atomic_number, nuclear_coord, charge=0.0_wp, uhf=0, &
      & lattice=lattice, periodic=periodic)
end subroutine build_structure

!> Read optional [Cell] section for 3 lattice vectors from buffered lines.
subroutine read_cell(buffer, spec, periodic, lattice, error)
   !> Buffered [Cell] section
   type(line_buffer_type), intent(in) :: buffer
   !> Unit specifier from the [Cell] header
   character(len=*), intent(in) :: spec
   !> Periodic directions
   logical, allocatable, intent(out) :: periodic(:)
   !> Lattice vectors in Bohr
   real(wp), allocatable, intent(out) :: lattice(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=512) :: line
   integer :: ilat, iline, stat
   real(wp) :: conversion

   if (buffer%count == 0) return

   ! Read optional unit specifier from the header
   conversion = 1.0_wp
   if (len_trim(spec) > 0) then
      select case(trim(spec))
      case("au", "bohr")
         conversion = 1.0_wp
      case("angs")
         conversion = aatoau
      case default
         call fatal_error(error, "Unrecognized unit specifier in [Cell] header: "//trim(spec))
         return
      end select
   end if

   ! If there is a [Cell] section, we require three dimensional periodicity
   allocate(periodic(3), source=.true.)
   allocate(lattice(3, 3), source=0.0_wp)
   ilat = 0
   do iline = 1, buffer%count
      line = buffer%lines(iline)
      ilat = ilat + 1
      if (ilat > 3) then
         call fatal_error(error, "[Cell] has more than three lattice vectors")
         return
      end if
      read(line, *, iostat=stat) lattice(:, ilat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to read [Cell] vector: "//trim(line))
         return
      end if
      lattice(:, ilat) = conversion * lattice(:, ilat)
   end do

   if (ilat < 3) then
      call fatal_error(error, "[Cell] requires three lattice vectors")
      return
   end if
end subroutine read_cell

!> Read required [Atoms] section for atomic numbers and coordinates from buffered lines.
subroutine read_atoms(buffer, spec, nat, atomic_number, nuclear_coord, error)
   !> Buffered [Atoms] section
   type(line_buffer_type), intent(in) :: buffer
   !> Unit specifier from the [Atoms] header
   character(len=*), intent(in) :: spec
   !> Number of atoms
   integer, intent(out) :: nat
   !> Atomic numbers
   integer, allocatable, intent(out) :: atomic_number(:)
   !> Nuclear coordinates in Bohr
   real(wp), allocatable, intent(out) :: nuclear_coord(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=512) :: line
   character(len=16) :: sym
   integer :: iline, iat, stat, tmp_z
   real(wp) :: conversion, tmp_xyz(3)

   if (buffer%count == 0) then
      call fatal_error(error, "Molden input is missing required [Atoms] section")
      return
   end if

   ! Read optional unit specifier from the header
   conversion = 1.0_wp
   if (len_trim(spec) > 0) then
      select case(trim(spec))
      case("au", "bohr")
         conversion = 1.0_wp
      case("angs")
         conversion = aatoau
      case default
         call fatal_error(error, "Unrecognized unit specifier in [Atoms] header: "//trim(spec))
         return
      end select
   end if

   nat = buffer%count
   allocate(atomic_number(nat), nuclear_coord(3, nat))
   do iline = 1, buffer%count
      line = buffer%lines(iline)
      read(line, *, iostat=stat) sym, iat, tmp_z, tmp_xyz
      if (stat /= 0) then
         call fatal_error(error, "Failed to read [Atoms] entry: "//trim(line))
         return
      end if
      if (iat < 1 .or. iat > nat) then
         call fatal_error(error, "Atom index out of range in [Atoms] entry: "//trim(line))
         return
      end if
      atomic_number(iat) = tmp_z
      nuclear_coord(:, iat) = conversion * tmp_xyz
   end do
end subroutine read_atoms


!> Build basis set from buffered [GTO] section.
subroutine build_basis(gto, mol, bas, error)
   !> Buffered [GTO] section
   type(line_buffer_type), intent(in) :: gto
   !> Molecular structure data loaded from Molden file
   type(structure_type), intent(in) :: mol
   !> Basis set information loaded from molden File
   type(basis_type), intent(out) :: bas
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: nsh
   integer, allocatable :: nsh_id(:)
   type(cgto_type), allocatable :: cgto(:, :)

   call read_gto(gto, mol, nsh, nsh_id, cgto, error)
   if (allocated(error)) return

   ! Set up new basis set object with dummy accuracy setting
   call new_basis(bas, mol, nsh_id, cgto, 1.0_wp)
end subroutine build_basis

!> Read [GTO] section for CGTO information from buffered lines.
subroutine read_gto(buffer, mol, nsh, nsh_id, cgto, error)
   !> Buffered [GTO] section
   type(line_buffer_type), intent(in) :: buffer
   !> Molecular structure data loaded from Molden file
   type(structure_type), intent(in) :: mol
   !> Number of shells
   integer, intent(out) :: nsh
   !> Number of shells per species
   integer, allocatable, intent(out) :: nsh_id(:)
   !> CGTO for each shell and species
   type(cgto_type), allocatable, intent(out) :: cgto(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=512) :: line
   integer :: iat, isp, ish, iprim, l, nprim, iline, stat
   integer, allocatable :: nsh_at(:), nsh_work(:)
   real(wp) :: alpha(maxg), coeff(maxg)
   logical, allocatable :: seen_cgto(:, :)
   type(cgto_type) :: tmp_cgto

   if (buffer%count == 0) then
      call fatal_error(error, "Molden input is missing required [GTO] section")
      return
   end if

   ! First pass to count and map shells
   allocate(nsh_at(mol%nat), source=0)
   nsh = 0
   iat = 0
   iline = 1
   do while (iline <= buffer%count)
      line = buffer%lines(iline)
      ! Check for for atom index line
      call parse_atom_index(line, iat, stat)
      if (stat == 0) then
         if (iat < 1 .or. iat > mol%nat) then
            call fatal_error(error, "[GTO] atom index out of range: "//trim(line))
            return
         end if
         iline = iline + 1
         cycle
      end if
      ! Check for shell header
      call parse_shell_header(line, l, nprim, stat)
      if (stat == 0) then
         if (iat < 1 .or. iat > mol%nat) then
            call fatal_error(error, "[GTO] shell before valid atom index")
            return
         end if
         if (nprim <= 0 .or. nprim > maxg) then
            call fatal_error(error, "[GTO] shell has zero primitives or more primitives than supported")
            return
         end if
         nsh = nsh + 1
         nsh_at(iat) = nsh_at(iat) + 1
         iline = iline + nprim + 1
         cycle
      end if

      call fatal_error(error, "Failed to parse [GTO] line: "//trim(line))
      return
   end do

   ! Check for identical shell counts for all atoms of a given species
   allocate(nsh_id(mol%nid))
   nsh_id(:) = -1
   do iat = 1, mol%nat
      isp = mol%id(iat)
      if (nsh_id(isp) < 0) then
         nsh_id(isp) = nsh_at(iat)
      else if (nsh_id(isp) /= nsh_at(iat)) then
         call fatal_error(error, "Molden input has inconsistent shell counts for one species")
         return
      end if
   end do
   if (any(nsh_id <= 0)) then
      call fatal_error(error, "Molden input missing shells for at least one species")
      return
   end if

   ! Build CGTOs for each species
   allocate(cgto(maxval(nsh_id), mol%nid))
   allocate(seen_cgto(maxval(nsh_id), mol%nid), source=.false.)
   allocate(nsh_work(mol%nat), source=0)
   iline = 1
   do while (iline <= buffer%count)
      line = buffer%lines(iline)
      call parse_atom_index(line, iat, stat)
      if (stat == 0) then
         iline = iline + 1
         cycle
      end if
      call parse_shell_header(line, l, nprim, stat)
      if (stat /= 0) cycle

      ! Read primitive GTO information
      alpha(:) = 0.0_wp
      coeff(:) = 0.0_wp
      do iprim = 1, nprim
         iline = iline + 1
         line = buffer%lines(iline)
         read(line, *, iostat=stat) alpha(iprim), coeff(iprim)
         if (stat /= 0) then
            call fatal_error(error, "Failed to parse [GTO] primitive: "//trim(line))
            return
         end if
      end do

      ! Assume normalization is already included in the coefficients
      call new_cgto(tmp_cgto, nprim, l, alpha, coeff, .false.)

      ! Check if current CGTO matches earlier CGTOs for the same species
      isp = mol%id(iat)
      nsh_work(iat) = nsh_work(iat) + 1
      ish = nsh_work(iat)
      if (.not. seen_cgto(ish, isp)) then
         cgto(ish, isp) = tmp_cgto
         seen_cgto(ish, isp) = .true.
      else if (.not.cgto(ish, isp)%compare(tmp_cgto)) then
         call fatal_error(error, "Molden basis differs between atoms of the same species")
         return
      end if
      iline = iline + 1
   end do
end subroutine read_gto

subroutine parse_atom_index(line, iat, stat)
   !> Line to parse for an atom index
   character(len=*), intent(in) :: line
   !> Atom index if the line contains an atom header
   integer, intent(inout) :: iat
   !> Status of the parsing
   integer, intent(out) :: stat

   integer :: tmp_iat, dummy

   ! Check for Molden atom header: atom_index 0
   read(line, *, iostat=stat) tmp_iat, dummy
   if (stat == 0 .and. dummy == 0) then
      iat = tmp_iat
   else
      stat = 1
   end if
end subroutine parse_atom_index

subroutine parse_shell_header(line, l, nprim, stat)
   !> Line to parse for an shell header
   character(len=*), intent(in) :: line
   !> Angular momentum of the shell
   integer, intent(out) :: l
   !> Number of primitives in the shell
   integer, intent(out) :: nprim
   !> Status of the parsing
   integer, intent(out) :: stat

   character(len=16) :: label
   real(wp) :: scale

   l = -1
   nprim = 0
   scale = 1.0_wp

   ! Check for Molden shell header (ignoring scale): label nprim [scale]
   read(line, *, iostat=stat) label, nprim, scale
   if (stat /= 0) then
      read(line, *, iostat=stat) label, nprim
      scale = 1.0_wp
   end if
   if (stat /= 0) return

   select case(trim(lower_string(adjustl(label))))
   case("s")
      l = 0
   case("p")
      l = 1
   case("d")
      l = 2
   case("f")
      l = 3
   case("g")
      l = 4
   case default
      stat = 1
      return
   end select
end subroutine parse_shell_header


!> Build wavefunction from buffered [Core], [Pseudo], and [Nval] sections
!> while reading the [MO] section on-the-fly.
subroutine build_wavefunction(unit, mo_position, core, pseudo, nval, &
   & use_6d, use_10f, use_15g, mol, bas, wfn, error)
   !> Open Molden file unit
   integer, intent(in) :: unit
   !> Stream position of first line after [MO] header
   integer(i8), intent(in) :: mo_position
   !> Buffered [Core] section
   type(line_buffer_type), intent(in) :: core
   !> Buffered [Pseudo] section
   type(line_buffer_type), intent(in) :: pseudo
   !> Buffered [Nval] section
   type(line_buffer_type), intent(in) :: nval
   !> Use Cartesian d coefficients in [MO]
   logical, intent(in) :: use_6d
   !> Use Cartesian f coefficients in [MO]
   logical, intent(in) :: use_10f
   !> Use Cartesian g coefficients in [MO]
   logical, intent(in) :: use_15g
   !> Molecular structure, updated with charge/uhf
   type(structure_type), intent(inout) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Wavefunction data
   type(wavefunction_type), intent(out) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: spin, nspin
   integer, allocatable :: nvalence(:)
   real(wp), allocatable :: coeff(:, :, :), focc(:, :), emo(:,:), tmp_focc(:)
   logical :: cartesian

   ! Parse optional core electron sections checking for consistency.
   call read_core(core, mol, nvalence, error)
   if (.not. allocated(error)) call read_pseudo(pseudo, mol, nvalence, error)
   if (.not. allocated(error)) call read_nval(nval, mol, nvalence, error)
   if (.not. allocated(error) .and. allocated(nvalence)) then
      if (any(nvalence <= 0)) then
         call fatal_error(error, "No valence electrons assigned for at least one atom")
         return
      end if
   end if
   if (allocated(error)) return

   ! Select Cartesian (default) or spherical AO basis
   select case(bas%maxl)
   case(2)
      cartesian = use_6d
   case(3)
      cartesian = use_6d .and. use_10f
      if (.not.cartesian .and. (use_6d .or. use_10f)) then
         call fatal_error(error, "Molden file must use consistent Cartesian or spherical ordering")
         return
      end if
   case(4)
      cartesian = use_6d .and. use_10f .and. use_15g
      if (.not.cartesian .and. (use_6d .or. use_10f .or. use_15g)) then
         call fatal_error(error, "Molden file must use consistent Cartesian or spherical ordering")
         return
      end if
   case default
      cartesian = .false.
   end select

   ! Read the MO section directly from the file stream
   call read_mo(unit, mo_position, mol, bas, cartesian, nspin, coeff, focc, emo, error)
   if (allocated(error)) return

   ! Setup new wavefunction object (by default 0K electronic temperature)
   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, nspin, 0.0_wp)

   ! Set MO occupation and unpaired electrons
   wfn%nocc = sum(focc)
   wfn%nel(1) = sum(focc(:, 1))
   wfn%nel(2) = sum(focc(:, 2))
   wfn%nuhf = abs(wfn%nel(1) - wfn%nel(2))
   wfn%emo(:, :) = emo(:, 1:nspin)
   wfn%focc(:, :) = focc

   ! Set atomic reference occupation from valence electron count
   if (allocated(nvalence)) then
      wfn%n0at(:) = nvalence
   else
      wfn%n0at(:) = mol%num(mol%id)
   end if

   ! Update unpaired electrons and charge also for the structure type
   mol%uhf = nint(wfn%nuhf)
   mol%charge = sum(real(wfn%n0at, wp)) - sum(wfn%nel)

   ! Transform reordered MO coefficients to spherical basis if necessary
   if (cartesian) then
      do spin = 1, nspin
         call bas%cartesian_to_spherical_trafo(mol, coeff(:, :, spin), &
            & wfn%coeff(:, :, spin))
      end do
   else
      do spin = 1, nspin
         wfn%coeff(:, :, spin) = coeff(:, :, spin)
      end do
   end if

   ! Calculate the density matrix from MO coefficients and occupation numbers
   if (nspin == 1) then
      allocate(tmp_focc(size(wfn%focc, 1)))
      tmp_focc(:) = wfn%focc(:, 1) + wfn%focc(:, 2)
      call get_density_matrix(tmp_focc, wfn%coeff(:, :, 1), wfn%density(:, :, 1))
   else
      do spin = 1, nspin
         call get_density_matrix(wfn%focc(:, spin), wfn%coeff(:, :, spin), &
            & wfn%density(:, :, spin))
      end do
   end if
end subroutine build_wavefunction

!> Read optional [Core] section for valence electron count from buffered lines.
subroutine read_core(buffer, mol, nvalence, error)
   !> Buffered [Core] section
   type(line_buffer_type), intent(in) :: buffer
   !> Molecular structure data loaded from Molden file
   type(structure_type), intent(in) :: mol
   !> Number of valence electrons per atom
   integer, allocatable, intent(inout) :: nvalence(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=512) :: line
   character(len=16) :: sym, sep
   integer :: iline, iat, isp, stat, val, tmp_valence
   logical :: check, matched, by_index

   if (buffer%count == 0) return

   ! Check for consistency, if valence electrons were also read in a different section
   if (.not. allocated(nvalence)) then
      check = .false.
      allocate(nvalence(mol%nat), source=0)
   else
      check = .true.
   end if

   do iline = 1, buffer%count
      line = buffer%lines(iline)
      iat = -1
      by_index = .false.
      ! Read core entry: iat [:] val
      read(line, *, iostat=stat) iat, sep, val
      if (stat == 0) by_index = .true.
      if (stat /= 0) read(line, *, iostat=stat) iat, val
      if (stat == 0) by_index = .true.
      ! Read core entry: sym [:] val
      if (stat /= 0) read(line, *, iostat=stat) sym, sep, val
      if (stat /= 0) read(line, *, iostat=stat) sym, val
      if (stat /= 0) then
         call fatal_error(error, "Failed to read [Core] entry: "//trim(line))
         return
      end if

      if (by_index) then
         ! Core electrons provided for atom index
         if (iat < 1 .or. iat > mol%nat) then
            call fatal_error(error, "[Core] atom index out of range: "//trim(line))
            return
         end if
         isp = mol%id(iat)
         tmp_valence = mol%num(isp) - val
         if (check .and. nvalence(iat) /= tmp_valence) then
            call fatal_error(error, "Inconsistent core electron counts in [Core]")
            return
         end if
         nvalence(iat) = tmp_valence
         cycle
      else
         ! Core electrons provided for atom label
         matched = .false.
         do iat = 1, mol%nat
            isp = mol%id(iat)
            if (trim(lower_string(adjustl(sym))) == &
               & trim(lower_string(adjustl(mol%sym(isp))))) then
               matched = .true.
               tmp_valence = mol%num(isp) - val
               if (check .and. nvalence(iat) /= tmp_valence) then
                  call fatal_error(error, "Inconsistent core electron counts in [Core]")
                  return
               end if
               nvalence(iat) = tmp_valence
            end if
         end do
         if (.not. matched) then
            call fatal_error(error, "No matching atom for [Core] entry: "//trim(line))
            return
         end if
      end if
   end do
end subroutine read_core

!> Read optional [Pseudo] section for valence electron count from buffered lines.
subroutine read_pseudo(buffer, mol, nvalence, error)
   !> Buffered [Pseudo] section
   type(line_buffer_type), intent(in) :: buffer
   !> Molecular structure data loaded from Molden file
   type(structure_type), intent(in) :: mol
   !> Number of valence electrons per atom
   integer, allocatable, intent(inout) :: nvalence(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=512) :: line
   character(len=16) :: sym
   integer :: iline, iat, isp, stat, val
   logical :: check

   if (buffer%count == 0) return

   ! Check for consistency, if valence electrons were also read in a different section
   if (.not. allocated(nvalence)) then
      check = .false.
      allocate(nvalence(mol%nat), source=0)
   else
      check = .true.
   end if

   do iline = 1, buffer%count
      line = buffer%lines(iline)
      ! Read pseudo entry: sym iat val
      read(line, *, iostat=stat) sym, iat, val
      if (stat /= 0) then
         call fatal_error(error, "Failed to read [Pseudo] entry: "//trim(line))
         return
      end if
      if (iat < 1 .or. iat > mol%nat) then
         call fatal_error(error, "[Pseudo] atom index out of range: "//trim(line))
         return
      end if
      ! Valence electrons provided for atom index and atom label
      isp = mol%id(iat)
      if (trim(lower_string(adjustl(sym))) /= &
         & trim(lower_string(adjustl(mol%sym(isp))))) then
         call fatal_error(error, "Atom index and symbol do not match in [Pseudo] entry: "//trim(line))
         return
      end if
      if (check .and. nvalence(iat) /= val) then
         call fatal_error(error, "Inconsistent valence electron counts in [Pseudo]")
         return
      end if
      nvalence(iat) = val
   end do
end subroutine read_pseudo

!> Read optional [Nval] section for valence electron count from buffered lines.
subroutine read_nval(buffer, mol, nvalence, error)
   !> Buffered [Nval] section
   type(line_buffer_type), intent(in) :: buffer
   !> Molecular structure data loaded from Molden file
   type(structure_type), intent(in) :: mol
   !> Number of valence electrons per atom
   integer, allocatable, intent(inout) :: nvalence(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=512) :: line
   character(len=16) :: sym
   integer :: iline, iat, isp, stat, val
   logical :: check, matched

   if (buffer%count == 0) return

   ! Check for consistency, if valence electrons were also read in a different section
   if (.not. allocated(nvalence)) then
      check = .false.
      allocate(nvalence(mol%nat), source=0)
   else
      check = .true.
   end if

   do iline = 1, buffer%count
      line = buffer%lines(iline)
      ! Read nval entry: sym val
      read(line, *, iostat=stat) sym, val
      if (stat /= 0) then
         call fatal_error(error, "Failed to read [Nval] entry: "//trim(line))
         return
      end if
      ! Valence electrons provided for atom label
      matched = .false.
      do iat = 1, mol%nat
         isp = mol%id(iat)
         if (trim(lower_string(adjustl(sym))) /= &
            & trim(lower_string(adjustl(mol%sym(isp))))) then
            cycle
         end if
         matched = .true.
         if (check .and. nvalence(iat) /= val) then
            call fatal_error(error, "Inconsistent valence electron counts in [Nval]")
            return
         end if
         nvalence(iat) = val
      end do
      if (.not. matched) then
         call fatal_error(error, "No matching atom for [Nval] entry: "//trim(line))
         return
      end if
   end do
end subroutine read_nval

!> Read [MO] section for MO coefficients, energies, and occupations.
subroutine read_mo(unit, mo_position, mol, bas, cartesian, &
   & nspin, coeff, focc, emo, error)
   !> Open Molden file unit
   integer, intent(in) :: unit
   !> Stream position of first line after [MO] header
   integer(i8), intent(in) :: mo_position
   !> Molecular structure data loaded from Molden file
   type(structure_type), intent(in) :: mol
   !> Basis set information loaded from Molden File
   type(basis_type), intent(in) :: bas
   !> Use Cartesian (default) or spherical AO basis
   logical, intent(in) :: cartesian
   !> Spin channel count found in [MO]
   integer, intent(out) :: nspin
   !> MO coefficients, shape coeff(nao, nmo, 2)
   real(wp), allocatable, intent(out) :: coeff(:, :, :)
   !> MO occupation numbers, shape focc(nmo, 2)
   real(wp), allocatable, intent(out) :: focc(:, :)
   !> MO energies, shape emo(nmo, 2)
   real(wp), allocatable, intent(out) :: emo(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=512) :: line
   character(len=16) :: key
   integer :: stat, iat, isp, is, ish, ii, iicart, iao, imo, nao, nmo, l, tmp_nao
   integer :: idx_eq, tmp_spin, spin_count(2), perm_cart(15), perm_sphr(9)
   integer, allocatable :: ao_map(:)
   real(wp) :: tmp, tmp_energy, tmp_occup
   real(wp), allocatable :: tmp_coeff(:)
   logical :: have_block

   if (cartesian) then
      nao = bas%nao_cart
   else
      nao = bas%nao
   end if
   nmo = bas%nao
   nspin = 1

   ! Build map from Molden AO index to tblite AO index.
   allocate(ao_map(nao))
   do iat = 1, mol%nat
      isp = mol%id(iat)
      is = bas%ish_at(iat)
      do ish = 1, bas%nsh_id(isp)
         ii = bas%iao_sh(is + ish)
         iicart = bas%iao_cart_sh(is + ish)
         l = bas%cgto(ish, isp)%ang
         if (cartesian) then
            call get_molden_to_tblite_cart_perm(l, perm_cart, error)
            if (allocated(error)) return
            do iao = 1, bas%nao_cart_sh(is + ish)
               ao_map(iicart + iao) = iicart + perm_cart(iao)
            end do
         else
            call get_molden_to_tblite_sphr_perm(l, perm_sphr, error)
            if (allocated(error)) return
            do iao = 1, bas%nao_sh(is + ish)
               ao_map(ii + iao) = ii + perm_sphr(iao)
            end do
         end if
      end do
   end do

   ! Since the number of spin channels is unknown we allocated the maximum.
   allocate(coeff(nao, nmo, 2), focc(nmo, 2), emo(nmo, 2), source=0.0_wp)

   allocate(tmp_coeff(nao), source=0.0_wp)
   spin_count(:) = 0
   have_block = .false.
   tmp_spin = 1
   tmp_energy = 0.0_wp
   tmp_occup = 0.0_wp
   tmp_nao = 0

   ! Read the [MO] section starting from the first line after the header.
   read(unit, "(A)", pos=mo_position, iostat=stat) line
   do
      ! Check for end of file
      if (stat == iostat_end) exit
      ! Check for other I/O errors
      if (stat /= 0) then
         call fatal_error(error, "I/O error reading [MO]")
         return
      end if
      ! Check for start of new section
      if (len_trim(section_id(line)) > 0) exit
      ! Ignore blank lines.
      if (len_trim(line) == 0) then
         ! allow(C181): stat is checked at the top of the enclosing loop after `cycle`
         read(unit, "(A)", iostat=stat) line
         cycle
      end if

      ! Check for start of new MO and store previous block
      if (is_mo_key(line, "sym")) then
         if (have_block) then
            spin_count(tmp_spin) = spin_count(tmp_spin) + 1
            imo = spin_count(tmp_spin)
            if (imo > nmo) then
               call fatal_error(error, "[MO] has too many orbitals for one spin channel")
               return
            end if
            if (tmp_nao /= nao) then
               call fatal_error(error, "[MO] has inconsistent number of coefficients for one orbital")
               return
            end if
            ! Store temporary data in the final position
            coeff(:, imo, tmp_spin) = tmp_coeff(:)
            focc(imo, tmp_spin) = tmp_occup
            emo(imo, tmp_spin) = tmp_energy
         end if
         ! Reset state for the new MO block
         have_block = .true.
         tmp_spin = 1
         tmp_energy = 0.0_wp
         tmp_occup = 0.0_wp
         tmp_coeff(:) = 0.0_wp
         tmp_nao = 0

      ! Check for MO energy
      else if (is_mo_key(line, "ene")) then
         if (.not. have_block) then
            call fatal_error(error, "[MO] energy before orbital header")
            return
         end if
         idx_eq = index(line, "=")
         read(line(idx_eq+1:), *, iostat=stat) tmp_energy
         if (stat /= 0) then
            call fatal_error(error, "Failed to read MO energy: "//trim(line))
            return
         end if
      ! Check for MO spin
      else if (is_mo_key(line, "spin")) then
         if (.not. have_block) then
            call fatal_error(error, "[MO] spin before orbital header")
            return
         end if
         idx_eq = index(line, "=")
         key = ""
         if (idx_eq > 0 .and. idx_eq < len_trim(line)) then
            key = adjustl(line(idx_eq+1:len_trim(line)))
         end if
         if (index(lower_string(key), "beta") > 0) then
            tmp_spin = 2
            nspin = 2
         else if (index(lower_string(key), "alpha") > 0) then
            tmp_spin = 1
         else
            call fatal_error(error, "Failed to parse MO spin: "//trim(line))
            return
         end if
      ! Check for MO occupation
      else if (is_mo_key(line, "occup")) then
         if (.not. have_block) then
            call fatal_error(error, "[MO] occupation before orbital header")
            return
         end if
         idx_eq = index(line, "=")
         read(line(idx_eq+1:), *, iostat=stat) tmp_occup
         if (stat /= 0) then
            call fatal_error(error, "Failed to read MO occupation: "//trim(line))
            return
         end if
      ! Check for MO coefficients
      else
         if (.not. have_block) then
            call fatal_error(error, "[MO] coefficient before orbital header")
            return
         end if

         read(line, *, iostat=stat) iao, tmp
         if (stat /= 0) then
            call fatal_error(error, "Failed to read MO coefficient: "//trim(line))
            return
         end if
         if (iao < 1 .or. iao > nao) then
            call fatal_error(error, "[MO] coefficient index out of range")
            return
         end if
         tmp_nao = tmp_nao + 1
         tmp_coeff(ao_map(iao)) = tmp
      end if
      ! Read the next line
      ! allow(C181): stat is checked at the top of the enclosing loop on the next iteration
      read(unit, "(A)", iostat=stat) line
   end do

   ! Store the last block if it exists and check for consistent orbital counts
   if (have_block) then
      spin_count(tmp_spin) = spin_count(tmp_spin) + 1
      imo = spin_count(tmp_spin)
      if (imo > nmo) then
         call fatal_error(error, "[MO] has too many orbitals for one spin channel")
         return
      end if
      if (tmp_nao /= nao) then
         call fatal_error(error, "[MO] has inconsistent number of coefficients for one orbital")
         return
      end if
      ! Store temporary data in the final position
      coeff(:, imo, tmp_spin) = tmp_coeff(:)
      focc(imo, tmp_spin) = tmp_occup
      emo(imo, tmp_spin) = tmp_energy
   end if

   ! Transform restricted occupation into separate alpha and beta occupation
   if (nspin == 1) then
      call magnet_to_updown(focc)
   end if

   if (spin_count(1) /= nmo) then
      call fatal_error(error, "[MO] has inconsistent alpha orbital count")
      return
   end if
   if (nspin == 2 .and. spin_count(2) /= nmo) then
      call fatal_error(error, "[MO] has inconsistent beta orbital count")
      return
   end if
end subroutine read_mo

!> Spherical AO order permutation from Molden (0, ... l, -l) to tblite (-l, ..., 0, ..., l)
subroutine get_molden_to_tblite_sphr_perm(l, perm, error)
   !> Angular momentum quantum number
   integer, intent(in)  :: l
   !> Permutation array
   integer, intent(out) :: perm(9)
   !> Error handling
   type(error_type), intent(out), allocatable :: error

   perm = 0
   select case(l)
   case(0)
      perm(1) = 1
   case(1)
      ! Exception in the order following the Cartesian convention
      perm(1:3) = [3, 1, 2]
   case(2)
      perm(1:5) = [3, 4, 2, 5, 1]
   case(3)
      perm(1:7) = [4, 5, 3, 6, 2, 7, 1]
   case(4)
      perm(1:9) = [5, 6, 4, 7, 3, 8, 2, 9, 1]
   case default
      call fatal_error(error, "Molden reader only supports angular momenta up to g")
      return
   end select
end subroutine get_molden_to_tblite_sphr_perm

!> Cartesian AO order permutation from Molden to tblite.
subroutine get_molden_to_tblite_cart_perm(l, perm, error)
   !> Angular momentum of the shell
   integer, intent(in) :: l
   !> Index permutation
   integer, intent(out) :: perm(15)
   !> Error handling
   type(error_type), intent(out), allocatable :: error

   perm = 0

   select case(l)
   case(0)
      perm(1) = 1
   case(1)
      perm(1:3) = [3, 1, 2]
   case(2)
      perm(1:6) = [1, 2, 3, 4, 5, 6]
   case(3)
      perm(1:10) = [1, 2, 3, 6, 4, 5, 8, 9, 7, 10]
   case(4)
      perm(1:15) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
   case default
      call fatal_error(error, "Molden reader only supports angular momenta up to g")
      return
   end select
end subroutine get_molden_to_tblite_cart_perm


!> Write optional [Cell] section in atomic units
subroutine write_cell(unit, mol, error)
   !> Output Molden file unit
   integer, intent(in) :: unit
   !> Molecular structure
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ilat

   if (.not. allocated(mol%periodic) .or. .not. allocated(mol%lattice)) return

   if (any(mol%periodic)) then
      write(unit,"(A)") "[Cell] AU"
      do ilat = 1, 3
         write(unit,"(3(1x,ES24.16))") mol%lattice(:, ilat)
      end do
   end if
end subroutine write_cell

!> Write [Atoms] section in atomic units
subroutine write_atoms(unit, mol, error)
   !> Output Molden file unit
   integer, intent(in) :: unit
   !> Molecular structure
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=2) :: sym
   integer :: iat, isp

   write(unit,"(A)") "[Atoms] AU"
   do iat = 1, mol%nat
      isp = mol%id(iat)
      sym = trim(mol%sym(isp))
      write(unit,"(a2,1x,i6,1x,i6,3(1x,ES24.16))") sym, iat, mol%num(isp), &
         & mol%xyz(1, iat), mol%xyz(2, iat), mol%xyz(3, iat)
   end do
end subroutine write_atoms

!> Write [Core] section: core electrons per species.
subroutine write_core(unit, mol, wfn, error)
   !> Output Molden file unit
   integer, intent(in) :: unit
   !> Molecular structure
   type(structure_type), intent(in) :: mol
   !> Converged wavefunction
   type(wavefunction_type), intent(in) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=2) :: sym
   integer :: iat, isp

   write(unit,"(A)") "[Core]"
   do isp = 1, mol%nid
      sym = trim(mol%sym(isp))
      do iat = 1, mol%nat
         if (mol%id(iat) == isp) exit
      end do
      write(unit,"(a2,1x,a,1x,i0)") sym, ":", mol%num(isp) - nint(wfn%n0at(iat))
   end do
end subroutine write_core

!> Write [Pseudo] section: valence electrons per atom.
subroutine write_pseudo(unit, mol, wfn, error)
   !> Output Molden file unit
   integer, intent(in) :: unit
   !> Molecular structure
   type(structure_type), intent(in) :: mol
   !> Converged wavefunction
   type(wavefunction_type), intent(in) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=2) :: sym
   integer :: iat, isp

   write(unit,"(A)") "[Pseudo]"
   do iat = 1, mol%nat
      isp = mol%id(iat)
      sym = trim(mol%sym(isp))
      write(unit,"(a2,1x,i6,1x,i6)") sym, iat, nint(wfn%n0at(iat))
   end do
end subroutine write_pseudo

!> Write [Nval] section: valence electrons per species.
subroutine write_nval(unit, mol, wfn, error)
   !> Output Molden file unit
   integer, intent(in) :: unit
   !> Molecular structure
   type(structure_type), intent(in) :: mol
   !> Converged wavefunction
   type(wavefunction_type), intent(in) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=2) :: sym
   integer :: iat, isp

   write(unit,"(A)") "[Nval]"
   do isp = 1, mol%nid
      sym = trim(mol%sym(isp))
      do iat = 1, mol%nat
         if (mol%id(iat) == isp) exit
      end do
      write(unit,"(a2,1x,i6)") sym, nint(wfn%n0at(iat))
   end do
end subroutine write_nval

!> Write [GTO] section.
subroutine write_gto(unit, mol, bas, error)
   !> Output Molden file unit
   integer, intent(in) :: unit
   !> Molecular structure
   type(structure_type), intent(in) :: mol
   !> Basis metadata
   class(basis_type), intent(in) :: bas
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=1) :: aang
   integer :: iat, iprim, ish

   write(unit,"(A)") "[GTO]"
   do iat = 1, mol%nat
      write(unit,"(i6,a)") iat, " 0"
      do ish = 1, bas%nsh_at(iat)
         associate(p_cgto => &
               & bas%cgto(ish, mol%id(iat)))
            select case(p_cgto%ang)
            case(0)
               aang = "s"
            case(1)
               aang = "p"
            case(2)
               aang = "d"
            case(3)
               aang = "f"
            case(4)
               aang = "g"
            case default
               call fatal_error(error, "Molden writer only supports angular momenta up to g")
               return
            end select

            write(unit,"(a,1x,i6,1x,f8.2)") aang, p_cgto%nprim, 1.00_wp
            do iprim = 1, p_cgto%nprim
               write(unit,"(2(1x,ES24.16))") p_cgto%alpha(iprim), p_cgto%coeff(iprim)
            end do
         end associate
      end do
      write(unit,*)
   end do
end subroutine write_gto

!> Write [MO] section with spherical-to-Cartesian transformation.
subroutine write_mo(unit, mol, bas, wfn, error)
   !> Output Molden file unit
   integer, intent(in) :: unit
   !> Molecular structure
   type(structure_type), intent(in) :: mol
   !> Basis metadata
   class(basis_type), intent(in) :: bas
   !> Converged wavefunction
   type(wavefunction_type), intent(in) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iao_cart, iat, isp, ii, imo, is, ish, li, jao
   integer :: spin, cspin, nspin_mo, perm(15)
   logical :: spin_resolved_mo
   real(wp), allocatable :: coeff_cart(:, :, :)

   ! Represent only restricted closed-shell MOs as one spin-summed MO list
   ! and both restricted open-shell and unrestricted MOs with two spin channels
   spin_resolved_mo = wfn%nspin > 1 .or. abs(wfn%nuhf) > epsilon(1.0_wp)
   if (spin_resolved_mo) then
      nspin_mo = 2
   else
      nspin_mo = 1
   end if

   ! Transform MO coefficients to the Cartesian AO basis
   allocate(coeff_cart(bas%nao_cart, bas%nao, wfn%nspin))
   do spin = 1, wfn%nspin
      call bas%spherical_to_cartesian_trafo(mol, wfn%coeff(:, :, spin), &
         & coeff_cart(:, :, spin))
   end do

   write(unit,"(A)") "[MO]"
   do spin = 1, nspin_mo
      cspin = min(spin, wfn%nspin)

      do imo = 1, bas%nao
         write(unit,"(A)", advance="no") "Sym= "
         if (spin_resolved_mo) then
            if (spin == 1) then
               write(unit,"(i5,a)") imo, "a (alpha)"
            else
               write(unit,"(i5,a)") imo, "a (beta)"
            end if
         else
            write(unit,"(i5,a)") imo, "a"
         end if

         write(unit,"(A)", advance="no") "Ene= "
         write(unit,*) wfn%emo(imo, cspin)

         write(unit,"(A)", advance="no") "Spin= "
         if (spin_resolved_mo .and. spin == 2) then
            write(unit,"(A)") "Beta"
         else
            write(unit,"(A)") "Alpha"
         end if

         write(unit,"(A)", advance="no") "Occup= "
         if (spin_resolved_mo) then
            write(unit,"(F14.8)") wfn%focc(imo, spin)
         else
            write(unit,"(F14.8)") wfn%focc(imo, 1) + wfn%focc(imo, 2)
         end if

         jao = 0
         do iat = 1, mol%nat
            isp = mol%id(iat)
            is = bas%ish_at(iat)
            do ish = 1, bas%nsh_at(iat)
               li = bas%cgto(ish, isp)%ang
               ii = bas%iao_cart_sh(is + ish)

               call get_molden_to_tblite_cart_perm(li, perm, error)
               if (allocated(error)) return

               do iao_cart = 1, bas%nao_cart_sh(is + ish)
                  jao = jao + 1
                  write(unit,"(i6,1x,ES24.16)") jao, coeff_cart(ii + perm(iao_cart), &
                     & imo, cspin)
               end do
            end do
         end do
      end do
   end do
end subroutine write_mo

!> Extract section name
pure function section_id(line) result(id)
   !> Input line
   character(len=*), intent(in) :: line
   !> Lower-case section id without brackets
   character(len=32) :: id

   character(len=len(line)) :: work
   integer :: right

   id = ""
   work = adjustl(line)
   if (len_trim(work) < 3) return
   if (work(1:1) /= "[") return

   right = index(work, "]")
   if (right <= 2) return
   id = lower_string(adjustl(work(2:right-1)))
end function section_id

!> Extract additional section specification
pure function section_spec(line) result(spec)
   !> Input line
   character(len=*), intent(in) :: line
   !> Specification after section id
   character(len=32) :: spec

   character(len=len(line)) :: work
   integer :: right, stat

   spec = ""
   work = adjustl(line)
   right = index(work, "]")
   if (right <= 0 .or. right >= len_trim(work)) return
   spec = adjustl(work(right+1:))
   ! allow(C181): a failed list-directed parse intentionally leaves spec unchanged
   read(spec, *, iostat=stat) spec
   spec = lower_string(spec)
end function section_spec

!> Convert a string to lowercase.
pure function lower_string(input) result(output)
   !> Input string
   character(len=*), intent(in) :: input
   !> Lower-case string
   character(len=len(input)) :: output

   integer :: i, ic

   output = input
   do i = 1, len(input)
      ic = iachar(input(i:i))
      if (ic >= iachar("A") .and. ic <= iachar("Z")) then
         output(i:i) = achar(ic + iachar("a") - iachar("A"))
      end if
   end do
end function lower_string

!> Check if a line has a given MO key before '='.
pure function is_mo_key(line, key) result(is_key)
   !> Input line
   character(len=*), intent(in) :: line
   !> Key to check for
   character(len=*), intent(in) :: key
   !> Whether the line matches
   logical :: is_key

   character(len=len(line)) :: left, work
   integer :: eq

   work = adjustl(line)
   eq = index(work, "=")
   if (eq <= 0) then
      is_key = .false.
      return
   end if
   left = lower_string(adjustl(work(:eq-1)))
   is_key = trim(left) == lower_string(key)
end function is_mo_key

end module tblite_io_molden
