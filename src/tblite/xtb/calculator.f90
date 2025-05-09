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

!> @file tblite/xtb/calculator.f90
!> Provides the calculator type for holding xTB Hamiltonian parametrization.

!> Implementation of calculator type for the extended-tight binding Hamiltonian.
!> The #tblite_xtb_calculator::xtb_calculator collects the basic interactions
!> required to perform a tight-binding calculation.
module tblite_xtb_calculator
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_basis_ortho, only : orthogonalize
   use tblite_basis_type, only : basis_type, new_basis, cgto_type
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_classical_halogen, only : halogen_correction, new_halogen_correction
   use tblite_container, only : container_type, container_list
   use tblite_coulomb_charge, only : coulomb_kernel, new_gamma_coulomb, gamma_coulomb, &
      & new_effective_coulomb, effective_coulomb, average_interface, &
      & harmonic_average, arithmetic_average, geometric_average
   use tblite_coulomb_multipole, only : new_damped_multipole
   use tblite_coulomb_thirdorder, only : new_onsite_thirdorder
   use tblite_disp, only : dispersion_type, d4_dispersion, new_d4_dispersion, &
      & d3_dispersion, new_d3_dispersion
   use tblite_ncoord, only : ncoord_type, new_ncoord
   use tblite_param, only : param_record
   use tblite_repulsion, only : new_repulsion
   use tblite_repulsion_effective, only : tb_repulsion
   use tblite_xtb_coulomb, only : tb_coulomb
   use tblite_xtb_h0, only : tb_hamiltonian, new_hamiltonian
   use tblite_xtb_spec, only : tb_h0spec
   implicit none
   private

   public :: new_xtb_calculator
   public :: param_h0spec

   !> Default value for self-consistent iteration mixing
   real(wp), parameter :: mixer_damping_default = 0.4_wp

   !> Default maximum number of self-consistent iterations
   integer, parameter :: max_iter_default = 250

   !> Extended tight-binding calculator
   type, public :: xtb_calculator
      !> Basis set definition
      type(basis_type) :: bas
      !> Core Hamiltonian
      type(tb_hamiltonian) :: h0
      !> Coordination number for modifying the self-energies
      class(ncoord_type), allocatable :: ncoord
      !> Electronegativity-weighted coordination number for modifying the self-energies
      class(ncoord_type), allocatable :: ncoord_en
      !> Repulsion energy interactions
      type(tb_repulsion), allocatable :: repulsion
      !> Collection of all Coulombic interactions
      type(tb_coulomb), allocatable :: coulomb
      !> Halogen bonding correction
      type(halogen_correction), allocatable :: halogen
      !> London-dispersion interaction
      class(dispersion_type), allocatable :: dispersion
      !> Parameter for self-consistent iteration mixing
      real(wp) :: mixer_damping = mixer_damping_default
      !> Maximum number of self-consistent iteractions
      integer :: max_iter = max_iter_default
      !> Store calculated integral intermediates
      logical :: save_integrals = .false.
      !> List of additional interaction containers
      type(container_list), allocatable :: interactions
      !> string with method or "custom"
      character(len=:), allocatable :: method
   contains
      !> Get information about density dependent quantities used in the energy
      procedure :: variable_info
      !> Information on calculator
      procedure :: info
      !> Add an interaction container
      procedure :: push_back
      !> Remove an interaction container
      procedure :: pop
   end type xtb_calculator


   !> Specification of the Hamiltonian
   type, extends(tb_h0spec) :: param_h0spec
      type(param_record), pointer :: param => null()
      integer, pointer :: irc(:) => null()
      logical, allocatable :: valence(:, :)
   contains
      !> Generator for the self energy / atomic levels of the Hamiltonian
      procedure :: get_selfenergy
      !> Generator for the coordination number dependent shift of the self energy
      procedure :: get_cnshift
      !> Generator for the enhancement factor to for scaling Hamiltonian elements
      procedure :: get_hscale
      !> Generator for the polynomial parameters for the distant dependent scaling
      procedure :: get_shpoly
      !> Generator for the reference occupation numbers of the atoms
      procedure :: get_reference_occ
   end type param_h0spec

   !> Constructor for Hamiltonian specification
   interface param_h0spec
      module procedure :: new_param_h0spec
   end interface param_h0spec


contains


!> Create new xTB Hamiltonian calculator from parametrization data
subroutine new_xtb_calculator(calc, mol, param, error)
   !> Instance of the xTB calculator
   type(xtb_calculator), intent(out) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: isp
   integer, allocatable :: irc(:)

   allocate(irc(mol%nid))

   do isp = 1, mol%nid
      call param%get(mol%sym(isp), mol%num(isp), irc(isp))
      if (irc(isp) == 0) then
         call fatal_error(error, "No entry in parametrization for element "//mol%sym(isp))
         exit
      end if
   end do
   if (allocated(error)) return

   call add_basis(calc, mol, param, irc)
   call add_ncoord(calc, mol, param)
   call add_hamiltonian(calc, mol, param, irc)
   call add_repulsion(calc, mol, param, irc)
   call add_halogen(calc, mol, param, irc)
   call add_dispersion(calc, mol, param)
   call add_coulomb(calc, mol, param, irc)

   calc%method = "custom"

end subroutine new_xtb_calculator


subroutine add_basis(calc, mol, param, irc)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)

   integer :: isp, ir, ish, stat, ng, il
   integer, allocatable :: nsh_id(:)
   integer :: ang_idx(0:4), ortho(10)
   type(cgto_type), allocatable :: cgto(:, :)

   nsh_id = param%record(irc)%nsh
   allocate(cgto(maxval(nsh_id), mol%nid))
   do isp = 1, mol%nid
      ang_idx = 0
      ortho = 0
      ir = irc(isp)
      do ish = 1, nsh_id(isp)
         il = param%record(ir)%lsh(ish)
         ng = param%record(ir)%ngauss(ish)
         if (ang_idx(il) > 0) then
            ortho(ish) = ang_idx(il)
         else
            ang_idx(il) = ish
         end if
         call slater_to_gauss(ng, param%record(ir)%pqn(ish), il, &
            & param%record(ir)%slater(ish), cgto(ish, isp), .true., stat)
      end do

      do ish = 1, nsh_id(isp)
         if (ortho(ish) > 0) then
            call orthogonalize(cgto(ortho(ish), isp), cgto(ish, isp))
         end if
      end do
   end do

   call new_basis(calc%bas, mol, nsh_id, cgto, 1.0_wp)

end subroutine add_basis


subroutine add_ncoord(calc, mol, param)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param

   if (allocated(param%hamiltonian%cn)) then
      call new_ncoord(calc%ncoord, mol, cn_type=param%hamiltonian%cn)
   end if
end subroutine add_ncoord


subroutine add_hamiltonian(calc, mol, param, irc)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)

   call new_hamiltonian(calc%h0, mol, calc%bas, new_param_h0spec(mol, param, irc))
end subroutine add_hamiltonian


subroutine add_dispersion(calc, mol, param)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param

   type(d4_dispersion), allocatable :: d4
   type(d3_dispersion), allocatable :: d3

   if (.not.allocated(param%dispersion)) return
   associate(par => param%dispersion)
      if (par%d3) then
         allocate(d3)
         call new_d3_dispersion(d3, mol, s6=par%s6, s8=par%s8, a1=par%a1, a2=par%a2, s9=par%s9)
         call move_alloc(d3, calc%dispersion)
      else
         allocate(d4)
         call new_d4_dispersion(d4, mol, s6=par%s6, s8=par%s8, a1=par%a1, a2=par%a2, s9=par%s9)
         call move_alloc(d4, calc%dispersion)
      end if
   end associate
end subroutine add_dispersion

subroutine add_repulsion(calc, mol, param, irc)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)

   real(wp), allocatable :: alpha(:), zeff(:)
   real(wp), parameter :: rep_rexp = 1.0_wp

   if (.not.allocated(param%repulsion)) return
   allocate(calc%repulsion)
   alpha = param%record(irc)%alpha
   zeff = param%record(irc)%zeff
   call new_repulsion(calc%repulsion, mol, alpha, zeff, param%repulsion%kexp, &
      & param%repulsion%klight, rep_rexp)

end subroutine add_repulsion


subroutine add_halogen(calc, mol, param, irc)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)

   real(wp), allocatable :: bond_strength(:)

   if (.not.allocated(param%halogen)) return
   allocate(calc%halogen)
   bond_strength = param%record(irc)%xbond
   associate(par => param%halogen)
      call new_halogen_correction(calc%halogen, mol, par%damping, par%rscale, &
         & bond_strength)
   end associate
end subroutine add_halogen

subroutine add_coulomb(calc, mol, param, irc)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)

   real(wp), allocatable :: hardness(:, :), hubbard_derivs(:, :)
   real(wp), allocatable :: dkernel(:), qkernel(:), rad(:), vcn(:)
   procedure(average_interface), pointer :: average

   if (.not.any([allocated(param%charge), allocated(param%thirdorder), &
      & allocated(param%multipole)])) return
   allocate(calc%coulomb)

   if (allocated(param%charge)) then
      call get_shell_hardness(mol, param, irc, hardness)
      select case(param%charge%kernel)
      case(coulomb_kernel%effective)
         block
            type(effective_coulomb), allocatable :: es2
            allocate(es2)
            call get_average(param%charge%average, average)
            call new_effective_coulomb(es2, mol, param%charge%gexp, hardness, &
               & average, calc%bas%nsh_id)
            call move_alloc(es2, calc%coulomb%es2)
         end block
      case(coulomb_kernel%dftb_gamma)
         block
            type(gamma_coulomb), allocatable :: es2
            allocate(es2)
            call new_gamma_coulomb(es2, mol, hardness, calc%bas%nsh_id)
            call move_alloc(es2, calc%coulomb%es2)
         end block
      end select
   end if

   if (allocated(param%thirdorder)) then
      allocate(calc%coulomb%es3)
      if (param%thirdorder%shell) then
         call get_hubbard_derivs(mol, param, irc, hubbard_derivs)
         call new_onsite_thirdorder(calc%coulomb%es3, mol, hubbard_derivs, calc%bas%nsh_id)
      else
         allocate(hubbard_derivs(1, mol%nid))
         hubbard_derivs(1, :) = param%record(irc)%gam3
         call new_onsite_thirdorder(calc%coulomb%es3, mol, hubbard_derivs)
      end if
   end if

   if (allocated(param%multipole)) then
      allocate(calc%coulomb%aes2)
      dkernel = param%record(irc)%dkernel
      qkernel = param%record(irc)%qkernel
      rad = param%record(irc)%mprad
      vcn = param%record(irc)%mpvcn
      associate(par => param%multipole)
         call new_damped_multipole(calc%coulomb%aes2, mol, par%dmp3, par%dmp5, &
            & dkernel, qkernel, par%shift, par%kexp, par%rmax, rad, vcn)
      end associate
   end if

end subroutine add_coulomb


subroutine get_average(average_type, averager)
   character(len=*), intent(in) :: average_type
   procedure(average_interface), pointer, intent(out) :: averager
   select case(average_type)
   case default
      nullify(averager)
   case("harmonic")
      averager => harmonic_average
   case("geometric")
      averager => geometric_average
   case("arithmetic")
      averager => arithmetic_average
   end select
end subroutine get_average


subroutine get_shell_hardness(mol, param, irc, hardness)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)
   !> Shell resolved hardness parameters
   real(wp), allocatable, intent(out) :: hardness(:, :)

   integer :: isp, ir, ish, il

   allocate(hardness(maxval(param%record(irc)%nsh), mol%nid))
   hardness(:, :) = 0.0_wp
   do isp = 1, mol%nid
      ir = irc(isp)
      do ish = 1, param%record(ir)%nsh
         il = param%record(ir)%lsh(ish)
         hardness(ish, isp) = param%record(ir)%gam * param%record(ir)%lgam(ish)
      end do
   end do
end subroutine get_shell_hardness


subroutine get_hubbard_derivs(mol, param, irc, hubbard_derivs)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)
   !> Shell resolved Hubbard derivatives
   real(wp), allocatable, intent(out) :: hubbard_derivs(:, :)

   integer :: isp, ir, ish, il

   allocate(hubbard_derivs(maxval(param%record(irc)%nsh), mol%nid))
   hubbard_derivs(:, :) = 0.0_wp
   do isp = 1, mol%nid
      ir = irc(isp)
      do ish = 1, param%record(ir)%nsh
         il = param%record(ir)%lsh(ish)
         hubbard_derivs(ish, isp) = param%record(ir)%gam3 * param%thirdorder%ksh(il)
      end do
   end do
end subroutine get_hubbard_derivs


function new_param_h0spec(mol, param, irc) result(self)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), target, intent(in) :: param
   !> Record identifiers
   integer, target, intent(in) :: irc(:)
   !> Instance of the Hamiltonian specification
   type(param_h0spec) :: self

   integer :: isp, il, ir, ish
   integer :: ang_idx(0:4)

   self%param => param
   self%irc => irc

   allocate(self%valence(maxval(param%record(irc)%nsh), mol%nid))
   do isp = 1, mol%nid
      ang_idx = 0
      ir = irc(isp)
      do ish = 1, param%record(ir)%nsh
         il = param%record(ir)%lsh(ish)
         self%valence(ish, isp) = ang_idx(il) == 0
         if (self%valence(ish, isp)) ang_idx(il) = ish
      end do
   end do
end function new_param_h0spec


!> Generator for the enhancement factor to for scaling Hamiltonian elements
subroutine get_hscale(self, mol, bas, hscale)
   !> Instance of the Hamiltonian specification
   class(param_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Scaling parameters for the Hamiltonian elements
   real(wp), intent(out) :: hscale(:, :, :, :)

   integer :: isp, jsp, izp, jzp, ish, jsh, il, jl, ir, jr
   real(wp) :: zi, zj, zij, den, enp, km

   hscale(:, :, :, :) = 0.0_wp

   associate(par => self%param%hamiltonian, record => self%param%record, irc => self%irc)
      do isp = 1, mol%nid
         izp = mol%num(isp)
         ir = irc(isp)
         do jsp = 1, mol%nid
            jzp = mol%num(jsp)
            jr = irc(jsp)
            den = (record(ir)%en - record(jr)%en)**2
            do ish = 1, bas%nsh_id(isp)
               il = bas%cgto(ish, isp)%ang
               do jsh = 1, bas%nsh_id(jsp)
                  jl = bas%cgto(jsh, jsp)%ang
                  zi = record(ir)%slater(ish)
                  zj = record(jr)%slater(jsh)
                  zij = (2*sqrt(zi*zj)/(zi+zj))**par%wexp
                  if (self%valence(ish, isp) .and. self%valence(jsh, jsp)) then
                     enp = 1.0_wp + par%enscale * den
                     km = par%kpair(jr, ir) * par%ksh(jl, il) * enp
                  else if (self%valence(ish, isp) .and. .not.self%valence(jsh, jsp)) then
                     km = 0.5_wp * (par%ksh(il, il) + par%kpol)
                  else if (.not.self%valence(ish, isp) .and. self%valence(jsh, jsp)) then
                     km = 0.5_wp * (par%ksh(jl, jl) + par%kpol)
                  else
                     km = par%kpol
                  end if
                  hscale(jsh, ish, jsp, isp) = zij * km
               end do
            end do
         end do
      end do
   end associate
end subroutine get_hscale


!> Generator for the self energy / atomic levels of the Hamiltonian
subroutine get_selfenergy(self, mol, bas, selfenergy)
   !> Instance of the Hamiltonian specification
   class(param_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Self energy / atomic levels
   real(wp), intent(out) :: selfenergy(:, :)

   integer :: isp, ir, ish

   selfenergy(:, :) = 0.0_wp

   associate(record => self%param%record, irc => self%irc)
      do isp = 1, mol%nid
         ir = irc(isp)
         do ish = 1, bas%nsh_id(isp)
            selfenergy(ish, isp) = record(ir)%levels(ish)
         end do
      end do
   end associate
end subroutine get_selfenergy


!> Generator of the coordination number dependent shift of the self energy
subroutine get_cnshift(self, mol, bas, kcn)
   !> Instance of the Hamiltonian specification
   class(param_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Coordination number dependent shift
   real(wp), intent(out) :: kcn(:, :)

   integer :: isp, ir, ish

   kcn(:, :) = 0.0_wp

   associate(record => self%param%record, irc => self%irc)
      do isp = 1, mol%nid
         ir = irc(isp)
         do ish = 1, bas%nsh_id(isp)
            kcn(ish, isp) = record(ir)%kcn(ish)
         end do
      end do
   end associate
end subroutine get_cnshift


!> Generator for the polynomial parameters for the distant dependent scaling
subroutine get_shpoly(self, mol, bas, shpoly)
   !> Instance of the Hamiltonian specification
   class(param_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Polynomial parameters for distant dependent scaleing
   real(wp), intent(out) :: shpoly(:, :)

   integer :: isp, ir, ish

   shpoly(:, :) = 0.0_wp

   associate(record => self%param%record, irc => self%irc)
      do isp = 1, mol%nid
         ir = irc(isp)
         do ish = 1, bas%nsh_id(isp)
            shpoly(ish, isp) = record(ir)%shpoly(ish)
         end do
      end do
   end associate
end subroutine get_shpoly


subroutine get_reference_occ(self, mol, bas, refocc)
   !> Instance of the Hamiltonian specification
   class(param_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Reference occupation numbers
   real(wp), intent(out) :: refocc(:, :)

   integer :: isp, ir, ish

   refocc(:, :) = 0.0_wp

   associate(record => self%param%record, irc => self%irc)
      do isp = 1, mol%nid
         ir = irc(isp)
         do ish = 1, bas%nsh_id(isp)
            refocc(ish, isp) = merge(record(ir)%refocc(ish), 0.0_wp, self%valence(ish, isp))
         end do
      end do
   end associate
end subroutine get_reference_occ


subroutine update(self, mol)
   class(xtb_calculator), intent(inout) :: self
   type(structure_type), intent(in) :: mol
end subroutine update


!> Add an interaction container
subroutine push_back(self, cont)
   !> Instance of the tight-binding calculator
   class(xtb_calculator), intent(inout) :: self
   !> Container to be added
   class(container_type), allocatable, intent(inout) :: cont

   if (.not.allocated(self%interactions)) allocate(self%interactions)

   call self%interactions%push_back(cont)
end subroutine push_back


!> Add a container
subroutine pop(self, cont, idx)
   !> Instance of the tight-binding calculator
   class(xtb_calculator), intent(inout) :: self
   !> Container to be removed
   class(container_type), allocatable, intent(out) :: cont
   !> Index to remove container from
   integer, intent(in), optional :: idx

   if (.not.allocated(self%interactions)) return

   call self%interactions%pop(cont, idx)
end subroutine pop


pure function variable_info(self) result(info)
   use tblite_scf_info, only : scf_info, max
   !> Instance of the electrostatic container
   class(xtb_calculator), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info()

   if (allocated(self%coulomb)) then
      info = max(info, self%coulomb%variable_info())
   end if

   if (allocated(self%dispersion)) then
      info = max(info, self%dispersion%variable_info())
   end if

   if (allocated(self%interactions)) then
      info = max(info, self%interactions%variable_info())
   end if

end function variable_info


!> Information on container
pure function info(self, verbosity, indent) result(str)
   !> Instance of the interaction container
   class(xtb_calculator), intent(in) :: self
   !> Verbosity level
   integer, intent(in) :: verbosity
   !> Indentation level
   character(len=*), intent(in) :: indent
   !> Information on the container
   character(len=:), allocatable :: str

   character(len=*), parameter :: nl = new_line('a')

   str = "xTB calculator"

   if (allocated(self%repulsion)) then
      str = str // nl // indent // self%repulsion%info(verbosity, indent)
   end if

   if (allocated(self%coulomb)) then
      str = str // nl // indent // self%coulomb%info(verbosity, indent)
   end if

   if (allocated(self%dispersion)) then
      str = str // nl // indent // self%dispersion%info(verbosity, indent)
   end if

   if (allocated(self%interactions)) then
      str = str // nl // indent // self%interactions%info(verbosity, indent)
   end if
end function info

end module tblite_xtb_calculator
