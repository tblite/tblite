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

!> @file tblite/api/container.f90
!> Provides API exports for the #tblite_container handle.

!> API export for managing interaction containers
module tblite_api_container
   use, intrinsic :: iso_c_binding
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type
   use tblite_api_version, only : namespace
   use tblite_api_calculator, only : vp_calculator
   use tblite_api_context, only : vp_context
   use tblite_api_structure, only : vp_structure
   use tblite_basis, only : basis_type
   use tblite_container, only : container_type, container_list
   use tblite_data_spin, only : get_spin_constant
   use tblite_external_field, only : electric_field
   use tblite_spin, only : spin_polarization, new_spin_polarization
   use tblite_solvation, only : solvation_input, cpcm_input, alpb_input, &
      & solvent_data, get_solvent_data, solvation_type, new_solvation, solutionState, &
      & new_solvation_cds, new_solvation_shift, cds_input, shift_input
   use tblite_api_utils, only: c_f_character
   implicit none
   private

   public :: vp_container, delete_container_api


   !> Void pointer to a container instance
   type :: vp_container
      !> Actual container
      class(container_type), allocatable :: ptr
   end type vp_container

   logical, parameter :: debug = .false.

contains


subroutine push_back_api(vctx, vcalc, vcont) &
      & bind(C, name=namespace//"calculator_push_back")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   type(c_ptr), intent(inout) :: vcont
   type(vp_container), pointer :: cont

   type(error_type), allocatable :: error

   if (debug) print '("[Info]", 1x, a)', "calculator_push_back"

   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vcalc)) then
      call fatal_error(error, "Calculator object is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vcalc, calc)

   if (.not.c_associated(vcont)) then
      call fatal_error(error, "Container object is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vcont, cont)

   ! Hacky way to propagate spin channel information to calculator
   select type(tcont => cont%ptr)
   type is(spin_polarization)
      calc%nspin = 2
   end select

   call calc%ptr%push_back(cont%ptr)

   deallocate(cont)
   vcont = c_null_ptr
end subroutine push_back_api


function new_electric_field_api(efield) result(vcont) &
      & bind(C, name=namespace//"new_electric_field")
   real(c_double), intent(in) :: efield(3)
   type(c_ptr) :: vcont
   type(vp_container), pointer :: cont

   allocate(cont)
   cont%ptr = electric_field(efield)
   vcont = c_loc(cont)
end function new_electric_field_api


function new_spin_polarization_api(vctx, vmol, vcalc, wscale) result(vcont) &
      & bind(C, name=namespace//"new_spin_polarization")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   type(c_ptr) :: vcont
   type(vp_container), pointer :: cont
   real(c_double), value :: wscale

   type(error_type), allocatable :: error
   type(spin_polarization), allocatable :: spin
   real(wp), allocatable :: wll(:, :, :)

   if (debug) print '("[Info]", 1x, a)', "new_spin_polarization"
   vcont = c_null_ptr

   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vmol)) then
      call fatal_error(error, "Molecular structure data is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vmol, mol)

   if (.not.c_associated(vcalc)) then
      call fatal_error(error, "Calculator object is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vcalc, calc)

   allocate(spin)
   call get_spin_constants(wll, mol%ptr, calc%ptr%bas)
   wll(:, :, :) = wscale * wll
   call new_spin_polarization(spin, mol%ptr, wll, calc%ptr%bas%nsh_id)

   allocate(cont)
   call move_alloc(spin, cont%ptr)
   vcont = c_loc(cont)
end function new_spin_polarization_api

function new_cpcm_solvation_solvent_api(vctx, vmol, vcalc, solvstr) result(vcont) &
   & bind(C, name=namespace//"new_cpcm_solvation_solvent")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   type(c_ptr) :: vcont
   type(vp_container), pointer :: cont
   character(kind=c_char), intent(in) :: solvstr(*)
   character(len=:), allocatable :: solvinp
   type(solvation_input) :: solvmodel
   type(solvent_data) :: solvent
   class(solvation_type), allocatable :: solv
   type(error_type), allocatable :: error
   integer :: stat
   logical :: ok

   if (debug) print '("[Info]", 1x, a)', "new_cpcm_solvation"
   vcont = c_null_ptr

   call resolve_ptr_input(vctx, vmol, vcalc, ctx, mol, calc, ok)
   if (.not.ok) return

   call c_f_character(solvstr, solvinp)

   solvent = get_solvent_data(solvinp)
   if (solvent%eps <= 0.0_wp) then
      call fatal_error(error, "String value for epsilon was not found among database of solvents")
      call ctx%ptr%set_error(error)
      return
   end if
   solvmodel%cpcm = cpcm_input(solvent%eps)
   call new_solvation(solv, mol%ptr, solvmodel, error)
   if (allocated(error)) then 
      call ctx%ptr%set_error(error)
      return
   end if
   
   allocate(cont)
   call move_alloc(solv, cont%ptr)
   
   vcont = c_loc(cont)

end function

function new_cpcm_solvation_epsilon_api(vctx, vmol, vcalc, eps) result(vcont) &
   & bind(C, name=namespace//"new_cpcm_solvation_epsilon")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   type(c_ptr) :: vcont
   type(vp_container), pointer :: cont
   real(kind=c_double), value :: eps
   type(solvation_input) :: solvmodel
   class(solvation_type), allocatable :: solv
   type(error_type), allocatable :: error
   integer :: stat
   logical :: ok
   if (debug) print '("[Info]", 1x, a)', "new_cpcm_solvation"
   vcont = c_null_ptr

   call resolve_ptr_input(vctx, vmol, vcalc, ctx, mol, calc, ok)
   if (.not.ok) return

   solvmodel%cpcm = cpcm_input(eps)
   call new_solvation(solv, mol%ptr, solvmodel, error)
   if (allocated(error)) then
      call ctx%ptr%set_error(error)
      return
   end if
   
   allocate(cont)
   call move_alloc(solv, cont%ptr)
   
   vcont = c_loc(cont)

end function

function new_alpb_solvation_solvent_api(vctx, vmol, vcalc, solvstr, c_refstate) result(vcont) &
   & bind(C, name=namespace//"new_alpb_solvation_solvent")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   type(c_ptr) :: vcont
   type(vp_container), pointer :: cont
   character(kind=c_char), intent(in) :: solvstr(*)
   integer(c_int), optional, intent(in) :: c_refstate(1)
   character(len=:), allocatable :: solvinp, refstr
   type(container_list), allocatable :: cont_list
   type(error_type), allocatable :: error
   integer :: stat, kernel = 2, refstate
   logical :: ok, alpb = .true.
   
   if (debug) print '("[Info]", 1x, a)', "new_alpb_solvation"
   vcont = c_null_ptr

   call resolve_ptr_input(vctx, vmol, vcalc, ctx, mol, calc, ok)
   if (.not.ok) return

   call c_f_character(solvstr, solvinp)
   if (present(c_refstate)) then 
      refstate = c_refstate(1)
   else
      refstate = 1
   end if
   call get_ref_state_from_enum(refstate, refstr, error)

   call setup_gbsa_alpb_solvent_model(solvinp, refstr, calc%ptr%method, kernel, alpb, mol%ptr, cont_list, error)
   if (allocated(error)) then
      call ctx%ptr%set_error(error)
      return
   end if

   allocate(cont)
   call move_alloc(cont_list, cont%ptr)

   vcont = c_loc(cont)
   
end function


function new_alpb_solvation_epsilon_api(vctx, vmol, vcalc, eps, c_refstate) result(vcont) &
   & bind(C, name=namespace//"new_alpb_solvation_epsilon")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   type(c_ptr) :: vcont
   type(vp_container), pointer :: cont
   real(c_double), value :: eps
   integer(c_int), optional, intent(in) :: c_refstate(1)
   type(solvation_input) :: solvmodel
   class(solvation_type), allocatable :: solv
   type(error_type), allocatable :: error
   character(len=:), allocatable :: refstr
   integer :: stat, refstate
   logical :: ok, alpb = .true.

   if (debug) print '("[Info]", 1x, a)', "new_alpb_solvation float input"
   vcont = c_null_ptr

   call resolve_ptr_input(vctx, vmol, vcalc, ctx, mol, calc, ok)
   if (.not.ok) return
   if (present(c_refstate)) then 
      refstate = c_refstate(1)
   else
      refstate = 1
   end if
   call get_ref_state_from_enum(refstate, refstr, error)
   
   if (refstr /= "gsolv") then
      call fatal_error(error, "Solution state shift is only supported for named solvents")
      call ctx%ptr%set_error(error)
      return
   end if

   solvmodel%alpb = alpb_input(eps, alpb=alpb)
   call new_solvation(solv, mol%ptr, solvmodel, error)
   if (allocated(error)) then 
      call ctx%ptr%set_error(error)
      return
   end if
   
   allocate(cont)
   call move_alloc(solv, cont%ptr)

   vcont = c_loc(cont)
   
end function


function new_gbsa_solvation_solvent_api(vctx, vmol, vcalc, solvstr, c_refstate) result(vcont) &
   & bind(C, name=namespace//"new_gbsa_solvation_solvent")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   type(c_ptr) :: vcont
   type(vp_container), pointer :: cont
   character(kind=c_char), intent(in) :: solvstr(*)
   integer(c_int), optional, intent(in) :: c_refstate(1)
   character(len=:), allocatable :: solvinp, refstr
   type(container_list), allocatable :: cont_list
   type(error_type), allocatable :: error
   integer :: stat, kernel = 1, refstate
   logical :: ok, alpb = .false.
   
   if (debug) print '("[Info]", 1x, a)', "new_gbsa_solvation"
   vcont = c_null_ptr

   call resolve_ptr_input(vctx, vmol, vcalc, ctx, mol, calc, ok)
   if (.not.ok) return

   call c_f_character(solvstr, solvinp)
   if (present(c_refstate)) then 
      refstate = c_refstate(1)
   else
      refstate = 1
   end if
   call get_ref_state_from_enum(refstate, refstr, error)
   

   call setup_gbsa_alpb_solvent_model(solvinp, refstr, calc%ptr%method, kernel, alpb, mol%ptr, cont_list, error)
   if (allocated(error)) then
      call ctx%ptr%set_error(error)
      return
   end if

   allocate(cont)
   call move_alloc(cont_list, cont%ptr)

   vcont = c_loc(cont)
   
end function


function new_gbsa_solvation_epsilon_api(vctx, vmol, vcalc, eps, c_refstate) result(vcont) &
   & bind(C, name=namespace//"new_gbsa_solvation_epsilon")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   type(c_ptr) :: vcont
   type(vp_container), pointer :: cont
   real(c_double), value :: eps
   integer(c_int), optional, intent(in) :: c_refstate(1)
   type(solvation_input) :: solvmodel
   class(solvation_type), allocatable :: solv
   type(error_type), allocatable :: error
   character(len=:), allocatable :: refstr
   integer :: stat, refstate
   logical :: ok, alpb = .false.

   if (debug) print '("[Info]", 1x, a)', "new_gbsa_solvation float input"
   vcont = c_null_ptr

   call resolve_ptr_input(vctx, vmol, vcalc, ctx, mol, calc, ok)
   if (.not.ok) return
   if (present(c_refstate)) then 
      refstate = c_refstate(1)
   else
      refstate = 1
   end if
   call get_ref_state_from_enum(refstate, refstr, error)
   
   if (refstr /= "gsolv") then
      call fatal_error(error, "Solution state shift is only supported for named solvents")
      call ctx%ptr%set_error(error)
      return
   end if

   solvmodel%alpb = alpb_input(eps, alpb=alpb)
   call new_solvation(solv, mol%ptr, solvmodel, error)
   if (allocated(error)) then
      call ctx%ptr%set_error(error)
      return
   end if

   allocate(cont)
   call move_alloc(solv, cont%ptr)

   vcont = c_loc(cont)
   
end function


subroutine resolve_ptr_input(vctx, vmol, vcalc, ctx, mol, calc, ok)
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   type(error_type), allocatable :: error
   logical :: ok
   ok = .false.
   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vmol)) then
      call fatal_error(error, "Molecular structure data is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vmol, mol)

   if (.not.c_associated(vcalc)) then
      call fatal_error(error, "Calculator object is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vcalc, calc)
   ok = .true.
end subroutine

subroutine get_spin_constants(wll, mol, bas)
   real(wp), allocatable, intent(out) :: wll(:, :, :)
   type(structure_type), intent(in) :: mol
   type(basis_type), intent(in) :: bas

   integer :: izp, ish, jsh, il, jl

   allocate(wll(bas%nsh, bas%nsh, mol%nid), source=0.0_wp)

   do izp = 1, mol%nid
      do ish = 1, bas%nsh_id(izp)
         il = bas%cgto(ish, izp)%ang
         do jsh = 1, bas%nsh_id(izp)
            jl = bas%cgto(jsh, izp)%ang
            wll(jsh, ish, izp) = get_spin_constant(jl, il, mol%num(izp))
         end do
      end do
   end do
end subroutine get_spin_constants


subroutine delete_container_api(vcont) &
      & bind(C, name=namespace//"delete_container")
   type(c_ptr), intent(inout) :: vcont
   type(vp_container), pointer :: cont

   if (debug) print '("[Info]", 1x, a)', "delete_container"

   if (c_associated(vcont)) then
      call c_f_pointer(vcont, cont)

      deallocate(cont)
      vcont = c_null_ptr
   end if
end subroutine delete_container_api


subroutine setup_gbsa_alpb_solvent_model(solvstr, refstr, method, kernel, alpb, mol, cont_list, error)
   !> String containing the solvent name
   character(len=:), allocatable, intent(in) :: solvstr
   !> String for the reference state
   character(len=:), allocatable, intent(in) :: refstr
   !> String with the short name of the GFN Hamiltoniann in use
   character(len=:), allocatable, intent(inout) :: method
   !> which kernel are we using
   integer, intent(in) :: kernel
   !> do we use gbsa or alpb
   logical, intent(in) :: alpb
   type(structure_type), intent(in) :: mol
   !> Container list passed back, containing all solvent contributions
   type(container_list), allocatable, intent(out) :: cont_list
   !> Error object returned for bad input
   type(error_type), allocatable, intent(out) :: error
   integer :: sol_state
   
   type(solvation_input) :: solvmodel
   type(solvent_data) :: solvent
   
   class(container_type), allocatable :: tmp_cont
   class(solvation_type), allocatable :: solv
   class(solvation_type), allocatable :: cds
   class(solvation_type), allocatable :: shift
   allocate(cont_list)
   solvent = get_solvent_data(solvstr)
   if (solvent%eps <= 0.0_wp) then
      call fatal_error(error, "String value for epsilon was not found among database of solvents")
      return
   end if

   sol_state = solutionState%gsolv
   select case(refstr)
   case("gsolv")
      sol_state = solutionState%gsolv
   case("bar1mol")
      sol_state = solutionState%bar1mol
   case("reference")
      sol_state = solutionState%reference
   case default
      sol_state = solutionState%gsolv
   end select

   solvmodel%alpb = alpb_input(solvent%eps, solvent=solvent%solvent, kernel=kernel, alpb=alpb)
   solvmodel%cds = cds_input(alpb=alpb, solvent=solvent%solvent)
   solvmodel%shift = shift_input(alpb=alpb, solvent=solvent%solvent, state=sol_state)
   if (method == "custom") method = "gfn2"
   call new_solvation(solv, mol, solvmodel, error, method)
   call move_alloc(solv, tmp_cont)
   if (allocated(error)) return
   call cont_list%push_back(tmp_cont)
   
   call new_solvation_cds(cds, mol, solvmodel, error, method)
   call move_alloc(cds, tmp_cont)
   if (allocated(error)) return
   call cont_list%push_back(tmp_cont)
   call new_solvation_shift(shift, solvmodel, error, method)
   if (allocated(error)) return
   call move_alloc(shift, tmp_cont)
   call cont_list%push_back(tmp_cont)

end subroutine

subroutine get_ref_state_from_enum(refstate, refstr, error)
   integer(c_int), intent(in) :: refstate
   character(len=:), allocatable, intent(out) :: refstr
   type(error_type), allocatable, intent(inout) :: error

   select case(refstate)
   case(1)
      refstr = "gsolv"
   case(2)
      refstr = "bar1mol"
   case(3)
      refstr = "reference"
   case default
      call fatal_error(error, "Given reference state is not known check enumerator in container.h for valid options.")
   end select 

end subroutine

end module tblite_api_container
