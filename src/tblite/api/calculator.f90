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

!> API export for managing tight-binding parameters and calculators
module tblite_api_calculator
   use, intrinsic :: iso_c_binding
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_api_context, only : vp_context
   use tblite_api_result, only : vp_result
   use tblite_api_structure, only : vp_structure
   use tblite_api_version, only : namespace
   use tblite_wavefunction_mulliken, only : get_molecular_dipole_moment, &
      & get_molecular_quadrupole_moment
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_ipea1, only : new_ipea1_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: vp_calculator, delete_calculator_api
   public :: new_gfn2_calculator_api, new_ipea1_calculator_api, new_gfn1_calculator_api
   public :: set_calculator_mixer_damping_api, set_calculator_max_iter_api, &
      & set_calculator_accuracy_api, set_calculator_temperature_api
   public :: get_singlepoint_api


   !> Void pointer to calculator type
   type :: vp_calculator
      !> Actual payload
      type(xtb_calculator) :: ptr
      !> Calculation accuracy
      real(wp) :: accuracy = 1.0_wp
      !> Electronic temperature
      real(wp) :: etemp = 300.0_wp * 3.166808578545117e-06_wp
   end type vp_calculator


   logical, parameter :: debug = .false.


contains


function new_gfn2_calculator_api(vctx, vmol) result(vcalc) &
      & bind(C, name=namespace//"new_gfn2_calculator")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr) :: vcalc
   type(vp_calculator), pointer :: calc

   if (debug) print'("[Info]", 1x, a)', "new_gfn2_calculator"

   vcalc = c_null_ptr
   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vmol)) return
   call c_f_pointer(vmol, mol)

   allocate(calc)
   call new_gfn2_calculator(calc%ptr, mol%ptr)
   vcalc = c_loc(calc)

end function new_gfn2_calculator_api


function new_ipea1_calculator_api(vctx, vmol) result(vcalc) &
      & bind(C, name=namespace//"new_ipea1_calculator")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr) :: vcalc
   type(vp_calculator), pointer :: calc

   if (debug) print'("[Info]", 1x, a)', "new_ipea1_calculator"

   vcalc = c_null_ptr
   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vmol)) return
   call c_f_pointer(vmol, mol)

   allocate(calc)
   call new_ipea1_calculator(calc%ptr, mol%ptr)
   vcalc = c_loc(calc)

end function new_ipea1_calculator_api


function new_gfn1_calculator_api(vctx, vmol) result(vcalc) &
      & bind(C, name=namespace//"new_gfn1_calculator")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr) :: vcalc
   type(vp_calculator), pointer :: calc

   if (debug) print'("[Info]", 1x, a)', "new_gfn1_calculator"

   vcalc = c_null_ptr
   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vmol)) return
   call c_f_pointer(vmol, mol)

   allocate(calc)
   call new_gfn1_calculator(calc%ptr, mol%ptr)
   vcalc = c_loc(calc)

end function new_gfn1_calculator_api


subroutine delete_calculator_api(vcalc) &
      & bind(C, name=namespace//"delete_calculator")
   type(c_ptr), intent(inout) :: vcalc
   type(vp_calculator), pointer :: calc

   if (debug) print'("[Info]", 1x, a)', "delete_context"

   if (c_associated(vcalc)) then
      call c_f_pointer(vcalc, calc)

      deallocate(calc)
      vcalc = c_null_ptr
   end if
end subroutine delete_calculator_api


subroutine set_calculator_mixer_damping_api(vctx, vcalc, damping) &
      & bind(C, name=namespace//"set_calculator_mixer_damping")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   real(c_double), value :: damping
   type(error_type), allocatable :: error

   if (debug) print'("[Info]", 1x, a)', "set_calculator_mixer_damping"

   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vcalc)) then
      call fatal_error(error, "Calculator object is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vcalc, calc)

   calc%ptr%mixer_damping = damping
end subroutine set_calculator_mixer_damping_api


subroutine set_calculator_max_iter_api(vctx, vcalc, max_iter) &
      & bind(C, name=namespace//"set_calculator_max_iter")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   integer(c_int), value :: max_iter
   type(error_type), allocatable :: error

   if (debug) print'("[Info]", 1x, a)', "get_calculator_max_iter"

   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vcalc)) then
      call fatal_error(error, "Calculator object is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vcalc, calc)

   calc%ptr%max_iter = max_iter
end subroutine set_calculator_max_iter_api


subroutine set_calculator_accuracy_api(vctx, vcalc, accuracy) &
      & bind(C, name=namespace//"set_calculator_accuracy")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   real(c_double), value :: accuracy
   type(error_type), allocatable :: error

   if (debug) print'("[Info]", 1x, a)', "set_calculator_accuracy"

   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vcalc)) then
      call fatal_error(error, "Calculator object is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vcalc, calc)

   calc%accuracy = accuracy
end subroutine set_calculator_accuracy_api


subroutine set_calculator_temperature_api(vctx, vcalc, etemp) &
      & bind(C, name=namespace//"set_calculator_temperature")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   real(c_double), value :: etemp
   type(error_type), allocatable :: error

   if (debug) print'("[Info]", 1x, a)', "get_calculator_temperature"

   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vcalc)) then
      call fatal_error(error, "Calculator object is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vcalc, calc)

   calc%etemp = etemp
end subroutine set_calculator_temperature_api


subroutine get_singlepoint_api(vctx, vmol, vcalc, vres) &
      & bind(C, name=namespace//"get_singlepoint")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   type(error_type), allocatable :: error

   if (debug) print'("[Info]", 1x, a)', "get_singlepoint"

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

   if (.not.c_associated(vres)) then
      call fatal_error(error, "Result container is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vres, res)

   res%energy = 0.0_wp
   res%gradient = spread([0.0_wp, 0.0_wp, 0.0_wp], 2, mol%ptr%nat)
   res%sigma = spread([0.0_wp, 0.0_wp, 0.0_wp], 2, 3)
   res%dipole = spread(0.0_wp, 1, 3)
   res%quadrupole = spread(0.0_wp, 1, 6)
   call check_wavefunction(res%wfn, mol%ptr, calc%ptr, calc%etemp)

   call xtb_singlepoint(ctx%ptr, mol%ptr, calc%ptr, res%wfn, calc%accuracy, res%energy, &
      & res%gradient, res%sigma)

   call get_molecular_dipole_moment(mol%ptr, res%wfn%qat, res%wfn%dpat, res%dipole)
   call get_molecular_quadrupole_moment(mol%ptr, res%wfn%qat, res%wfn%dpat, res%wfn%qpat, &
      & res%quadrupole)

end subroutine get_singlepoint_api


subroutine check_wavefunction(wfn, mol, calc, etemp)
   type(wavefunction_type), allocatable, intent(inout) :: wfn
   type(structure_type), intent(in) :: mol
   type(xtb_calculator), intent(in) :: calc
   real(wp), intent(in) :: etemp

   if (allocated(wfn)) then
      wfn%kt = etemp

      if (size(wfn%qat) /= mol%nat .or. size(wfn%emo) /= calc%bas%nao &
         & .or. size(wfn%qsh) /= calc%bas%nsh) then
         deallocate(wfn)
      end if
   end if

   if (.not.allocated(wfn)) then
      allocate(wfn)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, etemp)
   end if
end subroutine check_wavefunction


end module tblite_api_calculator
