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

!> @file tblite/api/calculator.f90
!> Provides API exports for the #tblite_calculator handle.

!> API export for managing tight-binding parameters and calculators
module tblite_api_calculator
   use, intrinsic :: iso_c_binding
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_api_container, only : vp_container
   use tblite_api_context, only : vp_context
   use tblite_api_param, only : vp_param
   use tblite_api_result, only : vp_result
   use tblite_api_structure, only : vp_structure
   use tblite_api_version, only : namespace
   use tblite_results, only : results_type
   use tblite_wavefunction_mulliken, only : get_molecular_dipole_moment, &
      & get_molecular_quadrupole_moment
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction, sad_guess, eeq_guess
   use tblite_xtb_calculator, only : xtb_calculator, new_xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_ipea1, only : new_ipea1_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: vp_calculator, delete_calculator_api
   public :: new_gfn2_calculator_api, new_ipea1_calculator_api, new_gfn1_calculator_api, &
      & new_xtb_calculator_api
   public :: set_calculator_mixer_damping_api, set_calculator_max_iter_api, &
      & set_calculator_accuracy_api, set_calculator_temperature_api, &
      & set_calculator_save_integrals_api
   public :: get_singlepoint_api


   enum, bind(c)
      enumerator :: &
         tblite_guess_sad = 0_c_int, &
         tblite_guess_eeq = 1_c_int
   end enum


   !> Void pointer to calculator type
   type :: vp_calculator
      !> Actual payload
      type(xtb_calculator) :: ptr
      !> Calculation accuracy
      real(wp) :: accuracy = 1.0_wp
      !> Electronic temperature
      real(wp) :: etemp = 300.0_wp * 3.166808578545117e-06_wp
      !> Wavefunction guess
      integer :: guess = tblite_guess_sad
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

   if (debug) print '("[Info]", 1x, a)', "new_gfn2_calculator"

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

   if (debug) print '("[Info]", 1x, a)', "new_ipea1_calculator"

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

   if (debug) print '("[Info]", 1x, a)', "new_gfn1_calculator"

   vcalc = c_null_ptr
   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vmol)) return
   call c_f_pointer(vmol, mol)

   allocate(calc)
   call new_gfn1_calculator(calc%ptr, mol%ptr)
   vcalc = c_loc(calc)

end function new_gfn1_calculator_api


function new_xtb_calculator_api(vctx, vmol, vparam) result(vcalc) &
      & bind(C, name=namespace//"new_xtb_calculator")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr), value :: vparam
   type(vp_param), pointer :: param
   type(c_ptr) :: vcalc
   type(vp_calculator), pointer :: calc

   type(error_type), allocatable :: error

   if (debug) print '("[Info]", 1x, a)', "new_xtb_calculator"

   vcalc = c_null_ptr
   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vmol)) return
   call c_f_pointer(vmol, mol)

   if (.not.c_associated(vparam)) return
   call c_f_pointer(vparam, param)

   allocate(calc)
   call new_xtb_calculator(calc%ptr, mol%ptr, param%ptr, error)
   if (allocated(error)) then
      deallocate(calc)
   else
      vcalc = c_loc(calc)
   end if

end function new_xtb_calculator_api


subroutine delete_calculator_api(vcalc) &
      & bind(C, name=namespace//"delete_calculator")
   type(c_ptr), intent(inout) :: vcalc
   type(vp_calculator), pointer :: calc

   if (debug) print '("[Info]", 1x, a)', "delete_calculator"

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

   if (debug) print '("[Info]", 1x, a)', "set_calculator_mixer_damping"

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

   if (debug) print '("[Info]", 1x, a)', "get_calculator_max_iter"

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


subroutine set_calculator_guess_api(vctx, vcalc, guess) &
      & bind(C, name=namespace//"set_calculator_guess")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   integer(c_int), value :: guess
   type(error_type), allocatable :: error

   if (debug) print '("[Info]", 1x, a)', "set_calculator_guess"

   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vcalc)) then
      call fatal_error(error, "Calculator object is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vcalc, calc)

   calc%guess = guess
end subroutine set_calculator_guess_api


subroutine set_calculator_accuracy_api(vctx, vcalc, accuracy) &
      & bind(C, name=namespace//"set_calculator_accuracy")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   real(c_double), value :: accuracy
   type(error_type), allocatable :: error

   if (debug) print '("[Info]", 1x, a)', "set_calculator_accuracy"

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

   if (debug) print '("[Info]", 1x, a)', "get_calculator_temperature"

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


subroutine set_calculator_save_integrals_api(vctx, vcalc, save_integrals) &
      & bind(C, name=namespace//"set_calculator_save_integrals")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   integer(c_int), value :: save_integrals
   type(error_type), allocatable :: error

   if (debug) print '("[Info]", 1x, a)', "set_calculator_save_integrals"

   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vcalc)) then
      call fatal_error(error, "Calculator object is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vcalc, calc)

   calc%ptr%save_integrals = save_integrals /= 0
end subroutine set_calculator_save_integrals_api


subroutine get_calculator_shell_count(vctx, vcalc, nsh) &
      & bind(C, name=namespace//"get_calculator_shell_count")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   integer(c_int), intent(inout) :: nsh
   type(error_type), allocatable :: error

   if (debug) print '("[Info]", 1x, a)', "get_calculator_shell_count"

   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vcalc)) then
      call fatal_error(error, "Calculator object is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vcalc, calc)

   associate(sh2at => calc%ptr%bas%sh2at)
      nsh = size(sh2at)
   end associate
end subroutine get_calculator_shell_count

subroutine get_calculator_shell_map(vctx, vcalc, imap) &
      & bind(C, name=namespace//"get_calculator_shell_map")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   integer(c_int), intent(inout) :: imap(*)
   type(error_type), allocatable :: error

   if (debug) print '("[Info]", 1x, a)', "get_calculator_shell_map"

   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vcalc)) then
      call fatal_error(error, "Calculator object is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vcalc, calc)

   associate(sh2at => calc%ptr%bas%sh2at)
      imap(:size(sh2at)) = sh2at - 1_c_int
   end associate
end subroutine get_calculator_shell_map

subroutine get_calculator_angular_momenta(vctx, vcalc, am) &
      & bind(C, name=namespace//"get_calculator_angular_momenta")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   integer(c_int), intent(inout) :: am(*)
   type(error_type), allocatable :: error
   integer :: ish

   if (debug) print '("[Info]", 1x, a)', "get_calculator_angular_momenta"

   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vcalc)) then
      call fatal_error(error, "Calculator object is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vcalc, calc)

   associate(bas => calc%ptr%bas)
      do ish = 1, bas%nsh
         ! Angular momentum is only available for each species,
         ! which is unavailable in this API call.
         ! Therefore, we calculate it by using 2*l + 1 == nao
         am(ish) = (bas%nao_sh(ish)-1)/2
      end do
   end associate
end subroutine get_calculator_angular_momenta

subroutine get_calculator_orbital_count(vctx, vcalc, nao) &
      & bind(C, name=namespace//"get_calculator_orbital_count")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   integer(c_int), intent(inout) :: nao
   type(error_type), allocatable :: error

   if (debug) print '("[Info]", 1x, a)', "get_calculator_orbital_count"

   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vcalc)) then
      call fatal_error(error, "Calculator object is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vcalc, calc)

   associate(ao2sh => calc%ptr%bas%ao2sh)
      nao = size(ao2sh)
   end associate
end subroutine get_calculator_orbital_count

subroutine get_calculator_orbital_map(vctx, vcalc, imap) &
      & bind(C, name=namespace//"get_calculator_orbital_map")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   integer(c_int), intent(inout) :: imap(*)
   type(error_type), allocatable :: error

   if (debug) print '("[Info]", 1x, a)', "get_calculator_orbital_map"

   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vcalc)) then
      call fatal_error(error, "Calculator object is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vcalc, calc)

   associate(ao2sh => calc%ptr%bas%ao2sh)
      imap(:size(ao2sh)) = ao2sh - 1_c_int
   end associate
end subroutine get_calculator_orbital_map

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

   if (debug) print '("[Info]", 1x, a)', "get_singlepoint"

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
   res%results = results_type()
   call check_wavefunction(res%wfn, mol%ptr, calc%ptr, calc%etemp, calc%guess)

   call xtb_singlepoint(ctx%ptr, mol%ptr, calc%ptr, res%wfn, calc%accuracy, res%energy, &
      & res%gradient, res%sigma, results=res%results)

   call get_molecular_dipole_moment(mol%ptr, res%wfn%qat(:, 1), res%wfn%dpat(:, :, 1), &
      & res%dipole)
   call get_molecular_quadrupole_moment(mol%ptr, res%wfn%qat(:, 1), res%wfn%dpat(:, :, 1), &
      & res%wfn%qpat(:, :, 1), res%quadrupole)

end subroutine get_singlepoint_api


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
      call fatal_error(error, "Molecular structure data is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vcont, cont)

   call calc%ptr%push_back(cont%ptr)

   deallocate(cont)
   vcont = c_null_ptr
end subroutine push_back_api


subroutine check_wavefunction(wfn, mol, calc, etemp, guess)
   type(wavefunction_type), allocatable, intent(inout) :: wfn
   type(structure_type), intent(in) :: mol
   type(xtb_calculator), intent(in) :: calc
   real(wp), intent(in) :: etemp
   integer, intent(in) :: guess

   integer, parameter :: nspin = 1

   if (allocated(wfn)) then
      wfn%kt = etemp

      if (size(wfn%qat) /= mol%nat .or. size(wfn%emo) /= calc%bas%nao &
         & .or. size(wfn%qsh) /= calc%bas%nsh) then
         deallocate(wfn)
      end if
   end if

   if (.not.allocated(wfn)) then
      allocate(wfn)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, nspin, etemp)

      select case(guess)
      case default
         call sad_guess(mol, calc, wfn)
      case(tblite_guess_eeq)
         call eeq_guess(mol, calc, wfn)
      end select
   end if
end subroutine check_wavefunction


end module tblite_api_calculator
