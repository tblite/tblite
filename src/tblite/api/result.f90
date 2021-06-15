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

module tblite_api_result
   use, intrinsic :: iso_c_binding
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_api_error, only : vp_error
   use tblite_api_version, only : namespace
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: vp_result


   !> Void pointer holding results of a calculation
   type :: vp_result
      !> Single point energy
      real(wp), allocatable :: energy

      !> Molecular gradient
      real(wp), allocatable :: gradient(:, :)

      !> Virial
      real(wp), allocatable :: sigma(:, :)

      !> Wavefunction
      type(wavefunction_type), allocatable :: wfn
   end type vp_result


   logical, parameter :: debug = .false.


contains


!> Create new result container
function new_result_api() result(vres) &
      & bind(C, name=namespace//"new_result")
   type(vp_result), pointer :: res
   type(c_ptr) :: vres

   if (debug) print'("[Info]", 1x, a)', "new_result"

   allocate(res)
   vres = c_loc(res)

end function new_result_api


!> Delete result container
subroutine delete_result_api(vres) &
      & bind(C, name=namespace//"delete_result")
   type(c_ptr), intent(inout) :: vres
   type(vp_result), pointer :: res

   if (debug) print'("[Info]", 1x, a)', "delete_result"

   if (c_associated(vres)) then
      call c_f_pointer(vres, res)

      deallocate(res)
      vres = c_null_ptr
   end if
end subroutine delete_result_api


subroutine get_result_energy(verror, vres, energy) &
      & bind(C, name=namespace//"get_result_energy")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   real(c_double), intent(out) :: energy

   if (debug) print'("[Info]", 1x, a)', "get_result_energy"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vres)) then
      call fatal_error(error%ptr, "Result container is missing")
      return
   end if
   call c_f_pointer(vres, res)

   if (.not.allocated(res%energy)) then
      call fatal_error(error%ptr, "Result does not contain energy")
      return
   end if

   energy = res%energy
end subroutine get_result_energy


subroutine get_result_gradient(verror, vres, gradient) &
      & bind(C, name=namespace//"get_result_gradient")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   real(c_double), intent(out) :: gradient(*)

   if (debug) print'("[Info]", 1x, a)', "get_result_gradient"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vres)) then
      call fatal_error(error%ptr, "Result container is missing")
      return
   end if
   call c_f_pointer(vres, res)

   if (.not.allocated(res%gradient)) then
      call fatal_error(error%ptr, "Result does not contain gradient")
      return
   end if

   gradient(:size(res%gradient)) = reshape(res%gradient, [size(res%gradient)])
end subroutine get_result_gradient


subroutine get_result_virial(verror, vres, sigma) &
      & bind(C, name=namespace//"get_result_virial")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   real(c_double), intent(out) :: sigma(*)

   if (debug) print'("[Info]", 1x, a)', "get_result_virial"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vres)) then
      call fatal_error(error%ptr, "Result container is missing")
      return
   end if
   call c_f_pointer(vres, res)

   if (.not.allocated(res%sigma)) then
      call fatal_error(error%ptr, "Result does not contain virial")
      return
   end if

   sigma(:size(res%sigma)) = reshape(res%sigma, [size(res%sigma)])
end subroutine get_result_virial


end module tblite_api_result
