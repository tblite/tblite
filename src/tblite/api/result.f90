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

!> @file tblite/api/result.f90
!> Provides API exports for the #tblite_result handle.

!> API export for managing calculation results
module tblite_api_result
   use, intrinsic :: iso_c_binding
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_api_error, only : vp_error
   use tblite_api_version, only : namespace
   use tblite_results, only : results_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: vp_result, new_result_api, copy_result_api, delete_result_api
   public :: get_result_number_of_atoms_api, get_result_number_of_shells_api, &
      & get_result_number_of_orbitals_api, get_result_energy_api, get_result_gradient_api, &
      & get_result_virial_api, get_result_charges_api, get_result_dipole_api, &
      & get_result_quadrupole_api, get_result_orbital_energies_api, &
      & get_result_orbital_occupations_api, get_result_orbital_coefficients_api, &
      & get_result_energies_api, get_result_density_matrix_api, &
      & get_result_overlap_matrix_api, get_result_hamiltonian_matrix_api, &
      & get_result_bond_orders_api


   !> Void pointer holding results of a calculation
   type :: vp_result
      !> Single point energy
      real(wp), allocatable :: energy
      !> Molecular gradient
      real(wp), allocatable :: gradient(:, :)
      !> Virial
      real(wp), allocatable :: sigma(:, :)
      !> Dipole moment
      real(wp), allocatable :: dipole(:)
      !> Quadrupole moment
      real(wp), allocatable :: quadrupole(:)
      !> Wavefunction
      type(wavefunction_type), allocatable :: wfn
      !> Additional results
      type(results_type), allocatable :: results
   end type vp_result


   logical, parameter :: debug = .false.


contains


!> Create new result container
function new_result_api() result(vres) &
      & bind(C, name=namespace//"new_result")
   type(vp_result), pointer :: res
   type(c_ptr) :: vres

   if (debug) print '("[Info]", 1x, a)', "new_result"

   allocate(res)
   vres = c_loc(res)

end function new_result_api


!> Create copy result container
function copy_result_api(vold) result(vnew) &
      & bind(C, name=namespace//"copy_result")
   type(vp_result), pointer :: new
   type(c_ptr), value :: vold
   type(vp_result), pointer :: old
   type(c_ptr) :: vnew

   if (debug) print '("[Info]", 1x, a)', "copy_result"

   vnew = c_null_ptr
   if (c_associated(vold)) then
      call c_f_pointer(vold, old)
      allocate(new, source=old)
      vnew = c_loc(new)
   end if

end function copy_result_api


!> Delete result container
subroutine delete_result_api(vres) &
      & bind(C, name=namespace//"delete_result")
   type(c_ptr), intent(inout) :: vres
   type(vp_result), pointer :: res

   if (debug) print '("[Info]", 1x, a)', "delete_result"

   if (c_associated(vres)) then
      call c_f_pointer(vres, res)

      deallocate(res)
      vres = c_null_ptr
   end if
end subroutine delete_result_api


subroutine get_result_number_of_atoms_api(verror, vres, natoms) &
      & bind(C, name=namespace//"get_result_number_of_atoms")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   integer(c_int), intent(out) :: natoms
   logical :: ok

   if (debug) print '("[Info]", 1x, a)', "get_result_number_of_atoms"

   call get_result(verror, vres, error, res, ok)
   if (.not.ok) return

   if (.not.allocated(res%wfn)) then
      call fatal_error(error%ptr, "Result does not contain number of atoms")
      return
   end if

   natoms = size(res%wfn%qat)
end subroutine get_result_number_of_atoms_api


subroutine get_result_number_of_shells_api(verror, vres, nshells) &
      & bind(C, name=namespace//"get_result_number_of_shells")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   integer(c_int), intent(out) :: nshells
   logical :: ok

   if (debug) print '("[Info]", 1x, a)', "get_result_number_of_shells"

   call get_result(verror, vres, error, res, ok)
   if (.not.ok) return

   if (.not.allocated(res%wfn)) then
      call fatal_error(error%ptr, "Result does not contain number of shells")
      return
   end if

   nshells = size(res%wfn%qsh)
end subroutine get_result_number_of_shells_api


subroutine get_result_number_of_orbitals_api(verror, vres, norb) &
      & bind(C, name=namespace//"get_result_number_of_orbitals")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   integer(c_int), intent(out) :: norb
   logical :: ok

   if (debug) print '("[Info]", 1x, a)', "get_result_number_of_orbitals"

   call get_result(verror, vres, error, res, ok)
   if (.not.ok) return

   if (.not.allocated(res%wfn)) then
      call fatal_error(error%ptr, "Result does not contain number of orbitals")
      return
   end if

   norb = size(res%wfn%emo)
end subroutine get_result_number_of_orbitals_api


subroutine get_result_energy_api(verror, vres, energy) &
      & bind(C, name=namespace//"get_result_energy")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   real(c_double), intent(out) :: energy
   logical :: ok

   if (debug) print '("[Info]", 1x, a)', "get_result_energy"

   call get_result(verror, vres, error, res, ok)
   if (.not.ok) return

   if (.not.allocated(res%energy)) then
      call fatal_error(error%ptr, "Result does not contain energy")
      return
   end if

   energy = res%energy
end subroutine get_result_energy_api


subroutine get_result_energies_api(verror, vres, energies) &
      & bind(C, name=namespace//"get_result_energies")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   real(c_double), intent(out) :: energies(*)
   logical :: ok

   if (debug) print '("[Info]", 1x, a)', "get_result_energies"

   call get_result(verror, vres, error, res, ok)
   if (.not.ok) return

   if (.not.allocated(res%results)) then
      call fatal_error(error%ptr, "Result does not contain energies")
      return
   end if

   if (.not.allocated(res%results%energies)) then
      call fatal_error(error%ptr, "Result does not contain energies")
      return
   end if

   energies(:size(res%results%energies)) = res%results%energies
end subroutine get_result_energies_api


subroutine get_result_gradient_api(verror, vres, gradient) &
      & bind(C, name=namespace//"get_result_gradient")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   real(c_double), intent(out) :: gradient(*)
   logical :: ok

   if (debug) print '("[Info]", 1x, a)', "get_result_gradient"

   call get_result(verror, vres, error, res, ok)
   if (.not.ok) return

   if (.not.allocated(res%gradient)) then
      call fatal_error(error%ptr, "Result does not contain gradient")
      return
   end if

   gradient(:size(res%gradient)) = reshape(res%gradient, [size(res%gradient)])
end subroutine get_result_gradient_api


subroutine get_result_virial_api(verror, vres, sigma) &
      & bind(C, name=namespace//"get_result_virial")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   real(c_double), intent(out) :: sigma(*)
   logical :: ok

   if (debug) print '("[Info]", 1x, a)', "get_result_virial"

   call get_result(verror, vres, error, res, ok)
   if (.not.ok) return

   if (.not.allocated(res%sigma)) then
      call fatal_error(error%ptr, "Result does not contain virial")
      return
   end if

   sigma(:size(res%sigma)) = reshape(res%sigma, [size(res%sigma)])
end subroutine get_result_virial_api


subroutine get_result_charges_api(verror, vres, charges) &
      & bind(C, name=namespace//"get_result_charges")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   real(c_double), intent(out) :: charges(*)
   logical :: ok

   if (debug) print '("[Info]", 1x, a)', "get_result_charges"

   call get_result(verror, vres, error, res, ok)
   if (.not.ok) return

   if (.not.allocated(res%wfn)) then
      call fatal_error(error%ptr, "Result does not contain atomic charges")
      return
   end if

   charges(:size(res%wfn%qat)) = res%wfn%qat(:, 1)
end subroutine get_result_charges_api


subroutine get_result_dipole_api(verror, vres, dipole) &
      & bind(C, name=namespace//"get_result_dipole")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   real(c_double), intent(out) :: dipole(*)
   logical :: ok

   if (debug) print '("[Info]", 1x, a)', "get_result_dipole"

   call get_result(verror, vres, error, res, ok)
   if (.not.ok) return

   if (.not.allocated(res%dipole)) then
      call fatal_error(error%ptr, "Result does not contain dipole moment")
      return
   end if

   dipole(:size(res%dipole)) = res%dipole
end subroutine get_result_dipole_api


subroutine get_result_quadrupole_api(verror, vres, quadrupole) &
      & bind(C, name=namespace//"get_result_quadrupole")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   real(c_double), intent(out) :: quadrupole(*)
   logical :: ok

   if (debug) print '("[Info]", 1x, a)', "get_result_quadrupole"

   call get_result(verror, vres, error, res, ok)
   if (.not.ok) return

   if (.not.allocated(res%quadrupole)) then
      call fatal_error(error%ptr, "Result does not contain quadrupole moment")
      return
   end if

   quadrupole(:size(res%quadrupole)) = res%quadrupole
end subroutine get_result_quadrupole_api


subroutine get_result_orbital_energies_api(verror, vres, emo) &
      & bind(C, name=namespace//"get_result_orbital_energies")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   real(c_double), intent(out) :: emo(*)
   logical :: ok

   if (debug) print '("[Info]", 1x, a)', "get_result_orbital_energies"

   call get_result(verror, vres, error, res, ok)
   if (.not.ok) return

   if (.not.allocated(res%wfn)) then
      call fatal_error(error%ptr, "Result does not contain orbital energies")
      return
   end if

   emo(:size(res%wfn%emo)) = res%wfn%emo(:, 1)
end subroutine get_result_orbital_energies_api


subroutine get_result_orbital_occupations_api(verror, vres, occ) &
      & bind(C, name=namespace//"get_result_orbital_occupations")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   real(c_double), intent(out) :: occ(*)
   logical :: ok

   if (debug) print '("[Info]", 1x, a)', "get_result_orbital_occupations"

   call get_result(verror, vres, error, res, ok)
   if (.not.ok) return

   if (.not.allocated(res%wfn)) then
      call fatal_error(error%ptr, "Result does not contain orbital occupations")
      return
   end if

   occ(:size(res%wfn%focc)) = res%wfn%focc(:, 1)
end subroutine get_result_orbital_occupations_api


subroutine get_result_orbital_coefficients_api(verror, vres, cmo) &
      & bind(C, name=namespace//"get_result_orbital_coefficients")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   real(c_double), intent(out) :: cmo(*)
   logical :: ok

   if (debug) print '("[Info]", 1x, a)', "get_result_orbital_coefficients"

   call get_result(verror, vres, error, res, ok)
   if (.not.ok) return

   if (.not.allocated(res%wfn)) then
      call fatal_error(error%ptr, "Result does not contain orbital coefficients")
      return
   end if

   cmo(:size(res%wfn%coeff)) = reshape(res%wfn%coeff(:, :, 1), [size(res%wfn%coeff)])
end subroutine get_result_orbital_coefficients_api


subroutine get_result_density_matrix_api(verror, vres, pmat) &
      & bind(C, name=namespace//"get_result_density_matrix")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   real(c_double), intent(out) :: pmat(*)
   logical :: ok

   if (debug) print '("[Info]", 1x, a)', "get_result_density_matrix"

   call get_result(verror, vres, error, res, ok)
   if (.not.ok) return

   if (.not.allocated(res%wfn)) then
      call fatal_error(error%ptr, "Result does not contain density matrix")
      return
   end if

   pmat(:size(res%wfn%density)) = reshape(res%wfn%density(:, :, 1), [size(res%wfn%density)])
end subroutine get_result_density_matrix_api


subroutine get_result_overlap_matrix_api(verror, vres, smat) &
      & bind(C, name=namespace//"get_result_overlap_matrix")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   real(c_double), intent(out) :: smat(*)
   logical :: ok

   if (debug) print '("[Info]", 1x, a)', "get_result_overlap_matrix"

   call get_result(verror, vres, error, res, ok)
   if (.not.ok) return

   if (.not.allocated(res%results)) then
      call fatal_error(error%ptr, "Result does not contain overlap matrix")
      return
   end if

   if (.not.allocated(res%results%overlap)) then
      call fatal_error(error%ptr, "Result does not contain overlap matrix")
      return
   end if

   smat(:size(res%results%overlap)) = &
      & reshape(res%results%overlap, [size(res%results%overlap)])
end subroutine get_result_overlap_matrix_api


subroutine get_result_hamiltonian_matrix_api(verror, vres, hmat) &
      & bind(C, name=namespace//"get_result_hamiltonian_matrix")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   real(c_double), intent(out) :: hmat(*)
   logical :: ok

   if (debug) print '("[Info]", 1x, a)', "get_result_hamiltonian_matrix"

   call get_result(verror, vres, error, res, ok)
   if (.not.ok) return

   if (.not.allocated(res%results)) then
      call fatal_error(error%ptr, "Result does not contain Hamiltonian matrix")
      return
   end if

   if (.not.allocated(res%results%hamiltonian)) then
      call fatal_error(error%ptr, "Result does not contain Hamiltonian matrix")
      return
   end if

   hmat(:size(res%results%hamiltonian)) = &
      & reshape(res%results%hamiltonian, [size(res%results%hamiltonian)])
end subroutine get_result_hamiltonian_matrix_api


subroutine get_result_bond_orders_api(verror, vres, mbo) &
      & bind(C, name=namespace//"get_result_bond_orders")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vres
   type(vp_result), pointer :: res
   real(c_double), intent(out) :: mbo(*)
   logical :: ok

   if (debug) print '("[Info]", 1x, a)', "get_result_bond_orders"

   call get_result(verror, vres, error, res, ok)
   if (.not.ok) return

   if (.not.allocated(res%results)) then
      call fatal_error(error%ptr, "Result does not contain bond orders")
      return
   end if

   if (.not.allocated(res%results%bond_orders)) then
      call fatal_error(error%ptr, "Result does not contain bond orders")
      return
   end if

   mbo(:size(res%results%bond_orders)) = &
      & reshape(res%results%bond_orders, [size(res%results%bond_orders)])
end subroutine get_result_bond_orders_api


subroutine get_result(verror, vres, error, res, ok)
   type(c_ptr), intent(in) :: verror
   type(c_ptr), intent(in) :: vres
   type(vp_error), pointer, intent(out) :: error
   type(vp_result), pointer, intent(out) :: res
   logical, intent(out) :: ok

   ok = .false.
   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vres)) then
      call fatal_error(error%ptr, "Result container is missing")
      return
   end if
   call c_f_pointer(vres, res)
   ok = .true.
end subroutine get_result


end module tblite_api_result
