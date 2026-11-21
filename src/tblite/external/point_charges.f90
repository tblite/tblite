! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! External point-charge embedding for the xTB interaction model.
module tblite_external_point_charges
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_container_cache, only : container_cache
   use tblite_container_type, only : container_type
   use tblite_partition, only : owns_index
   use tblite_scf_info, only : scf_info, shell_resolved
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private
   public :: point_charges, new_point_charges

   type, public, extends(container_type) :: point_charges
      real(wp), allocatable :: charge(:), xyz(:, :), hardness(:), shell_hardness(:)
      integer, allocatable :: shell_atom(:)
      real(wp) :: gexp = 2.0_wp
      logical :: harmonic = .false.
   contains
      procedure :: get_potential
      procedure :: variable_info
      procedure :: get_energy
      procedure :: get_gradient
      procedure :: point_gradient
   end type point_charges

contains

pure function variable_info(self) result(info)
   class(point_charges), intent(in) :: self
   type(scf_info) :: info

   info = scf_info(charge=shell_resolved)
end function variable_info

subroutine new_point_charges(self, charge, xyz, hardness, shell_hardness, shell_atom, gexp, harmonic)
   type(point_charges), intent(out) :: self
   real(wp), intent(in) :: charge(:), xyz(:, :), hardness(:), shell_hardness(:), gexp
   integer, intent(in) :: shell_atom(:)
   logical, intent(in) :: harmonic
   self%charge = charge
   self%xyz = xyz
   self%hardness = hardness
   self%shell_hardness = shell_hardness
   self%shell_atom = shell_atom
   self%gexp = gexp
   self%harmonic = harmonic
   self%label = "external point charges"
end subroutine new_point_charges

pure real(wp) function pair_hardness(self, ish, ipc) result(gam)
   class(point_charges), intent(in) :: self
   integer, intent(in) :: ish, ipc
   if (self%harmonic) then
      gam = 2.0_wp/(1.0_wp/self%shell_hardness(ish) + 1.0_wp/self%hardness(ipc))
   else
      gam = 0.5_wp*(self%shell_hardness(ish) + self%hardness(ipc))
   end if
end function pair_hardness

subroutine get_potential(self, mol, cache, wfn, pot)
   class(point_charges), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: cache
   type(wavefunction_type), intent(in) :: wfn
   type(potential_type), intent(inout) :: pot
   integer :: ish, ipc, iat
   real(wp) :: dr(3), r, gam, value
   if (.not.owns_index(self%partition, 1)) return
   do ish = 1, size(self%shell_hardness)
      iat = self%shell_atom(ish)
      do ipc = 1, size(self%charge)
         dr = mol%xyz(:, iat) - self%xyz(:, ipc)
         r = norm2(dr)
         gam = pair_hardness(self, ish, ipc)
         value = self%charge(ipc)/(r**self%gexp + gam**(-self%gexp))**(1.0_wp/self%gexp)
         pot%vsh(ish, 1) = pot%vsh(ish, 1) + value
      end do
   end do
end subroutine get_potential

subroutine get_energy(self, mol, cache, wfn, energies)
   class(point_charges), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: cache
   type(wavefunction_type), intent(in) :: wfn
   real(wp), intent(inout) :: energies(:)
   integer :: iat, ish, ipc
   real(wp) :: dr(3), r, gam, value
   if (.not.owns_index(self%partition, 1)) return
   do iat = 1, mol%nat
      do ish = 1, size(self%shell_hardness)
         if (self%shell_atom(ish) /= iat) cycle
         do ipc = 1, size(self%charge)
            dr = mol%xyz(:, iat) - self%xyz(:, ipc)
            r = norm2(dr)
            gam = pair_hardness(self, ish, ipc)
            value = self%charge(ipc)/(r**self%gexp + gam**(-self%gexp))**(1.0_wp/self%gexp)
            energies(iat) = energies(iat) + wfn%qsh(ish, 1)*value
         end do
      end do
   end do
end subroutine get_energy

subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   class(point_charges), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: cache
   type(wavefunction_type), intent(in) :: wfn
   real(wp), contiguous, intent(inout) :: gradient(:, :), sigma(:, :)
   integer :: iat, ish, ipc
   real(wp) :: dr(3), r, gam, fac
   if (.not.owns_index(self%partition, 1)) return
   do iat = 1, mol%nat
      do ish = 1, size(self%shell_hardness)
         if (self%shell_atom(ish) /= iat) cycle
         do ipc = 1, size(self%charge)
            dr = mol%xyz(:, iat) - self%xyz(:, ipc)
            r = norm2(dr)
            gam = pair_hardness(self, ish, ipc)
            fac = -wfn%qsh(ish, 1)*self%charge(ipc)*r**(self%gexp-2.0_wp) &
               & /(r**self%gexp + gam**(-self%gexp))**(1.0_wp/self%gexp+1.0_wp)
            gradient(:, iat) = gradient(:, iat) + fac*dr
         end do
      end do
   end do
end subroutine get_gradient

subroutine point_gradient(self, mol, wfn, gradient)
   class(point_charges), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wavefunction_type), intent(in) :: wfn
   real(wp), intent(out) :: gradient(:, :)
   integer :: iat, ish, ipc
   real(wp) :: dr(3), r, gam, fac
   gradient = 0.0_wp
   if (.not.owns_index(self%partition, 1)) return
   do iat = 1, mol%nat
      do ish = 1, size(self%shell_hardness)
         if (self%shell_atom(ish) /= iat) cycle
         do ipc = 1, size(self%charge)
            dr = mol%xyz(:, iat) - self%xyz(:, ipc)
            r = norm2(dr)
            gam = pair_hardness(self, ish, ipc)
            fac = wfn%qsh(ish, 1)*self%charge(ipc)*r**(self%gexp-2.0_wp) &
               & /(r**self%gexp + gam**(-self%gexp))**(1.0_wp/self%gexp+1.0_wp)
            gradient(:, ipc) = gradient(:, ipc) + fac*dr
         end do
      end do
   end do
end subroutine point_gradient

end module tblite_external_point_charges
