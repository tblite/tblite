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

!> @dir tblite/coulomb
!> Contains the Coulomb related interactions.

!> @file tblite/coulomb.f90
!> Reexports of the Coulombic interaction containers.

!> Proxy module for handling Coulomb-type interactions
module tblite_coulomb
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_charge, only : coulomb_charge_type, coulomb_kernel, &
      & effective_coulomb, new_effective_coulomb, &
      & average_interface, harmonic_average, arithmetic_average, geometric_average, &
      & gamma_coulomb, new_gamma_coulomb
   use tblite_coulomb_multipole, only : damped_multipole, new_damped_multipole
   use tblite_coulomb_thirdorder, only : onsite_thirdorder, new_onsite_thirdorder
   use tblite_coulomb_type, only : coulomb_type
   implicit none

   public :: coulomb_type
   public :: coulomb_cache
   public :: coulomb_charge_type, coulomb_kernel
   public :: effective_coulomb, new_effective_coulomb
   public :: average_interface, harmonic_average, arithmetic_average, geometric_average
   public :: gamma_coulomb, new_gamma_coulomb
   public :: damped_multipole, new_damped_multipole
   public :: onsite_thirdorder, new_onsite_thirdorder
end module tblite_coulomb
