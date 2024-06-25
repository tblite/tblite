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

!> @dir tblite/ncoord
!> Contains the implementation for the coordination number evaluators.

!> @file tblite/ncoord.f90
!> Reexports the coordination number evaluation modules.

!> Proxy module to expose coordination number containers
module tblite_ncoord
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_ncoord_type, only : ncoord_type
   use tblite_ncoord_gfn, only : gfn_ncoord_type, new_gfn_ncoord
   use tblite_ncoord_exp, only : exp_ncoord_type, new_exp_ncoord
   use tblite_ncoord_erf, only : erf_ncoord_type, new_erf_ncoord
   use tblite_ncoord_erf_en, only : erf_en_ncoord_type, new_erf_en_ncoord
   implicit none
   private

   public :: ncoord_type, new_ncoord

contains

!> Create a new generic coordination number container
subroutine new_ncoord(self, mol, cn_type, kcn, rcov, en)
   !> Instance of the coordination number container
   class(ncoord_type), allocatable, intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Coordination number type
   character(len=*), intent(in) :: cn_type
   !> Steepness of counting function
   real(wp), optional :: kcn
   !> Optional set of covalent radii to be used in CN
   real(wp), intent(in), optional :: rcov(:)
   !> Optional set of electronegativity to be use din CN
   real(wp), intent(in), optional :: en(:)


   select case(cn_type)
   case("exp")
      block
         type(exp_ncoord_type), allocatable :: tmp
         allocate(tmp)
         call new_exp_ncoord(tmp, mol, kcn=kcn, rcov=rcov)
         call move_alloc(tmp, self)
      end block
   case("gfn")
      block
         type(gfn_ncoord_type), allocatable :: tmp
         allocate(tmp)
         call new_gfn_ncoord(tmp, mol, rcov=rcov)
         call move_alloc(tmp, self)
      end block
   case("erf")
      block
         type(erf_ncoord_type), allocatable :: tmp
         allocate(tmp)
         call new_erf_ncoord(tmp, mol, kcn=kcn, rcov=rcov)
         call move_alloc(tmp, self)
      end block
   case("erf_en")
      block
         type(erf_en_ncoord_type), allocatable :: tmp
         allocate(tmp)
         call new_erf_en_ncoord(tmp, mol, kcn=kcn, rcov=rcov, en=en)
         call move_alloc(tmp, self)
      end block
   end select
end subroutine new_ncoord

end module tblite_ncoord
