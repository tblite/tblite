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

!> @file tblite/ncoord/ceh_en.f90
!> Provides an implementation for the electronegativity-weighted CN as used in the CEH method

!> Coordination number implementation with single error function and EN-weighting for the CEH method.
module tblite_ncoord_ceh_en
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use tblite_data_covrad_ceh, only : get_covalent_cehrad
   use tblite_data_paulingen_ceh, only : get_pauling_en_ceh
   use tblite_ncoord_type, only : ncoord_type
   implicit none
   private

   public :: new_ceh_en_ncoord

   !> Coordination number evaluator
   type, public, extends(ncoord_type) :: ceh_en_ncoord_type
      real(wp), allocatable :: rcov(:)
      real(wp), allocatable :: en(:)
   contains
      !> Evaluates the EN-weighted error counting function 
      procedure :: ncoord_count
      !> Evaluates the derivative of the EN-weighted error counting function 
      procedure :: ncoord_dcount
   end type ceh_en_ncoord_type

   !> Steepness of counting function
   real(wp), parameter :: kcn = 3.09_wp

   real(wp), parameter :: default_cutoff = 25.0_wp

  !  ! Pauling EN normalized to EN(F)=1
  !  ! TM and group 1/2 (from K on) hand optimized
  !  ! also adjusted: Rn,Xe,Kr,He,LNs
  !  real(wp), parameter :: en(86)    = (1.0_wp/3.98_wp) * [ &
  !  &         2.200,3.100&
  !  &        ,0.980,1.570&                          !Li
  !  &        ,2.040,2.550,3.040,3.440,3.980,4.500 & !  -Ne
  !  &        ,0.930,1.310&                          !Na-
  !  &        ,1.610,1.900,2.190,2.580,3.160,3.500 & !   Ar
  !  &        ,0.700,1.050 &                         !K-
  !  !             Sc    Ti     V     Cr   Mn     Fe    Co    Ni   Cu     Zn
  !  &        ,1.280,1.350,1.620,1.710,1.800,1.850,1.930,1.870,1.870,1.600 &
  !  &        ,1.810,2.010,2.180,2.550,2.960,3.200 & !   Kr
  !  &        ,0.700,0.900 &                         !Rb-
  !  !             Y     Zr     Nb    Mo   Tc     Ru    Rh    Pd   Ag     Cd
  !  &        ,1.320,1.380,1.570,1.800,1.900,2.180,2.300,2.100,1.800,1.600 &
  !  &        ,1.780,1.960,2.050,2.100,2.660,2.750 & !   Xe
  !  &        ,0.700,0.800 &                         !Cs-
  !  !             La
  !     ,1.050 &
  !     ,1.350,1.350,1.350,1.350,1.350,1.350,1.350,1.350,1.350,1.350,1.350,1.350,1.350,1.350  & ! Ce-Lu
  !  !                   Hf     Ta    W    Re     Os    Ir    Pt   Au     Hg
  !  &        ,      1.350,1.530,1.610,1.730,1.920,2.150,2.010,2.000,1.600 &
  !  &        ,1.620,2.330,2.020,2.000,2.200,2.600 ] !   Rn


contains


   subroutine new_ceh_en_ncoord(self, mol, cutoff, rcov, en)
      !> Coordination number container
      type(ceh_en_ncoord_type), intent(out) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Real space cutoff
      real(wp), intent(in), optional :: cutoff
      !> Covalent radii
      real(wp), intent(in), optional :: rcov(:)
      !> Covalent radii
      real(wp), intent(in), optional :: en(:)

      if (present(cutoff)) then
         self%cutoff = cutoff
      else
         self%cutoff = default_cutoff
      end if

      allocate(self%rcov(mol%nid))
      if (present(rcov)) then
         self%rcov(:) = rcov
      else
         self%rcov(:) = get_covalent_cehrad(mol%num)
      end if

      allocate(self%en(mol%nid))
      if (present(en)) then
         self%en(:) = en
      else
         self%en(:) = get_pauling_en_ceh(mol%num)
      end if

      !> CN is directed due to the EN contribution 
      !> i.e. contribution added to higher EN and removed from lower EN partner
      self%directed_factor = -1.0_wp

   end subroutine new_ceh_en_ncoord


!> Error counting function for coordination number contributions.
elemental function ncoord_count(self, mol, izp, jzp, r) result(count)

   !> Coordination number container
   class(ceh_en_ncoord_type), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Atom i index
   integer, intent(in)  :: izp

   !> Atom j index
   integer, intent(in)  :: jzp

   !> Current distance.
   real(wp), intent(in) :: r

   real(wp) :: rc, diff_en, count

   rc = self%rcov(izp) + self%rcov(jzp)
  
   diff_en = self%en(jzp) - self%en(izp)
   
   count = 0.5_wp * diff_en * (1.0_wp + erf(-kcn*(r-rc)/rc))

end function ncoord_count


!> Derivative of the error counting function w.r.t. the distance.
elemental function ncoord_dcount(self, mol, izp, jzp, r) result(count)

   !> Coordination number container
   class(ceh_en_ncoord_type), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Atom i index
   integer, intent(in)  :: izp

   !> Atom j index
   integer, intent(in)  :: jzp

   !> Current distance.
   real(wp), intent(in) :: r

   real(wp) :: rc, diff_en, exponent, expterm, count
   
   rc = self%rcov(izp) + self%rcov(jzp)
   
   diff_en = self%en(jzp) - self%en(izp)

   exponent = kcn*(r-rc)/rc
      
   expterm = exp(-exponent**2)
   
   count = - diff_en * (kcn*expterm)/(rc*sqrt(pi))

end function ncoord_dcount


end module tblite_ncoord_ceh_en
