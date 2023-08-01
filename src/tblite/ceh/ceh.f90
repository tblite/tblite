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

!> @dir tblite/ceh/ceh.f90
!> Contains the implementation of the Charge Extended Hückel (CEH) method.


module tblite_ceh_ceh
   use mctc_env, only : error_type, wp
   use mctc_io, only: structure_type
   use tblite_basis_ortho, only : orthogonalize
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_basis_type, only : cgto_type, new_basis
   use tblite_ceh_coordination, only : ceh_ncoord
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: new_ceh_calculator

   integer, parameter, private :: max_elem = 86
   integer, parameter, private :: max_shell = 3

   !> Number of shells # MM, August 01, 2023
   integer, parameter :: nshell(max_elem) = [&
   & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 2, 3, & ! 1-20
   & 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, & ! 21-40
   & 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, & ! 41-60
   & 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, & ! 61-80
   & 3, 3, 3, 3, 3, 3]                                             ! 81-86

   !> Angular momentum of each shell # MM, August 01, 2023
   ! 0 = s, 1 = p, 2 = d
   ! CAUTION: Ordering from original CEH model is taken for consistency with the parameterization
   ! I.e., the ordering of the shells is always: "s", "p", "d"
   integer, parameter :: ang_shell(max_shell, max_elem) = reshape([&
   & 0, 0, 0,  0, 0, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0, & ! 1-7
   & 0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 2,  0, 1, 2, & ! 8-14
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 0,  0, 1, 2,  0, 1, 2, & ! 15-21
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 22-28
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 29-35
   & 0, 1, 2,  0, 1, 0,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 36-42
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 43-49
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 0,  0, 1, 2, & ! 50-56
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 57-63
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 64-70
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 71-77
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 78-84
   & 0, 1, 2,  0, 1, 2], shape(ang_shell))                                  ! 85-86

   !> Principal quantum number of each shell
   integer, parameter :: principal_quantum_number(max_shell, max_elem) = reshape([&
   & 1, 2, 0,  1, 0, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0, & ! 1-7
   & 2, 2, 0,  2, 2, 0,  2, 2, 3,  3, 3, 0,  3, 3, 0,  3, 3, 3,  3, 3, 3, & ! 8-14
   & 3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  4, 4, 0,  4, 4, 3,  4, 4, 3, & ! 15-21
   & 4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3, & ! 22-28
   & 4, 4, 3,  4, 4, 3,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4, & ! 29-35
   & 4, 4, 4,  5, 5, 0,  5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 4, & ! 36-42
   & 5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 5, & ! 43-49
   & 5, 5, 5,  5, 5, 5,  5, 5, 5,  5, 5, 5,  5, 5, 5,  6, 6, 0,  6, 6, 5, & ! 50-56
   & 6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5, & ! 57-63
   & 6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5, & ! 64-70
   & 6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5, & ! 71-77
   & 6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5, & ! 78-84
   & 6, 6, 5,  6, 6, 5], shape(principal_quantum_number))                   ! 85-86

   !> Number of primitive gaussians per shell # MM, August 01, 2023
   integer, parameter :: number_of_primitives(max_shell, max_elem) = reshape([&
   & 4, 0, 0,  4, 0, 0,  4, 4, 0,  4, 4, 0,  4, 4, 0,  4, 4, 0,  4, 4, 0, & ! 1-7
   & 4, 4, 0,  4, 4, 0,  4, 4, 0,  4, 4, 0,  4, 4, 0,  4, 4, 4,  4, 4, 4, & ! 8-14
   & 4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 0,  4, 4, 4,  4, 4, 4, & ! 15-21
   & 4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4, & ! 22-28
   & 4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4, & ! 29-35
   & 4, 4, 4,  4, 4, 0,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4, & ! 36-42
   & 4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4, & ! 43-49
   & 4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  6, 6, 0,  6, 6, 4, & ! 50-56
   & 6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4, & ! 57-63
   & 6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4, & ! 64-70
   & 6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4, & ! 71-77
   & 6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4, & ! 78-84
   & 6, 6, 4,  6, 6, 4], shape(number_of_primitives))                       ! 85-86

   !> Exponent of the Slater function # MM, August 01, 2023
   real(wp), parameter :: slater_exponent(max_shell, max_elem) = reshape([&
   &  1.13056314_wp,  0.00000000_wp,  0.00000000_wp,  1.50065094_wp,  0.00000000_wp,  0.00000000_wp, &
   &  1.19413208_wp,  1.00435986_wp,  0.00000000_wp,  1.04160365_wp,  1.58776712_wp,  0.00000000_wp, &
   &  1.80615826_wp,  1.68397018_wp,  0.00000000_wp,  1.92578557_wp,  1.64086778_wp,  0.00000000_wp, &
   &  2.16223922_wp,  1.79669019_wp,  0.00000000_wp,  2.31036294_wp,  2.07689830_wp,  0.00000000_wp, &
   &  2.63689506_wp,  2.39472851_wp,  0.00000000_wp,  3.02968365_wp,  2.22048798_wp,  0.00000000_wp, &
   &  1.75555631_wp,  1.33838074_wp,  0.00000000_wp,  1.70435065_wp,  1.39174584_wp,  0.00000000_wp, &
   &  1.80712664_wp,  1.62636355_wp,  2.28258151_wp,  2.02930097_wp,  1.63990927_wp,  1.85910992_wp, &
   &  2.24859618_wp,  1.77685356_wp,  1.76755789_wp,  2.72470838_wp,  2.01797698_wp,  2.04253078_wp, &
   &  2.79738729_wp,  2.29520297_wp,  1.61932628_wp,  2.53272473_wp,  1.81071325_wp,  1.69070808_wp, &
   &  1.98719413_wp,  1.22379806_wp,  0.00000000_wp,  1.62021149_wp,  1.54648435_wp,  2.66580857_wp, &
   &  2.14684072_wp,  1.36068696_wp,  1.74878015_wp,  0.69639491_wp,  1.03767289_wp,  1.86366428_wp, &
   &  1.69851374_wp,  2.17779533_wp,  1.83191727_wp,  1.12031062_wp,  1.41500509_wp,  2.16118253_wp, &
   &  1.76160033_wp,  1.98316214_wp,  2.06245601_wp,  1.92205103_wp,  1.30492118_wp,  2.50946021_wp, &
   &  2.07573975_wp,  1.46306532_wp,  2.63099323_wp,  2.44652113_wp,  1.68710423_wp,  3.04820511_wp, &
   &  1.80368223_wp,  1.71760742_wp,  3.02263698_wp,  1.92342980_wp,  1.82453676_wp,  3.49908671_wp, &
   &  2.24151618_wp,  2.08921540_wp,  2.12567962_wp,  2.23205133_wp,  1.89430927_wp,  1.97513779_wp, &
   &  2.69088900_wp,  2.12267272_wp,  1.94145540_wp,  2.94651193_wp,  2.35079302_wp,  1.79682240_wp, &
   &  2.62665945_wp,  2.50566210_wp,  1.75613478_wp,  2.72057897_wp,  2.16290669_wp,  1.88632737_wp, &
   &  1.70116961_wp,  1.64545292_wp,  0.00000000_wp,  1.58807973_wp,  1.80661603_wp,  2.99671810_wp, &
   &  1.08395599_wp,  1.42679219_wp,  2.19662852_wp,  1.47634306_wp,  1.29489380_wp,  2.37913368_wp, &
   &  2.69077953_wp,  1.95585989_wp,  2.11756166_wp,  2.44347621_wp,  1.73442204_wp,  2.02065708_wp, &
   &  2.23908639_wp,  1.95016590_wp,  2.18967157_wp,  2.40381905_wp,  1.44282665_wp,  2.59962263_wp, &
   &  1.96295898_wp,  1.58471984_wp,  2.95895238_wp,  2.38179585_wp,  2.23204400_wp,  3.22819136_wp, &
   &  2.51271048_wp,  1.85519318_wp,  3.24389973_wp,  2.11457529_wp,  1.96275823_wp,  3.51477490_wp, &
   &  2.38640874_wp,  2.12255171_wp,  2.27363232_wp,  2.62089386_wp,  2.06375224_wp,  2.75185414_wp, &
   &  2.79382121_wp,  2.28410073_wp,  2.04070104_wp,  2.97205755_wp,  2.32287719_wp,  1.74734495_wp, &
   &  2.72692370_wp,  2.65571946_wp,  1.85299047_wp,  2.73000797_wp,  2.51650283_wp,  2.24007230_wp, &
   &  1.68663815_wp,  2.16354335_wp,  0.00000000_wp,  1.71010093_wp,  1.58444852_wp,  3.02035513_wp, &
   &  1.11292580_wp,  1.67649145_wp,  2.19394682_wp,  2.04931361_wp,  2.02737301_wp,  2.57681508_wp, &
   &  2.01436702_wp,  2.02975950_wp,  2.59448588_wp,  1.97942043_wp,  2.03214600_wp,  2.61215668_wp, &
   &  1.94447383_wp,  2.03453249_wp,  2.62982748_wp,  1.90952724_wp,  2.03691898_wp,  2.64749828_wp, &
   &  1.87458065_wp,  2.03930547_wp,  2.66516908_wp,  1.83963406_wp,  2.04169196_wp,  2.68283988_wp, &
   &  1.80468747_wp,  2.04407845_wp,  2.70051068_wp,  1.76974088_wp,  2.04646494_wp,  2.71818148_wp, &
   &  1.73479429_wp,  2.04885143_wp,  2.73585228_wp,  1.69984770_wp,  2.05123792_wp,  2.75352308_wp, &
   &  1.66490111_wp,  2.05362441_wp,  2.77119388_wp,  1.62995451_wp,  2.05601091_wp,  2.78886467_wp, &
   &  1.59500792_wp,  2.05839740_wp,  2.80653547_wp,  1.23998866_wp,  1.55902021_wp,  2.56002100_wp, &
   &  2.21476417_wp,  1.21067246_wp,  2.27168445_wp,  2.35183881_wp,  2.17165636_wp,  2.38705146_wp, &
   &  2.29642436_wp,  2.38947110_wp,  2.49180058_wp,  2.78850456_wp,  1.82182727_wp,  3.11161203_wp, &
   &  2.18906814_wp,  1.84759137_wp,  3.45359549_wp,  2.08875168_wp,  2.56532973_wp,  3.59828518_wp, &
   &  2.39586959_wp,  2.60766526_wp,  3.87307779_wp,  2.39636916_wp,  2.53950059_wp,  3.95389569_wp, &
   &  2.72134635_wp,  2.34645674_wp,  2.67736087_wp,  2.99323675_wp,  2.17128412_wp,  5.38232746_wp, &
   &  3.05653080_wp,  2.49075153_wp,  2.67638930_wp,  3.02905756_wp,  2.63479560_wp,  1.94913437_wp, &
   &  2.54746694_wp,  2.83550170_wp,  1.88029428_wp,  2.26386287_wp,  2.46706218_wp,  2.09966650_wp],&
   & shape(slater_exponent))

contains

   subroutine new_ceh_calculator(mol, error)

      type(structure_type), intent(in) :: mol
      type(error_type), intent(out)   :: error
      type(xtb_calculator)            :: calc
      !  local variables
      real(wp),allocatable :: F(:), eps(:)    ! Fock and eigenvalues
      real(wp),allocatable :: S(:), P(:)      ! overlap and density
      real(wp),allocatable :: D(:,:)          ! dipole integrals
      real(wp),allocatable :: cn1(:),cn2(:)   ! CN for H0
      real(wp),allocatable :: norm(:)         ! AO norm
      real(wp),allocatable :: psh(:,:)        ! shell populations
      real(wp),allocatable :: wbo(:,:)        ! BO

      real(wp)             :: tmp, hav, hdii, hdjj, ps, dum, dcal1, dcal2
      real(wp)             :: klli, kllj
      real(wp)             :: dip(3)
      integer              :: i, j, k, l, m, ii, jj, iat, jat
      integer              :: ati, atj, ish, jsh, iish, li, lj, ij
      integer              :: ierr
      integer, parameter   :: llao(0:3) = [1,3,5,7] ! # of AO in shell
      logical              :: EX, FIT
      character(len=80)    :: atmp

      ! PARAMETER SECTION
      real(wp), parameter   :: kll(1:3) = [0.6366_wp, 0.9584_wp, 1.2320_wp] ! H0 spd-scaling 1.23
      real(wp), parameter   :: kcn  = -3.09d0 ! CN erf expo

      !allocate(F(ndim*(ndim+1)/2), S(ndim*(ndim+1)/2), P(ndim*(ndim+1)/2), wbo(n,n), &
      !&         cn1(n), cn2(n), D(ndim*(ndim+1)/2,3), norm(ndim), eps(ndim), psh(maxsh,n), &
      !&         source=0.0_wp, stat=ierr)

      write(*,*) "CEH: Initializing"
      write(*,*) "Hello world!"
      call add_ceh_basis(calc, mol)

      ! call ceh_ncoord(n,at,kcn,rab,cn1,cn2) ! routine with standard radii-> CN, EN weigthed CN

      ! call ncoord_basq(n,at,    rab,cn_eeq)  ! CN for q-vSZP basis compression term

      ! call sint_ceh(norm,S,F)                ! sig,pi,del scaled overlap on F, standard overlap on S
      ! call dipint_ceh(norm,D) ! dipole ints


   end subroutine new_ceh_calculator

   subroutine add_ceh_basis(calc, mol)
      !> Instance of the xTB evaluator
      type(xtb_calculator), intent(inout) :: calc
      !> Molecular structure data
      type(structure_type), intent(in) :: mol

      integer :: isp, izp, ish, stat, ng, il
      integer, allocatable :: nsh_id(:)
      integer :: ang_idx(0:2), ortho(max_shell)
      type(cgto_type), allocatable :: cgto(:, :)

      nsh_id = nshell(mol%num)
      allocate(cgto(maxval(nsh_id), mol%nid))
      do isp = 1, mol%nid
         ang_idx = 0
         ortho = 0
         izp = mol%num(isp)
         do ish = 1, nsh_id(isp)
            il = ang_shell(ish, izp)
            ng = number_of_primitives(ish, izp)
            if (ang_idx(il) > 0) then
               ortho(ish) = ang_idx(il)
            else
               ang_idx(il) = ish
            end if
            call slater_to_gauss(ng, principal_quantum_number(ish, izp), il, &
            & slater_exponent(ish, izp), cgto(ish, isp), .true., stat)
         end do

         do ish = 1, nsh_id(isp)
            if (ortho(ish) > 0) then
               call orthogonalize(cgto(ortho(ish), isp), cgto(ish, isp))
            end if
         end do
      end do

      call new_basis(calc%bas, mol, nsh_id, cgto, 1.0_wp)

   end subroutine add_ceh_basis

end module tblite_ceh_ceh
