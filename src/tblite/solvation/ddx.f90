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

!> @file tblite/solvation/ddx.f90
!> Provides a polarizable continuum model

!> Implicit solvation model based on a polarizable dielectric continuum
module tblite_solvation_ddx
   use ddx, only: allocate_state, check_error, ddinit, ddrun, ddx_error_type, &
      & ddx_state_type, ddx_type, fill_guess, fill_guess_adjoint, setup, &
      & solvation_force_terms, solve, solve_adjoint
   use ddx_core, only: ddx_electrostatics_type
   use ddx_multipolar_solutes, only: multipole_electrostatics, multipole_force_terms, multipole_psi
   use mctc_env, only: error_type, fatal_error, wp
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use omp_lib, only: omp_get_max_threads
   use tblite_blas, only: dot, gemv
   use tblite_container_cache, only: container_cache
   use tblite_mesh_lebedev, only: grid_size
   use tblite_scf_info, only: atom_resolved, scf_info
   use tblite_scf_potential, only: potential_type
   use tblite_solvation_data, only: get_vdw_rad_cosmo
   use tblite_solvation_type, only: solvation_type
   use tblite_wavefunction_type, only: wavefunction_type

   implicit none
   private

   public :: ddx_solvation, new_ddx, ddx_input, ddx_cache, ddx_solvation_model

   !> Possible solvation models to be used within the dd framework
   type :: enum_ddx_solvation_model
      ! COSMO, CPCM, and PCM are internally defined as 100, 101, and 200 to get correct feps in the first step
      ! and to avoid collisions with other aliase elsewhere in the code. 
      ! Labels are dumbed to ddX-specific values 1 and 2 before passing the model input to the ddX routine.
      !> Conductor like screening model
      integer :: cosmo = 100
      !> Conductor-like polarizable continuum model
      integer :: cpcm = 101
      !> Polarizable continuum model
      integer :: pcm = 200
   end type enum_ddx_solvation_model

   !> Actual enumerator for the dd solvation models
   type(enum_ddx_solvation_model), parameter :: ddx_solvation_model = enum_ddx_solvation_model()


   !> Input for ddX solvation
   type :: ddx_input
      !> Dielectric constant
      real(wp) :: dielectric_const
      !> ddx model
      integer :: ddx_model = ddx_solvation_model%cosmo
      !> Scaling of van-der-Waals radii
      real(wp) :: rscale = 1.0_wp
      !> Accuracy for iterative solver
      real(wp) :: conv = 1.0e-10_wp
      !> Regularization parameter
      real(wp) :: eta = 0.1_wp
      !> Number of grid points for each atom (=110)
      integer :: nang = grid_size(8)
      !> Maximum angular momentum of basis functions
      integer :: lmax = 1
      !> Van-der-Waals radii for all species
      real(wp), allocatable :: rvdw(:)
      !> Number of OMP threads
      integer :: nproc = 1
      !> Shift of the characteristic function 
      ! (default value depends on the model, for COSMO/CPCM it is -1)
      real(wp) :: shift = -1.0_wp
      !> Maximum number of iterations for the iterative solver
      integer :: max_iter = 100
      !> Number of extrapolation points for the Jacobi/DIIS solver
      integer :: jacobi_ndiis = 20
      !> Maximal degree of multipole spherical harmonics
      integer :: pm = 8
      !> Maximal degree of local spherical harmonics
      integer :: pl = 8
      !> Handling of the sparse matrices 
      integer :: incore = 0
      !> 1 to use FMM acceleration and 0 otherwise
      integer :: enable_fmm = 1
   end type ddx_input

   !> Definition of polarizable continuum model
   type, extends(solvation_type) :: ddx_solvation
      !> ddX instance
      type(ddx_input) :: ddx_input
      !> Dielectric function
      real(wp) :: feps
      !> Dielctric constant
      real(wp) :: dielectric_const
      !> Van-der-Waal radii for all atoms
      real(wp), allocatable :: rvdw(:)
   contains
      !> Update cache from container
      procedure :: update
      !> Return dependency on density
      procedure :: variable_info
      !> Get electric field energy
      procedure :: get_energy
      !> Get electric field potential
      procedure :: get_potential
      !> Get electric field gradient
      procedure :: get_gradient
   end type ddx_solvation

   !> Provide constructor for ddX solvation
   interface ddx_solvation
      module procedure :: create_ddx
   end interface ddX_solvation

   !> Restart data for ddX calculation
   type :: ddx_cache
      !> ddX instance
      type(ddx_type) :: ddx 
      !> ddX container with quantities common to all models
      type(ddx_state_type) :: ddx_state
      !> ddX container for the electrostatic properties  
      type(ddx_electrostatics_type) :: ddx_electrostatics
      !> ddX error handling
      type(ddx_error_type) :: ddx_error
      !> Interaction matrix with surface charges jmat(ncav, nat)
      real(wp), allocatable :: jmat(:, :)
      !> ddX potential
      real(wp), allocatable :: ddx_pot(:)
      !> Solvation energy as returned by ddx
      real(wp) :: esolv
      !> ddx multipole, (1, mol%nat)
      real(wp), allocatable :: multipoles(:, :)
      !> ddx forces (i.e. gradient of the solvation energy)
      real(wp), allocatable :: force(:, :)
   end type ddx_cache

contains

!> Create new electric field container
subroutine new_ddx(self, mol, input, error)
   !> Instance of the solvation model
   type(ddx_solvation), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input for ddX solvation
   type(ddx_input), intent(in) :: input
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, izp
   real(wp) :: feps_param 

   ! Set label
   if (input%ddx_model == ddx_solvation_model%cosmo) then
      self%label = "ddcosmo solvation model"
   else if (input%ddx_model == ddx_solvation_model%cpcm) then
      self%label = "ddcpcm solvation model"
   else if (input%ddx_model == ddx_solvation_model%pcm) then
      self%label = "ddpcm solvation model"
   end if

   ! Set model 
   self%ddx_input%ddx_model = input%ddx_model

   ! Get number of OMP threads
   self%ddx_input%nproc = omp_get_max_threads()

   ! Get radii for all atoms
   allocate(self%rvdw(mol%nat), source=0.0_wp)
   if (allocated(input%rvdw)) then
      self%rvdw(:) = input%rscale * input%rvdw(mol%id)
   else
      do iat = 1, mol%nat
         izp = mol%num(mol%id(iat))
         self%rvdw(iat) = input%rscale * get_vdw_rad_cosmo(izp)
      end do
   end if

   ! Get epsilon and calculate feps
   self%dielectric_const = input%dielectric_const
   if (input%ddx_model == ddx_solvation_model%cosmo) then
      feps_param = 0.5_wp
      self%feps = (self%dielectric_const - 1.0_wp) / (self%dielectric_const + feps_param)
   else if (input%ddx_model == ddx_solvation_model%cpcm) then
      feps_param = 0.0_wp
      self%feps = (self%dielectric_const - 1.0_wp) / (self%dielectric_const + feps_param)
   else 
      self%feps = 1.0_wp
   end if

   ! Initialize the shift of the switching function depending on the model: 
   ! ddCOSMO/ddCPCM has an internal shift, ddPCM has a symmetric shift
   if (input%ddx_model == ddx_solvation_model%cosmo .or. &
      & input%ddx_model == ddx_solvation_model%cpcm) then
      self%ddx_input%shift = -1.0_wp
   else 
      self%ddx_input%shift = 0.0_wp
   end if

   ! Initialize the rest of the input
   self%ddx_input%conv = input%conv
   self%ddx_input%eta = input%eta
   self%ddx_input%nang = input%nang
   self%ddx_input%lmax = input%lmax
   self%ddx_input%max_iter = input%max_iter
   self%ddx_input%jacobi_ndiis = input%jacobi_ndiis
   self%ddx_input%pm = input%pm
   self%ddx_input%pl = input%pl
   self%ddx_input%incore = input%incore
   self%ddx_input%enable_fmm = input%enable_fmm

end subroutine new_ddx

!> Type constructor for ddX splvation
function create_ddx(mol, input) result(self)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input for ddX solvation
   type(ddx_input), intent(in) :: input
   !> Instance of the solvation model
   type(ddx_solvation) :: self
   !> Error handling
   type(error_type), allocatable :: error

   ! Create new instance of the solvation model
   call new_ddx(self, mol, input, error)
   if (allocated(error)) then
      call fatal_error(error)
   end if

end function create_ddx


!> Update cache from container
subroutine update(self, mol, cache)
   !> Instance of the solvation model
   class(ddx_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   type(ddx_cache), pointer :: ptr

   integer :: model

   call taint(cache, ptr)

   ! Electrostatics at the cavity points
   ! Electric potential
   ptr%ddx_electrostatics%do_phi = .true.
   ! Electric field
   ptr%ddx_electrostatics%do_e = .true.
   ! Electric field gradient
   ptr%ddx_electrostatics%do_g = .true.

   ! Adjust the model to what ddX expects
   if (self%ddx_input%ddx_model == ddx_solvation_model%cosmo .or. &
      & self%ddx_input%ddx_model == ddx_solvation_model%cpcm) then
      model = 1 ! ddCOSMO and ddCPCM are handeled the same way in ddX
   else
      model = 2 ! ddPCM 
   end if

   call ddinit(model, mol%nat, mol%xyz, self%rvdw, self%dielectric_const, ptr%ddx, ptr%ddx_error, &
      & force=1, ngrid=self%ddx_input%nang, &
      & lmax=self%ddx_input%lmax, nproc=self%ddx_input%nproc, &
      & eta=self%ddx_input%eta, shift=self%ddx_input%shift, &
      & maxiter=self%ddx_input%max_iter, jacobi_ndiis=self%ddx_input%jacobi_ndiis, &
      & pm=self%ddx_input%pm, pl=self%ddx_input%pl, &
      & incore=self%ddx_input%incore, enable_fmm=self%ddx_input%enable_fmm)
   call check_error(ptr%ddx_error)

   call allocate_state(ptr%ddx%params, ptr%ddx%constants, ptr%ddx_state, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   if (allocated(ptr%multipoles)) then
      deallocate(ptr%multipoles)
   end if
   allocate(ptr%multipoles(1, mol%nat), source=0.0_wp)

   if (allocated(ptr%jmat)) then
      deallocate(ptr%jmat)
   end if
   allocate(ptr%jmat(ptr%ddx%constants%ncav, mol%nat), source=0.0_wp)
   call get_coulomb_matrix(mol%xyz, ptr%ddx%constants%ccav, ptr%jmat)

   if (allocated(ptr%force)) then
      deallocate(ptr%force)
   end if
   allocate(ptr%force(3, mol%nat), source=0.0_wp)

   if (allocated(ptr%ddx_pot)) then
      deallocate(ptr%ddx_pot)
   end if
   allocate(ptr%ddx_pot(mol%nat), source=0.0_wp)

   call multipole_electrostatics(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%multipoles, 0, ptr%ddx_electrostatics, ptr%ddx_error)
   call multipole_psi(ptr%ddx%params, ptr%multipoles, 0, ptr%ddx_state%psi)
   
   call setup(ptr%ddx%params,ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, ptr%ddx_electrostatics, &
      & ptr%ddx_state%psi, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   call fill_guess(ptr%ddx%params, ptr%ddx%constants, &
         & ptr%ddx%workspace, ptr%ddx_state, self%ddx_input%conv, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   call fill_guess_adjoint(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, self%ddx_input%conv, ptr%ddx_error)
   call check_error(ptr%ddx_error)

end subroutine update

!> Get electric field energy
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the solvation model
   class(ddx_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Solvation free energy
   real(wp), intent(inout) :: energies(:)
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   type(ddx_cache), pointer :: ptr

   call view(cache, ptr)

   ! Recalculate the solution of the ddX system with the new charges after diagonalization
   ! This solution cannot be reused in the potential due to intermediate mixing

   ptr%multipoles(1, :) = wfn%qat(:, 1) / sqrt(4.0_wp*pi)
   call multipole_electrostatics(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%multipoles, 0, ptr%ddx_electrostatics, ptr%ddx_error)
   call multipole_psi(ptr%ddx%params, ptr%multipoles, 0, ptr%ddx_state%psi)
   
   call setup(ptr%ddx%params,ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, ptr%ddx_electrostatics, &
      & ptr%ddx_state%psi, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   call solve(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, self%ddx_input%conv, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   ! Add solvation energy to total energy
   energies(:) = energies + self%feps * 0.5_wp * sum(ptr%ddx_state%xs * ptr%ddx_state%psi, 1) 

end subroutine get_energy

!> Get electric field potential
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the solvation model
   class(ddx_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   type(ddx_cache), pointer :: ptr

   call view(cache, ptr)

   ! Solution of the ddX system (direct and adjoint) with the mixed charges
   ! This solution cannot be reused in the energy calculation due to intermediate diagonalization

   ptr%multipoles(1, :) = wfn%qat(:, 1) / sqrt(4.0_wp*pi)
   call multipole_electrostatics(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%multipoles, 0, ptr%ddx_electrostatics, ptr%ddx_error)
   call multipole_psi(ptr%ddx%params, ptr%multipoles, 0, ptr%ddx_state%psi)

   call setup(ptr%ddx%params,ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, ptr%ddx_electrostatics, &
      & ptr%ddx_state%psi, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   call solve(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, self%ddx_input%conv, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   call solve_adjoint(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, self%ddx_input%conv, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   ! Contract with the Coulomb matrix
   ptr%ddx_pot = 0.0_wp
   call gemv(ptr%jmat, ptr%ddx_state%zeta, ptr%ddx_pot(:), alpha=-1.0_wp, beta=1.0_wp, trans='t')   
   ! Scale with 0.5 and feps, and get second contribution to potential
   ptr%ddx_pot(:) = 0.5_wp * self%feps * (ptr%ddx_pot(:) + sqrt(4.0_wp*pi) * ptr%ddx_state%xs(1, :))
 
   ! Add potential to overall potential for new SCF step 
   pot%vat(:,1) = pot%vat(:,1) + ptr%ddx_pot(:)

end subroutine get_potential

!> Get electric field gradient
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the solvation model
   class(ddx_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular gradient of the solvation free energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the solvation free energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   type(ddx_cache), pointer :: ptr
   
   call view(cache, ptr)

   ptr%force = 0.0_wp

   ptr%multipoles(1, :) = wfn%qat(:, 1) / sqrt(4.0_wp*pi)
   call multipole_electrostatics(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%multipoles, 0, ptr%ddx_electrostatics, ptr%ddx_error)
   call multipole_psi(ptr%ddx%params, ptr%multipoles, 0, ptr%ddx_state%psi)

   call setup(ptr%ddx%params,ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, ptr%ddx_electrostatics, &
      & ptr%ddx_state%psi, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   call solve(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, self%ddx_input%conv, ptr%ddx_error)
      
   call check_error(ptr%ddx_error)

   call solve_adjoint(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, self%ddx_input%conv, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   call solvation_force_terms(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, ptr%ddx_electrostatics, ptr%force, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   call multipole_force_terms(ptr%ddx%params, ptr%ddx%constants, ptr%ddx%workspace, &
      ptr%ddx_state, 0, ptr%multipoles, ptr%force, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   ! Calculate the gradient of the solvation energy
   ptr%force = self%feps * ptr%force 

   ! Add the gradient of the solvation energy to the total gradient
   gradient =  gradient + ptr%force

end subroutine get_gradient


!> Return dependency on density
pure function variable_info(self) result(info)
   !> Instance of the solvation model
   class(ddx_solvation), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=atom_resolved)
end function variable_info


subroutine taint(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(ddx_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(ddx_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

subroutine view(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(ddx_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(ddx_cache)
      ptr => target
   end select
end subroutine view

!> Evaluate the Coulomb interactions between the atomic sides (xyz) and the
!> surface elements of the cavity (ccav).
subroutine get_coulomb_matrix(xyz, ccav, jmat)
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: ccav(:, :)
   real(wp), intent(inout) :: jmat(:, :)

   integer :: ic, j
   real(wp) :: vec(3), d2, d

   jmat(:, :) = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp shared(ccav, xyz, jmat) private(ic, j, vec, d2, d)
   do ic = 1, size(ccav, 2)
      do jat = 1, size(xyz, 2)
         vec(:) = ccav(:, ic) - xyz(:, j)
         d2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         d = sqrt(d2)
         jmat(ic, j) = 1.0_wp / d
      end do
   end do

end subroutine get_coulomb_matrix

end module tblite_solvation_ddx
