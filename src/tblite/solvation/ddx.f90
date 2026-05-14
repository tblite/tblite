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
!> Provides a polarizable continuum model in domain decomposition (dd) framework of ddX

!> Implicit solvation model based on a polarizable dielectric continuum in domain
!> decomposition framework of ddX
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
      ! COSMO, CPCM, and PCM are internally defined as 100, 101, and 200 to get the correct
      ! dielectric factor (feps) during initialization and to avoid collisions with other aliases
      ! elsewhere in the code
      ! Labels are mapped to ddX-specific values 1 and 2 before passing the model input to the ddX routine
      !> Conductor-like screening model
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
      !> ddX model
      integer :: ddx_model 
      !> Dielectric constant
      real(wp) :: dielectric_const
      !> van der Waals radii for all atoms
      real(wp), allocatable :: rvdw(:)
      !> Scaling of van-der-Waals radii
      real(wp) :: rscale = 1.0_wp
      !> Number of grid points for each atom (default=110)
      integer :: nang = grid_size(8)
      !> Regularization parameter / width of the switching function 
      real(wp) :: eta = 0.1_wp
      !> Maximum angular momentum of basis functions
      integer :: lmax = 1
   end type ddx_input

   !> Provide constructor for ddX input
   interface ddx_input
      module procedure :: create_ddx_input
   end interface ddx_input

   !> Definition of polarizable continuum model
   type, extends(solvation_type) :: ddx_solvation
      !> ddX model 
      integer :: ddx_model
      !> Dielectric function
      real(wp) :: feps
      !> Dielectric constant
      real(wp) :: dielectric_const
      !> van der Waals radii for all atoms
      real(wp), allocatable :: rvdw(:)
      !> Accuracy for iterative solver
      real(wp) :: conv = 1.0e-10_wp
      !> Regularization parameter
      real(wp) :: eta 
      !> Number of grid points for each atom (=110)
      integer :: nang 
      !> Maximum angular momentum of basis functions
      integer :: lmax 
      !> Number of OMP threads
      integer :: nproc = 1
      !> Shift of the switching function 
      ! (default value depends on the model, for COSMO/CPCM it is -1)
      real(wp) :: shift
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
      !> ddX multipoles, dim=(1, mol%nat)
      real(wp), allocatable :: multipoles(:, :)
   end type ddx_cache

contains

!> Constructor for ddX input
function create_ddx_input(ddx_model, dielectric_const, rvdw, rscale, nang, eta, lmax) result(self)
   !> ddX model
   integer, intent(in), optional :: ddx_model
   !> Dielectric constant
   real(wp), intent(in) :: dielectric_const
   !> van der Waals radii for all atoms
   real(wp), intent(in), optional :: rvdw(:)
   !> Scaling of van-der-Waals radii
   real(wp), intent(in), optional :: rscale
   !> Number of grid points for each atom
   integer, intent(in), optional :: nang
   !> Regularization parameter
   real(wp), intent(in), optional :: eta
   !> Maximum angular momentum of basis functions
   integer, intent(in), optional :: lmax

   type(ddx_input) :: self


   if (present(ddx_model)) then
      self%ddx_model = ddx_model
   end if

   self%dielectric_const = dielectric_const

   if (present(rvdw)) then
      self%rvdw = rvdw
   end if

   if (present(rscale)) then
      self%rscale = rscale
   end if

   if (present(nang)) then
      self%nang = nang
   end if

   if (present(eta)) then
      self%eta = eta
   end if

   if (present(lmax)) then
      self%lmax = lmax
   end if

end function create_ddx_input

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
   self%ddx_model = input%ddx_model

   ! Get number of OMP threads
   self%nproc = omp_get_max_threads()

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

   ! Get epsilon and calculate dielectric function depending on the model
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
   ! In ddX, ddCOSMO/ddCPCM uses an internal shift by default, ddPCM has a symmetric shift
   if (input%ddx_model == ddx_solvation_model%cosmo .or. &
      & input%ddx_model == ddx_solvation_model%cpcm) then
      self%shift = -1.0_wp
   else 
      self%shift = 0.0_wp
   end if

   ! Initialize the rest of the adjustable input parameters
   self%nang = input%nang
   self%eta = input%eta
   self%lmax = input%lmax

end subroutine new_ddx

!> Type constructor for ddX solvation
function create_ddx(mol, input, error) result(self)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input for ddX solvation
   type(ddx_input), intent(in) :: input
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Instance of the solvation model
   type(ddx_solvation) :: self

   ! Create new instance of the solvation model
   call new_ddx(self, mol, input, error)
   if (allocated(error)) return

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

   ! Request all electrostatic quantities that may be needed later
   ! The ddX multipole routine allocates and fills only the selected arrays
   ! Electric potential at the cavity points
   ptr%ddx_electrostatics%do_phi = .true.
   ! Electric field at the cavity points
   ptr%ddx_electrostatics%do_e = .true.
   ! Electric field gradient at the cavity points
   ptr%ddx_electrostatics%do_g = .true.

   ! Adjust the tblite model labels to the public ddX API:
   ! model 1 is the COSMO-type equation, model 2 is the PCM equation
   if (self%ddx_model == ddx_solvation_model%cosmo .or. &
      & self%ddx_model == ddx_solvation_model%cpcm) then
      model = 1 ! ddCOSMO and ddCPCM are handled the same way in ddX
   else
      model = 2 ! ddPCM
   end if

   ! Initialize the ddX model:
   ! ddX uses two main kinds of data, labeled as model- and state-specific.
   ! In ddinit, we initialize the model class, where we specify the desired model,
   ! dielectric permittivity, general information on the discretization (such as grid size),
   ! and cavity parameters (including number of spheres, coordinates, and radii)
   ! Moreover, geometry-dependent constants (for instance the number of solvent-exposed
   ! Lebedev grid points) are pre-computed and stored.
   call ddinit(model, mol%nat, mol%xyz, self%rvdw, self%dielectric_const, ptr%ddx, ptr%ddx_error, &
      & force=1, ngrid=self%nang, &
      & lmax=self%lmax, nproc=self%nproc, &
      & eta=self%eta, shift=self%shift, &
      & maxiter=self%max_iter, jacobi_ndiis=self%jacobi_ndiis, &
      & pm=self%pm, pl=self%pl, &
      & incore=self%incore, enable_fmm=self%enable_fmm)
   call check_error(ptr%ddx_error)

   ! The state class contains all quantities that are explicitly solute-dependent. That will be the 
   ! right-hand sides of the primal & adjoint linear systems (the spherical harmonics expansion of
   ! the solute potential and the spherical harmonics representation of the solute charges)
   ! as well as the relevant intermediates.
   ! First we start by allocating the state-specific arrays using the dimensions now defined via ddinit
   call allocate_state(ptr%ddx%params, ptr%ddx%constants, ptr%ddx_state, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   ! Allocate the multipole array that later contains the xTB Mulliken charges
   ! (We only pass monopoles)
   if (.not.allocated(ptr%multipoles))then
         allocate(ptr%multipoles(1, mol%nat), source=0.0_wp) 
   endif

   ! Allocate and compute the Coulomb matrix
   if (allocated(ptr%jmat))then
      deallocate(ptr%jmat)
   endif
   allocate(ptr%jmat(ptr%ddx%constants%ncav, mol%nat), source=0.0_wp)
   call get_coulomb_matrix(mol%xyz, ptr%ddx%constants%ccav, ptr%jmat)

   ! Allocate atom-resolved solvation potential contribution produced by get_potential
   if (.not.allocated(ptr%ddx_pot))then
      allocate(ptr%ddx_pot(mol%nat), source=0.0_wp) 
   endif

   ! Now we initialize the RHSs (they are zero in the first iteration)
   ! Given the monopole distribution...
   ! ...we compute the electrostatic potential on the cavity surface
   call multipole_electrostatics(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%multipoles, 0, ptr%ddx_electrostatics, ptr%ddx_error)
   ! ...we convert the monopole charges into the spherical harmonics representation
   call multipole_psi(ptr%ddx%params, ptr%multipoles, 0, ptr%ddx_state%psi)
   ! ...and we transform the potential into the spherical harmonics representation,
   !    marking the RHSs as ready for the solver
   call setup(ptr%ddx%params,ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, ptr%ddx_electrostatics, &
      & ptr%ddx_state%psi, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   ! Construct model-specific initial guesses for the solutions of the
   ! primal and adjoint linear systems
   call fill_guess(ptr%ddx%params, ptr%ddx%constants, &
         & ptr%ddx%workspace, ptr%ddx_state, self%conv, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   call fill_guess_adjoint(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, self%conv, ptr%ddx_error)
   call check_error(ptr%ddx_error)

end subroutine update

!> Get solvation energy
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

   ! Assemble the multipole array based on the new xTB charges and
   ! apply a normalization for the spherical harmonics representation
   ptr%multipoles(1, :) = wfn%qat(:, 1) / sqrt(4.0_wp*pi)
   ! Compute the electrostatic potential on the cavity surface
   call multipole_electrostatics(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%multipoles, 0, ptr%ddx_electrostatics, ptr%ddx_error)
   ! Convert the normalized monopole charges into the spherical harmonics representation
   call multipole_psi(ptr%ddx%params, ptr%multipoles, 0, ptr%ddx_state%psi)
   ! Transform the potential into the spherical harmonics representation
   call setup(ptr%ddx%params,ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, ptr%ddx_electrostatics, &
      & ptr%ddx_state%psi, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   ! Solve the primal linear system to get the surface charge distribution
   ! for the current solute potential
   call solve(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, self%conv, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   ! Add solvation energy to total energy
   energies(:) = energies + self%feps * 0.5_wp * sum(ptr%ddx_state%xs * ptr%ddx_state%psi, 1) 

end subroutine get_energy

!> Get solvation potential
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

   ! Due to intermediate mixing, the xTB charges have changed since the last call
   ! to get_energy, and we need to solve the ddX systems again

   ! Assemble the multipole array based on the new xTB charges and
   ! apply a normalization for the spherical harmonics representation
   ptr%multipoles(1, :) = wfn%qat(:, 1) / sqrt(4.0_wp*pi)
   ! Compute the electrostatic potential on the cavity surface
   call multipole_electrostatics(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%multipoles, 0, ptr%ddx_electrostatics, ptr%ddx_error)
   ! Convert the normalized monopole charges into the spherical harmonics representation
   call multipole_psi(ptr%ddx%params, ptr%multipoles, 0, ptr%ddx_state%psi)
   ! Transform the potential into the spherical harmonics representation
   call setup(ptr%ddx%params,ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, ptr%ddx_electrostatics, &
      & ptr%ddx_state%psi, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   ! Solve the primal linear system to get the surface charge distribution
   ! for the current solute potential
   call solve(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, self%conv, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   ! Solve the adjoint linear system:
   ! We need to do this any time we compute a derivative of the energy (here wrt to the basis
   ! function coefficients), because this lets us avoid computing the derivative of
   ! the surface charge distribution.
   ! For reference, see J. Chem. Theory Comput. 2013, 9, 3637−3648
   call solve_adjoint(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, self%conv, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   ptr%ddx_pot = 0.0_wp
   ! The solvation potential comes from a contribution of the derivative of Psi as well as
   ! from a contribution of the derivative of the cavity potential & coupling matrix.
   ! (Again, see J. Chem. Theory Comput. 2013, 9, 3637−3648)
   ! Zeta is an intermediate in the computation of the latter and needs to be
   ! contracted with the Coulomb matrix in a last step.
   call gemv(ptr%jmat, ptr%ddx_state%zeta, ptr%ddx_pot(:), alpha=-1.0_wp, beta=1.0_wp, trans='t')  

   ! Scale with 0.5 and feps, and get the Psi contribution to potential
   ptr%ddx_pot(:) = 0.5_wp * self%feps * (ptr%ddx_pot(:) + sqrt(4.0_wp*pi) * ptr%ddx_state%xs(1, :))
   ! Add potential to overall potential for new SCC step
   pot%vat(:,1) = pot%vat(:,1) + ptr%ddx_pot(:)

end subroutine get_potential

!> Get solvation gradient
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

   !> Temporary variable for the ddX force/gradient
   real(wp), allocatable :: force(:,:)

   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   type(ddx_cache), pointer :: ptr
   
   call view(cache, ptr)

   allocate(force(3, mol%nat), source=0.0_wp)

   ! Compute all the solute-aspecific solvation force terms
   call solvation_force_terms(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, ptr%ddx_electrostatics, force, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   ! Compute all the solute-specific solvation force terms
   call multipole_force_terms(ptr%ddx%params, ptr%ddx%constants, ptr%ddx%workspace, &
      ptr%ddx_state, 0, ptr%multipoles, force, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   ! Add the dielectric factor to the gradient of the solvation energy
   ! (ddX computes the gradient without it)
   force = self%feps * force 

   ! Add the gradient of the solvation energy to the total gradient
   gradient =  gradient + force

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

!> Evaluate the Coulomb interactions between the atomic sites (xyz) and the
!> surface elements of the cavity (ccav)
subroutine get_coulomb_matrix(xyz, ccav, jmat)
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: ccav(:, :)
   real(wp), intent(inout) :: jmat(:, :)

   integer :: ic, jat
   real(wp) :: vec(3), d2, d

   jmat(:, :) = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp shared(ccav, xyz, jmat) private(ic, jat, vec, d2, d)
   do ic = 1, size(ccav, 2)
      do jat = 1, size(xyz, 2)
         vec(:) = ccav(:, ic) - xyz(:, jat)
         d2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         d = sqrt(d2)
         jmat(ic, jat) = 1.0_wp / d
      end do
   end do

end subroutine get_coulomb_matrix

end module tblite_solvation_ddx
