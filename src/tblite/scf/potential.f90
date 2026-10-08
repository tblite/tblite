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

!> @file tblite/scf/potential.f90
!> Defines a container for holding density dependent potentials

!> Implementation of the density dependent potential and its contribution
!> to the effective Hamiltonian
module tblite_scf_potential
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_integral_type, only : integral_type
   use tblite_mpi_utils, only : mpi_allreduce_sum, new_mpi_work_partition, &
      & mpi_transpose_columns, mpi_gather_columns
   use tblite_partition, only : work_partition
   use tblite_wavefunction_spin, only : magnet_to_updown
   implicit none
   private

   public :: new_potential, add_pot_to_h1, mpi_add_pot_to_h1, reduce_potential


   !> Container for density dependent potential-shifts
   type, public :: potential_type
      !> Flag if potential includes potential gradients
      logical :: grad = .false.
      !> Atom-resolved charge-dependent potential shift
      real(wp), allocatable :: vat(:, :)
      !> Shell-resolved charge-dependent potential shift
      real(wp), allocatable :: vsh(:, :)
      !> Orbital-resolved charge-dependent potential shift
      real(wp), allocatable :: vao(:, :)

      !> Atom-resolved dipolar potential
      real(wp), allocatable :: vdp(:, :, :)
      !> Atom-resolved quadrupolar potential
      real(wp), allocatable :: vqp(:, :, :)

      !> Reusable workspace for assembling all potential components together
      real(wp), allocatable :: reduction_buffer(:)

      !> Position derivative of atom-resolved charge-dependent potential shift
      real(wp), allocatable :: dvatdr(:, :, :, :)
      !> Lattice vector derivative of atom-resolved charge-dependent potential shift
      real(wp), allocatable :: dvatdL(:, :, :, :)

      !> Position derivative of shell-resolved charge-dependent potential shift
      real(wp), allocatable :: dvshdr(:, :, :, :)
      !> Lattice vector derivative of shell-resolved charge-dependent potential shift
      real(wp), allocatable :: dvshdL(:, :, :, :)
   contains
      !> Reset the density dependent potential
      procedure :: reset
   end type potential_type


contains


!> Create a new potential object
subroutine new_potential(self, mol, bas, nspin, grad)
   !> Instance of the density dependent potential
   type(potential_type), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Description of the basis set
   type(basis_type), intent(in) :: bas
   !> Number of spin channels
   integer, intent(in) :: nspin
   !> Flag to indicate if potential gradients are requested
   logical, intent(in), optional :: grad

   allocate(self%vat(mol%nat, nspin))
   allocate(self%vsh(bas%nsh, nspin))
   allocate(self%vao(bas%nao, nspin))

   allocate(self%vdp(3, mol%nat, nspin))
   allocate(self%vqp(6, mol%nat, nspin))

   if(present(grad)) then
      if(grad) then
         self%grad = .true.
         allocate(self%dvatdr(3, mol%nat, mol%nat, nspin))
         allocate(self%dvatdL(3, 3, mol%nat, nspin))

         allocate(self%dvshdr(3, mol%nat, bas%nsh, nspin))
         allocate(self%dvshdL(3, 3, bas%nsh, nspin))
      end if
   end if

end subroutine new_potential

!> Reset the density dependent potential
subroutine reset(self)
   !> Instance of the density dependent potential
   class(potential_type), intent(inout) :: self

   self%vat(:, :) = 0.0_wp
   self%vsh(:, :) = 0.0_wp
   self%vao(:, :) = 0.0_wp
   self%vdp(:, :, :) = 0.0_wp
   self%vqp(:, :, :) = 0.0_wp

   if(self%grad) then
      self%dvatdr(:, :, :, :) = 0.0_wp
      self%dvatdL(:, :, :, :) = 0.0_wp

      self%dvshdr(:, :, :, :) = 0.0_wp
      self%dvshdL(:, :, :, :) = 0.0_wp
   end if
end subroutine reset

!> Assemble the SCF potential before expanding atom/shell shifts into vao.
subroutine reduce_potential(error, pot, comm)
   !> Error handling; a pending error skips the reduction
   type(error_type), allocatable, intent(inout) :: error
   !> Partial atomic and shell potentials, replaced by their sum over all ranks
   type(potential_type), intent(inout) :: pot
   !> Communicator to reduce over, absent leaves the potential unchanged
   integer, intent(in), optional :: comm
   integer :: na, ns, nd, nq, n

   if (.not.present(comm) .or. allocated(error)) return
   na = size(pot%vat)
   ns = size(pot%vsh)
   nd = size(pot%vdp)
   nq = size(pot%vqp)
   n = na + ns + nd + nq
   if (allocated(pot%reduction_buffer)) then
      if (size(pot%reduction_buffer) /= n) deallocate(pot%reduction_buffer)
   end if
   if (.not.allocated(pot%reduction_buffer)) allocate(pot%reduction_buffer(n))
   pot%reduction_buffer(:na) = reshape(pot%vat, [na])
   pot%reduction_buffer(na+1:na+ns) = reshape(pot%vsh, [ns])
   pot%reduction_buffer(na+ns+1:na+ns+nd) = reshape(pot%vdp, [nd])
   pot%reduction_buffer(na+ns+nd+1:) = reshape(pot%vqp, [nq])
   call mpi_allreduce_sum(error, pot%reduction_buffer, comm)
   if (allocated(error)) return
   pot%vat = reshape(pot%reduction_buffer(:na), shape(pot%vat))
   pot%vsh = reshape(pot%reduction_buffer(na+1:na+ns), shape(pot%vsh))
   pot%vdp = reshape(pot%reduction_buffer(na+ns+1:na+ns+nd), shape(pot%vdp))
   pot%vqp = reshape(pot%reduction_buffer(na+ns+nd+1:), shape(pot%vqp))
end subroutine reduce_potential

!> Expand the replicated potential into disjoint Hamiltonian columns.
!> Local storage needs a distributed transpose of one-centre contributions;
!> the compatibility path gathers complete columns for a replicated solver.
subroutine mpi_add_pot_to_h1(error, bas, ints, pot, h1, comm)
   !> Error handling
   type(error_type), allocatable, intent(inout) :: error
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Integrals stored as complete matrices or complete local columns
   type(integral_type), intent(in) :: ints
   !> Replicated potential shifts, expanded to shell and orbital resolution
   type(potential_type), intent(inout) :: pot
   !> Effective Hamiltonian to overwrite, using the same column layout as ints
   real(wp), contiguous, intent(inout) :: h1(:, :, :)
   !> Communicator for distributed assembly, absent selects serial execution
   integer, intent(in), optional :: comm

   type(work_partition) :: partition
   real(wp), allocatable :: transposed(:, :)
   integer, allocatable :: columns(:)
   integer :: spin

   if (allocated(error)) return
   if (present(comm)) call new_mpi_work_partition(error, partition, comm)
   if (allocated(error)) return
   if (.not.ints%local .and. partition%get_nparts() > 1) columns = partition%get_columns(size(h1, 2))
   call add_pot_to_h1(bas, ints, pot, h1, columns)
   if (ints%local) then
      allocate(transposed(size(h1, 1), size(h1, 2)))
      do spin = 1, size(h1, 3)
         call mpi_transpose_columns(error, h1(:, :, spin), transposed, comm)
         if (allocated(error)) return
         h1(:, :, spin) = h1(:, :, spin) + transposed
      end do
   else if (allocated(columns)) then
      do spin = 1, size(h1, 3)
         call mpi_gather_columns(error, b=h1(:, :, spin), n=size(h1, 2), comm=comm)
      end do
   end if
end subroutine mpi_add_pot_to_h1


!> Add the collected potential shifts to complete or locally stored columns.
!> Local storage returns one-centre contributions; the frontend adds their
!> distributed transpose. Complete storage includes both centres directly.
subroutine add_pot_to_h1(bas, ints, pot, h1, columns)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Integrals stored as complete matrices or complete local columns
   type(integral_type), intent(in) :: ints
   !> Potential shifts, expanded to shell and orbital resolution
   type(potential_type), intent(inout) :: pot
   !> Hamiltonian columns to overwrite; local storage returns one-centre terms
   real(wp), intent(out) :: h1(:, :, :)
   !> Optional range to compute within replicated storage
   integer, intent(in), optional :: columns(2)

   integer :: spin, col, iao, jao, iat, jat, first, last, offset
   real(wp) :: hij, scale
   logical :: symmetric

   call add_vat_to_vsh(bas, pot%vat, pot%vsh)
   call add_vsh_to_vao(bas, pot%vsh, pot%vao)
   first = 1
   last = size(h1, 2)
   if (present(columns)) then
      first = columns(1)
      last = columns(2)
   end if
   offset = ints%columns(1) - 1
   scale = merge(0.5_wp, 1.0_wp, ints%local)
   symmetric = .not.ints%local .and. .not.present(columns)

   !$omp parallel do collapse(2) schedule(runtime) default(none) &
   !$omp shared(bas, ints, pot, h1, first, last, offset, scale, symmetric) &
   !$omp private(spin, col, iao, jao, iat, jat, hij)
   do spin = 1, size(h1, 3)
      do col = first, last
         iao = col + offset
         iat = bas%ao2at(iao)
         do jao = 1, merge(iao, bas%nao, symmetric)
            hij = 0.0_wp
            if (spin == 1) hij = ints%hamiltonian(jao, col)
            hij = scale*(hij - 0.5_wp*ints%overlap(jao, col)*(pot%vao(jao, spin) + pot%vao(iao, spin))) &
               & - 0.5_wp*dot_product(ints%dipole(:, jao, col), pot%vdp(:, iat, spin)) &
               & - 0.5_wp*dot_product(ints%quadrupole(:, jao, col), pot%vqp(:, iat, spin))
            if (.not.ints%local) then
               jat = bas%ao2at(jao)
               hij = hij - 0.5_wp*dot_product(ints%dipole(:, iao, jao), pot%vdp(:, jat, spin)) &
                  & - 0.5_wp*dot_product(ints%quadrupole(:, iao, jao), pot%vqp(:, jat, spin))
            end if
            h1(jao, col, spin) = hij
            if (symmetric) h1(iao, jao, spin) = hij
         end do
      end do
   end do
   call magnet_to_updown(h1(:, first:last, :))
end subroutine add_pot_to_h1

!> Expand an atom-resolved potential shift to a shell-resolved potential shift
subroutine add_vat_to_vsh(bas, vat, vsh)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Atom-resolved charge-dependent potential shift
   real(wp), intent(in) :: vat(:, :)
   !> Shell-resolved charge-dependent potential shift
   real(wp), intent(inout) :: vsh(:, :)

   integer :: iat, ish, ii, spin

   !$omp parallel do schedule(runtime) collapse(2) default(none) &
   !$omp shared(bas, vat, vsh) private(spin, ii, ish, iat)
   do spin = 1, size(vat, 2)
      do iat = 1, size(vat, 1)
         ii = bas%ish_at(iat)
         do ish = 1, bas%nsh_at(iat)
            vsh(ii+ish, spin) = vsh(ii+ish, spin) + vat(iat, spin)
         end do
      end do
   end do
end subroutine add_vat_to_vsh

!> Expand a shell-resolved potential shift to an orbital-resolved potential shift
subroutine add_vsh_to_vao(bas, vsh, vao)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Shell-resolved charge-dependent potential shift
   real(wp), intent(in) :: vsh(:, :)
   !> Orbital-resolved charge-dependent potential shift
   real(wp), intent(inout) :: vao(:, :)

   integer :: ish, iao, ii, spin

   !$omp parallel do schedule(runtime) collapse(2) default(none) &
   !$omp shared(bas, vsh, vao) private(ii, iao, ish)
   do spin = 1, size(vsh, 2)
      do ish = 1, size(vsh, 1)
         ii = bas%iao_sh(ish)
         do iao = 1, bas%nao_sh(ish)
            vao(ii+iao, spin) = vao(ii+iao, spin) + vsh(ish, spin)
         end do
      end do
   end do
end subroutine add_vsh_to_vao


end module tblite_scf_potential
