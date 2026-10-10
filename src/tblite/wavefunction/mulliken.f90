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

!> @file tblite/wavefunction/mulliken.f90
!> Provides Mulliken population analysis

!> Wavefunction analysis via Mulliken populations
module tblite_wavefunction_mulliken
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_blas, only: gemm
   use tblite_partition, only : work_partition, owns_index
   use tblite_wavefunction_spin, only : updown_to_magnet
   implicit none
   private

   public :: contract_mayer_columns
   public :: get_mulliken_shell_charges, get_mulliken_atomic_multipoles
   public :: get_molecular_dipole_moment, get_molecular_quadrupole_moment
   public :: get_mayer_bond_orders, get_mulliken_shell_multipoles

contains


subroutine get_mulliken_shell_charges(bas, smat, pmat, n0sh, qsh, partition, columns)
   type(basis_type), intent(in) :: bas
   real(wp), intent(in) :: smat(:, :)
   real(wp), intent(in) :: pmat(:, :, :)
   real(wp), intent(in) :: n0sh(:)
   real(wp), intent(out) :: qsh(:, :)
   !> Optional ownership of output shells; contributions add across partitions
   type(work_partition), intent(in), optional :: partition
   integer, intent(in), optional :: columns(2)
   integer :: first, last, offset

   integer :: ish, ii, iao, jao, spin
   real(wp) :: psh

   first = 1
   last = bas%nao
   if (present(columns)) then
      first = columns(1)
      last = columns(2)
   end if
   offset = first - 1

   qsh(:, :) = 0.0_wp
   !$omp parallel do default(none) collapse(2) schedule(runtime) &
   !$omp shared(qsh, bas, pmat, smat, partition, columns, first, last, offset) private(spin, ish, ii, iao, jao, psh)
   do spin = 1, size(pmat, 3)
      do ish = 1, bas%nsh
         if (.not.present(columns)) then
            if (.not.owns_index(partition, ish)) cycle
         end if
         psh = 0.0_wp
         ii = bas%iao_sh(ish)
         do iao = max(1, first-ii), min(bas%nao_sh(ish), last-ii)
            do jao = 1, bas%nao
               psh = psh + pmat(jao, ii+iao-offset, spin) * smat(jao, ii+iao-offset)
            end do
         end do
         qsh(ish, spin) = -psh
      end do
   end do

   call updown_to_magnet(qsh)
   do ish = 1, bas%nsh
      if (present(columns)) then
         ! Assign the neutral shell occupation once, even when a shell is split.
         if (bas%iao_sh(ish)+1 < first .or. bas%iao_sh(ish)+1 > last) cycle
      else
         if (.not.owns_index(partition, ish)) cycle
      end if
      qsh(ish, 1) = qsh(ish, 1) + n0sh(ish)
   end do

end subroutine get_mulliken_shell_charges


subroutine get_mulliken_atomic_multipoles(bas, mpmat, pmat, mpat, partition, columns)
   type(basis_type), intent(in) :: bas
   real(wp), intent(in) :: mpmat(:, :, :)
   real(wp), intent(in) :: pmat(:, :, :)
   real(wp), intent(out) :: mpat(:, :, :)
   !> Optional ownership of output atoms; contributions add across partitions
   type(work_partition), intent(in), optional :: partition
   integer, intent(in), optional :: columns(2)
   integer :: first, last, offset

   integer ::  iat, ish, is, ii, iao, jao, spin
   real(wp) :: pat(size(mpmat, 1))

   first = 1
   last = bas%nao
   if (present(columns)) then
      first = columns(1)
      last = columns(2)
   end if
   offset = first - 1

   mpat(:, :, :) = 0.0_wp
   !$omp parallel do default(none) collapse(2) schedule(runtime) &
   !$omp shared(bas, pmat, mpmat, mpat, partition, columns, first, last, offset) &
   !$omp private(spin, iat, ish, is, ii, iao, jao, pat)
   do spin = 1, size(pmat, 3)
      do iat = 1, size(mpat, 2)
         if (.not.present(columns)) then
            if (.not.owns_index(partition, iat)) cycle
         end if
         pat(:) = 0.0_wp
         is = bas%ish_at(iat)
         do ish = 1, bas%nsh_at(iat)
            ii = bas%iao_sh(is + ish)
            do iao = max(1, first-ii), min(bas%nao_sh(is + ish), last-ii)
               do jao = 1, bas%nao
                  pat(:) = pat(:) + pmat(jao, ii+iao-offset, spin) * mpmat(:, jao, ii+iao-offset)
               end do
            end do
         end do
         mpat(:, iat, spin) = -pat(:)
      end do
   end do

   call updown_to_magnet(mpat)

end subroutine get_mulliken_atomic_multipoles

subroutine get_mulliken_shell_multipoles(bas, mpmat, pmat, mpsh)
   type(basis_type), intent(in) :: bas
   real(wp), intent(in) :: mpmat(:, :, :)
   real(wp), intent(in) :: pmat(:, :, :)
   real(wp), intent(out) :: mpsh(:, :, :)

   integer :: iao, jao, spin, ish, ii
   real(wp) :: psh(size(mpmat, 1))

   mpsh(:, :, :) = 0.0_wp
   !$omp parallel do default(none) collapse(2) schedule(runtime) &
   !$omp shared(bas, pmat, mpmat, mpsh) &
   !$omp private(spin, ish, iao, jao, ii, psh)
   do spin = 1, size(pmat, 3)
      do ish = 1, bas%nsh
         psh(:) = 0.0_wp
         ii = bas%iao_sh(ish)
         do iao = 1, bas%nao_sh(ish)
            do jao = 1, bas%nao
               psh(:) = psh + pmat(jao, ii+iao, spin) * mpmat(:, jao, ii+iao)
            end do
         end do
         mpsh(:, ish, spin) = -psh(:)
      end do
   end do

   call updown_to_magnet(mpsh)

end subroutine get_mulliken_shell_multipoles

subroutine get_molecular_dipole_moment(mol, qat, dpat, dpmom)
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: qat(:)
   real(wp), intent(in) :: dpat(:, :)
   real(wp), intent(out) :: dpmom(:)

   integer :: iat

   dpmom(:) = 0.0_wp
   do iat = 1, mol%nat
      dpmom(:) = dpmom + mol%xyz(:, iat) * qat(iat) + dpat(:, iat)
   end do
end subroutine get_molecular_dipole_moment

subroutine get_molecular_quadrupole_moment(mol, qat, dpat, qpat, qpmom)
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: qat(:)
   real(wp), intent(in) :: dpat(:, :)
   real(wp), intent(in) :: qpat(:, :)
   real(wp), intent(out) :: qpmom(:)

   integer :: iat
   real(wp) :: vec(3), cart(6), tr

   qpmom(:) = 0.0_wp
   do iat = 1, mol%nat
      vec(:) = mol%xyz(:, iat)*qat(iat)
      cart([1, 3, 6]) = mol%xyz(:, iat) * (vec + 2*dpat(:, iat))
      cart(2) = mol%xyz(1, iat) * (vec(2) + dpat(2, iat)) + dpat(1, iat)*mol%xyz(2, iat)
      cart(4) = mol%xyz(1, iat) * (vec(3) + dpat(3, iat)) + dpat(1, iat)*mol%xyz(3, iat)
      cart(5) = mol%xyz(2, iat) * (vec(3) + dpat(3, iat)) + dpat(2, iat)*mol%xyz(3, iat)
      tr = 0.5_wp * (cart(1) + cart(3) + cart(6))
      cart(1) = 1.5_wp * cart(1) - tr
      cart(2) = 3.0_wp * cart(2)
      cart(3) = 1.5_wp * cart(3) - tr
      cart(4) = 3.0_wp * cart(4)
      cart(5) = 3.0_wp * cart(5)
      cart(6) = 1.5_wp * cart(6) - tr
      qpmom(:) = qpmom(:) + qpat(:, iat) + cart
   end do
end subroutine get_molecular_quadrupole_moment


!> Evaluate Wiberg/Mayer bond orders including factor 2 scaling for open-shell cases
subroutine get_mayer_bond_orders(mol, bas, smat, pmat, mbo, partition)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Overlap matrix
   real(wp), intent(in) :: smat(:, :)
   !> Density matrix
   real(wp), intent(in) :: pmat(:, :, :)
   !> Wiberg/Mayer bond orders
   real(wp), intent(out) :: mbo(:, :, :)
   !> Optional ownership of the outer atom of each unordered pair
   type(work_partition), intent(in), optional :: partition

   integer :: i, nlocal, spin
   integer, allocatable :: columns(:)
   real(wp), allocatable, target :: pscols(:, :), row_storage(:, :)
   real(wp), contiguous, pointer :: psrows(:, :)

   mbo = 0.0_wp
   columns = pack([(i, i=1, bas%nao)], owns_index(partition, bas%ao2at))
   nlocal = size(columns)
   if (nlocal == 0) return
   allocate(pscols(bas%nao, nlocal))
   if (nlocal < bas%nao) then
      allocate(row_storage(nlocal, bas%nao))
      psrows => row_storage
   else
      ! Complete storage needs only one product for both orientations.
      psrows => pscols
   end if
   do spin = 1, size(pmat, 3)
      if (nlocal < bas%nao) then
         call gemm(pmat(:, :, spin), smat(:, columns), pscols)
         call gemm(pmat(columns, :, spin), smat, psrows)
      else
         call gemm(pmat(:, :, spin), smat, pscols)
      end if
      call contract_mayer_columns(bas, pscols, psrows, mbo(:, :, spin), columns)
   end do
   if (size(pmat, 3) == 2) mbo = 2*mbo
   call updown_to_magnet(mbo)
end subroutine get_mayer_bond_orders


!> Contract owned PS columns with the corresponding rows into atom pairs.
!> Distributed callers supply the rows as columns of the transposed product.
subroutine contract_mayer_columns(bas, ps, pst, mbo, columns, transposed)
   type(basis_type), intent(in) :: bas
   real(wp), intent(in) :: ps(:, :), pst(:, :)
   real(wp), intent(out) :: mbo(:, :)
   integer, intent(in) :: columns(:)
   logical, intent(in), optional :: transposed
   integer :: iat, jat, col, jao
   real(wp) :: row_value
   logical :: rows_as_columns

   rows_as_columns = .false.
   if (present(transposed)) rows_as_columns = transposed
   mbo = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(bas, ps, pst, mbo, columns, rows_as_columns) private(iat, jat, col, jao, row_value)
   do iat = 1, size(mbo, 2)
      do col = 1, size(columns)
         if (bas%ao2at(columns(col)) /= iat) cycle
         do jao = 1, bas%nao
            jat = bas%ao2at(jao)
            if (jat >= iat) cycle
            if (rows_as_columns) then
               row_value = pst(jao, col)
            else
               row_value = pst(col, jao)
            end if
            mbo(jat, iat) = mbo(jat, iat) + ps(jao, col)*row_value
         end do
      end do
      mbo(iat, :iat-1) = mbo(:iat-1, iat)
   end do
end subroutine contract_mayer_columns

end module tblite_wavefunction_mulliken
