! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
module test_gpu_solver
   use mctc_env, only : sp, dp, error_type
   use mctc_env_testing, only : new_unittest, unittest_type, check
   use tblite_cusolver_sygvd, only : sygvd_cusolver
   use tblite_features, only : tblite_has_cuda
   implicit none
   private
   public :: collect_gpu_solver
contains

subroutine collect_gpu_solver(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [new_unittest("generalized-eigensystem", test_sygvd)]
end subroutine collect_gpu_solver

subroutine test_sygvd(error)
   type(error_type), allocatable, intent(out) :: error
   type(error_type), allocatable :: solver_error
   type(sygvd_cusolver) :: solver
   real(dp), parameter :: h(2,2) = reshape([2.0_dp, 0.3_dp, 0.3_dp, 4.0_dp], [2,2])
   real(dp), parameter :: weight(2) = [-0.8_dp, 0.4_dp]
   real(dp) :: s(2,2), coeff(2,2), eval(2), density(2,2), expected(2,2)
   real(sp) :: cs(2,2), es(2)

   s = reshape([1.5_dp, 0.1_dp, 0.1_dp, 1.0_dp], [2,2])
   coeff = h
   call solver%solve(coeff, s, eval, solver_error)
   if (.not.tblite_has_cuda) then
      call check(error, allocated(solver_error), "CPU build must report an unavailable GPU solver")
      return
   end if
   if (allocated(solver_error)) then
      call move_alloc(solver_error, error)
      return
   end if
   call check(error, maxval(abs(matmul(h, coeff) - matmul(s, coeff)*spread(eval, 1, 2))) < 1.0e-12_dp)
   if (allocated(error)) return
   ! A different host coefficient matrix catches reuse of a previous spin channel.
   coeff(:, 1) = [0.4_dp, 0.6_dp]
   expected = matmul(coeff*spread(weight, 1, 2), transpose(coeff))
   call solver%get_density_matrix(weight, coeff, density, error)
   if (allocated(error)) return
   call check(error, maxval(abs(density - expected)) < 1.0e-12_dp)
   if (allocated(error)) return

   cs = real(h, sp)
   call solver%solve(cs, real(s, sp), es, error)
   if (allocated(error)) return
   call check(error, maxval(abs(real(es,dp) - eval)) < 2.0e-6_dp)
   if (allocated(error)) return

   s(1,1) = -1.0_dp
   coeff = h
   call solver%solve(coeff, s, eval, solver_error)
   call check(error, allocated(solver_error), "Expected a non-positive-definite overlap error")
end subroutine test_sygvd
end module test_gpu_solver
