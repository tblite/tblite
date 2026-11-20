/* This file is part of tblite. SPDX-License-Identifier: LGPL-3.0-or-later */
#ifdef USE_CUDA
#include <cublas_v2.h>
#include <cuda_runtime_api.h>
#include <cusolverDn.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

struct gpu_workspace {
    int n, lwork;
    bool single;
    cusolverDnHandle_t solver;
    cublasHandle_t blas;
    void *a, *b, *values, *work;
    int *info;
};

void tblite_cuda_delete(struct gpu_workspace *w)
{
    if (!w) return;
    cudaFree(w->work);
    cudaFree(w->info);
    cudaFree(w->values);
    cudaFree(w->b);
    cudaFree(w->a);
    if (w->blas) cublasDestroy(w->blas);
    if (w->solver) cusolverDnDestroy(w->solver);
    free(w);
}

/* Keep API failures distinct from LAPACK's numerical info. */
#define CUDA_OK(expr) do { int rc = (expr); if (rc) return -1000 - rc; } while (0)
#define SOLVER_OK(expr) do { int rc = (expr); if (rc) return -2000 - rc; } while (0)
#define BLAS_OK(expr) do { int rc = (expr); if (rc) return -3000 - rc; } while (0)

static int prepare(struct gpu_workspace *w, int n, bool single)
{
    const size_t item = single ? sizeof(float) : sizeof(double);
    if (n <= 0 || (size_t)n > SIZE_MAX / item / (size_t)n) return -1;
    const size_t square = (size_t)n * n;
    w->n = n;
    w->single = single;
    SOLVER_OK(cusolverDnCreate(&w->solver));
    BLAS_OK(cublasCreate(&w->blas));
    CUDA_OK(cudaMalloc(&w->a, square * item));
    CUDA_OK(cudaMalloc(&w->b, square * item));
    CUDA_OK(cudaMalloc(&w->values, n * item));
    CUDA_OK(cudaMalloc((void **)&w->info, sizeof(int)));
    if (single) {
        SOLVER_OK(cusolverDnSsygvd_bufferSize(w->solver, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR,
            CUBLAS_FILL_MODE_UPPER, n, w->a, n, w->b, n, w->values, &w->lwork));
    } else {
        SOLVER_OK(cusolverDnDsygvd_bufferSize(w->solver, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR,
            CUBLAS_FILL_MODE_UPPER, n, w->a, n, w->b, n, w->values, &w->lwork));
    }
    if (w->lwork <= 0) return -2;
    /* Density construction reuses this buffer for the weighted coefficients. */
    size_t count = (size_t)w->lwork;
    if (!single && count < square) count = square;
    CUDA_OK(cudaMalloc(&w->work, count * item));
    return 0;
}

int tblite_cuda_sygvd(struct gpu_workspace **handle, int n, int single,
                      void *h, const void *s, void *eval)
{
    struct gpu_workspace *w = *handle;
    if (!w || w->n != n || w->single != (bool)single) {
        tblite_cuda_delete(w);
        *handle = w = calloc(1, sizeof(*w));
        if (!w) return -3;
        int status = prepare(w, n, single);
        if (status) {
            tblite_cuda_delete(w);
            *handle = NULL;
            return status;
        }
    }
    const size_t item = single ? sizeof(float) : sizeof(double);
    const size_t bytes = (size_t)n * n * item;
    CUDA_OK(cudaMemcpy(w->a, h, bytes, cudaMemcpyHostToDevice));
    CUDA_OK(cudaMemcpy(w->b, s, bytes, cudaMemcpyHostToDevice));
    if (single) {
        SOLVER_OK(cusolverDnSsygvd(w->solver, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR,
            CUBLAS_FILL_MODE_UPPER, n, w->a, n, w->b, n, w->values, w->work, w->lwork, w->info));
    } else {
        SOLVER_OK(cusolverDnDsygvd(w->solver, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR,
            CUBLAS_FILL_MODE_UPPER, n, w->a, n, w->b, n, w->values, w->work, w->lwork, w->info));
    }
    int info;
    CUDA_OK(cudaMemcpy(&info, w->info, sizeof(info), cudaMemcpyDeviceToHost));
    if (info) return info;
    CUDA_OK(cudaMemcpy(eval, w->values, n * item, cudaMemcpyDeviceToHost));
    CUDA_OK(cudaMemcpy(h, w->a, bytes, cudaMemcpyDeviceToHost));
    return 0;
}

int tblite_cuda_density(struct gpu_workspace *w, const double *coeff,
                        const double *weight, double *density)
{
    if (!w || w->single) return -1;
    const int n = w->n;
    const size_t bytes = (size_t)n * n * sizeof(double);
    const double one = 1.0, zero = 0.0;
    /* Weighted densities may refer to either spin channel after both solves. */
    CUDA_OK(cudaMemcpy(w->a, coeff, bytes, cudaMemcpyHostToDevice));
    CUDA_OK(cudaMemcpy(w->values, weight, n * sizeof(double), cudaMemcpyHostToDevice));
    BLAS_OK(cublasDdgmm(w->blas, CUBLAS_SIDE_RIGHT, n, n, w->a, n, w->values, 1, w->work, n));
    BLAS_OK(cublasDgemm(w->blas, CUBLAS_OP_N, CUBLAS_OP_T, n, n, n, &one,
        w->work, n, w->a, n, &zero, w->b, n));
    CUDA_OK(cudaMemcpy(density, w->b, bytes, cudaMemcpyDeviceToHost));
    return 0;
}
#endif
