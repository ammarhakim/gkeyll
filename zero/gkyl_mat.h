#pragma once

#include <gkyl_ref_count.h>
#include <gkyl_util.h>

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

/** Flags for indicating (conjugate) transpose */
enum gkyl_mat_trans { GKYL_NO_TRANS, GKYL_TRANS, GKYL_CONJ_TRANS };

/**
 * Matrix object. Stored in column major order.
 */
struct gkyl_mat {
  size_t nr, nc; // Number of rows, columns
  double *data; // Pointer to data
};

/**
 * Multiple matrices, each stored in column major order
 */
struct gkyl_nmat {
  size_t num; // Number of matrices
  size_t nr, nc; // Number of rows, columns
  double *data; // Pointer to data
  double **mptr; // pointers to start of each sub-matrix

  uint32_t flags;  
  struct gkyl_ref_count ref_count;
  struct  gkyl_nmat *on_dev; // pointer to itself or device data
};

// Type for storing preallocating memory needed in various batch
// operations
typedef struct gkyl_nmat_mem gkyl_nmat_mem;

/**
 * Construct new matrix with all elements initialized to @a
 * val. Delete using gkyl_mat_release method.
 *
 * @param nr Number of rows
 * @param nc Number of cols
 * @param val Initial value
 * @return Pointer to new matrix.
 */
struct gkyl_mat* gkyl_mat_new(size_t nr, size_t nc, double val);

/**
 * Clone matrix.
 */
struct gkyl_mat* gkyl_mat_clone(const struct gkyl_mat *in);

/**
 * Set value in matrix.
 */
GKYL_CU_DH
static inline void
gkyl_mat_set(struct gkyl_mat *mat, size_t r, size_t c, double val)
{
  mat->data[c*mat->nr+r] = val;
}

/**
 * Set value in matrix.
 */
GKYL_CU_DH
static inline void
gkyl_mat_inc(struct gkyl_mat *mat, size_t r, size_t c, double val)
{
  mat->data[c*mat->nr+r] += val;
}

/**
 * Get value from matrix.
 */
GKYL_CU_DH
static inline double
gkyl_mat_get(const struct gkyl_mat *mat, size_t r, size_t c)
{
  return mat->data[c*mat->nr+r];
}

/**
 * Get column of matrix as const pointer.
 */
GKYL_CU_DH
static inline const double*
gkyl_mat_get_ccol(const struct gkyl_mat *mat, size_t c)
{
  return mat->data+c*mat->nr;
}

/**
 * Get column of matrix as pointer.
 */
GKYL_CU_DH
static inline double*
gkyl_mat_get_col(struct gkyl_mat *mat, size_t c)
{
  return mat->data+c*mat->nr;
}

/**
 * Set all elements of matrix to specified value. Returns pointer to @a mat.
 */
GKYL_CU_DH
static inline struct gkyl_mat*
gkyl_mat_clear(struct gkyl_mat *mat, double val)
{
  for (size_t i=0; i<mat->nr*mat->nc; ++i) mat->data[i] = val;
  return mat;
}

/**
 * Set all elements on diagonal to specified value. All other elements
 * are set to 0.0. Returns pointer to @a mat.
 */ 
struct gkyl_mat* gkyl_mat_diag(struct gkyl_mat *mat, double val);

/**
 * Write matrix to file. Output is in Maxima matrix format
 */
void gkyl_mat_show(const char *name, FILE *fp, const struct gkyl_mat *mat);

/**
 * Computes matrix-matrix product:
 *
 * C = alpha*OP(A)*OP(B) + beta*C
 *
 * where OP(A) indicates transpose/no-transpose based on the
 * transa/transb flags.
 *
 * C is returned
 */
struct gkyl_mat* gkyl_mat_mm(double alpha, double beta,
  enum gkyl_mat_trans transa, const struct gkyl_mat *A,
  enum gkyl_mat_trans transb, const struct gkyl_mat *B, struct gkyl_mat *C);

/**
 * Solve system of linear equations using LU decomposition. On input
 * the RHS must be in the "x" matrix (each column represents a RHS
 * vector) and on output "x" is replaced with the solution(s). Returns
 * true on success, false otherwise. Note that on output A is replaced
 * by its LU factors.
 *
 * The ipiv input is an chunk of memory that is sizeof(lapack_int[N]),
 * where N is the number of equations. It is safest to assume
 * lapack_int is long (it may be smaller). You must
 * allocate/deallocate ipiv yourself! Use:
 *
 * ipiv = gkyl_mem_buff_new(sizeof(long[N]));
 * gkyl_mat_linsolve_lu(A, x, gkyl_mem_buff_data(ipiv));
 * gkyl_mem_buff_release(ipiv);
 *
 * The reason for passing ipiv to this function is that it avoids
 * allocations inside this function.
 */
bool gkyl_mat_linsolve_lu(struct gkyl_mat *A, struct gkyl_mat *x, void* ipiv);

/**
 * Release matrix
 *
 * @param mat Pointer to matrix to release
 */
void gkyl_mat_release(struct gkyl_mat *mat);

//////////////// gkyl_nmat API

/**
 * Construct new multi-matrix (batch of matrices). Delete using
 * gkyl_nmat_release method. Each matrix has the same shape.
 *
 * @param num Number of matrices
 * @param nr Number of rows
 * @param nc Number of cols
 * @return Pointer to new multi-matrix.
 */
struct gkyl_nmat *gkyl_nmat_new(size_t num, size_t nr, size_t nc);

/**
 * Construct new multi-matrix (batch of matrices). Delete using
 * gkyl_nmat_release method. Each matrix has the same shape.
 *
 * CAUTION: The nmat returned by this method lives on the GPU. You
 * CAN'T modify it directly on the host! If you try, it will crash the
 * program.
 *
 * NOTE: the data member lives on GPU, but the struct lives on the
 * host.  However, the on_dev member for this cal is set to a device
 * clone of the host struct, and is what should be used to pass to
 * CUDA kernels which require the entire array struct on device.
 * 
 * @param num Number of matrices
 * @param nr Number of rows
 * @param nc Number of cols
 * @return Pointer to new multi-matrix.
 */
struct gkyl_nmat* gkyl_nmat_cu_dev_new(size_t num, size_t nr, size_t nc);

/**
 * Copy into nmat: pointer to dest nmat is returned. 'dest' and 'src'
 * must not point to same data.
 *
 * @param dest Destination for copy.
 * @param src Srouce to copy from.
 * @return dest is returned
 */
struct gkyl_nmat* gkyl_nmat_copy(struct gkyl_nmat *dest, const struct gkyl_nmat *src);

/**
 * Get a matrix from multi-matrix. DO NOT free the returned matrix!
 *
 * @param n Matrix to fetch
 * @return Matrix (DO NOT free/release this)
 */
GKYL_CU_DH
static inline struct gkyl_mat
gkyl_nmat_get(struct gkyl_nmat *mat, size_t num)
{
  return (struct gkyl_mat) {
    .nr = mat->nr,
    .nc = mat->nc,
    .data = mat->data+num*mat->nr*mat->nc
  };
}

/**
 * Returns true if multi-matrix lives on NV-GPU.
 *
 * @param mat Multi-matrix to check
 * @return true if on NV-GPU, false otherwise
 */
bool gkyl_nmat_is_cu_dev(const struct gkyl_nmat *mat);

/**
 * Acquire pointer to multi-matrix. The pointer must be released using
 * gkyl_nmar_release method.
 *
 * @param mat Multi-matrix to which a pointer is needed
 * @return Pointer to acquired multi-matrix.
 */
struct gkyl_nmat* gkyl_nmat_acquire(const struct gkyl_nmat *mat);

/**
 * Allocate memory needed in batched LU solves on host. Free using the
 * release method.
 *
 * @param num Number of matrices in batch.
 * @param nrow Number of rows in each matrix
 * @return Preallocated memory
 */
gkyl_nmat_mem *gkyl_nmat_linsolve_lu_new(size_t num, size_t nrow);
// Same as above, except for GPUs
gkyl_nmat_mem *gkyl_nmat_linsolve_lu_cu_dev_new(size_t num, size_t nrow);

/**
 * Release memory allocated for batched LU solves.
 *
 * @param mem Memory to release
 */
void gkyl_nmat_linsolve_lu_release(gkyl_nmat_mem *mem);

/**
 * Solve a batched system of linear equations using LU
 * decomposition. On input the RHSs must be in the "x" multi-matrix
 * (each column represents a RHS vector) and on output "x" is replaced
 * with the solution(s). Returns true on success, false
 * otherwise. Note that on output each of the As is replaced by its LU
 * factors.
 */
bool gkyl_nmat_linsolve_lu(struct gkyl_nmat *A, struct gkyl_nmat *x);

/**
 * Solve a batched system of linear equations using LU
 * decomposition. On input the RHSs must be in the "x" multi-matrix
 * (each column represents a RHS vector) and on output "x" is replaced
 * with the solution(s). Returns true on success, false
 * otherwise. Note that on output each of the As is replaced by its LU
 * factors.
 *
 * The memory required in this call must be pre-allocated.
 *
 * @param mem Preallocated memory needed in the solve
 * @param A list of LHS matrices, replaced by LU factors on return
 * @param x list of RHS vectors, replace by solution in exit.
 */
bool gkyl_nmat_linsolve_lu_pa(gkyl_nmat_mem *mem, struct gkyl_nmat *A, struct gkyl_nmat *x);

/**
 * Release multi-matrix
 *
 * @param mat Pointer to multi-matrix to release
 */
void gkyl_nmat_release(struct gkyl_nmat *nmat);
