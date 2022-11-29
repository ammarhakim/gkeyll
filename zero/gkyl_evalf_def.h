#pragma once

#include <stddef.h>

/**
 * Type of function to project.
 *
 * @param t Time to evaluate function
 * @param xn Coordinates for evaluation
 * @param fout Output vector of 'num_ret_vals'
 * @param ctx Context for function evaluation. Can be NULL
 */
typedef void (*evalf_t)(double t, const double *xn, double *fout, void *ctx);

/**
 * Type of function to apply BC
 *
 * @param t Time at which BC is applied
 * @param xc Coordinate of skin cell
 * @param ncomp Number of compontents (size of skin and ghost arrays)
 * @param skin Pointer to data in skin-cell
 * @param ghost Pointer to data in ghost-cell
 * @param ctx Context for function evaluation. Can be NULL
 */
typedef void (*wv_bc_func_t)(double t, double xc[3], int ncomp, const double *skin, double *ghost, void *ctx);

/**
 * Type of function for use in array copy op.
 *
 * @param nc Number of elements in @a out and @a inp
 * @param out Output buffer
 * @param inp Input buffer
 * @param ctx Context for function evaluation. Can be NULL
 */
typedef void (*array_copy_func_t)(size_t nc, double *out, const double *inp, void *ctx);
