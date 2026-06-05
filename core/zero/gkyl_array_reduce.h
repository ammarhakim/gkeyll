#pragma once

#include <float.h>
#include <gkyl_array.h>
#include <gkyl_range.h>

/**
 * Perform an "reduce" operation of data in the array.
 *
 * @param res On output, reduces values (ncomp size).
 * @param arr Array to perform reduction on.
 * @param op Reduction operators.
 */
void gkyl_array_reduce(double *res, const struct gkyl_array *arr, enum gkyl_array_op op);

/**
 * Perform an "reduce" operation of data in the array in a specified range.
 * Data is reduced component-wise.
 *
 * @param res On output, reduced values (ncomp size).
 * @param arr Array to perform reduction on.
 * @param op Reduction operators.
 * @param range Range specifying region.
 */
void gkyl_array_reduce_range(double *res,
  const struct gkyl_array *arr, enum gkyl_array_op op, const struct gkyl_range *range);

/**
 * Perform a weighted reduction of data in the array.
 *
 * @param res On output, reduces values (ncomp size).
 * @param arr Array to perform reduction on.
 * @param wgt Weight to include in the reduction.
 * @param op Reduction operators.
 */
void gkyl_array_reduce_weighted(double *out, const struct gkyl_array *arr,
  const struct gkyl_array *wgt, enum gkyl_array_op op);

/**
 * Perform a weighted reduction of data in the array in a specified range.
 *
 * @param res On output, reduced values (ncomp size).
 * @param arr Array to perform reduction on.
 * @param op Reduction operators.
 * @param range Range specifying region.
 */
void gkyl_array_reduce_weighted_range(double *res, const struct gkyl_array *arr,
  const struct gkyl_array *wgt, enum gkyl_array_op op, const struct gkyl_range *range);

