#pragma once

#include <gkyl_array_reduce.h>

#ifdef GKYL_HAVE_CUDA

/**
 * Reduce a gkyl_array component-wise, using the max operator.
 *
 * @param inp A gkyl_array to be reduced.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_max_cu(double *out_d, const struct gkyl_array* inp);

/**
 * Reduce a gkyl_array component-wise over a specified range, using the max operator.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_range_max_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_range *range);

/**
 * Reduce a gkyl_array component-wise, using the min operator.
 *
 * @param inp A gkyl_array to be reduced.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_min_cu(double *out_d, const struct gkyl_array* inp);

/**
 * Reduce a gkyl_array component-wise over a specified range, using the min operator.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_range_min_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_range *range);

/**
 * Reduce a gkyl_array component-wise, using the sum operator.
 *
 * @param inp A gkyl_array to be reduced.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_sum_cu(double *out_d, const struct gkyl_array* inp);

/**
 * Reduce a gkyl_array component-wise over a specified range, using the sum operator.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_range_sum_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_range *range);

/**
 * Reduce a gkyl_array component-wise, using the abs and max operators.
 *
 * @param inp A gkyl_array to be reduced.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_abs_max_cu(double *out_d, const struct gkyl_array* inp);

/**
 * Reduce a gkyl_array component-wise over a specified range, using the abs and max operators.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_range_abs_max_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_range *range);

/**
 * Reduce a gkyl_array component-wise, using the sq sum operator.
 *
 * @param inp A gkyl_array to be reduced.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_sq_sum_cu(double *out_d, const struct gkyl_array* inp);

/**
 * Reduce a gkyl_array component-wise over a specified range, using the sq sum operator.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_range_sq_sum_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_range *range);

/**
 * Reduce a gkyl_array component-wise, using the rms operator.
 *
 * @param inp A gkyl_array to be reduced.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_rms_cu(double *out_d, const struct gkyl_array* inp);

/**
 * Reduce a gkyl_array component-wise over a specified range, using the rms operator.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_range_rms_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_range *range);

/**
 * Reduce a weighted gkyl_array component-wise, using the max operator.
 *
 * @param inp A gkyl_array to be reduced.
 * @param wgt A gkyl_array weight.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_weighted_max_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_array* wgt);

/**
 * Reduce a weighted gkyl_array component-wise over a specified range, using the max operator.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param wgt A gkyl_array weight.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_weighted_range_max_cu(double *out_d, const struct gkyl_array* inp,
  const struct gkyl_array* wgt, const struct gkyl_range *range);

/**
 * Reduce a weighted gkyl_array component-wise, using the min operator.
 *
 * @param inp A gkyl_array to be reduced.
 * @param wgt A gkyl_array weight.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_weighted_min_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_array* wgt);

/**
 * Reduce a weighted gkyl_array component-wise over a specified range, using the min operator.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param wgt A gkyl_array weight.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_weighted_range_min_cu(double *out_d, const struct gkyl_array* inp,
  const struct gkyl_array* wgt, const struct gkyl_range *range);

/**
 * Reduce a weighted gkyl_array component-wise, using the sum operator.
 *
 * @param inp A gkyl_array to be reduced.
 * @param wgt A gkyl_array weight.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_weighted_sum_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_array* wgt);

/**
 * Reduce a weighted gkyl_array component-wise over a specified range, using the sum operator.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param wgt A gkyl_array weight.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_weighted_range_sum_cu(double *out_d, const struct gkyl_array* inp,
  const struct gkyl_array* wgt, const struct gkyl_range *range);

/**
 * Reduce a weighted gkyl_array component-wise, using the abs and max operators.
 *
 * @param inp A gkyl_array to be reduced.
 * @param wgt A gkyl_array weight.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_weighted_abs_max_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_array* wgt);

/**
 * Reduce a weighted gkyl_array component-wise over a specified range, using the abs and max operators.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param wgt A gkyl_array weight.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_weighted_range_abs_max_cu(double *out_d, const struct gkyl_array* inp,
  const struct gkyl_array* wgt, const struct gkyl_range *range);

/**
 * Reduce a weighted gkyl_array component-wise, using the sq sum operator.
 *
 * @param inp A gkyl_array to be reduced.
 * @param wgt A gkyl_array weight.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_weighted_sq_sum_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_array* wgt);

/**
 * Reduce a weighted gkyl_array component-wise over a specified range, using the sq sum operator.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param wgt A gkyl_array weight.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_weighted_range_sq_sum_cu(double *out_d, const struct gkyl_array* inp,
  const struct gkyl_array* wgt, const struct gkyl_range *range);

/**
 * Reduce a weighted gkyl_array component-wise, using the RMS operator.
 *
 * @param inp A gkyl_array to be reduced.
 * @param wgt A gkyl_array weight.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_weighted_rms_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_array* wgt);

/**
 * Reduce a weighted gkyl_array component-wise over a specified range, using the RMS operator.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param wgt A gkyl_array weight.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_weighted_range_rms_cu(double *out_d, const struct gkyl_array* inp,
  const struct gkyl_array* wgt, const struct gkyl_range *range);

#endif
