#include <gkyl_array_reduce.h>
#include <gkyl_array_reduce_priv.h>
#include <assert.h>
#include <float.h>

void 
gkyl_array_reduce(double *out, const struct gkyl_array *arr, enum gkyl_array_op op)
{
  assert(arr->type == GKYL_DOUBLE);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) {
    switch (op) {
      case GKYL_MAX:
        gkyl_array_reduce_max_cu(out, arr);
        break;
      case GKYL_MIN:
        gkyl_array_reduce_min_cu(out, arr);
        break;
      case GKYL_SUM:
        gkyl_array_reduce_sum_cu(out, arr);
        break;
      case GKYL_ABS_MAX:
        gkyl_array_reduce_abs_max_cu(out, arr);
        break;
      case GKYL_SQ_SUM:
        gkyl_array_reduce_sq_sum_cu(out, arr);
        break;
      case GKYL_ABS:
      case GKYL_INV:
      case GKYL_PROD:
      case GKYL_DIV:
      case GKYL_AXPBY:
      case GKYL_RMS:
        assert(false);
        break;
    }
    return;
  }
#endif

  long nc = arr->ncomp;

  switch (op) {
    case GKYL_MIN:
      for (long k=0; k<nc; ++k) out[k] = DBL_MAX;
      for (size_t i=0; i<arr->size; ++i) {
        const double *d = gkyl_array_cfetch(arr, i);
        for (long k=0; k<nc; ++k)
          out[k] = fmin(out[k], d[k]);
      }
      break;

    case GKYL_MAX:
      for (long k=0; k<nc; ++k) out[k] = -DBL_MAX;
      for (size_t i=0; i<arr->size; ++i) {
        const double *d = gkyl_array_cfetch(arr, i);
        for (long k=0; k<nc; ++k)
          out[k] = fmax(out[k], d[k]);
      }
      break;

    case GKYL_SUM:
      for (long k=0; k<nc; ++k) out[k] = 0;
      for (size_t i=0; i<arr->size; ++i) {
        const double *d = gkyl_array_cfetch(arr, i);
        for (long k=0; k<nc; ++k)
          out[k] += d[k];
      }
      break;

    case GKYL_ABS_MAX:
      for (long k=0; k<nc; ++k) out[k] = -DBL_MAX;
      for (size_t i=0; i<arr->size; ++i) {
        const double *d = gkyl_array_cfetch(arr, i);
        for (long k=0; k<nc; ++k)
          out[k] = fmax(out[k], fabs(d[k]));
      }
      break;

    case GKYL_SQ_SUM:
      for (long k=0; k<nc; ++k) out[k] = 0;
      for (size_t i=0; i<arr->size; ++i) {
        const double *d = gkyl_array_cfetch(arr, i);
        for (long k=0; k<nc; ++k)
          out[k] += d[k]*d[k];
      }
      break;

    case GKYL_ABS:
    case GKYL_INV:
    case GKYL_PROD:
    case GKYL_DIV:
    case GKYL_AXPBY:
    case GKYL_RMS:
      assert(false);
      break;
  }
}

void
gkyl_array_reduce_range(double *res,
  const struct gkyl_array *arr, enum gkyl_array_op op, const struct gkyl_range *range)
{
  assert(arr->type == GKYL_DOUBLE);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) {
    switch (op) {
      case GKYL_MAX:
        gkyl_array_reduce_range_max_cu(res, arr, range);
        break;
      case GKYL_MIN:
        gkyl_array_reduce_range_min_cu(res, arr, range);
        break;
      case GKYL_SUM:
        gkyl_array_reduce_range_sum_cu(res, arr, range);
        break;
      case GKYL_ABS_MAX:
        gkyl_array_reduce_range_abs_max_cu(res, arr, range);
        break;
      case GKYL_SQ_SUM:
        gkyl_array_reduce_range_sq_sum_cu(res, arr, range);
        break;
      case GKYL_ABS:
      case GKYL_INV:
      case GKYL_PROD:
      case GKYL_DIV:
      case GKYL_AXPBY:
      case GKYL_RMS:
        assert(false);
        break;
    }
    return;
  }
#endif

  long n = arr->ncomp;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  switch (op) {
    case GKYL_MIN:
      for (long i=0; i<n; ++i) res[i] = DBL_MAX;

      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(range, iter.idx);
        const double *d = gkyl_array_cfetch(arr, start);
        for (long i=0; i<n; ++i)
          res[i] = fmin(res[i], d[i]);
      }
      break;
    case GKYL_MAX:
      for (long i=0; i<n; ++i) res[i] = -DBL_MAX;

      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(range, iter.idx);
        const double *d = gkyl_array_cfetch(arr, start);
        for (long i=0; i<n; ++i)
          res[i] = fmax(res[i], d[i]);
      }
      break;
    case GKYL_SUM:
      for (long i=0; i<n; ++i) res[i] = 0;

      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(range, iter.idx);
        const double *d = gkyl_array_cfetch(arr, start);
        for (long i=0; i<n; ++i)
          res[i] += d[i];
      }
      break;
    case GKYL_ABS_MAX:
      for (long i=0; i<n; ++i) res[i] = -DBL_MAX;

      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(range, iter.idx);
        const double *d = gkyl_array_cfetch(arr, start);
        for (long i=0; i<n; ++i)
          res[i] = fmax(res[i], fabs(d[i]));
      }
      break;
    case GKYL_SQ_SUM:
      for (long i=0; i<n; ++i) res[i] = 0;

      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(range, iter.idx);
        const double *d = gkyl_array_cfetch(arr, start);
        for (long i=0; i<n; ++i)
          res[i] += d[i]*d[i];
      }
      break;
    case GKYL_ABS:
    case GKYL_INV:
    case GKYL_PROD:
    case GKYL_DIV:
    case GKYL_AXPBY:
    case GKYL_RMS:
      assert(false);
      break;
  }
}

void 
gkyl_array_reduce_weighted(double *out, const struct gkyl_array *arr,
  const struct gkyl_array *wgt, enum gkyl_array_op op)
{
  assert(arr->type == GKYL_DOUBLE);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) {
    switch (op) {
      case GKYL_MAX:
        gkyl_array_reduce_weighted_max_cu(out, arr, wgt);
        break;
      case GKYL_MIN:
        gkyl_array_reduce_weighted_min_cu(out, arr, wgt);
        break;
      case GKYL_SUM:
        gkyl_array_reduce_weighted_sum_cu(out, arr, wgt);
        break;
      case GKYL_ABS_MAX:
        gkyl_array_reduce_weighted_abs_max_cu(out, arr, wgt);
        break;
      case GKYL_SQ_SUM:
        gkyl_array_reduce_weighted_sq_sum_cu(out, arr, wgt);
        break;
      case GKYL_RMS:
        gkyl_array_reduce_weighted_rms_cu(out, arr, wgt);
        break;
      case GKYL_ABS:
      case GKYL_INV:
      case GKYL_PROD:
      case GKYL_DIV:
      case GKYL_AXPBY:
        assert(false);
        break;
    }
    return;
  }
#endif

  long nc = arr->ncomp;

  switch (op) {
    case GKYL_MIN:
      for (long k=0; k<nc; ++k) out[k] = DBL_MAX;
      for (size_t i=0; i<arr->size; ++i) {
        const double *arr_c = gkyl_array_cfetch(arr, i);
        const double *wgt_c = gkyl_array_cfetch(wgt, i);
        for (long k=0; k<nc; ++k)
          out[k] = fmin(out[k], wgt_c[k]*arr_c[k]);
      }
      break;

    case GKYL_MAX:
      for (long k=0; k<nc; ++k) out[k] = -DBL_MAX;
      for (size_t i=0; i<arr->size; ++i) {
        const double *arr_c = gkyl_array_cfetch(arr, i);
        const double *wgt_c = gkyl_array_cfetch(wgt, i);
        for (long k=0; k<nc; ++k)
          out[k] = fmax(out[k], wgt_c[k]*arr_c[k]);
      }
      break;

    case GKYL_SUM:
      for (long k=0; k<nc; ++k) out[k] = 0;
      for (size_t i=0; i<arr->size; ++i) {
        const double *arr_c = gkyl_array_cfetch(arr, i);
        const double *wgt_c = gkyl_array_cfetch(wgt, i);
        for (long k=0; k<nc; ++k)
          out[k] += wgt_c[k]*arr_c[k];
      }
      break;
    case GKYL_ABS_MAX:
      for (long k=0; k<nc; ++k) out[k] = -DBL_MAX;
      for (size_t i=0; i<arr->size; ++i) {
        const double *arr_c = gkyl_array_cfetch(arr, i);
        const double *wgt_c = gkyl_array_cfetch(wgt, i);
        for (long k=0; k<nc; ++k)
          out[k] = fmax(out[k], fabs(wgt_c[k]*arr_c[k]));
      }
      break;

    case GKYL_SQ_SUM:
      for (long k=0; k<nc; ++k) out[k] = 0;
      for (size_t i=0; i<arr->size; ++i) {
        const double *arr_c = gkyl_array_cfetch(arr, i);
        const double *wgt_c = gkyl_array_cfetch(wgt, i);
        for (long k=0; k<nc; ++k)
          out[k] += wgt_c[k]*arr_c[k]*arr_c[k];
      }
      break;

    case GKYL_RMS:
      for (long k=0; k<nc; ++k) out[k] = 0;
      for (size_t i=0; i<arr->size; ++i) {
        const double *arr_c = gkyl_array_cfetch(arr, i);
        const double *wgt_c = gkyl_array_cfetch(wgt, i);
	double rms = 0.0;
        for (long k=0; k<nc; ++k)
          rms += arr_c[k]*arr_c[k];

	rms = sqrt(rms/nc);
        for (long k=0; k<nc; ++k)
          out[k] += wgt_c[k]*rms;
      }
      break;
    case GKYL_ABS:
    case GKYL_INV:
    case GKYL_PROD:
    case GKYL_DIV:
    case GKYL_AXPBY:
      assert(false);
      break;
  }
}

void
gkyl_array_reduce_weighted_range(double *res, const struct gkyl_array *arr,
  const struct gkyl_array *wgt, enum gkyl_array_op op, const struct gkyl_range *range)
{
  assert(arr->type == GKYL_DOUBLE);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) {
    switch (op) {
      case GKYL_MAX:
        gkyl_array_reduce_weighted_range_max_cu(res, arr, wgt, range);
        break;
      case GKYL_MIN:
        gkyl_array_reduce_weighted_range_min_cu(res, arr, wgt, range);
        break;
      case GKYL_SUM:
        gkyl_array_reduce_weighted_range_sum_cu(res, arr, wgt, range);
        break;
      case GKYL_ABS_MAX:
        gkyl_array_reduce_weighted_range_abs_max_cu(res, arr, wgt, range);
        break;
      case GKYL_SQ_SUM:
        gkyl_array_reduce_weighted_range_sq_sum_cu(res, arr, wgt, range);
        break;
      case GKYL_RMS:
        gkyl_array_reduce_weighted_range_rms_cu(res, arr, wgt, range);
        break;
      case GKYL_ABS:
      case GKYL_INV:
      case GKYL_PROD:
      case GKYL_DIV:
      case GKYL_AXPBY:
        assert(false);
        break;
    }
    return;
  }
#endif

  long n = arr->ncomp;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  switch (op) {
    case GKYL_MIN:
      for (long i=0; i<n; ++i) res[i] = DBL_MAX;

      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(range, iter.idx);
        const double *arr_c = gkyl_array_cfetch(arr, start);
        const double *wgt_c = gkyl_array_cfetch(wgt, start);
        for (long i=0; i<n; ++i)
          res[i] = fmin(res[i], wgt_c[i]*arr_c[i]);
      }
      break;
    case GKYL_MAX:
      for (long i=0; i<n; ++i) res[i] = -DBL_MAX;

      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(range, iter.idx);
        const double *arr_c = gkyl_array_cfetch(arr, start);
        const double *wgt_c = gkyl_array_cfetch(wgt, start);
        for (long i=0; i<n; ++i)
          res[i] = fmax(res[i], wgt_c[i]*arr_c[i]);
      }
      break;
    case GKYL_SUM:
      for (long i=0; i<n; ++i) res[i] = 0;

      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(range, iter.idx);
        const double *arr_c = gkyl_array_cfetch(arr, start);
        const double *wgt_c = gkyl_array_cfetch(wgt, start);
        for (long i=0; i<n; ++i)
          res[i] += wgt_c[i]*arr_c[i];
      }
      break;
    case GKYL_ABS_MAX:
      for (long i=0; i<n; ++i) res[i] = -DBL_MAX;

      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(range, iter.idx);
        const double *arr_c = gkyl_array_cfetch(arr, start);
        const double *wgt_c = gkyl_array_cfetch(wgt, start);
        for (long i=0; i<n; ++i)
          res[i] = fmax(res[i], fabs(wgt_c[i]*arr_c[i]));
      }
      break;
    case GKYL_SQ_SUM:
      for (long i=0; i<n; ++i) res[i] = 0;

      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(range, iter.idx);
        const double *arr_c = gkyl_array_cfetch(arr, start);
        const double *wgt_c = gkyl_array_cfetch(wgt, start);
        for (long i=0; i<n; ++i)
          res[i] += wgt_c[i]*arr_c[i]*arr_c[i];
      }
      break;
    case GKYL_RMS:
      for (long i=0; i<n; ++i) res[i] = 0;

      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(range, iter.idx);
        const double *arr_c = gkyl_array_cfetch(arr, start);
        const double *wgt_c = gkyl_array_cfetch(wgt, start);
	double rms = 0.0;
        for (long i=0; i<n; ++i)
          rms += pow(arr_c[i],2.0);

	rms = sqrt(rms/n);
        for (long i=0; i<n; ++i)
          res[i] += wgt_c[i]*rms*rms;
      }
      break;
    case GKYL_ABS:
    case GKYL_INV:
    case GKYL_PROD:
    case GKYL_DIV:
    case GKYL_AXPBY:
      assert(false);
      break;
  }
}

