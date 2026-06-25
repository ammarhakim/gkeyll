#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

bool
gkyl_array_copy_func_is_cu_dev(const struct gkyl_array_copy_func *bc)
{
  return GKYL_IS_CU_ALLOC(bc->flags);
}

struct gkyl_array*
gkyl_array_clear(struct gkyl_array* out, double val)
{
  assert(out->type == GKYL_DOUBLE);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {gkyl_array_clear_cu(out, val); return out; }
#endif

  double *out_d = out->data;
  for (size_t i=0; i<NELM(out); ++i)
    out_d[i] = val;
  return out;
}

struct gkyl_array*
gkyl_array_randomize(struct gkyl_array* out, double lo, double hi)
{
  assert(out->type == GKYL_DOUBLE);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_randomize_cu(out, lo, hi); return out; }
#endif

  pcg64_random_t rng = gkyl_pcg64_init(false);
  double *out_d = out->data;
  double range = hi - lo;
  for (size_t i=0; i<NELM(out); ++i)
    out_d[i] = lo + range * gkyl_pcg64_rand_double(&rng);
  return out;
}

struct gkyl_array*
gkyl_array_accumulate(struct gkyl_array* out, double a,
  const struct gkyl_array* inp)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->size == inp->size && out->elemsz == inp->elemsz);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out) && gkyl_array_is_cu_dev(inp)) { gkyl_array_accumulate_cu(out, a, inp); return out; }
#endif

  double *out_d = out->data;
  const double *inp_d = inp->data;
  for (size_t i=0; i<NELM(out); ++i)
    out_d[i] += a*inp_d[i];
  return out;
}

struct gkyl_array*
gkyl_array_accumulate_offset(struct gkyl_array* out, double a,
  const struct gkyl_array* inp, int coff)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->size == inp->size);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out) && gkyl_array_is_cu_dev(inp)) { gkyl_array_accumulate_offset_cu(out, a, inp, coff); return out; }
#endif

  double *out_d = out->data;
  const double *inp_d = inp->data;
  if (NCOM(out) < NCOM(inp)) {
    // Interpret offset as offset in input components.
    for (size_t i=0; i<out->size; ++i)
      for (size_t c=0; c<NCOM(out); ++c)
        out_d[i*NCOM(out)+c] += a*inp_d[i*NCOM(inp)+c+coff];
  } else {
    // Interpret offset as offset in output components.
    for (size_t i=0; i<out->size; ++i)
      for (size_t c=0; c<NCOM(inp); ++c)
        out_d[i*NCOM(out)+c+coff] += a*inp_d[i*NCOM(inp)+c];
  }
  return out;
}

struct gkyl_array*
gkyl_array_set(struct gkyl_array* out, double a,
  const struct gkyl_array* inp)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->size == inp->size && out->elemsz == inp->elemsz);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_set_cu(out, a, inp); return out; }
#endif

  double *out_d = out->data;
  const double *inp_d = inp->data;
  for (size_t i=0; i<NELM(out); ++i)
    out_d[i] = a*inp_d[i];
  return out;
}

struct gkyl_array*
gkyl_array_set_offset(struct gkyl_array* out, double a,
  const struct gkyl_array* inp, int coff)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->size == inp->size);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_set_offset_cu(out, a, inp, coff); return out; }
#endif

  double *out_d = out->data;
  const double *inp_d = inp->data;
  if (NCOM(out) < NCOM(inp)) {
    // Interpret offset as offset in input components.
    for (size_t i=0; i<out->size; ++i)
      for (size_t c=0; c<NCOM(out); ++c)
        out_d[i*NCOM(out)+c] = a*inp_d[i*NCOM(inp)+c+coff];
  } else {
    // Interpret offset as offset in output components.
    for (size_t i=0; i<out->size; ++i)
      for (size_t c=0; c<NCOM(inp); ++c)
        out_d[i*NCOM(out)+c+coff] = a*inp_d[i*NCOM(inp)+c];
  }
  return out;
}

struct gkyl_array*
gkyl_array_scale(struct gkyl_array* out, double a)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_scale_cu(out, a); return out; }
#endif

  return gkyl_array_set(out, a, out);
}

struct gkyl_array*
gkyl_array_scale_by_cell(struct gkyl_array* out, const struct gkyl_array* a)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->size == a->size && NCOM(a) == 1);
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_scale_by_cell_cu(out, a); return out; }
#endif

  double *out_d = out->data;
  const double *a_d = a->data;
  for (size_t i=0; i<out->size; ++i)
    for (size_t c=0; c<NCOM(out); ++c)
      out_d[i*NCOM(out)+c] = a_d[i]*out_d[i*NCOM(out)+c];
  return out;
}

struct gkyl_array*
gkyl_array_divide_by_cell(struct gkyl_array* out, const struct gkyl_array* a)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->size == a->size && NCOM(a) == 1);
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_divide_by_cell_cu(out, a); return out; }
#endif

  double *out_d = out->data;
  const double *a_d = a->data;
  for (size_t i=0; i<out->size; ++i)
    for (size_t c=0; c<NCOM(out); ++c)
      out_d[i*NCOM(out)+c] = out_d[i*NCOM(out)+c]/a_d[i];
  return out;
}

struct gkyl_array*
gkyl_array_shiftc(struct gkyl_array* out, double a, unsigned k)
{
  assert(out->type == GKYL_DOUBLE);
  assert(k < NCOM(out));
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_shiftc_cu(out, a, k); return out; }
#endif

  double *out_d = out->data;
  for (size_t i=0; i<out->size; ++i)
    out_d[i*NCOM(out)+k] = a+out_d[i*NCOM(out)+k];
  return out;
}

struct gkyl_array*
gkyl_array_comp_op(struct gkyl_array *out, enum gkyl_array_op op, double a,
 const struct gkyl_array *in1, double b, const struct gkyl_array *in2)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->ncomp == in1->ncomp);
  assert(out->size == in1->size);
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_comp_op_cu(out, op, a, in1, b, in2); return out; }
#endif

  double *out_d = out->data;
  int nc = out->ncomp;
  long sz = out->size;

  switch (op) {
    case GKYL_ABS:
      for (size_t i=0; i<sz; ++i) {
        double *out_c = gkyl_array_fetch(out, i);
        const double *in1_c = gkyl_array_cfetch(in1, i);
        for (size_t k=0; k<nc; ++k)
          out_c[k] = fabs(a*in1_c[k]);
      }
      break;                                            
    case GKYL_INV:
      for (size_t i=0; i<sz; ++i) {
        double *out_c = gkyl_array_fetch(out, i);
        const double *in1_c = gkyl_array_cfetch(in1, i);
        for (size_t k=0; k<nc; ++k)
          out_c[k] = a/in1_c[k];
      }
      break;                                            
    case GKYL_PROD:
      assert(out->ncomp == in2->ncomp);
      assert(out->size == in2->size);
      for (size_t i=0; i<sz; ++i) {
        double *out_c = gkyl_array_fetch(out, i);
        const double *in1_c = gkyl_array_cfetch(in1, i);
        const double *in2_c = gkyl_array_cfetch(in2, i);
        for (size_t k=0; k<nc; ++k)
          out_c[k] = a*in1_c[k]*in2_c[k]+b;
      }
      break;                                            
    case GKYL_DIV:
      assert(out->ncomp == in2->ncomp);
      assert(out->size == in2->size);
      for (size_t i=0; i<sz; ++i) {
        double *out_c = gkyl_array_fetch(out, i);
        const double *in1_c = gkyl_array_cfetch(in1, i);
        const double *in2_c = gkyl_array_cfetch(in2, i);
        for (size_t k=0; k<nc; ++k)
          out_c[k] = a*in1_c[k]/in2_c[k]+b;
      }
      break;                                            
    case GKYL_AXPBY:
      assert(out->ncomp == in2->ncomp);
      assert(out->size == in2->size);
      for (size_t i=0; i<sz; ++i) {
        double *out_c = gkyl_array_fetch(out, i);
        const double *in1_c = gkyl_array_cfetch(in1, i);
        const double *in2_c = gkyl_array_cfetch(in2, i);
        for (size_t k=0; k<nc; ++k)
          out_c[k] = a*in1_c[k]+b*in2_c[k];
      }
      break;                                            
    case GKYL_MAX:
    case GKYL_MIN:
    case GKYL_SUM:
    case GKYL_ABS_MAX:
    case GKYL_SQ_SUM:
    case GKYL_RMS:
      assert(false);
      break;
  }
  return out;
}

// range based methods
struct gkyl_array*
gkyl_array_clear_range(struct gkyl_array *out, double val, const struct gkyl_range *range)
{
  assert(out->type == GKYL_DOUBLE);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_clear_range_cu(out, val, range); return out; }
#endif

  long n = NCOM(out);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    array_clear1(n, gkyl_array_fetch(out, start), val);
  }

  return out;
}

struct gkyl_array*
gkyl_array_accumulate_range(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, const struct gkyl_range *range)
{
  assert(out->type == GKYL_DOUBLE);

  assert(out->size == inp->size);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_accumulate_range_cu(out, a, inp, range); return out; }
#endif

  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n = outnc<inpnc ? outnc : inpnc;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    array_acc1(n,
      gkyl_array_fetch(out, start), a, gkyl_array_cfetch(inp, start));
  }

  return out;
}

struct gkyl_array*
gkyl_array_accumulate_offset_range(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, int coff, const struct gkyl_range *range)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->size == inp->size);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_accumulate_offset_range_cu(out, a, inp, coff, range); return out; }
#endif

  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n;
  int outoff, inoff;
  if (outnc < inpnc) {
    n = outnc;
    outoff = 0;
    inoff = coff;
  } else {
    n = inpnc;
    outoff = coff;
    inoff = 0;
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    double *out_d = gkyl_array_fetch(out, start);
    const double *inp_d = gkyl_array_cfetch(inp, start);
    array_acc1(n, out_d+outoff, a, inp_d+inoff);
  }

  return out;
}

struct gkyl_array*
gkyl_array_set_range(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, const struct gkyl_range *range)
{
  assert(out->type == GKYL_DOUBLE && inp->type == GKYL_DOUBLE);
  assert(out->size == inp->size);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_set_range_cu(out, a, inp, range); return out; }
#endif

  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n = outnc<inpnc ? outnc : inpnc;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    array_set1(n,
      gkyl_array_fetch(out, start), a, gkyl_array_cfetch(inp, start));
  }

  return out;
}

struct gkyl_array*
gkyl_array_set_range_to_range(struct gkyl_array *out, double a,
  const struct gkyl_array *inp, const struct gkyl_range *out_range, const struct gkyl_range *inp_range)
{
  assert(out->elemsz == inp->elemsz);
  assert((inp_range->volume < 1) || (out_range->volume == inp_range->volume));

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_set_range_to_range_cu(out, a, inp, out_range, inp_range); return out; }
#endif

  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n = outnc<inpnc ? outnc : inpnc;

  // Setup linear counter offset for output range/array.
  int iloLocal_out[GKYL_MAX_DIM], iloLocal_inp[GKYL_MAX_DIM];
  for (int d=0; d<out_range->ndim; ++d){
    iloLocal_out[d] = out_range->lower[d];
    iloLocal_inp[d] = inp_range->lower[d];
  }

  int idx_out[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, inp_range);
  while (gkyl_range_iter_next(&iter)) {
    for (int d=0; d<out_range->ndim; ++d)
      idx_out[d] = iloLocal_out[d] + (iter.idx[d] - iloLocal_inp[d]);

    long linidx_inp = gkyl_range_idx(inp_range, iter.idx);
    long linidx_out = gkyl_range_idx(out_range, idx_out);
    array_set1(n,
      gkyl_array_fetch(out, linidx_out), a, gkyl_array_cfetch(inp, linidx_inp));
  }

  return out;
}

struct gkyl_array*
gkyl_array_set_offset_range(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, int coff, const struct gkyl_range *range)
{
  assert(out->type == GKYL_DOUBLE && inp->type == GKYL_DOUBLE);
  assert(out->size == inp->size);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_set_offset_range_cu(out, a, inp, coff, range); return out; }
#endif

  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n;
  int outoff, inoff;
  if (outnc < inpnc) {
    n = outnc;
    outoff = 0;
    inoff = coff;
  } else {
    n = inpnc;
    outoff = coff;
    inoff = 0;
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    double *out_d = gkyl_array_fetch(out, start);
    const double *inp_d = gkyl_array_cfetch(inp, start);
    array_set1(n, out_d+outoff, a, inp_d+inoff);
  }

  return out;
}

struct gkyl_array*
gkyl_array_scale_range(struct gkyl_array *out,
  double a, const struct gkyl_range *range)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_scale_range_cu(out, a, range); return out; }
#endif

  return gkyl_array_set_range(out, a, out, range);
}

struct gkyl_array*
gkyl_array_shiftc_range(struct gkyl_array* out, double a, unsigned k, const struct gkyl_range *range)
{
  assert(out->type == GKYL_DOUBLE);
  assert(k < NCOM(out));
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_shiftc_range_cu(out, a, k, range); return out; }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    double *out_d = gkyl_array_fetch(out, start);
    out_d[k] += a;
  }
  return out;
}

struct gkyl_array*
gkyl_array_comp_op_range(struct gkyl_array *out, enum gkyl_array_op op, double a,
  const struct gkyl_array *in1, double b, const struct gkyl_array *in2,
  const struct gkyl_range *range)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->ncomp == in1->ncomp);
  assert(out->size == in1->size);
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_comp_op_range_cu(out, op, a, in1, b, in2, range); return out; }
#endif

  double *out_d = out->data;
  int nc = out->ncomp;
  long sz = out->size;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  switch (op) {
    case GKYL_MAX:
    case GKYL_MIN:
    case GKYL_SUM:
    case GKYL_ABS_MAX:
    case GKYL_SQ_SUM:
    case GKYL_RMS:
      assert(false);
      break;
    case GKYL_ABS:
      while (gkyl_range_iter_next(&iter)) {
        long linidx = gkyl_range_idx(range, iter.idx);
        double *out_c = gkyl_array_fetch(out, linidx);
        const double *in1_c = gkyl_array_cfetch(in1, linidx);
        for (size_t k=0; k<nc; ++k)
          out_c[k] = fabs(a*in1_c[k]);
      }
      break;                                            
    case GKYL_INV:
      while (gkyl_range_iter_next(&iter)) {
        long linidx = gkyl_range_idx(range, iter.idx);
        double *out_c = gkyl_array_fetch(out, linidx);
        const double *in1_c = gkyl_array_cfetch(in1, linidx);
        for (size_t k=0; k<nc; ++k)
          out_c[k] = a/in1_c[k];
      }
      break;                                            
    case GKYL_PROD:
      assert(out->ncomp == in2->ncomp);
      assert(out->size == in2->size);
      while (gkyl_range_iter_next(&iter)) {
        long linidx = gkyl_range_idx(range, iter.idx);
        double *out_c = gkyl_array_fetch(out, linidx);
        const double *in1_c = gkyl_array_cfetch(in1, linidx);
        const double *in2_c = gkyl_array_cfetch(in2, linidx);
        for (size_t k=0; k<nc; ++k)
          out_c[k] = a*in1_c[k]*in2_c[k]+b;
      }
      break;                                            
    case GKYL_DIV:
      assert(out->ncomp == in2->ncomp);
      assert(out->size == in2->size);
      while (gkyl_range_iter_next(&iter)) {
        long linidx = gkyl_range_idx(range, iter.idx);
        double *out_c = gkyl_array_fetch(out, linidx);
        const double *in1_c = gkyl_array_cfetch(in1, linidx);
        const double *in2_c = gkyl_array_cfetch(in2, linidx);
        for (size_t k=0; k<nc; ++k)
          out_c[k] = a*in1_c[k]/in2_c[k]+b;
      }
      break;                                            
    case GKYL_AXPBY:
      assert(out->ncomp == in2->ncomp);
      assert(out->size == in2->size);
      while (gkyl_range_iter_next(&iter)) {
        long linidx = gkyl_range_idx(range, iter.idx);
        double *out_c = gkyl_array_fetch(out, linidx);
        const double *in1_c = gkyl_array_cfetch(in1, linidx);
        const double *in2_c = gkyl_array_cfetch(in2, linidx);
        for (size_t k=0; k<nc; ++k)
          out_c[k] = a*in1_c[k]+b*in2_c[k];
      }
      break;                                            
  }
  return out;
}

struct gkyl_array*
gkyl_array_copy_range(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  assert(out->size == inp->size && out->elemsz == inp->elemsz);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_copy_range_cu(out, inp, range); return out; }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    memcpy(gkyl_array_fetch(out, start), gkyl_array_cfetch(inp, start), inp->esznc);
  }
  return out;
}

struct gkyl_array*
gkyl_array_copy_range_to_range(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *out_range, const struct gkyl_range *inp_range)
{
  assert(out->elemsz == inp->elemsz);
  assert((inp_range->volume < 1) || (out_range->volume == inp_range->volume));

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_copy_range_to_range_cu(out, inp, out_range, inp_range); return out; }
#endif

  // Setup linear counter offset for output range/array.
  int iloLocal_out[GKYL_MAX_DIM], iloLocal_inp[GKYL_MAX_DIM];
  for (int d=0; d<out_range->ndim; ++d){
    iloLocal_out[d] = out_range->lower[d];
    iloLocal_inp[d] = inp_range->lower[d];
  }

  int idx_out[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, inp_range);
  while (gkyl_range_iter_next(&iter)) {
    for (int d=0; d<out_range->ndim; ++d)
      idx_out[d] = iloLocal_out[d] + (iter.idx[d] - iloLocal_inp[d]);

    long linidx_inp = gkyl_range_idx(inp_range, iter.idx);
    long linidx_out = gkyl_range_idx(out_range, idx_out);
    memcpy(gkyl_array_fetch(out, linidx_out), gkyl_array_cfetch(inp, linidx_inp), inp->esznc);
  }
  return out;
}

void
gkyl_array_copy_to_buffer(void *data, const struct gkyl_array *arr,
  const struct gkyl_range *range)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) { gkyl_array_copy_to_buffer_cu(data, arr, range); return; }
#endif

#define _F(loc) gkyl_array_cfetch(arr, loc)

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    memcpy(((char*) data) + arr->esznc*count++, _F(start), arr->esznc);
  }

#undef _F
}

void
gkyl_array_copy_from_buffer(struct gkyl_array *arr,
  const void *data, const struct gkyl_range *range)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) { gkyl_array_copy_from_buffer_cu(arr, data, range); return; }
#endif

#define _F(loc) gkyl_array_fetch(arr, loc)

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    memcpy(_F(start), ((char*) data) + arr->esznc*count++, arr->esznc);
  }

#undef _F
}

void
gkyl_array_copy_to_buffer_fn(void *data, const struct gkyl_array *arr,
  const struct gkyl_range *range, struct gkyl_array_copy_func *cf)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) { gkyl_array_copy_to_buffer_fn_cu(data, arr, range, cf); return; }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *inp = gkyl_array_cfetch(arr, loc);
    double *out = gkyl_flat_fetch(data, arr->esznc*count);
    cf->func(NCOM(arr), out, inp, cf->ctx);
    count += 1;
  }
}

void
gkyl_array_flip_copy_to_buffer_fn(void *data, const struct gkyl_array *arr,
  int dir, const struct gkyl_range *range, struct gkyl_array_copy_func *cf)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) {
    if (gkyl_array_is_cu_dev(arr)) { gkyl_array_flip_copy_to_buffer_fn_cu(data, arr, dir, range, cf); return; }
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  int fidx[GKYL_MAX_DIM]; // flipped index
  struct gkyl_range buff_range;
  gkyl_range_init(&buff_range, range->ndim, range->lower, range->upper);

  int uplo = range->upper[dir]+range->lower[dir];

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    gkyl_copy_int_arr(range->ndim, iter.idx, fidx);
    fidx[dir] = uplo - iter.idx[dir];
    
    long count = gkyl_range_idx(&buff_range, fidx);

    const double *inp = gkyl_array_cfetch(arr, loc);
    double *out = gkyl_flat_fetch(data, arr->esznc*count);
    cf->func(NCOM(arr), out, inp, cf->ctx);
  }
}

static double
calc_rel_diff(double a, double b)
{
  if (isnan(a) || isnan(b)) return DBL_MAX;
  
  double absa = fabs(a), absb = fabs(b), diff = fabs(a-b);
  if (a == b) return 0;
  if (a == 0 || b == 0 || (absa+absb < DBL_MIN)) return diff;
  return diff/fmin(absa+absb, DBL_MAX);
}

struct gkyl_array_diff
gkyl_array_diff(const struct gkyl_array *arr1, const struct gkyl_array *arr2, const struct gkyl_range *range)
{
  struct gkyl_array_diff incompat = {
    .is_compatible = false,
    .max_abs_diff = DBL_MAX,
    .min_abs_diff = DBL_MAX,
    .max_rel_diff = DBL_MAX,
    .min_rel_diff = DBL_MAX
  };

  if ((arr1->type != GKYL_DOUBLE) && (arr2->type != GKYL_DOUBLE))
    return incompat;

  if (gkyl_array_is_cu_dev(arr1) || gkyl_array_is_cu_dev(arr2))
    return incompat;    

  if (arr1->elemsz != arr2->elemsz)
    return incompat;    

  if (arr1->ncomp != arr2->ncomp)
    return incompat;

  if (arr1->size != arr2->size)
    return incompat;

  double max_abs_diff = -DBL_MAX, max_rel_diff = -DBL_MAX;
  double min_abs_diff = DBL_MAX, min_rel_diff = DBL_MAX;
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    
    long loc = gkyl_range_idx(range, iter.idx);
    const double *a1 = gkyl_array_cfetch(arr1, loc);
    const double *a2 = gkyl_array_cfetch(arr2, loc);

    for (int c=0; c<arr1->ncomp; ++c) {
      max_abs_diff = fmax(max_abs_diff, a1[c]-a2[c]);
      min_abs_diff = fmin(min_abs_diff, a1[c]-a2[c]);
      double rel_diff = calc_rel_diff(a1[c], a2[c]);
      max_rel_diff = fmax(max_rel_diff, rel_diff);
      min_rel_diff = fmin(min_rel_diff, rel_diff);
    }
  }

  return (struct gkyl_array_diff) {
    .is_compatible = true,
    .max_abs_diff = max_abs_diff,
    .min_abs_diff = min_abs_diff,
    .max_rel_diff = max_rel_diff,
    .min_rel_diff = min_rel_diff
  };
}

struct gkyl_array*
gkyl_array_error_denom_fac(struct gkyl_array* out, double eps_rel, double eps_abs, const struct gkyl_array* inp)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->size == inp->size);
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_error_denom_fac_cu(out, eps_rel, eps_abs, inp); return out; }
#endif

  double *out_d = out->data;
  double *inp_d = inp->data;
  for (size_t i=0; i<out->size; ++i) {
    double sqsum = 0.0;
    for (size_t k=0; k<inp->ncomp; ++k)
      sqsum += pow(inp_d[i*inp->ncomp+k],2);

    double error_denom_fac = 1.0/(eps_rel*sqrt(sqsum/inp->ncomp) + eps_abs);

    for (size_t k=0; k<out->ncomp; ++k)
      out_d[i*out->ncomp+k] = error_denom_fac;
  }
  return out;
}

struct gkyl_array*
gkyl_array_error_denom_fac_range(struct gkyl_array* out, double eps_rel, double eps_abs,
  const struct gkyl_array* inp, const struct gkyl_range *range)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->size == inp->size);
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    gkyl_array_error_denom_fac_range_cu(out, eps_rel, eps_abs, inp, range);
    return out;
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    double *out_d = gkyl_array_fetch(out, start);
    const double *inp_d = gkyl_array_cfetch(inp, start);

    double sqsum = 0.0;
    for (size_t k=0; k<inp->ncomp; ++k)
      sqsum += pow(inp_d[k],2);

    double error_denom_fac = pow(1.0/(eps_rel*sqrt(sqsum/inp->ncomp) + eps_abs), 2.0);

    for (size_t k=0; k<out->ncomp; ++k)
      out_d[k] = error_denom_fac;
  }
  return out;
}

struct gkyl_array*
gkyl_array_max_by_cell_per_cell_avg_range(struct gkyl_array* out, const struct gkyl_array* inp, struct gkyl_range *range)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->size == inp->size && out->elemsz == inp->elemsz);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out) && gkyl_array_is_cu_dev(inp)) {
    gkyl_array_max_by_cell_per_cell_avg_range_cu(out, inp, range);
    return out;
  }
#endif

  long nc = NCOM(out);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    double *out_d = gkyl_array_fetch(out, start);
    const double *inp_d = gkyl_array_cfetch(inp, start);
    if (out_d[0] < inp_d[0]) {
      for (int c=0; c<nc; ++c) {
        out_d[c] = inp_d[c];
      }
    }
  }
  return out;
}
