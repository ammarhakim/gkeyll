/* -*- c++ -*- */

extern "C" {
#include <gkyl_mom_weighted_gyrokinetic.h>
#include <gkyl_mom_weighted_gyrokinetic_priv.h>
#include <gkyl_array_ops.h>
#include <float.h>
}

void
mom_weighted_choose_kernel_cu(struct gkyl_mom_weighted_gyrokinetic_kernels *kernels,
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, enum gkyl_distribution_moments mom_type,
  enum gkyl_mom_weight_type wgt_type, bool is_integrated)
{
}

