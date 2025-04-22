#pragma once

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_mom_vlasov_sr_auxfields { 
  const struct gkyl_array *gamma; // gamma = sqrt(1 + p^2)
  const struct gkyl_array *vmap; // velocity space mapping for mapped velocity grids
  const struct gkyl_array *jacob_vel_inv; // inverse velocity space Jacobian for mapped velocity grids
};

/**
 * Create new special relativistic Vlasov moment type object. 
 * Valid 'mom' strings are "M0", "M1i", "M2", "M3i"
 * Note: M2 is the integral(gamma*f) velocity moment and M3i the integral(p*f) velocity moment in relativity
 * "Ni" = (M0, M1i) (1+vdim components)
 * "Tij" = stress-energy tensor (M2, M3i (vdim components), Stress tensor (vdim*(vdim+1))/2 components))
 * Also permits M0_upper and M0_lower, which integrate the zeroth moment over a sub-domain of the upper
 * or lower half-plane respectively with integration limits set by the v_thresh input parameter. 
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param vel_range Velocity space range
 * @param use_vmap bool to determine if we are using mapped velocity grid kernels
 * @param v_thresh Threshold velocity for integrations over only the upper or lower half-plane. 
 * @param mom Name of moment to compute
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_mom_type* gkyl_mom_vlasov_sr_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range, 
  bool use_vmap, double v_thresh, const char *mom, bool use_gpu);

/**
 * Create new special relativistic Vlasov moment type object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_mom_type* gkyl_mom_vlasov_sr_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range, 
  bool use_vmap, double v_thresh, const char *mom);

/**
 * Create new special relativistic Vlasov integrated moment type
 * object. Integrates (M0, M2, M3i) (2 + vdim) components.
 * Note the different order from non-relativistic since we compute the
 * integrated energy flux instead of the integrated mass flux. 
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param vel_range Velocity space range
 * @param use_vmap bool to determine if we are using mapped velocity grid kernels
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_mom_type* gkyl_int_mom_vlasov_sr_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range, 
  bool use_vmap, bool use_gpu);

/**
 * Create new special relativistic Vlasov integrated moment type
 * object on NV-GPU: see new() method above for documentation.
 */
struct gkyl_mom_type* gkyl_int_mom_vlasov_sr_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range, 
  bool use_vmap);

/**
 * Set the auxiliary fields needed in computing moments.
 * 
 * @param momt moment type.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_mom_vlasov_sr_set_auxfields(const struct gkyl_mom_type *momt,
  struct gkyl_mom_vlasov_sr_auxfields auxin);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set auxiliary fields needed in computing moments.
 * 
 * @param momt moment type.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_mom_vlasov_sr_set_auxfields_cu(const struct gkyl_mom_type *momt,
  struct gkyl_mom_vlasov_sr_auxfields auxin);

#endif
