#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_mom_gyrokinetic_priv.h>
#include <gkyl_util.h>

void
gkyl_gk_mom_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_mom_type *base = container_of(ref, struct gkyl_mom_type, ref_count);
  struct mom_type_gyrokinetic *mom_gk = container_of(base, struct mom_type_gyrokinetic, momt);
  gkyl_velocity_map_release(mom_gk->vel_map);
  gkyl_gk_geometry_release(mom_gk->gk_geom);
  if (mom_gk->phi != 0)
    gkyl_array_release(mom_gk->phi);

  if (gkyl_mom_type_is_cu_dev(base)) {
    // free inner on_dev object
    struct mom_type_gyrokinetic *mom_gk_cu = container_of(base->on_dev, struct mom_type_gyrokinetic, momt);
    gkyl_cu_free(mom_gk_cu);
  }

  gkyl_free(mom_gk);
}

struct gkyl_mom_type*
gkyl_mom_gyrokinetic_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, double mass, double charge, const struct gkyl_velocity_map* vel_map,
  const struct gk_geometry *gk_geom, struct gkyl_array *phi, enum gkyl_distribution_moments mom_type, bool use_gpu)
{
  assert(cbasis->poly_order == pbasis->poly_order);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_mom_gyrokinetic_cu_dev_new(cbasis, pbasis, conf_range, mass, charge, vel_map, gk_geom, phi, mom_type);
#endif    

  struct mom_type_gyrokinetic *mom_gk = gkyl_malloc(sizeof(struct mom_type_gyrokinetic));
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_gk->momt.cdim = cdim;
  mom_gk->momt.pdim = pdim;
  mom_gk->momt.poly_order = poly_order;
  mom_gk->momt.num_config = cbasis->num_basis;
  mom_gk->momt.num_phase = pbasis->num_basis;

  // choose kernel tables based on basis-function type
  const gkyl_gyrokinetic_mom_kern_list *m0_kernels, *m1_kernels, *m2_kernels, *m3_kernels, 
    *m2par_kernels, *m2perp_kernels, *m3par_kernels, *m3perp_kernels,
    *m0m1m2_kernels, *m0m1m2m3_kernels,
    *m0m1m2parm2perp_kernels, *m0m1m2parm2perpm3parm3perp_kernels,
    *hamiltonian_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      m0_kernels = ser_m0_kernels;
      m1_kernels = ser_m1_kernels;
      m2_kernels = ser_m2_kernels;
      m3_kernels = ser_m3_kernels;
      m2par_kernels = ser_m2par_kernels;
      m2perp_kernels = ser_m2perp_kernels;
      m3par_kernels = ser_m3par_kernels;
      m3perp_kernels = ser_m3perp_kernels;
      m0m1m2_kernels = ser_m0m1m2_kernels;
      m0m1m2m3_kernels = ser_m0m1m2m3_kernels;
      m0m1m2parm2perp_kernels = ser_m0m1m2parm2perp_kernels;
      m0m1m2parm2perpm3parm3perp_kernels = ser_m0m1m2parm2perpm3parm3perp_kernels;
      hamiltonian_kernels = ser_hamiltonian_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  if (mom_type == GKYL_F_MOMENT_M0) { // Density.
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m0_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = m0_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M1) { // Parallel momentum.
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m1_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = m1_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M2) { // Total kinetic energy.
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = m2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M3) { // Total heat flux.
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m3_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = m3_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M2PAR) { // Parallel energy.
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2par_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = m2par_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M2PERP) { // Perpendicular energy.
    assert(vdim == 2);
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = m2perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M3PAR) { // Parallel heat flux.
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m3par_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = m3par_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M3PERP) { // Perpendicular heat flux.
    assert(vdim == 2);
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m3perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = m3perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M0M1M2) { 
    // Density, parallel momentum, and total energy computed together.
    assert(cv_index[cdim].vdim[vdim] != -1);   
    assert(NULL != m0m1m2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = m0m1m2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 3;
  }
  else if (mom_type == GKYL_F_MOMENT_M0M1M2M3) { 
    // Density, parallel momentum, total energy and total heat flux computed together.
    assert(cv_index[cdim].vdim[vdim] != -1);   
    assert(NULL != m0m1m2m3_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = m0m1m2m3_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 4;
  }
  else if (mom_type == GKYL_F_MOMENT_M0M1M2PARM2PERP) {
    // Density, parallel momentum, parallel and perpendicular
    // kinetic energy computed together.
    assert(vdim == 2);
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m0m1m2parm2perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = m0m1m2parm2perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = vdim+2;
  }
  else if (mom_type == GKYL_F_MOMENT_M0M1M2PARM2PERPM3PARM3PERP) {
    // Density, parallel momentum, parallel and perpendicular
    // kinetic energy, and parallel and perpendicular heat flux computed together.
    assert(vdim == 2);
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m0m1m2parm2perpm3parm3perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = m0m1m2parm2perpm3parm3perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 2*vdim+2;
  }
  else if (mom_type == GKYL_F_MOMENT_HAMILTONIAN) {
    // M0, mass*M0 and total particle energy
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != hamiltonian_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = hamiltonian_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 3;
  }
  else {
    // string not recognized
    printf("Error: requested moment %d.\n", mom_type);
    gkyl_exit("gkyl_mom_type_gyrokinetic: Unrecognized moment requested!");
  }

  mom_gk->conf_range = *conf_range;
  mom_gk->mass = mass;
  mom_gk->charge = charge;
  mom_gk->vel_map = gkyl_velocity_map_acquire(vel_map);
  mom_gk->gk_geom = gkyl_gk_geometry_acquire(gk_geom);
  mom_gk->phi = 0;
  if (phi)
    mom_gk->phi = gkyl_array_acquire(phi);
  
  mom_gk->momt.flags = 0;
  GKYL_CLEAR_CU_ALLOC(mom_gk->momt.flags);
  mom_gk->momt.ref_count = gkyl_ref_count_init(gkyl_gk_mom_free);
  
  mom_gk->momt.on_dev = &mom_gk->momt; // on host, self-reference
    
  return &mom_gk->momt;
}

struct gkyl_mom_type*
gkyl_int_mom_gyrokinetic_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, double mass, double charge, const struct gkyl_velocity_map* vel_map,
  const struct gk_geometry *gk_geom, struct gkyl_array *phi, enum gkyl_distribution_moments mom_type, bool use_gpu)
{
  assert(cbasis->poly_order == pbasis->poly_order);

#ifdef GKYL_HAVE_CUDA
  if(use_gpu)
    return gkyl_int_mom_gyrokinetic_cu_dev_new(cbasis, pbasis, conf_range, mass, charge, vel_map, gk_geom, phi, mom_type);
#endif

  struct mom_type_gyrokinetic *mom_gk = gkyl_malloc(sizeof(struct mom_type_gyrokinetic));
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_gk->momt.cdim = cdim;
  mom_gk->momt.pdim = pdim;
  mom_gk->momt.poly_order = poly_order;
  mom_gk->momt.num_config = cbasis->num_basis;
  mom_gk->momt.num_phase = pbasis->num_basis;
  
  // Choose kernel tables based on basis-function type.
  const gkyl_gyrokinetic_mom_kern_list *int_m0_kernels, *int_m1_kernels, *int_m2_kernels, *int_m3_kernels, 
    *int_m2par_kernels, *int_m2perp_kernels, *int_m3par_kernels, *int_m3perp_kernels,
    *int_m0m1m2_kernels, *int_m0m1m2m3_kernels,
    *int_m0m1m2parm2perp_kernels, *int_m0m1m2parm2perpm3parm3perp_kernels,
    *int_hamiltonian_kernels;

  // Set kernel pointer.
  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      int_m0_kernels = ser_int_m0_kernels;
      int_m1_kernels = ser_int_m1_kernels;
      int_m2_kernels = ser_int_m2_kernels;
      int_m3_kernels = ser_int_m3_kernels;
      int_m2par_kernels = ser_int_m2par_kernels;
      int_m2perp_kernels = ser_int_m2perp_kernels;
      int_m3par_kernels = ser_int_m3par_kernels;
      int_m3perp_kernels = ser_int_m3perp_kernels;
      int_m0m1m2_kernels = ser_int_m0m1m2_kernels;
      int_m0m1m2m3_kernels = ser_int_m0m1m2m3_kernels;
      int_m0m1m2parm2perp_kernels = ser_int_m0m1m2parm2perp_kernels;
      int_m0m1m2parm2perpm3parm3perp_kernels = ser_int_m0m1m2parm2perpm3parm3perp_kernels;
      int_hamiltonian_kernels = ser_int_hamiltonian_kernels;
      break;

    default:
      assert(false);
      break;    
  }  

  assert(cv_index[cdim].vdim[vdim] != -1);   

  if (mom_type == GKYL_F_MOMENT_M0) {
    // Density moment only.
    assert(NULL != int_m0_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = int_m0_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M1) {
    // Parallel momentum moment only.
    assert(NULL != int_m1_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = int_m1_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M2PAR) {
    // Parallel kinetic energy moment only.
    assert(NULL != int_m2par_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = int_m2par_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M2PERP) {
    // Perpendicular kinetic energy moment only.
    assert(vdim == 2);
    assert(NULL != int_m2perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = int_m2perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M2) {
    // Kinetic energy moment only.
    assert(NULL != int_m2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = int_m2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M3) {
    // Heat flux moment only.
    assert(NULL != int_m3_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = int_m3_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M3PAR) {
    // Parallel heat flux moment only.
    assert(NULL != int_m3par_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = int_m3par_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M3PERP) {
    // Perpendicular heat flux moment only.
    assert(vdim == 2);
    assert(NULL != int_m3perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = int_m3perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M0M1M2) {
    // Density, parallel momentum, and kinetic energy computed together.
    assert(NULL != int_m0m1m2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = int_m0m1m2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 3;
  }
  else if (mom_type == GKYL_F_MOMENT_M0M1M2) {
    // Density, parallel momentum, kinetic energy and heat flux computed together.
    assert(NULL != int_m0m1m2m3_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = int_m0m1m2m3_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 4;
  }
  else if (mom_type == GKYL_F_MOMENT_M0M1M2PARM2PERP) {
    // Density, parallel momentum, and parallel and perpendicular kinetic energy.
    assert(vdim == 2);
    assert(NULL != int_m0m1m2parm2perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = int_m0m1m2parm2perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = vdim+2;
  }
  else if (mom_type == GKYL_F_MOMENT_M0M1M2PARM2PERPM3PARM3PERP) {
    // Density, parallel momentum, parallel and perpendicular kinetic energy,
    // and parallel and perpendicular heat flux.
    assert(vdim == 2);
    assert(NULL != int_m0m1m2parm2perpm3parm3perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = int_m0m1m2parm2perpm3parm3perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 2*vdim+2;
  }
  else if (mom_type == GKYL_F_MOMENT_HAMILTONIAN) {
    // Density), parallel momentum, and total energy computed together.
    assert(NULL != int_hamiltonian_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gk->momt.kernel = int_hamiltonian_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gk->momt.num_mom = 3;
  }
  else {
    fprintf(stderr,"Moment option %d not available.\n",mom_type);
    assert(false);
  }

  mom_gk->conf_range = *conf_range;
  mom_gk->mass = mass;
  mom_gk->charge = charge;
  mom_gk->vel_map = gkyl_velocity_map_acquire(vel_map);
  mom_gk->gk_geom = gkyl_gk_geometry_acquire(gk_geom);
  mom_gk->phi = 0;
  if (phi)
    mom_gk->phi = gkyl_array_acquire(phi);

  mom_gk->momt.flags = 0;
  GKYL_CLEAR_CU_ALLOC(mom_gk->momt.flags);
  mom_gk->momt.ref_count = gkyl_ref_count_init(gkyl_gk_mom_free);
  
  mom_gk->momt.on_dev = &mom_gk->momt; // on host, self-reference
    
  return &mom_gk->momt;  
}

