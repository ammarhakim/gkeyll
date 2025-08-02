/* -*- c++ -*- */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_mom_vlasov_priv.h>
#include <gkyl_util.h>
}

static int
v_num_mom(int vdim, enum gkyl_distribution_moments mom_type)
{
  int m3ijk_count[] = { 1, 4, 10 };
  int num_mom = 0;
  
  switch (mom_type) {
    case GKYL_F_MOMENT_M0:
    case GKYL_F_MOMENT_ENERGY:
    case GKYL_F_MOMENT_M2:
    case GKYL_F_MOMENT_M0_UPPER:
    case GKYL_F_MOMENT_M0_LOWER:
      num_mom = 1;
      break;

    case GKYL_F_MOMENT_M1:
    case GKYL_F_MOMENT_M1_FROM_H:
    case GKYL_F_MOMENT_M3:
      num_mom = vdim;
      break;

    case GKYL_F_MOMENT_M2IJ:
      num_mom = vdim*(vdim+1)/2;
      break;

    case GKYL_F_MOMENT_M3IJK:
      num_mom = m3ijk_count[vdim-1];
      break;

    case GKYL_F_MOMENT_M0M1M2:
      num_mom = vdim+2;
      break;      
      
    default: // Can't happen.
      fprintf(stderr,"Moment option %d not available.\n",mom_type);
      assert(false);
      break;
  }

  return num_mom;
}

__global__
static void
set_cu_ptrs(struct mom_type_vlasov* mom_vlasov, enum gkyl_distribution_moments mom_type,
  enum gkyl_basis_type b_type, int cdim, int vdim, int poly_order, 
  enum gkyl_model_id model_id, const struct gkyl_array *hamil)
{
  int m3ijk_count[] = { 1, 4, 10 };

  mom_vlasov->momt.kernel = kernel;
  
  // choose kernel tables based on basis-function type
  const gkyl_mom_kern_list *m0_kernels, *m1i_hamil_vel_kernels, *m1i_hamil_gen_kernels,
    *m2_hamil_vel_kernels, *m2_hamil_gen_kernels, *m3i_hamil_vel_kernels, *m3i_hamil_gen_kernels, 
    *m2ij_kernels, *m3ijk_kernels, *five_moments_hamil_vel_kernels, *five_moments_hamil_gen_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      m0_kernels = ser_m0_kernels;
      m1i_hamil_vel_kernels = ser_m1i_hamil_vel_kernels;
      m2_hamil_vel_kernels = ser_m2_hamil_vel_kernels;
      // m3i_hamil_vel_kernels = ser_m3i_hamil_vel_kernels;
      m1i_hamil_gen_kernels = ser_m1i_hamil_gen_kernels;
      m2_hamil_gen_kernels = ser_m2_hamil_gen_kernels;
      // m3i_hamil_gen_kernels = ser_m3i_hamil_gen_kernels;
      m2ij_kernels = ser_m2ij_kernels;
      m3ijk_kernels = ser_m3ijk_kernels;
      five_moments_hamil_vel_kernels = ser_five_moments_hamil_vel_kernels;
      five_moments_hamil_gen_kernels = ser_five_moments_hamil_gen_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      m0_kernels = tensor_m0_kernels;
      m1i_hamil_vel_kernels = tensor_m1i_hamil_vel_kernels;
      m2_hamil_vel_kernels = tensor_m2_hamil_vel_kernels;
      // m3i_hamil_vel_kernels = tensor_m3i_hamil_vel_kernels;
      m1i_hamil_gen_kernels = tensor_m1i_hamil_gen_kernels;
      m2_hamil_gen_kernels = tensor_m2_hamil_gen_kernels;
      // m3i_hamil_gen_kernels = tensor_m3i_hamil_gen_kernels;
      m2ij_kernels = tensor_m2ij_kernels;
      m3ijk_kernels = tensor_m3ijk_kernels;
      five_moments_hamil_vel_kernels = tensor_five_moments_hamil_vel_kernels;
      five_moments_hamil_gen_kernels = tensor_five_moments_hamil_gen_kernels;
      break;

    default:
      assert(false);
      break;    
  }
  
  switch (mom_type) {
    case GKYL_F_MOMENT_M0:
      mom_vlasov->kernel = m0_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      mom_vlasov->momt.num_mom = 1;
      break;

    case GKYL_F_MOMENT_M1:
    case GKYL_F_MOMENT_M1_FROM_H:
      if (model_id == GKYL_MODEL_DEFAULT || model_id == GKYL_MODEL_SR) {
        mom_vlasov->kernel = m1i_hamil_vel_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      }
      else {
        mom_vlasov->kernel = m1i_hamil_gen_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      }
      mom_vlasov->momt.num_mom = vdim;
      break;

    case GKYL_F_MOMENT_M2:
    case GKYL_F_MOMENT_ENERGY:
      if (model_id == GKYL_MODEL_DEFAULT || model_id == GKYL_MODEL_SR) {
        mom_vlasov->kernel = m2_hamil_vel_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      }
      else {
        mom_vlasov->kernel = m2_hamil_gen_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      }
      mom_vlasov->momt.num_mom = 1;
      break;

    case GKYL_F_MOMENT_M3:
      // if (model_id == GKYL_MODEL_DEFAULT || model_id == GKYL_MODEL_SR) {
      //   mom_vlasov->kernel = m3i_hamil_vel_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      // }
      // else {
      //   mom_vlasov->kernel = m3i_hamil_gen_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      // }
      mom_vlasov->momt.num_mom = vdim;
      break;

    case GKYL_F_MOMENT_M2IJ:
      mom_vlasov->kernel = m2ij_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      mom_vlasov->momt.num_mom = vdim*(vdim+1)/2;
      break;

    case GKYL_F_MOMENT_M3IJK:
      mom_vlasov->kernel = m3ijk_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      mom_vlasov->momt.num_mom = m3ijk_count[vdim-1];
      break;

    case GKYL_F_MOMENT_M0M1M2:
      if (model_id == GKYL_MODEL_DEFAULT || model_id == GKYL_MODEL_SR) {
        mom_vlasov->kernel = five_moments_hamil_vel_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      }
      else {
        mom_vlasov->kernel = five_moments_hamil_gen_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      }
      mom_vlasov->momt.num_mom = vdim+2;
      break;

    case GKYL_F_MOMENT_M0_UPPER:
      mom_vlasov->momt.kernel = tensor_m0_upper_kernels[cdim-1].kernels[poly_order];
      mom_vlasov->momt.num_mom = 1;
      break;

    case GKYL_F_MOMENT_M0_LOWER:
      mom_vlasov->momt.kernel = tensor_m0_lower_kernels[cdim-1].kernels[poly_order];
      mom_vlasov->momt.num_mom = 1;
      break;   

    default: // can't happen
      break;
  }
}

struct gkyl_mom_type*
gkyl_mom_vlasov_cu_dev_inew(const struct gkyl_mom_vlasov_inp *inp)
{
  assert(inp->conf_basis->poly_order == inp->phase_basis->poly_order);

  struct mom_type_vlasov *mom_vlasov = (struct mom_type_vlasov*) gkyl_malloc(sizeof(*mom_vlasov));
  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim, vdim = pdim-cdim;
  int poly_order = inp->conf_basis->poly_order;

  mom_vlasov->momt.cdim = cdim;
  mom_vlasov->momt.pdim = pdim;
  mom_vlasov->momt.poly_order = poly_order;
  mom_vlasov->momt.num_config = inp->conf_basis->num_basis;
  mom_vlasov->momt.num_phase = inp->phase_basis->num_basis;

  // Determine Hamiltonian dimensionality and index offset for indexing Hamiltonian
  // from an input phase space index. 
  if (inp->model_id == GKYL_MODEL_DEFAULT || inp->model_id == GKYL_MODEL_SR) {
    mom_vlasov->hamil_dim = vdim; 
    mom_vlasov->hamil_offset = cdim; 
  }
  else {
    mom_vlasov->hamil_dim = pdim; 
    mom_vlasov->hamil_offset = 0; 
  }
  mom_vlasov->hamil_range = *inp->hamil_range;
  struct gkyl_array *hamil_ho = gkyl_array_acquire(inp->hamil); 
  mom_vlasov->hamil = hamil_ho->on_dev; // store pointer to on_dev for copying over to device. 

  mom_vlasov->vel_range = *inp->vel_range;
  struct gkyl_array *vmap_ho = 0; 
  struct gkyl_array *jacob_vel_ho = 0; 
  if (inp->use_vmap) {
    vmap_ho = gkyl_array_acquire(inp->vmap); 
    jacob_vel_ho = gkyl_array_acquire(inp->jacob_vel); 
    mom_vlasov->vmap = vmap_ho->on_dev;
    mom_vlasov->jacob_vel = jacob_vel_ho->on_dev; 
  }

  // Threshold velocity for integration of moments over a subset of the domain. 
  mom_vlasov->v_thresh = inp->v_thresh; 

  mom_vlasov->momt.num_mom = v_num_mom(vdim, inp->mom_type); // Number of moments.

  mom_vlasov->momt.flags = 0;
  GKYL_SET_CU_ALLOC(mom_vlasov->momt.flags);
  mom_vlasov->momt.ref_count = gkyl_ref_count_init(gkyl_mom_free);
  
  // copy struct to device
  struct mom_type_vlasov *mom_vlasov_cu = (struct mom_type_vlasov*) gkyl_cu_malloc(sizeof(*mom_vlasov_cu));
  gkyl_cu_memcpy(mom_vlasov_cu, mom_vlasov, sizeof(struct mom_type_vlasov), GKYL_CU_MEMCPY_H2D);

  set_cu_ptrs<<<1,1>>>(mom_vlasov_cu, inp->mom_type, inp->conf_basis->b_type, 
    cdim, vdim, poly_order, inp->model_id, inp->hamil->on_dev);

  mom_vlasov->momt.on_dev = &mom_vlasov_cu->momt;

  // Host-side moment type object should store host pointers.
  mom_vlasov->vmap = vmap_ho; 
  mom_vlasov->jacob_vel = jacob_vel_ho; 
  mom_vlasov->hamil = hamil_ho; 
  
  return &mom_vlasov->momt;
}

__global__
static void
set_int_cu_ptrs(struct mom_type_vlasov* mom_vlasov, enum gkyl_distribution_moments mom_type,
  enum gkyl_basis_type b_type, int cdim, int vdim, int poly_order, 
  enum gkyl_model_id model_id, const struct gkyl_array *hamil)
{
  mom_vlasov->momt.kernel = kernel;

  // Choose kernel tables based on basis-function type.
  const gkyl_mom_kern_list *int_five_moments_hamil_vel_kernels, *int_five_moments_hamil_gen_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      int_five_moments_hamil_vel_kernels = ser_int_five_moments_hamil_vel_kernels;
      int_five_moments_hamil_gen_kernels = ser_int_five_moments_hamil_gen_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      int_five_moments_hamil_vel_kernels = tensor_int_five_moments_hamil_vel_kernels;
      int_five_moments_hamil_gen_kernels = tensor_int_five_moments_hamil_gen_kernels;
      break;

    default:
      assert(false);
      break;    
  }  

  switch (mom_type) {
    case GKYL_F_MOMENT_M0M1M2:
      if (model_id == GKYL_MODEL_DEFAULT || model_id == GKYL_MODEL_SR) {
        mom_vlasov->kernel = int_five_moments_hamil_vel_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      }
      else {
        mom_vlasov->kernel = int_five_moments_hamil_gen_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      } 
      mom_vlasov->momt.num_mom = 2+vdim;
      break;

    default:
      assert(false);
      break;
  }
}

struct gkyl_mom_type*
gkyl_int_mom_vlasov_cu_dev_inew(const struct gkyl_mom_vlasov_inp *inp)
{
  assert(inp->conf_basis->poly_order == inp->phase_basis->poly_order);

  struct mom_type_vlasov *mom_vlasov = (struct mom_type_vlasov*) gkyl_malloc(sizeof(*mom_vlasov));
  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim, vdim = pdim-cdim;
  int poly_order = inp->conf_basis->poly_order;

  mom_vlasov->momt.cdim = cdim;
  mom_vlasov->momt.pdim = pdim;
  mom_vlasov->momt.poly_order = poly_order;
  mom_vlasov->momt.num_config = inp->conf_basis->num_basis;
  mom_vlasov->momt.num_phase = inp->phase_basis->num_basis;
  // Determine Hamiltonian dimensionality and index offset for indexing Hamiltonian
  // from an input phase space index. 
  if (inp->model_id == GKYL_MODEL_DEFAULT || inp->model_id == GKYL_MODEL_SR) {
    mom_vlasov->hamil_dim = vdim; 
    mom_vlasov->hamil_offset = cdim; 
  }
  else {
    mom_vlasov->hamil_dim = pdim; 
    mom_vlasov->hamil_offset = 0; 
  }
  mom_vlasov->hamil_range = *inp->hamil_range;
  struct gkyl_array *hamil_ho = gkyl_array_acquire(inp->hamil); 
  mom_vlasov->hamil = hamil_ho->on_dev; // store pointer to on_dev for copying over to device. 

  mom_vlasov->momt.num_mom = v_num_mom(vdim, inp->mom_type); // Number of moments.

  mom_vlasov->momt.flags = 0;
  GKYL_SET_CU_ALLOC(mom_vlasov->momt.flags);
  mom_vlasov->momt.ref_count = gkyl_ref_count_init(gkyl_mom_free);
  
  // copy struct to device
  struct mom_type_vlasov *mom_vlasov_cu = (struct mom_type_vlasov*) gkyl_cu_malloc(sizeof(*mom_vlasov_cu));
  gkyl_cu_memcpy(mom_vlasov_cu, mom_vlasov, sizeof(struct mom_type_vlasov), GKYL_CU_MEMCPY_H2D);

  set_int_cu_ptrs<<<1,1>>>(mom_vlasov_cu, inp->mom_type, inp->conf_basis->b_type, 
    cdim, vdim, poly_order, inp->model_id, inp->hamil->on_dev);

  mom_vlasov->momt.on_dev = &mom_vlasov_cu->momt;

  // Host-side moment type object should store host pointers.
  mom_vlasov->hamil = hamil_ho;   

  return &mom_vlasov->momt;
}
