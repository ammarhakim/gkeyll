#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_mom_vlasov_priv.h>
#include <gkyl_util.h>

void
gkyl_mom_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_mom_type *base = container_of(ref, struct gkyl_mom_type, ref_count);

  if (gkyl_mom_type_is_cu_dev(base)) {
    // free inner on_dev object
    struct mom_type_vlasov *mom_vlasov = container_of(base->on_dev, struct mom_type_vlasov, momt);
    gkyl_cu_free(mom_vlasov);
  }

  struct mom_type_vlasov *mom_vlasov = container_of(base, struct mom_type_vlasov, momt);
  gkyl_array_release(mom_vlasov->hamil);
  gkyl_free(mom_vlasov);
}

struct gkyl_mom_type*
gkyl_mom_vlasov_inew(const struct gkyl_mom_vlasov_inp *inp)
{
  assert(inp->conf_basis->poly_order == inp->phase_basis->poly_order);

#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_mom_vlasov_cu_dev_inew(inp);
  } 
#endif
  struct mom_type_vlasov *mom_vlasov = gkyl_malloc(sizeof(*mom_vlasov));
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
  mom_vlasov->hamil = gkyl_array_acquire(inp->hamil); 

  mom_vlasov->vel_range = *inp->vel_range;
  mom_vlasov->vmap = 0;
  mom_vlasov->jacob_vel = 0;
  if (inp->use_vmap) {
    mom_vlasov->vmap = gkyl_array_acquire(inp->vmap); 
    mom_vlasov->jacob_vel = gkyl_array_acquire(inp->jacob_vel); 
  }

  // Threshold velocity for integration of moments over a subset of the domain. 
  mom_vlasov->v_thresh = inp->v_thresh; 

  // choose kernel tables based on basis-function type
  const gkyl_vlasov_mom_kern_list *m0_kernels, *m1i_hamil_vel_kernels, *m1i_hamil_gen_kernels,
    *m2_hamil_vel_kernels, *m2_hamil_gen_kernels, 
    *m2ij_kernels, *m3ijk_kernels, *five_moments_hamil_vel_kernels, *five_moments_hamil_gen_kernels;

  switch (inp->conf_basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      m0_kernels = ser_m0_kernels;
      m2ij_kernels = ser_m2ij_kernels;
      m3ijk_kernels = ser_m3ijk_kernels;
      m1i_hamil_vel_kernels = ser_hamil_vel_m1i_kernels;
      m2_hamil_vel_kernels = ser_hamil_vel_m2_kernels;
      five_moments_hamil_vel_kernels = ser_hamil_vel_five_moments_kernels;
      m1i_hamil_gen_kernels = ser_hamil_gen_m1i_kernels;
      m2_hamil_gen_kernels = ser_hamil_gen_m2_kernels;
      five_moments_hamil_gen_kernels = ser_hamil_gen_five_moments_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      m0_kernels = tensor_m0_kernels;
      m2ij_kernels = tensor_m2ij_kernels;
      m3ijk_kernels = tensor_m3ijk_kernels;
      m1i_hamil_vel_kernels = tensor_hamil_vel_m1i_kernels;
      m2_hamil_vel_kernels = tensor_hamil_vel_m2_kernels;
      five_moments_hamil_vel_kernels = tensor_hamil_vel_five_moments_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  if (inp->mom_type == GKYL_F_MOMENT_M0) { // density
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m0_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_vlasov->momt.kernel = m0_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_vlasov->momt.num_mom = 1;
  }
  else if (inp->mom_type == GKYL_F_MOMENT_M1 || inp->mom_type == GKYL_F_MOMENT_M1_FROM_H) { 
    // As part of Hamiltonian Vlasov refactor, assume user wants dH/dv moment when they
    // request M1 in some form. JJ 07/25/25
    assert(cv_index[cdim].vdim[vdim] != -1);    
    if (inp->model_id == GKYL_MODEL_DEFAULT || inp->model_id == GKYL_MODEL_SR) {
      assert(NULL != m1i_hamil_vel_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
      mom_vlasov->momt.kernel = m1i_hamil_vel_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    }
    else {
      if (inp->conf_basis->b_type == GKYL_BASIS_MODAL_TENSOR) {
        gkyl_exit("mom_vlasov M1i hamil_gen: Not currently supported in tensor basis.");
      }
      assert(NULL != m1i_hamil_gen_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
      mom_vlasov->momt.kernel = m1i_hamil_gen_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    }
    mom_vlasov->momt.num_mom = vdim;
  }
  else if (inp->mom_type == GKYL_F_MOMENT_M2 || inp->mom_type == GKYL_F_MOMENT_ENERGY) { 
    // As part of Hamiltonian Vlasov refactor, assume user wants H moment when they
    // request M2 in some form. JJ 07/25/25
    assert(cv_index[cdim].vdim[vdim] != -1);
    if (inp->model_id == GKYL_MODEL_DEFAULT || inp->model_id == GKYL_MODEL_SR) {
      assert(NULL != m2_hamil_vel_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
      mom_vlasov->momt.kernel = m2_hamil_vel_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    }
    else {
      if (inp->conf_basis->b_type == GKYL_BASIS_MODAL_TENSOR) {
        gkyl_exit("mom_vlasov M2 hamil_gen: Not currently supported in tensor basis.");
      }
      assert(NULL != m2_hamil_gen_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
      mom_vlasov->momt.kernel = m2_hamil_gen_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    } 
    mom_vlasov->momt.num_mom = 1;
  }
  else if (inp->mom_type == GKYL_F_MOMENT_M3 || inp->mom_type == GKYL_F_MOMENT_ENERGY_FLUX) { 
    // As part of Hamiltonian Vlasov refactor, assume user wants H*dH/dv moment when they
    // request M3 in some form. JJ 07/25/25
    gkyl_exit("mom_vlasov M3: Not currently supported in new Hamiltonian formulation.");
    // assert(cv_index[cdim].vdim[vdim] != -1);
    // if (inp->model_id == GKYL_MODEL_DEFAULT || inp->model_id == GKYL_MODEL_SR) {
    //   assert(NULL != m3i_hamil_vel_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    //   mom_vlasov->momt.kernel = m3i_hamil_vel_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    // }
    // else {
    //   assert(NULL != m3i_hamil_gen_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    //   mom_vlasov->momt.kernel = m3i_hamil_gen_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    // }    
    // mom_vlasov->momt.num_mom = vdim;
  }
  else if (inp->mom_type == GKYL_F_MOMENT_M2IJ) { // pressure tensor in lab-frame
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2ij_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_vlasov->momt.kernel = m2ij_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_vlasov->momt.num_mom = vdim*(vdim+1)/2;
  }
  else if (inp->mom_type == GKYL_F_MOMENT_M3IJK) { // heat-flux tensor in lab-frame
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m3ijk_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_vlasov->momt.kernel = m3ijk_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];

    int m3ijk_count[] = { 1, 4, 10 };
    mom_vlasov->momt.num_mom = m3ijk_count[vdim-1];
  }
  else if (inp->mom_type == GKYL_F_MOMENT_M0M1M2) { // Zeroth, First, and Second moment computed together
    // As part of Hamiltonian Vlasov refactor, assume user wants {1, dH/dv, H} moments when they
    // request five_moments/M0M1M2. JJ 07/25/25
    assert(cv_index[cdim].vdim[vdim] != -1);
    if (inp->model_id == GKYL_MODEL_DEFAULT || inp->model_id == GKYL_MODEL_SR) {
      assert(NULL != five_moments_hamil_vel_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
      mom_vlasov->momt.kernel = five_moments_hamil_vel_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    }
    else {
      if (inp->conf_basis->b_type == GKYL_BASIS_MODAL_TENSOR) {
        gkyl_exit("mom_vlasov five moments hamil_gen: Not currently supported in tensor basis.");
      }
      assert(NULL != five_moments_hamil_gen_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
      mom_vlasov->momt.kernel = five_moments_hamil_gen_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    } 
    mom_vlasov->momt.num_mom = 2+vdim;
  }
  else if (inp->mom_type == GKYL_F_MOMENT_M0_UPPER) {
    if (inp->conf_basis->b_type != GKYL_BASIS_MODAL_TENSOR && vdim !=1) {
      gkyl_exit("mom_vlasov: M0 upper only defined for tensor basis and vdim = 1!");
    }
    mom_vlasov->momt.kernel = tensor_m0_upper_kernels[cdim-1].kernels[poly_order];
    mom_vlasov->momt.num_mom = 1;
  }
  else if (inp->mom_type == GKYL_F_MOMENT_M0_LOWER) {
    if (inp->conf_basis->b_type != GKYL_BASIS_MODAL_TENSOR && vdim !=1) {
      gkyl_exit("mom_vlasov: M0 lower only defined for tensor basis and vdim = 1!");
    }
    mom_vlasov->momt.kernel = tensor_m0_lower_kernels[cdim-1].kernels[poly_order];
    mom_vlasov->momt.num_mom = 1;
  }
  else {
    // string not recognized
    gkyl_exit("mom_vlasov: Unrecognized moment requested!");
  }

  mom_vlasov->momt.flags = 0;
  GKYL_CLEAR_CU_ALLOC(mom_vlasov->momt.flags);
  mom_vlasov->momt.ref_count = gkyl_ref_count_init(gkyl_mom_free);
  
  mom_vlasov->momt.on_dev = &mom_vlasov->momt; // on host, self-reference
    
  return &mom_vlasov->momt;
}

struct gkyl_mom_type*
gkyl_int_mom_vlasov_inew(const struct gkyl_mom_vlasov_inp *inp)
{
  assert(inp->conf_basis->poly_order == inp->phase_basis->poly_order);

#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_int_mom_vlasov_cu_dev_inew(inp);
  } 
#endif
  struct mom_type_vlasov *mom_vlasov = gkyl_malloc(sizeof(*mom_vlasov));
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
  mom_vlasov->hamil = gkyl_array_acquire(inp->hamil); 

  mom_vlasov->vel_range = *inp->vel_range;
  mom_vlasov->vmap = 0;
  mom_vlasov->jacob_vel = 0;
  if (inp->use_vmap) {
    mom_vlasov->vmap = gkyl_array_acquire(inp->vmap); 
    mom_vlasov->jacob_vel = gkyl_array_acquire(inp->jacob_vel); 
  }

  // Choose kernel tables based on basis-function type.
  const gkyl_vlasov_mom_kern_list *int_five_moments_hamil_vel_kernels, *int_five_moments_hamil_gen_kernels;

  switch (inp->conf_basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      int_five_moments_hamil_vel_kernels = ser_hamil_vel_int_five_moments_kernels;
      int_five_moments_hamil_gen_kernels = ser_hamil_gen_int_five_moments_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      int_five_moments_hamil_vel_kernels = tensor_hamil_vel_int_five_moments_kernels;
      break;

    default:
      assert(false);
      break;    
  }  

  assert(cv_index[cdim].vdim[vdim] != -1);   

  if (inp->mom_type == GKYL_F_MOMENT_M0M1M2) { // Zeroth, First, and Second moment computed together
    if (inp->model_id == GKYL_MODEL_DEFAULT || inp->model_id == GKYL_MODEL_SR) {
      assert(NULL != int_five_moments_hamil_vel_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
      mom_vlasov->momt.kernel = int_five_moments_hamil_vel_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    }
    else {
      if (inp->conf_basis->b_type == GKYL_BASIS_MODAL_TENSOR) {
        gkyl_exit("mom_vlasov integrated moments hamil_gen: Not currently supported in tensor basis.");
      }
      assert(NULL != int_five_moments_hamil_gen_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
      mom_vlasov->momt.kernel = int_five_moments_hamil_gen_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    } 
    mom_vlasov->momt.num_mom = 2+vdim;
  }
  else {
    gkyl_exit("int_mom_vlasov: Unrecognized moment requested!");
  }

  mom_vlasov->momt.flags = 0;
  GKYL_CLEAR_CU_ALLOC(mom_vlasov->momt.flags);
  mom_vlasov->momt.ref_count = gkyl_ref_count_init(gkyl_mom_free);
  
  mom_vlasov->momt.on_dev = &mom_vlasov->momt; // on host, self-reference
    
  return &mom_vlasov->momt;  
}
