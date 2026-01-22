#include <acutest.h>

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_triad_geom.h>
#include <gkyl_vlasov_triad_geom_priv.h>
#include <gkyl_vlasov_triad_geom_preset_priv.h>
#include <gkyl_util.h>

#include <gkyl_vlasov_kernels.h>

void
test_triad_2x2v_rphi_ks_pnt_alpha(int poly_order)
{

  // Black hole parameters
  struct vm_geom geom;
  geom.mass_bh = 0.7;
  geom.spin_bh = 0.2;

  double lower[] = {1.2, 1.0, 1.0, 1.0}, upper[] = {1.201, 1.001, 1.001, 1.001};
  int cells[] = {2, 2, 2, 2};
  int vdim = 2, cdim = 2;
  int ndim = cdim+vdim;

  double confLower[] = {lower[0], lower[1]}, confUpper[] = {upper[0], upper[1]};
  int confCells[] = {cells[0], cells[1]};
  double velLower[] = { lower[2], lower[3]}, velUpper[] = { upper[2], upper[3]};
  int velCells[] = {cells[2], cells[3]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

    struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = { 0, 0 };
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = { 0, 0 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = { confGhost[0], confGhost[1], 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Construct the input map
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom;
  inp_triad_geom.use_vierbein = true;
  inp_triad_geom.use_preset_geom = false;
  inp_triad_geom.eval_vierbein = eval_ks_rphi_vierbein_2v; 
  inp_triad_geom.eval_vierbein_gradient = eval_ks_rphi_vierbein_gradient_2v; 
  inp_triad_geom.eval_vierbein_ctx = &geom; 
  inp_triad_geom.eval_vierbein_gradient_ctx = &geom;  

  // Make the memory for arrays (modal)
  struct gkyl_array *conf_poisson_tensor; // Configuration space Poisson tensor representation

  // Size of the PT
  int num_pt_indices[3] = { 1 , 6, 18 }; 

  // Allocate arrays for covariant tangent basis 
  conf_poisson_tensor = mkarr(false, basis.num_basis*num_pt_indices[vdim-1], local_ext.volume);

  // Construct the Geometry for this configuration
  gkyl_vlasov_triad_geom_new(&confGrid, &confLocal, confBasis, 
    &grid, &local, basis, inp_triad_geom, conf_poisson_tensor);

  // Iterate over the grid, conf space, checking output
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *conf_poisson_tensor_d = gkyl_array_cfetch(conf_poisson_tensor, gkyl_range_idx(&confLocal, iter.idx));
    
    int NC = confBasis.num_basis;

    // index size of the symmetric vierbein
    int n_sym = vdim*(vdim+1)/2; 

    // Double precision values of the expected output
    // pt at the first and second points of conflocal
    double conf_pt_pnt1_vals[48] = {
      1.3588086393098029e+00,  4.3989488197243091e-05,  0.0000000000000000e+00,  0.0000000000000000e+00,
     -2.8209032673875084e-09,  1.5884597076836006e-16,  0.0000000000000000e+00, -4.4892812507152318e-21,

      3.3326390817298623e-01, -4.0077069914578266e-05,  0.0000000000000000e+00,  4.6259292692714869e-18,
      4.3107081734112447e-09,  4.0163355468972541e-17,  0.0000000000000000e+00,  3.3710617091548934e-21,

      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,

      1.6663195408649309e+00, -2.0038534957297145e-04,  0.0000000000000000e+00,  0.0000000000000000e+00,
      2.1553540958340271e-08,  1.8375096125123659e-16,  0.0000000000000000e+00,  1.7028157575878893e-20,

      7.4747990717562870e-02, -2.2127022055885313e-05,  0.0000000000000000e+00,  0.0000000000000000e+00,
      3.9991364151008178e-09,  7.6698926571441654e-18,  0.0000000000000000e+00,  1.9701838207745808e-21,

      1.1321046895815821e+00, -9.9492398515537507e-05,  0.0000000000000000e+00, -1.8503717077085947e-17,
      8.3511852151311742e-09,  1.3368141403477171e-16,  0.0000000000000000e+00,  9.9229411160645424e-21
    };

    double conf_pt_pnt2_vals[48] = {
      1.3588086393098029e+00,  4.3989488197243091e-05,  0.0000000000000000e+00,  0.0000000000000000e+00,
     -2.8209032673875084e-09,  1.5884597076836006e-16,  0.0000000000000000e+00, -4.4892812507152318e-21,

      3.3326390817298623e-01, -4.0077069914578266e-05,  0.0000000000000000e+00,  4.6259292692714869e-18,
      4.3107081734112447e-09,  4.0163355468972541e-17,  0.0000000000000000e+00,  3.3710617091548934e-21,

      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,

      1.6663195408649309e+00, -2.0038534957297145e-04,  0.0000000000000000e+00,  0.0000000000000000e+00,
      2.1553540958340271e-08,  1.8375096125123659e-16,  0.0000000000000000e+00,  1.7028157575878893e-20,

      7.4747990717562870e-02, -2.2127022055885313e-05,  0.0000000000000000e+00,  0.0000000000000000e+00,
      3.9991364151008178e-09,  7.6698926571441654e-18,  0.0000000000000000e+00,  1.9701838207745808e-21,

      1.1321046895815821e+00, -9.9492398515537507e-05,  0.0000000000000000e+00, -1.8503717077085947e-17,
      8.3511852151311742e-09,  1.3368141403477171e-16,  0.0000000000000000e+00,  9.9229411160645424e-21
    };

    double conf_pt_pnt3_vals[48] = {
      1.3589609855304008e+00,  4.3967646125575739e-05,  0.0000000000000000e+00,  0.0000000000000000e+00,
     -2.8186954269208275e-09,  1.2945958434809000e-16,  0.0000000000000000e+00, -3.5482011072808203e-21,

      3.3312513494050400e-01, -4.0043700164429930e-05,  0.0000000000000000e+00,  4.6259292692714869e-18,
      4.3053254227575864e-09,  3.0521154308908952e-17,  0.0000000000000000e+00,  4.6283179257499702e-21,

      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,

      1.6656256747025202e+00, -2.0021850082239002e-04,  0.0000000000000000e+00, -1.8503717077085947e-17,
      2.1526627144731539e-08,  2.0492889448999093e-16,  0.0000000000000000e+00,  2.1966070824554527e-20,

      7.4671394086855941e-02, -2.2096072509741311e-05,  0.0000000000000000e+00,  0.0000000000000000e+00,
      3.9920039782016186e-09,  8.0867145517377461e-18,  0.0000000000000000e+00,  2.4846788383995600e-21,

      1.1317601498076699e+00, -9.9427743788891502e-05,  0.0000000000000000e+00,  1.8503717077085947e-17,
      8.3425950860921977e-09,  9.8703706655908211e-17,  0.0000000000000000e+00,  9.1561817091435301e-21
    };

    double conf_pt_pnt4_vals[48] = {
      1.3589609855304008e+00,  4.3967646125575739e-05,  0.0000000000000000e+00,  0.0000000000000000e+00,
     -2.8186954269208275e-09,  1.2945958434809000e-16,  0.0000000000000000e+00, -3.5482011072808203e-21,

      3.3312513494050400e-01, -4.0043700164429930e-05,  0.0000000000000000e+00,  4.6259292692714869e-18,
      4.3053254227575864e-09,  3.0521154308908952e-17,  0.0000000000000000e+00,  4.6283179257499702e-21,

      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
      0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,

      1.6656256747025202e+00, -2.0021850082239002e-04,  0.0000000000000000e+00, -1.8503717077085947e-17,
      2.1526627144731539e-08,  2.0492889448999093e-16,  0.0000000000000000e+00,  2.1966070824554527e-20,

      7.4671394086855941e-02, -2.2096072509741311e-05,  0.0000000000000000e+00,  0.0000000000000000e+00,
      3.9920039782016186e-09,  8.0867145517377461e-18,  0.0000000000000000e+00,  2.4846788383995600e-21,

      1.1317601498076699e+00, -9.9427743788891502e-05,  0.0000000000000000e+00,  1.8503717077085947e-17,
      8.3425950860921977e-09,  9.8703706655908211e-17,  0.0000000000000000e+00,  9.1561817091435301e-21
    };


    for (int k = 0; k<num_pt_indices[vdim-1]; ++k) {
      for (int m = 0; m<NC; ++m) {
        int test_idx = k*NC + m;
        double expected = 0;
        if (iter.idx[0] == 1) expected = conf_pt_pnt1_vals[test_idx];
        if (iter.idx[0] == 3) expected = conf_pt_pnt2_vals[test_idx];
        if (iter.idx[0] == 2) expected = conf_pt_pnt3_vals[test_idx];
        if (iter.idx[0] == 4) expected = conf_pt_pnt4_vals[test_idx];
        //printf("conf_poisson_tensor_d[%d]: %1.16e\n",test_idx,conf_poisson_tensor_d[test_idx]);
        TEST_CHECK( gkyl_compare_double(conf_poisson_tensor_d[test_idx], expected, 1e-12) );
      }
    }
  }

  // Build the hamiltonian 
  // Hamiltonian (GR) is a full phase-space array. 
  struct gkyl_range hamil_range = local; 
  struct gkyl_array *hamil = mkarr(false, basis.num_basis, local_ext.volume);


  // Evaluate specified hamiltonian function at nodes to ensure continuity of hamiltonian
  struct gkyl_eval_on_nodes* hamil_proj = gkyl_eval_on_nodes_new(&grid, &basis, 1, gkyl_vlasov_triad_preset_hamil(vdim,GKYL_TRIAD_GR_KERR_SCHILD_RPHI), &geom);
  gkyl_eval_on_nodes_advance(hamil_proj, 0.0, &local_ext, hamil);
  gkyl_eval_on_nodes_release(hamil_proj);

  
  // Now build the charateristics components in Velocity and Configuration Space
  gkyl_range_iter_init(&iter, &local);
  int pdim = cdim + vdim;
  int idx[GKYL_MAX_DIM], conf_idx[GKYL_MAX_DIM];

  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idx);
    gkyl_copy_int_arr(cdim, iter.idx, conf_idx);
    
    // printf("idx[0]: %d, idx[1]: %d, idx[2]: %d, idx[3]: %d\n",idx[0], idx[1], idx[2], idx[3]);
    // Only select the first cell which doesn't require ghost cells to compute
    if (idx[0] == 1 && idx[1] == 1 && idx[2] == 1 && idx[3] == 1) {

      // Make room for charateristics:
      double alpha_quad_r[8] = {0.0};
      double alpha_quad_phi[8] = {0.0};
      double alpha_quad_vr[8] = {0.0};
      double alpha_quad_vphi[8] = {0.0};

      // Grab the cell center location for NC bracket calculation 
      double xcC[GKYL_MAX_DIM], xcR[GKYL_MAX_DIM];
      long cidx = gkyl_range_idx(&confLocal, conf_idx);

      // Assuming the hamil is dependant on the phase space
      long hidx = gkyl_range_idx(&hamil_range, idx);
      const double *hamil_d = gkyl_array_cfetch(hamil, hidx); 
      const double *conf_poisson_tensor_d = gkyl_array_cfetch(conf_poisson_tensor, cidx);

      gkyl_rect_grid_cell_center(&grid, idx, xcC);
      hamil_phase_ho_alpha_quad_x_2x2v_ser_p2(xcC, grid.dx, -1, conf_poisson_tensor_d, hamil_d, alpha_quad_r);
      hamil_phase_ho_alpha_quad_y_2x2v_ser_p2(xcC, grid.dx, -1, conf_poisson_tensor_d, hamil_d, alpha_quad_phi);
      nc_hamil_phase_ho_alpha_quad_vx_2x2v_ser_p2(xcC, grid.dx, conf_poisson_tensor_d, hamil_d, alpha_quad_vr);
      nc_hamil_phase_ho_alpha_quad_vy_2x2v_ser_p2(xcC, grid.dx, conf_poisson_tensor_d, hamil_d, alpha_quad_vphi);

      printf("\n");
      for (int m = 0; m<1; ++m) {
        printf("alpha(r): %1.16e, alpha(phi): %1.16e, alpha(vr): %1.16e, alpha(vphi): %1.16e\n", alpha_quad_r[m], alpha_quad_phi[m], alpha_quad_vr[m], alpha_quad_vphi[m]);
      }
    }
  }


  // Release the memory
  gkyl_array_release(conf_poisson_tensor);
  gkyl_array_release(hamil);
}

void test_triad_2x2v_rphi_ks_alpha() { test_triad_2x2v_rphi_ks_pnt_alpha(2); }

TEST_LIST = {
  { "test_triad_2x2v_rphi_ks_alpha", test_triad_2x2v_rphi_ks_alpha},
  {NULL, NULL}
};
