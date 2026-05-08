#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_const.h>
#include <gkyl_positivity_fdot_restrict.h>
#include <gkyl_positivity_fdot_restrict_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

#include <acutest.h>

#include <float.h>
#include <math.h>

enum test_basis_type {
  SERENDIPITY,
  GKHYBRID
};

static double
cell_min_from_modal(const struct gkyl_basis *basis, const double *modal_vals)
{
  double fmin = DBL_MAX;
  const int nbasis = basis->num_quad;
  double *fquad = gkyl_malloc(nbasis * sizeof(double));

  for (int k = 0; k < nbasis; ++k) {
    fquad[k] = 0.0;
    basis->modal_to_quad_nodal(modal_vals, fquad, k);
    if (fquad[k] < fmin)
      fmin = fquad[k];
  }

  gkyl_free(fquad);
  return fmin;
}

static void
test_ndim_poly_order(int ndim, int poly_order, enum test_basis_type basis_type,
  bool use_gpu, enum gkyl_positivity_fdot_restrict_type lim_type)
{
  struct gkyl_basis basis, *basis_on_dev;
  switch (basis_type) {
    case SERENDIPITY:
      gkyl_cart_modal_serendip(&basis, ndim, poly_order);
      break;
    case GKHYBRID:
      if (ndim == 2) {
        gkyl_cart_modal_gkhybrid(&basis, 1, 1);
      }
      else {
        gkyl_cart_modal_gkhybrid(&basis, ndim - 2, 2);
      }
      break;
  }
  if (use_gpu) {
    switch (basis_type) {
      case SERENDIPITY:
        basis_on_dev = gkyl_cart_modal_serendip_cu_dev_new(ndim, poly_order);
        break;
      case GKHYBRID:
        if (ndim == 2) {
          basis_on_dev = gkyl_cart_modal_gkhybrid_cu_dev_new(1, 1);
        }
        else {
          basis_on_dev = gkyl_cart_modal_gkhybrid_cu_dev_new(ndim - 2, 2);
        }
        break;
    }
  }
  else {
    basis_on_dev = &basis;
  }

  double *lower = gkyl_malloc(ndim * sizeof(double));
  double *upper = gkyl_malloc(ndim * sizeof(double));
  int *cells = gkyl_malloc(ndim * sizeof(int));
  int *ghost = gkyl_malloc(ndim * sizeof(int));
  int num_cells = 4;
  if (use_gpu) {
    num_cells = 8;
  }
  for (int d = 0; d < ndim; ++d) {
    lower[d] = 0.0;
    upper[d] = 1.0;
    cells[d] = num_cells;
    ghost[d] = 0;
  }

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  struct gkyl_range range, range_ext;
  gkyl_create_grid_ranges(&grid, ghost, &range_ext, &range);

  // Host arrays: used for initialization and result verification.
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, range.volume);
  struct gkyl_array *dfdt = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, range.volume);
  struct gkyl_array *f_adv = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, range.volume);
  struct gkyl_array *dfdt_adv = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, range.volume);

  // Device arrays: input/output for the GPU operator.
  struct gkyl_array *f_dev = NULL, *dfdt_adv_dev = NULL;
  if (use_gpu) {
    f_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, range.volume);
    dfdt_adv_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, range.volume);
  }

  const double c0 = pow(sqrt(2.0), ndim);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&range, iter.idx);
    double *fc = gkyl_array_fetch(f, lidx);
    double *dfdtc = gkyl_array_fetch(dfdt, lidx);

    if (lim_type == GKYL_POSITIVITY_FDOT_RESTRICT_QUAD ||
      lim_type == GKYL_POSITIVITY_FDOT_RESTRICT_AVG) {
      fc[0] = c0;
      dfdtc[0] = -c0 * (0.9 + 0.32 * (double)iter.idx[0] / num_cells);

      double slope = -0.64 * (double)iter.idx[0] / num_cells;
      fc[1] = slope;
      for (int k = 2; k < basis.num_basis; ++k) {
        fc[k] = 0.0;
        dfdtc[k] = 0.0;
      }
    }
    else if (lim_type == GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD ||
      lim_type == GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG) {
      fc[0] = c0 * (1.0 - (double)iter.idx[0] / 16);
      dfdtc[0] = -c0* cos(4 * GKYL_PI * (double)iter.idx[0] / num_cells);

      double slope = -0.64 * (double)iter.idx[0] / num_cells;
      fc[1] = slope;
      for (int k = 2; k < basis.num_basis; ++k) {
        fc[k] = 0.0;
        dfdtc[k] = 0.0;
      }
    }
  }

  struct gkyl_positivity_fdot_restrict_inp inp = {
    .basis = basis_on_dev,
    .type = lim_type,
    .safety_factor = 0.0,
    .use_gpu = use_gpu,
  };

  struct gkyl_positivity_fdot_restrict *up = gkyl_positivity_fdot_restrict_new(inp);
  double dt = 1.0;

  if (use_gpu) {
    // Copy initialized host data to device, run restriction on GPU, copy result back.
    gkyl_array_copy(f_dev, f);
    gkyl_array_copy(dfdt_adv_dev, dfdt);
    gkyl_positivity_fdot_restrict_advance(up, &range, f_dev, dfdt_adv_dev, dt);
    gkyl_array_copy(dfdt_adv, dfdt_adv_dev);
  }
  else {
    gkyl_array_copy(dfdt_adv, dfdt);
    gkyl_positivity_fdot_restrict_advance(up, &range, f, dfdt_adv, dt);
  }

  // Compute the Euler stepped solution (f + dt*dfdt) to check positivity.
  gkyl_array_copy(f_adv, f);
  gkyl_array_accumulate(f_adv, 1.0, gkyl_array_scale(dfdt_adv, dt));

  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&range, iter.idx);
    const double *f_before = gkyl_array_cfetch(f, lidx);
    const double *f_after = gkyl_array_cfetch(f_adv, lidx);
    const double *dfdt_before = gkyl_array_cfetch(dfdt, lidx);
    if (lim_type == GKYL_POSITIVITY_FDOT_RESTRICT_AVG) {
      TEST_CHECK(f_after[0] >= -1e-13);
    }
    else if (lim_type == GKYL_POSITIVITY_FDOT_RESTRICT_QUAD) {
      double fquad_after[basis.num_quad];
      for (int k = 0; k < basis.num_quad; ++k) {
        basis.modal_to_quad_nodal(f_after, fquad_after, k);
        TEST_CHECK(fquad_after[k] >= -1e-13);
      }
    }
    else if (lim_type == GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG) {
      if (f_before[0] >= 0.0) {
        // If f_before[0] was positive, it should have been updated maximally
        TEST_CHECK(gkyl_compare_double(f_after[0], f_before[0] + dfdt_before[0] * dt, 1e-14));
      }
      else {
        // If f_before was negative, then it should not have been made more negative
        TEST_CHECK(f_after[0] >= f_before[0]);
      }
    }
    else if (lim_type == GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD) {
      // Loop over each quadrature point
      double fquad_before[basis.num_quad];
      double fquad_after[basis.num_quad];
      double dfdt_quad_before[basis.num_quad];
      for (int k = 0; k < basis.num_quad; ++k) {
        // Get quadrature point k
        basis.modal_to_quad_nodal(f_before, fquad_before, k);
        basis.modal_to_quad_nodal(f_after, fquad_after, k);
        basis.modal_to_quad_nodal(dfdt_before, dfdt_quad_before, k);
        if (fquad_before[k] >= 0.0) {
          // If f was positive at this quadrature point before, it should have been updated maximally
          TEST_CHECK(gkyl_compare_double(fquad_after[k], fquad_before[k] + dfdt_quad_before[k] * dt,
            1e-13));
        }
        else {
          // If f was negative at this quadrature point, then it should not have been made more negative
          TEST_CHECK(fquad_after[k] >= fquad_before[k] - 1e-13); // allow small tolerance for numerical errors
        }
      }
    }
  }

  // char file_name[256];
  // sprintf(file_name, "ctest_positivity_fdot_restrict_f_input.gkyl");
  // gkyl_grid_sub_array_write(&grid, &range, 0, f, file_name);
  // sprintf(file_name, "ctest_positivity_fdot_restrict_f_output.gkyl");
  // gkyl_grid_sub_array_write(&grid, &range, 0, f_adv, file_name);
  // sprintf(file_name, "ctest_positivity_fdot_restrict_dfdt_input.gkyl");
  // gkyl_grid_sub_array_write(&grid, &range, 0, dfdt, file_name);
  // sprintf(file_name, "ctest_positivity_fdot_restrict_dfdt_output.gkyl");
  // gkyl_grid_sub_array_write(&grid, &range, 0, dfdt_adv, file_name);

  if (use_gpu) {
    gkyl_array_release(f_dev);
    gkyl_array_release(dfdt_adv_dev);
    gkyl_cart_modal_basis_release_cu(basis_on_dev);
  }
  gkyl_positivity_fdot_restrict_release(up);
  gkyl_array_release(f_adv);
  gkyl_array_release(f);
  gkyl_array_release(dfdt);
  gkyl_array_release(dfdt_adv);

  gkyl_free(lower);
  gkyl_free(upper);
  gkyl_free(cells);
  gkyl_free(ghost);
}

// All combinations of
// 1d, 2d, 3d, 4d, 5d, 6d
// p1, p2
// serendipity, gkhybrid (only p1, 2d-5d)
// positivity_fdot_restrict modes: avg, quad, diode_avg, diode_quad

void test_1d_p1_ser_fdot_restrict_avg(void){test_ndim_poly_order(1, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_1d_p1_ser_fdot_restrict_quad(void){test_ndim_poly_order(1, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_1d_p1_ser_fdot_restrict_diode_avg(void){test_ndim_poly_order(1, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_1d_p1_ser_fdot_restrict_diode_quad(void){test_ndim_poly_order(1, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_1d_p2_ser_fdot_restrict_avg(void){test_ndim_poly_order(1, 2, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_1d_p2_ser_fdot_restrict_quad(void){test_ndim_poly_order(1, 2, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_1d_p2_ser_fdot_restrict_diode_avg(void){test_ndim_poly_order(1, 2, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_1d_p2_ser_fdot_restrict_diode_quad(void){test_ndim_poly_order(1, 2, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_2d_p1_ser_fdot_restrict_avg(void){test_ndim_poly_order(2, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_2d_p1_ser_fdot_restrict_quad(void){test_ndim_poly_order(2, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_2d_p1_ser_fdot_restrict_diode_avg(void){test_ndim_poly_order(2, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_2d_p1_ser_fdot_restrict_diode_quad(void){test_ndim_poly_order(2, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_2d_p2_ser_fdot_restrict_avg(void){test_ndim_poly_order(2, 2, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_2d_p2_ser_fdot_restrict_quad(void){test_ndim_poly_order(2, 2, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_2d_p2_ser_fdot_restrict_diode_avg(void){test_ndim_poly_order(2, 2, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_2d_p2_ser_fdot_restrict_diode_quad(void){test_ndim_poly_order(2, 2, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_3d_p1_ser_fdot_restrict_avg(void){test_ndim_poly_order(3, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_3d_p1_ser_fdot_restrict_quad(void){test_ndim_poly_order(3, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_3d_p1_ser_fdot_restrict_diode_avg(void){test_ndim_poly_order(3, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_3d_p1_ser_fdot_restrict_diode_quad(void){test_ndim_poly_order(3, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_3d_p2_ser_fdot_restrict_avg(void){test_ndim_poly_order(3, 2, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_3d_p2_ser_fdot_restrict_quad(void){test_ndim_poly_order(3, 2, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_3d_p2_ser_fdot_restrict_diode_avg(void){test_ndim_poly_order(3, 2, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_3d_p2_ser_fdot_restrict_diode_quad(void){test_ndim_poly_order(3, 2, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_4d_p1_ser_fdot_restrict_avg(void){test_ndim_poly_order(4, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_4d_p1_ser_fdot_restrict_quad(void){test_ndim_poly_order(4, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_4d_p1_ser_fdot_restrict_diode_avg(void){test_ndim_poly_order(4, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_4d_p1_ser_fdot_restrict_diode_quad(void){test_ndim_poly_order(4, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_4d_p2_ser_fdot_restrict_avg(void){test_ndim_poly_order(4, 2, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_4d_p2_ser_fdot_restrict_quad(void){test_ndim_poly_order(4, 2, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_4d_p2_ser_fdot_restrict_diode_avg(void){test_ndim_poly_order(4, 2, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_4d_p2_ser_fdot_restrict_diode_quad(void){test_ndim_poly_order(4, 2, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_5d_p1_ser_fdot_restrict_avg(void){test_ndim_poly_order(5, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_5d_p1_ser_fdot_restrict_quad(void){test_ndim_poly_order(5, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_5d_p1_ser_fdot_restrict_diode_avg(void){test_ndim_poly_order(5, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_5d_p1_ser_fdot_restrict_diode_quad(void){test_ndim_poly_order(5, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_5d_p2_ser_fdot_restrict_avg(void){test_ndim_poly_order(5, 2, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_5d_p2_ser_fdot_restrict_quad(void){test_ndim_poly_order(5, 2, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_5d_p2_ser_fdot_restrict_diode_avg(void){test_ndim_poly_order(5, 2, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_5d_p2_ser_fdot_restrict_diode_quad(void){test_ndim_poly_order(5, 2, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_6d_p1_ser_fdot_restrict_avg(void){test_ndim_poly_order(6, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_6d_p1_ser_fdot_restrict_quad(void){test_ndim_poly_order(6, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_6d_p1_ser_fdot_restrict_diode_avg(void){test_ndim_poly_order(6, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_6d_p1_ser_fdot_restrict_diode_quad(void){test_ndim_poly_order(6, 1, SERENDIPITY, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_2d_p1_gkhybrid_fdot_restrict_avg(void){test_ndim_poly_order(2, 1, GKHYBRID, false, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_2d_p1_gkhybrid_fdot_restrict_quad(void){test_ndim_poly_order(2, 1, GKHYBRID, false, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_2d_p1_gkhybrid_fdot_restrict_diode_avg(void){test_ndim_poly_order(2, 1, GKHYBRID, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_2d_p1_gkhybrid_fdot_restrict_diode_quad(void){test_ndim_poly_order(2, 1, GKHYBRID, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_3d_p1_gkhybrid_fdot_restrict_avg(void){test_ndim_poly_order(3, 1, GKHYBRID, false, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_3d_p1_gkhybrid_fdot_restrict_quad(void){test_ndim_poly_order(3, 1, GKHYBRID, false, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_3d_p1_gkhybrid_fdot_restrict_diode_avg(void){test_ndim_poly_order(3, 1, GKHYBRID, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_3d_p1_gkhybrid_fdot_restrict_diode_quad(void){test_ndim_poly_order(3, 1, GKHYBRID, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_4d_p1_gkhybrid_fdot_restrict_avg(void){test_ndim_poly_order(4, 1, GKHYBRID, false, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_4d_p1_gkhybrid_fdot_restrict_quad(void){test_ndim_poly_order(4, 1, GKHYBRID, false, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_4d_p1_gkhybrid_fdot_restrict_diode_avg(void){test_ndim_poly_order(4, 1, GKHYBRID, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_4d_p1_gkhybrid_fdot_restrict_diode_quad(void){test_ndim_poly_order(4, 1, GKHYBRID, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_5d_p1_gkhybrid_fdot_restrict_avg(void){test_ndim_poly_order(5, 1, GKHYBRID, false, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_5d_p1_gkhybrid_fdot_restrict_quad(void){test_ndim_poly_order(5, 1, GKHYBRID, false, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_5d_p1_gkhybrid_fdot_restrict_diode_avg(void){test_ndim_poly_order(5, 1, GKHYBRID, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_5d_p1_gkhybrid_fdot_restrict_diode_quad(void){test_ndim_poly_order(5, 1, GKHYBRID, false, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

#ifdef GKYL_HAVE_CUDA

void test_1d_p1_ser_fdot_restrict_avg_cu(void){test_ndim_poly_order(1, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_1d_p1_ser_fdot_restrict_quad_cu(void){test_ndim_poly_order(1, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_1d_p1_ser_fdot_restrict_diode_avg_cu(void){test_ndim_poly_order(1, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_1d_p1_ser_fdot_restrict_diode_quad_cu(void){test_ndim_poly_order(1, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_1d_p2_ser_fdot_restrict_avg_cu(void){test_ndim_poly_order(1, 2, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_1d_p2_ser_fdot_restrict_quad_cu(void){test_ndim_poly_order(1, 2, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_1d_p2_ser_fdot_restrict_diode_avg_cu(void){test_ndim_poly_order(1, 2, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_1d_p2_ser_fdot_restrict_diode_quad_cu(void){test_ndim_poly_order(1, 2, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_2d_p1_ser_fdot_restrict_avg_cu(void){test_ndim_poly_order(2, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_2d_p1_ser_fdot_restrict_quad_cu(void){test_ndim_poly_order(2, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_2d_p1_ser_fdot_restrict_diode_avg_cu(void){test_ndim_poly_order(2, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_2d_p1_ser_fdot_restrict_diode_quad_cu(void){test_ndim_poly_order(2, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_2d_p2_ser_fdot_restrict_avg_cu(void){test_ndim_poly_order(2, 2, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_2d_p2_ser_fdot_restrict_quad_cu(void){test_ndim_poly_order(2, 2, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_2d_p2_ser_fdot_restrict_diode_avg_cu(void){test_ndim_poly_order(2, 2, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_2d_p2_ser_fdot_restrict_diode_quad_cu(void){test_ndim_poly_order(2, 2, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_3d_p1_ser_fdot_restrict_avg_cu(void){test_ndim_poly_order(3, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_3d_p1_ser_fdot_restrict_quad_cu(void){test_ndim_poly_order(3, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_3d_p1_ser_fdot_restrict_diode_avg_cu(void){test_ndim_poly_order(3, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_3d_p1_ser_fdot_restrict_diode_quad_cu(void){test_ndim_poly_order(3, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_3d_p2_ser_fdot_restrict_avg_cu(void){test_ndim_poly_order(3, 2, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_3d_p2_ser_fdot_restrict_quad_cu(void){test_ndim_poly_order(3, 2, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_3d_p2_ser_fdot_restrict_diode_avg_cu(void){test_ndim_poly_order(3, 2, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_3d_p2_ser_fdot_restrict_diode_quad_cu(void){test_ndim_poly_order(3, 2, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_4d_p1_ser_fdot_restrict_avg_cu(void){test_ndim_poly_order(4, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_4d_p1_ser_fdot_restrict_quad_cu(void){test_ndim_poly_order(4, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_4d_p1_ser_fdot_restrict_diode_avg_cu(void){test_ndim_poly_order(4, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_4d_p1_ser_fdot_restrict_diode_quad_cu(void){test_ndim_poly_order(4, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_4d_p2_ser_fdot_restrict_avg_cu(void){test_ndim_poly_order(4, 2, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_4d_p2_ser_fdot_restrict_quad_cu(void){test_ndim_poly_order(4, 2, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_4d_p2_ser_fdot_restrict_diode_avg_cu(void){test_ndim_poly_order(4, 2, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_4d_p2_ser_fdot_restrict_diode_quad_cu(void){test_ndim_poly_order(4, 2, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_5d_p1_ser_fdot_restrict_avg_cu(void){test_ndim_poly_order(5, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_5d_p1_ser_fdot_restrict_quad_cu(void){test_ndim_poly_order(5, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_5d_p1_ser_fdot_restrict_diode_avg_cu(void){test_ndim_poly_order(5, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_5d_p1_ser_fdot_restrict_diode_quad_cu(void){test_ndim_poly_order(5, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_5d_p2_ser_fdot_restrict_avg_cu(void){test_ndim_poly_order(5, 2, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_5d_p2_ser_fdot_restrict_quad_cu(void){test_ndim_poly_order(5, 2, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_5d_p2_ser_fdot_restrict_diode_avg_cu(void){test_ndim_poly_order(5, 2, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_5d_p2_ser_fdot_restrict_diode_quad_cu(void){test_ndim_poly_order(5, 2, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_6d_p1_ser_fdot_restrict_avg_cu(void){test_ndim_poly_order(6, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_6d_p1_ser_fdot_restrict_quad_cu(void){test_ndim_poly_order(6, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_6d_p1_ser_fdot_restrict_diode_avg_cu(void){test_ndim_poly_order(6, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_6d_p1_ser_fdot_restrict_diode_quad_cu(void){test_ndim_poly_order(6, 1, SERENDIPITY, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_2d_p1_gkhybrid_fdot_restrict_avg_cu(void){test_ndim_poly_order(2, 1, GKHYBRID, true, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_2d_p1_gkhybrid_fdot_restrict_quad_cu(void){test_ndim_poly_order(2, 1, GKHYBRID, true, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_2d_p1_gkhybrid_fdot_restrict_diode_avg_cu(void){test_ndim_poly_order(2, 1, GKHYBRID, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_2d_p1_gkhybrid_fdot_restrict_diode_quad_cu(void){test_ndim_poly_order(2, 1, GKHYBRID, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_3d_p1_gkhybrid_fdot_restrict_avg_cu(void){test_ndim_poly_order(3, 1, GKHYBRID, true, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_3d_p1_gkhybrid_fdot_restrict_quad_cu(void){test_ndim_poly_order(3, 1, GKHYBRID, true, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_3d_p1_gkhybrid_fdot_restrict_diode_avg_cu(void){test_ndim_poly_order(3, 1, GKHYBRID, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_3d_p1_gkhybrid_fdot_restrict_diode_quad_cu(void){test_ndim_poly_order(3, 1, GKHYBRID, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_4d_p1_gkhybrid_fdot_restrict_avg_cu(void){test_ndim_poly_order(4, 1, GKHYBRID, true, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_4d_p1_gkhybrid_fdot_restrict_quad_cu(void){test_ndim_poly_order(4, 1, GKHYBRID, true, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_4d_p1_gkhybrid_fdot_restrict_diode_avg_cu(void){test_ndim_poly_order(4, 1, GKHYBRID, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_4d_p1_gkhybrid_fdot_restrict_diode_quad_cu(void){test_ndim_poly_order(4, 1, GKHYBRID, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

void test_5d_p1_gkhybrid_fdot_restrict_avg_cu(void){test_ndim_poly_order(5, 1, GKHYBRID, true, GKYL_POSITIVITY_FDOT_RESTRICT_AVG);}
void test_5d_p1_gkhybrid_fdot_restrict_quad_cu(void){test_ndim_poly_order(5, 1, GKHYBRID, true, GKYL_POSITIVITY_FDOT_RESTRICT_QUAD);}
void test_5d_p1_gkhybrid_fdot_restrict_diode_avg_cu(void){test_ndim_poly_order(5, 1, GKHYBRID, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG);}
void test_5d_p1_gkhybrid_fdot_restrict_diode_quad_cu(void){test_ndim_poly_order(5, 1, GKHYBRID, true, GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD);}

#endif // GKYL_HAVE_CUDA

TEST_LIST = {
  { "test_1d_p1_ser_fdot_restrict_avg", test_1d_p1_ser_fdot_restrict_avg },
  { "test_1d_p1_ser_fdot_restrict_quad", test_1d_p1_ser_fdot_restrict_quad },
  { "test_1d_p1_ser_fdot_restrict_diode_avg", test_1d_p1_ser_fdot_restrict_diode_avg },
  { "test_1d_p1_ser_fdot_restrict_diode_quad", test_1d_p1_ser_fdot_restrict_diode_quad },
  { "test_1d_p2_ser_fdot_restrict_avg", test_1d_p2_ser_fdot_restrict_avg },
  { "test_1d_p2_ser_fdot_restrict_quad", test_1d_p2_ser_fdot_restrict_quad },
  { "test_1d_p2_ser_fdot_restrict_diode_avg", test_1d_p2_ser_fdot_restrict_diode_avg },
  { "test_1d_p2_ser_fdot_restrict_diode_quad", test_1d_p2_ser_fdot_restrict_diode_quad },
  { "test_2d_p1_ser_fdot_restrict_avg", test_2d_p1_ser_fdot_restrict_avg },
  { "test_2d_p1_ser_fdot_restrict_quad", test_2d_p1_ser_fdot_restrict_quad },
  { "test_2d_p1_ser_fdot_restrict_diode_avg", test_2d_p1_ser_fdot_restrict_diode_avg },
  { "test_2d_p1_ser_fdot_restrict_diode_quad", test_2d_p1_ser_fdot_restrict_diode_quad },
  { "test_2d_p2_ser_fdot_restrict_avg", test_2d_p2_ser_fdot_restrict_avg },
  { "test_2d_p2_ser_fdot_restrict_quad", test_2d_p2_ser_fdot_restrict_quad },
  { "test_2d_p2_ser_fdot_restrict_diode_avg", test_2d_p2_ser_fdot_restrict_diode_avg },
  { "test_2d_p2_ser_fdot_restrict_diode_quad", test_2d_p2_ser_fdot_restrict_diode_quad },
  { "test_3d_p1_ser_fdot_restrict_avg", test_3d_p1_ser_fdot_restrict_avg },
  { "test_3d_p1_ser_fdot_restrict_quad", test_3d_p1_ser_fdot_restrict_quad },
  { "test_3d_p1_ser_fdot_restrict_diode_avg", test_3d_p1_ser_fdot_restrict_diode_avg },
  { "test_3d_p1_ser_fdot_restrict_diode_quad", test_3d_p1_ser_fdot_restrict_diode_quad },
  { "test_3d_p2_ser_fdot_restrict_avg", test_3d_p2_ser_fdot_restrict_avg },
  { "test_3d_p2_ser_fdot_restrict_quad", test_3d_p2_ser_fdot_restrict_quad },
  { "test_3d_p2_ser_fdot_restrict_diode_avg", test_3d_p2_ser_fdot_restrict_diode_avg },
  { "test_3d_p2_ser_fdot_restrict_diode_quad", test_3d_p2_ser_fdot_restrict_diode_quad },
  { "test_4d_p1_ser_fdot_restrict_avg", test_4d_p1_ser_fdot_restrict_avg },
  { "test_4d_p1_ser_fdot_restrict_quad", test_4d_p1_ser_fdot_restrict_quad },
  { "test_4d_p1_ser_fdot_restrict_diode_avg", test_4d_p1_ser_fdot_restrict_diode_avg },
  { "test_4d_p1_ser_fdot_restrict_diode_quad", test_4d_p1_ser_fdot_restrict_diode_quad },
  { "test_4d_p2_ser_fdot_restrict_avg", test_4d_p2_ser_fdot_restrict_avg },
  { "test_4d_p2_ser_fdot_restrict_quad", test_4d_p2_ser_fdot_restrict_quad },
  { "test_4d_p2_ser_fdot_restrict_diode_avg", test_4d_p2_ser_fdot_restrict_diode_avg },
  { "test_4d_p2_ser_fdot_restrict_diode_quad", test_4d_p2_ser_fdot_restrict_diode_quad },
  { "test_5d_p1_ser_fdot_restrict_avg", test_5d_p1_ser_fdot_restrict_avg },
  { "test_5d_p1_ser_fdot_restrict_quad", test_5d_p1_ser_fdot_restrict_quad },
  { "test_5d_p1_ser_fdot_restrict_diode_avg", test_5d_p1_ser_fdot_restrict_diode_avg },
  { "test_5d_p1_ser_fdot_restrict_diode_quad", test_5d_p1_ser_fdot_restrict_diode_quad },
  { "test_5d_p2_ser_fdot_restrict_avg", test_5d_p2_ser_fdot_restrict_avg },
  { "test_5d_p2_ser_fdot_restrict_quad", test_5d_p2_ser_fdot_restrict_quad },
  { "test_5d_p2_ser_fdot_restrict_diode_avg", test_5d_p2_ser_fdot_restrict_diode_avg },
  { "test_5d_p2_ser_fdot_restrict_diode_quad", test_5d_p2_ser_fdot_restrict_diode_quad },
  { "test_6d_p1_ser_fdot_restrict_avg", test_6d_p1_ser_fdot_restrict_avg },
  { "test_6d_p1_ser_fdot_restrict_quad", test_6d_p1_ser_fdot_restrict_quad },
  { "test_6d_p1_ser_fdot_restrict_diode_avg", test_6d_p1_ser_fdot_restrict_diode_avg },
  { "test_6d_p1_ser_fdot_restrict_diode_quad", test_6d_p1_ser_fdot_restrict_diode_quad },
  { "test_2d_p1_gkhybrid_fdot_restrict_avg", test_2d_p1_gkhybrid_fdot_restrict_avg },
  { "test_2d_p1_gkhybrid_fdot_restrict_quad", test_2d_p1_gkhybrid_fdot_restrict_quad },
  { "test_2d_p1_gkhybrid_fdot_restrict_diode_avg", test_2d_p1_gkhybrid_fdot_restrict_diode_avg },
  { "test_2d_p1_gkhybrid_fdot_restrict_diode_quad", test_2d_p1_gkhybrid_fdot_restrict_diode_quad },
  { "test_3d_p1_gkhybrid_fdot_restrict_avg", test_3d_p1_gkhybrid_fdot_restrict_avg },
  { "test_3d_p1_gkhybrid_fdot_restrict_quad", test_3d_p1_gkhybrid_fdot_restrict_quad },
  { "test_3d_p1_gkhybrid_fdot_restrict_diode_avg", test_3d_p1_gkhybrid_fdot_restrict_diode_avg },
  { "test_3d_p1_gkhybrid_fdot_restrict_diode_quad", test_3d_p1_gkhybrid_fdot_restrict_diode_quad },
  { "test_4d_p1_gkhybrid_fdot_restrict_avg", test_4d_p1_gkhybrid_fdot_restrict_avg },
  { "test_4d_p1_gkhybrid_fdot_restrict_quad", test_4d_p1_gkhybrid_fdot_restrict_quad },
  { "test_4d_p1_gkhybrid_fdot_restrict_diode_avg", test_4d_p1_gkhybrid_fdot_restrict_diode_avg },
  { "test_4d_p1_gkhybrid_fdot_restrict_diode_quad", test_4d_p1_gkhybrid_fdot_restrict_diode_quad },
  { "test_5d_p1_gkhybrid_fdot_restrict_avg", test_5d_p1_gkhybrid_fdot_restrict_avg },
  { "test_5d_p1_gkhybrid_fdot_restrict_quad", test_5d_p1_gkhybrid_fdot_restrict_quad },
  { "test_5d_p1_gkhybrid_fdot_restrict_diode_avg", test_5d_p1_gkhybrid_fdot_restrict_diode_avg },
  { "test_5d_p1_gkhybrid_fdot_restrict_diode_quad", test_5d_p1_gkhybrid_fdot_restrict_diode_quad },
#ifdef GKYL_HAVE_CUDA
  { "test_1d_p1_ser_fdot_restrict_avg_cu", test_1d_p1_ser_fdot_restrict_avg_cu },
  { "test_1d_p1_ser_fdot_restrict_quad_cu", test_1d_p1_ser_fdot_restrict_quad_cu },
  { "test_1d_p1_ser_fdot_restrict_diode_avg_cu", test_1d_p1_ser_fdot_restrict_diode_avg_cu },
  { "test_1d_p1_ser_fdot_restrict_diode_quad_cu", test_1d_p1_ser_fdot_restrict_diode_quad_cu },
  { "test_1d_p2_ser_fdot_restrict_avg_cu", test_1d_p2_ser_fdot_restrict_avg_cu },
  { "test_1d_p2_ser_fdot_restrict_quad_cu", test_1d_p2_ser_fdot_restrict_quad_cu },
  { "test_1d_p2_ser_fdot_restrict_diode_avg_cu", test_1d_p2_ser_fdot_restrict_diode_avg_cu },
  { "test_1d_p2_ser_fdot_restrict_diode_quad_cu", test_1d_p2_ser_fdot_restrict_diode_quad_cu },
  { "test_2d_p1_ser_fdot_restrict_avg_cu", test_2d_p1_ser_fdot_restrict_avg_cu },
  { "test_2d_p1_ser_fdot_restrict_quad_cu", test_2d_p1_ser_fdot_restrict_quad_cu },
  { "test_2d_p1_ser_fdot_restrict_diode_avg_cu", test_2d_p1_ser_fdot_restrict_diode_avg_cu },
  { "test_2d_p1_ser_fdot_restrict_diode_quad_cu", test_2d_p1_ser_fdot_restrict_diode_quad_cu },
  { "test_2d_p2_ser_fdot_restrict_avg_cu", test_2d_p2_ser_fdot_restrict_avg_cu },
  { "test_2d_p2_ser_fdot_restrict_quad_cu", test_2d_p2_ser_fdot_restrict_quad_cu },
  { "test_2d_p2_ser_fdot_restrict_diode_avg_cu", test_2d_p2_ser_fdot_restrict_diode_avg_cu },
  { "test_2d_p2_ser_fdot_restrict_diode_quad_cu", test_2d_p2_ser_fdot_restrict_diode_quad_cu },
  { "test_3d_p1_ser_fdot_restrict_avg_cu", test_3d_p1_ser_fdot_restrict_avg_cu },
  { "test_3d_p1_ser_fdot_restrict_quad_cu", test_3d_p1_ser_fdot_restrict_quad_cu },
  { "test_3d_p1_ser_fdot_restrict_diode_avg_cu", test_3d_p1_ser_fdot_restrict_diode_avg_cu },
  { "test_3d_p1_ser_fdot_restrict_diode_quad_cu", test_3d_p1_ser_fdot_restrict_diode_quad_cu },
  { "test_3d_p2_ser_fdot_restrict_avg_cu", test_3d_p2_ser_fdot_restrict_avg_cu },
  { "test_3d_p2_ser_fdot_restrict_quad_cu", test_3d_p2_ser_fdot_restrict_quad_cu },
  { "test_3d_p2_ser_fdot_restrict_diode_avg_cu", test_3d_p2_ser_fdot_restrict_diode_avg_cu },
  { "test_3d_p2_ser_fdot_restrict_diode_quad_cu", test_3d_p2_ser_fdot_restrict_diode_quad_cu },
  { "test_4d_p1_ser_fdot_restrict_avg_cu", test_4d_p1_ser_fdot_restrict_avg_cu },
  { "test_4d_p1_ser_fdot_restrict_quad_cu", test_4d_p1_ser_fdot_restrict_quad_cu },
  { "test_4d_p1_ser_fdot_restrict_diode_avg_cu", test_4d_p1_ser_fdot_restrict_diode_avg_cu },
  { "test_4d_p1_ser_fdot_restrict_diode_quad_cu", test_4d_p1_ser_fdot_restrict_diode_quad_cu },
  { "test_4d_p2_ser_fdot_restrict_avg_cu", test_4d_p2_ser_fdot_restrict_avg_cu },
  { "test_4d_p2_ser_fdot_restrict_quad_cu", test_4d_p2_ser_fdot_restrict_quad_cu },
  { "test_4d_p2_ser_fdot_restrict_diode_avg_cu", test_4d_p2_ser_fdot_restrict_diode_avg_cu },
  { "test_4d_p2_ser_fdot_restrict_diode_quad_cu", test_4d_p2_ser_fdot_restrict_diode_quad_cu },
  { "test_5d_p1_ser_fdot_restrict_avg_cu", test_5d_p1_ser_fdot_restrict_avg_cu },
  { "test_5d_p1_ser_fdot_restrict_quad_cu", test_5d_p1_ser_fdot_restrict_quad_cu },
  { "test_5d_p1_ser_fdot_restrict_diode_avg_cu", test_5d_p1_ser_fdot_restrict_diode_avg_cu },
  { "test_5d_p1_ser_fdot_restrict_diode_quad_cu", test_5d_p1_ser_fdot_restrict_diode_quad_cu },
  { "test_5d_p2_ser_fdot_restrict_avg_cu", test_5d_p2_ser_fdot_restrict_avg_cu },
  { "test_5d_p2_ser_fdot_restrict_quad_cu", test_5d_p2_ser_fdot_restrict_quad_cu },
  { "test_5d_p2_ser_fdot_restrict_diode_avg_cu", test_5d_p2_ser_fdot_restrict_diode_avg_cu },
  { "test_5d_p2_ser_fdot_restrict_diode_quad_cu", test_5d_p2_ser_fdot_restrict_diode_quad_cu },
  { "test_6d_p1_ser_fdot_restrict_avg_cu", test_6d_p1_ser_fdot_restrict_avg_cu },
  { "test_6d_p1_ser_fdot_restrict_quad_cu", test_6d_p1_ser_fdot_restrict_quad_cu },
  { "test_6d_p1_ser_fdot_restrict_diode_avg_cu", test_6d_p1_ser_fdot_restrict_diode_avg_cu },
  { "test_6d_p1_ser_fdot_restrict_diode_quad_cu", test_6d_p1_ser_fdot_restrict_diode_quad_cu },
  { "test_2d_p1_gkhybrid_fdot_restrict_avg_cu", test_2d_p1_gkhybrid_fdot_restrict_avg_cu },
  { "test_2d_p1_gkhybrid_fdot_restrict_quad_cu", test_2d_p1_gkhybrid_fdot_restrict_quad_cu },
  { "test_2d_p1_gkhybrid_fdot_restrict_diode_avg_cu",
    test_2d_p1_gkhybrid_fdot_restrict_diode_avg_cu },
  { "test_2d_p1_gkhybrid_fdot_restrict_diode_quad_cu",
    test_2d_p1_gkhybrid_fdot_restrict_diode_quad_cu },
  { "test_3d_p1_gkhybrid_fdot_restrict_avg_cu", test_3d_p1_gkhybrid_fdot_restrict_avg_cu },
  { "test_3d_p1_gkhybrid_fdot_restrict_quad_cu", test_3d_p1_gkhybrid_fdot_restrict_quad_cu },
  { "test_3d_p1_gkhybrid_fdot_restrict_diode_avg_cu",
    test_3d_p1_gkhybrid_fdot_restrict_diode_avg_cu },
  { "test_3d_p1_gkhybrid_fdot_restrict_diode_quad_cu",
    test_3d_p1_gkhybrid_fdot_restrict_diode_quad_cu },
  { "test_4d_p1_gkhybrid_fdot_restrict_avg_cu", test_4d_p1_gkhybrid_fdot_restrict_avg_cu },
  { "test_4d_p1_gkhybrid_fdot_restrict_quad_cu", test_4d_p1_gkhybrid_fdot_restrict_quad_cu },
  { "test_4d_p1_gkhybrid_fdot_restrict_diode_avg_cu",
    test_4d_p1_gkhybrid_fdot_restrict_diode_avg_cu },
  { "test_4d_p1_gkhybrid_fdot_restrict_diode_quad_cu",
    test_4d_p1_gkhybrid_fdot_restrict_diode_quad_cu },
  { "test_5d_p1_gkhybrid_fdot_restrict_avg_cu", test_5d_p1_gkhybrid_fdot_restrict_avg_cu },
  { "test_5d_p1_gkhybrid_fdot_restrict_quad_cu", test_5d_p1_gkhybrid_fdot_restrict_quad_cu },
  { "test_5d_p1_gkhybrid_fdot_restrict_diode_avg_cu",
    test_5d_p1_gkhybrid_fdot_restrict_diode_avg_cu },
  { "test_5d_p1_gkhybrid_fdot_restrict_diode_quad_cu",
    test_5d_p1_gkhybrid_fdot_restrict_diode_quad_cu },
#endif // GKYL_HAVE_CUDA
  { NULL, NULL },
};
