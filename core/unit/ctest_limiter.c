#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_limiter.h>
#include <gkyl_limiter_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

#include <acutest.h>

#include <float.h>
#include <math.h>

static double
cell_min_from_modal(const struct gkyl_basis *basis, const double *modal_vals)
{
  double fmin = DBL_MAX;
  const int nbasis = basis->num_basis;
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
test_cdim_poly_order(int ndim, int poly_order, bool ser,
  bool gkhyb, bool use_gpu, enum gkyl_limiter_type lim_type)
{
  (void)use_gpu;

  struct gkyl_basis basis;
  if (ser) {
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  }
  else if (gkhyb) {
    if (ndim == 2) {
      gkyl_cart_modal_gkhybrid(&basis, 1, 1);
    }
    else {
      gkyl_cart_modal_gkhybrid(&basis, ndim - 2, 2);
    }
  }

  double lower[] = { 0.0 };
  double upper[] = { 1.0 };
  int cells[] = { 32 };

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  int ghost[] = { 0 };
  struct gkyl_range range, range_ext;
  gkyl_create_grid_ranges(&grid, ghost, &range_ext, &range);

  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, range.volume);

  const double c0 = sqrt(2.0);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&range, iter.idx);
    double *fc = gkyl_array_fetch(f, lidx);

    fc[0] = c0;

    double slope = -1.1 - 0.1 * lidx;
    for (int k = 1; k < basis.num_basis; ++k) {
      fc[k] = slope;
    }
  }

  struct gkyl_array *dfdt = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, range.volume);
  gkyl_array_clear(dfdt, -c0);

  struct gkyl_limiter_inp inp = {
    .basis = basis,
    .type = lim_type,
    .dt_factor = 0.9
  };

  struct gkyl_limiter *up = gkyl_limiter_new(inp);
  double dt;
  gkyl_limiter_advance(up, &range, f, dfdt, &dt);

  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&range, iter.idx);
    const double *after = gkyl_array_cfetch(f, lidx);
    TEST_CHECK(cell_min_from_modal(&basis, after) >= -1e-13);
  }
  TEST_CHECK(gkyl_compare_double(dt, 0.9, 1e-13));

  gkyl_limiter_release(up);
  gkyl_array_release(f);
}

void test_1d_p1_ser_zs(void){test_cdim_poly_order(1, 1, true, false, false, GKYL_LIMITER_ZS);}
void test_1d_p2_ser_zs(void){test_cdim_poly_order(1, 2, true, false, false, GKYL_LIMITER_ZS);}
void test_2d_p1_ser_zs(void){test_cdim_poly_order(2, 1, true, false, false, GKYL_LIMITER_ZS);}
void test_2d_p2_ser_zs(void){test_cdim_poly_order(2, 2, true, false, false, GKYL_LIMITER_ZS);}
void test_3d_p1_ser_zs(void){test_cdim_poly_order(3, 1, true, false, false, GKYL_LIMITER_ZS);}
void test_3d_p2_ser_zs(void){test_cdim_poly_order(3, 2, true, false, false, GKYL_LIMITER_ZS);}
void test_4d_p1_ser_zs(void){test_cdim_poly_order(4, 1, true, false, false, GKYL_LIMITER_ZS);}
void test_4d_p2_ser_zs(void){test_cdim_poly_order(4, 2, true, false, false, GKYL_LIMITER_ZS);}
void test_5d_p1_ser_zs(void){test_cdim_poly_order(5, 1, true, false, false, GKYL_LIMITER_ZS);}
void test_5d_p2_ser_zs(void){test_cdim_poly_order(5, 2, true, false, false, GKYL_LIMITER_ZS);}
void test_6d_p1_ser_zs(void){test_cdim_poly_order(6, 1, true, false, false, GKYL_LIMITER_ZS);}

void test_2d_p1_gkhyb_zs(void){test_cdim_poly_order(2, 1, false, true, false, GKYL_LIMITER_ZS);}
void test_3d_p1_gkhyb_zs(void){test_cdim_poly_order(3, 1, false, true, false, GKYL_LIMITER_ZS);}
void test_4d_p1_gkhyb_zs(void){test_cdim_poly_order(4, 1, false, true, false, GKYL_LIMITER_ZS);}
void test_5d_p1_gkhyb_zs(void){test_cdim_poly_order(5, 1, false, true, false, GKYL_LIMITER_ZS);}


void test_1d_p1_ser_mrs(void){test_cdim_poly_order(1, 1, true, false, false, GKYL_LIMITER_MRS);}
void test_1d_p2_ser_mrs(void){test_cdim_poly_order(1, 2, true, false, false, GKYL_LIMITER_MRS);}
void test_2d_p1_ser_mrs(void){test_cdim_poly_order(2, 1, true, false, false, GKYL_LIMITER_MRS);}
void test_2d_p2_ser_mrs(void){test_cdim_poly_order(2, 2, true, false, false, GKYL_LIMITER_MRS);}
void test_3d_p1_ser_mrs(void){test_cdim_poly_order(3, 1, true, false, false, GKYL_LIMITER_MRS);}
void test_3d_p2_ser_mrs(void){test_cdim_poly_order(3, 2, true, false, false, GKYL_LIMITER_MRS);}
void test_4d_p1_ser_mrs(void){test_cdim_poly_order(4, 1, true, false, false, GKYL_LIMITER_MRS);}
void test_4d_p2_ser_mrs(void){test_cdim_poly_order(4, 2, true, false, false, GKYL_LIMITER_MRS);}
void test_5d_p1_ser_mrs(void){test_cdim_poly_order(5, 1, true, false, false, GKYL_LIMITER_MRS);}
void test_5d_p2_ser_mrs(void){test_cdim_poly_order(5, 2, true, false, false, GKYL_LIMITER_MRS);}
void test_6d_p1_ser_mrs(void){test_cdim_poly_order(6, 1, true, false, false, GKYL_LIMITER_MRS);}

void test_2d_p1_gkhyb_mrs(void){test_cdim_poly_order(2, 1, false, true, false, GKYL_LIMITER_MRS);}
void test_3d_p1_gkhyb_mrs(void){test_cdim_poly_order(3, 1, false, true, false, GKYL_LIMITER_MRS);}
void test_4d_p1_gkhyb_mrs(void){test_cdim_poly_order(4, 1, false, true, false, GKYL_LIMITER_MRS);}
void test_5d_p1_gkhyb_mrs(void){test_cdim_poly_order(5, 1, false, true, false, GKYL_LIMITER_MRS);}

TEST_LIST = {
  { "test_1d_p1_ser_zs", test_1d_p1_ser_zs },
  { "test_1d_p2_ser_zs", test_1d_p2_ser_zs },
  { "test_2d_p1_ser_zs", test_2d_p1_ser_zs },
  { "test_2d_p2_ser_zs", test_2d_p2_ser_zs },
  { "test_3d_p1_ser_zs", test_3d_p1_ser_zs },
  { "test_3d_p2_ser_zs", test_3d_p2_ser_zs },
  { "test_4d_p1_ser_zs", test_4d_p1_ser_zs },
  { "test_4d_p2_ser_zs", test_4d_p2_ser_zs },
  { "test_5d_p1_ser_zs", test_5d_p1_ser_zs },
  { "test_5d_p2_ser_zs", test_5d_p2_ser_zs },
  { "test_6d_p1_ser_zs", test_6d_p1_ser_zs },
  { "test_2d_p1_gkhyb_zs", test_2d_p1_gkhyb_zs },
  { "test_3d_p1_gkhyb_zs", test_3d_p1_gkhyb_zs },
  { "test_4d_p1_gkhyb_zs", test_4d_p1_gkhyb_zs },
  { "test_5d_p1_gkhyb_zs", test_5d_p1_gkhyb_zs },
  { "test_1d_p1_ser_mrs", test_1d_p1_ser_mrs },
  { "test_1d_p2_ser_mrs", test_1d_p2_ser_mrs },
  { "test_2d_p1_ser_mrs", test_2d_p1_ser_mrs },
  { "test_2d_p2_ser_mrs", test_2d_p2_ser_mrs },
  { "test_3d_p1_ser_mrs", test_3d_p1_ser_mrs },
  { "test_3d_p2_ser_mrs", test_3d_p2_ser_mrs },
  { "test_4d_p1_ser_mrs", test_4d_p1_ser_mrs },
  { "test_4d_p2_ser_mrs", test_4d_p2_ser_mrs },
  { "test_5d_p1_ser_mrs", test_5d_p1_ser_mrs },
  { "test_5d_p2_ser_mrs", test_5d_p2_ser_mrs },
  { "test_6d_p1_ser_mrs", test_6d_p1_ser_mrs },
  { "test_2d_p1_gkhyb_mrs", test_2d_p1_gkhyb_mrs },
  { "test_3d_p1_gkhyb_mrs", test_3d_p1_gkhyb_mrs },
  { "test_4d_p1_gkhyb_mrs", test_4d_p1_gkhyb_mrs },
  { "test_5d_p1_gkhyb_mrs", test_5d_p1_gkhyb_mrs },
  { NULL, NULL },
};