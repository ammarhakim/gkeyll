#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>


static inline double sq(double x) { return x*x; }
static inline double cub(double x) { return x*x*x; }

void
test_cubic_1d(void)
{
  double val[2] = { 1.0, 2.0 };
  double grad[2] = { -1.0, -2.0 };
  double coeff[4] = { 0.0 };
  
  gkyl_dg_calc_cubic_1d(val, grad, coeff);

  struct gkyl_basis b3;
  gkyl_cart_modal_serendip(&b3, 1, 3);

  TEST_CHECK( gkyl_compare_double(val[0], b3.eval_expand((double[1]) { -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(val[1], b3.eval_expand((double[1]) { 1.0 }, coeff), 1.0e-15) );

  TEST_CHECK( gkyl_compare_double(grad[0], b3.eval_grad_expand(0, (double[1]) { -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(grad[1], b3.eval_grad_expand(0, (double[1]) { 1.0 }, coeff), 1.0e-15) );
}

void
test_cubic_2d(void)
{
  double val[4] = { 1.0, 2.0, 3.0, 4.0 };
  double gradx[4] = { -1.0, -2.0, -3.0, -4.0 };
  double grady[4] = { 1.0, 2.0, 3.0, 4.0 };
  double gradxy[4] = { 1.5, 2.5, 3.5, 4.5 };
  double coeff[16] = { 0.0 };
  
  gkyl_dg_calc_cubic_2d(val, gradx, grady, gradxy, coeff);

  struct gkyl_basis b3;
  gkyl_cart_modal_tensor(&b3, 2, 3);

  TEST_CHECK( gkyl_compare_double(val[0], b3.eval_expand((double[2]) { -1.0, -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(val[1], b3.eval_expand((double[2]) { -1.0, 1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(val[2], b3.eval_expand((double[2]) { 1.0, -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(val[3], b3.eval_expand((double[2]) { 1.0, 1.0 }, coeff), 1.0e-15) );

  TEST_CHECK( gkyl_compare_double(gradx[0], b3.eval_grad_expand(0, (double[2]) { -1.0, -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(gradx[1], b3.eval_grad_expand(0, (double[2]) { -1.0, 1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(gradx[2], b3.eval_grad_expand(0, (double[2]) { 1.0, -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(gradx[3], b3.eval_grad_expand(0, (double[2]) { 1.0, 1.0 }, coeff), 1.0e-15) );

  TEST_CHECK( gkyl_compare_double(grady[0], b3.eval_grad_expand(1, (double[2]) { -1.0, -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(grady[1], b3.eval_grad_expand(1, (double[2]) { -1.0, 1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(grady[2], b3.eval_grad_expand(1, (double[2]) { 1.0, -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(grady[3], b3.eval_grad_expand(1, (double[2]) { 1.0, 1.0 }, coeff), 1.0e-15) );
}

void
test_cubic_evalf_2d(void)
{
  double lower[] = { 0.0, 0.0 }, upper[] = { 5.0, 5.0 };
  int cells[] = { 8, 8 };

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  // nodal grid used in IO so we can plot things
  double nc_lower[] = { lower[0] - 0.5*grid.dx[0], lower[1] - 0.5*grid.dx[1] };
  double nc_upper[] = { upper[0] + 0.5*grid.dx[0], upper[1] + 0.5*grid.dx[1] };
  int nc_cells[] = { cells[0] + 1, cells[1] + 1 };
  struct gkyl_rect_grid nc_grid;
  gkyl_rect_grid_init(&nc_grid, 2, nc_lower, nc_upper, nc_cells);

  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 0, 0 };  
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  struct gkyl_range nc_local, nc_local_ext;
  gkyl_create_grid_ranges(&nc_grid, nghost, &nc_local_ext, &nc_local);

  struct gkyl_basis basis;
  gkyl_cart_modal_tensor(&basis, 2, 3);

  struct gkyl_array *psi_nodal = gkyl_array_new(GKYL_DOUBLE, 1, (cells[0]+1)*(cells[1]+1));
  struct gkyl_array *psi_cubic = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_dg_basis_op_mem *mem = gkyl_dg_alloc_cubic_2d(cells);

  double xn[2];
  
  // initialize 2D nodal values f = x^y^2 which we can represent exactly
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &nc_local);
  while (gkyl_range_iter_next(&iter)) {
    long nidx = gkyl_range_idx(&nc_local, iter.idx);
    
    gkyl_rect_grid_ll_node(&grid, iter.idx, xn);
    
    double *pn = gkyl_array_fetch(psi_nodal, nidx);
    pn[0] = sq(xn[0])*sq(xn[1]);
  }
  // compute cubic expansion
  gkyl_dg_calc_cubic_2d_from_nodal_vals(mem, cells, grid.dx,
    psi_nodal, psi_cubic);
    
  // compute evalf function from nodal values
  struct gkyl_basis_ops_evalf *evf = gkyl_dg_basis_ops_evalf_new(&grid, psi_nodal);

  // Test wgrad function: f, df/dx, and df/dy
  xn[0] = 2.5;
  xn[1] = 2.5;
  double fout[4];
  evf->eval_cubic_wgrad(0.0, xn, fout, evf->ctx);
  TEST_CHECK( gkyl_compare_double(fout[0], sq(xn[0])*sq(xn[1]), 1.0e-14) );
  TEST_CHECK( gkyl_compare_double(fout[1], 2*xn[0]*sq(xn[1]), 1.0e-14) );
  TEST_CHECK( gkyl_compare_double(fout[2], 2*sq(xn[0])*xn[1], 1.0e-14) );

  // Test wgrad2 function: f, d2f/dx2, d2f/dy2, and d2f/dxdy
  evf->eval_cubic_wgrad2(0.0, xn, fout, evf->ctx);
  TEST_CHECK( gkyl_compare_double(fout[0], sq(xn[0])*sq(xn[1]), 1.0e-14) );
  TEST_CHECK( gkyl_compare_double(fout[1], 2*sq(xn[1]), 1.0e-14) );
  TEST_CHECK( gkyl_compare_double(fout[2], 2*sq(xn[0]), 1.0e-14) );
  TEST_CHECK( gkyl_compare_double(fout[3], 4*xn[0]*xn[1], 1.0e-14) );

  gkyl_array_release(psi_nodal);
  gkyl_array_release(psi_cubic);
  
  gkyl_dg_basis_op_mem_release(mem);
  gkyl_dg_basis_ops_evalf_release(evf);
}


TEST_LIST = {
  { "cubic_1d", test_cubic_1d },
  { "cubic_2d", test_cubic_2d },
  { "cubic_evalf_2d", test_cubic_evalf_2d },
  { NULL, NULL },
};

