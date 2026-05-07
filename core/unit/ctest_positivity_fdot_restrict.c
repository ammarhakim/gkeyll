#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_positivity_fdot_restrict.h>
#include <gkyl_range.h>

#include <acutest.h>

#include <math.h>
#include <stdio.h>

void
test_positivity_fdot_restrict_avg_1x1v()
{
  // Test cell-average mode for 1x1v with DG1 (linear basis).
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, 1); // 1x1v, p=1
  
  struct gkyl_positivity_fdot_restrict_inp inp = {
    .basis = basis,
    .mode = GKYL_POSITIVITY_FDOT_RESTRICT_AVG,
    .safety_factor = 0.9,
  };
  
  struct gkyl_positivity_fdot_restrict *up = gkyl_positivity_fdot_restrict_new(inp);
  
  // Create simple arrays for f and dfdt.
  struct gkyl_array *f = gkyl_array_new(basis.num_basis, 1);
  struct gkyl_array *dfdt = gkyl_array_new(basis.num_basis, 1);
  
  // Test case: f = [1.0, 0.1], dfdt = [-0.5, 0.01]
  // Cell average: f_avg = 1.0, dfdt_avg = -0.5
  // With dt = 2.0: f_pred = 1.0 + 2.0 * (-0.5) = 0.0 (boundary case)
  // With dt = 3.0: f_pred = 1.0 + 3.0 * (-0.5) = -0.5 (would be negative)
  // After restriction: dfdt_avg should be clamped to reach 0.9 * 1.0 = 0.9
  // dfdt_new = (0.9 - 1.0) / 3.0 = -0.0333...
  
  double *f_ptr = gkyl_array_fetch(f, 0);
  double *dfdt_ptr = gkyl_array_fetch(dfdt, 0);
  
  f_ptr[0] = 1.0; // Cell average
  f_ptr[1] = 0.1; // Slope
  
  dfdt_ptr[0] = -0.5; // Cell average of dfdt
  dfdt_ptr[1] = 0.01; // Slope of dfdt
  
  double dt = 3.0;
  
  struct gkyl_range range;
  gkyl_range_init_1d(&range, 0, 0); // Single cell range.
  
  // Apply restriction.
  gkyl_positivity_fdot_restrict_advance(up, &range, f, dfdt, dt);
  
  // Check that dfdt[0] was restricted.
  double expected_dfdt_avg = (0.9 - 1.0) / dt; // = -0.0333...
  TEST_CHECK(fabs(dfdt_ptr[0] - expected_dfdt_avg) < 1e-10);
  TEST_CHECK(fabs(dfdt_ptr[1] - 0.01) < 1e-10); // Slope unchanged in avg mode.
  
  gkyl_array_release(f);
  gkyl_array_release(dfdt);
  gkyl_positivity_fdot_restrict_release(up);
}

void
test_positivity_fdot_restrict_avg_safe()
{
  // Test case where prediction is already safe (no restriction needed).
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, 1); // 1x1v, p=1
  
  struct gkyl_positivity_fdot_restrict_inp inp = {
    .basis = basis,
    .mode = GKYL_POSITIVITY_FDOT_RESTRICT_AVG,
    .safety_factor = 0.9,
  };
  
  struct gkyl_positivity_fdot_restrict *up = gkyl_positivity_fdot_restrict_new(inp);
  
  struct gkyl_array *f = gkyl_array_new(basis.num_basis, 1);
  struct gkyl_array *dfdt = gkyl_array_new(basis.num_basis, 1);
  
  double *f_ptr = gkyl_array_fetch(f, 0);
  double *dfdt_ptr = gkyl_array_fetch(dfdt, 0);
  
  f_ptr[0] = 1.0;
  f_ptr[1] = 0.1;
  
  dfdt_ptr[0] = 0.1; // Positive dfdt, safe.
  dfdt_ptr[1] = 0.01;
  
  double dt = 1.0;
  double dfdt_avg_orig = dfdt_ptr[0];
  
  struct gkyl_range range;
  gkyl_range_init_1d(&range, 0, 0);
  
  gkyl_positivity_fdot_restrict_advance(up, &range, f, dfdt, dt);
  
  // dfdt should remain unchanged since prediction is safe.
  TEST_CHECK(fabs(dfdt_ptr[0] - dfdt_avg_orig) < 1e-10);
  
  gkyl_array_release(f);
  gkyl_array_release(dfdt);
  gkyl_positivity_fdot_restrict_release(up);
}

int
main(int argc, char* argv[])
{
  acutest_suite_begin();
  acutest_run_test(test_positivity_fdot_restrict_avg_1x1v);
  acutest_run_test(test_positivity_fdot_restrict_avg_safe);
  acutest_suite_end();
  
  return acutest_suite_retval();
}
