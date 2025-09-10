#include <acutest.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_util.h>

void test_g1()
{
  double w[gkyl_gauss_max], x[gkyl_gauss_max];
  
  for (unsigned n=1; n<=gkyl_gauss_max; ++n) {
    // use the generic routine to get the ordinates and weights
    gkyl_gauleg(-1, 1, &x[0], &w[0], n);

    // compare with pre-computed values
    const double *xp = gkyl_gauss_ordinates[n];
    const double *wp = gkyl_gauss_weights[n];
    
    for (unsigned i=0; i<n; ++i) {
      TEST_CHECK( gkyl_compare(x[i], xp[i], 1e-12) );
      TEST_CHECK( gkyl_compare(w[i], wp[i], 1e-12) );
    }
  }
}

TEST_LIST = {
  { "basic", test_g1 },
  { NULL, NULL },
};
