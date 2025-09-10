#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gkgeom.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

#include <math.h>

// Elliptical "equilibrium"



static inline double sq(double x) { return x*x; }
static inline double cub(double x) { return x*x*x; }
static inline double qad(double x) { return x*x*x*x; }
static inline double pen(double x) { return x*x*x*x*x; }
static inline double hex(double x) { return x*x*x*x*x*x; }

void
psi_ellip(double t, const double *xn, double *fout, void *ctx)
{
  double R = xn[0], Z = xn[1];
  fout[0] = (R-2)*(R-2) + Z*Z/4;
}

void
ellip_unit(void)
{
  // create RZ grid
  double lower[] = { 0.5, -4.0 }, upper[] = { 6.0, 4.0 };
  // as ellipitical surfaces are exact, we only need 1 cell in each
  // direction
  int cells[] = { 1, 1 };

  struct gkyl_rect_grid rzgrid;
  gkyl_rect_grid_init(&rzgrid, 2, lower, upper, cells);

  // RZ ranges
  struct gkyl_range rzlocal, rzlocal_ext;
  int nghost[GKYL_MAX_CDIM] = { 0, 0 };
  gkyl_create_grid_ranges(&rzgrid, nghost, &rzlocal_ext, &rzlocal);

  // RZ basis function
  int rz_poly_order = 2;
  struct gkyl_basis rzbasis;
  gkyl_cart_modal_serendip(&rzbasis, 2, rz_poly_order);

  // allocate psiRZ array, initialize and write it to file
  struct gkyl_array *psiRZ = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  
  gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&rzgrid,
    &rzbasis, 1, &psi_ellip, 0);
  gkyl_eval_on_nodes_advance(eon, 0.0, &rzlocal, psiRZ);
  gkyl_eval_on_nodes_release(eon);

  gkyl_grid_sub_array_write(&rzgrid, &rzlocal, 0, psiRZ, "ellip_psi.gkyl");

  gkyl_gkgeom *geo = gkyl_gkgeom_new(&(struct gkyl_gkgeom_inp) {
      // psiRZ and related inputs
      .rzgrid = &rzgrid,
      .rzbasis = &rzbasis,
      .psiRZ = psiRZ,
      .rzlocal = &rzlocal,

    }
  );

  // exact values computed with the following Maxima code
  /*
    Zlo : -4.0$
    Zup : 4.0$
    psi : (R-2)^2 + Z^2/4$

    psi0 : 10.1$
    
    Zmin : max(-2*sqrt(psi0)+1e-12, Zlo)$
    Zmax : min(2*sqrt(psi0)-1e-12, Zup)$
    
    R : 2 + sqrt(psi0-Z^2/4)$
    fR : sqrt(1+diff(R,Z)^2)$
    Ipsi : quad_qag(fR,Z,Zmin,Zmax, 3, 'epsrel=1e-12)$
  */
  
  do {
    double psi_ref = 6.0;
    double arcL = gkyl_gkgeom_integrate_psi_contour(geo, psi_ref,
      lower[1], upper[1], upper[0]);

    TEST_CHECK( gkyl_compare(8.382428377712543, arcL, 1e-12) );
    
  } while(0);

  do {
    double psi_ref = 10.1;
    double arcL = gkyl_gkgeom_integrate_psi_contour(geo, psi_ref,
      lower[1], upper[1], upper[0]);

    TEST_CHECK( gkyl_compare(8.172574228918158, arcL, 1e-12) );
    
  } while(0);    
  
  gkyl_gkgeom_release(geo);
  gkyl_array_release(psiRZ);
}

TEST_LIST = {
  { "ellip", ellip_unit },
  { NULL, NULL }
};
