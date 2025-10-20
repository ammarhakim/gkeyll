#include <gkyl_efit.h>
#include <gkyl_dg_basis_ops.h>

bool 
newton_raphson(struct gkyl_efit *up, const double *coeffs, double *xsol, bool cubics);

int 
find_xpts(gkyl_efit* up, double *Rxpt, double *Zxpt);

int 
find_xpts_cubic(gkyl_efit* up, double *Rxpt, double *Zxpt);

