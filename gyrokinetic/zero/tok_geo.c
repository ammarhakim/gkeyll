#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <assert.h>
#include <gkyl_basis.h>
#include <gkyl_math.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_position_map.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_priv.h>
#include <gkyl_tok_geo_priv.h>
#include <gkyl_dg_bin_ops.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static bool
tok_geo_same_flux(double psi_a, double psi_b)
{
  double scale = fmax(1.0, fmax(fabs(psi_a), fabs(psi_b)));
  return fabs(psi_a-psi_b) <= 64.0*DBL_EPSILON*scale;
}

static bool
tok_geo_trace_enabled(void)
{
  const char *trace = getenv("GKYL_TOK_GEO_TRACE");
  return trace && trace[0] != '\0' && trace[0] != '0';
}

static bool
tok_geo_trace_block_selected(void)
{
  const char *filter = getenv("GKYL_TOK_GEO_TRACE_BLOCK_FILTER");
  const char *block = getenv("GKYL_TOK_GEO_TRACE_BLOCK");
  return !(filter && filter[0] != '\0') || (block && strcmp(filter, block) == 0);
}

static int
tok_geo_trace_int_env(const char *name, int fallback)
{
  const char *val = getenv(name);
  return val && val[0] != '\0' ? atoi(val) : fallback;
}

static double
tok_geo_trace_double_env(const char *name, double fallback)
{
  const char *val = getenv(name);
  return val && val[0] != '\0' ? atof(val) : fallback;
}

static bool
tok_ordered_map_diag_enabled(void)
{
  const char *diag = getenv("GKYL_TOK_ORDERED_MAP_DIAG");
  return diag && diag[0] != '\0' && diag[0] != '0';
}

static bool
tok_xpt_mapping_requested(const struct gkyl_tok_geo_grid_inp *inp)
{
  if (inp->straight_xpt_ray)
    return true;
  return inp->straight_core_xpt_ray &&
    (inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R ||
     inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L);
}

static void
tok_init_xpt_ray_target(const struct gkyl_tok_geo_grid_inp *inp,
  const struct gkyl_position_map *position_map, struct arc_length_ctx *arc_ctx)
{
  if (!tok_xpt_mapping_requested(inp))
    return;
  bool inner_radial_boundary =
    inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R ||
    inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L ||
    inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R ||
    inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_L;
  double psi_comp = inner_radial_boundary
    ? inp->cgrid.lower[0] : inp->cgrid.upper[0];
  position_map->maps[0](0.0, &psi_comp, &arc_ctx->xpt_ray_psi0,
    position_map->ctxs[0]);
}

static double
tok_xpt_theta_to_arc(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double theta)
{
  if (!arc_ctx->xpt_map_valid)
    return (theta+M_PI)*arc_ctx->arcL_tot/(2.0*M_PI);
  double frac = (theta-inp->cgrid.lower[2])
    /(inp->cgrid.upper[2]-inp->cgrid.lower[2]);
  double arc = arc_ctx->xpt_map_arc_lo
    + frac*(arc_ctx->xpt_map_arc_hi-arc_ctx->xpt_map_arc_lo);
  if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R ||
      inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L) {
    arc = fmod(arc, arc_ctx->arcL_tot);
    if (arc < 0.0)
      arc += arc_ctx->arcL_tot;
  }
  return arc;
}

static bool
tok_xpt_at_seam(const struct gkyl_tok_geo_grid_inp *inp, int it,
  const struct gkyl_range *nrange, bool lower_is_global, bool upper_is_global)
{
  bool lower = it == nrange->lower[2] && lower_is_global;
  bool upper = it == nrange->upper[2] && upper_is_global;
  switch (inp->ftype) {
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_R:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID:
    case GKYL_GEOMETRY_TOKAMAK_CORE_L:
      return upper;
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_L:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO:
    case GKYL_GEOMETRY_TOKAMAK_CORE_R:
      return lower;
    default:
      return false;
  }
}

static bool
tok_xpt_at_fixed_edge(const struct gkyl_tok_geo_grid_inp *inp,
  const struct arc_length_ctx *arc_ctx, int it,
  const struct gkyl_range *nrange, bool lower_is_global,
  bool upper_is_global, double *z, double *root_reference)
{
  if (!arc_ctx->xpt_map_valid)
    return false;
  bool lower = it == nrange->lower[2] && lower_is_global;
  bool upper = it == nrange->upper[2] && upper_is_global;
  switch (inp->ftype) {
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_R:
      if (lower) {
        *z = arc_ctx->zmin_right; *root_reference = inp->rright;
        return true;
      }
      break;
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_L:
      if (upper) {
        *z = arc_ctx->zmin_left; *root_reference = inp->rleft;
        return true;
      }
      break;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO:
      if (lower) {
        *z = arc_ctx->zmin; *root_reference = inp->rright;
        return true;
      }
      break;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID:
    case GKYL_GEOMETRY_TOKAMAK_CORE_R:
      if (upper) {
        *z = arc_ctx->geo->zmaxis; *root_reference = inp->rright;
        return true;
      }
      break;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID:
    case GKYL_GEOMETRY_TOKAMAK_CORE_L:
      if (lower) {
        *z = arc_ctx->geo->zmaxis; *root_reference = inp->rleft;
        return true;
      }
      break;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO:
      if (upper) {
        *z = arc_ctx->zmin; *root_reference = inp->rleft;
        return true;
      }
      break;
    default:
      break;
  }
  return false;
}

static void
tok_geo_trace_surface(FILE *fp, const char *stage, int block, int ftype,
  int dir, int ip, int it, int ia, int ip_delta, double psi, double alpha,
  double theta, double arcL, double zmin, double zmax, double rclose,
  double ridders_min, double ridders_max, double res, int nevals, double elapsed)
{
  fprintf(fp, "%s,%d,%d,%d,%d,%d,%d,%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%d,%.17e\n",
    stage, block, ftype, dir, ip, it, ia, ip_delta, psi, alpha, theta, arcL,
    zmin, zmax, rclose, ridders_min, ridders_max, res, nevals, elapsed);
}

static void
tok_geo_trace_surface_row(const char *stage, int ftype, int dir, int ip,
  int it, int ia, int ip_delta, double psi, double alpha, double theta,
  double arcL, double zmin, double zmax, double rclose, double ridders_min,
  double ridders_max, double res, int nevals, double elapsed)
{
  if (!tok_geo_trace_enabled() || !tok_geo_trace_block_selected())
    return;

  const char *block_env = getenv("GKYL_TOK_GEO_TRACE_BLOCK");
  int block = block_env ? atoi(block_env) : -1;
  const char *fname = getenv("GKYL_TOK_GEO_TRACE_FILE");

  if (fname && fname[0] != '\0') {
    static bool wrote_header = false;
    FILE *fp = fopen(fname, "a");
    if (fp) {
      if (!wrote_header) {
        fprintf(fp, "stage,block,ftype,dir,ip,it,ia,ip_delta,psi,alpha,theta,arcL,zmin,zmax,rclose,ridders_min,ridders_max,res,nevals,elapsed_sec\n");
        wrote_header = true;
      }
      tok_geo_trace_surface(fp, stage, block, ftype, dir, ip, it, ia, ip_delta,
        psi, alpha, theta, arcL, zmin, zmax, rclose, ridders_min, ridders_max,
        res, nevals, elapsed);
      fclose(fp);
    }
  }

  if (tok_geo_trace_int_env("GKYL_TOK_GEO_TRACE_STDERR", 1)) {
    tok_geo_trace_surface(stderr, stage, block, ftype, dir, ip, it, ia, ip_delta,
      psi, alpha, theta, arcL, zmin, zmax, rclose, ridders_min, ridders_max,
      res, nevals, elapsed);
    fflush(stderr);
  }
}

double
tok_plate_psi_func(double s, void *ctx){
  // uses a pointer to the plate function to get R(s), Z(s)
  // Then calculates psi(R, Z)
  // will be used by ridders later
  
  struct plate_ctx *gc = ctx;
  double RZ[2];
  if(gc->lower==true)
    gc->geo->plate_func_lower(s, RZ);
  else
    gc->geo->plate_func_upper(s, RZ);

  double R = RZ[0];
  double Z = RZ[1];

  // Now find the cell where this R and Z is
  if (gc->geo->use_cubics) {
    double xn[2] = {R, Z};
    double psi;
    gc->geo->efit->evf->eval_cubic(0.0, xn, &psi, gc->geo->efit->evf->ctx);
    return psi - gc->psi_curr;
  }
  else {
    int rzidx[2];
    rzidx[0] = fmin(gc->geo->rzlocal.lower[0] + (int) floor((R - gc->geo->rzgrid.lower[0])/gc->geo->rzgrid.dx[0]), gc->geo->rzlocal.upper[0]);
    rzidx[1] = fmin(gc->geo->rzlocal.lower[1] + (int) floor((Z - gc->geo->rzgrid.lower[1])/gc->geo->rzgrid.dx[1]), gc->geo->rzlocal.upper[1]);
    long loc = gkyl_range_idx(&gc->geo->rzlocal, rzidx);
    const double *coeffs = gkyl_array_cfetch(gc->geo->psiRZ,loc);

    double xc[2];
    gkyl_rect_grid_cell_center(&gc->geo->rzgrid, rzidx, xc);
    double xy[2];
    xy[0] = (R-xc[0])/(gc->geo->rzgrid.dx[0]*0.5);
    xy[1] = (Z-xc[1])/(gc->geo->rzgrid.dx[1]*0.5);
    double psi = gc->geo->rzbasis.eval_expand(xy, coeffs);
    return psi - gc->psi_curr;
  }
}



// Function to pass to root-finder to find Z location for given arc-length
static inline double
arc_length_func(double Z, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL;
  double zmax = actx->zmax;
  double ival = 0.0;

  if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_CORE){
    if(actx->right==true){
      double *arc_memo = actx->arc_memo_right;
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, false, arc_memo) - arcL;
    }
    else{
      double *arc_memo = actx->arc_memo_left;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, true, false, arc_memo)  - arcL + actx->arcL_right;
    }
  }

  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_CORE_L || actx->ftype==GKYL_GEOMETRY_TOKAMAK_CORE_R){
    if (actx->xpt_map_valid && actx->right) {
      double *arc_memo = actx->arc_memo_right;
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose,
        true, false, arc_memo)-arcL;
    }
    else if (actx->xpt_map_valid) {
      double *arc_memo = actx->arc_memo_left;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose,
        true, false, arc_memo)-arcL+actx->arcL_right;
    }
    else if(actx->pre==true){
      ival = actx->arcL_start - integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo) - arcL ;
    }
    else if(actx->right==true){
      double *arc_memo = actx->arc_memo_right;
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, false, arc_memo) - arcL + actx->arcL_start;
    }
    else{
      double *arc_memo = actx->arc_memo_left;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, true, false, arc_memo)  - arcL + actx->arcL_right + actx->arcL_start;
    }
  }

  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_PF_LO_L || actx->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R){
    if(actx->right==true){
      double *arc_memo = actx->arc_memo_right;
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, false, arc_memo) - arcL;
    }
    else{
      double *arc_memo = actx->arc_memo_left;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, true, false, arc_memo)  - arcL + actx->arcL_right;
    }
  }

  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_PF_UP_L || actx->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_R){
    if(actx->right==false){
      double *arc_memo = actx->arc_memo_left;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, true, false, arc_memo) - arcL;
    }
    else{
      double *arc_memo = actx->arc_memo_right;
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, false, arc_memo)  - arcL + actx->arcL_left;
    }
  }

  else if( (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT) || (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO) || (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID) || (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP) ){
    double *arc_memo = actx->arc_memo;
    ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, false, arc_memo) - arcL;
  }
  else if( (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN) || (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN) || (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO) || (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID) || (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP) ){
    double *arc_memo = actx->arc_memo;
    ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, true, false, arc_memo) - arcL;
  }

  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_LSN_SOL || actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO || actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID || actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP){
    if(actx->right==true){
      double *arc_memo = actx->arc_memo_right;
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo) - arcL;
    }
    else{
      double *arc_memo = actx->arc_memo_left;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo)  - arcL + actx->arcL_right;
    }
  }

  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_IWL){
    if(actx->q3) {
      double *arc_memo = actx->arc_memo;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo) - arcL;
    }
    else if (actx->q4) {
      double *arc_memo = actx->arc_memo;
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo)  - arcL + actx->arcL_q3;
    }
    else if (actx->q1) {
      double *arc_memo = actx->arc_memo;
      ival = integrate_psi_contour_memo(actx->geo, psi, actx->geo->zmaxis, Z, rclose, false, false, arc_memo)  - arcL + actx->arcL_q3 +actx->arcL_q4;
    }
    else {
      double *arc_memo = actx->arc_memo;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo)  - arcL + actx->arcL_q3 +actx->arcL_q4 + actx->arcL_q1;
    }

  }

  return ival;
}

static double
tok_eval_psi_rz_local(const struct gkyl_tok_geo *geo, double R, double Z)
{
  if (geo->use_cubics) {
    double xn[2] = { R, Z }, psi = 0.0;
    geo->efit->evf->eval_cubic(0.0, xn, &psi, geo->efit->evf->ctx);
    return psi;
  }

  int idx[2];
  idx[0] = GKYL_MIN2(geo->rzlocal.upper[0], GKYL_MAX2(geo->rzlocal.lower[0],
    geo->rzlocal.lower[0]
      +(int) floor((R-geo->rzgrid.lower[0])/geo->rzgrid.dx[0])));
  idx[1] = GKYL_MIN2(geo->rzlocal.upper[1], GKYL_MAX2(geo->rzlocal.lower[1],
    geo->rzlocal.lower[1]
      +(int) floor((Z-geo->rzgrid.lower[1])/geo->rzgrid.dx[1])));
  long loc = gkyl_range_idx(&geo->rzlocal, idx);
  const double *coeffs = gkyl_array_cfetch(geo->psiRZ, loc);
  double xc[2], eta[2];
  gkyl_rect_grid_cell_center(&geo->rzgrid, idx, xc);
  eta[0] = (R-xc[0])/(0.5*geo->rzgrid.dx[0]);
  eta[1] = (Z-xc[1])/(0.5*geo->rzgrid.dx[1]);
  return geo->rzbasis.eval_expand(eta, coeffs);
}

static void
tok_append_unique_root(double root, double *roots, int *nr, int nmax)
{
  for (int i=0; i<*nr; ++i)
    if (fabs(root-roots[i]) <= 2e-11*fmax(1.0, fabs(root)))
      return;
  if (*nr < nmax)
    roots[(*nr)++] = root;
}

// Symmetric companion to gkyl_tok_geo_R_psiZ: find every Z root at fixed R.
// The quadratic tensor representation is solved analytically in each DG cell,
// so horizontal contour segments and Z turning points are not discarded.
static int
tok_geo_Z_psiR(const struct gkyl_tok_geo *geo, double psi, double R,
  int nmax, double *Z)
{
  int nr = 0;
  if (!geo->use_cubics && geo->efit->rzbasis.poly_order <= 2) {
    int ridx = GKYL_MIN2(geo->rzlocal.upper[0],
      GKYL_MAX2(geo->rzlocal.lower[0], geo->rzlocal.lower[0]
        +(int) floor((R-geo->rzgrid.lower[0])/geo->rzgrid.dx[0])));
    for (int iz=geo->rzlocal.lower[1]; iz<=geo->rzlocal.upper[1]; ++iz) {
      int idx[2] = { ridx, iz };
      long loc = gkyl_range_idx(&geo->rzlocal, idx);
      const double *p = gkyl_array_cfetch(geo->psiRZ, loc);
      double xc[2];
      gkyl_rect_grid_cell_center(&geo->rzgrid, idx, xc);
      double x = (R-xc[0])/(0.5*geo->rzgrid.dx[0]);
      double roots[2] = { 0.0, 0.0 };
      int ncell = 0;

      if (geo->efit->rzbasis.poly_order == 1) {
        double den = 3.0*p[3]*x+1.732050807568877*p[2];
        double num = -1.732050807568877*p[1]*x+2.0*psi-p[0];
        double scale = fmax(1.0, fabs(num));
        if (fabs(den) > 64.0*DBL_EPSILON*scale)
          roots[ncell++] = num/den;
      }
      else {
        double a = 0.125*(45.0*p[8]*x*x+23.2379000772445*p[7]*x
          -15.0*p[8]+13.41640786499874*p[5]);
        double b = 0.125*(23.2379000772445*p[6]*x*x+12.0*p[3]*x
          -7.745966692414834*p[6]+6.928203230275509*p[2]);
        double c = 0.125*((13.41640786499874*p[4]-15.0*p[8])*x*x
          +(6.928203230275509*p[1]-7.745966692414834*p[7])*x
          +5.0*p[8]-4.47213595499958*p[4]-4.47213595499958*p[5]
          +4.0*p[0])-psi;
        double coeff_scale = fmax(1.0, fmax(fabs(a), fmax(fabs(b), fabs(c))));
        if (fabs(a) <= 64.0*DBL_EPSILON*coeff_scale) {
          if (fabs(b) > 64.0*DBL_EPSILON*coeff_scale)
            roots[ncell++] = -c/b;
        }
        else {
          double disc = b*b-4.0*a*c;
          double disc_tol = 256.0*DBL_EPSILON*
            fmax(1.0, b*b+fabs(4.0*a*c));
          if (disc >= -disc_tol) {
            disc = fmax(0.0, disc);
            double sd = sqrt(disc);
            if (sd == 0.0)
              roots[ncell++] = -0.5*b/a;
            else {
              double q = -0.5*(b+copysign(sd, b));
              roots[ncell++] = q/a;
              roots[ncell++] = c/q;
            }
          }
        }
      }

      for (int k=0; k<ncell; ++k) {
        double y = roots[k];
        if (isfinite(y) && y >= -1.0-2e-12 && y <= 1.0+2e-12) {
          double z = xc[1]+0.5*geo->rzgrid.dx[1]*
            fmin(1.0, fmax(-1.0, y));
          tok_append_unique_root(z, Z, &nr, nmax);
        }
      }
    }
    return nr;
  }

  // Cubic fallback: bracket sign-changing roots on a fine Z scan.  Exact
  // X-point endpoints are supplied explicitly by the trace builder.
  const int nsamp = 8*geo->rzgrid_cubic.cells[1];
  double z0 = geo->rzgrid_cubic.lower[1];
  double f0 = tok_eval_psi_rz_local(geo, R, z0)-psi;
  for (int i=1; i<=nsamp; ++i) {
    double z1 = geo->rzgrid_cubic.lower[1]
      +(geo->rzgrid_cubic.upper[1]-geo->rzgrid_cubic.lower[1])
        *i/(double) nsamp;
    double f1 = tok_eval_psi_rz_local(geo, R, z1)-psi;
    if (isfinite(f0) && isfinite(f1) && f0*f1 <= 0.0) {
      double lo = z0, hi = z1, flo = f0;
      for (int n=0; n<64; ++n) {
        double mid = 0.5*(lo+hi);
        double fm = tok_eval_psi_rz_local(geo, R, mid)-psi;
        if (flo*fm <= 0.0)
          hi = mid;
        else { lo = mid; flo = fm; }
      }
      tok_append_unique_root(0.5*(lo+hi), Z, &nr, nmax);
    }
    z0 = z1; f0 = f1;
  }
  return nr;
}

static double
tok_nearest_value(double ref, const double *values, int n)
{
  int ibest = 0;
  for (int i=1; i<n; ++i)
    if (fabs(values[i]-ref) < fabs(values[ibest]-ref))
      ibest = i;
  return values[ibest];
}

static bool tok_fixed_edge_is_midplane(enum gkyl_tok_geo_type ftype);

static bool
tok_build_contour_candidate(const struct gkyl_tok_geo *geo, double psi,
  bool param_is_r, double rfixed, double zfixed, double rx, double zx, int n,
  double *r, double *z, double *score, const char *name,
  bool cluster_fixed_endpoint)
{
  r[0] = rfixed; z[0] = zfixed;
  double total = 0.0, max_step = 0.0, final_step = 0.0;
  for (int i=1; i<n-1; ++i) {
    double t = i/(double) (n-1);
    // A midplane endpoint is an R turning point of the contour.  There,
    // Z-Z_fixed scales as sqrt(|R-R_fixed|), so uniform R sampling creates an
    // artificially large first physical step.  Quadratic R spacing recovers
    // approximately uniform arc-length spacing without weakening the checks
    // that reject a genuine root-branch jump.
    double f = param_is_r && cluster_fixed_endpoint ? t*t : t;
    if (param_is_r) {
      r[i] = rfixed+f*(rx-rfixed);
      double roots[16] = { 0.0 };
      int nr = tok_geo_Z_psiR(geo, psi, r[i], 16, roots);
      if (nr == 0)
        return false;
      z[i] = tok_nearest_value(z[i-1], roots, nr);
    }
    else {
      z[i] = zfixed+f*(zx-zfixed);
      double roots[8] = { 0.0 }, dRdZ[8] = { 0.0 };
      double dR[8] = { 0.0 }, dZ[8] = { 0.0 };
      int nr = gkyl_tok_geo_R_psiZ(geo, psi, z[i], 8,
        roots, dRdZ, dR, dZ);
      if (nr == 0)
        return false;
      r[i] = tok_nearest_value(r[i-1], roots, nr);
    }
    if (!isfinite(r[i]) || !isfinite(z[i]))
      return false;
    double residual = tok_eval_psi_rz_local(geo, r[i], z[i])-psi;
    if (!isfinite(residual) || fabs(residual) >
        1e-9*fmax(1.0, fabs(psi)))
      return false;
    double ds = hypot(r[i]-r[i-1], z[i]-z[i-1]);
    if (!(ds > 0.0) || !isfinite(ds))
      return false;
    total += ds; max_step = fmax(max_step, ds);
  }
  r[n-1] = rx; z[n-1] = zx;
  final_step = hypot(rx-r[n-2], zx-z[n-2]);
  if (!(final_step > 0.0) || !isfinite(final_step))
    return false;
  total += final_step; max_step = fmax(max_step, final_step);
  double mean = total/(n-1);
  double cell_diag = hypot(geo->rzgrid.dx[0], geo->rzgrid.dx[1]);
  if (max_step > 2.0*cell_diag || max_step > 16.0*mean ||
      final_step > 2.0*cell_diag) {
    if (tok_ordered_map_diag_enabled())
      fprintf(stderr,
        "TOK_ORDERED_MAP_DIAG reason=discontinuous_trace trace=%s param=%c n=%d psi=%.17g max_step=%.17g final_step=%.17g mean_step=%.17g cell_diag=%.17g endpoints=(%.17g,%.17g)->(%.17g,%.17g)\n",
        name, param_is_r ? 'R' : 'Z', n, psi, max_step, final_step, mean,
        cell_diag, rfixed, zfixed, rx, zx);
    return false;
  }
  *score = max_step/mean+4.0*final_step/mean;
  return isfinite(*score);
}

static bool
tok_sep_fixed_edge_is_first(enum gkyl_tok_geo_type ftype)
{
  return ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L ||
    ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID ||
    ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R ||
    ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO;
}

static void
tok_build_sep_trace(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double zfixed, double rfixed_ref)
{
  if (arc_ctx->sep_trace_initialized)
    return;
  if (!arc_ctx->sep_trace_r || !arc_ctx->sep_trace_z ||
      !arc_ctx->sep_trace_s || arc_ctx->sep_trace_capacity < 17) {
    fprintf(stderr, "TOK_SEP_TRACE missing storage ftype=%d\n", inp->ftype);
    abort();
  }
  double R[8] = { 0.0 }, dRdZ[8] = { 0.0 };
  double dR[8] = { 0.0 }, dZ[8] = { 0.0 };
  int nr = gkyl_tok_geo_R_psiZ(arc_ctx->geo, arc_ctx->geo->psisep,
    zfixed, 8, R, dRdZ, dR, dZ);
  if (nr == 0) {
    fprintf(stderr,
      "TOK_SEP_TRACE no fixed-edge root ftype=%d z=%.17g psi=%.17g\n",
      inp->ftype, zfixed, arc_ctx->geo->psisep);
    abort();
  }
  double rfixed = tok_nearest_value(rfixed_ref, R, nr);
  double rx = arc_ctx->geo->use_cubics
    ? arc_ctx->geo->efit->Rxpt_cubic[0] : arc_ctx->geo->efit->Rxpt[0];
  double zx = arc_ctx->geo->use_cubics
    ? arc_ctx->geo->efit->Zxpt_cubic[0] : arc_ctx->geo->efit->Zxpt[0];
  // A 257-point trace resolves a quadratic EFIT cell much more finely than
  // the geometry grid while avoiding the large startup penalty of rebuilding
  // 1025-point candidate traces in each geometry pass.
  int n = GKYL_MIN2(257, arc_ctx->sep_trace_capacity);
  double *rr = gkyl_malloc(sizeof(double[n]));
  double *zr = gkyl_malloc(sizeof(double[n]));
  double *rz = gkyl_malloc(sizeof(double[n]));
  double *zz = gkyl_malloc(sizeof(double[n]));
  double score_r = DBL_MAX, score_z = DBL_MAX;
  bool ok_r = tok_build_contour_candidate(arc_ctx->geo,
    arc_ctx->geo->psisep, true, rfixed, zfixed, rx, zx, n,
    rr, zr, &score_r, "separatrix", false);
  bool ok_z = tok_build_contour_candidate(arc_ctx->geo,
    arc_ctx->geo->psisep, false, rfixed, zfixed, rx, zx, n,
    rz, zz, &score_z, "separatrix", false);
  if (!ok_r && !ok_z) {
    fprintf(stderr,
      "TOK_SEP_TRACE both parameterizations failed ftype=%d fixed=(%.17g,%.17g) xpt=(%.17g,%.17g)\n",
      inp->ftype, rfixed, zfixed, rx, zx);
    abort();
  }
  bool use_r = ok_r && (!ok_z || score_r <= score_z);
  const double *src_r = use_r ? rr : rz;
  const double *src_z = use_r ? zr : zz;
  bool fixed_first = tok_sep_fixed_edge_is_first(inp->ftype);
  for (int i=0; i<n; ++i) {
    int src = fixed_first ? i : n-1-i;
    arc_ctx->sep_trace_r[i] = src_r[src];
    arc_ctx->sep_trace_z[i] = src_z[src];
  }
  arc_ctx->sep_trace_s[0] = 0.0;
  for (int i=1; i<n; ++i)
    arc_ctx->sep_trace_s[i] = arc_ctx->sep_trace_s[i-1]
      +hypot(arc_ctx->sep_trace_r[i]-arc_ctx->sep_trace_r[i-1],
        arc_ctx->sep_trace_z[i]-arc_ctx->sep_trace_z[i-1]);
  if (!(arc_ctx->sep_trace_s[n-1] > 0.0) ||
      !isfinite(arc_ctx->sep_trace_s[n-1])) {
    fprintf(stderr, "TOK_SEP_TRACE invalid arc length ftype=%d\n", inp->ftype);
    abort();
  }
  arc_ctx->sep_trace_n = n;
  arc_ctx->sep_trace_param_is_r = use_r;
  arc_ctx->sep_trace_initialized = true;
  gkyl_free(rr); gkyl_free(zr); gkyl_free(rz); gkyl_free(zz);
}

static bool
tok_half_domain_sep_rz(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double theta, double *r, double *z)
{
  if (!inp->half_domain ||
      !tok_geo_same_flux(arc_ctx->psi, arc_ctx->geo->psisep))
    return false;

  double zx = arc_ctx->geo->use_cubics
    ? arc_ctx->geo->efit->Zxpt_cubic[0] : arc_ctx->geo->efit->Zxpt[0];
  double za = arc_ctx->geo->zmaxis, z0 = 0.0, z1 = 0.0, rclose = 0.0;
  switch (inp->ftype) {
    case GKYL_GEOMETRY_TOKAMAK_CORE_R:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID:
      z0 = zx; z1 = za; rclose = inp->rright;
      break;
    case GKYL_GEOMETRY_TOKAMAK_CORE_L:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID:
      z0 = za; z1 = zx; rclose = inp->rleft;
      break;
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_R:
      z0 = arc_ctx->zmin_right; z1 = zx; rclose = inp->rright;
      break;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO:
      z0 = arc_ctx->zmin; z1 = zx; rclose = inp->rright;
      break;
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_L:
      z0 = zx; z1 = arc_ctx->zmin_left; rclose = inp->rleft;
      break;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO:
      z0 = zx; z1 = arc_ctx->zmin; rclose = inp->rleft;
      break;
    default:
      return false;
  }

  // Use the actual (position-mapped) computational coordinate.  Corner
  // nodes lie at frac=0 or 1, while interior and radial/alpha-surface
  // quadrature nodes do not.  Inferring frac from the nodal-array index
  // incorrectly maps the first and last Gauss points onto the segment ends.
  double frac = (theta-inp->cgrid.lower[2])
    /(inp->cgrid.upper[2]-inp->cgrid.lower[2]);
  double zfixed = tok_sep_fixed_edge_is_first(inp->ftype) ? z0 : z1;
  tok_build_sep_trace(inp, arc_ctx, zfixed, rclose);
  int n = arc_ctx->sep_trace_n;
  double total = arc_ctx->sep_trace_s[n-1];
  arc_ctx->xpt_map_darc_dtheta = total/
    (inp->cgrid.upper[2]-inp->cgrid.lower[2]);
  if (frac <= 0.0) {
    *r = arc_ctx->sep_trace_r[0]; *z = arc_ctx->sep_trace_z[0];
    return true;
  }
  if (frac >= 1.0) {
    *r = arc_ctx->sep_trace_r[n-1]; *z = arc_ctx->sep_trace_z[n-1];
    return true;
  }
  double target = frac*total;
  int lo = 0, hi = n-1;
  while (hi-lo > 1) {
    int mid = (lo+hi)/2;
    if (arc_ctx->sep_trace_s[mid] < target) lo = mid;
    else hi = mid;
  }
  double ds = arc_ctx->sep_trace_s[hi]-arc_ctx->sep_trace_s[lo];
  double w = ds > 0.0 ? (target-arc_ctx->sep_trace_s[lo])/ds : 0.0;
  double rlin = arc_ctx->sep_trace_r[lo]
    +w*(arc_ctx->sep_trace_r[hi]-arc_ctx->sep_trace_r[lo]);
  double zlin = arc_ctx->sep_trace_z[lo]
    +w*(arc_ctx->sep_trace_z[hi]-arc_ctx->sep_trace_z[lo]);
  if (arc_ctx->sep_trace_param_is_r) {
    double roots[16] = { 0.0 };
    int nr = tok_geo_Z_psiR(arc_ctx->geo, arc_ctx->geo->psisep,
      rlin, 16, roots);
    if (nr == 0) {
      fprintf(stderr,
        "TOK_SEP_TRACE no Z root during resampling ftype=%d frac=%.17g R=%.17g\n",
        inp->ftype, frac, rlin);
      abort();
    }
    *r = rlin; *z = tok_nearest_value(zlin, roots, nr);
  }
  else {
    double roots[8] = { 0.0 }, dRdZ[8] = { 0.0 };
    double dR[8] = { 0.0 }, dZ[8] = { 0.0 };
    int nr = gkyl_tok_geo_R_psiZ(arc_ctx->geo, arc_ctx->geo->psisep,
      zlin, 8, roots, dRdZ, dR, dZ);
    if (nr == 0) {
      fprintf(stderr,
        "TOK_SEP_TRACE no R root during resampling ftype=%d frac=%.17g Z=%.17g\n",
        inp->ftype, frac, zlin);
      abort();
    }
    *r = tok_nearest_value(rlin, roots, nr); *z = zlin;
  }
  double residual = tok_eval_psi_rz_local(arc_ctx->geo, *r, *z)
    -arc_ctx->geo->psisep;
  if (!isfinite(*r) || !isfinite(*z) || !isfinite(residual) ||
      fabs(residual) > 1e-9*fmax(1.0, fabs(arc_ctx->geo->psisep))) {
    fprintf(stderr,
      "TOK_SEP_TRACE invalid resampled point ftype=%d frac=%.17g R=%.17g Z=%.17g residual=%.17g\n",
      inp->ftype, frac, *r, *z, residual);
    abort();
  }
  return true;
}

struct tok_ordered_point {
  double r, z, phi;
  double dr_dtheta, dz_dtheta, dphi_dtheta;
};

static bool
tok_fixed_edge_is_midplane(enum gkyl_tok_geo_type ftype)
{
  return ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R ||
    ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L ||
    ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID ||
    ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID;
}

static bool
tok_fixed_edge_uses_lower_plate(enum gkyl_tok_geo_type ftype)
{
  return ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R ||
    ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO ||
    ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO;
}

// Locate the physical (plate or midplane) end of a block at a specified psi.
// Plate intersections are scanned before refinement so an exact endpoint is
// returned directly and an ambiguous multi-intersection plate is rejected.
static bool
tok_fixed_edge_point(const struct gkyl_tok_geo_grid_inp *inp,
  const struct gkyl_tok_geo *geo, double psi, double *r, double *z)
{
  if (tok_fixed_edge_is_midplane(inp->ftype)) {
    double R[16] = { 0.0 }, dRdZ[16] = { 0.0 };
    double dR[16] = { 0.0 }, dZ[16] = { 0.0 };
    *z = geo->zmaxis;
    int nr = gkyl_tok_geo_R_psiZ(geo, psi, *z, 8, R, dRdZ, dR, dZ);
    if (nr <= 0)
      return false;
    bool outboard = inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R ||
      inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID;
    *r = tok_nearest_value(outboard ? inp->rright : inp->rleft, R, nr);
    return isfinite(*r);
  }

  plate_func plate = 0;
  if (tok_fixed_edge_uses_lower_plate(inp->ftype))
    plate = geo->plate_func_lower;
  else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_L)
    plate = geo->plate_func_upper;
  if (!plate)
    return false;

  const int nsamp = 512;
  double roots_s[16] = { 0.0 }, roots_r[16] = { 0.0 };
  double roots_z[16] = { 0.0 };
  int nroots = 0;
  double rz0[2] = { 0.0 }, rz1[2] = { 0.0 };
  plate(0.0, rz0);
  double f0 = tok_eval_psi_rz_local(geo, rz0[0], rz0[1])-psi;
  double flux_tol = 1e-10*fmax(1.0, fabs(psi));
  if (isfinite(f0) && fabs(f0) <= flux_tol) {
    roots_s[nroots] = 0.0; roots_r[nroots] = rz0[0];
    roots_z[nroots++] = rz0[1];
  }
  for (int i=1; i<=nsamp; ++i) {
    double s1 = i/(double) nsamp;
    plate(s1, rz1);
    double f1 = tok_eval_psi_rz_local(geo, rz1[0], rz1[1])-psi;
    bool endpoint_root = isfinite(f1) && fabs(f1) <= flux_tol;
    bool bracket = isfinite(f0) && isfinite(f1) && f0*f1 < 0.0;
    if ((endpoint_root || bracket) && nroots < 16) {
      double sr = s1;
      if (!endpoint_root) {
        double slo = (i-1)/(double) nsamp, shi = s1, flo = f0;
        for (int k=0; k<64; ++k) {
          double smid = 0.5*(slo+shi), rzm[2];
          plate(smid, rzm);
          double fm = tok_eval_psi_rz_local(geo, rzm[0], rzm[1])-psi;
          if (!isfinite(fm))
            return false;
          if (flo*fm <= 0.0)
            shi = smid;
          else { slo = smid; flo = fm; }
        }
        sr = 0.5*(slo+shi);
      }
      double rzr[2];
      plate(sr, rzr);
      if (nroots == 0 || fabs(sr-roots_s[nroots-1]) > 2e-10) {
        roots_s[nroots] = sr; roots_r[nroots] = rzr[0];
        roots_z[nroots++] = rzr[1];
      }
    }
    rz0[0] = rz1[0]; rz0[1] = rz1[1]; f0 = f1;
  }
  if (nroots != 1) {
    fprintf(stderr,
      "TOK_ORDERED_MAP plate root count=%d ftype=%d psi=%.17g\n",
      nroots, inp->ftype, psi);
    return false;
  }
  *r = roots_r[0]; *z = roots_z[0];
  double residual = tok_eval_psi_rz_local(geo, *r, *z)-psi;
  return isfinite(*r) && isfinite(*z) && isfinite(residual) &&
    fabs(residual) <= 1e-9*fmax(1.0, fabs(psi));
}

static bool
tok_trace_sample(const struct gkyl_tok_geo *geo, double psi,
  const double *tr, const double *tz, const double *ts, int n,
  bool param_is_r, double u, double *r, double *z)
{
  if (n < 2 || !(ts[n-1] > 0.0))
    return false;
  // Logical edge coordinates acquire a few ulps of roundoff when reconstructed
  // from the arc-length coordinate.  Snap those values before resampling: near
  // an X point, asking the polynomial root finder for an infinitesimally
  // interior point can select the other root of the saddle instead of the
  // explicitly stored, shared endpoint.
  const double endpoint_tol = 256.0*DBL_EPSILON;
  if (u <= endpoint_tol) { *r = tr[0]; *z = tz[0]; return true; }
  if (u >= 1.0-endpoint_tol) {
    *r = tr[n-1]; *z = tz[n-1]; return true;
  }
  double target = u*ts[n-1];
  int lo = 0, hi = n-1;
  while (hi-lo > 1) {
    int mid = (lo+hi)/2;
    if (ts[mid] < target) lo = mid;
    else hi = mid;
  }
  double ds = ts[hi]-ts[lo];
  double w = ds > 0.0 ? (target-ts[lo])/ds : 0.0;
  double rlin = tr[lo]+w*(tr[hi]-tr[lo]);
  double zlin = tz[lo]+w*(tz[hi]-tz[lo]);
  if (param_is_r) {
    double roots[32] = { 0.0 };
    int nr = tok_geo_Z_psiR(geo, psi, rlin, 16, roots);
    if (nr <= 0) return false;
    *r = rlin; *z = tok_nearest_value(zlin, roots, nr);
  }
  else {
    double roots[16] = { 0.0 }, dRdZ[16] = { 0.0 };
    double dR[16] = { 0.0 }, dZ[16] = { 0.0 };
    int nr = gkyl_tok_geo_R_psiZ(geo, psi, zlin, 8,
      roots, dRdZ, dR, dZ);
    if (nr <= 0) return false;
    *r = tok_nearest_value(rlin, roots, nr); *z = zlin;
  }
  double residual = tok_eval_psi_rz_local(geo, *r, *z)-psi;
  return isfinite(*r) && isfinite(*z) && isfinite(residual) &&
    fabs(residual) <= 1e-9*fmax(1.0, fabs(psi));
}

static bool
tok_build_far_trace(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx)
{
  if (arc_ctx->far_trace_initialized)
    return true;
  if (!arc_ctx->xpt_ray_initialized)
    return false;
  double rf = 0.0, zf = 0.0;
  if (!tok_fixed_edge_point(inp, arc_ctx->geo, arc_ctx->xpt_ray_psi0,
      &rf, &zf))
    return false;
  int n = GKYL_MIN2(257, arc_ctx->sep_trace_capacity);
  double *rr = gkyl_malloc(sizeof(double[n]));
  double *zr = gkyl_malloc(sizeof(double[n]));
  double *rz = gkyl_malloc(sizeof(double[n]));
  double *zz = gkyl_malloc(sizeof(double[n]));
  double score_r = DBL_MAX, score_z = DBL_MAX;
  bool ok_r = tok_build_contour_candidate(arc_ctx->geo,
    arc_ctx->xpt_ray_psi0, true, rf, zf, arc_ctx->xpt_ray_r0,
    arc_ctx->xpt_ray_z0, n, rr, zr, &score_r, "far_boundary", false);
  bool ok_z = tok_build_contour_candidate(arc_ctx->geo,
    arc_ctx->xpt_ray_psi0, false, rf, zf, arc_ctx->xpt_ray_r0,
    arc_ctx->xpt_ray_z0, n, rz, zz, &score_z, "far_boundary", false);
  // Preserve both original parameterizations whenever either works.  Only a
  // trace for which both were rejected and whose fixed edge is a known
  // midplane turning point receives the endpoint-clustered R retry.
  if (!ok_r && !ok_z && tok_fixed_edge_is_midplane(inp->ftype))
    ok_r = tok_build_contour_candidate(arc_ctx->geo,
      arc_ctx->xpt_ray_psi0, true, rf, zf, arc_ctx->xpt_ray_r0,
      arc_ctx->xpt_ray_z0, n, rr, zr, &score_r, "far_boundary", true);
  if (!ok_r && !ok_z) {
    fprintf(stderr,
      "TOK_ORDERED_MAP far trace failed ftype=%d psi=%.17g fixed=(%.17g,%.17g) ray=(%.17g,%.17g)\n",
      inp->ftype, arc_ctx->xpt_ray_psi0, rf, zf, arc_ctx->xpt_ray_r0,
      arc_ctx->xpt_ray_z0);
    gkyl_free(rr); gkyl_free(zr); gkyl_free(rz); gkyl_free(zz);
    return false;
  }
  bool use_r = ok_r && (!ok_z || score_r <= score_z);
  const double *src_r = use_r ? rr : rz;
  const double *src_z = use_r ? zr : zz;
  bool fixed_first = tok_sep_fixed_edge_is_first(inp->ftype);
  arc_ctx->far_trace_s[0] = 0.0;
  for (int i=0; i<n; ++i) {
    int src = fixed_first ? i : n-1-i;
    arc_ctx->far_trace_r[i] = src_r[src];
    arc_ctx->far_trace_z[i] = src_z[src];
    if (i > 0)
      arc_ctx->far_trace_s[i] = arc_ctx->far_trace_s[i-1]
        +hypot(arc_ctx->far_trace_r[i]-arc_ctx->far_trace_r[i-1],
          arc_ctx->far_trace_z[i]-arc_ctx->far_trace_z[i-1]);
  }
  arc_ctx->far_trace_n = n;
  arc_ctx->far_trace_param_is_r = use_r;
  arc_ctx->far_trace_initialized =
    isfinite(arc_ctx->far_trace_s[n-1]) && arc_ctx->far_trace_s[n-1] > 0.0;
  gkyl_free(rr); gkyl_free(zr); gkyl_free(rz); gkyl_free(zz);
  return arc_ctx->far_trace_initialized;
}

static bool
tok_build_trace_correspondence(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx)
{
  if (arc_ctx->ordered_boundaries_initialized)
    return true;
  if (!arc_ctx->sep_trace_initialized) {
    double rf = 0.0, zf = 0.0;
    if (!tok_fixed_edge_point(inp, arc_ctx->geo, arc_ctx->geo->psisep,
        &rf, &zf))
      return false;
    tok_build_sep_trace(inp, arc_ctx, zf, rf);
  }
  if (!tok_build_far_trace(inp, arc_ctx))
    return false;

  int n = GKYL_MIN2(arc_ctx->sep_trace_n, arc_ctx->far_trace_n);
  arc_ctx->trace_corr_v[0] = 0.0;
  for (int i=1; i<n-1; ++i) {
    double u = i/(double) (n-1), rs = 0.0, zs = 0.0;
    if (!tok_trace_sample(arc_ctx->geo, arc_ctx->geo->psisep,
        arc_ctx->sep_trace_r, arc_ctx->sep_trace_z,
        arc_ctx->sep_trace_s, arc_ctx->sep_trace_n,
        arc_ctx->sep_trace_param_is_r, u, &rs, &zs))
      return false;
    double best_d2 = DBL_MAX, best_v = 0.0;
    for (int j=0; j<arc_ctx->far_trace_n-1; ++j) {
      double r0 = arc_ctx->far_trace_r[j], z0 = arc_ctx->far_trace_z[j];
      double dr = arc_ctx->far_trace_r[j+1]-r0;
      double dz = arc_ctx->far_trace_z[j+1]-z0;
      double den = dr*dr+dz*dz;
      double t = den > 0.0 ? ((rs-r0)*dr+(zs-z0)*dz)/den : 0.0;
      t = fmin(1.0, fmax(0.0, t));
      double rp = r0+t*dr, zp = z0+t*dz;
      double d2 = SQ(rs-rp)+SQ(zs-zp);
      if (d2 < best_d2) {
        best_d2 = d2;
        best_v = (arc_ctx->far_trace_s[j]
          +t*(arc_ctx->far_trace_s[j+1]-arc_ctx->far_trace_s[j]))
          /arc_ctx->far_trace_s[arc_ctx->far_trace_n-1];
      }
    }
    arc_ctx->trace_corr_v[i] = best_v;
  }
  arc_ctx->trace_corr_v[n-1] = 1.0;

  // Isotonic regression (increasing PAVA) finds the least-squares monotone
  // correspondence, rather than greedily flattening every value after the
  // first local reversal.  Endpoints are already bounded by every projection,
  // so applying PAVA to the interior preserves the exact 0 and 1 anchors.
  int nblock = 0;
  double *block_mean = gkyl_malloc(sizeof(double[n]));
  int *block_count = gkyl_malloc(sizeof(int[n]));
  for (int i=1; i<n-1; ++i) {
    block_mean[nblock] = arc_ctx->trace_corr_v[i];
    block_count[nblock] = 1;
    ++nblock;
    while (nblock > 1 && block_mean[nblock-2] > block_mean[nblock-1]) {
      int count = block_count[nblock-2]+block_count[nblock-1];
      block_mean[nblock-2] =
        (block_count[nblock-2]*block_mean[nblock-2]
          +block_count[nblock-1]*block_mean[nblock-1])/count;
      block_count[nblock-2] = count;
      --nblock;
    }
  }
  int out = 1;
  for (int b=0; b<nblock; ++b)
    for (int j=0; j<block_count[b]; ++j)
      arc_ctx->trace_corr_v[out++] = block_mean[b];
  gkyl_free(block_mean);
  gkyl_free(block_count);
  // Nearest-point projection can legitimately collapse a finite interval of
  // the separatrix onto one far-boundary location.  Blend in a small identity
  // map so g remains well-conditioned and strictly order preserving without
  // discarding the nearest-normal correspondence.  Shot 203730 remains
  // fold-free over a broad range around this conservative one-percent blend.
  const double identity_fraction = 0.01;
  for (int i=1; i<n-1; ++i) {
    double u = i/(double) (n-1);
    arc_ctx->trace_corr_v[i] = identity_fraction*u
      +(1.0-identity_fraction)*arc_ctx->trace_corr_v[i];
  }
  for (int i=1; i<n; ++i)
    if (!(arc_ctx->trace_corr_v[i] > arc_ctx->trace_corr_v[i-1]))
      return false;
  arc_ctx->trace_corr_n = n;
  arc_ctx->ordered_boundaries_initialized = true;
  return true;
}

static double
tok_trace_correspondence(const struct arc_length_ctx *arc_ctx, double u)
{
  const double endpoint_tol = 256.0*DBL_EPSILON;
  if (u <= endpoint_tol) return 0.0;
  if (u >= 1.0-endpoint_tol) return 1.0;
  double x = u*(arc_ctx->trace_corr_n-1);
  int i = GKYL_MIN2(arc_ctx->trace_corr_n-2, (int) floor(x));
  double w = x-i;
  return arc_ctx->trace_corr_v[i]
    +w*(arc_ctx->trace_corr_v[i+1]-arc_ctx->trace_corr_v[i]);
}

static bool
tok_ordered_chord_point(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double psi, double u, double *r, double *z)
{
  if (!tok_build_trace_correspondence(inp, arc_ctx))
    return false;
  u = fmin(1.0, fmax(0.0, u));
  double v = tok_trace_correspondence(arc_ctx, u);
  double ra = 0.0, za = 0.0, rb = 0.0, zb = 0.0;
  if (!tok_trace_sample(arc_ctx->geo, arc_ctx->geo->psisep,
      arc_ctx->sep_trace_r, arc_ctx->sep_trace_z,
      arc_ctx->sep_trace_s, arc_ctx->sep_trace_n,
      arc_ctx->sep_trace_param_is_r, u, &ra, &za) ||
      !tok_trace_sample(arc_ctx->geo, arc_ctx->xpt_ray_psi0,
        arc_ctx->far_trace_r, arc_ctx->far_trace_z,
        arc_ctx->far_trace_s, arc_ctx->far_trace_n,
        arc_ctx->far_trace_param_is_r, v, &rb, &zb))
    return false;
  double scale = fmax(1.0, fmax(fabs(psi), fabs(arc_ctx->geo->psisep)));
  if (tok_geo_same_flux(psi, arc_ctx->geo->psisep)) {
    *r = ra; *z = za; return true;
  }
  if (tok_geo_same_flux(psi, arc_ctx->xpt_ray_psi0)) {
    *r = rb; *z = zb; return true;
  }
  double fa = tok_eval_psi_rz_local(arc_ctx->geo, ra, za)-psi;
  double fb = tok_eval_psi_rz_local(arc_ctx->geo, rb, zb)-psi;
  if (!isfinite(fa) || !isfinite(fb) || fa*fb > 0.0) {
    fprintf(stderr,
      "TOK_ORDERED_MAP unbracketed chord ftype=%d psi=%.17g u=%.17g fa=%.17g fb=%.17g\n",
      inp->ftype, psi, u, fa, fb);
    return false;
  }
  double slo = 0.0, shi = 1.0, flo = fa;
  for (int k=0; k<44; ++k) {
    double smid = 0.5*(slo+shi);
    double rm = ra+smid*(rb-ra), zm = za+smid*(zb-za);
    double fm = tok_eval_psi_rz_local(arc_ctx->geo, rm, zm)-psi;
    if (!isfinite(fm)) return false;
    if (flo*fm <= 0.0)
      shi = smid;
    else { slo = smid; flo = fm; }
  }
  double s = 0.5*(slo+shi);
  *r = ra+s*(rb-ra); *z = za+s*(zb-za);
  double residual = tok_eval_psi_rz_local(arc_ctx->geo, *r, *z)-psi;
  return isfinite(*r) && isfinite(*z) && isfinite(residual) &&
    fabs(residual) <= 1e-9*scale;
}

static bool
tok_xpt_seam_endpoint(enum gkyl_tok_geo_type ftype, double u)
{
  const double endpoint_tol = 256.0*DBL_EPSILON;
  return tok_sep_fixed_edge_is_first(ftype)
    ? u >= 1.0-endpoint_tol : u <= endpoint_tol;
}

static bool
tok_eval_psi_grad_rz_local(const struct gkyl_tok_geo *geo,
  double R, double Z, double *dpsidR, double *dpsidZ)
{
  if (geo->use_cubics) {
    double xn[2] = { R, Z }, out[3] = { 0.0 };
    geo->efit->evf->eval_cubic_wgrad(0.0, xn, out, geo->efit->evf->ctx);
    *dpsidR = out[1]; *dpsidZ = out[2];
    return isfinite(*dpsidR) && isfinite(*dpsidZ);
  }
  int idx[2];
  idx[0] = GKYL_MIN2(geo->rzlocal.upper[0], GKYL_MAX2(geo->rzlocal.lower[0],
    geo->rzlocal.lower[0]
      +(int) floor((R-geo->rzgrid.lower[0])/geo->rzgrid.dx[0])));
  idx[1] = GKYL_MIN2(geo->rzlocal.upper[1], GKYL_MAX2(geo->rzlocal.lower[1],
    geo->rzlocal.lower[1]
      +(int) floor((Z-geo->rzgrid.lower[1])/geo->rzgrid.dx[1])));
  long loc = gkyl_range_idx(&geo->rzlocal, idx);
  const double *p = gkyl_array_cfetch(geo->psiRZ, loc);
  double xc[2];
  gkyl_rect_grid_cell_center(&geo->rzgrid, idx, xc);
  double x = (R-xc[0])/(0.5*geo->rzgrid.dx[0]);
  double y = (Z-xc[1])/(0.5*geo->rzgrid.dx[1]);
  if (geo->efit->rzbasis.poly_order == 1) {
    *dpsidR = (1.5*p[3]*y+0.8660254037844386*p[1])
      *2.0/geo->rzgrid.dx[0];
    *dpsidZ = (1.5*p[3]*x+0.8660254037844386*p[2])
      *2.0/geo->rzgrid.dx[1];
  }
  else {
    *dpsidR = (5.625*p[8]*(2.0*x*y*y-0.6666666666666666*x)
      +2.904737509655563*p[7]*(y*y-0.3333333333333333)
      +5.809475019311126*p[6]*x*y+1.5*p[3]*y
      +3.354101966249684*p[4]*x+0.8660254037844386*p[1])
      *2.0/geo->rzgrid.dx[0];
    *dpsidZ = (5.625*p[8]*(2.0*x*x*y-0.6666666666666666*y)
      +5.809475019311126*p[7]*x*y+3.354101966249684*p[5]*y
      +2.904737509655563*p[6]*(x*x-0.3333333333333333)
      +1.5*p[3]*x+0.8660254037844386*p[2])
      *2.0/geo->rzgrid.dx[1];
  }
  return isfinite(*dpsidR) && isfinite(*dpsidZ);
}

// Project a continuation predictor to the exact requested DG/cubic flux
// surface. Holding the coordinate that changes most along the contour keeps
// the root solve well conditioned at either an R or Z turning point.
static bool
tok_project_xpt_seam_candidate(const struct gkyl_tok_geo *geo, double psi,
  double rprev, double zprev, double rpred, double zpred,
  double tangent_r, double tangent_z, double step,
  double *r, double *z)
{
  const struct gkyl_rect_grid *grid = geo->use_cubics
    ? &geo->rzgrid_cubic : &geo->rzgrid;
  if (rpred < grid->lower[0] || rpred > grid->upper[0] ||
      zpred < grid->lower[1] || zpred > grid->upper[1])
    return false;

  double best_distance = DBL_MAX, second_distance = DBL_MAX;
  if (fabs(tangent_r) >= fabs(tangent_z)) {
    double roots[32] = { 0.0 };
    int nr = tok_geo_Z_psiR(geo, psi, rpred, 32, roots);
    if (nr <= 0)
      return false;
    int best = 0;
    for (int i=0; i<nr; ++i) {
      double distance = fabs(roots[i]-zpred);
      if (distance < best_distance) {
        second_distance = best_distance;
        best_distance = distance;
        best = i;
      }
      else if (distance < second_distance)
        second_distance = distance;
    }
    double ambiguity_scale = fmax(grid->dx[1], fabs(step));
    if (nr > 1 && second_distance-best_distance <=
        1e-8*fmax(1.0, ambiguity_scale))
      return false;
    *r = rpred; *z = roots[best];
  }
  else {
    double roots[32] = { 0.0 }, dRdZ[32] = { 0.0 };
    double dR[32] = { 0.0 }, dZ[32] = { 0.0 };
    int nr = gkyl_tok_geo_R_psiZ(geo, psi, zpred, 32,
      roots, dRdZ, dR, dZ);
    if (nr <= 0)
      return false;
    int best = 0;
    for (int i=0; i<nr; ++i) {
      double distance = fabs(roots[i]-rpred);
      if (distance < best_distance) {
        second_distance = best_distance;
        best_distance = distance;
        best = i;
      }
      else if (distance < second_distance)
        second_distance = distance;
    }
    double ambiguity_scale = fmax(grid->dx[0], fabs(step));
    if (nr > 1 && second_distance-best_distance <=
        1e-8*fmax(1.0, ambiguity_scale))
      return false;
    *r = roots[best]; *z = zpred;
  }

  double dr = *r-rprev, dz = *z-zprev;
  double forward = step*(dr*tangent_r+dz*tangent_z);
  double distance = hypot(dr, dz);
  double residual = tok_eval_psi_rz_local(geo, *r, *z)-psi;
  double flux_scale = fmax(1.0, fabs(psi));
  return isfinite(*r) && isfinite(*z) && isfinite(distance) &&
    isfinite(residual) && forward > 0.0 &&
    distance <= 2.5*fabs(step)+1e-12 &&
    fabs(residual) <= 1e-10*flux_scale;
}

static bool
tok_displace_xpt_seam_on_flux(const struct gkyl_tok_geo *geo, double psi,
  double delta_s, double *r, double *z)
{
  if (delta_s == 0.0)
    return true;
  const struct gkyl_rect_grid *grid = geo->use_cubics
    ? &geo->rzgrid_cubic : &geo->rzgrid;
  double max_step = 0.05*fmin(grid->dx[0], grid->dx[1]);
  if (!(max_step > 0.0) || !isfinite(max_step))
    return false;
  double required_steps = ceil(fabs(delta_s)/max_step);
  if (!isfinite(required_steps) || required_steps > 4096.0)
    return false;
  int nstep = GKYL_MAX2(16, (int) required_steps);
  double rstart = *r, zstart = *z;
  double relative_arc_error = DBL_MAX, realized = 0.0;
  for (;;) {
    *r = rstart; *z = zstart;
    double step = delta_s/nstep;
    realized = 0.0;

    for (int i=0; i<nstep; ++i) {
      double grad_r = 0.0, grad_z = 0.0;
      if (!tok_eval_psi_grad_rz_local(geo, *r, *z, &grad_r, &grad_z))
        return false;
      double grad = hypot(grad_r, grad_z);
      if (!(grad > 1e-14))
        return false;
      double tangent_r = -grad_z/grad, tangent_z = grad_r/grad;

      double rmid = *r+0.5*step*tangent_r;
      double zmid = *z+0.5*step*tangent_z;
      if (!tok_project_xpt_seam_candidate(geo, psi, *r, *z, rmid, zmid,
          tangent_r, tangent_z, 0.5*step, &rmid, &zmid)) {
        if (tok_ordered_map_diag_enabled())
          fprintf(stderr,
            "TOK_XPT_SEAM_DIAG reason=midpoint_projection_failed psi=%.17g step=%d/%d R=%.17g Z=%.17g delta_s=%.17g\n",
            psi, i, nstep, *r, *z, delta_s);
        return false;
      }
      if (!tok_eval_psi_grad_rz_local(geo, rmid, zmid, &grad_r, &grad_z))
        return false;
      grad = hypot(grad_r, grad_z);
      if (!(grad > 1e-14))
        return false;
      tangent_r = -grad_z/grad; tangent_z = grad_r/grad;

      double rnext = *r+step*tangent_r;
      double znext = *z+step*tangent_z;
      if (!tok_project_xpt_seam_candidate(geo, psi, *r, *z, rnext, znext,
          tangent_r, tangent_z, step, &rnext, &znext)) {
        if (tok_ordered_map_diag_enabled())
          fprintf(stderr,
            "TOK_XPT_SEAM_DIAG reason=full_projection_failed psi=%.17g step=%d/%d R=%.17g Z=%.17g delta_s=%.17g\n",
            psi, i, nstep, *r, *z, delta_s);
        return false;
      }
      realized += hypot(rnext-*r, znext-*z);
      *r = rnext; *z = znext;
    }
    relative_arc_error = fabs(realized-fabs(delta_s))
      /fmax(fabs(delta_s), 1e-14);
    if (isfinite(relative_arc_error) && relative_arc_error <= 5e-3)
      return true;
    if (tok_ordered_map_diag_enabled())
      fprintf(stderr,
        "TOK_XPT_SEAM_DIAG reason=arc_length_refine psi=%.17g requested=%.17g realized=%.17g signed_relative_error=%.17g nstep=%d step_size=%.17g\n",
        psi, delta_s, realized,
        (realized-fabs(delta_s))/fmax(fabs(delta_s), 1e-14), nstep,
        fabs(step));
    if (nstep > 4096/2)
      break;
    nstep *= 2;
  }
  *r = rstart; *z = zstart;
  if (tok_ordered_map_diag_enabled())
    fprintf(stderr,
      "TOK_XPT_SEAM_DIAG reason=arc_length_not_converged psi=%.17g requested=%.17g realized=%.17g relative_error=%.17g nstep=%d\n",
      psi, delta_s, realized, relative_arc_error, nstep);
  return false;
}

// B1(q)=4q(1-q) is smooth, has unit peak, and uses explicit endpoint
// branches so both physical anchors remain exactly fixed.
static bool
tok_xpt_seam_delta_s(const struct gkyl_tok_geo_grid_inp *inp,
  const struct arc_length_ctx *arc_ctx, double psi, double *q,
  double *delta_s)
{
  double span = arc_ctx->xpt_ray_psi0-arc_ctx->geo->psisep;
  double scale = fmax(1.0, fmax(fabs(arc_ctx->xpt_ray_psi0),
    fabs(arc_ctx->geo->psisep)));
  if (!isfinite(span) || fabs(span) <= 256.0*DBL_EPSILON*scale)
    return false;
  if (tok_geo_same_flux(psi, arc_ctx->geo->psisep)) {
    *q = 0.0; *delta_s = 0.0;
    return true;
  }
  if (tok_geo_same_flux(psi, arc_ctx->xpt_ray_psi0)) {
    *q = 1.0; *delta_s = 0.0;
    return true;
  }
  *q = (psi-arc_ctx->geo->psisep)/span;
  if (*q < -1e-8 || *q > 1.0+1e-8)
    return false;
  *q = fmin(1.0, fmax(0.0, *q));
  if (*q == 0.0 || *q == 1.0) {
    *delta_s = 0.0;
    return true;
  }
  double coefficient = inp->relaxed_xpt_seam_delta_s_coeff;
  double bound = inp->relaxed_xpt_seam_delta_s_bound;
  if (coefficient == 0.0) {
    *delta_s = 0.0;
    return true;
  }
  if (!isfinite(coefficient) || !isfinite(bound) || !(bound > 0.0) ||
      fabs(coefficient) > bound*(1.0+64.0*DBL_EPSILON))
    return false;
  *delta_s = coefficient*4.0*(*q)*(1.0-*q);
  return isfinite(*delta_s) && fabs(*delta_s) <=
    bound*(1.0+64.0*DBL_EPSILON);
}

static bool
tok_parameterized_xpt_seam_point(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double psi, double u, double *r, double *z)
{
  // Retain the validated straight-ray point before considering relaxation.
  if (!tok_ordered_chord_point(inp, arc_ctx, psi, u, r, z))
    return false;
  if (!inp->relaxed_xpt_seam || !inp->half_domain)
    return true;

  double q = 0.0, delta_s = 0.0;
  if (!tok_xpt_seam_delta_s(inp, arc_ctx, psi, &q, &delta_s)) {
    fprintf(stderr,
      "TOK_XPT_SEAM_DIAG mode=rejected ftype=%d psi=%.17g q=%.17g coefficient=%.17g bound=%.17g reason=invalid_delta_s action=%s\n",
      inp->ftype, psi, q, inp->relaxed_xpt_seam_delta_s_coeff,
      inp->relaxed_xpt_seam_delta_s_bound,
      inp->relaxed_xpt_seam_sweep ? "reject_candidate" : "retain_straight");
    return !inp->relaxed_xpt_seam_sweep;
  }

  // Use signed local contour arc length, with s0=0 at the straight-ray
  // intersection. Return before doing coordinate arithmetic so zero mode is
  // bitwise identical to the existing construction.
  const double s0 = 0.0;
  if (delta_s == 0.0) {
    if (tok_ordered_map_diag_enabled() &&
        tok_xpt_seam_endpoint(inp->ftype, u))
      fprintf(stderr,
        "TOK_XPT_SEAM_DIAG mode=zero ftype=%d psi=%.17g q=%.17g s0=%.17g delta_s=%.17g seam_s=%.17g R=%.17g Z=%.17g construction=straight\n",
        inp->ftype, psi, q, s0, delta_s, s0, *r, *z);
    return true;
  }

  if (!inp->relaxed_xpt_seam_sweep || !inp->straight_xpt_ray) {
    fprintf(stderr,
      "TOK_XPT_SEAM_DIAG mode=rejected ftype=%d psi=%.17g q=%.17g s0=%.17g delta_s=%.17g reason=diagnostic_sweep_not_enabled fallback=straight\n",
      inp->ftype, psi, q, s0, delta_s);
    return true;
  }

  // The only free data are on the seam. A fixed linear arc weight carries
  // that boundary condition to the unchanged opposite edge, giving a smooth
  // diagnostic shadow map without introducing interior degrees of freedom.
  double seam_weight = tok_sep_fixed_edge_is_first(inp->ftype) ? u : 1.0-u;
  double point_delta_s = seam_weight*delta_s;
  if (point_delta_s == 0.0)
    return true;
  double candidate_r = *r, candidate_z = *z;
  if (!tok_displace_xpt_seam_on_flux(arc_ctx->geo, psi, point_delta_s,
      &candidate_r, &candidate_z)) {
    fprintf(stderr,
      "TOK_XPT_SEAM_DIAG mode=rejected ftype=%d psi=%.17g q=%.17g u=%.17g delta_s=%.17g point_delta_s=%.17g reason=contour_tracking_failed action=reject_candidate fallback=none\n",
      inp->ftype, psi, q, u, delta_s, point_delta_s);
    return false;
  }
  *r = candidate_r; *z = candidate_z;
  if (tok_ordered_map_diag_enabled() &&
      tok_xpt_seam_endpoint(inp->ftype, u))
    fprintf(stderr,
      "TOK_XPT_SEAM_DIAG mode=candidate ftype=%d psi=%.17g q=%.17g s0=%.17g delta_s=%.17g seam_s=%.17g R=%.17g Z=%.17g extension=fixed_edge_linear_arc_blend\n",
      inp->ftype, psi, q, s0, delta_s, s0+delta_s, *r, *z);
  return true;
}

static double
tok_fpol_at_psi(const struct gkyl_tok_geo *geo, double psi)
{
  double p = psi;
  if (p < geo->fgrid.lower[0] || p > geo->fgrid.upper[0])
    p = geo->sibry;
  int idx = GKYL_MIN2(geo->frange.upper[0], GKYL_MAX2(geo->frange.lower[0],
    geo->frange.lower[0]
      +(int) floor((p-geo->fgrid.lower[0])/geo->fgrid.dx[0])));
  long loc = gkyl_range_idx(&geo->frange, &idx);
  const double *coeffs = gkyl_array_cfetch(geo->fpoldg, loc);
  double xc;
  gkyl_rect_grid_cell_center(&geo->fgrid, &idx, &xc);
  double x = (p-xc)/(0.5*geo->fgrid.dx[0]);
  return geo->fbasis.eval_expand(&x, coeffs);
}

static bool
tok_build_current_ordered_trace(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx)
{
  double psi = arc_ctx->psi;
  if (arc_ctx->map_trace_initialized && tok_geo_same_flux(
      psi, arc_ctx->map_trace_psi))
    return true;
  if (!tok_build_trace_correspondence(inp, arc_ctx))
    return false;
  int n = GKYL_MIN2(arc_ctx->sep_trace_capacity,
    GKYL_MAX2(49, 4*inp->cgrid.cells[2]+1));
  double fpol = tok_fpol_at_psi(arc_ctx->geo, psi);
  arc_ctx->map_trace_s[0] = 0.0;
  arc_ctx->map_trace_phi[0] = 0.0;
  for (int i=0; i<n; ++i) {
    double u = i/(double) (n-1);
    if (!tok_parameterized_xpt_seam_point(inp, arc_ctx, psi, u,
        &arc_ctx->map_trace_r[i], &arc_ctx->map_trace_z[i]))
      return false;
    if (i > 0) {
      double ds = hypot(arc_ctx->map_trace_r[i]-arc_ctx->map_trace_r[i-1],
        arc_ctx->map_trace_z[i]-arc_ctx->map_trace_z[i-1]);
      if (!(ds > 0.0) || !isfinite(ds)) return false;
      arc_ctx->map_trace_s[i] = arc_ctx->map_trace_s[i-1]+ds;
      double rm = 0.5*(arc_ctx->map_trace_r[i]
        +arc_ctx->map_trace_r[i-1]);
      double zm = 0.5*(arc_ctx->map_trace_z[i]
        +arc_ctx->map_trace_z[i-1]);
      double gr = 0.0, gz = 0.0;
      if (!tok_eval_psi_grad_rz_local(arc_ctx->geo, rm, zm, &gr, &gz))
        return false;
      double grad = hypot(gr, gz);
      if (!(grad > 1e-14) || !(rm > 0.0)) return false;
      arc_ctx->map_trace_phi[i] = arc_ctx->map_trace_phi[i-1]
        +fpol*ds/(rm*grad);
    }
  }
  arc_ctx->map_trace_n = n;
  arc_ctx->map_trace_psi = psi;
  arc_ctx->map_trace_initialized = true;
  return true;
}

static bool
tok_ordered_map_lookup(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double theta, double alpha,
  struct tok_ordered_point *out)
{
  if (!tok_xpt_mapping_requested(inp))
    return false;
  if (!tok_build_current_ordered_trace(inp, arc_ctx)) {
    fprintf(stderr,
      "TOK_ORDERED_MAP initialization failed ftype=%d psi=%.17g\n",
      inp->ftype, arc_ctx->psi);
    abort();
  }
  double dtheta = inp->cgrid.upper[2]-inp->cgrid.lower[2];
  double u = (theta-inp->cgrid.lower[2])/dtheta;
  if (u < -1e-10 || u > 1.0+1e-10) return false;
  if (u <= 256.0*DBL_EPSILON) u = 0.0;
  else if (u >= 1.0-256.0*DBL_EPSILON) u = 1.0;
  else u = fmin(1.0, fmax(0.0, u));
  if (!tok_parameterized_xpt_seam_point(inp, arc_ctx, arc_ctx->psi, u,
      &out->r, &out->z))
    return false;
  double x = u*(arc_ctx->map_trace_n-1);
  int i = GKYL_MIN2(arc_ctx->map_trace_n-2, GKYL_MAX2(0, (int) floor(x)));
  double du = 1.0/(arc_ctx->map_trace_n-1);
  double dr = arc_ctx->map_trace_r[i+1]-arc_ctx->map_trace_r[i];
  double dz = arc_ctx->map_trace_z[i+1]-arc_ctx->map_trace_z[i];
  double speed_u = hypot(dr, dz)/du;
  double gr = 0.0, gz = 0.0;
  if (!tok_eval_psi_grad_rz_local(arc_ctx->geo, out->r, out->z,
      &gr, &gz))
    return false;
  double tr = dr, tz = dz, grad = hypot(gr, gz);
  if (grad > 1e-14) {
    tr = -gz/grad; tz = gr/grad;
    if (tr*dr+tz*dz < 0.0) { tr = -tr; tz = -tz; }
  }
  else {
    double tmag = hypot(tr, tz);
    if (!(tmag > 0.0)) return false;
    tr /= tmag; tz /= tmag;
  }
  out->dr_dtheta = tr*speed_u/dtheta;
  out->dz_dtheta = tz*speed_u/dtheta;
  double w = x-i;
  double path_phi = arc_ctx->map_trace_phi[i]
    +w*(arc_ctx->map_trace_phi[i+1]-arc_ctx->map_trace_phi[i]);
  double ref_phi = tok_sep_fixed_edge_is_first(inp->ftype) ? 0.0
    : arc_ctx->map_trace_phi[arc_ctx->map_trace_n-1];
  out->phi = alpha+path_phi-ref_phi;
  out->dphi_dtheta = (arc_ctx->map_trace_phi[i+1]
    -arc_ctx->map_trace_phi[i])/du/dtheta;
  return isfinite(out->r) && isfinite(out->z) && isfinite(out->phi) &&
    isfinite(out->dr_dtheta) && isfinite(out->dz_dtheta) &&
    isfinite(out->dphi_dtheta);
}

// Function to calculate phi given alpha
double
phi_func(double alpha_curr, double Z, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo = actx->arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL, zmax = actx->zmax;

  // Here we will abandon conventions about alpha and phi except for full core and full SN cases
  // The convention for phi only affects b_x - it does not affect any quantities used in axisymmetric simulations
  // I have not quite figured out full 3D yet. b_x presents a serious problem as of now. Akash Shukla 1/20/2024
  // The idea for axisymmetry is that I am avoiding starting integrals at the x-point to minimize issues
  double ival = 0;
  double phi_ref = 0.0;
  if (actx->ftype==GKYL_GEOMETRY_TOKAMAK_CORE){ // phi = alpha at outboard midplane
    if(actx->right==true){
      if(Z<actx->zmaxis)
        ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
      else
        ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
    }
    else{
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmax, rclose, false, false, arc_memo);
      phi_ref = actx->phi_right;
    }
  }
  else if (actx->ftype==GKYL_GEOMETRY_TOKAMAK_CORE_L){ // alpha = phi at inboard midplane
    if(Z<actx->zmaxis)
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
    else
      ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
  }

  else if (actx->ftype==GKYL_GEOMETRY_TOKAMAK_CORE_R){ // alpha = phi at outboard midplane
    if(Z<actx->zmaxis)
      ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
    else
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
  }

  else if( (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT) || (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID)){ // alpha = phi at outboard midplane
    if (Z<actx->zmaxis)
      ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
    else
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
  }
  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO){ // alpha = phi at lower plate and increases towards xpt
    ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo);
  }
  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP){ //alpha = phi at upper plate and decreases towards xpt
    ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo);
  }
  if( (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN) || (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID) ){ // alpha = phi at inboard midplane
    if (Z<actx->zmaxis)
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
    else
      ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
  }
  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO){ // alpha = phi at lower plate and decreases towards xpt
    ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo);
  }
  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP){// alpha = phi at upper plate and increases towards xpt
    ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo);
  }
  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_LSN_SOL  || actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID){
    // alpha = phi at outboard midplane
    if (actx->right==true){
      if (Z<actx->zmaxis)
        ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
      else
        ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
    }
    else{
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmax, rclose, false, false, arc_memo);
      phi_ref = actx->phi_right;
    }
  }
  else if(actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO){ //alpha = phi at outer plate and increases towards xpt
    ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo);
  }
  else if(actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP){ //alpha = phi at inner plate and decreases towards xpt
    ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo);
  }
  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_PF_LO_R){ // alpha = phi at outer plate and increases towards xpt
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo);
  }
  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_PF_LO_L){ //alpha = phi at inner plate and decreases towards xpt
      ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo);
  }
  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_PF_UP_R){ // alpha = phi at outer plate and decreases towards Xpt
      ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo);
  }
  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_PF_UP_L){ // alpha = phi at inner plate and increases towards xpt
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo);
  }
  else if (actx->ftype==GKYL_GEOMETRY_TOKAMAK_IWL) {
    // phi = alpha at outboard midplane
    if(actx->right==true){
      if(Z<actx->zmaxis)
        ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
      else
        ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
    }
    else{
      if (Z<actx->zmaxis) {
        ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmin, Z, rclose, false, false, arc_memo) ;
        phi_ref  = -actx->phi_right;
      }
      else {
        ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmax, rclose, false, false, arc_memo) ;
        phi_ref  = actx->phi_right;
      }
    }
  }
  // Now multiply by fpol
  double R[4] = {0};
  double dRdZ[4] = {0};
  double dR[4] = {0};
  double dZ[4] = {0};
  int nr = gkyl_tok_geo_R_psiZ(actx->geo, psi, Z, 4, R, dRdZ, dR, dZ);
  double r_curr = nr == 1 ? R[0] : choose_closest(rclose, R, R, nr);
  double psi_fpol = psi;
  if ( (psi_fpol < actx->geo->fgrid.lower[0]) || (psi_fpol > actx->geo->fgrid.upper[0]) ) // F = F(psi_sep) in the SOL.
    psi_fpol = actx->geo->sibry;
  int idx = fmin(actx->geo->frange.lower[0] + (int) floor((psi_fpol - actx->geo->fgrid.lower[0])/actx->geo->fgrid.dx[0]), actx->geo->frange.upper[0]);
  long loc = gkyl_range_idx(&actx->geo->frange, &idx);
  const double *coeffs = gkyl_array_cfetch(actx->geo->fpoldg,loc);
  double fxc;
  gkyl_rect_grid_cell_center(&actx->geo->fgrid, &idx, &fxc);
  double fx = (psi_fpol-fxc)/(actx->geo->fgrid.dx[0]*0.5);
  double fpol = actx->geo->fbasis.eval_expand(&fx, coeffs);
  ival = ival*fpol;

  //while(ival < -M_PI){
  //  ival +=2*M_PI;
  //}
  //while(ival > M_PI){
  //  ival -=2*M_PI;
  //}
  return alpha_curr + ival + phi_ref;
}

double
qprofile_func(void *ctx)
{
  // Function to calculate the flux surface averaged q profile.

  struct arc_length_ctx *actx = ctx;
  double *arc_memo = actx->arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL, zmax = actx->zmax;
  double rleft =  actx->rleft, rright = actx->rright;

  // Calculate q(psi) = -F(psi)/2pi * integral_zmin^zmax 1/Rgrad(psi).

  double ival = 0;
  double phi_ref = 0.0;
  if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_CORE) {
    double ival1 = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmin, actx->zmax, rright, false, false, arc_memo);
    double ival2 = -integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmin, actx->zmax, rleft, false, false, arc_memo);
    ival = ival1 + ival2;
  }

  if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_IWL) {
    double ival1 = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmin, actx->zmax, rright, false, false, arc_memo);
    double ival2 = -integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmin, actx->zmax, rleft, false, false, arc_memo);
    ival = ival1 + ival2;
  }

  if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL) {
    double ival1 = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmin_right, actx->zmax, rright, false, false, arc_memo);
    double ival2 = -integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmin_left, actx->zmax, rleft, false, false, arc_memo);
    ival = ival1 + ival2;
  }

  if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT || actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID || actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO || actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP){
    ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmin, actx->zmax, rclose, false, false, arc_memo);
  }

  if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN || actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID || actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO || actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP) {
    ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmin, actx->zmax, rclose, false, false, arc_memo);
  }

  // Now multiply by fpol/2pi.
  double R[4] = {0};
  double dR[4] = {0};
  double psi_fpol = psi;
  if ( (psi_fpol < actx->geo->fgrid.lower[0]) || (psi_fpol > actx->geo->fgrid.upper[0]) ) // F = F(psi_sep) in the SOL.
    psi_fpol = actx->geo->sibry;
  int idx = fmin(actx->geo->frange.lower[0] + (int) floor((psi_fpol - actx->geo->fgrid.lower[0])/actx->geo->fgrid.dx[0]), actx->geo->frange.upper[0]);
  long loc = gkyl_range_idx(&actx->geo->frange, &idx);
  const double *coeffs = gkyl_array_cfetch(actx->geo->fpoldg,loc);
  double fxc;
  gkyl_rect_grid_cell_center(&actx->geo->fgrid, &idx, &fxc);
  double fx = (psi_fpol-fxc)/(actx->geo->fgrid.dx[0]*0.5);
  double fpol = actx->geo->fbasis.eval_expand(&fx, coeffs);
  double qout = -ival*fpol/M_PI;

  // AS 1/15/25: The 3 lines below are a useful check to compare against q from efit.
  //coeffs = gkyl_array_cfetch(actx->geo->qdg,loc);
  //double q_efit = actx->geo->fbasis.eval_expand(&fx, coeffs);
  // printf("psi_curr = %g, my q = %g, efit q = %g\n", psi_fpol, qout, q_efit);

  return qout;
}

static double
dphidtheta_func(double Z, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo = actx->arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL, zmax = actx->zmax;

  // Get the integrand
  double integrand = 0.0;
  struct contour_ctx cctx = {
    .geo = actx->geo,
    .psi = psi,
    .ncall = 0,
    .last_R = rclose
  };
  integrand = dphidtheta_integrand(Z, &cctx);
  // Now multiply by fpol
  double R[4] = {0};
  double dRdZ[4] = {0};
  double dR[4] = {0};
  double dZ[4] = {0};
  int nr = gkyl_tok_geo_R_psiZ(actx->geo, psi, Z, 4, R, dRdZ, dR, dZ);
  double r_curr = nr == 1 ? R[0] : choose_closest(rclose, R, R, nr);
  double psi_fpol = psi;
  if ( (psi_fpol < actx->geo->fgrid.lower[0]) || (psi_fpol > actx->geo->fgrid.upper[0]) ) // F = F(psi_sep) in the SOL.
    psi_fpol = actx->geo->sibry;
  int idx = fmin(actx->geo->frange.lower[0] + (int) floor((psi_fpol - actx->geo->fgrid.lower[0])/actx->geo->fgrid.dx[0]), actx->geo->frange.upper[0]);
  long loc = gkyl_range_idx(&actx->geo->frange, &idx);
  const double *coeffs = gkyl_array_cfetch(actx->geo->fpoldg,loc);
  double fxc;
  gkyl_rect_grid_cell_center(&actx->geo->fgrid, &idx, &fxc);
  double fx = (psi_fpol-fxc)/(actx->geo->fgrid.dx[0]*0.5);
  double fpol = actx->geo->fbasis.eval_expand(&fx, coeffs);
  integrand = integrand*fpol;
  double darc_dtheta = actx->xpt_map_valid
    ? actx->xpt_map_darc_dtheta : actx->arcL_tot/(2.0*M_PI);
  integrand = integrand*darc_dtheta;
  return integrand;
}

static double
bmag_func(double r_curr, double Z, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo = actx->arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL, zmax = actx->zmax;
  // Calculate fpol
  double psi_fpol = psi;
  if ( (psi_fpol < actx->geo->fgrid.lower[0]) || (psi_fpol > actx->geo->fgrid.upper[0]) ) // F = F(psi_sep) in the SOL.
    psi_fpol = actx->geo->sibry;
  int idx = fmin(actx->geo->frange.lower[0] + (int) floor((psi_fpol - actx->geo->fgrid.lower[0])/actx->geo->fgrid.dx[0]), actx->geo->frange.upper[0]);
  long loc = gkyl_range_idx(&actx->geo->frange, &idx);
  const double *coeffs = gkyl_array_cfetch(actx->geo->fpoldg,loc);
  double fxc;
  gkyl_rect_grid_cell_center(&actx->geo->fgrid, &idx, &fxc);
  double fx = (psi_fpol-fxc)/(actx->geo->fgrid.dx[0]*0.5);
  double fpol = actx->geo->fbasis.eval_expand(&fx, coeffs);
  double Bphi = fpol/r_curr;
  double Br = 0.0, Bz = 0.0, bmag = 0.0;

  if (actx->geo->use_cubics) {
  double xn[2] = {r_curr, Z};
  double fout[3];
  actx->geo->efit->evf->eval_cubic_wgrad(0.0, xn, fout, actx->geo->efit->evf->ctx);
  double dpsidR = fout[1];
  double dpsidZ = fout[2];

  Br = 1.0/r_curr*dpsidZ;
  Bz = -1.0/r_curr*dpsidR;
  }
  else {
    int rzidx[2];
    int idxtemp = actx->geo->rzlocal.lower[0] + (int) floor((r_curr - actx->geo->rzgrid.lower[0])/actx->geo->rzgrid.dx[0]);
    idxtemp = GKYL_MIN2(idxtemp, actx->geo->rzlocal.upper[0]);
    idxtemp = GKYL_MAX2(idxtemp, actx->geo->rzlocal.lower[0]);
    rzidx[0] = idxtemp;
    idxtemp = actx->geo->rzlocal.lower[1] + (int) floor((Z - actx->geo->rzgrid.lower[1])/actx->geo->rzgrid.dx[1]);
    idxtemp = GKYL_MIN2(idxtemp, actx->geo->rzlocal.upper[1]);
    idxtemp = GKYL_MAX2(idxtemp, actx->geo->rzlocal.lower[1]);
    rzidx[1] = idxtemp;

    long loc = gkyl_range_idx((&actx->geo->rzlocal), rzidx);
    const double *psih = gkyl_array_cfetch(actx->geo->psiRZ, loc);

    double xc[2];
    gkyl_rect_grid_cell_center((&actx->geo->rzgrid), rzidx, xc);
    double x = (r_curr-xc[0])/(actx->geo->rzgrid.dx[0]*0.5);
    double y = (Z-xc[1])/(actx->geo->rzgrid.dx[1]*0.5);

    double dpsidx = 5.625*psih[8]*(2.0*x*SQ(y)-0.6666666666666666*x)+2.904737509655563*psih[7]*(SQ(y)-0.3333333333333333)+5.809475019311126*psih[6]*x*y+1.5*psih[3]*y+3.354101966249684*psih[4]*x+0.8660254037844386*psih[1];
    double dpsidy = 5.625*psih[8]*(2.0*SQ(x)*y-0.6666666666666666*y)+5.809475019311126*psih[7]*x*y+3.354101966249684*psih[5]*y+2.904737509655563*psih[6]*(SQ(x)-0.3333333333333333)+1.5*psih[3]*x+0.8660254037844386*psih[2];
    double dpsidR = dpsidx*2.0/actx->geo->rzgrid.dx[0];
    double dpsidZ = dpsidy*2.0/actx->geo->rzgrid.dx[1];
    Br = 1.0/r_curr*dpsidZ;
    Bz = -1.0/r_curr*dpsidR;
  }
  bmag = sqrt(Br*Br+Bz*Bz+Bphi*Bphi);
  return bmag;
}

static void
curlbhat_func(double psi, double r_curr, double Z, double phi, double *curlbhat, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo = actx->arc_memo;
  // Calculate fpol and fpolprime
  // First get the location on the flux grid
  double psi_fpol = psi;
  if ( (psi_fpol < actx->geo->fgrid.lower[0]) || (psi_fpol > actx->geo->fgrid.upper[0]) ) // F = F(psi_sep) in the SOL.
    psi_fpol = actx->geo->sibry;
  int idx = fmin(actx->geo->frange.lower[0] + (int) floor((psi_fpol - actx->geo->fgrid.lower[0])/actx->geo->fgrid.dx[0]), actx->geo->frange.upper[0]);
  long loc = gkyl_range_idx(&actx->geo->frange, &idx);
  const double *coeffs = gkyl_array_cfetch(actx->geo->fpoldg,loc);
  double fxc;
  gkyl_rect_grid_cell_center(&actx->geo->fgrid, &idx, &fxc);
  double fx = (psi_fpol-fxc)/(actx->geo->fgrid.dx[0]*0.5);
  // Second calculate fpol and bphi 
  double fpol = actx->geo->fbasis.eval_expand(&fx, coeffs);
  double Bphi = fpol/r_curr;
  // Third calculate fpolprime
  coeffs = gkyl_array_cfetch(actx->geo->fpolprimedg,loc);
  double fpolprime = actx->geo->fbasis.eval_expand(&fx, coeffs);


  // Now calculate psi's various derivatives
  double Br = 0.0, Bz = 0.0, bmag = 0.0;
  double dpsidR = 0.0, dpsidZ = 0.0;
  double d2psidR2 = 0.0, d2psidZ2 = 0.0, d2psidRdZ = 0.0;
  double dBdR = 0.0, dBdZ = 0.0;
  double dBrdR = 0.0, dBrdZ = 0.0;
  double dBzdR = 0.0, dBzdZ = 0.0;

  if (actx->geo->use_cubics) {
    double xn[2] = {r_curr, Z};
    double fout[4];
    actx->geo->efit->evf->eval_cubic_wgrad(0.0, xn, fout, actx->geo->efit->evf->ctx);
    dpsidR = fout[1];
    dpsidZ = fout[2];
    actx->geo->efit->evf->eval_cubic_wgrad2(0.0, xn, fout, actx->geo->efit->evf->ctx);
    d2psidR2 = fout[1];
    d2psidZ2 = fout[2];
    d2psidRdZ = fout[3];
  }
  else {
    int rzidx[2];
    int idxtemp = actx->geo->rzlocal.lower[0] + (int) floor((r_curr - actx->geo->rzgrid.lower[0])/actx->geo->rzgrid.dx[0]);
    idxtemp = GKYL_MIN2(idxtemp, actx->geo->rzlocal.upper[0]);
    idxtemp = GKYL_MAX2(idxtemp, actx->geo->rzlocal.lower[0]);
    rzidx[0] = idxtemp;
    idxtemp = actx->geo->rzlocal.lower[1] + (int) floor((Z - actx->geo->rzgrid.lower[1])/actx->geo->rzgrid.dx[1]);
    idxtemp = GKYL_MIN2(idxtemp, actx->geo->rzlocal.upper[1]);
    idxtemp = GKYL_MAX2(idxtemp, actx->geo->rzlocal.lower[1]);
    rzidx[1] = idxtemp;

    long loc = gkyl_range_idx((&actx->geo->rzlocal), rzidx);
    const double *psih = gkyl_array_cfetch(actx->geo->psiRZ, loc);

    double xc[2];
    gkyl_rect_grid_cell_center((&actx->geo->rzgrid), rzidx, xc);
    double x = (r_curr-xc[0])/(actx->geo->rzgrid.dx[0]*0.5);
    double y = (Z-xc[1])/(actx->geo->rzgrid.dx[1]*0.5);

    double dpsidx = 5.625*psih[8]*(2.0*x*SQ(y)-0.6666666666666666*x)+2.904737509655563*psih[7]*(SQ(y)-0.3333333333333333)+5.809475019311126*psih[6]*x*y+1.5*psih[3]*y+3.354101966249684*psih[4]*x+0.8660254037844386*psih[1];
    double dpsidy = 5.625*psih[8]*(2.0*SQ(x)*y-0.6666666666666666*y)+5.809475019311126*psih[7]*x*y+3.354101966249684*psih[5]*y+2.904737509655563*psih[6]*(SQ(x)-0.3333333333333333)+1.5*psih[3]*x+0.8660254037844386*psih[2];
    dpsidR = dpsidx*2.0/actx->geo->rzgrid.dx[0];
    dpsidZ = dpsidy*2.0/actx->geo->rzgrid.dx[1];
    double d2psidx2 = 11.25*psih[8]*y*y+5.809475019311125*psih[6]*y-3.75*psih[8]+3.354101966249685*psih[4];
    double d2psidy2 = 11.25*psih[8]*x*x+5.809475019311125*psih[7]*x-3.75*psih[8]+3.354101966249685*psih[5];
    double d2psidxdy = 22.5*psih[8]*x*y+5.809475019311125*psih[7]*y+5.809475019311125*psih[6]*x+1.5*psih[3]; 
    d2psidR2 = d2psidx2*(2.0/actx->geo->rzgrid.dx[0])*(2.0/actx->geo->rzgrid.dx[0]);
    d2psidZ2 = d2psidy2*(2.0/actx->geo->rzgrid.dx[1])*(2.0/actx->geo->rzgrid.dx[1]);
    d2psidRdZ = d2psidxdy*(2.0/actx->geo->rzgrid.dx[0])*(2.0/actx->geo->rzgrid.dx[1]);
  }

  Br = 1.0/r_curr*dpsidZ;
  Bz = -1.0/r_curr*dpsidR;
  bmag = sqrt(Br*Br+Bz*Bz+Bphi*Bphi);

  dBrdR  = 1.0/r_curr*d2psidRdZ;
  dBrdZ  = 1.0/r_curr*d2psidZ2;
  dBzdR  = -1.0/r_curr*d2psidR2;
  dBzdZ  = -1.0/r_curr*d2psidRdZ;


  double dFdR = fpolprime*dpsidR;
  double dFdZ = fpolprime*dpsidZ;
  dBdR = 1/bmag*(Br*dBrdR + Bz*dBzdR + fpol/r_curr*(dFdR/r_curr - fpol/r_curr/r_curr));
  dBdZ = 1/bmag*(Br*dBrdZ + Bz*dBzdZ + fpol/r_curr*dFdZ);

  // Get the polar components (contravariant, upperscript components on tangent basis)
  double polar_comp[3] = {0.0};
  polar_comp[0] = 1.0/bmag*-1.0/r_curr*dFdZ - 1.0/bmag*1.0/bmag*dBdZ*fpol/r_curr; // R component ^1
  polar_comp[1] = 1.0/bmag*1.0/r_curr*(dBrdZ - dBzdR) + (-dBdR*Bz/r_curr + dBdZ*Br/r_curr); // Phi component ^2
  polar_comp[2] = 1.0/bmag*1.0/r_curr*dFdR + 1.0/bmag*1.0/bmag*dBdR*fpol/r_curr;

  // Convert to cartesian
  curlbhat[0] = polar_comp[0] * cos(phi) - polar_comp[1] * sin(phi) * r_curr;
  curlbhat[1] = polar_comp[0] * sin(phi) + polar_comp[1] * cos(phi) * r_curr;
  curlbhat[2] = polar_comp[2];
}






struct gkyl_tok_geo*
gkyl_tok_geo_new(const struct gkyl_efit_inp *inp, const struct gkyl_tok_geo_grid_inp *ginp)
{
  struct gkyl_tok_geo *geo = gkyl_malloc(sizeof(*geo));
  *geo = (struct gkyl_tok_geo) {};

  geo->efit = gkyl_efit_new(inp);

  geo->plate_spec = ginp->plate_spec;
  geo->plate_func_lower = ginp->plate_func_lower;
  geo->plate_func_upper = ginp->plate_func_upper;

  geo->rzbasis = geo->efit->rzbasis;
  geo->rzbasis_cubic = geo->efit->rzbasis_cubic;
  geo->rzgrid = geo->efit->rzgrid;
  geo->rzgrid_cubic = geo->efit->rzgrid_cubic;
  geo->psiRZ = gkyl_array_acquire(geo->efit->psizr);
  geo->psiRZ_cubic = gkyl_array_acquire(geo->efit->psizr_cubic);

  geo->num_rzbasis = geo->rzbasis.num_basis;
  geo->rzlocal = geo->efit->rzlocal;
  geo->rzlocal_ext = geo->efit->rzlocal_ext;
  geo->rzlocal_cubic = geo->efit->rzlocal_cubic;
  geo->rzlocal_cubic_ext = geo->efit->rzlocal_cubic_ext;
  geo->fgrid = geo->efit->fluxgrid;
  geo->fbasis = geo->efit->fluxbasis;
  geo->frange = geo->efit->fluxlocal;
  geo->frange_ext = geo->efit->fluxlocal_ext;
  geo->fpoldg= gkyl_array_acquire(geo->efit->fpolflux);
  geo->fpolprimedg= gkyl_array_acquire(geo->efit->fpolprimeflux);
  geo->qdg= gkyl_array_acquire(geo->efit->qflux);
  geo->sibry = geo->efit->sibry;
  geo->psisep = geo->efit->psisep;
  geo->zmaxis = geo->efit->zmaxis;

  geo->use_cubics = ginp->use_cubics;
  geo->use_hyperbolic_numbers = ginp->use_hyperbolic_numbers;
  geo->root_param.eps =
    ginp->root_param.eps > 0 ? ginp->root_param.eps : 1e-10;
  geo->root_param.max_iter =
    ginp->root_param.max_iter > 0 ? ginp->root_param.max_iter : 100;

  geo->quad_param.max_level =
    ginp->quad_param.max_levels > 0 ? ginp->quad_param.max_levels : 10;
  geo->quad_param.eps =
    ginp->quad_param.eps > 0 ? ginp->quad_param.eps : 1e-10;

  if (geo->use_cubics) {
    if(geo->use_hyperbolic_numbers)
      geo->calc_roots = calc_RdR_p3_hyperbolic;
    else
      geo->calc_roots = calc_RdR_p3;
    geo->calc_grad_psi = calc_grad_psi_p3;
  }
  else if (geo->efit->rzbasis.poly_order == 1) {
    geo->calc_roots = calc_RdR_p1;
    geo->calc_grad_psi = calc_grad_psi_p1;
  }
  else if (geo->efit->rzbasis.poly_order == 2){
    geo->calc_roots = calc_RdR_p2_tensor_nrc;
    geo->calc_grad_psi = calc_grad_psi_p2_tensor;
  }

  geo->stat = (struct gkyl_tok_geo_stat) { };

  
  return geo;
}

double
gkyl_tok_geo_integrate_psi_contour(const struct gkyl_tok_geo *geo, double psi,
  double zmin, double zmax, double rclose)
{
  return integrate_psi_contour_memo(geo, psi, zmin, zmax, rclose,
    false, false, 0);
}

int
gkyl_tok_geo_R_psiZ(const struct gkyl_tok_geo *geo, double psi, double Z, int nmaxroots,
  double *R, double *dRdZ, double *dR, double *dZ)
{
  if(geo->use_cubics)
    return R_psiZ_cubic(geo, psi, Z, nmaxroots, R, dRdZ, dR, dZ);
  else
    return R_psiZ(geo, psi, Z, nmaxroots, R, dRdZ, dR, dZ);
}

void gkyl_tok_geo_calc(struct gk_geometry* up, struct gkyl_range *nrange, struct gkyl_tok_geo *geo, 
  struct gkyl_tok_geo_grid_inp *inp, struct gkyl_position_map *position_map)
{

  geo->rleft = inp->rleft;
  geo->rright = inp->rright;

  geo->inexact_roots = inp->inexact_roots;

  geo->rmax = inp->rmax;
  geo->rmin = inp->rmin;

  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  
  double dtheta = inp->cgrid.dx[TH_IDX],
    dpsi = inp->cgrid.dx[PSI_IDX],
    dalpha = inp->cgrid.dx[AL_IDX];
  
  double theta_lo = up->grid.lower[TH_IDX] + (up->local.lower[TH_IDX] - up->global.lower[TH_IDX])*up->grid.dx[TH_IDX],
    psi_lo = up->grid.lower[PSI_IDX] + (up->local.lower[PSI_IDX] - up->global.lower[PSI_IDX])*up->grid.dx[PSI_IDX],
    alpha_lo = up->grid.lower[AL_IDX] + (up->local.lower[AL_IDX] - up->global.lower[AL_IDX])*up->grid.dx[AL_IDX];
    
  double dx_fact = up->basis.poly_order == 1.0/up->basis.poly_order;
  dtheta *= dx_fact; dpsi *= dx_fact; dalpha *= dx_fact;

  double rclose = inp->rclose;
  double rright = inp->rright;
  double rleft = inp->rleft;

  int nzcells;
  if(geo->use_cubics)
    nzcells = geo->rzgrid_cubic.cells[1];
  else
    nzcells = geo->rzgrid.cells[1];
  double *arc_memo = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_left = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_right = gkyl_malloc(sizeof(double[nzcells]));
  int sep_trace_capacity = 16*nzcells+1;
  double *sep_trace_r = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *sep_trace_z = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *sep_trace_s = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *ordered_trace_storage =
    gkyl_malloc(sizeof(double[8*sep_trace_capacity]));

  struct arc_length_ctx arc_ctx = {
    .geo = geo,
    .arc_memo = arc_memo,
    .arc_memo_right = arc_memo_right,
    .arc_memo_left = arc_memo_left,
    .sep_trace_r = sep_trace_r,
    .sep_trace_z = sep_trace_z,
    .sep_trace_s = sep_trace_s,
    .sep_trace_capacity = sep_trace_capacity,
    .far_trace_r = ordered_trace_storage,
    .far_trace_z = ordered_trace_storage+sep_trace_capacity,
    .far_trace_s = ordered_trace_storage+2*sep_trace_capacity,
    .trace_corr_v = ordered_trace_storage+3*sep_trace_capacity,
    .map_trace_r = ordered_trace_storage+4*sep_trace_capacity,
    .map_trace_z = ordered_trace_storage+5*sep_trace_capacity,
    .map_trace_s = ordered_trace_storage+6*sep_trace_capacity,
    .map_trace_phi = ordered_trace_storage+7*sep_trace_capacity,
    .ftype = inp->ftype,
    .zmaxis = geo->zmaxis
  };
  struct plate_ctx pctx = {
    .geo = geo
  };
  tok_init_xpt_ray_target(inp, position_map, &arc_ctx);

  int cidx[3] = { 0 };
  for (int ia=nrange->lower[AL_IDX]; ia<nrange->lower[AL_IDX]+1; ++ia){
    cidx[AL_IDX] = ia;
    double alpha_curr = alpha_lo + ia*dalpha;

    for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
      double psi_curr = psi_lo + ip*dpsi;

      // Non-uniform psi. Finite differences are calculated in calc_metric.c
      position_map->maps[0](0.0, &psi_curr,  &psi_curr,  position_map->ctxs[0]);

      double darcL, arcL_curr, arcL_lo;

      // For double null blocks this should set arc_ctx :
      // zmin, zmax, rclose, arcL_tot for all blocks. No left and right
      // For a full core case:
      // also set phi_right and arcL_right
      // For a single null case:
      // also set zmin_left and zmin_right 
      if (tok_xpt_mapping_requested(inp))
        tok_prepare_ordered_map(inp, &arc_ctx, psi_curr);
      else
        tok_find_endpoints(inp, geo, &arc_ctx, &pctx, psi_curr, alpha_curr,
          arc_memo, arc_memo_left, arc_memo_right);

      darcL = arc_ctx.arcL_tot/(up->basis.poly_order*inp->cgrid.cells[TH_IDX])
        * (inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX])/2/M_PI;
      // At the beginning of each theta loop we need to reset things.
      cidx[PSI_IDX] = ip;
      arcL_curr = 0.0;
      arcL_lo = (theta_lo + M_PI)/2/M_PI*arc_ctx.arcL_tot;
      double ridders_min, ridders_max;
      // Set node coordinates.
      for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
        int it_delta = 0;
        arcL_curr = arcL_lo + it*darcL;
        double theta_curr = arcL_curr*(2*M_PI/arc_ctx.arcL_tot) - M_PI ; 

        // Calculate derivatives using finite difference for ddtheta,
        // as well as transform the computational coordiante to the non-uniform field-aligned value

        // We cannot do non-uniform alpha because we are modeling axisymmetric systems
        // Non-uniform theta
        double Theta_curr;
        position_map->maps[2](0.0, &theta_curr,  &Theta_curr,  position_map->ctxs[2]);
        theta_curr = Theta_curr;

        struct tok_ordered_point ordered = { 0.0 };
        bool ordered_mapping = tok_ordered_map_lookup(inp, &arc_ctx,
          theta_curr, alpha_curr, &ordered);
        if (tok_xpt_mapping_requested(inp) && !ordered_mapping) {
          fprintf(stderr,
            "TOK_ORDERED_MAP lookup failed ftype=%d psi=%.17g theta=%.17g\n",
            inp->ftype, psi_curr, theta_curr);
          abort();
        }
        double r_curr = 0.0, z_curr = 0.0, phi_curr = 0.0;
        if (ordered_mapping) {
          r_curr = ordered.r; z_curr = ordered.z; phi_curr = ordered.phi;
        }
        else {
          arcL_curr = tok_xpt_theta_to_arc(inp, &arc_ctx, theta_curr);

          tok_set_ridders(inp, &arc_ctx, psi_curr, arcL_curr, &rclose, &ridders_min, &ridders_max);

          struct gkyl_qr_res res = gkyl_ridders(arc_length_func, &arc_ctx,
            arc_ctx.zmin, arc_ctx.zmax, ridders_min, ridders_max,
            geo->root_param.max_iter, 1e-10);
          z_curr = res.res;
          ((struct gkyl_tok_geo *)geo)->stat.nroot_cont_calls += res.nevals;

        if (inp->half_domain) { //Alternative for half domain
          if (tok_geo_same_flux(psi_curr, geo->psisep)) {
            if (it == nrange->upper[TH_IDX] && (up->local.upper[TH_IDX]== up->global.upper[TH_IDX])) {
              if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L || inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO|| inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID)
                z_curr = geo->efit->Zxpt[0];
              else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO || inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE)
                z_curr = geo->efit->Zxpt[0];
            }
            if (it == nrange->lower[TH_IDX] && (up->local.lower[TH_IDX]== up->global.lower[TH_IDX])) {
              if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_L || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R || inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID|| inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO)
                z_curr = geo->efit->Zxpt[0];
              else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP || inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE)
                z_curr = geo->efit->Zxpt[0];
            }
          }
        }
        else { // For full domain
          if (tok_geo_same_flux(psi_curr, geo->psisep)) {
            if (it == nrange->upper[TH_IDX] && (up->local.upper[TH_IDX]== up->global.upper[TH_IDX])) {
              if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_L || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R || inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID || inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP)
                z_curr = geo->efit->Zxpt[1];
              else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L || inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO|| inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID)
                z_curr = geo->efit->Zxpt[0];
              else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO || inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE)
                z_curr = geo->efit->Zxpt[0];
            }
            if (it == nrange->lower[TH_IDX] && (up->local.lower[TH_IDX]== up->global.lower[TH_IDX])) {
              if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_R || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L || inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP|| inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID)
                z_curr = geo->efit->Zxpt[1];
              else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_L || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R || inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID|| inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO)
                z_curr = geo->efit->Zxpt[0];
              else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP || inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE)
                z_curr = geo->efit->Zxpt[0];
            }
          }
        }

        double sep_r_curr = 0.0;
        bool at_sep_trace = tok_half_domain_sep_rz(inp, &arc_ctx,
          theta_curr, &sep_r_curr, &z_curr);

        double fixed_root_ref = rclose;
        bool at_fixed_edge = tok_xpt_at_fixed_edge(inp, &arc_ctx, it,
          nrange,
          up->local.lower[TH_IDX] == up->global.lower[TH_IDX],
          up->local.upper[TH_IDX] == up->global.upper[TH_IDX],
          &z_curr, &fixed_root_ref);
        bool at_xpt_anchor = arc_ctx.xpt_anchor_valid &&
          tok_xpt_at_seam(inp, it, nrange,
            up->local.lower[TH_IDX] == up->global.lower[TH_IDX],
            up->local.upper[TH_IDX] == up->global.upper[TH_IDX]);
        if (at_xpt_anchor)
          z_curr = arc_ctx.xpt_anchor_z;

        double R[4] = { 0 }, dRdZ[4] = { 0 };
        double dR[4] = { 0 }, dZ[4] = { 0 };
        int nr = gkyl_tok_geo_R_psiZ(geo, psi_curr, z_curr, 4, R, dRdZ, dR, dZ);
        double root_ref = at_sep_trace ? sep_r_curr
          : (at_xpt_anchor ? arc_ctx.xpt_anchor_r
            : (at_fixed_edge ? fixed_root_ref : rclose));
        r_curr = choose_closest(root_ref, R, R, nr);
        double drdz_curr = choose_closest(root_ref, R, dRdZ, nr);
        double dr_curr = choose_closest(root_ref, R, dR, nr);
        double dz_curr = choose_closest(root_ref, R, dZ, nr);
        if (at_xpt_anchor)
          r_curr = arc_ctx.xpt_anchor_r;
        if (at_sep_trace)
          r_curr = sep_r_curr;

        if (tok_geo_same_flux(psi_curr, geo->psisep)) {
          if (z_curr == geo->efit->Zxpt[0]) {
            nr = 1;
            r_curr = geo->efit->Rxpt[0];
          }
          if (z_curr == geo->efit->Zxpt[1]) {
            nr = 1;
            r_curr = geo->efit->Rxpt[1];
          }
        }

        if (nr==0) {
          printf(" ip = %d, it = %d, ia = %d\n", ip, it, ia);
          printf("Block Type = %d | Failed to find a root at psi = %g, Z = %1.16f\n", inp->ftype, psi_curr, z_curr);
          assert(false);
        }

          phi_curr = phi_func(alpha_curr, z_curr, &arc_ctx);
        }
        cidx[TH_IDX] = it;
        double *mc2p_n = gkyl_array_fetch(up->geo_corn.mc2p_nodal, gkyl_range_idx(nrange, cidx));
        double *mc2nu_n = gkyl_array_fetch(up->geo_corn.mc2nu_pos_nodal, gkyl_range_idx(nrange, cidx));
        double *bmag_n = gkyl_array_fetch(up->geo_corn.bmag_nodal, gkyl_range_idx(nrange, cidx));

        mc2p_n[X_IDX] = r_curr;
        mc2p_n[Y_IDX] = z_curr;
        mc2p_n[Z_IDX] = phi_curr;
        mc2nu_n[X_IDX] = psi_curr;
        mc2nu_n[Y_IDX] = alpha_curr;
        mc2nu_n[Z_IDX] = theta_curr;
        bmag_n[0] = bmag_func(r_curr, z_curr, &arc_ctx);
      }
    }
  }

  // Populate other alpha indices by using axisymmetry
  for (int ia=nrange->lower[AL_IDX]+1; ia<=nrange->upper[AL_IDX]; ++ia){
    cidx[AL_IDX] = ia;
    double alpha_curr = alpha_lo + ia*dalpha;
    double alpha_donor = alpha_lo + nrange->lower[AL_IDX]*dalpha;
    double alpha_diff = alpha_curr -  alpha_donor;
    for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
      cidx[PSI_IDX] = ip;
      for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
        cidx[TH_IDX] = it;

        double *mc2p_n = gkyl_array_fetch(up->geo_corn.mc2p_nodal, gkyl_range_idx(nrange, cidx));
        double *mc2nu_n = gkyl_array_fetch(up->geo_corn.mc2nu_pos_nodal, gkyl_range_idx(nrange, cidx));
        double *bmag_n = gkyl_array_fetch(up->geo_corn.bmag_nodal, gkyl_range_idx(nrange, cidx));

        int donor_cidx[3] ;
        donor_cidx[AL_IDX] = nrange->lower[AL_IDX];
        donor_cidx[PSI_IDX] = ip;
        donor_cidx[TH_IDX] = it;

        double *donor_mc2p_n = gkyl_array_fetch(up->geo_corn.mc2p_nodal, gkyl_range_idx(nrange, donor_cidx));
        double *donor_mc2nu_n = gkyl_array_fetch(up->geo_corn.mc2nu_pos_nodal, gkyl_range_idx(nrange, donor_cidx));
        double *donor_bmag_n = gkyl_array_fetch(up->geo_corn.bmag_nodal, gkyl_range_idx(nrange, donor_cidx));

        mc2p_n[X_IDX] = donor_mc2p_n[X_IDX];
        mc2p_n[Y_IDX] = donor_mc2p_n[Y_IDX];
        mc2p_n[Z_IDX] = donor_mc2p_n[Z_IDX] + alpha_diff;
        mc2nu_n[X_IDX] = donor_mc2nu_n[X_IDX];
        mc2nu_n[Y_IDX] = donor_mc2nu_n[AL_IDX] + alpha_diff;
        mc2nu_n[Z_IDX] = donor_mc2nu_n[Z_IDX];
        bmag_n[0] = donor_bmag_n[0];
      }
    }
  }

  struct gkyl_nodal_ops *n2m =  gkyl_nodal_ops_new(&inp->cbasis, &inp->cgrid, false);
  gkyl_nodal_ops_n2m(n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 3, up->geo_corn.mc2p_nodal, up->geo_corn.mc2p, false);
  gkyl_nodal_ops_n2m(n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 3, up->geo_corn.mc2nu_pos_nodal, up->geo_corn.mc2nu_pos, false);
  gkyl_nodal_ops_n2m(n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 1, up->geo_corn.bmag_nodal, up->geo_corn.bmag, false);
  gkyl_nodal_ops_release(n2m);

  // Need 1/B for LBO collisions, computed weakly.
  gkyl_dg_inv_op_range(&inp->cbasis, 0, up->geo_corn.bmag_inv, 0, up->geo_corn.bmag, &up->local);

  gkyl_free(arc_memo);
  gkyl_free(arc_memo_left);
  gkyl_free(arc_memo_right);
  gkyl_free(sep_trace_r);
  gkyl_free(sep_trace_z);
  gkyl_free(sep_trace_s);
  gkyl_free(ordered_trace_storage);
}

void gkyl_tok_geo_calc_interior(struct gk_geometry* up, struct gkyl_range *nrange, double dzc[3], 
    struct gkyl_tok_geo *geo, struct gkyl_tok_geo_grid_inp *inp, struct gkyl_position_map *position_map)
{

  geo->rleft = inp->rleft;
  geo->rright = inp->rright;

  geo->inexact_roots = inp->inexact_roots;

  geo->rmax = inp->rmax;
  geo->rmin = inp->rmin;

  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  
  double dtheta = inp->cgrid.dx[TH_IDX],
    dpsi = inp->cgrid.dx[PSI_IDX],
    dalpha = inp->cgrid.dx[AL_IDX];
  
  double theta_lo = up->grid.lower[TH_IDX] + (up->local.lower[TH_IDX] - up->global.lower[TH_IDX])*up->grid.dx[TH_IDX],
    psi_lo = up->grid.lower[PSI_IDX] + (up->local.lower[PSI_IDX] - up->global.lower[PSI_IDX])*up->grid.dx[PSI_IDX],
    alpha_lo = up->grid.lower[AL_IDX] + (up->local.lower[AL_IDX] - up->global.lower[AL_IDX])*up->grid.dx[AL_IDX];

  double dels[2] = {1.0/sqrt(3), 1.0-1.0/sqrt(3) };
  theta_lo = theta_lo + dels[1]*dtheta/2.0;
  psi_lo = psi_lo + dels[1]*dpsi/2.0;
  alpha_lo = alpha_lo + dels[1]*dalpha/2.0;
    
  double dx_fact = up->basis.poly_order == 1.0/up->basis.poly_order;
  dtheta *= dx_fact; dpsi *= dx_fact; dalpha *= dx_fact;

  // used for finite differences 
  double delta_alpha = dalpha*1e-2;
  double delta_psi = dpsi*1e-2;
  double delta_theta = dtheta*1e-2;
  dzc[0] = delta_psi;
  dzc[1] = delta_alpha;
  dzc[2] = delta_theta;
  int modifiers[5] = {0, -1, 1, -2, 2};

  double rclose = inp->rclose;
  double rright = inp->rright;
  double rleft = inp->rleft;

  int nzcells;
  if(geo->use_cubics)
    nzcells = geo->rzgrid_cubic.cells[1];
  else
    nzcells = geo->rzgrid.cells[1];
  double *arc_memo = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_left = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_right = gkyl_malloc(sizeof(double[nzcells]));
  int sep_trace_capacity = 16*nzcells+1;
  double *sep_trace_r = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *sep_trace_z = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *sep_trace_s = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *ordered_trace_storage =
    gkyl_malloc(sizeof(double[8*sep_trace_capacity]));

  struct arc_length_ctx arc_ctx = {
    .geo = geo,
    .arc_memo = arc_memo,
    .arc_memo_right = arc_memo_right,
    .arc_memo_left = arc_memo_left,
    .sep_trace_r = sep_trace_r,
    .sep_trace_z = sep_trace_z,
    .sep_trace_s = sep_trace_s,
    .sep_trace_capacity = sep_trace_capacity,
    .far_trace_r = ordered_trace_storage,
    .far_trace_z = ordered_trace_storage+sep_trace_capacity,
    .far_trace_s = ordered_trace_storage+2*sep_trace_capacity,
    .trace_corr_v = ordered_trace_storage+3*sep_trace_capacity,
    .map_trace_r = ordered_trace_storage+4*sep_trace_capacity,
    .map_trace_z = ordered_trace_storage+5*sep_trace_capacity,
    .map_trace_s = ordered_trace_storage+6*sep_trace_capacity,
    .map_trace_phi = ordered_trace_storage+7*sep_trace_capacity,
    .ftype = inp->ftype,
    .zmaxis = geo->zmaxis
  };
  struct plate_ctx pctx = {
    .geo = geo
  };
  tok_init_xpt_ray_target(inp, position_map, &arc_ctx);

  // Temporary array to store nodal q profile.
  struct gkyl_array *qprofile_nodal = gkyl_array_new(GKYL_DOUBLE, up->geo_int.bmag_nodal->ncomp, up->geo_int.bmag_nodal->size);

  int cidx[3] = { 0 };
  for(int ia=nrange->lower[AL_IDX]; ia<nrange->lower[AL_IDX]+1; ++ia){
    cidx[AL_IDX] = ia;
    double alpha_curr = calc_running_coord(alpha_lo, ia-nrange->lower[AL_IDX], dalpha);

    for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
      int ip_delta_max = 3;
      for (int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++) {

        double psi_curr = calc_running_coord(psi_lo, ip-nrange->lower[PSI_IDX], dpsi) + modifiers[ip_delta]*delta_psi;

        // Non-uniform psi. Finite differences are calculated in calc_metric.c
        double Psi_curr;
        position_map->maps[0](0.0, &psi_curr,  &Psi_curr,  position_map->ctxs[0]);
        double dPsi_dpsi = gkyl_position_map_slope(position_map, 0, psi_curr, delta_psi, ip, nrange);
        psi_curr = Psi_curr;
        
        double darcL, arcL_curr, arcL_lo;

        // For double null blocks this should set arc_ctx :
        // zmin, zmax, rclose, arcL_tot for all blocks. No left and right
        // For a full core case:
        // also set phi_right and arcL_right
        // For a single null case:
        // also set zmin_left and zmin_right 
        if (tok_xpt_mapping_requested(inp))
          tok_prepare_ordered_map(inp, &arc_ctx, psi_curr);
        else
          tok_find_endpoints(inp, geo, &arc_ctx, &pctx, psi_curr, alpha_curr,
            arc_memo, arc_memo_left, arc_memo_right);

        // Calculate the q profile
        // qhat = - F(psi) * s(psi) / (R * grad(psi))
        // q = integral_0^2pi qhat dtheta
        //   = F(psi)*s(psi) * integral 1/(R*grad(psi)) dl
        //   = 1/s(psi) * integral (dphidtheta) ; dphidtheta = F(psi)/(R*grad(psi))
        // The half-domain block types handled by the ordered map have always
        // returned zero from qprofile_func; avoid invoking its legacy contour
        // integration path.  Full-domain geometries retain the original path.
        double qprofile = tok_xpt_mapping_requested(inp)
          ? 0.0 : qprofile_func(&arc_ctx);

        darcL = arc_ctx.arcL_tot/(up->basis.poly_order*inp->cgrid.cells[TH_IDX])
          * (inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX])/2/M_PI;
        // at the beginning of each theta loop we need to reset things
        cidx[PSI_IDX] = ip;
        arcL_curr = 0.0;
        arcL_lo = (theta_lo + M_PI)/2/M_PI*arc_ctx.arcL_tot;
        double ridders_min, ridders_max;

        for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
          arcL_curr = calc_running_coord(arcL_lo, it-nrange->lower[TH_IDX], darcL);
          double theta_curr = arcL_curr*(2*M_PI/arc_ctx.arcL_tot) - M_PI ; 

          // Calculate derivatives using finite difference for ddtheta,
          // as well as transform the computational coordiante to the non-uniform field-aligned value

          // We cannot do non-uniform alpha because we are modeling axisymmetric systems
          // Non-uniform theta
          double Theta_curr;
          position_map->maps[2](0.0, &theta_curr,  &Theta_curr,  position_map->ctxs[2]);
          double dTheta_dtheta = gkyl_position_map_slope(position_map, 2, theta_curr,\
            delta_theta, it, nrange);
          theta_curr = Theta_curr;

          struct tok_ordered_point ordered = { 0.0 };
          bool ordered_mapping = tok_ordered_map_lookup(inp, &arc_ctx,
            theta_curr, alpha_curr, &ordered);
          if (tok_xpt_mapping_requested(inp) && !ordered_mapping) {
            fprintf(stderr,
              "TOK_ORDERED_MAP lookup failed ftype=%d psi=%.17g theta=%.17g\n",
              inp->ftype, psi_curr, theta_curr);
            abort();
          }
          double r_curr = 0.0, z_curr = 0.0, phi_curr = 0.0;
          double drdz_curr = 0.0, dr_curr = 0.0, dz_curr = 0.0;
          if (ordered_mapping) {
            r_curr = ordered.r; z_curr = ordered.z; phi_curr = ordered.phi;
          }
          else {
            arcL_curr = tok_xpt_theta_to_arc(inp, &arc_ctx, theta_curr);

            tok_set_ridders(inp, &arc_ctx, psi_curr, arcL_curr, &rclose, &ridders_min, &ridders_max);

          bool trace_this_block = tok_geo_trace_enabled() && tok_geo_trace_block_selected();
          if (trace_this_block && tok_geo_trace_int_env("GKYL_TOK_GEO_TRACE_BEFORE", 1))
            tok_geo_trace_surface_row("before_ridders", inp->ftype, -1, ip,
              it, ia, ip_delta, psi_curr, alpha_curr, theta_curr, arcL_curr,
              arc_ctx.zmin, arc_ctx.zmax, rclose, ridders_min, ridders_max,
              0.0, -1, 0.0);

          clock_t trace_t0 = trace_this_block ? clock() : 0;
          struct gkyl_qr_res res = gkyl_ridders(arc_length_func, &arc_ctx,
            arc_ctx.zmin, arc_ctx.zmax, ridders_min, ridders_max,
            geo->root_param.max_iter, 1e-10);
          double trace_elapsed = trace_this_block ?
            ((double) (clock() - trace_t0))/CLOCKS_PER_SEC : 0.0;
          z_curr = res.res;
          ((struct gkyl_tok_geo *)geo)->stat.nroot_cont_calls += res.nevals;

          if (trace_this_block) {
            int neval_threshold = tok_geo_trace_int_env("GKYL_TOK_GEO_TRACE_NEVAL_THRESHOLD", 25);
            double time_threshold = tok_geo_trace_double_env("GKYL_TOK_GEO_TRACE_TIME_THRESHOLD", 0.01);
            if (res.nevals >= neval_threshold || trace_elapsed >= time_threshold ||
              tok_geo_trace_int_env("GKYL_TOK_GEO_TRACE_ALL_AFTER", 0))
              tok_geo_trace_surface_row("after_ridders", inp->ftype, -1, ip,
                it, ia, ip_delta, psi_curr, alpha_curr, theta_curr, arcL_curr,
                arc_ctx.zmin, arc_ctx.zmax, rclose, ridders_min, ridders_max,
                z_curr, res.nevals, trace_elapsed);
          }

          double sep_r_curr = 0.0;
          bool at_sep_trace = tok_half_domain_sep_rz(inp, &arc_ctx,
            theta_curr, &sep_r_curr, &z_curr);

          // These are volume-interior Gauss nodes, even when their nodal
          // indices are the first or last in the local range.  The affine
          // arc map already places them correctly; only a true block-edge
          // node may be snapped exactly to the shared ray.
          bool at_xpt_anchor = false;

          double R[4] = { 0 }, dRdZ[4] = { 0 };
          double dR[4] = { 0 }, dZ[4] = { 0 };
          int nr = gkyl_tok_geo_R_psiZ(geo, psi_curr, z_curr, 4, R, dRdZ, dR, dZ);
          double root_ref = at_sep_trace ? sep_r_curr
            : (at_xpt_anchor ? arc_ctx.xpt_anchor_r : rclose);
          r_curr = choose_closest(root_ref, R, R, nr);
          drdz_curr = choose_closest(root_ref, R, dRdZ, nr);
          dr_curr = choose_closest(root_ref, R, dR, nr);
          dz_curr = choose_closest(root_ref, R, dZ, nr);
          if (at_xpt_anchor)
            r_curr = arc_ctx.xpt_anchor_r;
          if (at_sep_trace)
            r_curr = sep_r_curr;

          if (tok_geo_same_flux(psi_curr, geo->psisep) && ip_delta==0) {
            if (z_curr == geo->efit->Zxpt[0]) {
              nr = 1;
              r_curr = geo->efit->Rxpt[0];
            }
            if (z_curr == geo->efit->Zxpt[1]) {
              nr = 1;
              r_curr = geo->efit->Rxpt[1];
            }
          }

          if (nr==0) {
            printf("ip = %d, it = %d, ia = %d, ip_delta = %d\n", ip, it, ia, ip_delta);
            printf("Block Type = %d | Failed to find a root at psi = %g, Z = %1.16f\n", inp->ftype, psi_curr, z_curr);
            assert(false);
          }

            phi_curr = phi_func(alpha_curr, z_curr, &arc_ctx);
          }
          cidx[TH_IDX] = it;
          int lidx = 0;
          if (ip_delta != 0)
            lidx = 3 + 3*(ip_delta-1);
          double *mc2p_fd_n = gkyl_array_fetch(up->geo_int.mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
          double *ddtheta_n = gkyl_array_fetch(up->geo_int.ddtheta_nodal, gkyl_range_idx(nrange, cidx));
          double *ddpsi_n = gkyl_array_fetch(up->geo_int.ddpsi_nodal, gkyl_range_idx(nrange, cidx));
          double *mc2p_n = gkyl_array_fetch(up->geo_int.mc2p_nodal, gkyl_range_idx(nrange, cidx));
          double *bmag_n = gkyl_array_fetch(up->geo_int.bmag_nodal, gkyl_range_idx(nrange, cidx));
          double *curlbhat_n = gkyl_array_fetch(up->geo_int.curlbhat_nodal, gkyl_range_idx(nrange, cidx));
          double *qprofile_n = gkyl_array_fetch(qprofile_nodal, gkyl_range_idx(nrange, cidx));

          mc2p_fd_n[lidx+X_IDX] = r_curr;
          mc2p_fd_n[lidx+Y_IDX] = z_curr;
          mc2p_fd_n[lidx+Z_IDX] = phi_curr;

          if(ip_delta==0){
            if (ordered_mapping) {
              ddtheta_n[0] = ordered.dr_dtheta*dTheta_dtheta;
              ddtheta_n[1] = ordered.dz_dtheta*dTheta_dtheta;
              ddtheta_n[2] = ordered.dphi_dtheta*dTheta_dtheta;
            }
            else {
              double darc_dtheta = arc_ctx.xpt_map_valid
                ? arc_ctx.xpt_map_darc_dtheta : arc_ctx.arcL_tot/(2.0*M_PI);
              ddtheta_n[0] = sin(atan2(dr_curr, dz_curr))*darc_dtheta*dTheta_dtheta;
              ddtheta_n[1] = cos(atan2(dr_curr, dz_curr))*darc_dtheta*dTheta_dtheta;
              ddtheta_n[2] = dphidtheta_func(z_curr, &arc_ctx)*dTheta_dtheta;
            }
            ddpsi_n[0] = dPsi_dpsi;
            mc2p_n[lidx+X_IDX] = r_curr;
            mc2p_n[lidx+Y_IDX] = z_curr;
            mc2p_n[lidx+Z_IDX] = phi_curr;
            bmag_n[0] = bmag_func(r_curr, z_curr, &arc_ctx);
            qprofile_n[0] = qprofile;
            curlbhat_func(psi_curr, r_curr, z_curr, phi_curr, curlbhat_n, &arc_ctx);
          }
        }
      }
    }
  }

  // Populate other alpha indices by using axisymmetry
  for (int ia=nrange->lower[AL_IDX]+1; ia<=nrange->upper[AL_IDX]; ++ia){
    cidx[AL_IDX] = ia;
    double alpha_curr = calc_running_coord(alpha_lo, ia-nrange->lower[AL_IDX], dalpha);
    double alpha_donor = calc_running_coord(alpha_lo, 0, dalpha);
    double alpha_diff = alpha_curr -  alpha_donor;
    for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
      cidx[PSI_IDX] = ip;
      int ip_delta_max = 3;
      for(int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++){
        double psi_curr = calc_running_coord(psi_lo, ip-nrange->lower[PSI_IDX], dpsi) + modifiers[ip_delta]*delta_psi;
        for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
          cidx[TH_IDX] = it;
            int lidx = 0;
            if (ip_delta != 0)
              lidx = 3 + 3*(ip_delta-1);

          double *mc2p_fd_n = gkyl_array_fetch(up->geo_int.mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
          double *ddtheta_n = gkyl_array_fetch(up->geo_int.ddtheta_nodal, gkyl_range_idx(nrange, cidx));
          double *ddpsi_n = gkyl_array_fetch(up->geo_int.ddpsi_nodal, gkyl_range_idx(nrange, cidx));
          double *mc2p_n = gkyl_array_fetch(up->geo_int.mc2p_nodal, gkyl_range_idx(nrange, cidx));
          double *bmag_n = gkyl_array_fetch(up->geo_int.bmag_nodal, gkyl_range_idx(nrange, cidx));
          double *curlbhat_n = gkyl_array_fetch(up->geo_int.curlbhat_nodal, gkyl_range_idx(nrange, cidx));

          int donor_cidx[3] ;
          donor_cidx[AL_IDX] = nrange->lower[AL_IDX];
          donor_cidx[PSI_IDX] = ip;
          donor_cidx[TH_IDX] = it;

          double *donor_mc2p_fd_n = gkyl_array_fetch(up->geo_int.mc2p_nodal_fd, gkyl_range_idx(nrange, donor_cidx));
          double *donor_ddtheta_n = gkyl_array_fetch(up->geo_int.ddtheta_nodal, gkyl_range_idx(nrange, donor_cidx));
          double *donor_ddpsi_n = gkyl_array_fetch(up->geo_int.ddpsi_nodal, gkyl_range_idx(nrange, donor_cidx));
          double *donor_mc2p_n = gkyl_array_fetch(up->geo_int.mc2p_nodal, gkyl_range_idx(nrange, donor_cidx));
          double *donor_bmag_n = gkyl_array_fetch(up->geo_int.bmag_nodal, gkyl_range_idx(nrange, donor_cidx));
          double *donor_curlbhat_n = gkyl_array_fetch(up->geo_int.curlbhat_nodal, gkyl_range_idx(nrange, donor_cidx));

          mc2p_fd_n[lidx+X_IDX] =  donor_mc2p_fd_n[lidx+X_IDX];
          mc2p_fd_n[lidx+Y_IDX] =  donor_mc2p_fd_n[lidx+Y_IDX];
          mc2p_fd_n[lidx+Z_IDX] =  donor_mc2p_fd_n[lidx+Z_IDX] + alpha_diff;
          if(ip_delta==0){
            ddtheta_n[0] = donor_ddtheta_n[0];
            ddtheta_n[1] = donor_ddtheta_n[1];
            ddtheta_n[2] = donor_ddtheta_n[2];
            ddpsi_n[0] = donor_ddpsi_n[0];
            mc2p_n[lidx+X_IDX] = donor_mc2p_n[lidx+X_IDX];
            mc2p_n[lidx+Y_IDX] = donor_mc2p_n[lidx+Y_IDX];
            mc2p_n[lidx+Z_IDX] = donor_mc2p_n[lidx+Z_IDX] + alpha_diff;
            bmag_n[0] = donor_bmag_n[0];
            curlbhat_func(psi_curr, mc2p_n[lidx+X_IDX], mc2p_n[lidx+Y_IDX], mc2p_n[lidx+Z_IDX], curlbhat_n, &arc_ctx);
          }
        }
      }
    }
  }

  struct gkyl_nodal_ops *n2m =  gkyl_nodal_ops_new(&inp->cbasis, &inp->cgrid, false);
  gkyl_nodal_ops_n2m(n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 3, up->geo_int.mc2p_nodal, up->geo_int.mc2p, true);
  gkyl_nodal_ops_n2m(n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 1, up->geo_int.bmag_nodal, up->geo_int.bmag, true);
  gkyl_nodal_ops_n2m(n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 1, qprofile_nodal, up->geo_int.qprofile, true);
  gkyl_nodal_ops_release(n2m);
  gkyl_array_release(qprofile_nodal);

  gkyl_free(arc_memo);
  gkyl_free(arc_memo_left);
  gkyl_free(arc_memo_right);
  gkyl_free(sep_trace_r);
  gkyl_free(sep_trace_z);
  gkyl_free(sep_trace_s);
  gkyl_free(ordered_trace_storage);
}

void gkyl_tok_geo_calc_surface(struct gk_geometry* up, int dir, struct gkyl_range *nrange, double dzc[3], 
    struct gkyl_tok_geo *geo, struct gkyl_tok_geo_grid_inp *inp, struct gkyl_position_map *position_map)
{

  geo->rleft = inp->rleft;
  geo->rright = inp->rright;

  geo->inexact_roots = inp->inexact_roots;

  geo->rmax = inp->rmax;
  geo->rmin = inp->rmin;

  enum { PSI_IDX, AL_IDX, TH_IDX }; // Arrangement of computational coordinates.
  enum { X_IDX, Y_IDX, Z_IDX }; // Arrangement of cartesian coordinates.
  
  double dtheta = inp->cgrid.dx[TH_IDX],
    dpsi = inp->cgrid.dx[PSI_IDX],
    dalpha = inp->cgrid.dx[AL_IDX];
  
  double theta_lo = up->grid.lower[TH_IDX] + (up->local.lower[TH_IDX] - up->global.lower[TH_IDX])*up->grid.dx[TH_IDX],
    psi_lo = up->grid.lower[PSI_IDX] + (up->local.lower[PSI_IDX] - up->global.lower[PSI_IDX])*up->grid.dx[PSI_IDX],
    alpha_lo = up->grid.lower[AL_IDX] + (up->local.lower[AL_IDX] - up->global.lower[AL_IDX])*up->grid.dx[AL_IDX];

  double dels[2] = {1.0/sqrt(3), 1.0-1.0/sqrt(3) };
  theta_lo += dir == 2 ? 0.0 : dels[1]*dtheta/2.0;
  psi_lo += dir == 0 ? 0.0 : dels[1]*dpsi/2.0;
  alpha_lo += dir == 1 ? 0. : dels[1]*dalpha/2.0;

    
  double dx_fact = up->basis.poly_order == 1.0/up->basis.poly_order;
  dtheta *= dx_fact; dpsi *= dx_fact; dalpha *= dx_fact;

  // Used for finite differences.
  double delta_alpha = dalpha*1e-2;
  double delta_psi = dpsi*1e-2;
  double delta_theta = dtheta*1e-2;
  dzc[0] = delta_psi;
  dzc[1] = delta_alpha;
  dzc[2] = delta_theta;
  int modifiers[5] = {0, -1, 1, -2, 2};

  double rclose = inp->rclose;
  double rright = inp->rright;
  double rleft = inp->rleft;


  int nzcells;
  if (geo->use_cubics)
    nzcells = geo->rzgrid_cubic.cells[1];
  else
    nzcells = geo->rzgrid.cells[1];
  double *arc_memo = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_left = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_right = gkyl_malloc(sizeof(double[nzcells]));
  int sep_trace_capacity = 16*nzcells+1;
  double *sep_trace_r = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *sep_trace_z = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *sep_trace_s = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *ordered_trace_storage =
    gkyl_malloc(sizeof(double[8*sep_trace_capacity]));

  struct arc_length_ctx arc_ctx = {
    .geo = geo,
    .arc_memo = arc_memo,
    .arc_memo_right = arc_memo_right,
    .arc_memo_left = arc_memo_left,
    .sep_trace_r = sep_trace_r,
    .sep_trace_z = sep_trace_z,
    .sep_trace_s = sep_trace_s,
    .sep_trace_capacity = sep_trace_capacity,
    .far_trace_r = ordered_trace_storage,
    .far_trace_z = ordered_trace_storage+sep_trace_capacity,
    .far_trace_s = ordered_trace_storage+2*sep_trace_capacity,
    .trace_corr_v = ordered_trace_storage+3*sep_trace_capacity,
    .map_trace_r = ordered_trace_storage+4*sep_trace_capacity,
    .map_trace_z = ordered_trace_storage+5*sep_trace_capacity,
    .map_trace_s = ordered_trace_storage+6*sep_trace_capacity,
    .map_trace_phi = ordered_trace_storage+7*sep_trace_capacity,
    .ftype = inp->ftype,
    .zmaxis = geo->zmaxis
  };
  struct plate_ctx pctx = {
    .geo = geo
  };
  tok_init_xpt_ray_target(inp, position_map, &arc_ctx);

  int cidx[3] = { 0 };
  for(int ia=nrange->lower[AL_IDX]; ia<nrange->lower[AL_IDX]+1; ++ia){
    cidx[AL_IDX] = ia;
    double alpha_curr = dir==1 ? alpha_lo + ia*dalpha : calc_running_coord(alpha_lo, ia-nrange->lower[AL_IDX], dalpha);

    for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
      int ip_delta_max = 5;
      for(int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++){
        if((ip == nrange->lower[PSI_IDX]) && (up->local.lower[PSI_IDX]== up->global.lower[PSI_IDX]) && dir==0){
          if(ip_delta == 1 || ip_delta == 3)
            continue; // one sided stencils at edge
        }
        else if((ip == nrange->upper[PSI_IDX]) && (up->local.upper[PSI_IDX]== up->global.upper[PSI_IDX]) && dir==0){
          if(ip_delta == 2 || ip_delta == 4)
            continue; // one sided stencils at edge
        }
        else{ // interior 
          if( ip_delta == 3 || ip_delta == 4)
            continue;
        }

        double psi_curr = dir == 0 ? psi_lo + ip*dpsi : calc_running_coord(psi_lo, ip-nrange->lower[PSI_IDX], dpsi) ;
        psi_curr += modifiers[ip_delta]*delta_psi;
        // Non-uniform psi. Finite differences are calculated in calc_metric.c
        double Psi_curr;
        position_map->maps[0](0.0, &psi_curr,  &Psi_curr,  position_map->ctxs[0]);
        double dPsi_dpsi = gkyl_position_map_slope(position_map, 0, psi_curr, delta_psi, ip, nrange);
        psi_curr = Psi_curr;
        
        double darcL, arcL_curr, arcL_lo;

        // For double null blocks this should set arc_ctx :
        // zmin, zmax, rclose, arcL_tot for all blocks. No left and right
        // For a full core case:
        // also set phi_right and arcL_right
        // For a single null case:
        // also set zmin_left and zmin_right 
        if (tok_xpt_mapping_requested(inp))
          tok_prepare_ordered_map(inp, &arc_ctx, psi_curr);
        else
          tok_find_endpoints(inp, geo, &arc_ctx, &pctx, psi_curr, alpha_curr,
            arc_memo, arc_memo_left, arc_memo_right);

        darcL = arc_ctx.arcL_tot/(up->basis.poly_order*inp->cgrid.cells[TH_IDX]) * (inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX])/2/M_PI;
        // at the beginning of each theta loop we need to reset things
        cidx[PSI_IDX] = ip;
        arcL_curr = 0.0;
        arcL_lo = (theta_lo + M_PI)/2/M_PI*arc_ctx.arcL_tot;
        double ridders_min, ridders_max;

        for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
          arcL_curr = dir==2 ? arcL_lo + it*darcL: calc_running_coord(arcL_lo, it-nrange->lower[TH_IDX], darcL);
          double theta_curr = arcL_curr*(2*M_PI/arc_ctx.arcL_tot) - M_PI ; 

          // Calculate derivatives using finite difference for ddtheta,
          // as well as transform the computational coordiante to the non-uniform field-aligned value

          // We cannot do non-uniform alpha because we are modeling axisymmetric systems
          // Non-uniform theta
          double Theta_curr;
          position_map->maps[2](0.0, &theta_curr,  &Theta_curr,  position_map->ctxs[2]);
          double dTheta_dtheta = gkyl_position_map_slope(position_map, 2, theta_curr,\
            delta_theta, it, nrange);
          theta_curr = Theta_curr;

          struct tok_ordered_point ordered = { 0.0 };
          bool ordered_mapping = tok_ordered_map_lookup(inp, &arc_ctx,
            theta_curr, alpha_curr, &ordered);
          if (tok_xpt_mapping_requested(inp) && !ordered_mapping) {
            fprintf(stderr,
              "TOK_ORDERED_MAP lookup failed ftype=%d psi=%.17g theta=%.17g\n",
              inp->ftype, psi_curr, theta_curr);
            abort();
          }
          double r_curr = 0.0, z_curr = 0.0, phi_curr = 0.0;
          double drdz_curr = 0.0, dr_curr = 0.0, dz_curr = 0.0;
          if (ordered_mapping) {
            r_curr = ordered.r; z_curr = ordered.z; phi_curr = ordered.phi;
          }
          else {
            arcL_curr = tok_xpt_theta_to_arc(inp, &arc_ctx, theta_curr);

            tok_set_ridders(inp, &arc_ctx, psi_curr, arcL_curr, &rclose, &ridders_min, &ridders_max);

          bool trace_this_block = tok_geo_trace_enabled() && tok_geo_trace_block_selected();
          if (trace_this_block && tok_geo_trace_int_env("GKYL_TOK_GEO_TRACE_BEFORE", 1))
            tok_geo_trace_surface_row("before_ridders", inp->ftype, dir, ip,
              it, ia, ip_delta, psi_curr, alpha_curr, theta_curr, arcL_curr,
              arc_ctx.zmin, arc_ctx.zmax, rclose, ridders_min, ridders_max,
              0.0, -1, 0.0);

          clock_t trace_t0 = trace_this_block ? clock() : 0;
          struct gkyl_qr_res res = gkyl_ridders(arc_length_func, &arc_ctx,
            arc_ctx.zmin, arc_ctx.zmax, ridders_min, ridders_max,
            geo->root_param.max_iter, 1e-10);
          double trace_elapsed = trace_this_block ?
            ((double) (clock() - trace_t0))/CLOCKS_PER_SEC : 0.0;
          z_curr = res.res;
          ((struct gkyl_tok_geo *)geo)->stat.nroot_cont_calls += res.nevals;

          if (trace_this_block) {
            int neval_threshold = tok_geo_trace_int_env("GKYL_TOK_GEO_TRACE_NEVAL_THRESHOLD", 25);
            double time_threshold = tok_geo_trace_double_env("GKYL_TOK_GEO_TRACE_TIME_THRESHOLD", 0.01);
            if (res.nevals >= neval_threshold || trace_elapsed >= time_threshold ||
              tok_geo_trace_int_env("GKYL_TOK_GEO_TRACE_ALL_AFTER", 0))
              tok_geo_trace_surface_row("after_ridders", inp->ftype, dir, ip,
                it, ia, ip_delta, psi_curr, alpha_curr, theta_curr, arcL_curr,
                arc_ctx.zmin, arc_ctx.zmax, rclose, ridders_min, ridders_max,
                z_curr, res.nevals, trace_elapsed);
          }

          double sep_r_curr = 0.0;
          bool at_sep_trace = tok_half_domain_sep_rz(inp, &arc_ctx,
            theta_curr, &sep_r_curr, &z_curr);

          // Only a surface normal to theta contains true theta-boundary
          // nodes.  On radial and alpha surfaces the theta coordinates are
          // Gauss points and must not be snapped to the seam.
          double fixed_root_ref = rclose;
          bool at_fixed_edge = dir == TH_IDX &&
            tok_xpt_at_fixed_edge(inp, &arc_ctx, it, nrange,
              up->local.lower[TH_IDX] == up->global.lower[TH_IDX],
              up->local.upper[TH_IDX] == up->global.upper[TH_IDX],
              &z_curr, &fixed_root_ref);
          bool at_xpt_anchor = dir == TH_IDX && arc_ctx.xpt_anchor_valid &&
            tok_xpt_at_seam(inp, it, nrange,
              up->local.lower[TH_IDX] == up->global.lower[TH_IDX],
              up->local.upper[TH_IDX] == up->global.upper[TH_IDX]);
          if (at_xpt_anchor)
            z_curr = arc_ctx.xpt_anchor_z;

          double R[4] = { 0 }, dRdZ[4] = { 0 };
          double dR[4] = { 0 }, dZ[4] = { 0 };
          int nr = gkyl_tok_geo_R_psiZ(geo, psi_curr, z_curr, 4, R, dRdZ, dR, dZ);
          double root_ref = at_sep_trace ? sep_r_curr
            : (at_xpt_anchor ? arc_ctx.xpt_anchor_r
              : (at_fixed_edge ? fixed_root_ref : rclose));
          r_curr = choose_closest(root_ref, R, R, nr);
          drdz_curr = choose_closest(root_ref, R, dRdZ, nr);
          dr_curr = choose_closest(root_ref, R, dR, nr);
          dz_curr = choose_closest(root_ref, R, dZ, nr);
          if (at_xpt_anchor)
            r_curr = arc_ctx.xpt_anchor_r;
          if (at_sep_trace)
            r_curr = sep_r_curr;

          if (tok_geo_same_flux(psi_curr, geo->psisep) && ip_delta==0) {
            if (z_curr == geo->efit->Zxpt[0]) {
              nr = 1;
              r_curr = geo->efit->Rxpt[0];
            }
            if (z_curr == geo->efit->Zxpt[1]) {
              nr = 1;
              r_curr = geo->efit->Rxpt[1];
            }
          }

          if(nr==0){
            printf("ip = %d, it = %d, ia = %d, ip_delta = %d\n", ip, it, ia, ip_delta);
            printf("Block Type = %d | Failed to find a root at psi = %g, Z = %1.16f\n", inp->ftype, psi_curr, z_curr);
            assert(false);
          }

            phi_curr = phi_func(alpha_curr, z_curr, &arc_ctx);
          }
          cidx[TH_IDX] = it;
          int lidx = 0;
          if (ip_delta != 0)
            lidx = 3 + 3*(ip_delta-1);
          double *mc2p_fd_n = gkyl_array_fetch(up->geo_surf[dir].mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
          double *ddtheta_n = gkyl_array_fetch(up->geo_surf[dir].ddtheta_nodal, gkyl_range_idx(nrange, cidx));
          double *ddpsi_n = gkyl_array_fetch(up->geo_surf[dir].ddpsi_nodal, gkyl_range_idx(nrange, cidx));
          double *bmag_n = gkyl_array_fetch(up->geo_surf[dir].bmag_nodal, gkyl_range_idx(nrange, cidx));
          double *curlbhat_n = gkyl_array_fetch(up->geo_surf[dir].curlbhat_nodal, gkyl_range_idx(nrange, cidx));
          double *deltats_n = gkyl_array_fetch(up->geo_surf[dir].deltats_nodal, gkyl_range_idx(nrange, cidx));

          mc2p_fd_n[lidx+X_IDX] = r_curr;
          mc2p_fd_n[lidx+Y_IDX] = z_curr;
          mc2p_fd_n[lidx+Z_IDX] = phi_curr;

          if(ip_delta==0){
            if (ordered_mapping) {
              ddtheta_n[0] = ordered.dr_dtheta*dTheta_dtheta;
              ddtheta_n[1] = ordered.dz_dtheta*dTheta_dtheta;
              ddtheta_n[2] = ordered.dphi_dtheta*dTheta_dtheta;
            }
            else {
              double darc_dtheta = arc_ctx.xpt_map_valid
                ? arc_ctx.xpt_map_darc_dtheta : arc_ctx.arcL_tot/(2.0*M_PI);
              ddtheta_n[0] = sin(atan2(dr_curr,dz_curr))*darc_dtheta*dTheta_dtheta;
              ddtheta_n[1] = cos(atan2(dr_curr,dz_curr))*darc_dtheta*dTheta_dtheta;
              ddtheta_n[2] = dphidtheta_func(z_curr, &arc_ctx)*dTheta_dtheta;
            }
            ddpsi_n[0] = dPsi_dpsi;
            bmag_n[0] = bmag_func(r_curr, z_curr, &arc_ctx);
            curlbhat_func(psi_curr, r_curr, z_curr, phi_curr, curlbhat_n, &arc_ctx);
            deltats_n[0] = phi_curr - alpha_curr;
          }
        }
      }
    }
  }

  // Populate other alpha indices by using axisymmetry
  for (int ia=nrange->lower[AL_IDX]+1; ia<=nrange->upper[AL_IDX]; ++ia){
    cidx[AL_IDX] = ia;
    double alpha_curr = dir==1 ? alpha_lo + ia*dalpha : calc_running_coord(alpha_lo, ia-nrange->lower[AL_IDX], dalpha);
    double alpha_donor= dir==1 ? alpha_lo + nrange->lower[AL_IDX]*dalpha : calc_running_coord(alpha_lo, 0, dalpha);
    double alpha_diff = alpha_curr -  alpha_donor;

    for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
      cidx[PSI_IDX] = ip;
      int ip_delta_max = 5;
      for(int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++){
        if((ip == nrange->lower[PSI_IDX]) && (up->local.lower[PSI_IDX]== up->global.lower[PSI_IDX]) && dir==0){
          if(ip_delta == 1 || ip_delta == 3)
            continue; // one sided stencils at edge
        }
        else if((ip == nrange->upper[PSI_IDX]) && (up->local.upper[PSI_IDX]== up->global.upper[PSI_IDX]) && dir==0){
          if(ip_delta == 2 || ip_delta == 4)
            continue; // one sided stencils at edge
        }
        else{ // interior 
          if( ip_delta == 3 || ip_delta == 4)
            continue;
        }
        double psi_curr = dir == 0 ? psi_lo + ip*dpsi : calc_running_coord(psi_lo, ip-nrange->lower[PSI_IDX], dpsi) ;
        psi_curr += modifiers[ip_delta]*delta_psi;

        for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
          cidx[TH_IDX] = it;
            int lidx = 0;
            if (ip_delta != 0)
              lidx = 3 + 3*(ip_delta-1);

          double *mc2p_fd_n = gkyl_array_fetch(up->geo_surf[dir].mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
          double *ddtheta_n = gkyl_array_fetch(up->geo_surf[dir].ddtheta_nodal, gkyl_range_idx(nrange, cidx));
          double *ddpsi_n = gkyl_array_fetch(up->geo_surf[dir].ddpsi_nodal, gkyl_range_idx(nrange, cidx));
          double *bmag_n = gkyl_array_fetch(up->geo_surf[dir].bmag_nodal, gkyl_range_idx(nrange, cidx));
          double *curlbhat_n = gkyl_array_fetch(up->geo_surf[dir].curlbhat_nodal, gkyl_range_idx(nrange, cidx));
          double *deltats_n = gkyl_array_fetch(up->geo_surf[dir].deltats_nodal, gkyl_range_idx(nrange, cidx));

          int donor_cidx[3] ;
          donor_cidx[AL_IDX] = nrange->lower[AL_IDX];
          donor_cidx[PSI_IDX] = ip;
          donor_cidx[TH_IDX] = it;

          double *donor_mc2p_fd_n = gkyl_array_fetch(up->geo_surf[dir].mc2p_nodal_fd, gkyl_range_idx(nrange, donor_cidx));
          double *donor_ddtheta_n = gkyl_array_fetch(up->geo_surf[dir].ddtheta_nodal, gkyl_range_idx(nrange, donor_cidx));
          double *donor_ddpsi_n = gkyl_array_fetch(up->geo_surf[dir].ddpsi_nodal, gkyl_range_idx(nrange, donor_cidx));
          double *donor_bmag_n = gkyl_array_fetch(up->geo_surf[dir].bmag_nodal, gkyl_range_idx(nrange, donor_cidx));

          mc2p_fd_n[lidx+X_IDX] =  donor_mc2p_fd_n[lidx+X_IDX];
          mc2p_fd_n[lidx+Y_IDX] =  donor_mc2p_fd_n[lidx+Y_IDX];
          mc2p_fd_n[lidx+Z_IDX] =  donor_mc2p_fd_n[lidx+Z_IDX] + alpha_diff;
          if(ip_delta==0){
            ddtheta_n[0] = donor_ddtheta_n[0];
            ddtheta_n[1] = donor_ddtheta_n[1];
            ddtheta_n[2] = donor_ddtheta_n[2];
            ddpsi_n[0] = donor_ddpsi_n[0];
            bmag_n[0] = donor_bmag_n[0];
            curlbhat_func(psi_curr, mc2p_fd_n[X_IDX], mc2p_fd_n[Y_IDX], mc2p_fd_n[Z_IDX], curlbhat_n, &arc_ctx);
            deltats_n[0] = mc2p_fd_n[Z_IDX] - alpha_curr;
          }
        }
      }
    }
  }

  gkyl_free(arc_memo);
  gkyl_free(arc_memo_left);
  gkyl_free(arc_memo_right);
  gkyl_free(sep_trace_r);
  gkyl_free(sep_trace_z);
  gkyl_free(sep_trace_s);
  gkyl_free(ordered_trace_storage);
}


void
gkyl_tok_geo_set_extent(struct gkyl_tok_geo_grid_inp* inp, struct gkyl_tok_geo *geo, double *theta_lo, double *theta_up)
{
  tok_geo_set_extent(inp, geo, theta_lo, theta_up);
}

struct gkyl_tok_geo_stat
gkyl_tok_geo_get_stat(const struct gkyl_tok_geo *geo)
{
  return geo->stat;
}

void
gkyl_tok_geo_release(struct gkyl_tok_geo *geo)
{
  gkyl_array_release(geo->psiRZ);
  gkyl_array_release(geo->psiRZ_cubic);
  gkyl_array_release(geo->fpoldg);
  gkyl_array_release(geo->fpolprimedg);
  gkyl_array_release(geo->qdg);
  gkyl_efit_release(geo->efit);
  gkyl_free(geo);
}
