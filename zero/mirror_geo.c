#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <assert.h>
#include <gkyl_basis.h>
#include <gkyl_calc_bmag.h>
#include <gkyl_efit.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_math.h>
#include <gkyl_mirror_geo_priv.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_position_map.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

#include <math.h>
#include <string.h>

double
mirror_plate_psi_func(double s, void *ctx){
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
  int rzidx[2];

  if (gc->geo->use_cubics) {
    double xn[2] = {R, Z};
    double psi;
    gc->geo->efit->evf->eval_cubic(0.0, xn, &psi, gc->geo->efit->evf->ctx);
    return psi - gc->psi_curr;
  }
  else {
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
  double *arc_memo = actx->arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL;
  double zmax = actx->zmax;
  double ival = 0.0;
  ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, false, arc_memo) - arcL;
  return ival;
}

struct gkyl_mirror_geo*
gkyl_mirror_geo_new(const struct gkyl_efit_inp *inp, const struct gkyl_mirror_geo_grid_inp *ginp)
{
  struct gkyl_mirror_geo *geo = gkyl_malloc(sizeof(*geo));
  *geo = (struct gkyl_mirror_geo) {};

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
  geo->qdg= gkyl_array_acquire(geo->efit->qflux);
  geo->sibry= geo->efit->sibry;
  geo->psisep = geo->efit->psisep;
  geo->zmaxis = geo->efit->zmaxis;

  geo->use_cubics = ginp->use_cubics;
  geo->root_param.eps =
    ginp->root_param.eps > 0 ? ginp->root_param.eps : 1e-10;
  geo->root_param.max_iter =
    ginp->root_param.max_iter > 0 ? ginp->root_param.max_iter : 100;

  geo->quad_param.max_level =
    ginp->quad_param.max_levels > 0 ? ginp->quad_param.max_levels : 10;
  geo->quad_param.eps =
    ginp->quad_param.eps > 0 ? ginp->quad_param.eps : 1e-10;

  if (geo->use_cubics) {
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

  geo->stat = (struct gkyl_mirror_geo_stat) { };

  
  return geo;
}

double
gkyl_mirror_geo_integrate_psi_contour(const struct gkyl_mirror_geo *geo, double psi,
  double zmin, double zmax, double rclose)
{
  return integrate_psi_contour_memo(geo, psi, zmin, zmax, rclose,
    false, false, 0);
}

int
gkyl_mirror_geo_R_psiZ(const struct gkyl_mirror_geo *geo, double psi, double Z, int nmaxroots,
  double *R, double *dR)
{
  if(geo->use_cubics)
    return R_psiZ_cubic(geo, psi, Z, nmaxroots, R, dR);
  else
    return R_psiZ(geo, psi, Z, nmaxroots, R, dR);
}

void gkyl_mirror_geo_calc(struct gk_geometry* up, struct gkyl_range *nrange, double dzc[3], 
  struct gkyl_mirror_geo *geo, struct gkyl_mirror_geo_grid_inp *inp,
  struct gkyl_array *mc2p_nodal_fd, struct gkyl_array *mc2p_nodal, struct gkyl_array *mc2p, struct gkyl_array *ddtheta_nodal,
  struct gkyl_array *mc2nu_pos_nodal, struct gkyl_array *mc2nu_pos,
  struct gkyl_position_map *position_map)
{


  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  
  double dtheta = inp->cgrid.dx[TH_IDX],
    dpsi = inp->cgrid.dx[PSI_IDX],
    dalpha = inp->cgrid.dx[AL_IDX];
  
  double theta_lo = up->grid.lower[TH_IDX] + (up->local.lower[TH_IDX] - up->global.lower[TH_IDX])*up->grid.dx[TH_IDX],
    psi_lo = up->grid.lower[PSI_IDX] + (up->local.lower[PSI_IDX] - up->global.lower[PSI_IDX])*up->grid.dx[PSI_IDX],
    alpha_lo = up->grid.lower[AL_IDX] + (up->local.lower[AL_IDX] - up->global.lower[AL_IDX])*up->grid.dx[AL_IDX];

  double dx_fact = up->basis.poly_order == 1 ? 1 : 0.5;
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

  int nzcells;
  if(geo->use_cubics)
    nzcells = geo->rzgrid_cubic.cells[1];
  else
    nzcells = geo->rzgrid.cells[1];
  double *arc_memo = gkyl_malloc(sizeof(double[nzcells]));

  struct arc_length_ctx arc_ctx = {
    .geo = geo,
    .arc_memo = arc_memo,
    .zmaxis = geo->zmaxis,
  };
  struct plate_ctx pctx = {
    .geo = geo
  };

  position_map->constB_ctx->psi_max   = up->grid.upper[PSI_IDX];
  position_map->constB_ctx->psi_min   = up->grid.lower[PSI_IDX];
  position_map->constB_ctx->alpha_max = up->grid.upper[AL_IDX];
  position_map->constB_ctx->alpha_min = up->grid.lower[AL_IDX];
  position_map->constB_ctx->theta_max = up->grid.upper[TH_IDX];
  position_map->constB_ctx->theta_min = up->grid.lower[TH_IDX];
  position_map->constB_ctx->N_theta_boundaries = up->global.upper[TH_IDX] - up->global.lower[TH_IDX];
  gkyl_position_map_optimize(position_map);

  int cidx[3] = { 0 };
  for(int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia){
    cidx[AL_IDX] = ia;
    int ia_delta = 0;
    double alpha_curr = alpha_lo + ia*dalpha + modifiers[ia_delta]*delta_alpha;
    // This is the convention described in Noah Mandell's Thesis Eq 5.104. comp coord y = -alpha.
    alpha_curr*=-1.0;

    for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
      int ip_delta_max = 5;
      if(ia_delta != 0)
        ip_delta_max = 1;
      for(int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++){
        if((ip == nrange->lower[PSI_IDX]) && (up->local.lower[PSI_IDX]== up->global.lower[PSI_IDX]) ){
          if(ip_delta == 1 || ip_delta == 3)
            continue; // one sided stencils at edge
        }
        else if((ip == nrange->upper[PSI_IDX]) && (up->local.upper[PSI_IDX]== up->global.upper[PSI_IDX])){
          if(ip_delta == 2 || ip_delta == 4)
            continue; // one sided stencils at edge
        }
        else{ // interior 
          if( ip_delta == 3 || ip_delta == 4)
            continue;
        }

        double psi_curr = psi_lo + ip*dpsi + modifiers[ip_delta]*delta_psi;
        double zmin = inp->zmin, zmax = inp->zmax;
        double darcL, arcL_curr, arcL_lo;

        mirror_find_endpoints(inp, geo, &arc_ctx, &pctx, psi_curr, alpha_curr, &zmin, &zmax, arc_memo);
        double arcL = arc_ctx.arcL_tot;
        darcL = arcL/(up->basis.poly_order*inp->cgrid.cells[TH_IDX]) * (inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX])/2/M_PI;

        // at the beginning of each theta loop we need to reset things
        cidx[PSI_IDX] = ip;
        arcL_curr = 0.0;
        arcL_lo = (theta_lo + M_PI)/2/M_PI*arcL;
        double ridders_min, ridders_max;
        // set node coordinates
        for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
          int it_delta = 0;
          arcL_curr = arcL_lo + it*darcL + modifiers[it_delta]*delta_theta*(arcL/2/M_PI);
          double theta_curr = arcL_curr*(2*M_PI/arcL) - M_PI ; 

          position_map->maps[0](0.0, &psi_curr,   &psi_curr,   position_map->ctxs[0]);
          position_map->maps[1](0.0, &alpha_curr, &alpha_curr, position_map->ctxs[1]);
          position_map->maps[2](0.0, &theta_curr, &theta_curr, position_map->ctxs[2]);
          arcL_curr = (theta_curr + M_PI) / (2*M_PI/arc_ctx.arcL_tot);

          mirror_set_ridders(inp, &arc_ctx, psi_curr, arcL, arcL_curr, zmin, zmax, &rclose, &ridders_min, &ridders_max);

          struct gkyl_qr_res res = gkyl_ridders(arc_length_func, &arc_ctx,
            arc_ctx.zmin, arc_ctx.zmax, ridders_min, ridders_max,
            geo->root_param.max_iter, 1e-10);
          double z_curr = res.res;
          ((struct gkyl_mirror_geo *)geo)->stat.nroot_cont_calls += res.nevals;
          double R[4] = { 0 }, dR[4] = { 0 };
          int nr = gkyl_mirror_geo_R_psiZ(geo, psi_curr, z_curr, 4, R, dR);
          double r_curr = choose_closest(rclose, R, R, nr);
          double dr_curr = choose_closest(rclose, R, dR, nr);


          if(nr==0){
            printf(" ip = %d, it = %d, ia = %d, ip_delta = %d, it_delta = %d, ia_delta = %d\n", ip, it, ia, ip_delta, it_delta, ia_delta);
            printf("Failed to find a root at psi = %g, Z = %1.16f\n", psi_curr, z_curr);
            assert(false);
          }

          cidx[TH_IDX] = it;
          int lidx = 0;
          if (ip_delta != 0)
            lidx = 3 + 3*(ip_delta-1);
          if (ia_delta != 0)
            lidx = 15 + 3*(ia_delta-1);
          if (it_delta != 0)
            lidx = 27 + 3*(it_delta-1);

          if(ia_delta==0 && ip_delta==0 && it_delta==0)
            lidx = 0;

          double phi_curr = alpha_curr;
          double *mc2p_fd_n = gkyl_array_fetch(mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
          double *mc2p_n = gkyl_array_fetch(mc2p_nodal, gkyl_range_idx(nrange, cidx));
          double *mc2nu_n = gkyl_array_fetch(mc2nu_pos_nodal, gkyl_range_idx(nrange, cidx));
          double *ddtheta_n = gkyl_array_fetch(ddtheta_nodal, gkyl_range_idx(nrange, cidx));

          mc2p_fd_n[lidx+X_IDX] = r_curr;
          mc2p_fd_n[lidx+Y_IDX] = z_curr;
          mc2p_fd_n[lidx+Z_IDX] = phi_curr;

          if(ip_delta==0 && ia_delta==0 && it_delta==0){
            mc2p_n[X_IDX] = r_curr;
            mc2p_n[Y_IDX] = z_curr;
            mc2p_n[Z_IDX] = phi_curr;
            mc2nu_n[X_IDX] = psi_curr;
            mc2nu_n[Y_IDX] = -alpha_curr;
            mc2nu_n[Z_IDX] = theta_curr;
            ddtheta_n[0] = 0.0;
            ddtheta_n[1] = sin(atan(dr_curr))*arc_ctx.arcL_tot/2.0/M_PI;
            ddtheta_n[2] = cos(atan(dr_curr))*arc_ctx.arcL_tot/2.0/M_PI;
          }
        }
      }
    }
  }

  struct gkyl_nodal_ops *n2m =  gkyl_nodal_ops_new(&inp->cbasis, &inp->cgrid, false);
  gkyl_nodal_ops_n2m(n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 3, mc2p_nodal, mc2p);
  gkyl_nodal_ops_n2m(n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 3, mc2nu_pos_nodal, mc2nu_pos);
  gkyl_nodal_ops_release(n2m);

  gkyl_free(arc_memo);
}


struct gkyl_mirror_geo_stat
gkyl_mirror_geo_get_stat(const struct gkyl_mirror_geo *geo)
{
  return geo->stat;
}

void
gkyl_mirror_geo_release(struct gkyl_mirror_geo *geo)
{
  gkyl_array_release(geo->psiRZ);
  gkyl_array_release(geo->psiRZ_cubic);
  gkyl_array_release(geo->fpoldg);
  gkyl_array_release(geo->qdg);
  gkyl_efit_release(geo->efit);
  gkyl_free(geo);
}
