#include <gkyl_calc_metric.h>
#include <gkyl_calc_metric_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_array_ops_priv.h>

gkyl_calc_metric*
gkyl_calc_metric_new(const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid,
  const struct gkyl_range *global, const struct gkyl_range *global_ext,
  const struct gkyl_range *local, const struct gkyl_range *local_ext, bool use_gpu)
{
  gkyl_calc_metric *up = gkyl_malloc(sizeof(gkyl_calc_metric));
  up->cbasis = cbasis;
  up->cdim = cbasis->ndim;
  up->cnum_basis = cbasis->num_basis;
  up->poly_order = cbasis->poly_order;
  up->grid = grid;
  up->use_gpu = use_gpu;
  up->num_cells = up->grid->cells;
  up->n2m = gkyl_nodal_ops_new(up->cbasis, up->grid, up->use_gpu);

  up->global = *global;
  up->global_ext = *global_ext;

  up->local = *local;
  up->local_ext = *local_ext;
  return up;
}

static inline double calc_metric(double dxdz[3][3], int i, int j) 
{
  double sum = 0;
  for (int k = 0; k < 3; ++k)
    sum += dxdz[k][i - 1] * dxdz[k][j - 1];
  return sum;
} 

// Calculates e^1 = e_2 x e_3 /J
static inline void
calc_dual(double J, const double e_2[3], const double e_3[3], double e1[3])
{
  e1[0] = (e_2[1]*e_3[2] - e_2[2]*e_3[1] )/J;
  e1[1] = -(e_2[0]*e_3[2] - e_2[2]*e_3[0] )/J;
  e1[2] = (e_2[0]*e_3[1] - e_2[1]*e_3[0] )/J;
}

void gkyl_calc_metric_advance_rz(
  gkyl_calc_metric *up, struct gkyl_range *nrange,
  struct gkyl_array *mc2p_nodal_fd, struct gkyl_array *ddtheta_nodal,
  struct gkyl_array *bmag_nodal, double *dzc, struct gkyl_array *gFld,
  struct gkyl_array *tanvecFld,
  struct gkyl_array *dualFld,
  struct gkyl_array *dualmagFld,
  struct gkyl_array *normFld,
  struct gkyl_array *jFld, struct gkyl_array* bcartFld, const struct gkyl_range *update_range)
{
  struct gkyl_array* gFld_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nrange->volume);
  struct gkyl_array* jFld_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nrange->volume);
  struct gkyl_array* bcartFld_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange->volume);
  struct gkyl_array* tanvecFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  struct gkyl_array* dualFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  struct gkyl_array* dualmagFld_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange->volume);
  struct gkyl_array* normFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { R_IDX, Z_IDX, PHI_IDX }; // arrangement of cartesian coordinates
  int cidx[3];
  for(int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia){
      for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
          for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
              cidx[PSI_IDX] = ip;
              cidx[AL_IDX] = ia;
              cidx[TH_IDX] = it;
              const double *mc2p_n = gkyl_array_cfetch(mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
              double dxdz[3][3];

              if((ip == nrange->lower[PSI_IDX]) && (up->local.lower[PSI_IDX]== up->global.lower[PSI_IDX]) ) {
                dxdz[0][0] = (-3*mc2p_n[R_IDX] + 4*mc2p_n[6+R_IDX] - mc2p_n[12+R_IDX] )/dzc[0]/2;
                dxdz[1][0] = (-3*mc2p_n[Z_IDX] + 4*mc2p_n[6+Z_IDX] - mc2p_n[12+Z_IDX] )/dzc[0]/2;
                dxdz[2][0] = (-3*mc2p_n[PHI_IDX] + 4*mc2p_n[6+PHI_IDX] - mc2p_n[12+PHI_IDX] );
              }
              else if((ip == nrange->upper[PSI_IDX]) && (up->local.upper[PSI_IDX]== up->global.upper[PSI_IDX])) {
                dxdz[0][0] = (3*mc2p_n[R_IDX] - 4*mc2p_n[3+R_IDX] + mc2p_n[9+R_IDX] )/dzc[0]/2;
                dxdz[1][0] = (3*mc2p_n[Z_IDX] - 4*mc2p_n[3+Z_IDX] + mc2p_n[9+Z_IDX] )/dzc[0]/2;
                dxdz[2][0] = (3*mc2p_n[PHI_IDX] - 4*mc2p_n[3+PHI_IDX] + mc2p_n[9+PHI_IDX] );
              }
              else {
                dxdz[0][0] = -(mc2p_n[3 +R_IDX] -   mc2p_n[6+R_IDX])/2/dzc[0];
                dxdz[1][0] = -(mc2p_n[3 +Z_IDX] -   mc2p_n[6+Z_IDX])/2/dzc[0];
                dxdz[2][0] = -(mc2p_n[3 +PHI_IDX] -   mc2p_n[6+PHI_IDX]);
              }
              // Take into account wrapping of cyclic coordinate phi
              if (dxdz[2][0] < -M_PI) {
                dxdz[2][0] += 2*M_PI;
              } else if (dxdz[2][0] > M_PI) {
                dxdz[2][0] -= 2*M_PI;
              }
              dxdz[2][0] = dxdz[2][0]/2.0/dzc[0];

              if((ia == nrange->lower[AL_IDX]) && (up->local.lower[AL_IDX]== up->global.lower[AL_IDX]) ) {
                dxdz[0][1] = (-3*mc2p_n[R_IDX] +  4*mc2p_n[18+R_IDX] -  mc2p_n[24+R_IDX])/dzc[1]/2;
                dxdz[1][1] = (-3*mc2p_n[Z_IDX] +  4*mc2p_n[18+Z_IDX] -  mc2p_n[24+Z_IDX])/dzc[1]/2;
                dxdz[2][1] = (-3*mc2p_n[PHI_IDX] +  4*mc2p_n[18+PHI_IDX] -  mc2p_n[24+PHI_IDX])/dzc[1]/2;
              }
              else if((ia == nrange->upper[AL_IDX])  && (up->local.upper[AL_IDX]== up->global.upper[AL_IDX])){
                dxdz[0][1] = (3*mc2p_n[R_IDX] -  4*mc2p_n[15+R_IDX] +  mc2p_n[21+R_IDX] )/dzc[1]/2;
                dxdz[1][1] = (3*mc2p_n[Z_IDX] -  4*mc2p_n[15+Z_IDX] +  mc2p_n[21+Z_IDX] )/dzc[1]/2;
                dxdz[2][1] = (3*mc2p_n[PHI_IDX] -  4*mc2p_n[15+PHI_IDX] +  mc2p_n[21+PHI_IDX] )/dzc[1]/2;
              }
              else {
                dxdz[0][1] = -(mc2p_n[15 +R_IDX] -  mc2p_n[18 +R_IDX])/2/dzc[1];
                dxdz[1][1] = -(mc2p_n[15 +Z_IDX] -  mc2p_n[18 +Z_IDX])/2/dzc[1];
                dxdz[2][1] = -(mc2p_n[15 +PHI_IDX] -  mc2p_n[18 +PHI_IDX])/2/dzc[1];
              }

              if((it == nrange->lower[TH_IDX]) && (up->local.lower[TH_IDX]== up->global.lower[TH_IDX])){
                dxdz[0][2] = (-3*mc2p_n[R_IDX] + 4*mc2p_n[30+R_IDX] - mc2p_n[36+R_IDX])/dzc[2]/2;
                dxdz[1][2] = (-3*mc2p_n[Z_IDX] + 4*mc2p_n[30+Z_IDX] - mc2p_n[36+Z_IDX])/dzc[2]/2;
                dxdz[2][2] = (-3*mc2p_n[PHI_IDX] + 4*mc2p_n[30+PHI_IDX] - mc2p_n[36+PHI_IDX])/dzc[2]/2;
              }
              else if((it == nrange->upper[TH_IDX]) && (up->local.upper[TH_IDX]== up->global.upper[TH_IDX])){
                dxdz[0][2] = (3*mc2p_n[R_IDX] - 4*mc2p_n[27+R_IDX] + mc2p_n[33+R_IDX] )/dzc[2]/2;
                dxdz[1][2] = (3*mc2p_n[Z_IDX] - 4*mc2p_n[27+Z_IDX] + mc2p_n[33+Z_IDX] )/dzc[2]/2;
                dxdz[2][2] = (3*mc2p_n[PHI_IDX] - 4*mc2p_n[27+PHI_IDX] + mc2p_n[33+PHI_IDX] )/dzc[2]/2;
              }
              else {
                dxdz[0][2] = -(mc2p_n[27 +R_IDX] - mc2p_n[30 +R_IDX])/2/dzc[2];
                dxdz[1][2] = -(mc2p_n[27 +Z_IDX] - mc2p_n[30 +Z_IDX])/2/dzc[2];
                dxdz[2][2] = -(mc2p_n[27 +PHI_IDX] - mc2p_n[30 +PHI_IDX])/2/dzc[2];
              }

              // Use exact expressions for dR/dtheta and dZ/dtheta
              double *ddtheta_n = gkyl_array_fetch(ddtheta_nodal, gkyl_range_idx(nrange, cidx));
              dxdz[0][2] = ddtheta_n[0]; // dR/dtheta
              dxdz[1][2] = ddtheta_n[1]; // dZ/dtheta

              // use exact expressions for d/dalpha
              dxdz[0][1] = 0.0; // dR/dalpha
              dxdz[1][1] = 0.0; // dZ/dalpha
              dxdz[2][1] = -1.0; // dphi/dalpha

              // dxdz is in cylindrical coords, calculate J as
              // J = R(dR/dpsi*dZ/dtheta - dR/dtheta*dZ/dpsi)
              double *jFld_n= gkyl_array_fetch(jFld_nodal, gkyl_range_idx(nrange, cidx));
              double R = mc2p_n[R_IDX];
              jFld_n[0] = sqrt(R*R*(dxdz[0][0]*dxdz[0][0]*dxdz[1][2]*dxdz[1][2] + dxdz[0][2]*dxdz[0][2]*dxdz[1][0]*dxdz[1][0] - 2*dxdz[0][0]*dxdz[0][2]*dxdz[1][0]*dxdz[1][2])) ;

              // Calculate dphi/dtheta based on the divergence free condition
              // on B: 1 = J*B/sqrt(g_33)
              double *bmag_n = gkyl_array_fetch(bmag_nodal, gkyl_range_idx(nrange, cidx));
              double dphidtheta = (jFld_n[0]*jFld_n[0]*bmag_n[0]*bmag_n[0] - dxdz[0][2]*dxdz[0][2] - dxdz[1][2]*dxdz[1][2])/R/R;
              dphidtheta = sqrt(dphidtheta);
              // recover sign from exact dphidtheta = F(psi)/R/\grad(psi)Add commentMore actions
              if (ddtheta_n[2] < 0) {
                dphidtheta = -dphidtheta;
              }

              double *gFld_n= gkyl_array_fetch(gFld_nodal, gkyl_range_idx(nrange, cidx));
              gFld_n[0] = dxdz[0][0]*dxdz[0][0] + R*R*dxdz[2][0]*dxdz[2][0] + dxdz[1][0]*dxdz[1][0]; 
              gFld_n[1] = R*R*dxdz[2][0]; 
              gFld_n[2] = dxdz[0][0]*dxdz[0][2] + R*R*dxdz[2][0]*dphidtheta + dxdz[1][0]*dxdz[1][2];
              gFld_n[3] = R*R; 
              gFld_n[4] = R*R*dphidtheta;
              gFld_n[5] = dxdz[0][2]*dxdz[0][2] + R*R*dphidtheta*dphidtheta + dxdz[1][2]*dxdz[1][2]; 

              // Calculate cartesian components of bhat
              double *bcartFld_n= gkyl_array_fetch(bcartFld_nodal, gkyl_range_idx(nrange, cidx));
              double phi = mc2p_n[PHI_IDX];
              double b3 = 1/sqrt(gFld_n[5]);
              bcartFld_n[0] = b3*(dxdz[0][2]*cos(phi) - R*sin(phi)*dphidtheta);
              bcartFld_n[1] = b3*(dxdz[0][2]*sin(phi) + R*cos(phi)*dphidtheta);
              bcartFld_n[2] = b3*(dxdz[1][2]);

              // Set cartesian components of tangents and duals
              double Z = mc2p_n[Z_IDX];
              double J = jFld_n[0];
              double *tanvecFld_n= gkyl_array_fetch(tanvecFld_nodal, gkyl_range_idx(nrange, cidx));
              tanvecFld_n[0] = dxdz[0][0]*cos(phi) - R*sin(phi)*dxdz[2][0]; 
              tanvecFld_n[1] = dxdz[0][0]*sin(phi)  + R*cos(phi)*dxdz[2][0]; 
              tanvecFld_n[2] = dxdz[1][0];

              tanvecFld_n[3] = +R*sin(phi); 
              tanvecFld_n[4] = -R*cos(phi); 
              tanvecFld_n[5] = 0.0; 

              tanvecFld_n[6] = dxdz[0][2]*cos(phi) - R*sin(phi)*dphidtheta; 
              tanvecFld_n[7] = dxdz[0][2]*sin(phi)  + R*cos(phi)*dphidtheta; 
              tanvecFld_n[8] = dxdz[1][2];

              double *dualFld_n= gkyl_array_fetch(dualFld_nodal, gkyl_range_idx(nrange, cidx));
              dualFld_n[0] = -R/J*cos(phi)*dxdz[1][2];
              dualFld_n[1] = -R/J*sin(phi)*dxdz[1][2];
              dualFld_n[2] = +R/J*dxdz[0][2];

              dualFld_n[3] = 1/J * (dxdz[1][0]*dxdz[0][2]*sin(phi) + dxdz[1][0]*R*cos(phi)*dphidtheta - dxdz[1][2]*dxdz[0][0]*sin(phi) - dxdz[1][2]*R*cos(phi)*dxdz[2][0] );
              dualFld_n[4] = -1/J * (dxdz[1][0]*dxdz[0][2]*cos(phi) + dxdz[1][0]*R*sin(phi)*dphidtheta - dxdz[1][2]*dxdz[0][0]*cos(phi) - dxdz[1][2]*R*sin(phi)*dxdz[2][0] );
              dualFld_n[5] =  R/J * ( dxdz[0][2]*dxdz[2][0] - dxdz[0][0]*dphidtheta);

              dualFld_n[6] = +R/J*cos(phi)*dxdz[1][0];
              dualFld_n[7] = +R/J*sin(phi)*dxdz[1][0];
              dualFld_n[8] = -R/J*dxdz[0][0];

              double norm1 = sqrt(dualFld_n[0]*dualFld_n[0] + dualFld_n[1]*dualFld_n[1] + dualFld_n[2]*dualFld_n[2]);
              double norm2 = sqrt(dualFld_n[3]*dualFld_n[3] + dualFld_n[4]*dualFld_n[4] + dualFld_n[5]*dualFld_n[5]);
              double norm3 = sqrt(dualFld_n[6]*dualFld_n[6] + dualFld_n[7]*dualFld_n[7] + dualFld_n[8]*dualFld_n[8]);

              double *dualmagFld_n = gkyl_array_fetch(dualmagFld_nodal, gkyl_range_idx(nrange, cidx));
              dualmagFld_n[0] = norm1;
              dualmagFld_n[1] = norm2;
              dualmagFld_n[2] = norm3;
              
              // Set normal vectors
              double *normFld_n = gkyl_array_fetch(normFld_nodal, gkyl_range_idx(nrange, cidx));
              normFld_n[0] = dualFld_n[0]/norm1;
              normFld_n[1] = dualFld_n[1]/norm1;
              normFld_n[2] = dualFld_n[2]/norm1;

              normFld_n[3] = dualFld_n[3]/norm2;
              normFld_n[4] = dualFld_n[4]/norm2;
              normFld_n[5] = dualFld_n[5]/norm2;

              normFld_n[6] = dualFld_n[6]/norm3;
              normFld_n[7] = dualFld_n[7]/norm3;
              normFld_n[8] = dualFld_n[8]/norm3;
      }
    }
  }
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 6, gFld_nodal, gFld);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 1, jFld_nodal, jFld);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 3, bcartFld_nodal, bcartFld);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, tanvecFld_nodal, tanvecFld);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, dualFld_nodal, dualFld);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 3, dualmagFld_nodal, dualmagFld);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, normFld_nodal, normFld);
  gkyl_array_release(gFld_nodal);
  gkyl_array_release(jFld_nodal);
  gkyl_array_release(bcartFld_nodal);
  gkyl_array_release(tanvecFld_nodal);
  gkyl_array_release(dualFld_nodal);
  gkyl_array_release(dualmagFld_nodal);
  gkyl_array_release(normFld_nodal);
}

void gkyl_calc_metric_advance_rz_neut(
  gkyl_calc_metric *up, struct gkyl_range *nrange,
  struct gkyl_array *mc2p_nodal_fd, struct gkyl_array *ddtheta_nodal,
  double *dzc, struct gkyl_array *gFld, struct gkyl_array *grFld,
  const struct gkyl_range *update_range)
{
  struct gkyl_array* gFld_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nrange->volume);
  struct gkyl_array* grFld_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nrange->volume);
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { R_IDX, Z_IDX, PHI_IDX }; // arrangement of cartesian coordinates
  int cidx[3];
  for(int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia){
      for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
          for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
              cidx[PSI_IDX] = ip;
              cidx[AL_IDX] = ia;
              cidx[TH_IDX] = it;
              const double *mc2p_n = gkyl_array_cfetch(mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
              double dxdz[3][3];

              if((ip == nrange->lower[PSI_IDX]) && (up->local.lower[PSI_IDX]== up->global.lower[PSI_IDX]) ) {
                dxdz[0][0] = (-3*mc2p_n[R_IDX] + 4*mc2p_n[6+R_IDX] - mc2p_n[12+R_IDX] )/dzc[0]/2;
                dxdz[1][0] = (-3*mc2p_n[Z_IDX] + 4*mc2p_n[6+Z_IDX] - mc2p_n[12+Z_IDX] )/dzc[0]/2;
              }
              else if((ip == nrange->upper[PSI_IDX]) && (up->local.upper[PSI_IDX]== up->global.upper[PSI_IDX])) {
                dxdz[0][0] = (3*mc2p_n[R_IDX] - 4*mc2p_n[3+R_IDX] + mc2p_n[9+R_IDX] )/dzc[0]/2;
                dxdz[1][0] = (3*mc2p_n[Z_IDX] - 4*mc2p_n[3+Z_IDX] + mc2p_n[9+Z_IDX] )/dzc[0]/2;
              }
              else {
                dxdz[0][0] = -(mc2p_n[3 +R_IDX] -   mc2p_n[6+R_IDX])/2/dzc[0];
                dxdz[1][0] = -(mc2p_n[3 +Z_IDX] -   mc2p_n[6+Z_IDX])/2/dzc[0];
              }


              // dphi/dpsi =0
              dxdz[2][0] = 0.0;
              // Use exact expressions for dR/dtheta and dZ/dtheta, dphi/dtheta
              double *ddtheta_n = gkyl_array_fetch(ddtheta_nodal, gkyl_range_idx(nrange, cidx));
              dxdz[0][2] = ddtheta_n[0];
              dxdz[1][2] = ddtheta_n[1];
              dxdz[2][2] = 0.0;
              double dphidtheta = dxdz[2][2];

              // use exact expressions for d/dalpha
              dxdz[0][1] = 0.0;
              dxdz[1][1] = 0.0;
              dxdz[2][1] = -1.0;

              // dxdz is in cylindrical coords, calculate J as
              // J = R(dR/dpsi*dZ/dtheta - dR/dtheta*dZ/dpsi)
              double R = mc2p_n[R_IDX];
              double jac = sqrt(R*R*(dxdz[0][0]*dxdz[0][0]*dxdz[1][2]*dxdz[1][2] + dxdz[0][2]*dxdz[0][2]*dxdz[1][0]*dxdz[1][0] - 2*dxdz[0][0]*dxdz[0][2]*dxdz[1][0]*dxdz[1][2])) ;

              double *gFld_n= gkyl_array_fetch(gFld_nodal, gkyl_range_idx(nrange, cidx));
              gFld_n[0] = dxdz[0][0]*dxdz[0][0] + R*R*dxdz[2][0]*dxdz[2][0] + dxdz[1][0]*dxdz[1][0]; 
              gFld_n[1] = R*R*dxdz[2][0]; 
              gFld_n[2] = dxdz[0][0]*dxdz[0][2] + R*R*dxdz[2][0]*dphidtheta + dxdz[1][0]*dxdz[1][2];
              gFld_n[3] = R*R; 
              gFld_n[4] = R*R*dphidtheta;
              gFld_n[5] = dxdz[0][2]*dxdz[0][2] + R*R*dphidtheta*dphidtheta + dxdz[1][2]*dxdz[1][2]; 

              double *grFld_n= gkyl_array_fetch(grFld_nodal, gkyl_range_idx(nrange, cidx));
              grFld_n[0] = R*R/jac/jac*(dxdz[1][2]*dxdz[1][2] + dxdz[0][2]*dxdz[0][2] );
              grFld_n[1] = 0.0;
              grFld_n[2] = -R*R/jac/jac*(dxdz[0][0]*dxdz[0][2] + dxdz[1][0]*dxdz[1][2] );
              grFld_n[3] = 1/R/R;
              grFld_n[4] = 0.0;
              grFld_n[5] = R*R/jac/jac*(dxdz[0][0]*dxdz[0][0] + dxdz[1][0]*dxdz[1][0] );

      }
    }
  }
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 6, gFld_nodal, gFld);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 6, grFld_nodal, grFld);
  gkyl_array_release(gFld_nodal);
  gkyl_array_release(grFld_nodal);
}


void gkyl_calc_metric_advance_mirror(
  gkyl_calc_metric *up, struct gkyl_range *nrange,
  struct gkyl_array *mc2p_nodal_fd, struct gkyl_array *ddtheta_nodal,
  struct gkyl_array *bmag_nodal, double *dzc, struct gkyl_array *gFld,
  struct gkyl_array *tanvecFld,
  struct gkyl_array *dualFld,
  struct gkyl_array *dualmagFld,
  struct gkyl_array *normFld,
  struct gkyl_array *jFld, struct gkyl_array* bcartFld, const struct gkyl_range *update_range)
{
  struct gkyl_array* gFld_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nrange->volume);
  struct gkyl_array* jFld_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nrange->volume);
  struct gkyl_array* bcartFld_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange->volume);
  struct gkyl_array* tanvecFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  struct gkyl_array* dualFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  struct gkyl_array* dualmagFld_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange->volume);
  struct gkyl_array* normFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { R_IDX, Z_IDX, PHI_IDX }; // arrangement of cartesian coordinates
  int cidx[3];
  for(int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia){
      for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
          for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
              cidx[PSI_IDX] = ip;
              cidx[AL_IDX] = ia;
              cidx[TH_IDX] = it;
              const double *mc2p_n = gkyl_array_cfetch(mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
              double dxdz[3][3];

              if((ip == nrange->lower[PSI_IDX]) && (up->local.lower[PSI_IDX]== up->global.lower[PSI_IDX]) ) {
                dxdz[0][0] = (-3*mc2p_n[R_IDX] + 4*mc2p_n[6+R_IDX] - mc2p_n[12+R_IDX] )/dzc[0]/2;
                dxdz[1][0] = (-3*mc2p_n[Z_IDX] + 4*mc2p_n[6+Z_IDX] - mc2p_n[12+Z_IDX] )/dzc[0]/2;
                dxdz[2][0] = (-3*mc2p_n[PHI_IDX] + 4*mc2p_n[6+PHI_IDX] - mc2p_n[12+PHI_IDX] )/dzc[0]/2;
              }
              else if((ip == nrange->upper[PSI_IDX]) && (up->local.upper[PSI_IDX]== up->global.upper[PSI_IDX])) {
                dxdz[0][0] = (3*mc2p_n[R_IDX] - 4*mc2p_n[3+R_IDX] + mc2p_n[9+R_IDX] )/dzc[0]/2;
                dxdz[1][0] = (3*mc2p_n[Z_IDX] - 4*mc2p_n[3+Z_IDX] + mc2p_n[9+Z_IDX] )/dzc[0]/2;
                dxdz[2][0] = (3*mc2p_n[PHI_IDX] - 4*mc2p_n[3+PHI_IDX] + mc2p_n[9+PHI_IDX] )/dzc[0]/2;
              }
              else {
                dxdz[0][0] = -(mc2p_n[3 +R_IDX] -   mc2p_n[6+R_IDX])/2/dzc[0];
                dxdz[1][0] = -(mc2p_n[3 +Z_IDX] -   mc2p_n[6+Z_IDX])/2/dzc[0];
                dxdz[2][0] = -(mc2p_n[3 +PHI_IDX] -   mc2p_n[6+PHI_IDX])/2/dzc[0];
              }

              if((ia == nrange->lower[AL_IDX]) && (up->local.lower[AL_IDX]== up->global.lower[AL_IDX]) ) {
                dxdz[0][1] = (-3*mc2p_n[R_IDX] +  4*mc2p_n[18+R_IDX] -  mc2p_n[24+R_IDX])/dzc[1]/2;
                dxdz[1][1] = (-3*mc2p_n[Z_IDX] +  4*mc2p_n[18+Z_IDX] -  mc2p_n[24+Z_IDX])/dzc[1]/2;
                dxdz[2][1] = (-3*mc2p_n[PHI_IDX] +  4*mc2p_n[18+PHI_IDX] -  mc2p_n[24+PHI_IDX])/dzc[1]/2;
              }
              else if((ia == nrange->upper[AL_IDX])  && (up->local.upper[AL_IDX]== up->global.upper[AL_IDX])){
                dxdz[0][1] = (3*mc2p_n[R_IDX] -  4*mc2p_n[15+R_IDX] +  mc2p_n[21+R_IDX] )/dzc[1]/2;
                dxdz[1][1] = (3*mc2p_n[Z_IDX] -  4*mc2p_n[15+Z_IDX] +  mc2p_n[21+Z_IDX] )/dzc[1]/2;
                dxdz[2][1] = (3*mc2p_n[PHI_IDX] -  4*mc2p_n[15+PHI_IDX] +  mc2p_n[21+PHI_IDX] )/dzc[1]/2;
              }
              else {
                dxdz[0][1] = -(mc2p_n[15 +R_IDX] -  mc2p_n[18 +R_IDX])/2/dzc[1];
                dxdz[1][1] = -(mc2p_n[15 +Z_IDX] -  mc2p_n[18 +Z_IDX])/2/dzc[1];
                dxdz[2][1] = -(mc2p_n[15 +PHI_IDX] -  mc2p_n[18 +PHI_IDX])/2/dzc[1];
              }

              if((it == nrange->lower[TH_IDX]) && (up->local.lower[TH_IDX]== up->global.lower[TH_IDX])){
                dxdz[0][2] = (-3*mc2p_n[R_IDX] + 4*mc2p_n[30+R_IDX] - mc2p_n[36+R_IDX])/dzc[2]/2;
                dxdz[1][2] = (-3*mc2p_n[Z_IDX] + 4*mc2p_n[30+Z_IDX] - mc2p_n[36+Z_IDX])/dzc[2]/2;
                dxdz[2][2] = (-3*mc2p_n[PHI_IDX] + 4*mc2p_n[30+PHI_IDX] - mc2p_n[36+PHI_IDX])/dzc[2]/2;
              }
              else if((it == nrange->upper[TH_IDX]) && (up->local.upper[TH_IDX]== up->global.upper[TH_IDX])){
                dxdz[0][2] = (3*mc2p_n[R_IDX] - 4*mc2p_n[27+R_IDX] + mc2p_n[33+R_IDX] )/dzc[2]/2;
                dxdz[1][2] = (3*mc2p_n[Z_IDX] - 4*mc2p_n[27+Z_IDX] + mc2p_n[33+Z_IDX] )/dzc[2]/2;
                dxdz[2][2] = (3*mc2p_n[PHI_IDX] - 4*mc2p_n[27+PHI_IDX] + mc2p_n[33+PHI_IDX] )/dzc[2]/2;
              }
              else {
                dxdz[0][2] = -(mc2p_n[27 +R_IDX] - mc2p_n[30 +R_IDX])/2/dzc[2];
                dxdz[1][2] = -(mc2p_n[27 +Z_IDX] - mc2p_n[30 +Z_IDX])/2/dzc[2];
                dxdz[2][2] = -(mc2p_n[27 +PHI_IDX] - mc2p_n[30 +PHI_IDX])/2/dzc[2];
              }

              // Use exact expressions for dphidtheta, dR/dtheta, and dZ/dtheta
              double *ddtheta_n = gkyl_array_fetch(ddtheta_nodal, gkyl_range_idx(nrange, cidx));
              dxdz[0][2] = ddtheta_n[0]; // dR/dtheta
              dxdz[1][2] = ddtheta_n[1]; // dZ/dtheta
              dxdz[2][2] = ddtheta_n[2]; // dphi/dtheta

              // use exact expressions for d/dalpha
              dxdz[0][1] = 0.0; // dR/dalpha
              dxdz[1][1] = 0.0; // dZ/dalpha
              dxdz[2][1] = -1.0; // dphi/dalpha

              double R = mc2p_n[R_IDX];
              double *bmag_n = gkyl_array_fetch(bmag_nodal, gkyl_range_idx(nrange, cidx));

              double *gFld_n= gkyl_array_fetch(gFld_nodal, gkyl_range_idx(nrange, cidx));
              gFld_n[5] = dxdz[0][2]*dxdz[0][2] + dxdz[1][2]*dxdz[1][2];

              // dxdz is in cylindrical coordinates. Caculate dR/dpsi as
              // dR/dpsi = (1/(dZ/dtheta) ) * [sqrt(g_33)/RB + dR/dtheta*dZ/dtheta]
              // because this is more reliable near R=0 than using finite differences
              double dRdpsi = 1/dxdz[1][2]*(sqrt(gFld_n[5])/bmag_n[0]/R + dxdz[0][2]*dxdz[1][0]);
              // Calculate J as J = R(dR/dpsi*dZ/dtheta - dR/dtheta*dZ/dpsi)
              double *jFld_n= gkyl_array_fetch(jFld_nodal, gkyl_range_idx(nrange, cidx));
              jFld_n[0] = sqrt(R*R*(dRdpsi*dRdpsi*dxdz[1][2]*dxdz[1][2] + dxdz[0][2]*dxdz[0][2]*dxdz[1][0]*dxdz[1][0] - 2*dRdpsi*dxdz[0][2]*dxdz[1][0]*dxdz[1][2])) ;

              gFld_n[0] = dRdpsi*dRdpsi + dxdz[1][0]*dxdz[1][0]; 
              gFld_n[1] = 0.0; 
              gFld_n[2] = dRdpsi*dxdz[0][2] + dxdz[1][0]*dxdz[1][2];
              gFld_n[3] = R*R; 
              gFld_n[4] = 0.0; 

              // Now do bcart
              double *bcartFld_n= gkyl_array_fetch(bcartFld_nodal, gkyl_range_idx(nrange, cidx));
              double phi = mc2p_n[PHI_IDX];
              double b3 = 1/sqrt(gFld_n[5]);
              bcartFld_n[0] = b3*dxdz[0][2]*cos(phi);
              bcartFld_n[1] = b3*dxdz[0][2]*sin(phi);
              bcartFld_n[2] = b3*dxdz[1][2];

              // Set cartesian components of tangents and duals
              double Z = mc2p_n[Z_IDX];
              double J = jFld_n[0];
              double *tanvecFld_n= gkyl_array_fetch(tanvecFld_nodal, gkyl_range_idx(nrange, cidx));
              tanvecFld_n[0] = dxdz[0][0]*cos(phi) - R*sin(phi)*dxdz[2][0]; 
              tanvecFld_n[1] = dxdz[0][0]*sin(phi)  + R*cos(phi)*dxdz[2][0]; 
              tanvecFld_n[2] = dxdz[1][0];

              tanvecFld_n[3] = +R*sin(phi); 
              tanvecFld_n[4] = -R*cos(phi); 
              tanvecFld_n[5] = 0.0; 

              tanvecFld_n[6] = dxdz[0][2]*cos(phi); 
              tanvecFld_n[7] = dxdz[0][2]*sin(phi); 
              tanvecFld_n[8] = dxdz[1][2];

              double *dualFld_n= gkyl_array_fetch(dualFld_nodal, gkyl_range_idx(nrange, cidx));
              dualFld_n[0] = -R/J*cos(phi)*dxdz[1][2];
              dualFld_n[1] = -R/J*sin(phi)*dxdz[1][2];
              dualFld_n[2] = +R/J*dxdz[0][2];

              dualFld_n[3] = 1/J * (dxdz[1][0]*dxdz[0][2]*sin(phi) - dxdz[1][2]*dxdz[0][0]*sin(phi) - dxdz[1][2]*R*cos(phi)*dxdz[2][0] );
              dualFld_n[4] = -1/J * (dxdz[1][0]*dxdz[0][2]*cos(phi) - dxdz[1][2]*dxdz[0][0]*cos(phi) - dxdz[1][2]*R*sin(phi)*dxdz[2][0] );
              dualFld_n[5] =  R/J * dxdz[0][2]*dxdz[2][0];

              dualFld_n[6] = +R/J*cos(phi)*dxdz[1][0];
              dualFld_n[7] = +R/J*sin(phi)*dxdz[1][0];
              dualFld_n[8] = -R/J*dxdz[0][0];

              double norm1 = sqrt(dualFld_n[0]*dualFld_n[0] + dualFld_n[1]*dualFld_n[1] + dualFld_n[2]*dualFld_n[2]);
              double norm2 = sqrt(dualFld_n[3]*dualFld_n[3] + dualFld_n[4]*dualFld_n[4] + dualFld_n[5]*dualFld_n[5]);
              double norm3 = sqrt(dualFld_n[6]*dualFld_n[6] + dualFld_n[7]*dualFld_n[7] + dualFld_n[8]*dualFld_n[8]);

              double *dualmagFld_n = gkyl_array_fetch(dualmagFld_nodal, gkyl_range_idx(nrange, cidx));
              dualmagFld_n[0] = norm1;
              dualmagFld_n[1] = norm2;
              dualmagFld_n[2] = norm3;
              
              // Set normal vectors
              double *normFld_n = gkyl_array_fetch(normFld_nodal, gkyl_range_idx(nrange, cidx));
              normFld_n[0] = dualFld_n[0]/norm1;
              normFld_n[1] = dualFld_n[1]/norm1;
              normFld_n[2] = dualFld_n[2]/norm1;

              normFld_n[3] = dualFld_n[3]/norm2;
              normFld_n[4] = dualFld_n[4]/norm2;
              normFld_n[5] = dualFld_n[5]/norm2;

              normFld_n[6] = dualFld_n[6]/norm3;
              normFld_n[7] = dualFld_n[7]/norm3;
              normFld_n[8] = dualFld_n[8]/norm3;
      }
    }
  }
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 6, gFld_nodal, gFld);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 1, jFld_nodal, jFld);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 3, bcartFld_nodal, bcartFld);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, tanvecFld_nodal, tanvecFld);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, dualFld_nodal, dualFld);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 3, dualmagFld_nodal, dualmagFld);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, normFld_nodal, normFld);
  gkyl_array_release(gFld_nodal);
  gkyl_array_release(jFld_nodal);
  gkyl_array_release(bcartFld_nodal);
  gkyl_array_release(tanvecFld_nodal);
  gkyl_array_release(dualFld_nodal);
  gkyl_array_release(dualmagFld_nodal);
  gkyl_array_release(normFld_nodal);
}


void gkyl_calc_metric_advance(gkyl_calc_metric *up, struct gkyl_range *nrange,
  struct gkyl_array *mc2p_nodal_fd, double *dzc,
  struct gkyl_array *gFld,
  struct gkyl_array *tanvecFld,
  struct gkyl_array *dualFld,
  struct gkyl_array *dualmagFld,
  struct gkyl_array *normFld,
  const struct gkyl_range *update_range)
{
  struct gkyl_array* gFld_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nrange->volume);
  struct gkyl_array* tanvecFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  struct gkyl_array* dualFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  struct gkyl_array* dualmagFld_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange->volume);
  struct gkyl_array* normFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  int cidx[3];
  
  for(int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia) {
      for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
          for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
              cidx[PSI_IDX] = ip;
              cidx[AL_IDX] = ia;
              cidx[TH_IDX] = it;
              const double *mc2p_n = gkyl_array_cfetch(mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
              double dxdz[3][3]; // tan vecs at node
              double dzdx[3][3]; // duals at node

              if((ip == nrange->lower[PSI_IDX]) && (up->local.lower[PSI_IDX]== up->global.lower[PSI_IDX]) ) {
                dxdz[0][0] = (-3*mc2p_n[X_IDX] + 4*mc2p_n[6+X_IDX] - mc2p_n[12+X_IDX] )/dzc[0]/2;
                dxdz[1][0] = (-3*mc2p_n[Y_IDX] + 4*mc2p_n[6+Y_IDX] - mc2p_n[12+Y_IDX] )/dzc[0]/2;
                dxdz[2][0] = (-3*mc2p_n[Z_IDX] + 4*mc2p_n[6+Z_IDX] - mc2p_n[12+Z_IDX] )/dzc[0]/2;
              }
              else if((ip == nrange->upper[PSI_IDX]) && (up->local.upper[PSI_IDX]== up->global.upper[PSI_IDX])) {
                dxdz[0][0] = (3*mc2p_n[X_IDX] - 4*mc2p_n[3+X_IDX] + mc2p_n[9+X_IDX] )/dzc[0]/2;
                dxdz[1][0] = (3*mc2p_n[Y_IDX] - 4*mc2p_n[3+Y_IDX] + mc2p_n[9+Y_IDX] )/dzc[0]/2;
                dxdz[2][0] = (3*mc2p_n[Z_IDX] - 4*mc2p_n[3+Z_IDX] + mc2p_n[9+Z_IDX] )/dzc[0]/2;
              }
              else{
                dxdz[0][0] = -(mc2p_n[3 +X_IDX] -   mc2p_n[6+X_IDX])/2/dzc[0];
                dxdz[1][0] = -(mc2p_n[3 +Y_IDX] -   mc2p_n[6+Y_IDX])/2/dzc[0];
                dxdz[2][0] = -(mc2p_n[3 +Z_IDX] -   mc2p_n[6+Z_IDX])/2/dzc[0];
              }


              if((ia == nrange->lower[AL_IDX]) && (up->local.lower[AL_IDX]== up->global.lower[AL_IDX]) ) {
                dxdz[0][1] = (-3*mc2p_n[X_IDX] +  4*mc2p_n[18+X_IDX] -  mc2p_n[24+X_IDX])/dzc[1]/2;
                dxdz[1][1] = (-3*mc2p_n[Y_IDX] +  4*mc2p_n[18+Y_IDX] -  mc2p_n[24+Y_IDX])/dzc[1]/2;
                dxdz[2][1] = (-3*mc2p_n[Z_IDX] +  4*mc2p_n[18+Z_IDX] -  mc2p_n[24+Z_IDX])/dzc[1]/2;
              }
              else if((ia == nrange->upper[AL_IDX])  && (up->local.upper[AL_IDX]== up->global.upper[AL_IDX])) {
                dxdz[0][1] = (3*mc2p_n[X_IDX] -  4*mc2p_n[15+X_IDX] +  mc2p_n[21+X_IDX] )/dzc[1]/2;
                dxdz[1][1] = (3*mc2p_n[Y_IDX] -  4*mc2p_n[15+Y_IDX] +  mc2p_n[21+Y_IDX] )/dzc[1]/2;
                dxdz[2][1] = (3*mc2p_n[Z_IDX] -  4*mc2p_n[15+Z_IDX] +  mc2p_n[21+Z_IDX] )/dzc[1]/2;
              }
              else {
                dxdz[0][1] = -(mc2p_n[15 +X_IDX] -  mc2p_n[18 +X_IDX])/2/dzc[1];
                dxdz[1][1] = -(mc2p_n[15 +Y_IDX] -  mc2p_n[18 +Y_IDX])/2/dzc[1];
                dxdz[2][1] = -(mc2p_n[15 +Z_IDX] -  mc2p_n[18 +Z_IDX])/2/dzc[1];
              }

              if((it == nrange->lower[TH_IDX]) && (up->local.lower[TH_IDX]== up->global.lower[TH_IDX])) {
                dxdz[0][2] = (-3*mc2p_n[X_IDX] + 4*mc2p_n[30+X_IDX] - mc2p_n[36+X_IDX])/dzc[2]/2;
                dxdz[1][2] = (-3*mc2p_n[Y_IDX] + 4*mc2p_n[30+Y_IDX] - mc2p_n[36+Y_IDX])/dzc[2]/2;
                dxdz[2][2] = (-3*mc2p_n[Z_IDX] + 4*mc2p_n[30+Z_IDX] - mc2p_n[36+Z_IDX])/dzc[2]/2;
              }
              else if((it == nrange->upper[TH_IDX]) && (up->local.upper[TH_IDX]== up->global.upper[TH_IDX])) {
                dxdz[0][2] = (3*mc2p_n[X_IDX] - 4*mc2p_n[27+X_IDX] + mc2p_n[33+X_IDX] )/dzc[2]/2;
                dxdz[1][2] = (3*mc2p_n[Y_IDX] - 4*mc2p_n[27+Y_IDX] + mc2p_n[33+Y_IDX] )/dzc[2]/2;
                dxdz[2][2] = (3*mc2p_n[Z_IDX] - 4*mc2p_n[27+Z_IDX] + mc2p_n[33+Z_IDX] )/dzc[2]/2;
              }
              else{
                dxdz[0][2] = -(mc2p_n[27 +X_IDX] - mc2p_n[30 +X_IDX])/2/dzc[2];
                dxdz[1][2] = -(mc2p_n[27 +Y_IDX] - mc2p_n[30 +Y_IDX])/2/dzc[2];
                dxdz[2][2] = -(mc2p_n[27 +Z_IDX] - mc2p_n[30 +Z_IDX])/2/dzc[2];
              }

              double *gFld_n= gkyl_array_fetch(gFld_nodal, gkyl_range_idx(nrange, cidx));
              gFld_n[0] = calc_metric(dxdz, 1, 1); 
              gFld_n[1] = calc_metric(dxdz, 1, 2); 
              gFld_n[2] = calc_metric(dxdz, 1, 3); 
              gFld_n[3] = calc_metric(dxdz, 2, 2); 
              gFld_n[4] = calc_metric(dxdz, 2, 3); 
              gFld_n[5] = calc_metric(dxdz, 3, 3); 

              double Jsq = gFld_n[0]*(gFld_n[3]*gFld_n[5] - gFld_n[4]*gFld_n[4] ) - gFld_n[1]*(gFld_n[1]*gFld_n[5] - gFld_n[4]*gFld_n[2] ) + gFld_n[2]*(gFld_n[1]*gFld_n[4] - gFld_n[3]*gFld_n[2] );
              double J = sqrt(Jsq);
              double e_1[3], e_2[3], e_3[3];
              e_1[0] = dxdz[0][0]; e_1[1] = dxdz[1][0]; e_1[2] = dxdz[2][0];
              e_2[0] = dxdz[0][1]; e_2[1] = dxdz[1][1]; e_2[2] = dxdz[2][1];
              e_3[0] = dxdz[0][2]; e_3[1] = dxdz[1][2]; e_3[2] = dxdz[2][2];
              calc_dual(J, e_2, e_3, dzdx[0]);
              calc_dual(J, e_3, e_1, dzdx[1]);
              calc_dual(J, e_1, e_2, dzdx[2]);

              double *dualFld_n= gkyl_array_fetch(dualFld_nodal, gkyl_range_idx(nrange, cidx));
              dualFld_n[0] = dzdx[0][0];
              dualFld_n[1] = dzdx[0][1];
              dualFld_n[2] = dzdx[0][2];
              dualFld_n[3] = dzdx[1][0];
              dualFld_n[4] = dzdx[1][1];
              dualFld_n[5] = dzdx[1][2];
              dualFld_n[6] = dzdx[2][0];
              dualFld_n[7] = dzdx[2][1];
              dualFld_n[8] = dzdx[2][2];

              double *tanvecFld_n= gkyl_array_fetch(tanvecFld_nodal, gkyl_range_idx(nrange, cidx));
              tanvecFld_n[0] = dxdz[0][0]; 
              tanvecFld_n[1] = dxdz[1][0]; 
              tanvecFld_n[2] = dxdz[2][0]; 
              tanvecFld_n[3] = dxdz[0][1]; 
              tanvecFld_n[4] = dxdz[1][1]; 
              tanvecFld_n[5] = dxdz[2][1]; 
              tanvecFld_n[6] = dxdz[0][2]; 
              tanvecFld_n[7] = dxdz[1][2]; 
              tanvecFld_n[8] = dxdz[2][2]; 

              double norm1 = sqrt(dualFld_n[0]*dualFld_n[0] + dualFld_n[1]*dualFld_n[1] + dualFld_n[2]*dualFld_n[2]);
              double norm2 = sqrt(dualFld_n[3]*dualFld_n[3] + dualFld_n[4]*dualFld_n[4] + dualFld_n[5]*dualFld_n[5]);
              double norm3 = sqrt(dualFld_n[6]*dualFld_n[6] + dualFld_n[7]*dualFld_n[7] + dualFld_n[8]*dualFld_n[8]);

              double *dualmagFld_n = gkyl_array_fetch(dualmagFld_nodal, gkyl_range_idx(nrange, cidx));
              dualmagFld_n[0] = norm1;
              dualmagFld_n[1] = norm2;
              dualmagFld_n[2] = norm3;
              
              // Set normal vectors
              double *normFld_n = gkyl_array_fetch(normFld_nodal, gkyl_range_idx(nrange, cidx));
              normFld_n[0] = dualFld_n[0]/norm1;
              normFld_n[1] = dualFld_n[1]/norm1;
              normFld_n[2] = dualFld_n[2]/norm1;

              normFld_n[3] = dualFld_n[3]/norm2;
              normFld_n[4] = dualFld_n[4]/norm2;
              normFld_n[5] = dualFld_n[5]/norm2;

              normFld_n[6] = dualFld_n[6]/norm3;
              normFld_n[7] = dualFld_n[7]/norm3;
              normFld_n[8] = dualFld_n[8]/norm3;
      }
    }
  }
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 6, gFld_nodal, gFld);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, tanvecFld_nodal, tanvecFld);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, dualFld_nodal, dualFld);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 3, dualmagFld_nodal, dualmagFld);
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, normFld_nodal, normFld);
  gkyl_array_release(gFld_nodal);
  gkyl_array_release(tanvecFld_nodal);
  gkyl_array_release(dualFld_nodal);
  gkyl_array_release(dualmagFld_nodal);
  gkyl_array_release(normFld_nodal);
}

void gkyl_calc_metric_advance_bcart(gkyl_calc_metric *up, struct gkyl_range *nrange,
  struct gkyl_array *biFld, struct gkyl_array *dualFld, struct gkyl_array *bcartFld,
  const struct gkyl_range *update_range)
{


  struct gkyl_array* bcartFld_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange->volume);
  struct gkyl_array* biFld_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange->volume);
  struct gkyl_array* dualFld_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  int cidx[3];


  // Fill the inputs
  gkyl_nodal_ops_m2n(up->n2m, up->cbasis, up->grid, nrange, update_range, 9, dualFld_nodal, dualFld);
  gkyl_nodal_ops_m2n(up->n2m, up->cbasis, up->grid, nrange, update_range, 3, biFld_nodal, biFld);
  
  for(int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia) {
      for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
          for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
              cidx[PSI_IDX] = ip;
              cidx[AL_IDX] = ia;
              cidx[TH_IDX] = it;

              double *biFld_n= gkyl_array_fetch(biFld_nodal, gkyl_range_idx(nrange, cidx));
              double *dualFld_n= gkyl_array_fetch(dualFld_nodal, gkyl_range_idx(nrange, cidx));
              double dzdx[3][3]; // duals at node
              dzdx[0][0] = dualFld_n[0];
              dzdx[0][1] = dualFld_n[1];
              dzdx[0][2] = dualFld_n[2];
              dzdx[1][0] = dualFld_n[3];
              dzdx[1][1] = dualFld_n[4];
              dzdx[1][2] = dualFld_n[5];
              dzdx[2][0] = dualFld_n[6];
              dzdx[2][1] = dualFld_n[7];
              dzdx[2][2] = dualFld_n[8];
              double *bcartFld_n= gkyl_array_fetch(bcartFld_nodal, gkyl_range_idx(nrange, cidx));
              bcartFld_n[0] = dzdx[0][0]*biFld_n[0] + dzdx[1][0]*biFld_n[1] + dzdx[2][0]*biFld_n[2];
              bcartFld_n[1] = dzdx[0][1]*biFld_n[0] + dzdx[1][1]*biFld_n[1] + dzdx[2][1]*biFld_n[2];
              bcartFld_n[2] = dzdx[0][2]*biFld_n[0] + dzdx[1][2]*biFld_n[1] + dzdx[2][2]*biFld_n[2];
          }
      }
  }
  gkyl_nodal_ops_n2m(up->n2m, up->cbasis, up->grid, nrange, update_range, 3, bcartFld_nodal, bcartFld);
  gkyl_array_release(bcartFld_nodal);
  gkyl_array_release(biFld_nodal);
  gkyl_array_release(dualFld_nodal);
}

void
gkyl_calc_metric_release(gkyl_calc_metric* up)
{
  gkyl_nodal_ops_release(up->n2m);
  gkyl_free(up);
}
