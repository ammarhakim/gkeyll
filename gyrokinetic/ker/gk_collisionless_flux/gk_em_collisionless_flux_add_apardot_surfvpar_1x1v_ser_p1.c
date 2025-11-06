#include <gkyl_gk_collisionless_flux_kernels.h> 
GKYL_CU_DH double gk_em_collisionless_flux_add_apardot_surfvpar_1x1v_ser_p1(
    const double *w, const double *dxv, 
    const double *vmap_prime_l, const double *vmap_prime_r,
    const double *vmap, const double *vmapSq, const double q_, const double m_,
    const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv,
    const double *bmag, const double *phi, const double *apar, const double *JfL, const double *JfR,
    double* GKYL_RESTRICT flux_surf) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap_prime_l,vmap_prime_r: velocity space mapping derivative in left and right cells.
  // vmap: velocity space mapping.
  // vmapSq: velocity space mapping squared.
  // q_,m_: species charge and mass.
  // dgv: volume DG geometry.
  // gkdgv: gyrokinetic volume DG geometry.
  // bmag: magnetic field amplitude.
  // phi: electrostatic potential.
  // apar: parallel component of vector potential (here it will be apardot).
  // JfL: distribution times total jacobian in left cell.
  // JfR: distribution times total jacobian in right cell.
  // flux_surf: output surface phase space flux in each direction (cdim + 1 components).
  //            Note: Each cell owns their *lower* edge surface evaluation.

  double rdx2 = 2.0/dxv[0];
  double rdvpar2 = 2.0/dxv[1];


  double *flux_surf_nodal = &flux_surf[3]; 
  double cfl = 0.0; 
  double alpha_quad = 0.0; 
  double JfL_quad = 0.0; 
  double JfR_quad = 0.0; 
  double Jfavg_quad = 0.0; 
  double Jfjump_quad = 0.0; 



  alpha_quad = -q_/m_*(0.7071067811865475*apar[0]-0.7071067811865475*apar[1]); 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (-(1.1180339887498951*JfL[5])+1.118033988749895*JfL[4]-0.8660254037844386*JfL[3]+0.8660254037844386*JfL[2]-0.5*JfL[1]+0.5*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (-(1.1180339887498951*JfR[5])+1.118033988749895*JfR[4]+0.8660254037844386*JfR[3]-0.8660254037844386*JfR[2]-0.5*JfR[1]+0.5*JfR[0])/vmap_prime_r[0];
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[0] += alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;



  alpha_quad = -q_/m_*(0.7071067811865475*apar[1]+0.7071067811865475*apar[0]); 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (1.1180339887498951*JfL[5]+1.118033988749895*JfL[4]+0.8660254037844386*JfL[3]+0.8660254037844386*JfL[2]+0.5*JfL[1]+0.5*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (1.1180339887498951*JfR[5]+1.118033988749895*JfR[4]-0.8660254037844386*JfR[3]-0.8660254037844386*JfR[2]+0.5*JfR[1]+0.5*JfR[0])/vmap_prime_r[0];
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[1] += alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;


  double vmap_prime_min = fmin(fabs(vmap_prime_l[0]),fabs(vmap_prime_r[0]));

  return cfl/vmap_prime_min*2.5*rdvpar2; 

} 
