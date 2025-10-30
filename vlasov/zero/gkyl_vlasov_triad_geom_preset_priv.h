#include <gkyl_vlasov_triad_geom.h>

void
eval_cylindrical_vierbein_1v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Downstairs components of the Vierbeins
  double q_r = xn[0];

  // e_ij = e_i(x) . \sigma_j(x)
  double e_rr = 1.0; // Vierbein Coefficients (r-r coefficient).
  
  fout[0] = e_rr;
}

void
eval_cylindrical_vierbein_gradient_1v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Downstairs components of the Vierbein Gradients 
  double q_r = xn[0];

  // d (e_ij) / dr 
  double d_e_rr_dr = 0.0; // Vierbein Gradient Coefficients (r-r coefficient).
  
  fout[0] = d_e_rr_dr;
}

void
eval_cylindrical_vierbein_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Downstairs components of the Vierbeins
  double q_r = xn[0];

  // e_ij = e_i(x) . \sigma_j(x)
  double e_rr = 1.0; // Vierbein Coefficients (r-r coefficient).
  double e_rt = 0.0; // Vierbein Coefficients (r-theta coefficient).
  double e_tr = 0.0; // Vierbein Coefficients (theta-r coefficient).
  double e_tt = q_r; // Vierbein Coefficients (theta-theta coefficient).
  
  fout[0] = e_rr;
  fout[1] = e_rt;
  fout[2] = e_tr;
  fout[3] = e_tt;
}

void
eval_cylindrical_vierbein_gradient_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Downstairs components of the Vierbein Gradients 
  double q_r = xn[0];

  // d (e_ij) / dr 
  double d_e_rr_dr = 0.0; // Vierbein Gradient Coefficients (r-r coefficient).
  double d_e_rt_dr = 0.0; // Vierbein Gradient Coefficients (r-theta coefficient).
  double d_e_tr_dr = 0.0; // Vierbein Gradient Coefficients (theta-r coefficient).
  double d_e_tt_dr = 1.0; // Vierbein Gradient Coefficients (theta-theta coefficient).

  // d (e_ij) / dtheta 
  double d_e_rr_dt = 0.0; // Vierbein Gradient Coefficients (r-r coefficient).
  double d_e_rt_dt = 0.0; // Vierbein Gradient Coefficients (r-theta coefficient).
  double d_e_tr_dt = 0.0; // Vierbein Gradient Coefficients (theta-r coefficient).
  double d_e_tt_dt = 0.0; // Vierbein Gradient Coefficients (theta-theta coefficient).
  
  fout[0] = d_e_rr_dr;
  fout[1] = d_e_rt_dr;
  fout[2] = d_e_tr_dr;
  fout[3] = d_e_tt_dr;

  fout[4] = d_e_rr_dt;
  fout[5] = d_e_rt_dt;
  fout[6] = d_e_tr_dt;
  fout[7] = d_e_tt_dt;

}

void
eval_cylindrical_vierbein_3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Downstairs components of the Vierbeins
  double q_r = xn[0];

  // e_ij = e_i(x) . \sigma_j(x)
  double e_rr = 1.0; // Vierbein Coefficients (r-r coefficient).
  double e_rz = 0.0;
  double e_rt = 0.0; // Vierbein Coefficients (r-theta coefficient).
  
  double e_zr = 0.0;
  double e_zz = 1.0;
  double e_zt = 0.0; 

  double e_tr = 0.0; // Vierbein Coefficients (theta-r coefficient).
  double e_tz = 0.0; 
  double e_tt = q_r; // Vierbein Coefficients (theta-theta coefficient).
  
  fout[0] = e_rr;
  fout[1] = e_rz;
  fout[2] = e_rt;

  fout[3] = e_zr;
  fout[4] = e_zz;
  fout[5] = e_zt;

  fout[6] = e_tr;
  fout[7] = e_tz;
  fout[8] = e_tt;
}

void
eval_cylindrical_vierbein_gradient_3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Downstairs components of the Vierbein Gradients 
  double q_r = xn[0];

  // d (e_ij) / dr 
  double d_e_rr_dr = 0.0; // Vierbein Gradient Coefficients (r-r coefficient).
  double d_e_rz_dr = 0.0;
  double d_e_rt_dr = 0.0; // Vierbein Gradient Coefficients (r-theta coefficient).
  
  double d_e_zr_dr = 0.0;
  double d_e_zz_dr = 0.0;
  double d_e_zt_dr = 0.0;

  double d_e_tr_dr = 0.0; // Vierbein Gradient Coefficients (theta-r coefficient).
  double d_e_tz_dr = 0.0;
  double d_e_tt_dr = 1.0; // Vierbein Gradient Coefficients (theta-theta coefficient).

  // d(e_ij) / dz
  double d_e_rr_dz = 0.0; // Vierbein Gradient Coefficients (r-r coefficient).
  double d_e_rz_dz = 0.0;
  double d_e_rt_dz = 0.0; // Vierbein Gradient Coefficients (r-theta coefficient).
  
  double d_e_zr_dz = 0.0;
  double d_e_zz_dz = 0.0;
  double d_e_zt_dz = 0.0;

  double d_e_tr_dz = 0.0; // Vierbein Gradient Coefficients (theta-r coefficient).
  double d_e_tz_dz = 0.0;
  double d_e_tt_dz = 0.0; // Vierbein Gradient Coefficients (theta-theta coefficient).


  // d (e_ij) / dtheta 
  double d_e_rr_dt = 0.0; // Vierbein Gradient Coefficients (r-r coefficient).
  double d_e_rz_dt = 0.0;
  double d_e_rt_dt = 0.0; // Vierbein Gradient Coefficients (r-theta coefficient).

  double d_e_zr_dt = 0.0;
  double d_e_zz_dt = 0.0;
  double d_e_zt_dt = 0.0;

  double d_e_tr_dt = 0.0; // Vierbein Gradient Coefficients (theta-r coefficient).
  double d_e_tz_dt = 0.0;
  double d_e_tt_dt = 0.0; // Vierbein Gradient Coefficients (theta-theta coefficient).
  
  fout[0] = d_e_rr_dr;
  fout[1] = d_e_rz_dr;
  fout[2] = d_e_rt_dr;
  fout[3] = d_e_zr_dr;
  fout[4] = d_e_zz_dr;
  fout[5] = d_e_zt_dr;
  fout[6] = d_e_tr_dr;
  fout[7] = d_e_tz_dr;
  fout[8] = d_e_tt_dr;

  fout[9+0] = d_e_rr_dz;
  fout[9+1] = d_e_rz_dz;
  fout[9+2] = d_e_rt_dz;
  fout[9+3] = d_e_zr_dz;
  fout[9+4] = d_e_zz_dz;
  fout[9+5] = d_e_zt_dz;
  fout[9+6] = d_e_tr_dz;
  fout[9+7] = d_e_tz_dz;
  fout[9+8] = d_e_tt_dz;

  fout[18+0] = d_e_rr_dt;
  fout[18+1] = d_e_rz_dt;
  fout[18+2] = d_e_rt_dt;
  fout[18+3] = d_e_zr_dt;
  fout[18+4] = d_e_zz_dt;
  fout[18+5] = d_e_zt_dt;
  fout[18+6] = d_e_tr_dt;
  fout[18+7] = d_e_tz_dt;
  fout[18+8] = d_e_tt_dt;

}

// for use in kernel tables
typedef struct { evalf_t kernels[3]; } vierbein_kern_list;
typedef struct { evalf_t kernels[3]; } vierbein_gradient_kern_list;

static const vierbein_kern_list cylindrical_vierbein_list[] = {
  { eval_cylindrical_vierbein_1v },
  { eval_cylindrical_vierbein_2v },
  { eval_cylindrical_vierbein_3v }
};

static const vierbein_gradient_kern_list cylindrical_vierbein_gradient_list[] = {
  { eval_cylindrical_vierbein_gradient_1v },
  { eval_cylindrical_vierbein_gradient_2v },
  { eval_cylindrical_vierbein_gradient_3v }
};

static evalf_t
choose_vierbein_kern(enum gkyl_triad_vierbein_type type, int vdim)
{
  switch(type) {
    case GKYL_TRIAD_CYLINDRICAL:
      return cylindrical_vierbein_list[vdim-1].kernels[0];
      break;
    default:
      assert(false);
  }
}

static evalf_t
choose_vierbein_gradient_kern(enum gkyl_triad_vierbein_type type, int vdim)
{
  switch(type) {
    case GKYL_TRIAD_CYLINDRICAL:
      return cylindrical_vierbein_gradient_list[vdim-1].kernels[0];
      break;
    default:
      assert(false);
  }
}
