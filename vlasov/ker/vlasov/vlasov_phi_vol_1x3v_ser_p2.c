#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_phi_vol_1x3v_ser_p2(const double *w, const double *dxv,
  const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  const double *jacob_cx = &jacob_pos[0]; 
  const double jacob_cx_inv = 1.0/jacob_cx[0]; 
  out[2] += dv10*dx10*(-(4.743416490252569*f[1]*phi[2]*jacob_cx_inv*jacob_vx_inv)-2.1213203435596424*f[0]*phi[1]*jacob_cx_inv*jacob_vx_inv); 
  out[5] += dv10*dx10*(phi[2]*(-(4.242640687119286*f[11]*jacob_cx_inv*jacob_vx_inv)-4.743416490252569*f[0]*jacob_cx_inv*jacob_vx_inv)-2.1213203435596424*f[1]*phi[1]*jacob_cx_inv*jacob_vx_inv); 
  out[7] += dv10*dx10*(-(4.743416490252569*phi[2]*f[6]*jacob_cx_inv*jacob_vx_inv)-2.1213203435596424*phi[1]*f[3]*jacob_cx_inv*jacob_vx_inv); 
  out[9] += dv10*dx10*(-(4.743416490252569*phi[2]*f[8]*jacob_cx_inv*jacob_vx_inv)-2.1213203435596424*phi[1]*f[4]*jacob_cx_inv*jacob_vx_inv); 
  out[12] += dv10*dx10*(-(10.606601717798211*phi[2]*f[5]*jacob_cx_inv*jacob_vx_inv)-4.743416490252569*phi[1]*f[2]*jacob_cx_inv*jacob_vx_inv); 
  out[15] += dv10*dx10*(phi[2]*(-(4.242640687119286*f[21]*jacob_cx_inv*jacob_vx_inv)-4.743416490252569*f[3]*jacob_cx_inv*jacob_vx_inv)-2.1213203435596424*phi[1]*f[6]*jacob_cx_inv*jacob_vx_inv); 
  out[16] += dv10*dx10*(phi[2]*(-(4.242640687119286*f[25]*jacob_cx_inv*jacob_vx_inv)-4.743416490252569*f[4]*jacob_cx_inv*jacob_vx_inv)-2.1213203435596424*phi[1]*f[8]*jacob_cx_inv*jacob_vx_inv); 
  out[18] += dv10*dx10*(-(4.743416490252569*phi[2]*f[17]*jacob_cx_inv*jacob_vx_inv)-2.1213203435596424*phi[1]*f[10]*jacob_cx_inv*jacob_vx_inv); 
  out[19] += dv10*dx10*(-(2.1213203435596424*phi[1]*f[11]*jacob_cx_inv*jacob_vx_inv)-4.242640687119286*f[1]*phi[2]*jacob_cx_inv*jacob_vx_inv); 
  out[20] += dv10*dx10*(phi[2]*(-(9.48683298050514*f[19]*jacob_cx_inv*jacob_vx_inv)-10.606601717798213*f[2]*jacob_cx_inv*jacob_vx_inv)-4.743416490252569*phi[1]*f[5]*jacob_cx_inv*jacob_vx_inv); 
  out[22] += dv10*dx10*(-(10.606601717798213*phi[2]*f[15]*jacob_cx_inv*jacob_vx_inv)-4.743416490252569*phi[1]*f[7]*jacob_cx_inv*jacob_vx_inv); 
  out[24] += dv10*dx10*(-(4.743416490252569*phi[2]*f[23]*jacob_cx_inv*jacob_vx_inv)-2.1213203435596424*phi[1]*f[13]*jacob_cx_inv*jacob_vx_inv); 
  out[26] += dv10*dx10*(-(10.606601717798213*phi[2]*f[16]*jacob_cx_inv*jacob_vx_inv)-4.743416490252569*phi[1]*f[9]*jacob_cx_inv*jacob_vx_inv); 
  out[29] += dv10*dx10*(-(4.743416490252569*phi[2]*f[28]*jacob_cx_inv*jacob_vx_inv)-2.1213203435596424*phi[1]*f[14]*jacob_cx_inv*jacob_vx_inv); 
  out[31] += dv10*dx10*(phi[2]*(-(4.242640687119286*f[37]*jacob_cx_inv*jacob_vx_inv)-4.743416490252569*f[10]*jacob_cx_inv*jacob_vx_inv)-2.1213203435596424*phi[1]*f[17]*jacob_cx_inv*jacob_vx_inv); 
  out[32] += dv10*dx10*(-(2.1213203435596424*phi[1]*f[21]*jacob_cx_inv*jacob_vx_inv)-4.242640687119286*phi[2]*f[6]*jacob_cx_inv*jacob_vx_inv); 
  out[33] += dv10*dx10*(phi[2]*(-(9.48683298050514*f[32]*jacob_cx_inv*jacob_vx_inv)-10.606601717798211*f[7]*jacob_cx_inv*jacob_vx_inv)-4.743416490252569*phi[1]*f[15]*jacob_cx_inv*jacob_vx_inv); 
  out[34] += dv10*dx10*(-(2.1213203435596424*phi[1]*f[23]*jacob_cx_inv*jacob_vx_inv)-4.743416490252569*phi[2]*f[13]*jacob_cx_inv*jacob_vx_inv); 
  out[35] += dv10*dx10*(-(2.1213203435596424*phi[1]*f[25]*jacob_cx_inv*jacob_vx_inv)-4.242640687119286*phi[2]*f[8]*jacob_cx_inv*jacob_vx_inv); 
  out[36] += dv10*dx10*(phi[2]*(-(9.48683298050514*f[35]*jacob_cx_inv*jacob_vx_inv)-10.606601717798211*f[9]*jacob_cx_inv*jacob_vx_inv)-4.743416490252569*phi[1]*f[16]*jacob_cx_inv*jacob_vx_inv); 
  out[38] += dv10*dx10*(-(10.606601717798211*phi[2]*f[31]*jacob_cx_inv*jacob_vx_inv)-4.743416490252569*phi[1]*f[18]*jacob_cx_inv*jacob_vx_inv); 
  out[40] += dv10*dx10*(-(4.743416490252569*phi[2]*f[39]*jacob_cx_inv*jacob_vx_inv)-2.1213203435596424*phi[1]*f[27]*jacob_cx_inv*jacob_vx_inv); 
  out[41] += dv10*dx10*(-(2.1213203435596424*phi[1]*f[28]*jacob_cx_inv*jacob_vx_inv)-4.743416490252569*phi[2]*f[14]*jacob_cx_inv*jacob_vx_inv); 
  out[43] += dv10*dx10*(-(4.743416490252569*phi[2]*f[42]*jacob_cx_inv*jacob_vx_inv)-2.1213203435596424*phi[1]*f[30]*jacob_cx_inv*jacob_vx_inv); 
  out[44] += dv10*dx10*(-(2.1213203435596424*phi[1]*f[37]*jacob_cx_inv*jacob_vx_inv)-4.242640687119286*phi[2]*f[17]*jacob_cx_inv*jacob_vx_inv); 
  out[45] += dv10*dx10*(phi[2]*(-(9.48683298050514*f[44]*jacob_cx_inv*jacob_vx_inv)-10.606601717798213*f[18]*jacob_cx_inv*jacob_vx_inv)-4.743416490252569*phi[1]*f[31]*jacob_cx_inv*jacob_vx_inv); 
  out[46] += dv10*dx10*(-(2.1213203435596424*phi[1]*f[39]*jacob_cx_inv*jacob_vx_inv)-4.743416490252569*phi[2]*f[27]*jacob_cx_inv*jacob_vx_inv); 
  out[47] += dv10*dx10*(-(2.1213203435596424*phi[1]*f[42]*jacob_cx_inv*jacob_vx_inv)-4.743416490252569*phi[2]*f[30]*jacob_cx_inv*jacob_vx_inv); 

} 
