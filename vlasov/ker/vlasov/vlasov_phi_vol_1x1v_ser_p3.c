#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_phi_vol_1x1v_ser_p3(const double *w, const double *dxv,
  const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  const double *jacob_cx = &jacob_pos[0]; 
  const double jacob_cx_inv = 1.0/jacob_cx[0]; 
  out[2] += dv10*dx10*(phi[3]*(-(7.24568837309472*f[4]*jacob_cx_inv*jacob_vx_inv)-3.24037034920393*f[0]*jacob_cx_inv*jacob_vx_inv)-4.743416490252569*f[1]*phi[2]*jacob_cx_inv*jacob_vx_inv-2.1213203435596424*f[0]*phi[1]*jacob_cx_inv*jacob_vx_inv); 
  out[3] += dv10*dx10*(phi[3]*(-(6.363961030678928*f[8]*jacob_cx_inv*jacob_vx_inv)-9.721111047611789*f[1]*jacob_cx_inv*jacob_vx_inv)+phi[2]*(-(4.242640687119286*f[4]*jacob_cx_inv*jacob_vx_inv)-4.743416490252569*f[0]*jacob_cx_inv*jacob_vx_inv)-2.1213203435596424*f[1]*phi[1]*jacob_cx_inv*jacob_vx_inv); 
  out[5] += dv10*dx10*(phi[3]*(-(16.201851746019653*f[6]*jacob_cx_inv*jacob_vx_inv)-7.24568837309472*f[2]*jacob_cx_inv*jacob_vx_inv)-10.606601717798211*phi[2]*f[3]*jacob_cx_inv*jacob_vx_inv-4.743416490252569*phi[1]*f[2]*jacob_cx_inv*jacob_vx_inv); 
  out[6] += dv10*dx10*(phi[2]*(-(4.16619044897648*f[8]*jacob_cx_inv*jacob_vx_inv)-4.242640687119286*f[1]*jacob_cx_inv*jacob_vx_inv)+phi[3]*(-(7.869470848066687*f[4]*jacob_cx_inv*jacob_vx_inv)-7.245688373094719*f[0]*jacob_cx_inv*jacob_vx_inv)-2.1213203435596424*phi[1]*f[4]*jacob_cx_inv*jacob_vx_inv); 
  out[7] += dv10*dx10*(phi[3]*(-(14.230249470757704*f[10]*jacob_cx_inv*jacob_vx_inv)-21.73706511928416*f[3]*jacob_cx_inv*jacob_vx_inv)+phi[2]*(-(9.48683298050514*f[6]*jacob_cx_inv*jacob_vx_inv)-10.606601717798213*f[2]*jacob_cx_inv*jacob_vx_inv)-4.743416490252569*phi[1]*f[3]*jacob_cx_inv*jacob_vx_inv); 
  out[9] += dv10*dx10*(phi[2]*(-(16.201851746019653*f[7]*jacob_cx_inv*jacob_vx_inv)-7.24568837309472*f[1]*jacob_cx_inv*jacob_vx_inv)+phi[1]*(-(7.24568837309472*f[5]*jacob_cx_inv*jacob_vx_inv)-3.24037034920393*f[0]*jacob_cx_inv*jacob_vx_inv)+phi[3]*(-(11.067971810589325*f[5]*jacob_cx_inv*jacob_vx_inv)-11.067971810589325*f[4]*jacob_cx_inv*jacob_vx_inv-4.949747468305832*f[0]*jacob_cx_inv*jacob_vx_inv)); 
  out[10] += dv10*dx10*(phi[3]*(-(7.560864148142503*f[8]*jacob_cx_inv*jacob_vx_inv)-6.363961030678925*f[1]*jacob_cx_inv*jacob_vx_inv)-2.121320343559642*phi[1]*f[8]*jacob_cx_inv*jacob_vx_inv-4.166190448976481*phi[2]*f[4]*jacob_cx_inv*jacob_vx_inv); 
  out[11] += dv10*dx10*(phi[3]*(-(9.721111047611789*f[8]*jacob_cx_inv*jacob_vx_inv)-33.20391543176798*f[7]*jacob_cx_inv*jacob_vx_inv-14.849242404917497*f[1]*jacob_cx_inv*jacob_vx_inv)+phi[1]*(-(7.245688373094719*f[7]*jacob_cx_inv*jacob_vx_inv)-3.2403703492039297*f[1]*jacob_cx_inv*jacob_vx_inv)+phi[2]*(-(16.20185174601965*f[5]*jacob_cx_inv*jacob_vx_inv)-6.48074069840786*f[4]*jacob_cx_inv*jacob_vx_inv-7.245688373094718*f[0]*jacob_cx_inv*jacob_vx_inv)); 

} 
