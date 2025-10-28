#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_1x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  out[0] += 0.0125*(53.66563145999496*coeff[0]*qr[2]+53.66563145999496*coeff[0]*ql[2]-107.3312629199899*coeff[0]*qc[2]-95.26279441628824*coeff[0]*qr[1]+95.26279441628824*coeff[0]*ql[1]+75.0*coeff[0]*qr[0]+75.0*coeff[0]*ql[0]-150.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.003125*(236.2519841186524*coeff[0]*qr[2]-236.2519841186524*coeff[0]*ql[2]-465.0*coeff[0]*qr[1]-465.0*coeff[0]*ql[1]-1710.0*coeff[0]*qc[1]+381.051177665153*coeff[0]*qr[0]-381.051177665153*coeff[0]*ql[0])*Jfac; 
  out[2] += -0.015625*(9.0*coeff[0]*qr[2]+9.0*coeff[0]*ql[2]+402.0*coeff[0]*qc[2]+19.36491673103709*coeff[0]*qr[1]-19.36491673103709*coeff[0]*ql[1]-26.83281572999748*coeff[0]*qr[0]-26.83281572999748*coeff[0]*ql[0]+53.66563145999496*coeff[0]*qc[0])*Jfac; 

  return 0.;

}

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_1x_ser_p2_varcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  out[0] += 0.008838834764831839*((120.0*coeff[2]+92.951600308978*coeff[1]+53.66563145999496*coeff[0])*qr[2]+(120.0*coeff[2]-92.951600308978*coeff[1]+53.66563145999496*coeff[0])*ql[2]+((-240.0*coeff[2])-107.3312629199899*coeff[0])*qc[2]+((-213.014084041408*qr[1])+213.014084041408*ql[1]+167.7050983124843*qr[0]+167.7050983124843*ql[0]-335.4101966249685*qc[0])*coeff[2]+((-165.0*coeff[1])-95.26279441628824*coeff[0])*qr[1]+(95.26279441628824*coeff[0]-165.0*coeff[1])*ql[1]-330.0*coeff[1]*qc[1]+(129.9038105676658*qr[0]-129.9038105676658*ql[0])*coeff[1]+75.0*coeff[0]*qr[0]+75.0*coeff[0]*ql[0]-150.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.002209708691207959*((528.2754963085075*coeff[2]+409.2004398824615*coeff[1]+236.2519841186524*coeff[0])*qr[2]+((-528.2754963085075*coeff[2])+409.2004398824615*coeff[1]-236.2519841186524*coeff[0])*ql[2]-1757.549430314835*coeff[1]*qc[2]+((-1039.771609537402*qr[1])-1039.771609537402*ql[1]-1677.050983124843*qc[1]+852.0563361656318*qr[0]-852.0563361656318*ql[0])*coeff[2]+((-805.4036255195278*coeff[1])-465.0*coeff[0])*qr[1]+(805.4036255195278*coeff[1]-465.0*coeff[0])*ql[1]-1710.0*coeff[0]*qc[1]+(660.0*qr[0]+660.0*ql[0]-1320.0*qc[0])*coeff[1]+381.051177665153*coeff[0]*qr[0]-381.051177665153*coeff[0]*ql[0])*Jfac; 
  out[2] += -0.0110485434560398*((20.12461179749811*coeff[2]+15.58845726811989*coeff[1]+9.0*coeff[0])*qr[2]+(20.12461179749811*coeff[2]-15.58845726811989*coeff[1]+9.0*coeff[0])*ql[2]+(402.0*coeff[0]-389.0758280849634*coeff[2])*qc[2]+(43.30127018922193*qr[1]-43.30127018922193*ql[1]-60.0*qr[0]-60.0*ql[0]+120.0*qc[0])*coeff[2]+(33.54101966249685*coeff[1]+19.36491673103709*coeff[0])*qr[1]+(33.54101966249685*coeff[1]-19.36491673103709*coeff[0])*ql[1]+254.911749434976*coeff[1]*qc[1]+(46.475800154489*ql[0]-46.475800154489*qr[0])*coeff[1]-26.83281572999748*coeff[0]*qr[0]-26.83281572999748*coeff[0]*ql[0]+53.66563145999496*coeff[0]*qc[0])*Jfac; 

  return 0.;

}

