#include <gkyl_vlasov_poisson_kernels.h> 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_1x1v_ser_p1(const double *w, const double *dxv, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // edge: determines if the update is for the left edge (-1) or right edge (+1).
  // fskin,fedge: distribution function in skin cell/last edge cell.
  // out: output increment in center cell.

  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[1], wv = w[1]; 
  double Ghat[3]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = (1.224744871391589*fskin[1]+0.7071067811865475*fskin[0])*wv+(0.3535533905932737*fskin[3]+0.20412414523193148*fskin[2])*dv; 
  Ghat[1] = (1.224744871391589*fskin[3]+0.7071067811865475*fskin[2])*wv+(0.3162277660168379*fskin[5]+0.18257418583505536*fskin[4]+0.3535533905932737*fskin[1]+0.20412414523193148*fskin[0])*dv; 
  Ghat[2] = (1.224744871391589*fskin[5]+0.7071067811865475*fskin[4])*wv+(0.3162277660168379*fskin[3]+0.18257418583505536*fskin[2])*dv; 

  } else { 

  Ghat[0] = -(0.11785113019775789*((10.392304845413262*fedge[1]-6.0*fedge[0])*wv+(3.0*fedge[3]-1.7320508075688772*fedge[2])*dv)); 
  Ghat[1] = -(0.023570226039551577*((51.96152422706631*fedge[3]-30.0*fedge[2])*wv+(13.416407864998737*fedge[5]-7.745966692414834*fedge[4]+15.0*fedge[1]-8.660254037844386*fedge[0])*dv)); 
  Ghat[2] = -(0.04714045207910316*((25.98076211353316*fedge[5]-15.0*fedge[4])*wv+(6.708203932499369*fedge[3]-3.872983346207417*fedge[2])*dv)); 

  } 

  out[0] += -(0.7071067811865475*Ghat[0]*dx10); 
  out[1] += -(1.224744871391589*Ghat[0]*dx10); 
  out[2] += -(0.7071067811865475*Ghat[1]*dx10); 
  out[3] += -(1.224744871391589*Ghat[1]*dx10); 
  out[4] += -(0.7071067811865475*Ghat[2]*dx10); 
  out[5] += -(1.224744871391589*Ghat[2]*dx10); 

  } else { 

  if (wv>0) { 

  Ghat[0] = (1.224744871391589*fedge[1]+0.7071067811865475*fedge[0])*wv+(0.3535533905932737*fedge[3]+0.20412414523193148*fedge[2])*dv; 
  Ghat[1] = (1.224744871391589*fedge[3]+0.7071067811865475*fedge[2])*wv+(0.3162277660168379*fedge[5]+0.18257418583505536*fedge[4]+0.3535533905932737*fedge[1]+0.20412414523193148*fedge[0])*dv; 
  Ghat[2] = (1.224744871391589*fedge[5]+0.7071067811865475*fedge[4])*wv+(0.3162277660168379*fedge[3]+0.18257418583505536*fedge[2])*dv; 

  } else { 

  Ghat[0] = -(0.11785113019775789*((10.392304845413262*fskin[1]-6.0*fskin[0])*wv+(3.0*fskin[3]-1.7320508075688772*fskin[2])*dv)); 
  Ghat[1] = -(0.023570226039551577*((51.96152422706631*fskin[3]-30.0*fskin[2])*wv+(13.416407864998737*fskin[5]-7.745966692414834*fskin[4]+15.0*fskin[1]-8.660254037844386*fskin[0])*dv)); 
  Ghat[2] = -(0.04714045207910316*((25.98076211353316*fskin[5]-15.0*fskin[4])*wv+(6.708203932499369*fskin[3]-3.872983346207417*fskin[2])*dv)); 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -(1.224744871391589*Ghat[0]*dx10); 
  out[2] += 0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -(1.224744871391589*Ghat[1]*dx10); 
  out[4] += 0.7071067811865475*Ghat[2]*dx10; 
  out[5] += -(1.224744871391589*Ghat[2]*dx10); 

  } 
  return 0.;

} 
