#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_M0_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, 
  const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  out[0] += 1.4142135623730951*f[0]*volFact; 
  out[1] += 1.4142135623730951*f[1]*volFact; 
  out[2] += 1.4142135623730951*f[4]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M2ij_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double wx1 = w[1], dv1 = dxv[1]; 
  double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += volFact*(1.4142135623730951*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.10540925533894596*f[5]*dv1_sq+0.11785113019775789*f[0]*dv1_sq); 
  out[1] += volFact*(1.4142135623730951*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.105409255338946*f[7]*dv1_sq+0.11785113019775789*f[1]*dv1_sq); 
  out[2] += volFact*(1.4142135623730951*f[4]*wx1_sq+0.816496580927726*f[6]*dv1*wx1+0.10540925533894596*f[8]*dv1_sq+0.11785113019775789*f[4]*dv1_sq); 
} 
GKYL_CU_DH void mom_vlasov_M3ijk_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double wx1 = w[1], dv1 = dxv[1]; 
  double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  out[0] += volFact*(1.224744871391589*f[2]*dv1*wx1_sq+1.4142135623730951*f[0]*wx1_cu+0.3162277660168379*f[5]*dv1_sq*wx1+0.3535533905932737*f[0]*dv1_sq*wx1+0.06123724356957942*f[2]*dv1_cu); 
  out[1] += volFact*(1.224744871391589*f[3]*dv1*wx1_sq+1.4142135623730951*f[1]*wx1_cu+0.31622776601683794*f[7]*dv1_sq*wx1+0.3535533905932737*f[1]*dv1_sq*wx1+0.06123724356957942*f[3]*dv1_cu); 
  out[2] += volFact*(1.224744871391589*f[6]*dv1*wx1_sq+1.4142135623730951*f[4]*wx1_cu+0.3162277660168379*f[8]*dv1_sq*wx1+0.3535533905932737*f[4]*dv1_sq*wx1+0.06123724356957942*f[6]*dv1_cu); 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_M1i_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  out[0] += ((1.2909944487358058*hamil[2]*f[5]+0.4303314829119353*hamil[1]*f[5]+1.9364916731037085*f[2]*hamil[2]+1.4433756729740645*f[0]*hamil[2]+0.6454972243679029*hamil[1]*f[2]+0.48112522432468824*f[0]*hamil[1])*dv10*volFact)/jacob_vx[2]+((-(1.2909944487358058*hamil[2]*f[5])+0.4303314829119353*hamil[1]*f[5]+1.9364916731037085*f[2]*hamil[2]-1.4433756729740645*f[0]*hamil[2]-0.6454972243679029*hamil[1]*f[2]+0.48112522432468824*f[0]*hamil[1])*dv10*volFact)/jacob_vx[0]+((0.7698003589195012*f[0]*hamil[1]-0.8606629658238706*hamil[1]*f[5])*dv10*volFact)/jacob_vx[1]; 
  out[1] += ((1.2909944487358056*hamil[2]*f[7]+0.4303314829119352*hamil[1]*f[7]+1.9364916731037085*hamil[2]*f[3]+0.6454972243679029*hamil[1]*f[3]+1.4433756729740645*f[1]*hamil[2]+0.48112522432468824*f[1]*hamil[1])*dv10*volFact)/jacob_vx[2]+((-(1.2909944487358056*hamil[2]*f[7])+0.4303314829119352*hamil[1]*f[7]+1.9364916731037085*hamil[2]*f[3]-0.6454972243679029*hamil[1]*f[3]-1.4433756729740645*f[1]*hamil[2]+0.48112522432468824*f[1]*hamil[1])*dv10*volFact)/jacob_vx[0]+((0.7698003589195012*f[1]*hamil[1]-0.8606629658238704*hamil[1]*f[7])*dv10*volFact)/jacob_vx[1]; 
  out[2] += ((1.2909944487358058*hamil[2]*f[8]+0.4303314829119353*hamil[1]*f[8]+1.9364916731037085*hamil[2]*f[6]+0.6454972243679028*hamil[1]*f[6]+1.4433756729740645*hamil[2]*f[4]+0.48112522432468824*hamil[1]*f[4])*dv10*volFact)/jacob_vx[2]+((-(1.2909944487358058*hamil[2]*f[8])+0.4303314829119353*hamil[1]*f[8]+1.9364916731037085*hamil[2]*f[6]-0.6454972243679028*hamil[1]*f[6]-1.4433756729740645*hamil[2]*f[4]+0.48112522432468824*hamil[1]*f[4])*dv10*volFact)/jacob_vx[0]+((0.7698003589195012*hamil[1]*f[4]-0.8606629658238706*hamil[1]*f[8])*dv10*volFact)/jacob_vx[1]; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_M2_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  out[0] += (hamil[2]*f[5]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[1] += (1.0000000000000002*hamil[2]*f[7]+hamil[1]*f[3]+hamil[0]*f[1])*volFact; 
  out[2] += (hamil[2]*f[8]+1.0000000000000002*hamil[1]*f[6]+hamil[0]*f[4])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_five_moments_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  out[0] += 1.4142135623730951*f[0]*volFact; 
  out[1] += 1.4142135623730951*f[1]*volFact; 
  out[2] += 1.4142135623730951*f[4]*volFact; 
  out[3] += ((1.2909944487358058*hamil[2]*f[5]+0.4303314829119353*hamil[1]*f[5]+1.9364916731037085*f[2]*hamil[2]+1.4433756729740645*f[0]*hamil[2]+0.6454972243679029*hamil[1]*f[2]+0.48112522432468824*f[0]*hamil[1])*dv10*volFact)/jacob_vx[2]+((-(1.2909944487358058*hamil[2]*f[5])+0.4303314829119353*hamil[1]*f[5]+1.9364916731037085*f[2]*hamil[2]-1.4433756729740645*f[0]*hamil[2]-0.6454972243679029*hamil[1]*f[2]+0.48112522432468824*f[0]*hamil[1])*dv10*volFact)/jacob_vx[0]+((0.7698003589195012*f[0]*hamil[1]-0.8606629658238706*hamil[1]*f[5])*dv10*volFact)/jacob_vx[1]; 
  out[4] += ((1.2909944487358056*hamil[2]*f[7]+0.4303314829119352*hamil[1]*f[7]+1.9364916731037085*hamil[2]*f[3]+0.6454972243679029*hamil[1]*f[3]+1.4433756729740645*f[1]*hamil[2]+0.48112522432468824*f[1]*hamil[1])*dv10*volFact)/jacob_vx[2]+((-(1.2909944487358056*hamil[2]*f[7])+0.4303314829119352*hamil[1]*f[7]+1.9364916731037085*hamil[2]*f[3]-0.6454972243679029*hamil[1]*f[3]-1.4433756729740645*f[1]*hamil[2]+0.48112522432468824*f[1]*hamil[1])*dv10*volFact)/jacob_vx[0]+((0.7698003589195012*f[1]*hamil[1]-0.8606629658238704*hamil[1]*f[7])*dv10*volFact)/jacob_vx[1]; 
  out[5] += ((1.2909944487358058*hamil[2]*f[8]+0.4303314829119353*hamil[1]*f[8]+1.9364916731037085*hamil[2]*f[6]+0.6454972243679028*hamil[1]*f[6]+1.4433756729740645*hamil[2]*f[4]+0.48112522432468824*hamil[1]*f[4])*dv10*volFact)/jacob_vx[2]+((-(1.2909944487358058*hamil[2]*f[8])+0.4303314829119353*hamil[1]*f[8]+1.9364916731037085*hamil[2]*f[6]-0.6454972243679028*hamil[1]*f[6]-1.4433756729740645*hamil[2]*f[4]+0.48112522432468824*hamil[1]*f[4])*dv10*volFact)/jacob_vx[0]+((0.7698003589195012*hamil[1]*f[4]-0.8606629658238706*hamil[1]*f[8])*dv10*volFact)/jacob_vx[1]; 
  out[6] += (hamil[2]*f[5]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[7] += (1.0000000000000002*hamil[2]*f[7]+hamil[1]*f[3]+hamil[0]*f[1])*volFact; 
  out[8] += (hamil[2]*f[8]+1.0000000000000002*hamil[1]*f[6]+hamil[0]*f[4])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_int_five_moments_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*0.25; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += ((1.825741858350554*hamil[2]*f[5]+0.6085806194501847*hamil[1]*f[5]+2.7386127875258306*f[2]*hamil[2]+2.041241452319315*f[0]*hamil[2]+0.9128709291752769*hamil[1]*f[2]+0.6804138174397718*f[0]*hamil[1])*dv10*volFact)/jacob_vx[2]+((-(1.825741858350554*hamil[2]*f[5])+0.6085806194501847*hamil[1]*f[5]+2.7386127875258306*f[2]*hamil[2]-2.041241452319315*f[0]*hamil[2]-0.9128709291752769*hamil[1]*f[2]+0.6804138174397718*f[0]*hamil[1])*dv10*volFact)/jacob_vx[0]+((1.0886621079036354*f[0]*hamil[1]-1.2171612389003696*hamil[1]*f[5])*dv10*volFact)/jacob_vx[1]; 
  out[2] += (1.4142135623730951*hamil[2]*f[5]+1.4142135623730951*hamil[1]*f[2]+1.4142135623730951*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M0_upper_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, 
    const double *vmap, double v_thresh, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double *px = &vmap[0]; 
  double f_nodes[9] = {0.0};
  if ((0.6324555320336759*px[2]-0.9486832980505137*px[1]+0.7071067811865475*px[0]) > v_thresh) { 
    f_nodes[0] = 0.4*f[8]-0.5999999999999995*f[7]-0.5999999999999999*f[6]+0.4472135954999579*(f[5]+f[4])+0.9*f[3]-0.6708203932499369*(f[2]+f[1])+0.5*f[0]; 
    f_nodes[3] = -(0.5*f[8])+0.75*f[6]+0.4472135954999579*f[5]-0.5590169943749475*f[4]-0.6708203932499369*f[2]+0.5*f[0]; 
    f_nodes[6] = 0.4*f[8]+0.5999999999999995*f[7]-0.5999999999999999*f[6]+0.4472135954999579*(f[5]+f[4])-0.9*f[3]-0.6708203932499369*f[2]+0.6708203932499369*f[1]+0.5*f[0]; 
  } 
  if ((0.7071067811865475*px[0]-0.7905694150420947*px[2]) > v_thresh) { 
    f_nodes[1] = -(0.5*f[8])+0.75*f[7]-0.5590169943749475*f[5]+0.4472135954999579*f[4]-0.6708203932499369*f[1]+0.5*f[0]; 
    f_nodes[4] = 0.625*f[8]-0.5590169943749475*(f[5]+f[4])+0.5*f[0]; 
    f_nodes[7] = -(0.5*f[8])-0.75*f[7]-0.5590169943749475*f[5]+0.4472135954999579*f[4]+0.6708203932499369*f[1]+0.5*f[0]; 
  } 
  if ((0.6324555320336759*px[2]+0.9486832980505137*px[1]+0.7071067811865475*px[0]) > v_thresh) { 
    f_nodes[2] = 0.4*f[8]-0.5999999999999995*f[7]+0.5999999999999999*f[6]+0.4472135954999579*(f[5]+f[4])-0.9*f[3]+0.6708203932499369*f[2]-0.6708203932499369*f[1]+0.5*f[0]; 
    f_nodes[5] = -(0.5*f[8])-0.75*f[6]+0.4472135954999579*f[5]-0.5590169943749475*f[4]+0.6708203932499369*f[2]+0.5*f[0]; 
    f_nodes[8] = 0.4*f[8]+0.5999999999999995*f[7]+0.5999999999999999*f[6]+0.4472135954999579*(f[5]+f[4])+0.9*f[3]+0.6708203932499369*(f[2]+f[1])+0.5*f[0]; 
  } 
  out[0] += (0.21824283369955166*f_nodes[8]+0.34918853391928284*f_nodes[7]+0.21824283369955166*f_nodes[6]+0.34918853391928284*f_nodes[5]+0.5587016542708527*f_nodes[4]+0.34918853391928284*f_nodes[3]+0.21824283369955166*f_nodes[2]+0.34918853391928284*f_nodes[1]+0.21824283369955166*f_nodes[0])*volFact; 
  out[1] += (0.2928034870526277*f_nodes[8]+0.46848557928420453*f_nodes[7]+0.2928034870526277*f_nodes[6]-0.2928034870526277*f_nodes[2]-0.46848557928420453*f_nodes[1]-0.2928034870526277*f_nodes[0])*volFact; 
  out[2] += (0.19520232470175186*f_nodes[8]+0.312323719522803*f_nodes[7]+0.19520232470175186*f_nodes[6]-0.3904046494035038*f_nodes[5]-0.624647439045606*f_nodes[4]-0.3904046494035038*f_nodes[3]+0.19520232470175186*f_nodes[2]+0.312323719522803*f_nodes[1]+0.19520232470175186*f_nodes[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M0_lower_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, 
    const double *vmap, double v_thresh, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double *px = &vmap[0]; 
  double f_nodes[9] = {0.0};
  if ((0.6324555320336759*px[2]-0.9486832980505137*px[1]+0.7071067811865475*px[0]) < -v_thresh) { 
    f_nodes[0] = 0.4*f[8]-0.5999999999999995*f[7]-0.5999999999999999*f[6]+0.4472135954999579*(f[5]+f[4])+0.9*f[3]-0.6708203932499369*(f[2]+f[1])+0.5*f[0]; 
    f_nodes[3] = -(0.5*f[8])+0.75*f[6]+0.4472135954999579*f[5]-0.5590169943749475*f[4]-0.6708203932499369*f[2]+0.5*f[0]; 
    f_nodes[6] = 0.4*f[8]+0.5999999999999995*f[7]-0.5999999999999999*f[6]+0.4472135954999579*(f[5]+f[4])-0.9*f[3]-0.6708203932499369*f[2]+0.6708203932499369*f[1]+0.5*f[0]; 
  } 
  if ((0.7071067811865475*px[0]-0.7905694150420947*px[2]) < -v_thresh) { 
    f_nodes[1] = -(0.5*f[8])+0.75*f[7]-0.5590169943749475*f[5]+0.4472135954999579*f[4]-0.6708203932499369*f[1]+0.5*f[0]; 
    f_nodes[4] = 0.625*f[8]-0.5590169943749475*(f[5]+f[4])+0.5*f[0]; 
    f_nodes[7] = -(0.5*f[8])-0.75*f[7]-0.5590169943749475*f[5]+0.4472135954999579*f[4]+0.6708203932499369*f[1]+0.5*f[0]; 
  } 
  if ((0.6324555320336759*px[2]+0.9486832980505137*px[1]+0.7071067811865475*px[0]) < -v_thresh) { 
    f_nodes[2] = 0.4*f[8]-0.5999999999999995*f[7]+0.5999999999999999*f[6]+0.4472135954999579*(f[5]+f[4])-0.9*f[3]+0.6708203932499369*f[2]-0.6708203932499369*f[1]+0.5*f[0]; 
    f_nodes[5] = -(0.5*f[8])-0.75*f[6]+0.4472135954999579*f[5]-0.5590169943749475*f[4]+0.6708203932499369*f[2]+0.5*f[0]; 
    f_nodes[8] = 0.4*f[8]+0.5999999999999995*f[7]+0.5999999999999999*f[6]+0.4472135954999579*(f[5]+f[4])+0.9*f[3]+0.6708203932499369*(f[2]+f[1])+0.5*f[0]; 
  } 
  out[0] += (0.21824283369955166*f_nodes[8]+0.34918853391928284*f_nodes[7]+0.21824283369955166*f_nodes[6]+0.34918853391928284*f_nodes[5]+0.5587016542708527*f_nodes[4]+0.34918853391928284*f_nodes[3]+0.21824283369955166*f_nodes[2]+0.34918853391928284*f_nodes[1]+0.21824283369955166*f_nodes[0])*volFact; 
  out[1] += (0.2928034870526277*f_nodes[8]+0.46848557928420453*f_nodes[7]+0.2928034870526277*f_nodes[6]-0.2928034870526277*f_nodes[2]-0.46848557928420453*f_nodes[1]-0.2928034870526277*f_nodes[0])*volFact; 
  out[2] += (0.19520232470175186*f_nodes[8]+0.312323719522803*f_nodes[7]+0.19520232470175186*f_nodes[6]-0.3904046494035038*f_nodes[5]-0.624647439045606*f_nodes[4]-0.3904046494035038*f_nodes[3]+0.19520232470175186*f_nodes[2]+0.312323719522803*f_nodes[1]+0.19520232470175186*f_nodes[0])*volFact; 
} 
