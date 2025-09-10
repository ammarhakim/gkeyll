#include <gkyl_skin_surf_from_ghost_kernels.h> 

GKYL_CU_DH void skin_surf_from_ghost_lowerx_1x_ser_p1(const double *fghost, double *fskin) 
{ 
  // fghost: field in the ghost cell.
  // fskin: field in the skin cell.

  double fghostSurf[1];
  fghostSurf[0] = 1.224744871391589*fghost[1]+0.7071067811865475*fghost[0]; 

  double fskinNew[2];
  fskinNew[0] = 0.8660254037844386*fskin[1]+0.5*fskin[0]+0.7071067811865475*fghostSurf[0]; 
  fskinNew[1] = 0.5*fskin[1]+0.2886751345948129*fskin[0]-0.408248290463863*fghostSurf[0]; 

  fskin[0] = fskinNew[0]; 
  fskin[1] = fskinNew[1]; 

}

GKYL_CU_DH void skin_surf_from_ghost_upperx_1x_ser_p1(const double *fghost, double *fskin) 
{ 
  // fghost: field in the ghost cell.
  // fskin: field in the skin cell.

  double fghostSurf[1];
  fghostSurf[0] = 0.7071067811865475*fghost[0]-1.224744871391589*fghost[1]; 

  double fskinNew[2];
  fskinNew[0] = -(0.8660254037844386*fskin[1])+0.5*fskin[0]+0.7071067811865475*fghostSurf[0]; 
  fskinNew[1] = 0.5*fskin[1]-0.2886751345948129*fskin[0]+0.408248290463863*fghostSurf[0]; 

  fskin[0] = fskinNew[0]; 
  fskin[1] = fskinNew[1]; 

}

