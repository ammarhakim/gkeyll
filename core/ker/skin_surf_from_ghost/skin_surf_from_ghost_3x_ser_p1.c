#include <gkyl_skin_surf_from_ghost_kernels.h> 

GKYL_CU_DH void skin_surf_from_ghost_lowerx_3x_ser_p1(const double *fghost, double *fskin) 
{ 
  // fghost: field in the ghost cell.
  // fskin: field in the skin cell.

  double fghostSurf[4];
  fghostSurf[0] = 1.224744871391589*fghost[1]+0.7071067811865475*fghost[0]; 
  fghostSurf[1] = 1.224744871391589*fghost[4]+0.7071067811865475*fghost[2]; 
  fghostSurf[2] = 1.224744871391589*fghost[5]+0.7071067811865475*fghost[3]; 
  fghostSurf[3] = 1.224744871391589*fghost[7]+0.7071067811865475*fghost[6]; 

  double fskinNew[8];
  fskinNew[0] = 0.8660254037844386*fskin[1]+0.5*fskin[0]+0.7071067811865475*fghostSurf[0]; 
  fskinNew[1] = 0.5*fskin[1]+0.2886751345948129*fskin[0]-0.408248290463863*fghostSurf[0]; 
  fskinNew[2] = 0.8660254037844386*fskin[4]+0.5*fskin[2]+0.7071067811865475*fghostSurf[1]; 
  fskinNew[3] = 0.8660254037844386*fskin[5]+0.5*fskin[3]+0.7071067811865475*fghostSurf[2]; 
  fskinNew[4] = 0.5*fskin[4]+0.2886751345948129*fskin[2]-0.408248290463863*fghostSurf[1]; 
  fskinNew[5] = 0.5*fskin[5]+0.2886751345948129*fskin[3]-0.408248290463863*fghostSurf[2]; 
  fskinNew[6] = 0.8660254037844386*fskin[7]+0.5*fskin[6]+0.7071067811865475*fghostSurf[3]; 
  fskinNew[7] = 0.5*fskin[7]+0.2886751345948129*fskin[6]-0.408248290463863*fghostSurf[3]; 

  fskin[0] = fskinNew[0]; 
  fskin[1] = fskinNew[1]; 
  fskin[2] = fskinNew[2]; 
  fskin[3] = fskinNew[3]; 
  fskin[4] = fskinNew[4]; 
  fskin[5] = fskinNew[5]; 
  fskin[6] = fskinNew[6]; 
  fskin[7] = fskinNew[7]; 

}

GKYL_CU_DH void skin_surf_from_ghost_upperx_3x_ser_p1(const double *fghost, double *fskin) 
{ 
  // fghost: field in the ghost cell.
  // fskin: field in the skin cell.

  double fghostSurf[4];
  fghostSurf[0] = 0.7071067811865475*fghost[0]-1.224744871391589*fghost[1]; 
  fghostSurf[1] = 0.7071067811865475*fghost[2]-1.224744871391589*fghost[4]; 
  fghostSurf[2] = 0.7071067811865475*fghost[3]-1.224744871391589*fghost[5]; 
  fghostSurf[3] = 0.7071067811865475*fghost[6]-1.224744871391589*fghost[7]; 

  double fskinNew[8];
  fskinNew[0] = -(0.8660254037844386*fskin[1])+0.5*fskin[0]+0.7071067811865475*fghostSurf[0]; 
  fskinNew[1] = 0.5*fskin[1]-0.2886751345948129*fskin[0]+0.408248290463863*fghostSurf[0]; 
  fskinNew[2] = -(0.8660254037844386*fskin[4])+0.5*fskin[2]+0.7071067811865475*fghostSurf[1]; 
  fskinNew[3] = -(0.8660254037844386*fskin[5])+0.5*fskin[3]+0.7071067811865475*fghostSurf[2]; 
  fskinNew[4] = 0.5*fskin[4]-0.2886751345948129*fskin[2]+0.408248290463863*fghostSurf[1]; 
  fskinNew[5] = 0.5*fskin[5]-0.2886751345948129*fskin[3]+0.408248290463863*fghostSurf[2]; 
  fskinNew[6] = -(0.8660254037844386*fskin[7])+0.5*fskin[6]+0.7071067811865475*fghostSurf[3]; 
  fskinNew[7] = 0.5*fskin[7]-0.2886751345948129*fskin[6]+0.408248290463863*fghostSurf[3]; 

  fskin[0] = fskinNew[0]; 
  fskin[1] = fskinNew[1]; 
  fskin[2] = fskinNew[2]; 
  fskin[3] = fskinNew[3]; 
  fskin[4] = fskinNew[4]; 
  fskin[5] = fskinNew[5]; 
  fskin[6] = fskinNew[6]; 
  fskin[7] = fskinNew[7]; 

}

GKYL_CU_DH void skin_surf_from_ghost_lowery_3x_ser_p1(const double *fghost, double *fskin) 
{ 
  // fghost: field in the ghost cell.
  // fskin: field in the skin cell.

  double fghostSurf[4];
  fghostSurf[0] = 1.224744871391589*fghost[2]+0.7071067811865475*fghost[0]; 
  fghostSurf[1] = 1.224744871391589*fghost[4]+0.7071067811865475*fghost[1]; 
  fghostSurf[2] = 1.224744871391589*fghost[6]+0.7071067811865475*fghost[3]; 
  fghostSurf[3] = 1.224744871391589*fghost[7]+0.7071067811865475*fghost[5]; 

  double fskinNew[8];
  fskinNew[0] = 0.8660254037844386*fskin[2]+0.5*fskin[0]+0.7071067811865475*fghostSurf[0]; 
  fskinNew[1] = 0.8660254037844386*fskin[4]+0.5*fskin[1]+0.7071067811865475*fghostSurf[1]; 
  fskinNew[2] = 0.5*fskin[2]+0.2886751345948129*fskin[0]-0.408248290463863*fghostSurf[0]; 
  fskinNew[3] = 0.8660254037844386*fskin[6]+0.5*fskin[3]+0.7071067811865475*fghostSurf[2]; 
  fskinNew[4] = 0.5*fskin[4]+0.2886751345948129*fskin[1]-0.408248290463863*fghostSurf[1]; 
  fskinNew[5] = 0.8660254037844386*fskin[7]+0.5*fskin[5]+0.7071067811865475*fghostSurf[3]; 
  fskinNew[6] = 0.5*fskin[6]+0.2886751345948129*fskin[3]-0.408248290463863*fghostSurf[2]; 
  fskinNew[7] = 0.5*fskin[7]+0.2886751345948129*fskin[5]-0.408248290463863*fghostSurf[3]; 

  fskin[0] = fskinNew[0]; 
  fskin[1] = fskinNew[1]; 
  fskin[2] = fskinNew[2]; 
  fskin[3] = fskinNew[3]; 
  fskin[4] = fskinNew[4]; 
  fskin[5] = fskinNew[5]; 
  fskin[6] = fskinNew[6]; 
  fskin[7] = fskinNew[7]; 

}

GKYL_CU_DH void skin_surf_from_ghost_uppery_3x_ser_p1(const double *fghost, double *fskin) 
{ 
  // fghost: field in the ghost cell.
  // fskin: field in the skin cell.

  double fghostSurf[4];
  fghostSurf[0] = 0.7071067811865475*fghost[0]-1.224744871391589*fghost[2]; 
  fghostSurf[1] = 0.7071067811865475*fghost[1]-1.224744871391589*fghost[4]; 
  fghostSurf[2] = 0.7071067811865475*fghost[3]-1.224744871391589*fghost[6]; 
  fghostSurf[3] = 0.7071067811865475*fghost[5]-1.224744871391589*fghost[7]; 

  double fskinNew[8];
  fskinNew[0] = -(0.8660254037844386*fskin[2])+0.5*fskin[0]+0.7071067811865475*fghostSurf[0]; 
  fskinNew[1] = -(0.8660254037844386*fskin[4])+0.5*fskin[1]+0.7071067811865475*fghostSurf[1]; 
  fskinNew[2] = 0.5*fskin[2]-0.2886751345948129*fskin[0]+0.408248290463863*fghostSurf[0]; 
  fskinNew[3] = -(0.8660254037844386*fskin[6])+0.5*fskin[3]+0.7071067811865475*fghostSurf[2]; 
  fskinNew[4] = 0.5*fskin[4]-0.2886751345948129*fskin[1]+0.408248290463863*fghostSurf[1]; 
  fskinNew[5] = -(0.8660254037844386*fskin[7])+0.5*fskin[5]+0.7071067811865475*fghostSurf[3]; 
  fskinNew[6] = 0.5*fskin[6]-0.2886751345948129*fskin[3]+0.408248290463863*fghostSurf[2]; 
  fskinNew[7] = 0.5*fskin[7]-0.2886751345948129*fskin[5]+0.408248290463863*fghostSurf[3]; 

  fskin[0] = fskinNew[0]; 
  fskin[1] = fskinNew[1]; 
  fskin[2] = fskinNew[2]; 
  fskin[3] = fskinNew[3]; 
  fskin[4] = fskinNew[4]; 
  fskin[5] = fskinNew[5]; 
  fskin[6] = fskinNew[6]; 
  fskin[7] = fskinNew[7]; 

}

GKYL_CU_DH void skin_surf_from_ghost_lowerz_3x_ser_p1(const double *fghost, double *fskin) 
{ 
  // fghost: field in the ghost cell.
  // fskin: field in the skin cell.

  double fghostSurf[4];
  fghostSurf[0] = 1.224744871391589*fghost[3]+0.7071067811865475*fghost[0]; 
  fghostSurf[1] = 1.224744871391589*fghost[5]+0.7071067811865475*fghost[1]; 
  fghostSurf[2] = 1.224744871391589*fghost[6]+0.7071067811865475*fghost[2]; 
  fghostSurf[3] = 1.224744871391589*fghost[7]+0.7071067811865475*fghost[4]; 

  double fskinNew[8];
  fskinNew[0] = 0.8660254037844386*fskin[3]+0.5*fskin[0]+0.7071067811865475*fghostSurf[0]; 
  fskinNew[1] = 0.8660254037844386*fskin[5]+0.5*fskin[1]+0.7071067811865475*fghostSurf[1]; 
  fskinNew[2] = 0.8660254037844386*fskin[6]+0.5*fskin[2]+0.7071067811865475*fghostSurf[2]; 
  fskinNew[3] = 0.5*fskin[3]+0.2886751345948129*fskin[0]-0.408248290463863*fghostSurf[0]; 
  fskinNew[4] = 0.8660254037844386*fskin[7]+0.5*fskin[4]+0.7071067811865475*fghostSurf[3]; 
  fskinNew[5] = 0.5*fskin[5]+0.2886751345948129*fskin[1]-0.408248290463863*fghostSurf[1]; 
  fskinNew[6] = 0.5*fskin[6]+0.2886751345948129*fskin[2]-0.408248290463863*fghostSurf[2]; 
  fskinNew[7] = 0.5*fskin[7]+0.2886751345948129*fskin[4]-0.408248290463863*fghostSurf[3]; 

  fskin[0] = fskinNew[0]; 
  fskin[1] = fskinNew[1]; 
  fskin[2] = fskinNew[2]; 
  fskin[3] = fskinNew[3]; 
  fskin[4] = fskinNew[4]; 
  fskin[5] = fskinNew[5]; 
  fskin[6] = fskinNew[6]; 
  fskin[7] = fskinNew[7]; 

}

GKYL_CU_DH void skin_surf_from_ghost_upperz_3x_ser_p1(const double *fghost, double *fskin) 
{ 
  // fghost: field in the ghost cell.
  // fskin: field in the skin cell.

  double fghostSurf[4];
  fghostSurf[0] = 0.7071067811865475*fghost[0]-1.224744871391589*fghost[3]; 
  fghostSurf[1] = 0.7071067811865475*fghost[1]-1.224744871391589*fghost[5]; 
  fghostSurf[2] = 0.7071067811865475*fghost[2]-1.224744871391589*fghost[6]; 
  fghostSurf[3] = 0.7071067811865475*fghost[4]-1.224744871391589*fghost[7]; 

  double fskinNew[8];
  fskinNew[0] = -(0.8660254037844386*fskin[3])+0.5*fskin[0]+0.7071067811865475*fghostSurf[0]; 
  fskinNew[1] = -(0.8660254037844386*fskin[5])+0.5*fskin[1]+0.7071067811865475*fghostSurf[1]; 
  fskinNew[2] = -(0.8660254037844386*fskin[6])+0.5*fskin[2]+0.7071067811865475*fghostSurf[2]; 
  fskinNew[3] = 0.5*fskin[3]-0.2886751345948129*fskin[0]+0.408248290463863*fghostSurf[0]; 
  fskinNew[4] = -(0.8660254037844386*fskin[7])+0.5*fskin[4]+0.7071067811865475*fghostSurf[3]; 
  fskinNew[5] = 0.5*fskin[5]-0.2886751345948129*fskin[1]+0.408248290463863*fghostSurf[1]; 
  fskinNew[6] = 0.5*fskin[6]-0.2886751345948129*fskin[2]+0.408248290463863*fghostSurf[2]; 
  fskinNew[7] = 0.5*fskin[7]-0.2886751345948129*fskin[4]+0.408248290463863*fghostSurf[3]; 

  fskin[0] = fskinNew[0]; 
  fskin[1] = fskinNew[1]; 
  fskin[2] = fskinNew[2]; 
  fskin[3] = fskinNew[3]; 
  fskin[4] = fskinNew[4]; 
  fskin[5] = fskinNew[5]; 
  fskin[6] = fskinNew[6]; 
  fskin[7] = fskinNew[7]; 

}

