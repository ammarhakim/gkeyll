GKYL_CU_DH static inline void 
ser_2x_p4_upwind_quad_to_modal(const double* fUpwindQuad, double* GKYL_RESTRICT fUpwind) { 
  fUpwind[0] = 0.1675326070686369*fUpwindQuad[4]+0.3384415785804036*fUpwindQuad[3]+0.4022651910750141*fUpwindQuad[2]+0.3384415785804036*fUpwindQuad[1]+0.1675326070686369*fUpwindQuad[0]; 
  fUpwind[1] = 0.262950725347801*fUpwindQuad[4]+0.315649637758137*fUpwindQuad[3]-0.315649637758137*fUpwindQuad[1]-0.262950725347801*fUpwindQuad[0]; 
  fUpwind[2] = 0.2741213413710451*fUpwindQuad[4]-0.04924826331462695*fUpwindQuad[3]-0.4497461561128364*fUpwindQuad[2]-0.04924826331462695*fUpwindQuad[1]+0.2741213413710451*fUpwindQuad[0]; 
  fUpwind[3] = 0.2220818735672946*fUpwindQuad[4]-0.3737373963531613*fUpwindQuad[3]+0.3737373963531613*fUpwindQuad[1]-0.2220818735672946*fUpwindQuad[0]; 
  fUpwind[4] = 0.123506106334137*fUpwindQuad[4]-0.3497802763138322*fUpwindQuad[3]+0.4525483399593906*fUpwindQuad[2]-0.3497802763138322*fUpwindQuad[1]+0.123506106334137*fUpwindQuad[0]; 

} 
