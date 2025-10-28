// Thu Jul  7 14:00:20 2022
#include <gkyl_basis_gkhyb_kernels.h>
GKYL_CU_DH
void
flip_odd_sign_1x1v_gkhyb_p1(int dir, const double *f, double *fout )
{
  if (dir == 0) {
    fout[0] = 1*f[0];
    fout[1] = -1*f[1];
    fout[2] = 1*f[2];
    fout[3] = -1*f[3];
    fout[4] = 1*f[4];
    fout[5] = -1*f[5];
  }
  if (dir == 1) {
    fout[0] = 1*f[0];
    fout[1] = 1*f[1];
    fout[2] = -1*f[2];
    fout[3] = -1*f[3];
    fout[4] = 1*f[4];
    fout[5] = 1*f[5];
  }
}

GKYL_CU_DH
void
flip_even_sign_1x1v_gkhyb_p1(int dir, const double *f, double *fout )
{
  if (dir == 0) {
    fout[0] = -1*f[0];
    fout[1] = 1*f[1];
    fout[2] = -1*f[2];
    fout[3] = 1*f[3];
    fout[4] = -1*f[4];
    fout[5] = 1*f[5];
  }
  if (dir == 1) {
    fout[0] = -1*f[0];
    fout[1] = -1*f[1];
    fout[2] = 1*f[2];
    fout[3] = 1*f[3];
    fout[4] = -1*f[4];
    fout[5] = -1*f[5];
  }
}

GKYL_CU_DH
void
flip_odd_sign_1x2v_gkhyb_p1(int dir, const double *f, double *fout )
{
  if (dir == 0) {
    fout[0] = 1*f[0];
    fout[1] = -1*f[1];
    fout[2] = 1*f[2];
    fout[3] = 1*f[3];
    fout[4] = -1*f[4];
    fout[5] = -1*f[5];
    fout[6] = 1*f[6];
    fout[7] = -1*f[7];
    fout[8] = 1*f[8];
    fout[9] = -1*f[9];
    fout[10] = 1*f[10];
    fout[11] = -1*f[11];
  }
  if (dir == 1) {
    fout[0] = 1*f[0];
    fout[1] = 1*f[1];
    fout[2] = -1*f[2];
    fout[3] = 1*f[3];
    fout[4] = -1*f[4];
    fout[5] = 1*f[5];
    fout[6] = -1*f[6];
    fout[7] = -1*f[7];
    fout[8] = 1*f[8];
    fout[9] = 1*f[9];
    fout[10] = 1*f[10];
    fout[11] = 1*f[11];
  }
  if (dir == 2) {
    fout[0] = 1*f[0];
    fout[1] = 1*f[1];
    fout[2] = 1*f[2];
    fout[3] = -1*f[3];
    fout[4] = 1*f[4];
    fout[5] = -1*f[5];
    fout[6] = -1*f[6];
    fout[7] = -1*f[7];
    fout[8] = 1*f[8];
    fout[9] = 1*f[9];
    fout[10] = -1*f[10];
    fout[11] = -1*f[11];
  }
}

GKYL_CU_DH
void
flip_even_sign_1x2v_gkhyb_p1(int dir, const double *f, double *fout )
{
  if (dir == 0) {
    fout[0] = -1*f[0];
    fout[1] = 1*f[1];
    fout[2] = -1*f[2];
    fout[3] = -1*f[3];
    fout[4] = 1*f[4];
    fout[5] = 1*f[5];
    fout[6] = -1*f[6];
    fout[7] = 1*f[7];
    fout[8] = -1*f[8];
    fout[9] = 1*f[9];
    fout[10] = -1*f[10];
    fout[11] = 1*f[11];
  }
  if (dir == 1) {
    fout[0] = -1*f[0];
    fout[1] = -1*f[1];
    fout[2] = 1*f[2];
    fout[3] = -1*f[3];
    fout[4] = 1*f[4];
    fout[5] = -1*f[5];
    fout[6] = 1*f[6];
    fout[7] = 1*f[7];
    fout[8] = -1*f[8];
    fout[9] = -1*f[9];
    fout[10] = -1*f[10];
    fout[11] = -1*f[11];
  }
  if (dir == 2) {
    fout[0] = -1*f[0];
    fout[1] = -1*f[1];
    fout[2] = -1*f[2];
    fout[3] = 1*f[3];
    fout[4] = -1*f[4];
    fout[5] = 1*f[5];
    fout[6] = 1*f[6];
    fout[7] = 1*f[7];
    fout[8] = -1*f[8];
    fout[9] = -1*f[9];
    fout[10] = 1*f[10];
    fout[11] = 1*f[11];
  }
}

GKYL_CU_DH
void
flip_odd_sign_2x2v_gkhyb_p1(int dir, const double *f, double *fout )
{
  if (dir == 0) {
    fout[0] = 1*f[0];
    fout[1] = -1*f[1];
    fout[2] = 1*f[2];
    fout[3] = 1*f[3];
    fout[4] = 1*f[4];
    fout[5] = -1*f[5];
    fout[6] = -1*f[6];
    fout[7] = 1*f[7];
    fout[8] = -1*f[8];
    fout[9] = 1*f[9];
    fout[10] = 1*f[10];
    fout[11] = -1*f[11];
    fout[12] = -1*f[12];
    fout[13] = -1*f[13];
    fout[14] = 1*f[14];
    fout[15] = -1*f[15];
    fout[16] = 1*f[16];
    fout[17] = -1*f[17];
    fout[18] = 1*f[18];
    fout[19] = 1*f[19];
    fout[20] = -1*f[20];
    fout[21] = -1*f[21];
    fout[22] = 1*f[22];
    fout[23] = -1*f[23];
  }
  if (dir == 1) {
    fout[0] = 1*f[0];
    fout[1] = 1*f[1];
    fout[2] = -1*f[2];
    fout[3] = 1*f[3];
    fout[4] = 1*f[4];
    fout[5] = -1*f[5];
    fout[6] = 1*f[6];
    fout[7] = -1*f[7];
    fout[8] = 1*f[8];
    fout[9] = -1*f[9];
    fout[10] = 1*f[10];
    fout[11] = -1*f[11];
    fout[12] = -1*f[12];
    fout[13] = 1*f[13];
    fout[14] = -1*f[14];
    fout[15] = -1*f[15];
    fout[16] = 1*f[16];
    fout[17] = 1*f[17];
    fout[18] = -1*f[18];
    fout[19] = 1*f[19];
    fout[20] = -1*f[20];
    fout[21] = 1*f[21];
    fout[22] = -1*f[22];
    fout[23] = -1*f[23];
  }
  if (dir == 2) {
    fout[0] = 1*f[0];
    fout[1] = 1*f[1];
    fout[2] = 1*f[2];
    fout[3] = -1*f[3];
    fout[4] = 1*f[4];
    fout[5] = 1*f[5];
    fout[6] = -1*f[6];
    fout[7] = -1*f[7];
    fout[8] = 1*f[8];
    fout[9] = 1*f[9];
    fout[10] = -1*f[10];
    fout[11] = -1*f[11];
    fout[12] = 1*f[12];
    fout[13] = -1*f[13];
    fout[14] = -1*f[14];
    fout[15] = -1*f[15];
    fout[16] = 1*f[16];
    fout[17] = 1*f[17];
    fout[18] = 1*f[18];
    fout[19] = 1*f[19];
    fout[20] = 1*f[20];
    fout[21] = 1*f[21];
    fout[22] = 1*f[22];
    fout[23] = 1*f[23];
  }
  if (dir == 3) {
    fout[0] = 1*f[0];
    fout[1] = 1*f[1];
    fout[2] = 1*f[2];
    fout[3] = 1*f[3];
    fout[4] = -1*f[4];
    fout[5] = 1*f[5];
    fout[6] = 1*f[6];
    fout[7] = 1*f[7];
    fout[8] = -1*f[8];
    fout[9] = -1*f[9];
    fout[10] = -1*f[10];
    fout[11] = 1*f[11];
    fout[12] = -1*f[12];
    fout[13] = -1*f[13];
    fout[14] = -1*f[14];
    fout[15] = -1*f[15];
    fout[16] = 1*f[16];
    fout[17] = 1*f[17];
    fout[18] = 1*f[18];
    fout[19] = -1*f[19];
    fout[20] = 1*f[20];
    fout[21] = -1*f[21];
    fout[22] = -1*f[22];
    fout[23] = -1*f[23];
  }
}

GKYL_CU_DH
void
flip_even_sign_2x2v_gkhyb_p1(int dir, const double *f, double *fout )
{
  if (dir == 0) {
    fout[0] = -1*f[0];
    fout[1] = 1*f[1];
    fout[2] = -1*f[2];
    fout[3] = -1*f[3];
    fout[4] = -1*f[4];
    fout[5] = 1*f[5];
    fout[6] = 1*f[6];
    fout[7] = -1*f[7];
    fout[8] = 1*f[8];
    fout[9] = -1*f[9];
    fout[10] = -1*f[10];
    fout[11] = 1*f[11];
    fout[12] = 1*f[12];
    fout[13] = 1*f[13];
    fout[14] = -1*f[14];
    fout[15] = 1*f[15];
    fout[16] = -1*f[16];
    fout[17] = 1*f[17];
    fout[18] = -1*f[18];
    fout[19] = -1*f[19];
    fout[20] = 1*f[20];
    fout[21] = 1*f[21];
    fout[22] = -1*f[22];
    fout[23] = 1*f[23];
  }
  if (dir == 1) {
    fout[0] = -1*f[0];
    fout[1] = -1*f[1];
    fout[2] = 1*f[2];
    fout[3] = -1*f[3];
    fout[4] = -1*f[4];
    fout[5] = 1*f[5];
    fout[6] = -1*f[6];
    fout[7] = 1*f[7];
    fout[8] = -1*f[8];
    fout[9] = 1*f[9];
    fout[10] = -1*f[10];
    fout[11] = 1*f[11];
    fout[12] = 1*f[12];
    fout[13] = -1*f[13];
    fout[14] = 1*f[14];
    fout[15] = 1*f[15];
    fout[16] = -1*f[16];
    fout[17] = -1*f[17];
    fout[18] = 1*f[18];
    fout[19] = -1*f[19];
    fout[20] = 1*f[20];
    fout[21] = -1*f[21];
    fout[22] = 1*f[22];
    fout[23] = 1*f[23];
  }
  if (dir == 2) {
    fout[0] = -1*f[0];
    fout[1] = -1*f[1];
    fout[2] = -1*f[2];
    fout[3] = 1*f[3];
    fout[4] = -1*f[4];
    fout[5] = -1*f[5];
    fout[6] = 1*f[6];
    fout[7] = 1*f[7];
    fout[8] = -1*f[8];
    fout[9] = -1*f[9];
    fout[10] = 1*f[10];
    fout[11] = 1*f[11];
    fout[12] = -1*f[12];
    fout[13] = 1*f[13];
    fout[14] = 1*f[14];
    fout[15] = 1*f[15];
    fout[16] = -1*f[16];
    fout[17] = -1*f[17];
    fout[18] = -1*f[18];
    fout[19] = -1*f[19];
    fout[20] = -1*f[20];
    fout[21] = -1*f[21];
    fout[22] = -1*f[22];
    fout[23] = -1*f[23];
  }
  if (dir == 3) {
    fout[0] = -1*f[0];
    fout[1] = -1*f[1];
    fout[2] = -1*f[2];
    fout[3] = -1*f[3];
    fout[4] = 1*f[4];
    fout[5] = -1*f[5];
    fout[6] = -1*f[6];
    fout[7] = -1*f[7];
    fout[8] = 1*f[8];
    fout[9] = 1*f[9];
    fout[10] = 1*f[10];
    fout[11] = -1*f[11];
    fout[12] = 1*f[12];
    fout[13] = 1*f[13];
    fout[14] = 1*f[14];
    fout[15] = 1*f[15];
    fout[16] = -1*f[16];
    fout[17] = -1*f[17];
    fout[18] = -1*f[18];
    fout[19] = 1*f[19];
    fout[20] = -1*f[20];
    fout[21] = 1*f[21];
    fout[22] = 1*f[22];
    fout[23] = 1*f[23];
  }
}

GKYL_CU_DH
void
flip_odd_sign_3x2v_gkhyb_p1(int dir, const double *f, double *fout )
{
  if (dir == 0) {
    fout[0] = 1*f[0];
    fout[1] = -1*f[1];
    fout[2] = 1*f[2];
    fout[3] = 1*f[3];
    fout[4] = 1*f[4];
    fout[5] = 1*f[5];
    fout[6] = -1*f[6];
    fout[7] = -1*f[7];
    fout[8] = 1*f[8];
    fout[9] = -1*f[9];
    fout[10] = 1*f[10];
    fout[11] = 1*f[11];
    fout[12] = -1*f[12];
    fout[13] = 1*f[13];
    fout[14] = 1*f[14];
    fout[15] = 1*f[15];
    fout[16] = -1*f[16];
    fout[17] = -1*f[17];
    fout[18] = -1*f[18];
    fout[19] = 1*f[19];
    fout[20] = -1*f[20];
    fout[21] = -1*f[21];
    fout[22] = 1*f[22];
    fout[23] = -1*f[23];
    fout[24] = 1*f[24];
    fout[25] = 1*f[25];
    fout[26] = -1*f[26];
    fout[27] = -1*f[27];
    fout[28] = -1*f[28];
    fout[29] = -1*f[29];
    fout[30] = 1*f[30];
    fout[31] = -1*f[31];
    fout[32] = 1*f[32];
    fout[33] = -1*f[33];
    fout[34] = 1*f[34];
    fout[35] = 1*f[35];
    fout[36] = 1*f[36];
    fout[37] = -1*f[37];
    fout[38] = -1*f[38];
    fout[39] = 1*f[39];
    fout[40] = -1*f[40];
    fout[41] = 1*f[41];
    fout[42] = 1*f[42];
    fout[43] = -1*f[43];
    fout[44] = -1*f[44];
    fout[45] = -1*f[45];
    fout[46] = 1*f[46];
    fout[47] = -1*f[47];
  }
  if (dir == 1) {
    fout[0] = 1*f[0];
    fout[1] = 1*f[1];
    fout[2] = -1*f[2];
    fout[3] = 1*f[3];
    fout[4] = 1*f[4];
    fout[5] = 1*f[5];
    fout[6] = -1*f[6];
    fout[7] = 1*f[7];
    fout[8] = -1*f[8];
    fout[9] = 1*f[9];
    fout[10] = -1*f[10];
    fout[11] = 1*f[11];
    fout[12] = 1*f[12];
    fout[13] = -1*f[13];
    fout[14] = 1*f[14];
    fout[15] = 1*f[15];
    fout[16] = -1*f[16];
    fout[17] = -1*f[17];
    fout[18] = 1*f[18];
    fout[19] = -1*f[19];
    fout[20] = -1*f[20];
    fout[21] = 1*f[21];
    fout[22] = -1*f[22];
    fout[23] = 1*f[23];
    fout[24] = -1*f[24];
    fout[25] = 1*f[25];
    fout[26] = -1*f[26];
    fout[27] = -1*f[27];
    fout[28] = -1*f[28];
    fout[29] = 1*f[29];
    fout[30] = -1*f[30];
    fout[31] = -1*f[31];
    fout[32] = 1*f[32];
    fout[33] = 1*f[33];
    fout[34] = -1*f[34];
    fout[35] = 1*f[35];
    fout[36] = 1*f[36];
    fout[37] = -1*f[37];
    fout[38] = 1*f[38];
    fout[39] = -1*f[39];
    fout[40] = 1*f[40];
    fout[41] = -1*f[41];
    fout[42] = 1*f[42];
    fout[43] = -1*f[43];
    fout[44] = -1*f[44];
    fout[45] = 1*f[45];
    fout[46] = -1*f[46];
    fout[47] = -1*f[47];
  }
  if (dir == 2) {
    fout[0] = 1*f[0];
    fout[1] = 1*f[1];
    fout[2] = 1*f[2];
    fout[3] = -1*f[3];
    fout[4] = 1*f[4];
    fout[5] = 1*f[5];
    fout[6] = 1*f[6];
    fout[7] = -1*f[7];
    fout[8] = -1*f[8];
    fout[9] = 1*f[9];
    fout[10] = 1*f[10];
    fout[11] = -1*f[11];
    fout[12] = 1*f[12];
    fout[13] = 1*f[13];
    fout[14] = -1*f[14];
    fout[15] = 1*f[15];
    fout[16] = -1*f[16];
    fout[17] = 1*f[17];
    fout[18] = -1*f[18];
    fout[19] = -1*f[19];
    fout[20] = 1*f[20];
    fout[21] = -1*f[21];
    fout[22] = -1*f[22];
    fout[23] = 1*f[23];
    fout[24] = 1*f[24];
    fout[25] = -1*f[25];
    fout[26] = -1*f[26];
    fout[27] = -1*f[27];
    fout[28] = 1*f[28];
    fout[29] = -1*f[29];
    fout[30] = -1*f[30];
    fout[31] = -1*f[31];
    fout[32] = 1*f[32];
    fout[33] = 1*f[33];
    fout[34] = 1*f[34];
    fout[35] = -1*f[35];
    fout[36] = 1*f[36];
    fout[37] = 1*f[37];
    fout[38] = -1*f[38];
    fout[39] = -1*f[39];
    fout[40] = 1*f[40];
    fout[41] = 1*f[41];
    fout[42] = -1*f[42];
    fout[43] = -1*f[43];
    fout[44] = 1*f[44];
    fout[45] = -1*f[45];
    fout[46] = -1*f[46];
    fout[47] = -1*f[47];
  }
  if (dir == 3) {
    fout[0] = 1*f[0];
    fout[1] = 1*f[1];
    fout[2] = 1*f[2];
    fout[3] = 1*f[3];
    fout[4] = -1*f[4];
    fout[5] = 1*f[5];
    fout[6] = 1*f[6];
    fout[7] = 1*f[7];
    fout[8] = 1*f[8];
    fout[9] = -1*f[9];
    fout[10] = -1*f[10];
    fout[11] = -1*f[11];
    fout[12] = 1*f[12];
    fout[13] = 1*f[13];
    fout[14] = 1*f[14];
    fout[15] = -1*f[15];
    fout[16] = 1*f[16];
    fout[17] = -1*f[17];
    fout[18] = -1*f[18];
    fout[19] = -1*f[19];
    fout[20] = 1*f[20];
    fout[21] = 1*f[21];
    fout[22] = 1*f[22];
    fout[23] = -1*f[23];
    fout[24] = -1*f[24];
    fout[25] = -1*f[25];
    fout[26] = -1*f[26];
    fout[27] = 1*f[27];
    fout[28] = -1*f[28];
    fout[29] = -1*f[29];
    fout[30] = -1*f[30];
    fout[31] = -1*f[31];
    fout[32] = 1*f[32];
    fout[33] = 1*f[33];
    fout[34] = 1*f[34];
    fout[35] = 1*f[35];
    fout[36] = 1*f[36];
    fout[37] = 1*f[37];
    fout[38] = 1*f[38];
    fout[39] = 1*f[39];
    fout[40] = 1*f[40];
    fout[41] = 1*f[41];
    fout[42] = 1*f[42];
    fout[43] = 1*f[43];
    fout[44] = 1*f[44];
    fout[45] = 1*f[45];
    fout[46] = 1*f[46];
    fout[47] = 1*f[47];
  }
  if (dir == 4) {
    fout[0] = 1*f[0];
    fout[1] = 1*f[1];
    fout[2] = 1*f[2];
    fout[3] = 1*f[3];
    fout[4] = 1*f[4];
    fout[5] = -1*f[5];
    fout[6] = 1*f[6];
    fout[7] = 1*f[7];
    fout[8] = 1*f[8];
    fout[9] = 1*f[9];
    fout[10] = 1*f[10];
    fout[11] = 1*f[11];
    fout[12] = -1*f[12];
    fout[13] = -1*f[13];
    fout[14] = -1*f[14];
    fout[15] = -1*f[15];
    fout[16] = 1*f[16];
    fout[17] = 1*f[17];
    fout[18] = 1*f[18];
    fout[19] = 1*f[19];
    fout[20] = -1*f[20];
    fout[21] = -1*f[21];
    fout[22] = -1*f[22];
    fout[23] = -1*f[23];
    fout[24] = -1*f[24];
    fout[25] = -1*f[25];
    fout[26] = 1*f[26];
    fout[27] = -1*f[27];
    fout[28] = -1*f[28];
    fout[29] = -1*f[29];
    fout[30] = -1*f[30];
    fout[31] = -1*f[31];
    fout[32] = 1*f[32];
    fout[33] = 1*f[33];
    fout[34] = 1*f[34];
    fout[35] = 1*f[35];
    fout[36] = -1*f[36];
    fout[37] = 1*f[37];
    fout[38] = 1*f[38];
    fout[39] = 1*f[39];
    fout[40] = -1*f[40];
    fout[41] = -1*f[41];
    fout[42] = -1*f[42];
    fout[43] = 1*f[43];
    fout[44] = -1*f[44];
    fout[45] = -1*f[45];
    fout[46] = -1*f[46];
    fout[47] = -1*f[47];
  }
}

GKYL_CU_DH
void
flip_even_sign_3x2v_gkhyb_p1(int dir, const double *f, double *fout )
{
  if (dir == 0) {
    fout[0] = -1*f[0];
    fout[1] = 1*f[1];
    fout[2] = -1*f[2];
    fout[3] = -1*f[3];
    fout[4] = -1*f[4];
    fout[5] = -1*f[5];
    fout[6] = 1*f[6];
    fout[7] = 1*f[7];
    fout[8] = -1*f[8];
    fout[9] = 1*f[9];
    fout[10] = -1*f[10];
    fout[11] = -1*f[11];
    fout[12] = 1*f[12];
    fout[13] = -1*f[13];
    fout[14] = -1*f[14];
    fout[15] = -1*f[15];
    fout[16] = 1*f[16];
    fout[17] = 1*f[17];
    fout[18] = 1*f[18];
    fout[19] = -1*f[19];
    fout[20] = 1*f[20];
    fout[21] = 1*f[21];
    fout[22] = -1*f[22];
    fout[23] = 1*f[23];
    fout[24] = -1*f[24];
    fout[25] = -1*f[25];
    fout[26] = 1*f[26];
    fout[27] = 1*f[27];
    fout[28] = 1*f[28];
    fout[29] = 1*f[29];
    fout[30] = -1*f[30];
    fout[31] = 1*f[31];
    fout[32] = -1*f[32];
    fout[33] = 1*f[33];
    fout[34] = -1*f[34];
    fout[35] = -1*f[35];
    fout[36] = -1*f[36];
    fout[37] = 1*f[37];
    fout[38] = 1*f[38];
    fout[39] = -1*f[39];
    fout[40] = 1*f[40];
    fout[41] = -1*f[41];
    fout[42] = -1*f[42];
    fout[43] = 1*f[43];
    fout[44] = 1*f[44];
    fout[45] = 1*f[45];
    fout[46] = -1*f[46];
    fout[47] = 1*f[47];
  }
  if (dir == 1) {
    fout[0] = -1*f[0];
    fout[1] = -1*f[1];
    fout[2] = 1*f[2];
    fout[3] = -1*f[3];
    fout[4] = -1*f[4];
    fout[5] = -1*f[5];
    fout[6] = 1*f[6];
    fout[7] = -1*f[7];
    fout[8] = 1*f[8];
    fout[9] = -1*f[9];
    fout[10] = 1*f[10];
    fout[11] = -1*f[11];
    fout[12] = -1*f[12];
    fout[13] = 1*f[13];
    fout[14] = -1*f[14];
    fout[15] = -1*f[15];
    fout[16] = 1*f[16];
    fout[17] = 1*f[17];
    fout[18] = -1*f[18];
    fout[19] = 1*f[19];
    fout[20] = 1*f[20];
    fout[21] = -1*f[21];
    fout[22] = 1*f[22];
    fout[23] = -1*f[23];
    fout[24] = 1*f[24];
    fout[25] = -1*f[25];
    fout[26] = 1*f[26];
    fout[27] = 1*f[27];
    fout[28] = 1*f[28];
    fout[29] = -1*f[29];
    fout[30] = 1*f[30];
    fout[31] = 1*f[31];
    fout[32] = -1*f[32];
    fout[33] = -1*f[33];
    fout[34] = 1*f[34];
    fout[35] = -1*f[35];
    fout[36] = -1*f[36];
    fout[37] = 1*f[37];
    fout[38] = -1*f[38];
    fout[39] = 1*f[39];
    fout[40] = -1*f[40];
    fout[41] = 1*f[41];
    fout[42] = -1*f[42];
    fout[43] = 1*f[43];
    fout[44] = 1*f[44];
    fout[45] = -1*f[45];
    fout[46] = 1*f[46];
    fout[47] = 1*f[47];
  }
  if (dir == 2) {
    fout[0] = -1*f[0];
    fout[1] = -1*f[1];
    fout[2] = -1*f[2];
    fout[3] = 1*f[3];
    fout[4] = -1*f[4];
    fout[5] = -1*f[5];
    fout[6] = -1*f[6];
    fout[7] = 1*f[7];
    fout[8] = 1*f[8];
    fout[9] = -1*f[9];
    fout[10] = -1*f[10];
    fout[11] = 1*f[11];
    fout[12] = -1*f[12];
    fout[13] = -1*f[13];
    fout[14] = 1*f[14];
    fout[15] = -1*f[15];
    fout[16] = 1*f[16];
    fout[17] = -1*f[17];
    fout[18] = 1*f[18];
    fout[19] = 1*f[19];
    fout[20] = -1*f[20];
    fout[21] = 1*f[21];
    fout[22] = 1*f[22];
    fout[23] = -1*f[23];
    fout[24] = -1*f[24];
    fout[25] = 1*f[25];
    fout[26] = 1*f[26];
    fout[27] = 1*f[27];
    fout[28] = -1*f[28];
    fout[29] = 1*f[29];
    fout[30] = 1*f[30];
    fout[31] = 1*f[31];
    fout[32] = -1*f[32];
    fout[33] = -1*f[33];
    fout[34] = -1*f[34];
    fout[35] = 1*f[35];
    fout[36] = -1*f[36];
    fout[37] = -1*f[37];
    fout[38] = 1*f[38];
    fout[39] = 1*f[39];
    fout[40] = -1*f[40];
    fout[41] = -1*f[41];
    fout[42] = 1*f[42];
    fout[43] = 1*f[43];
    fout[44] = -1*f[44];
    fout[45] = 1*f[45];
    fout[46] = 1*f[46];
    fout[47] = 1*f[47];
  }
  if (dir == 3) {
    fout[0] = -1*f[0];
    fout[1] = -1*f[1];
    fout[2] = -1*f[2];
    fout[3] = -1*f[3];
    fout[4] = 1*f[4];
    fout[5] = -1*f[5];
    fout[6] = -1*f[6];
    fout[7] = -1*f[7];
    fout[8] = -1*f[8];
    fout[9] = 1*f[9];
    fout[10] = 1*f[10];
    fout[11] = 1*f[11];
    fout[12] = -1*f[12];
    fout[13] = -1*f[13];
    fout[14] = -1*f[14];
    fout[15] = 1*f[15];
    fout[16] = -1*f[16];
    fout[17] = 1*f[17];
    fout[18] = 1*f[18];
    fout[19] = 1*f[19];
    fout[20] = -1*f[20];
    fout[21] = -1*f[21];
    fout[22] = -1*f[22];
    fout[23] = 1*f[23];
    fout[24] = 1*f[24];
    fout[25] = 1*f[25];
    fout[26] = 1*f[26];
    fout[27] = -1*f[27];
    fout[28] = 1*f[28];
    fout[29] = 1*f[29];
    fout[30] = 1*f[30];
    fout[31] = 1*f[31];
    fout[32] = -1*f[32];
    fout[33] = -1*f[33];
    fout[34] = -1*f[34];
    fout[35] = -1*f[35];
    fout[36] = -1*f[36];
    fout[37] = -1*f[37];
    fout[38] = -1*f[38];
    fout[39] = -1*f[39];
    fout[40] = -1*f[40];
    fout[41] = -1*f[41];
    fout[42] = -1*f[42];
    fout[43] = -1*f[43];
    fout[44] = -1*f[44];
    fout[45] = -1*f[45];
    fout[46] = -1*f[46];
    fout[47] = -1*f[47];
  }
  if (dir == 4) {
    fout[0] = -1*f[0];
    fout[1] = -1*f[1];
    fout[2] = -1*f[2];
    fout[3] = -1*f[3];
    fout[4] = -1*f[4];
    fout[5] = 1*f[5];
    fout[6] = -1*f[6];
    fout[7] = -1*f[7];
    fout[8] = -1*f[8];
    fout[9] = -1*f[9];
    fout[10] = -1*f[10];
    fout[11] = -1*f[11];
    fout[12] = 1*f[12];
    fout[13] = 1*f[13];
    fout[14] = 1*f[14];
    fout[15] = 1*f[15];
    fout[16] = -1*f[16];
    fout[17] = -1*f[17];
    fout[18] = -1*f[18];
    fout[19] = -1*f[19];
    fout[20] = 1*f[20];
    fout[21] = 1*f[21];
    fout[22] = 1*f[22];
    fout[23] = 1*f[23];
    fout[24] = 1*f[24];
    fout[25] = 1*f[25];
    fout[26] = -1*f[26];
    fout[27] = 1*f[27];
    fout[28] = 1*f[28];
    fout[29] = 1*f[29];
    fout[30] = 1*f[30];
    fout[31] = 1*f[31];
    fout[32] = -1*f[32];
    fout[33] = -1*f[33];
    fout[34] = -1*f[34];
    fout[35] = -1*f[35];
    fout[36] = 1*f[36];
    fout[37] = -1*f[37];
    fout[38] = -1*f[38];
    fout[39] = -1*f[39];
    fout[40] = 1*f[40];
    fout[41] = 1*f[41];
    fout[42] = 1*f[42];
    fout[43] = -1*f[43];
    fout[44] = 1*f[44];
    fout[45] = 1*f[45];
    fout[46] = 1*f[46];
    fout[47] = 1*f[47];
  }
}

