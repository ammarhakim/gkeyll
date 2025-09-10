#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_ten_moment_grad_closure.h>
#include <gkyl_moment_non_ideal_priv.h>

// Makes indexing cleaner
static const unsigned T11 = 0;
static const unsigned T12 = 1;
static const unsigned T13 = 2;
static const unsigned T22 = 3;
static const unsigned T23 = 4;
static const unsigned T33 = 5;

static const unsigned Q111 = 0;
static const unsigned Q112 = 1;
static const unsigned Q113 = 2;
static const unsigned Q122 = 3;
static const unsigned Q123 = 4;
static const unsigned Q133 = 5;
static const unsigned Q222 = 6;
static const unsigned Q223 = 7;
static const unsigned Q233 = 8;
static const unsigned Q333 = 9;

// 1D stencil locations (L: lower, U: upper)
enum loc_1d {
  L_1D, U_1D
};

// 2D stencil locations (L: lower, U: upper)
enum loc_2d {
  LL_2D, LU_2D,
  UL_2D, UU_2D
};

// 3D stencil locations (L: lower, U: upper)
enum loc_3d {
  LLL_3D, LLU_3D,
  LUL_3D, LUU_3D,
  ULL_3D, ULU_3D,
  UUL_3D, UUU_3D
};

struct gkyl_ten_moment_grad_closure {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  double k0; // damping coefficient
};

static void
create_offsets_vertices(const struct gkyl_range *range, long offsets[])
{
  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, (int[]) { -1, -1, -1 }, (int[]) { 0, 0, 0 });

  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);

  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}

static void
create_offsets_centers(const struct gkyl_range *range, long offsets[])
{
  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, (int[]) { 0, 0, 0 }, (int[]) { 1, 1, 1 });

  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);

  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}

static void
var_setup(const gkyl_ten_moment_grad_closure *gces,
  int start, int end,
  const double *fluid_d[],
  double rho[], double p[], double Tij[][6])
{
  for (int j = start; j <= end; ++j) {
    rho[j] = fluid_d[j][RHO];
    p[j] = (fluid_d[j][P11] - fluid_d[j][MX] * fluid_d[j][MX] / fluid_d[j][RHO]
          + fluid_d[j][P22] - fluid_d[j][MY] * fluid_d[j][MY] / fluid_d[j][RHO]
          + fluid_d[j][P33] - fluid_d[j][MZ] * fluid_d[j][MZ] / fluid_d[j][RHO])/3.0;
    Tij[j][T11] = (fluid_d[j][P11] - fluid_d[j][MX] * fluid_d[j][MX] / fluid_d[j][RHO]) / fluid_d[j][RHO];
    Tij[j][T12] = (fluid_d[j][P12] - fluid_d[j][MX] * fluid_d[j][MY] / fluid_d[j][RHO]) / fluid_d[j][RHO];
    Tij[j][T13] = (fluid_d[j][P13] - fluid_d[j][MX] * fluid_d[j][MZ] / fluid_d[j][RHO]) / fluid_d[j][RHO];
    Tij[j][T22] = (fluid_d[j][P22] - fluid_d[j][MY] * fluid_d[j][MY] / fluid_d[j][RHO]) / fluid_d[j][RHO];
    Tij[j][T23] = (fluid_d[j][P23] - fluid_d[j][MY] * fluid_d[j][MZ] / fluid_d[j][RHO]) / fluid_d[j][RHO];
    Tij[j][T33] = (fluid_d[j][P33] - fluid_d[j][MZ] * fluid_d[j][MZ] / fluid_d[j][RHO]) / fluid_d[j][RHO];
  }
}

static void
calc_unmag_heat_flux(const gkyl_ten_moment_grad_closure *gces,
  const double *fluid_d[], double *cflrate, double *heat_flux_d)
{
  const int ndim = gces->ndim;
  double rho_avg = 0.0;
  double p_avg = 0.0;
  double dTdx[6] = {0.0};
  double dTdy[6] = {0.0};
  double dTdz[6] = {0.0};

  if (ndim == 1) {
    const double dx = gces->grid.dx[0];
    double Tij[2][6] = {0.0};
    double rho[2] = {0.0};
    double p[2] = {0.0};
    var_setup(gces, L_1D, U_1D, fluid_d, rho, p, Tij);

    rho_avg = calc_harmonic_avg_1D(rho[L_1D], rho[U_1D]);
    p_avg = calc_harmonic_avg_1D(p[L_1D], p[U_1D]);

    dTdx[T11] = calc_sym_grad_1D(dx, Tij[L_1D][T11], Tij[U_1D][T11]);
    dTdx[T12] = calc_sym_grad_1D(dx, Tij[L_1D][T12], Tij[U_1D][T12]);
    dTdx[T13] = calc_sym_grad_1D(dx, Tij[L_1D][T13], Tij[U_1D][T13]);
    dTdx[T22] = calc_sym_grad_1D(dx, Tij[L_1D][T22], Tij[U_1D][T22]);
    dTdx[T23] = calc_sym_grad_1D(dx, Tij[L_1D][T23], Tij[U_1D][T23]);
    dTdx[T33] = calc_sym_grad_1D(dx, Tij[L_1D][T33], Tij[U_1D][T33]);
  }
  else if (ndim == 2) {
    const double dx = gces->grid.dx[0];
    const double dy = gces->grid.dx[1];
    double Tij[4][6] = {0.0};
    double rho[4] = {0.0};
    double p[4] = {0.0};
    var_setup(gces, LL_2D, UU_2D, fluid_d, rho, p, Tij);

    rho_avg = calc_harmonic_avg_2D(rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
    p_avg = calc_harmonic_avg_2D(p[LL_2D], p[LU_2D], p[UL_2D], p[UU_2D]);

    dTdx[T11] = calc_sym_gradx_2D(dx, Tij[LL_2D][T11], Tij[LU_2D][T11], Tij[UL_2D][T11], Tij[UU_2D][T11]);
    dTdx[T12] = calc_sym_gradx_2D(dx, Tij[LL_2D][T12], Tij[LU_2D][T12], Tij[UL_2D][T12], Tij[UU_2D][T12]);
    dTdx[T13] = calc_sym_gradx_2D(dx, Tij[LL_2D][T13], Tij[LU_2D][T13], Tij[UL_2D][T13], Tij[UU_2D][T13]);
    dTdx[T22] = calc_sym_gradx_2D(dx, Tij[LL_2D][T22], Tij[LU_2D][T22], Tij[UL_2D][T22], Tij[UU_2D][T22]);
    dTdx[T23] = calc_sym_gradx_2D(dx, Tij[LL_2D][T23], Tij[LU_2D][T23], Tij[UL_2D][T23], Tij[UU_2D][T23]);
    dTdx[T33] = calc_sym_gradx_2D(dx, Tij[LL_2D][T33], Tij[LU_2D][T33], Tij[UL_2D][T33], Tij[UU_2D][T33]);

    dTdy[T11] = calc_sym_grady_2D(dy, Tij[LL_2D][T11], Tij[LU_2D][T11], Tij[UL_2D][T11], Tij[UU_2D][T11]);
    dTdy[T12] = calc_sym_grady_2D(dy, Tij[LL_2D][T12], Tij[LU_2D][T12], Tij[UL_2D][T12], Tij[UU_2D][T12]);
    dTdy[T13] = calc_sym_grady_2D(dy, Tij[LL_2D][T13], Tij[LU_2D][T13], Tij[UL_2D][T13], Tij[UU_2D][T13]);
    dTdy[T22] = calc_sym_grady_2D(dy, Tij[LL_2D][T22], Tij[LU_2D][T22], Tij[UL_2D][T22], Tij[UU_2D][T22]);
    dTdy[T23] = calc_sym_grady_2D(dy, Tij[LL_2D][T23], Tij[LU_2D][T23], Tij[UL_2D][T23], Tij[UU_2D][T23]);
    dTdy[T33] = calc_sym_grady_2D(dy, Tij[LL_2D][T33], Tij[LU_2D][T33], Tij[UL_2D][T33], Tij[UU_2D][T33]);
  }
  else if (ndim == 3) {
    const double dx = gces->grid.dx[0];
    const double dy = gces->grid.dx[1];
    const double dz = gces->grid.dx[2];

    double Tij[8][6] = {0.0};
    double rho[8] = {0.0};
    double p[8] = {0.0};
    var_setup(gces, LLL_3D, UUU_3D, fluid_d, rho, p, Tij);

    rho_avg = calc_harmonic_avg_3D(rho[LLL_3D], rho[LLU_3D], rho[LUL_3D], rho[LUU_3D],
                                   rho[ULL_3D], rho[ULU_3D], rho[UUL_3D], rho[UUU_3D]);
    p_avg = calc_harmonic_avg_3D(p[LLL_3D], p[LLU_3D], p[LUL_3D], p[LUU_3D],
                                 p[ULL_3D], p[ULU_3D], p[UUL_3D], p[UUU_3D]);

    dTdx[T11] = calc_sym_gradx_3D(dx, Tij[LLL_3D][T11], Tij[LLU_3D][T11], Tij[LUL_3D][T11], Tij[LUU_3D][T11],
                                      Tij[ULL_3D][T11], Tij[ULU_3D][T11], Tij[UUL_3D][T11], Tij[UUU_3D][T11]);
    dTdx[T12] = calc_sym_gradx_3D(dx, Tij[LLL_3D][T12], Tij[LLU_3D][T12], Tij[LUL_3D][T12], Tij[LUU_3D][T12],
                                      Tij[ULL_3D][T12], Tij[ULU_3D][T12], Tij[UUL_3D][T12], Tij[UUU_3D][T12]);
    dTdx[T13] = calc_sym_gradx_3D(dx, Tij[LLL_3D][T13], Tij[LLU_3D][T13], Tij[LUL_3D][T13], Tij[LUU_3D][T13],
                                      Tij[ULL_3D][T13], Tij[ULU_3D][T13], Tij[UUL_3D][T13], Tij[UUU_3D][T13]);
    dTdx[T22] = calc_sym_gradx_3D(dx, Tij[LLL_3D][T22], Tij[LLU_3D][T22], Tij[LUL_3D][T22], Tij[LUU_3D][T22],
                                      Tij[ULL_3D][T22], Tij[ULU_3D][T22], Tij[UUL_3D][T22], Tij[UUU_3D][T22]);
    dTdx[T23] = calc_sym_gradx_3D(dx, Tij[LLL_3D][T23], Tij[LLU_3D][T23], Tij[LUL_3D][T23], Tij[LUU_3D][T23],
                                      Tij[ULL_3D][T23], Tij[ULU_3D][T23], Tij[UUL_3D][T23], Tij[UUU_3D][T23]);
    dTdx[T33] = calc_sym_gradx_3D(dx, Tij[LLL_3D][T33], Tij[LLU_3D][T33], Tij[LUL_3D][T33], Tij[LUU_3D][T33],
                                      Tij[ULL_3D][T33], Tij[ULU_3D][T33], Tij[UUL_3D][T33], Tij[UUU_3D][T33]);

    dTdy[T11] = calc_sym_grady_3D(dy, Tij[LLL_3D][T11], Tij[LLU_3D][T11], Tij[LUL_3D][T11], Tij[LUU_3D][T11],
                                      Tij[ULL_3D][T11], Tij[ULU_3D][T11], Tij[UUL_3D][T11], Tij[UUU_3D][T11]);
    dTdy[T12] = calc_sym_grady_3D(dy, Tij[LLL_3D][T12], Tij[LLU_3D][T12], Tij[LUL_3D][T12], Tij[LUU_3D][T12],
                                      Tij[ULL_3D][T12], Tij[ULU_3D][T12], Tij[UUL_3D][T12], Tij[UUU_3D][T12]);
    dTdy[T13] = calc_sym_grady_3D(dy, Tij[LLL_3D][T13], Tij[LLU_3D][T13], Tij[LUL_3D][T13], Tij[LUU_3D][T13],
                                      Tij[ULL_3D][T13], Tij[ULU_3D][T13], Tij[UUL_3D][T13], Tij[UUU_3D][T13]);
    dTdy[T22] = calc_sym_grady_3D(dy, Tij[LLL_3D][T22], Tij[LLU_3D][T22], Tij[LUL_3D][T22], Tij[LUU_3D][T22],
                                      Tij[ULL_3D][T22], Tij[ULU_3D][T22], Tij[UUL_3D][T22], Tij[UUU_3D][T22]);
    dTdy[T23] = calc_sym_grady_3D(dy, Tij[LLL_3D][T23], Tij[LLU_3D][T23], Tij[LUL_3D][T23], Tij[LUU_3D][T23],
                                      Tij[ULL_3D][T23], Tij[ULU_3D][T23], Tij[UUL_3D][T23], Tij[UUU_3D][T23]);
    dTdy[T33] = calc_sym_grady_3D(dy, Tij[LLL_3D][T33], Tij[LLU_3D][T33], Tij[LUL_3D][T33], Tij[LUU_3D][T33],
                                      Tij[ULL_3D][T33], Tij[ULU_3D][T33], Tij[UUL_3D][T33], Tij[UUU_3D][T33]);

    dTdz[T11] = calc_sym_gradz_3D(dz, Tij[LLL_3D][T11], Tij[LLU_3D][T11], Tij[LUL_3D][T11], Tij[LUU_3D][T11],
                                      Tij[ULL_3D][T11], Tij[ULU_3D][T11], Tij[UUL_3D][T11], Tij[UUU_3D][T11]);
    dTdz[T12] = calc_sym_gradz_3D(dz, Tij[LLL_3D][T12], Tij[LLU_3D][T12], Tij[LUL_3D][T12], Tij[LUU_3D][T12],
                                      Tij[ULL_3D][T12], Tij[ULU_3D][T12], Tij[UUL_3D][T12], Tij[UUU_3D][T12]);
    dTdz[T13] = calc_sym_gradz_3D(dz, Tij[LLL_3D][T13], Tij[LLU_3D][T13], Tij[LUL_3D][T13], Tij[LUU_3D][T13],
                                      Tij[ULL_3D][T13], Tij[ULU_3D][T13], Tij[UUL_3D][T13], Tij[UUU_3D][T13]);
    dTdz[T22] = calc_sym_gradz_3D(dz, Tij[LLL_3D][T22], Tij[LLU_3D][T22], Tij[LUL_3D][T22], Tij[LUU_3D][T22],
                                      Tij[ULL_3D][T22], Tij[ULU_3D][T22], Tij[UUL_3D][T22], Tij[UUU_3D][T22]);
    dTdz[T23] = calc_sym_gradz_3D(dz, Tij[LLL_3D][T23], Tij[LLU_3D][T23], Tij[LUL_3D][T23], Tij[LUU_3D][T23],
                                      Tij[ULL_3D][T23], Tij[ULU_3D][T23], Tij[UUL_3D][T23], Tij[UUU_3D][T23]);
    dTdz[T33] = calc_sym_gradz_3D(dz, Tij[LLL_3D][T33], Tij[LLU_3D][T33], Tij[LUL_3D][T33], Tij[LUU_3D][T33],
                                      Tij[ULL_3D][T33], Tij[ULU_3D][T33], Tij[UUL_3D][T33], Tij[UUU_3D][T33]);
  }
  double alpha = 1.0/gces->k0;
  double vth_avg = sqrt(p_avg/rho_avg);
  heat_flux_d[Q111] = alpha*vth_avg*rho_avg*(dTdx[T11] + dTdx[T11] + dTdx[T11])/3.0;
  heat_flux_d[Q112] = alpha*vth_avg*rho_avg*(dTdx[T12] + dTdx[T12] + dTdy[T11])/3.0;
  heat_flux_d[Q113] = alpha*vth_avg*rho_avg*(dTdx[T13] + dTdx[T13] + dTdz[T11])/3.0;
  heat_flux_d[Q122] = alpha*vth_avg*rho_avg*(dTdx[T22] + dTdy[T12] + dTdy[T12])/3.0;
  heat_flux_d[Q123] = alpha*vth_avg*rho_avg*(dTdx[T23] + dTdy[T13] + dTdz[T12])/3.0;
  heat_flux_d[Q133] = alpha*vth_avg*rho_avg*(dTdx[T33] + dTdz[T13] + dTdz[T13])/3.0;
  heat_flux_d[Q222] = alpha*vth_avg*rho_avg*(dTdy[T22] + dTdy[T22] + dTdy[T22])/3.0;
  heat_flux_d[Q223] = alpha*vth_avg*rho_avg*(dTdy[T23] + dTdy[T23] + dTdz[T22])/3.0;
  heat_flux_d[Q233] = alpha*vth_avg*rho_avg*(dTdy[T33] + dTdz[T23] + dTdz[T23])/3.0;
  heat_flux_d[Q333] = alpha*vth_avg*rho_avg*(dTdz[T33] + dTdz[T33] + dTdz[T33])/3.0;
}

static void
calc_grad_closure_update(const gkyl_ten_moment_grad_closure *gces,
  const double *heat_flux_d[], double *rhs)
{
  const int ndim = gces->ndim;
  double divQx[6] = {0.0};
  double divQy[6] = {0.0};
  double divQz[6] = {0.0};

  if (ndim == 1) {
    const double dx = gces->grid.dx[0];
    double q[2][10] = {0.0};
    for (int j = L_1D; j <= U_1D; ++j)
      for (int k = 0; k < 10; ++k)
        q[j][k] = heat_flux_d[j][k];

    divQx[0] = calc_sym_grad_1D(dx, q[L_1D][Q111], q[U_1D][Q111]);
    divQx[1] = calc_sym_grad_1D(dx, q[L_1D][Q112], q[U_1D][Q112]);
    divQx[2] = calc_sym_grad_1D(dx, q[L_1D][Q113], q[U_1D][Q113]);
    divQx[3] = calc_sym_grad_1D(dx, q[L_1D][Q122], q[U_1D][Q122]);
    divQx[4] = calc_sym_grad_1D(dx, q[L_1D][Q123], q[U_1D][Q123]);
    divQx[5] = calc_sym_grad_1D(dx, q[L_1D][Q133], q[U_1D][Q133]);
  }
  else if (ndim == 2) {
    const double dx = gces->grid.dx[0];
    const double dy = gces->grid.dx[1];
    double q[4][10] = {0.0};
    for (int j = LL_2D; j <= UU_2D; ++j)
      for (int k = 0; k < 10; ++k)
        q[j][k] = heat_flux_d[j][k];

    divQx[0] = calc_sym_gradx_2D(dx, q[LL_2D][Q111], q[LU_2D][Q111], q[UL_2D][Q111], q[UU_2D][Q111]);
    divQx[1] = calc_sym_gradx_2D(dx, q[LL_2D][Q112], q[LU_2D][Q112], q[UL_2D][Q112], q[UU_2D][Q112]);
    divQx[2] = calc_sym_gradx_2D(dx, q[LL_2D][Q113], q[LU_2D][Q113], q[UL_2D][Q113], q[UU_2D][Q113]);
    divQx[3] = calc_sym_gradx_2D(dx, q[LL_2D][Q122], q[LU_2D][Q122], q[UL_2D][Q122], q[UU_2D][Q122]);
    divQx[4] = calc_sym_gradx_2D(dx, q[LL_2D][Q123], q[LU_2D][Q123], q[UL_2D][Q123], q[UU_2D][Q123]);
    divQx[5] = calc_sym_gradx_2D(dx, q[LL_2D][Q133], q[LU_2D][Q133], q[UL_2D][Q133], q[UU_2D][Q133]);

    divQy[0] = calc_sym_grady_2D(dy, q[LL_2D][Q112], q[LU_2D][Q112], q[UL_2D][Q112], q[UU_2D][Q112]);
    divQy[1] = calc_sym_grady_2D(dy, q[LL_2D][Q122], q[LU_2D][Q122], q[UL_2D][Q122], q[UU_2D][Q122]);
    divQy[2] = calc_sym_grady_2D(dy, q[LL_2D][Q123], q[LU_2D][Q123], q[UL_2D][Q123], q[UU_2D][Q123]);
    divQy[3] = calc_sym_grady_2D(dy, q[LL_2D][Q222], q[LU_2D][Q222], q[UL_2D][Q222], q[UU_2D][Q222]);
    divQy[4] = calc_sym_grady_2D(dy, q[LL_2D][Q223], q[LU_2D][Q223], q[UL_2D][Q223], q[UU_2D][Q223]);
    divQy[5] = calc_sym_grady_2D(dy, q[LL_2D][Q233], q[LU_2D][Q233], q[UL_2D][Q233], q[UU_2D][Q233]);
  }
  else if (ndim == 3) {
    const double dx = gces->grid.dx[0];
    const double dy = gces->grid.dx[1];
    const double dz = gces->grid.dx[2];
    double q[8][10] = {0.0};
    for (int j = LLL_3D; j <= UUU_3D; ++j)
      for (int k = 0; k < 10; ++k)
        q[j][k] = heat_flux_d[j][k];

    divQx[0] = calc_sym_gradx_3D(dx, q[LLL_3D][Q111], q[LLU_3D][Q111], q[LUL_3D][Q111], q[LUU_3D][Q111],
                                     q[ULL_3D][Q111], q[ULU_3D][Q111], q[UUL_3D][Q111], q[UUU_3D][Q111]);
    divQx[1] = calc_sym_gradx_3D(dx, q[LLL_3D][Q112], q[LLU_3D][Q112], q[LUL_3D][Q112], q[LUU_3D][Q112],
                                     q[ULL_3D][Q112], q[ULU_3D][Q112], q[UUL_3D][Q112], q[UUU_3D][Q112]);
    divQx[2] = calc_sym_gradx_3D(dx, q[LLL_3D][Q113], q[LLU_3D][Q113], q[LUL_3D][Q113], q[LUU_3D][Q113],
                                     q[ULL_3D][Q113], q[ULU_3D][Q113], q[UUL_3D][Q113], q[UUU_3D][Q113]);
    divQx[3] = calc_sym_gradx_3D(dx, q[LLL_3D][Q122], q[LLU_3D][Q122], q[LUL_3D][Q122], q[LUU_3D][Q122],
                                     q[ULL_3D][Q122], q[ULU_3D][Q122], q[UUL_3D][Q122], q[UUU_3D][Q122]);
    divQx[4] = calc_sym_gradx_3D(dx, q[LLL_3D][Q123], q[LLU_3D][Q123], q[LUL_3D][Q123], q[LUU_3D][Q123],
                                     q[ULL_3D][Q123], q[ULU_3D][Q123], q[UUL_3D][Q123], q[UUU_3D][Q123]);
    divQx[5] = calc_sym_gradx_3D(dx, q[LLL_3D][Q133], q[LLU_3D][Q133], q[LUL_3D][Q133], q[LUU_3D][Q133],
                                     q[ULL_3D][Q133], q[ULU_3D][Q133], q[UUL_3D][Q133], q[UUU_3D][Q133]);

    divQy[0] = calc_sym_grady_3D(dy, q[LLL_3D][Q112], q[LLU_3D][Q112], q[LUL_3D][Q112], q[LUU_3D][Q112],
                                     q[ULL_3D][Q112], q[ULU_3D][Q112], q[UUL_3D][Q112], q[UUU_3D][Q112]);
    divQy[1] = calc_sym_grady_3D(dy, q[LLL_3D][Q122], q[LLU_3D][Q122], q[LUL_3D][Q122], q[LUU_3D][Q122],
                                     q[ULL_3D][Q122], q[ULU_3D][Q122], q[UUL_3D][Q122], q[UUU_3D][Q122]);
    divQy[2] = calc_sym_grady_3D(dy, q[LLL_3D][Q123], q[LLU_3D][Q123], q[LUL_3D][Q123], q[LUU_3D][Q123],
                                     q[ULL_3D][Q123], q[ULU_3D][Q123], q[UUL_3D][Q123], q[UUU_3D][Q123]);
    divQy[3] = calc_sym_grady_3D(dy, q[LLL_3D][Q222], q[LLU_3D][Q222], q[LUL_3D][Q222], q[LUU_3D][Q222],
                                     q[ULL_3D][Q222], q[ULU_3D][Q222], q[UUL_3D][Q222], q[UUU_3D][Q222]);
    divQy[4] = calc_sym_grady_3D(dy, q[LLL_3D][Q223], q[LLU_3D][Q223], q[LUL_3D][Q223], q[LUU_3D][Q223],
                                     q[ULL_3D][Q223], q[ULU_3D][Q223], q[UUL_3D][Q223], q[UUU_3D][Q223]);
    divQy[5] = calc_sym_grady_3D(dy, q[LLL_3D][Q233], q[LLU_3D][Q233], q[LUL_3D][Q233], q[LUU_3D][Q233],
                                     q[ULL_3D][Q233], q[ULU_3D][Q233], q[UUL_3D][Q233], q[UUU_3D][Q233]);

    divQz[0] = calc_sym_gradz_3D(dz, q[LLL_3D][Q113], q[LLU_3D][Q113], q[LUL_3D][Q113], q[LUU_3D][Q113],
                                     q[ULL_3D][Q113], q[ULU_3D][Q113], q[UUL_3D][Q113], q[UUU_3D][Q113]);
    divQz[1] = calc_sym_gradz_3D(dz, q[LLL_3D][Q123], q[LLU_3D][Q123], q[LUL_3D][Q123], q[LUU_3D][Q123],
                                     q[ULL_3D][Q123], q[ULU_3D][Q123], q[UUL_3D][Q123], q[UUU_3D][Q123]);
    divQz[2] = calc_sym_gradz_3D(dz, q[LLL_3D][Q133], q[LLU_3D][Q133], q[LUL_3D][Q133], q[LUU_3D][Q133],
                                     q[ULL_3D][Q133], q[ULU_3D][Q133], q[UUL_3D][Q133], q[UUU_3D][Q133]);
    divQz[3] = calc_sym_gradz_3D(dz, q[LLL_3D][Q223], q[LLU_3D][Q223], q[LUL_3D][Q223], q[LUU_3D][Q223],
                                     q[ULL_3D][Q223], q[ULU_3D][Q223], q[UUL_3D][Q223], q[UUU_3D][Q223]);
    divQz[4] = calc_sym_gradz_3D(dz, q[LLL_3D][Q233], q[LLU_3D][Q233], q[LUL_3D][Q233], q[LUU_3D][Q233],
                                     q[ULL_3D][Q233], q[ULU_3D][Q233], q[UUL_3D][Q233], q[UUU_3D][Q233]);
    divQz[5] = calc_sym_gradz_3D(dz, q[LLL_3D][Q333], q[LLU_3D][Q333], q[LUL_3D][Q333], q[LUU_3D][Q333],
                                     q[ULL_3D][Q333], q[ULU_3D][Q333], q[UUL_3D][Q333], q[UUU_3D][Q333]);
  }
  rhs[RHO] = 0.0;
  rhs[MX] = 0.0;
  rhs[MY] = 0.0;
  rhs[MZ] = 0.0;
  rhs[P11] = divQx[0] + divQy[0] + divQz[0];
  rhs[P12] = divQx[1] + divQy[1] + divQz[1];
  rhs[P13] = divQx[2] + divQy[2] + divQz[2];
  rhs[P22] = divQx[3] + divQy[3] + divQz[3];
  rhs[P23] = divQx[4] + divQy[4] + divQz[4];
  rhs[P33] = divQx[5] + divQy[5] + divQz[5];
}

gkyl_ten_moment_grad_closure*
gkyl_ten_moment_grad_closure_new(struct gkyl_ten_moment_grad_closure_inp inp)
{
  gkyl_ten_moment_grad_closure *up = gkyl_malloc(sizeof(gkyl_ten_moment_grad_closure));

  up->grid = *(inp.grid);
  up->ndim = up->grid.ndim;
  up->k0 = inp.k0;

  return up;
}

void
gkyl_ten_moment_grad_closure_advance(const gkyl_ten_moment_grad_closure *gces,
  const struct gkyl_range *heat_flux_range, const struct gkyl_range *update_range,
  const struct gkyl_array *fluid, const struct gkyl_array *em_tot,
  struct gkyl_array *cflrate, struct gkyl_array *heat_flux,
  struct gkyl_array *rhs)
{
  int ndim = update_range->ndim;
  long sz[] = { 2, 4, 8 };

  long offsets_vertices[sz[ndim-1]];
  create_offsets_vertices(update_range, offsets_vertices);

  long offsets_centers[sz[ndim-1]];
  create_offsets_centers(heat_flux_range, offsets_centers);

  const double* fluid_d[sz[ndim-1]];
  const double* em_tot_d[sz[ndim-1]];
  double *heat_flux_d;
  const double* heat_flux_up[sz[ndim-1]];
  double *rhs_d;

  struct gkyl_range_iter iter_vertex;
  gkyl_range_iter_init(&iter_vertex, heat_flux_range);
  while (gkyl_range_iter_next(&iter_vertex)) {

    long linc_vertex = gkyl_range_idx(heat_flux_range, iter_vertex.idx);
    long linc_center = gkyl_range_idx(update_range, iter_vertex.idx);

    for (int i=0; i<sz[ndim-1]; ++i) {
      em_tot_d[i] =  gkyl_array_cfetch(em_tot, linc_center + offsets_vertices[i]);
      fluid_d[i] = gkyl_array_cfetch(fluid, linc_center + offsets_vertices[i]);
    }

    heat_flux_d = gkyl_array_fetch(heat_flux, linc_vertex);

    calc_unmag_heat_flux(gces, fluid_d, gkyl_array_fetch(cflrate, linc_center), heat_flux_d);
  }

  struct gkyl_range_iter iter_center;
  gkyl_range_iter_init(&iter_center, update_range);
  while (gkyl_range_iter_next(&iter_center)) {

    long linc_vertex = gkyl_range_idx(heat_flux_range, iter_center.idx);
    long linc_center = gkyl_range_idx(update_range, iter_center.idx);

    for (int i=0; i<sz[ndim-1]; ++i)
      heat_flux_up[i] = gkyl_array_fetch(heat_flux, linc_vertex + offsets_centers[i]);

    rhs_d = gkyl_array_fetch(rhs, linc_center);

    calc_grad_closure_update(gces, heat_flux_up, rhs_d);
  }
}

void
gkyl_ten_moment_grad_closure_release(gkyl_ten_moment_grad_closure* up)
{
  free(up);
}
