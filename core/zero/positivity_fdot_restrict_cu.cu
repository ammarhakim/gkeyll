extern "C" {
#include <gkyl_positivity_fdot_restrict.h>
#include <gkyl_positivity_fdot_restrict_priv.h>

#include <gkyl_alloc.h>
#include <gkyl_range.h>
}

#include <float.h>
#include <math.h>

// 1D kernel for AVG-type restrictions: one thread per cell.
__global__ static void
gkyl_positivity_fdot_restrict_avg_cu_kern(gkyl_positivity_fdot_restrict* up,
  const struct gkyl_range *range, const struct gkyl_array *f,
  struct gkyl_array *dfdt, double dt)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= range->volume) return;

  const double *fptr = (const double*) gkyl_array_cfetch(f, idx);
  double *dfdtptr = (double*) gkyl_array_fetch(dfdt, idx);

  up->fdot_restriction_func(up, fptr, dfdtptr, dt);
}

// 2D kernel for QUAD-type restriction: threadIdx.x = cell, threadIdx.y = quadrature node k.
// Parallelizes over both cells and basis functions simultaneously.
//
// Shared memory layout per block (blockDim.x = cells_per_block, blockDim.y = nquad):
//   [fquad     : cells_per_block * nquad doubles]
//   [dfdt_quad : cells_per_block * nquad doubles]
//   [reproject : cells_per_block ints            ]
__global__ static void
gkyl_positivity_fdot_restrict_quad_cu_kern(const gkyl_positivity_fdot_restrict* up,
  const struct gkyl_range *range, const struct gkyl_array *f,
  struct gkyl_array *dfdt, double dt)
{
  int cell_idx = blockIdx.x * blockDim.x + threadIdx.x;
  int k = threadIdx.y;
  int nquad = up->basis.num_quad;
  bool in_bounds = (cell_idx < range->volume);

  extern __shared__ char shmem[];
  double *fquad_sh     = (double*)shmem + threadIdx.x * nquad;
  double *dfdt_quad_sh = (double*)shmem + blockDim.x * nquad + threadIdx.x * nquad;
  int    *reproject    = (int*)(shmem + 2 * blockDim.x * nquad * sizeof(double)) + threadIdx.x;

  if (k == 0) *reproject = 0;
  __syncthreads();

  if (in_bounds) {
    const double *fptr    = (const double*) gkyl_array_cfetch(f,    cell_idx);
    double       *dfdtptr = (double*)       gkyl_array_fetch(dfdt, cell_idx);
    up->basis.modal_to_quad_nodal(fptr,    fquad_sh,    k);
    up->basis.modal_to_quad_nodal(dfdtptr, dfdt_quad_sh, k);
  }
  __syncthreads();

  if (in_bounds) {
    double f_pred = fquad_sh[k] + dt * dfdt_quad_sh[k];
    if (f_pred < 0.0) {
      dfdt_quad_sh[k] = (up->safety_factor - 1.0) * fquad_sh[k] / dt;
      atomicOr(reproject, 1);
    }
  }
  __syncthreads();

  if (in_bounds && *reproject)
    up->basis.quad_nodal_to_modal(dfdt_quad_sh, (double*)gkyl_array_fetch(dfdt, cell_idx), k);
}

// 2D kernel for DIODE_QUAD-type restriction: threadIdx.x = cell, threadIdx.y = quadrature node k.
// Same shared memory layout as gkyl_positivity_fdot_restrict_quad_cu_kern.
__global__ static void
gkyl_positivity_fdot_restrict_diode_quad_cu_kern(const gkyl_positivity_fdot_restrict* up,
  const struct gkyl_range *range, const struct gkyl_array *f,
  struct gkyl_array *dfdt, double dt)
{
  int cell_idx = blockIdx.x * blockDim.x + threadIdx.x;
  int k = threadIdx.y;
  int nquad = up->basis.num_quad;
  bool in_bounds = (cell_idx < range->volume);

  extern __shared__ char shmem[];
  double *fquad_sh     = (double*)shmem + threadIdx.x * nquad;
  double *dfdt_quad_sh = (double*)shmem + blockDim.x * nquad + threadIdx.x * nquad;
  int    *reproject    = (int*)(shmem + 2 * blockDim.x * nquad * sizeof(double)) + threadIdx.x;

  if (k == 0) *reproject = 0;
  __syncthreads();

  if (in_bounds) {
    const double *fptr    = (const double*) gkyl_array_cfetch(f,    cell_idx);
    double       *dfdtptr = (double*)       gkyl_array_fetch(dfdt, cell_idx);
    up->basis.modal_to_quad_nodal(fptr,    fquad_sh,    k);
    up->basis.modal_to_quad_nodal(dfdtptr, dfdt_quad_sh, k);
  }
  __syncthreads();

  if (in_bounds && fquad_sh[k] < 0.0 && dfdt_quad_sh[k] < 0.0) {
    dfdt_quad_sh[k] = 0.0;
    atomicOr(reproject, 1);
  }
  __syncthreads();

  if (in_bounds && *reproject)
    up->basis.quad_nodal_to_modal(dfdt_quad_sh, (double*)gkyl_array_fetch(dfdt, cell_idx), k);
}

// Per-type host-side launcher functions. Set as function pointers at init time.

static void
restrict_fdot_avg_cu_launch(struct gkyl_positivity_fdot_restrict *up,
  const struct gkyl_range *range, const struct gkyl_array *f,
  struct gkyl_array *dfdt, double dt)
{
  gkyl_positivity_fdot_restrict_avg_cu_kern<<<range->nblocks, range->nthreads>>>(
    up->on_dev, range, f->on_dev, dfdt->on_dev, dt);
}

static void
restrict_fdot_quad_cu_launch(struct gkyl_positivity_fdot_restrict *up,
  const struct gkyl_range *range, const struct gkyl_array *f,
  struct gkyl_array *dfdt, double dt)
{
  int nquad = up->basis.num_quad;
  int cells_per_block = max(1, 512 / nquad);
  dim3 block_dim(cells_per_block, nquad);
  dim3 grid_dim((range->volume + cells_per_block - 1) / cells_per_block, 1);
  size_t shmem_size = cells_per_block * (2 * nquad * sizeof(double) + sizeof(int));
  gkyl_positivity_fdot_restrict_quad_cu_kern<<<grid_dim, block_dim, shmem_size>>>(
    up->on_dev, range, f->on_dev, dfdt->on_dev, dt);
}

static void
restrict_fdot_diode_avg_cu_launch(struct gkyl_positivity_fdot_restrict *up,
  const struct gkyl_range *range, const struct gkyl_array *f,
  struct gkyl_array *dfdt, double dt)
{
  gkyl_positivity_fdot_restrict_avg_cu_kern<<<range->nblocks, range->nthreads>>>(
    up->on_dev, range, f->on_dev, dfdt->on_dev, dt);
}

static void
restrict_fdot_diode_quad_cu_launch(struct gkyl_positivity_fdot_restrict *up,
  const struct gkyl_range *range, const struct gkyl_array *f,
  struct gkyl_array *dfdt, double dt)
{
  int nquad = up->basis.num_quad;
  int cells_per_block = max(1, 512 / nquad);
  dim3 block_dim(cells_per_block, nquad);
  dim3 grid_dim((range->volume + cells_per_block - 1) / cells_per_block, 1);
  size_t shmem_size = cells_per_block * (2 * nquad * sizeof(double) + sizeof(int));
  gkyl_positivity_fdot_restrict_diode_quad_cu_kern<<<grid_dim, block_dim, shmem_size>>>(
    up->on_dev, range, f->on_dev, dfdt->on_dev, dt);
}

void
gkyl_positivity_fdot_restrict_cu_set_ptrs(struct gkyl_positivity_fdot_restrict *up)
{
  switch (up->type) {
    case GKYL_POSITIVITY_FDOT_RESTRICT_AVG:
      up->advance_cu_func = restrict_fdot_avg_cu_launch;
      break;
    case GKYL_POSITIVITY_FDOT_RESTRICT_QUAD:
      up->advance_cu_func = restrict_fdot_quad_cu_launch;
      break;
    case GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG:
      up->advance_cu_func = restrict_fdot_diode_avg_cu_launch;
      break;
    case GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD:
      up->advance_cu_func = restrict_fdot_diode_quad_cu_launch;
      break;
  }
}

void
gkyl_positivity_fdot_restrict_advance_cu(gkyl_positivity_fdot_restrict* up,
  const struct gkyl_range *range, const struct gkyl_array *f,
  struct gkyl_array *dfdt, double dt)
{
  up->advance_cu_func(up, range, f, dfdt, dt);
}
