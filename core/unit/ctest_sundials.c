#include <acutest.h>
#include <mpack.h>

#include <gkyl_sundials_nvector.h>
#include <gkyl_util.h>

static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  // Allocate array (filled with zeros)
  struct gkyl_array* a = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void sundials_init()
{
  struct gkyl_sundials *gk_sundials = gkyl_sundials_new();

  gkyl_sundials_release(gk_sundials);
}

void sundials_nvector_init(bool use_gpu)
{
  // Test creation and destruction.
  int num_cells = 10;
  int ncomp = 3;

  struct gkyl_comm *comm = 0;
  struct gkyl_range local;

  struct gkyl_sundials *gk_sundials = gkyl_sundials_new();

  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, ncomp, num_cells);

  struct gkyl_sundials_nvec *a1_snv = gkyl_sundials_nvec_new(gk_sundials, a1, comm, &local);

  gkyl_sundials_nvec_release(a1_snv);
  gkyl_array_release(a1);
  gkyl_sundials_release(gk_sundials);
}

void sundials_nvector_init_ho()
{
  sundials_nvector_init(false);
}

#ifdef GKYL_HAVE_CUDA

void sundials_nvector_init_dev()
{
  sundials_nvector_init(true);
}

#endif

TEST_LIST = {
  { "sundials_init", sundials_init },
  { "sundials_nvector_init_ho", sundials_nvector_init_ho },
#ifdef GKYL_HAVE_CUDA
  { "sundials_nvector_init_dev", sundials_nvector_init_dev},
#endif
  { NULL, NULL },
};
