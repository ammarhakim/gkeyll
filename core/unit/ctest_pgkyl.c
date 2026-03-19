#include <acutest.h>
#include <stdio.h>

#include <pgkyl_ifc.h>

struct cmd_info_inp {
  int val;
};

static bool
cmd_info(pgkyl_state ps, void *ctx)
{
  struct cmd_info_inp *inp = ctx;
  
//  fprintf(stdout, "Inside info (val = %d)!\n", inp->val);
  return true;
}

struct cmd_interp_inp {
  int ninterp;
};

static bool
cmd_interp(pgkyl_state ps, void *ctx)
{
  struct cmd_interp_inp *inp = ctx;
  
  //fprintf(stdout, "Inside interp (niterp = %d)!\n", inp->ninterp);
  return true;
}

struct cmd_write_inp {
  const char *fbase;
};

static bool
cmd_write(pgkyl_state ps, void *ctx)
{
  struct cmd_write_inp *inp = ctx;
  
  //fprintf(stdout, "Inside write (fbase = %s)!\n", inp->fbase);
  return true;
}

struct cmd_select_inp {
  int comp;
};

static bool
cmd_select(pgkyl_state ps, void *ctx)
{
  struct cmd_select_inp *inp = ctx;
  
  //fprintf(stdout, "Inside select (comp = %d)!\n", inp->comp);
  return true;
}

void
test_basic(void)
{
  pgkyl_state ps = pgkyl_state_new();

  int shape[] = { 10, 10 };
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);

  double lower[] = { 0.0, 0.0 };
  double upper[] = { 1.0, 1.0 };
  int cells[] = { shape[0], shape[1] };
  
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);
  
  int aid = pgkyl_start_appending(ps);
  pgkyl_append_array(ps, &range, &grid, arr);
  pgkyl_finish_appending(ps);

  TEST_CHECK( aid == pgkyl_get_active_dataset(ps) );

  int bid = pgkyl_start_appending(ps);
  pgkyl_append_array(ps, &range, &grid, arr);
  pgkyl_append_array(ps, &range, &grid, arr);  
  pgkyl_finish_appending(ps);

  TEST_CHECK( bid == pgkyl_get_active_dataset(ps) );
  
  TEST_CHECK( 2 == pgkyl_num_datasets(ps) );
  
  gkyl_array_release(arr);
  pgkyl_state_release(ps);
}

TEST_LIST = {
  { "test_basic", test_basic },
  { NULL, NULL },
};
