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
  
  fprintf(stdout, "Inside info (val = %d)!\n", inp->val);
  return true;
}

struct cmd_interp_inp {
  int ninterp;
};

static bool
cmd_interp(pgkyl_state ps, void *ctx)
{
  struct cmd_interp_inp *inp = ctx;
  
  fprintf(stdout, "Inside interp (niterp = %d)!\n", inp->ninterp);
  return true;
}

struct cmd_write_inp {
  const char *fbase;
};

static bool
cmd_write(pgkyl_state ps, void *ctx)
{
  struct cmd_write_inp *inp = ctx;
  
  fprintf(stdout, "Inside write (fbase = %s)!\n", inp->fbase);
  return true;
}

struct cmd_select_inp {
  int comp;
};

static bool
cmd_select(pgkyl_state ps, void *ctx)
{
  struct cmd_select_inp *inp = ctx;
  
  fprintf(stdout, "Inside select (comp = %d)!\n", inp->comp);
  return true;
}

void
test_basic(void)
{
  pgkyl_state ps = pgkyl_state_new();

  /* pgkyl_start_appending(ps); */
  /* pgkyl_append_array(ps, &range, &grid, array); */
  /* pgkyl_append_array(ps, &range, &grid, array); */
  /* pgkyl_stop_appending(ps); */

  struct pgkyl_cmd_ctx_pair cmd_list[] = {
    { cmd_info, &(struct cmd_info_inp){ .val = 0 } },
    { cmd_interp, &(struct cmd_interp_inp){ .ninterp = 2 } },
    { cmd_select, &(struct cmd_select_inp){ .comp = 1 } },
    { cmd_write, &(struct cmd_write_inp){ .fbase = "output" } },
    { 0, 0 },
  };
  
  pgkyl_execute_cmd_chain(ps, cmd_list);

  /* pgkyl_clear_all_datasets(ps); */
  
  pgkyl_state_release(ps);
}

TEST_LIST = {
  { "test_basic", test_basic },
  { NULL, NULL },
};
