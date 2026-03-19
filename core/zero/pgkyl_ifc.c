#include <gkyl_alloc.h>
#include <pgkyl_ifc.h>

// Maximum number of datasets: this should be arbitrary and should be
// changed in the future
#define PGKYL_MAX_DATASETS 124 

static int
array_cmp(struct pgkyl_array_data *a, struct pgkyl_array_data *b)
{
  return 0;
}

static void
array_drop(struct pgkyl_array_data *a)
{
  gkyl_array_release(a->array);
}

static struct pgkyl_array_data
array_valfrom(struct pgkyl_array_data a)
{
  return a;
}

// define vector of struct pgkyl_array_data
#define i_val struct pgkyl_array_data
#define i_cmp array_cmp
#define i_drop array_drop
#define i_valfrom array_valfrom
#define i_tag array_p
#include <stc/cvec.h>
// .. finish vector definition

struct pgkyl_dataset {
  cvec_array_p arrays;
};

struct pgkyl_state_tag {
  int last_dset_id; // last inserted dataset ID
  int active_dset_id; // active dataset ID
  
  struct pgkyl_dataset datasets[PGKYL_MAX_DATASETS];
};

pgkyl_state
pgkyl_state_new(void)
{
  struct pgkyl_state_tag *ps = gkyl_malloc(sizeof *ps);
  ps->last_dset_id = 0;
  ps->active_dset_id = 0;

  for (int i=0; i<PGKYL_MAX_DATASETS; ++i)
    ps->datasets[i].arrays = cvec_array_p_init();
  
  return ps;
}

int
pgkyl_start_appending(pgkyl_state ps)
{
  return ps->last_dset_id;
}

void
pgkyl_finish_appending(pgkyl_state ps)
{
  ps->active_dset_id = ps->last_dset_id;
  ps->last_dset_id += 1;  
}

int
pgkyl_num_datasets(const pgkyl_state ps)
{
  return ps->last_dset_id;
}

bool
pgkyl_set_active_dataset(pgkyl_state ps, int dsid)
{
  if (dsid < PGKYL_MAX_DATASETS) {
    ps->active_dset_id = dsid;
    return true;
  }
  return false;
}

int
pgkyl_get_active_dataset(const pgkyl_state ps)
{
  return ps->active_dset_id;
}

void
pgkyl_append_array(pgkyl_state ps, const struct gkyl_range *range,
  const struct gkyl_rect_grid *grid, const struct gkyl_array *arr)
{
  struct pgkyl_array_data ad = {
    .range = *range,
    .grid = *grid,
    .array = gkyl_array_acquire(arr)
  };
  int dsid = ps->last_dset_id;
  cvec_array_p_push_back(&ps->datasets[dsid].arrays, ad);
}

struct pgkyl_cmd_status
pgkyl_execute_cmd(pgkyl_state ps, struct pgkyl_cmd_ctx_pair cmd)
{
  bool s = cmd.cmd(ps, cmd.ctx);
  return (struct pgkyl_cmd_status) { .status = s, .outid = ps->last_dset_id };
}

struct pgkyl_cmd_status
pgkyl_execute_cmd_chain(pgkyl_state ps, const struct pgkyl_cmd_ctx_pair cmd_list[])
{
  bool s = true;
  int lastid = ps->last_dset_id;
  
  for (int i=0; cmd_list[i].cmd != 0; ++i) {
    struct pgkyl_cmd_status stid = pgkyl_execute_cmd(ps, cmd_list[i]);
    if (!stid.status) {
      s = false;
      break;
    }
    lastid = stid.outid; // next command in chain uses last command's output ID
  }
  
  return (struct pgkyl_cmd_status) { .status = s, .outid = lastid };
}

void
pgkyl_state_release(pgkyl_state ps)
{
  for (int i=0; i<PGKYL_MAX_DATASETS; ++i)
    cvec_array_p_drop(&ps->datasets[i].arrays);
  gkyl_free(ps);
}     
