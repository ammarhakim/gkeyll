#include <gkyl_alloc.h>
#include <pgkyl_ifc.h>

struct pgkyl_state_tag {
  int last_dset_id; // last inserted dataset ID
};

pgkyl_state
pgkyl_state_new(void)
{
  struct pgkyl_state_tag *ps = gkyl_malloc(sizeof *ps);
  ps->last_dset_id = 0;

  return ps;
}


struct pgkyl_cmd_status
pgkyl_execute_cmd(pgkyl_state ps, int dsetid, struct pgkyl_cmd_ctx_pair cmd)
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
    struct pgkyl_cmd_status stid = pgkyl_execute_cmd(ps, lastid, cmd_list[i]);
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
  gkyl_free(ps);
}     
