#include <math.h>
#include <assert.h>
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_coldfluid.h>
#include <gkyl_wv_coldfluid_priv.h>


static void
coldfluid_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_coldfluid *coldfluid = container_of(base, struct wv_coldfluid, eqn);
  gkyl_free(coldfluid);
}

struct gkyl_wv_eqn*
gkyl_wv_coldfluid_inew(const struct gkyl_wv_coldfluid_inp *inp)
{
  struct wv_coldfluid *coldfluid = gkyl_malloc(sizeof(struct wv_coldfluid));
  coldfluid->eqn.type = GKYL_EQN_COLDFLUID;
  coldfluid->eqn.num_equations = 4;
  coldfluid->eqn.num_waves = 2;
  coldfluid->eqn.num_diag = 5; // KE is final component
  
  coldfluid->eqn.waves_func = wave_roe_l;
  coldfluid->eqn.qfluct_func = qfluct_roe;
  coldfluid->eqn.ffluct_func = ffluct_roe;

  coldfluid->eqn.flux_jump = flux_jump;
  coldfluid->eqn.check_inv_func = check_inv;
  coldfluid->eqn.max_speed_func = max_speed;
  coldfluid->eqn.rotate_to_local_func = rot_to_local;
  coldfluid->eqn.rotate_to_global_func = rot_to_global;

  coldfluid->eqn.cons_to_riem = cons_to_riem;
  coldfluid->eqn.riem_to_cons = riem_to_cons;

  coldfluid->eqn.cons_to_diag = coldfluid_cons_to_diag;

  coldfluid->eqn.source_func = coldfluid_source;

  coldfluid->eqn.ref_count = gkyl_ref_count_init(coldfluid_free);

  coldfluid->eqn.embed_geo = inp->embed_geo;
  if (coldfluid->eqn.embed_geo) {
    switch (coldfluid->eqn.embed_geo->type) {
      case GKYL_EMBED_ABSORB:
        coldfluid->eqn.embed_geo->embed_func = wave_embed_absorb;
        break;

      case GKYL_EMBED_REFLECT:
        coldfluid->eqn.embed_geo->embed_func = wave_embed_reflect;
        break;

    case GKYL_EMBED_FUNC:
        break; // already set by gkyl_wv_embed_geo_new

      default:
        assert(false);
        break;
    }
  } 

  return &coldfluid->eqn;
}

struct gkyl_wv_eqn*
gkyl_wv_coldfluid_new(void)
{
  return gkyl_wv_coldfluid_inew( &(struct gkyl_wv_coldfluid_inp) {
    }
  );
}
