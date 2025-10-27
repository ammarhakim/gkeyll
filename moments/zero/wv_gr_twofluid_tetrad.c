#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_twofluid_tetrad.h>
#include <gkyl_wv_gr_twofluid_tetrad_priv.h>

void
gkyl_gr_twofluid_tetrad_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(base->on_dev, struct wv_gr_twofluid_tetrad, eqn);
    gkyl_cu_free(gr_twofluid_tetrad);
  }

  struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(base, struct wv_gr_twofluid_tetrad, eqn);
  gkyl_free(gr_twofluid_tetrad);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_twofluid_tetrad_new(double mass_elc, double mass_ion, double charge_elc, double charge_ion, double gas_gamma_elc, double gas_gamma_ion,
  double light_speed, double e_fact, double b_fact, enum gkyl_spacetime_gauge spacetime_gauge, int reinit_freq, struct gkyl_gr_spacetime* spacetime,
  bool use_gpu)
{
  return gkyl_wv_gr_twofluid_tetrad_inew(&(struct gkyl_wv_gr_twofluid_tetrad_inp) {
      .mass_elc = mass_elc,
      .mass_ion = mass_ion,
      .charge_elc = charge_elc,
      .charge_ion = charge_ion,
      .gas_gamma_elc = gas_gamma_elc,
      .gas_gamma_ion = gas_gamma_ion,
      .light_speed = light_speed,
      .e_fact = e_fact,
      .b_fact = b_fact,
      .spacetime_gauge = spacetime_gauge,
      .reinit_freq = reinit_freq,
      .spacetime = spacetime,
      .rp_type = WV_GR_TWOFLUID_TETRAD_RP_LAX,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_gr_twofluid_tetrad_inew(const struct gkyl_wv_gr_twofluid_tetrad_inp* inp)
{
  struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = gkyl_malloc(sizeof(struct wv_gr_twofluid_tetrad));

  gr_twofluid_tetrad->eqn.type = GKYL_EQN_GR_TWOFLUID_TETRAD;
  gr_twofluid_tetrad->eqn.num_equations = 84;
  gr_twofluid_tetrad->eqn.num_diag = 5;

  gr_twofluid_tetrad->mass_elc = inp->mass_elc;
  gr_twofluid_tetrad->mass_ion = inp->mass_ion;
  gr_twofluid_tetrad->charge_elc = inp->charge_elc;
  gr_twofluid_tetrad->charge_ion = inp->charge_ion;
  gr_twofluid_tetrad->gas_gamma_elc = inp->gas_gamma_elc;
  gr_twofluid_tetrad->gas_gamma_ion = inp->gas_gamma_ion;
  gr_twofluid_tetrad->light_speed = inp->light_speed;
  gr_twofluid_tetrad->e_fact = inp->e_fact;
  gr_twofluid_tetrad->b_fact = inp->b_fact;

  gr_twofluid_tetrad->spacetime_gauge = inp->spacetime_gauge;
  gr_twofluid_tetrad->reinit_freq = inp->reinit_freq;
  gr_twofluid_tetrad->spacetime = inp->spacetime;

  gr_twofluid_tetrad->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_twofluid_tetrad->eqn.flags);
  gr_twofluid_tetrad->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_twofluid_tetrad_free);
  gr_twofluid_tetrad->eqn.on_dev = &gr_twofluid_tetrad->eqn; // On the CPU, the equation object points to itself.

  return &gr_twofluid_tetrad->eqn;
}

double
gkyl_wv_gr_twofluid_tetrad_mass_elc(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double mass_elc = gr_twofluid_tetrad->mass_elc;

  return mass_elc;
}

double
gkyl_wv_gr_twofluid_tetrad_mass_ion(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double mass_ion = gr_twofluid_tetrad->mass_ion;

  return mass_ion;
}

double
gkyl_wv_gr_twofluid_tetrad_charge_elc(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double charge_elc = gr_twofluid_tetrad->charge_elc;

  return charge_elc;
}

double
gkyl_wv_gr_twofluid_tetrad_charge_ion(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double charge_ion = gr_twofluid_tetrad->charge_ion;

  return charge_ion;
}

double
gkyl_wv_gr_twofluid_tetrad_gas_gamma_elc(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double gas_gamma_elc = gr_twofluid_tetrad->gas_gamma_elc;

  return gas_gamma_elc;
}

double
gkyl_wv_gr_twofluid_tetrad_gas_gamma_ion(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double gas_gamma_ion = gr_twofluid_tetrad->gas_gamma_ion;

  return gas_gamma_ion;
}

double
gkyl_wv_gr_twofluid_tetrad_light_speed(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double light_speed = gr_twofluid_tetrad->light_speed;

  return light_speed;
}

double
gkyl_wv_gr_twofluid_tetrad_e_fact(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double e_fact = gr_twofluid_tetrad->e_fact;

  return e_fact;
}

double
gkyl_wv_gr_twofluid_tetrad_b_fact(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double b_fact = gr_twofluid_tetrad->b_fact;

  return b_fact;
}

enum gkyl_spacetime_gauge
gkyl_wv_gr_twofluid_tetrad_spacetime_gauge(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  enum gkyl_spacetime_gauge spacetime_gauge = gr_twofluid_tetrad->spacetime_gauge;

  return spacetime_gauge;
}

int
gkyl_wv_gr_twofluid_tetrad_reinit_freq(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  int reinit_freq = gr_twofluid_tetrad->reinit_freq;

  return reinit_freq;
}

struct gkyl_gr_spacetime*
gkyl_wv_gr_twofluid_tetrad_spacetime(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  struct gkyl_gr_spacetime *spacetime = gr_twofluid_tetrad->spacetime;

  return spacetime;
}
