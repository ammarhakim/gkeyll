#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Object type
typedef struct gkyl_spitzer_coll_freq gkyl_spitzer_coll_freq;


/**
 * Calculate normNu based on reference parameters for species s and r
 * @param ns density of species s
 * @param nr density of species r
 * @param ms mass of species s
 * @param mr mass of species r
 * @param qs charge of species s
 * @param qr charge of species r
 * @param Ts Temperature of species s
 * @param Tr Temperature of species r
 */
double gkyl_calc_norm_nu(double ns, double nr, double ms, double mr, double qs, double qr, double Ts, double Tr, double bmag_mid, double eps0, double hbar, double eV);

/**
 * Create new updater to either compute the Spitzer collision frequency from
 * scratch based on local parameters, or scale a normalized collision frequency
 * by the local n_r/(v_ts^2+v_tr^2)^(3/2).
 *
 * @param basis Basis object (configuration space).
 * @param num_quad Number of quadrature nodes.
 * @param use_gpu Boolean indicating whether to use the GPU.
 * @param nufrac Fraction to multiply the collision frequency by.
 * @param eps0 Permittivity of vacuum.
 * @param hbar Planck's constant divided by 2*pi.
 * @return New updater pointer.
 */
gkyl_spitzer_coll_freq* gkyl_spitzer_coll_freq_new(const struct gkyl_basis *basis,
  int num_quad, double nufrac, double eps0, double hbar, bool use_gpu);

/**
 * Scale the normalized collision frequency, normNu, by
 * n_r/(v_ts^2+v_tr^2)^(3/2) and project it on to the basis.
 *
 * @param up Spizer collision frequency updater object.
 * @param range Config-space rang.e
 * @param vtSqSelf Thermal speed squared of this species. 
 * @param vtSqMinSelf Minimum vtSq of this species supported by the grid.
 * @param m0Other Thermal speed squared of the other species. 
 * @param vtSqOther Thermal speed squared of the other species. 
 * @param vtSqMinOther Minimum vtSq of the other species supported by the grid.
 * @param normNu Normalized collision frequency to scale.
 * @param nuOut Output collision frequency.
 */
void gkyl_spitzer_coll_freq_advance_normnu(const gkyl_spitzer_coll_freq *up,
  const struct gkyl_range *range, const struct gkyl_array *vtSqSelf, double vtSqMinSelf,
  const struct gkyl_array *m0Other, const struct gkyl_array *vtSqOther, double vtSqMinOther,
  double normNu, struct gkyl_array *nuOut);

/**
 * Compute the Spitzer collision frequency from scratch. Coulomb Logarithm
 * is computed using cell averaged values.
 *
 * @param up Spizer collision frequency updater object.
 * @param range Config-space range.
 * @param bmag Magnetic field magnitude.
 * @param qSelf Charge of this species.
 * @param mSelf Mass of this species.
 * @param m0Self Thermal speed squared of the other species. 
 * @param vtSqSelf Thermal speed squared of this species. 
 * @param vtSqMinSelf Minimum vtSq of this species supported by the grid.
 * @param qOther Charge of this species.
 * @param mOther Mass of this species.
 * @param m0Other Thermal speed squared of the other species. 
 * @param vtSqOther Thermal speed squared of the other species. 
 * @param vtSqMinOther Minimum vtSq of the other species supported by the grid.
 * @param nuOut Output collision frequency.
 */
void gkyl_spitzer_coll_freq_advance(const gkyl_spitzer_coll_freq *up,
  const struct gkyl_range *range, const struct gkyl_array *bmag,
  double qSelf, double mSelf, const struct gkyl_array *m0Self, const struct gkyl_array *vtSqSelf, double vtSqMinSelf,
  double qOther, double mOther, const struct gkyl_array *m0Other, const struct gkyl_array *vtSqOther, double vtSqMinOther,
  struct gkyl_array *nuOut);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_spitzer_coll_freq_release(gkyl_spitzer_coll_freq* up);
