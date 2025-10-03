#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>

// Object type.
typedef struct gkyl_rescale_ghost_jacf gkyl_rescale_ghost_jacf;

/**
 * Create an updater which takes the phase-space quantity in the ghost cell,
 * divides it by the ghost cell of a conf-space quantity, and multiplies it
 * by the skin-cell value of the same ghost cell quantity but flipped in the
 * direction of the boundary (so that the conf-space quantity we multiply by
 * in the ghost cell has the same value at the boundary as it does in the skin
 * cell).
 *
 * @param dir Direction perpendicular to the boundary.
 * @param edge Boundary edge (lower/upper).
 * @param conf_basis Configuration space basis object.
 * @param phase_basis Phase space basis object.
 * @param use_gpu Whether to run on the GPU.
 * @return New rescale_ghost_jacf updater. 
 */
struct gkyl_rescale_ghost_jacf* gkyl_rescale_ghost_jacf_new(int dir,
  enum gkyl_edge_loc edge, const struct gkyl_basis *conf_basis,
  const struct gkyl_basis *phase_basis, bool use_gpu);

/**
 * Run the rescale_ghost_jacf updater.
 *
 * @param up Updater of type rescale_ghost_jacf. 
 * @param conf_skin_r Configuration space skin range.
 * @param conf_ghost_r Configuration space ghost range.
 * @param phase_ghost_r Phase space ghost range.
 * @param jac Conf-space field (e.g. the jacobian).
 * @param jf Phase-space field (e.g. jacobian*distribution).
 */
void gkyl_rescale_ghost_jacf_advance(const struct gkyl_rescale_ghost_jacf *up,
  const struct gkyl_range *conf_skin_r, const struct gkyl_range *conf_ghost_r,
  const struct gkyl_range *phase_ghost_r, const struct gkyl_array *jac_sync, struct gkyl_array *jf);

/**
 * Release the memory associated with rescale_ghost_jacf updater.
 *
 * @param up Updater of type rescale_ghost_jacf. 
 */
void gkyl_rescale_ghost_jacf_release(struct gkyl_rescale_ghost_jacf *up);
