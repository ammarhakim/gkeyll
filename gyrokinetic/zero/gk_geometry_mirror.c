#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <assert.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mirror.h>
#include <gkyl_gk_geometry_mirror_priv.h>
#include <gkyl_position_map.h>
#include <gkyl_mirror_grid_gen.h>
#include <gkyl_mirror_geo_gen.h>
#include <gkyl_mirror_geo_dg.h>

struct gk_geometry*
gk_geometry_mirror_init(struct gkyl_gk_geometry_inp *geometry_inp)
{
  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
  up->basis = geometry_inp->geo_basis;
  up->local = geometry_inp->geo_local;
  up->local_ext = geometry_inp->geo_local_ext;
  up->global = geometry_inp->geo_global;
  up->global_ext = geometry_inp->geo_global_ext;
  up->grid = geometry_inp->geo_grid;
  up->geqdsk_sign_convention = 0.0; // Hardcoded 0. Means psi increases from axis. Always true for mirror geometry.

  // read psi(R,Z) from file
  struct gkyl_rect_grid psi_grid;
  struct gkyl_array *psi = gkyl_grid_array_new_from_file(&psi_grid, geometry_inp->mirror_grid_info.filename_psi);

  // Check if file exists using gkyl_check_file_exists and handle error only on rank 0
  if (!gkyl_check_file_exists(geometry_inp->mirror_grid_info.filename_psi)) {
    fprintf(stderr, "Failed to open the eqdsk file: %s\n", geometry_inp->mirror_grid_info.filename_psi);
    assert(false);
  }

  // The mirror grid generators always require a computational grid of psi, phi, z.
  // We must convert the user giving sqrt(psi) to psi if the flag is set.
  struct gkyl_rect_grid comp_grid;
  memcpy(&comp_grid, &up->grid, sizeof(struct gkyl_rect_grid));
  if (geometry_inp->mirror_grid_info.fl_coord == GKYL_MIRROR_GRID_GEN_SQRT_PSI_CART_Z) {
    // Convert grid from sqrt(psi) to psi if the flag is set
    comp_grid.lower[0] = pow(comp_grid.lower[0], 2);
    comp_grid.upper[0] = pow(comp_grid.upper[0], 2);
    comp_grid.dx[0] = (comp_grid.upper[0] - comp_grid.lower[0]) / comp_grid.cells[0];

    comp_grid.cellVolume = 1.0;
    for (int i=0; i<3; ++i) // Geometry is always generated in 3D
      comp_grid.cellVolume *= comp_grid.dx[i];
  }

  // create mirror geometry
  struct gkyl_mirror_grid_gen *mirror_grid =
    gkyl_mirror_grid_gen_inew(&(struct gkyl_mirror_grid_gen_inp) {
        .comp_grid = &comp_grid,
        
        .R = { psi_grid.lower[0], psi_grid.upper[0] },
        .Z = { psi_grid.lower[1], psi_grid.upper[1] },
        
        // psi(R,Z) grid size
        .nrcells = psi_grid.cells[0]-1, // cells and not nodes
        .nzcells = psi_grid.cells[1]-1, // cells and not nodes

        .psiRZ = psi,
        .fl_coord = geometry_inp->mirror_grid_info.fl_coord,
        .include_axis = geometry_inp->mirror_grid_info.include_axis,
        .write_psi_cubic = false,
      }
    );

  // Generate the required geometry objects at each node
  struct gkyl_mirror_geo_gen *mirror_geo = 
    gkyl_mirror_geo_gen_inew(&(struct gkyl_mirror_geo_gen_inp) {
        .comp_grid = &comp_grid,
        .mirror_grid = mirror_grid,
        .range = up->global,
        .basis = up->basis,
      }
    );

  // Project the geometry objects from the nodes onto DG arrays
  struct gkyl_mirror_geo_dg *mirror_geo_dg = 
    gkyl_mirror_geo_dg_inew(&(struct gkyl_mirror_geo_dg_inp) {
        .comp_grid = &comp_grid,
        .mirror_geo = mirror_geo,
        .range = up->global,
        .range_ext = up->global_ext,
        .basis = up->basis,
      }
    );

  // Global geometry is generated on all cores. Copy the local geometry to each process
  if (up->local.volume == up->global.volume) // If local = global
  {
    gk_geometry_mirror_array_acquire(up, mirror_geo_dg);
  }
  else // Code is running in parallel and local != global
  {
    gk_geometry_mirror_array_new(up);
    gk_geometry_mirror_array_copy(up, mirror_geo_dg);
  }


  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself

  gkyl_mirror_grid_gen_release(mirror_grid);
  gkyl_mirror_geo_gen_release(mirror_geo);
  gkyl_mirror_geo_dg_release(mirror_geo_dg);
  gkyl_array_release(psi);
  return up;
}

struct gk_geometry*
gkyl_gk_geometry_mirror_new(struct gkyl_gk_geometry_inp *geometry_inp)
{
  struct gk_geometry* gk_geom_3d;
  struct gk_geometry* gk_geom;

  if (geometry_inp->position_map == 0){
    geometry_inp->position_map = gkyl_position_map_new((struct gkyl_position_map_inp) {}, \
      geometry_inp->grid, geometry_inp->local, geometry_inp->local_ext, geometry_inp->local, \
      geometry_inp->local_ext, geometry_inp->basis);
    gk_geom_3d = gk_geometry_mirror_init(geometry_inp);
    gkyl_position_map_release(geometry_inp->position_map);
  }
  else {
    // First construct the uniform 3d geometry
    gk_geom_3d = gk_geometry_mirror_init(geometry_inp);
    if (geometry_inp->position_map->id == GKYL_PMAP_CONSTANT_DB_POLYNOMIAL || \
        geometry_inp->position_map->id == GKYL_PMAP_CONSTANT_DB_NUMERIC) {
      // The array mc2nu is computed using the uniform geometry, so we need to deflate it
      // Must deflate the 3D uniform geometry in order for the allgather to work
      if(geometry_inp->grid.ndim < 3)
        gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, geometry_inp);
      else
        gk_geom = gkyl_gk_geometry_acquire(gk_geom_3d);

      gkyl_position_map_set_bmag(geometry_inp->position_map, geometry_inp->comm, \
        gk_geom->bmag);

      gkyl_gk_geometry_release(gk_geom_3d); // release temporary 3d geometry
      gkyl_gk_geometry_release(gk_geom); // release 3d geometry

      // Construct the non-uniform grid
      gk_geom_3d = gk_geometry_mirror_init(geometry_inp);
    }
  }
  return gk_geom_3d;
}