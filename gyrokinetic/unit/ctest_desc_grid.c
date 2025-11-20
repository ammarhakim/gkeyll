void test_3xp1(){
    // Create Nodal Range and Grid and Write Nodal Coordinates (For corner nodes)
    struct gkyl_range nrange;
    gkyl_gk_geometry_init_nodal_range(&nrange, &app->global, app->poly_order);
    struct gkyl_array* mc2p_nodal = mkarr(false, 3, nrange.volume);
    struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&app->basis, &app->grid, false);
    gkyl_nodal_ops_m2n(n2m, &app->basis, &app->grid, &nrange, &app->global, 3, mc2p_nodal, mc2p_global_ho, false);
    struct gkyl_rect_grid ngrid;
    gkyl_gk_geometry_init_nodal_grid(&ngrid, &app->grid, &nrange);
    struct gkylarray* tangentscorn = mkarr(......, nrange.volume*9)
    gkyl_grid_sub_array_read(nrange, ngrid, tangentscorn, 'tangetns.gkyl')


    // Do for interior
    // Then calculate derived quantities
    while(gkyl_range_iter_next(nrange_interior) {
        jacobgeo_n = e_1 x e_2 x e_3
    }


    // Do again for surfaces
    // Derived quantities


    // Do a n2m (use appropriate for interior,corner, and surface) on all neede quantities
    // write out expansions on the modal grid
}