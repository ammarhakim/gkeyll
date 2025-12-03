#include <acutest.h>
#include <stdio.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include <errno.h>

#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Private headers (for test files)
#include <gkyl_app_priv.h>
#include <gkyl_elem_type_priv.h>

// Static function to load a .gkyl file and return the data array as nodal grid
// (static means it doesn't need a header file - only visible in this file)
static struct gkyl_array*
load_gkyl_file_to_nodal(const char *filename)
{
    // Prepend relative path to data directory
    // Test is run from gkeyll/ directory, so path is relative to there
    char fullpath[512];
    snprintf(fullpath, sizeof(fullpath), "W7-X_field-aligned_coords_gkyl/%s", filename);

    // Step 1: Read header to get grid info and number of components
    struct gkyl_rect_grid grid;
    struct gkyl_array_header_info hdr;
    enum gkyl_array_rio_status status = gkyl_grid_sub_array_header_read(&grid, &hdr, fullpath);

    if (status != GKYL_ARRAY_RIO_SUCCESS) {
        fprintf(stderr, "Failed to load file: %s (status: %d)\n", fullpath, status);
        return NULL;
    }

    // Step 2: Extract number of components from header
    size_t ncomp = hdr.esznc / gkyl_elem_type_size[hdr.etype];

    // Step 3: Create computational range from grid's cell counts
    // The grid.cells[] array contains cell counts per dimension
    struct gkyl_range comp_range;
    gkyl_range_init_from_shape1(&comp_range, grid.ndim, grid.cells);
    // This uses grid.cells which was read from the file header

    // Step 4: Create nodal range from computational range
    struct gkyl_range nrange;
    int poly_order = 1;
    gkyl_gk_geometry_init_nodal_range(&nrange, &comp_range, poly_order);

    // Step 5: Allocate array with correct number of components
    struct gkyl_array *data_array = mkarr(false, ncomp, nrange.volume);

    // Step 6: Read the actual data
    struct gkyl_rect_grid ngrid;
    gkyl_gk_geometry_init_nodal_grid(&ngrid, &grid, &nrange);
    status = gkyl_grid_sub_array_read(&ngrid, &nrange, data_array, fullpath);

    // Clean up header if needed
    gkyl_array_header_info_release(&hdr);

    if (status != GKYL_ARRAY_RIO_SUCCESS) {
        gkyl_array_release(data_array);
        return NULL;
    }

    return data_array;
}

void test_3xp1(){
    //Load in all corner files and calculate derived quantities
    struct gkyl_array* bmag_corner_nodal = load_gkyl_file_to_nodal("Bmag_corner.gkyl");
    struct gkyl_array* mc2nu_pos_corner_nodal = load_gkyl_file_to_nodal("rza_corner.gkyl");
    struct gkyl_array* mc2p_corner_nodal = load_gkyl_file_to_nodal("rzp_corner.gkyl");
    // Check for NULL pointers (file loading failures)
    TEST_CHECK(bmag_corner_nodal != NULL);
    TEST_CHECK(mc2nu_pos_corner_nodal != NULL);
    TEST_CHECK(mc2p_corner_nodal != NULL);
    if (!bmag_corner_nodal || !mc2p_corner_nodal || !mc2nu_pos_corner_nodal) {
        return; // Exit early if files failed to load
    }
    // Determine number of nodes from loaded array size (bmag has 1 component, so size = number of nodes)
    int shape_corner[1] = {(int)bmag_corner_nodal->size};
    struct gkyl_range nrange_corner;
    gkyl_range_init_from_shape1(&nrange_corner, 1, shape_corner);

    // Do for interior
    struct gkyl_array* dxdz_interior_nodal = load_gkyl_file_to_nodal("drzp_drza_I.gkyl");
    struct gkyl_array* bmag_interior_nodal = load_gkyl_file_to_nodal("Bmag_I.gkyl");
    struct gkyl_array* curlbhat_interior_nodal = load_gkyl_file_to_nodal("curl_B_hat_xyz_I.gkyl");
    struct gkyl_array* bcart_interior_nodal = load_gkyl_file_to_nodal("B_xyz_I.gkyl");
    // LOAD IN bcart_interior_nodal from DESC data
    TEST_CHECK(dxdz_interior_nodal != NULL);
    TEST_CHECK(bmag_interior_nodal != NULL);
    TEST_CHECK(curlbhat_interior_nodal != NULL);
    TEST_CHECK(bcart_interior_nodal != NULL);
    if (!dxdz_interior_nodal || !bmag_interior_nodal || !curlbhat_interior_nodal || !bcart_interior_nodal) {
        return;
    }
    // Determine number of nodes from loaded array size (bmag has 1 component, so size = number of nodes)
    int shape_interior[1] = {(int)bmag_interior_nodal->size};
    struct gkyl_range nrange_interior;
    gkyl_range_init_from_shape1(&nrange_interior, 1, shape_interior);
    
    // Then calculate derived quantities
    struct gkyl_range_iter iter_interior;
    gkyl_range_iter_init(&iter_interior, &nrange_interior);  
    // initialize jacobian and other arrays which are to be calculated (mkarr(false, number of components, nrange_interior.volume))
    struct gkyl_array* B3_interior_nodal = mkarr(false, 1, nrange_interior.volume);
    struct gkyl_array* b_i_interior_nodal = mkarr(false, 3, nrange_interior.volume);
    struct gkyl_array* bmag_inv_interior_nodal = mkarr(false, 1, nrange_interior.volume);
    struct gkyl_array* bmag_inv_sq_interior_nodal = mkarr(false, 1, nrange_interior.volume);
    struct gkyl_array* bioverJB_interior_nodal = mkarr(false, 3, nrange_interior.volume);
    struct gkyl_array* cmag_interior_nodal = mkarr(false, 1, nrange_interior.volume);
    struct gkyl_array* dualcurlbhat_interior_nodal = mkarr(false, 3, nrange_interior.volume);
    struct gkyl_array* dualcurlbhatoverB_interior_nodal = mkarr(false, 3, nrange_interior.volume);
    struct gkyl_array* dzdx_interior_nodal = mkarr(false, 9, nrange_interior.volume);
    struct gkyl_array* eps2_interior_nodal = mkarr(false, 1, nrange_interior.volume);
    struct gkyl_array* g_ij_interior_nodal = mkarr(false, 6, nrange_interior.volume);
    struct gkyl_array* gij_interior_nodal = mkarr(false, 6, nrange_interior.volume);
    struct gkyl_array* gxxj_interior_nodal = mkarr(false, 1, nrange_interior.volume);
    struct gkyl_array* gxyj_interior_nodal = mkarr(false, 1, nrange_interior.volume);
    struct gkyl_array* gxzj_interior_nodal = mkarr(false, 1, nrange_interior.volume);
    struct gkyl_array* gyyj_interior_nodal = mkarr(false, 1, nrange_interior.volume);
    struct gkyl_array* jacobgeo_interior_nodal = mkarr(false, 1, nrange_interior.volume);
    struct gkyl_array* jacobgeo_inv_interior_nodal = mkarr(false, 1, nrange_interior.volume);
    struct gkyl_array* jacobtot_interior_nodal = mkarr(false, 1, nrange_interior.volume);
    struct gkyl_array* jacobtot_inv_interior_nodal = mkarr(false, 1, nrange_interior.volume);
    struct gkyl_array* normals_interior_nodal = mkarr(false, 9, nrange_interior.volume);
    struct gkyl_array* rtg33inv_interior_nodal = mkarr(false, 1, nrange_interior.volume);
    
    while (gkyl_range_iter_next(&iter_interior)) {
        // Get current node index
        long idx = gkyl_range_idx(&nrange_interior, iter_interior.idx);
        
        double *dxdz_data = (double*)dxdz_interior_nodal->data;
        double *bmag_data = (double*)bmag_interior_nodal->data;
        double *curlbhat_data = (double*)curlbhat_interior_nodal->data;
        double *B3_data = (double*)B3_interior_nodal->data;
        double *b_i_data = (double*)b_i_interior_nodal->data;
        double *bcart_data = (double*)bcart_interior_nodal->data;
        double *bmag_inv_data = (double*)bmag_inv_interior_nodal->data;
        double *bmag_inv_sq_data = (double*)bmag_inv_sq_interior_nodal->data;
        double *bioverJB_data = (double*)bioverJB_interior_nodal->data;
        double *cmag_data = (double*)cmag_interior_nodal->data;
        double *dualcurlbhat_data = (double*)dualcurlbhat_interior_nodal->data;
        double *dualcurlbhatoverB_data = (double*)dualcurlbhatoverB_interior_nodal->data;
        double *dzdx_data = (double*)dzdx_interior_nodal->data;
        double *eps2_data = (double*)eps2_interior_nodal->data;
        double *g_ij_data = (double*)g_ij_interior_nodal->data;
        double *gij_data = (double*)gij_interior_nodal->data;
        double *gxxj_data = (double*)gxxj_interior_nodal->data;
        double *gxyj_data = (double*)gxyj_interior_nodal->data;
        double *gxzj_data = (double*)gxzj_interior_nodal->data;
        double *gyyj_data = (double*)gyyj_interior_nodal->data;
        double *jacobgeo_data = (double*)jacobgeo_interior_nodal->data;
        double *jacobgeo_inv_data = (double*)jacobgeo_inv_interior_nodal->data;
        double *jacobtot_data = (double*)jacobtot_interior_nodal->data;
        double *jacobtot_inv_data = (double*)jacobtot_inv_interior_nodal->data;
        double *normals_data = (double*)normals_interior_nodal->data;
        double *rtg33inv_data = (double*)rtg33inv_interior_nodal->data;
        
        double e_1[3] = {dxdz_data[idx*9+0], dxdz_data[idx*9+1], dxdz_data[idx*9+2]};
        double e_2[3] = {dxdz_data[idx*9+3], dxdz_data[idx*9+4], dxdz_data[idx*9+5]};
        double e_3[3] = {dxdz_data[idx*9+6], dxdz_data[idx*9+7], dxdz_data[idx*9+8]};
        double curlb[3] = {curlbhat_data[idx*3+0], curlbhat_data[idx*3+1], curlbhat_data[idx*3+2]};
        double bmag = bmag_data[idx];
        
        // Jacobian: J = e_1 · (e_2 × e_3)
        double e2xe3[3] = {e_2[1]*e_3[2]-e_2[2]*e_3[1], e_2[2]*e_3[0]-e_2[0]*e_3[2], e_2[0]*e_3[1]-e_2[1]*e_3[0]};
        double J = e_1[0]*e2xe3[0] + e_1[1]*e2xe3[1] + e_1[2]*e2xe3[2];
        jacobgeo_data[idx] = J;
        if (fabs(J) < 1e-12) {
            // Jacobian is too small, skip this node or set to a small value
            J = 1e-12;
            jacobgeo_data[idx] = J;
        }
        double Jinv = 1.0/J;
        jacobgeo_inv_data[idx] = Jinv;
        
        // Dual vectors: e^i = (1/J) * (e_j × e_k)
        double e3xe1[3] = {e_3[1]*e_1[2]-e_3[2]*e_1[1], e_3[2]*e_1[0]-e_3[0]*e_1[2], e_3[0]*e_1[1]-e_3[1]*e_1[0]};
        double e1xe2[3] = {e_1[1]*e_2[2]-e_1[2]*e_2[1], e_1[2]*e_2[0]-e_1[0]*e_2[2], e_1[0]*e_2[1]-e_1[1]*e_2[0]};
        for (int i=0; i<3; i++) {
            dzdx_data[idx*9+0*3+i] = Jinv * e2xe3[i];  // e^1
            dzdx_data[idx*9+1*3+i] = Jinv * e3xe1[i];  // e^2
            dzdx_data[idx*9+2*3+i] = Jinv * e1xe2[i];  // e^3
        }
        // Extract dual vectors for convenience
        double e1[3] = {dzdx_data[idx*9+0], dzdx_data[idx*9+1], dzdx_data[idx*9+2]}; // e^1
        double e2[3] = {dzdx_data[idx*9+3], dzdx_data[idx*9+4], dzdx_data[idx*9+5]}; // e^2
        double e3[3] = {dzdx_data[idx*9+6], dzdx_data[idx*9+7], dzdx_data[idx*9+8]}; // e^3
        
        // Metric tensor: g_ij = e_i · e_j (Compute only the 6 upper triangle components)
        // g_ij has 6 components per node (upper triangle): [g_11, g_12, g_13, g_22, g_23, g_33]
        double g_11 = e_1[0]*e_1[0] + e_1[1]*e_1[1] + e_1[2]*e_1[2]; // g_11
        double g_12 = e_1[0]*e_2[0] + e_1[1]*e_2[1] + e_1[2]*e_2[2]; // g_12
        double g_13 = e_1[0]*e_3[0] + e_1[1]*e_3[1] + e_1[2]*e_3[2]; // g_13
        double g_22 = e_2[0]*e_2[0] + e_2[1]*e_2[1] + e_2[2]*e_2[2]; // g_22
        double g_23 = e_2[0]*e_3[0] + e_2[1]*e_3[1] + e_2[2]*e_3[2]; // g_23
        double g_33 = e_3[0]*e_3[0] + e_3[1]*e_3[1] + e_3[2]*e_3[2]; // g_33
        
        g_ij_data[idx*6 + 0] = g_11;
        g_ij_data[idx*6 + 1] = g_12;
        g_ij_data[idx*6 + 2] = g_13;
        g_ij_data[idx*6 + 3] = g_22;
        g_ij_data[idx*6 + 4] = g_23;
        g_ij_data[idx*6 + 5] = g_33;
        
        
        // Normal vectors: n^i = e^i / |e^i|
        double mage1 = sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
        double mage2 = sqrt(e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2]);
        double mage3 = sqrt(e3[0]*e3[0] + e3[1]*e3[1] + e3[2]*e3[2]);
        
        if (mage1 > 1e-12) {
            for (int i=0; i<3; i++) normals_data[idx*9+0*3+i] = e1[i] / mage1;
        } else {
            for (int i=0; i<3; i++) normals_data[idx*9+0*3+i] = 0.0;
        }
        if (mage2 > 1e-12) {
            for (int i=0; i<3; i++) normals_data[idx*9+1*3+i] = e2[i] / mage2;
        } else {
            for (int i=0; i<3; i++) normals_data[idx*9+1*3+i] = 0.0;
        }
        if (mage3 > 1e-12) {
            for (int i=0; i<3; i++) normals_data[idx*9+2*3+i] = e3[i] / mage3;
            } else {
            for (int i=0; i<3; i++) normals_data[idx*9+2*3+i] = 0.0;
        }
        
        // Contravariant metric tensor: g^ij = (g_ij)^(-1)
        // Compute cofactors first (without 1/det factor)
        double gr11 = g_22*g_33 - g_23*g_23;  // Note: g_23 = g_32 for symmetric matrix
        double gr12 = g_13*g_23 - g_12*g_33;
        double gr13 = g_12*g_23 - g_13*g_22;
        double gr22 = g_11*g_33 - g_13*g_13;
        double gr23 = g_12*g_13 - g_11*g_23;
        double gr33 = g_11*g_22 - g_12*g_12;
        
        // Determinant of metric tensor: det(g) = J^2
        double det_g = g_11*gr11 + g_12*gr12 + g_13*gr13;
        if (fabs(det_g) < 1e-12) det_g = 1e-12;
        
        // Divide by determinant to get contravariant metric
        double g11 = gr11 / det_g;
        double g12 = gr12 / det_g;
        double g13 = gr13 / det_g;
        double g22 = gr22 / det_g;
        double g23 = gr23 / det_g;
        double g33 = gr33 / det_g;
        
        gij_data[idx*6 + 0] = g11;
        gij_data[idx*6 + 1] = g12;
        gij_data[idx*6 + 2] = g13;
        gij_data[idx*6 + 3] = g22;
        gij_data[idx*6 + 4] = g23;
        gij_data[idx*6 + 5] = g33;
        
        // rtg33inv = 1/sqrt(g_33)
        double sqrt_g33 = sqrt(g_33);
        if (sqrt_g33 < 1e-12) sqrt_g33 = 1e-12;
        rtg33inv_data[idx] = 1.0 / sqrt_g33;
        
        // b_i = g_{i3} / sqrt(g_33) (covariant components of magnetic unit vector)
        b_i_data[idx*3 + 0] = g_13 / sqrt_g33;  // b_1
        b_i_data[idx*3 + 1] = g_23 / sqrt_g33;  // b_2
        b_i_data[idx*3 + 2] = sqrt_g33;         // b_3 = sqrt(g_33)
        
        // B3 = B / sqrt(g_33)
        B3_data[idx] = bmag / sqrt_g33;
        
        // dualcurlbhat = e^i · curl(b̂)
        dualcurlbhat_data[idx*3 + 0] = e1[0]*curlb[0] + e1[1]*curlb[1] + e1[2]*curlb[2];  // e^1 · curl(b̂)
        dualcurlbhat_data[idx*3 + 1] = e2[0]*curlb[0] + e2[1]*curlb[1] + e2[2]*curlb[2];  // e^2 · curl(b̂)
        dualcurlbhat_data[idx*3 + 2] = e3[0]*curlb[0] + e3[1]*curlb[1] + e3[2]*curlb[2];  // e^3 · curl(b̂)
        
        // dualcurlbhatoverB = (e^i · curl(b̂)) / B
        if (bmag > 1e-12) {
            dualcurlbhatoverB_data[idx*3 + 0] = dualcurlbhat_data[idx*3 + 0] / bmag;
            dualcurlbhatoverB_data[idx*3 + 1] = dualcurlbhat_data[idx*3 + 1] / bmag;
            dualcurlbhatoverB_data[idx*3 + 2] = dualcurlbhat_data[idx*3 + 2] / bmag;
        } else {
            dualcurlbhatoverB_data[idx*3 + 0] = 0.0;
            dualcurlbhatoverB_data[idx*3 + 1] = 0.0;
            dualcurlbhatoverB_data[idx*3 + 2] = 0.0;
        }
        
        // bioverJB = b_i / (J * B)
        double JB = J * bmag;
        if (JB > 1e-12) {
            bioverJB_data[idx*3 + 0] = b_i_data[idx*3 + 0] / JB;
            bioverJB_data[idx*3 + 1] = b_i_data[idx*3 + 1] / JB;
            bioverJB_data[idx*3 + 2] = b_i_data[idx*3 + 2] / JB;
        } else {
            bioverJB_data[idx*3 + 0] = 0.0;
            bioverJB_data[idx*3 + 1] = 0.0;
            bioverJB_data[idx*3 + 2] = 0.0;
        }
        
        // cmag = J * B / sqrt(g_33)
        cmag_data[idx] = JB / sqrt_g33;
        
        // jacobtot = J * B
        jacobtot_data[idx] = JB;
        
        // jacobtot_inv = 1/(J * B)
        if (JB > 1e-12) {
            jacobtot_inv_data[idx] = 1.0 / JB;
        } else {
            jacobtot_inv_data[idx] = 1.0 / 1e-12;
        }
        
        // bmag_inv = 1/B
        if (bmag > 1e-12) {
            bmag_inv_data[idx] = 1.0 / bmag;
            } else {
            bmag_inv_data[idx] = 1.0 / 1e-12;
        }
        
        // bmag_inv_sq = (1/B)^2
        double bmag_inv_val = bmag_inv_data[idx];
        bmag_inv_sq_data[idx] = bmag_inv_val * bmag_inv_val;
        
        // Poisson solver quantities: g^ij * J
        gxxj_data[idx] = g11 * J;  // g^11 * J
        gxyj_data[idx] = g12 * J;  // g^12 * J
        gyyj_data[idx] = g22 * J;  // g^22 * J
        gxzj_data[idx] = g13 * J;  // g^13 * J
        
        // eps2 = J*g^33 - J/g_33
        eps2_data[idx] = J * g33 - J / g_33;

    }


    // Do again for surfaces
    struct gkyl_array* dxdz_surface1_nodal = load_gkyl_file_to_nodal("drzp_drza_S1.gkyl");
    struct gkyl_array* bmag_surface1_nodal = load_gkyl_file_to_nodal("Bmag_S1.gkyl");
    struct gkyl_array* curlbhat_surface1_nodal = load_gkyl_file_to_nodal("curl_B_hat_xyz_S1.gkyl");
    struct gkyl_array* dxdz_surface2_nodal = load_gkyl_file_to_nodal("drzp_drza_S2.gkyl");
    struct gkyl_array* bmag_surface2_nodal = load_gkyl_file_to_nodal("Bmag_S2.gkyl");
    struct gkyl_array* curlbhat_surface2_nodal = load_gkyl_file_to_nodal("curl_B_hat_xyz_S2.gkyl");
    struct gkyl_array* dxdz_surface3_nodal = load_gkyl_file_to_nodal("drzp_drza_S3.gkyl");
    struct gkyl_array* bmag_surface3_nodal = load_gkyl_file_to_nodal("Bmag_S3.gkyl");
    struct gkyl_array* curlbhat_surface3_nodal = load_gkyl_file_to_nodal("curl_B_hat_xyz_S3.gkyl");
    
    // Check all surface arrays loaded successfully
    TEST_CHECK(dxdz_surface1_nodal != NULL && bmag_surface1_nodal != NULL && curlbhat_surface1_nodal != NULL);
    TEST_CHECK(dxdz_surface2_nodal != NULL && bmag_surface2_nodal != NULL && curlbhat_surface2_nodal != NULL);
    TEST_CHECK(dxdz_surface3_nodal != NULL && bmag_surface3_nodal != NULL && curlbhat_surface3_nodal != NULL);
    if (!dxdz_surface1_nodal || !bmag_surface1_nodal || !curlbhat_surface1_nodal ||
        !dxdz_surface2_nodal || !bmag_surface2_nodal || !curlbhat_surface2_nodal ||
        !dxdz_surface3_nodal || !bmag_surface3_nodal || !curlbhat_surface3_nodal) {
        return;
    }
    // Derived quantities
    // Determine number of nodes from loaded array size (bmag has 1 component, so size = number of nodes)
    // All surfaces have the same number of nodes, so use surface1
    int shape_surface[1] = {(int)bmag_surface1_nodal->size};
    struct gkyl_range nrange_surface;
    gkyl_range_init_from_shape1(&nrange_surface, 1, shape_surface);
    struct gkyl_range_iter iter_surface;
    gkyl_range_iter_init(&iter_surface, &nrange_surface);  

    // calculate quantities for surfaces
    struct gkyl_array* B3_surface1_nodal = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* b_i_surface1_nodal = mkarr(false, 3, nrange_surface.volume);
    struct gkyl_array* cmag_surface1_nodal = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* normcurlbhat_surface1_nodal = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* jacobgeo_surface1_nodal = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* jacobtot_inv_surface1_nodal = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* normals_surface1_nodal = mkarr(false, 9, nrange_surface.volume);
    struct gkyl_array* lenr_surface1_nodal = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* B3_surface2_nodal = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* b_i_surface2_nodal = mkarr(false, 3, nrange_surface.volume);
    struct gkyl_array* cmag_surface2_nodal = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* normcurlbhat_surface2_nodal = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* jacobgeo_surface2_nodal = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* jacobtot_inv_surface2_nodal = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* normals_surface2_nodal = mkarr(false, 9, nrange_surface.volume);
    struct gkyl_array* lenr_surface2_nodal = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* B3_surface3_nodal = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* b_i_surface3_nodal = mkarr(false, 3, nrange_surface.volume);
    struct gkyl_array* cmag_surface3_nodal = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* normcurlbhat_surface3_nodal = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* jacobgeo_surface3_nodal = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* jacobtot_inv_surface3_nodal = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* normals_surface3_nodal = mkarr(false, 9, nrange_surface.volume);
    struct gkyl_array* lenr_surface3_nodal = mkarr(false, 1, nrange_surface.volume);
    
    while (gkyl_range_iter_next(&iter_surface)) {
        // Get current node index
        long idx = gkyl_range_idx(&nrange_surface, iter_surface.idx);
        
        double *dxdz_data1 = (double*)dxdz_surface1_nodal->data;
        double *bmag_data1 = (double*)bmag_surface1_nodal->data;
        double *curlbhat_data1 = (double*)curlbhat_surface1_nodal->data;
        double *B3_data1 = (double*)B3_surface1_nodal->data;
        double *b_i_data1 = (double*)b_i_surface1_nodal->data;
        double *cmag_data1 = (double*)cmag_surface1_nodal->data;
        double *normcurlbhat_data1 = (double*)normcurlbhat_surface1_nodal->data;
        double *jacobgeo_data1 = (double*)jacobgeo_surface1_nodal->data;
        double *jacobtot_inv_data1 = (double*)jacobtot_inv_surface1_nodal->data;
        double *normals_data1 = (double*)normals_surface1_nodal->data;
        double *lenr_data1 = (double*)lenr_surface1_nodal->data;

        double *dxdz_data2 = (double*)dxdz_surface2_nodal->data;
        double *bmag_data2 = (double*)bmag_surface2_nodal->data;
        double *curlbhat_data2 = (double*)curlbhat_surface2_nodal->data;
        double *B3_data2 = (double*)B3_surface2_nodal->data;
        double *b_i_data2 = (double*)b_i_surface2_nodal->data;
        double *cmag_data2 = (double*)cmag_surface2_nodal->data;
        double *normcurlbhat_data2 = (double*)normcurlbhat_surface2_nodal->data;
        double *jacobgeo_data2 = (double*)jacobgeo_surface2_nodal->data;
        double *jacobtot_inv_data2 = (double*)jacobtot_inv_surface2_nodal->data;
        double *normals_data2 = (double*)normals_surface2_nodal->data;
        double *lenr_data2 = (double*)lenr_surface2_nodal->data;

        double *dxdz_data3 = (double*)dxdz_surface3_nodal->data;
        double *bmag_data3 = (double*)bmag_surface3_nodal->data;
        double *curlbhat_data3 = (double*)curlbhat_surface3_nodal->data;
        double *B3_data3 = (double*)B3_surface3_nodal->data;
        double *b_i_data3 = (double*)b_i_surface3_nodal->data;
        double *cmag_data3 = (double*)cmag_surface3_nodal->data;
        double *normcurlbhat_data3 = (double*)normcurlbhat_surface3_nodal->data;
        double *jacobgeo_data3 = (double*)jacobgeo_surface3_nodal->data;
        double *jacobtot_inv_data3 = (double*)jacobtot_inv_surface3_nodal->data;
        double *normals_data3 = (double*)normals_surface3_nodal->data;
        double *lenr_data3 = (double*)lenr_surface3_nodal->data;
        
        // ========== SURFACE 1 (dir=1, normal to rho direction) ==========
        double e_1_s1[3] = {dxdz_data1[idx*9+0], dxdz_data1[idx*9+1], dxdz_data1[idx*9+2]};
        double e_2_s1[3] = {dxdz_data1[idx*9+3], dxdz_data1[idx*9+4], dxdz_data1[idx*9+5]};
        double e_3_s1[3] = {dxdz_data1[idx*9+6], dxdz_data1[idx*9+7], dxdz_data1[idx*9+8]};
        double curlb_s1[3] = {curlbhat_data1[idx*3+0], curlbhat_data1[idx*3+1], curlbhat_data1[idx*3+2]};
        double bmag_s1 = bmag_data1[idx];
        
        // Jacobian: J = e_1 · (e_2 × e_3)
        double e2xe3_s1[3] = {e_2_s1[1]*e_3_s1[2]-e_2_s1[2]*e_3_s1[1], e_2_s1[2]*e_3_s1[0]-e_2_s1[0]*e_3_s1[2], e_2_s1[0]*e_3_s1[1]-e_2_s1[1]*e_3_s1[0]};
        double J_s1 = e_1_s1[0]*e2xe3_s1[0] + e_1_s1[1]*e2xe3_s1[1] + e_1_s1[2]*e2xe3_s1[2];
        if (fabs(J_s1) < 1e-12) J_s1 = 1e-12;
        jacobgeo_data1[idx] = J_s1;
        double Jinv_s1 = 1.0/J_s1;
        
        // Dual vectors: e^i = (1/J) * (e_j × e_k)
        double e3xe1_s1[3] = {e_3_s1[1]*e_1_s1[2]-e_3_s1[2]*e_1_s1[1], e_3_s1[2]*e_1_s1[0]-e_3_s1[0]*e_1_s1[2], e_3_s1[0]*e_1_s1[1]-e_3_s1[1]*e_1_s1[0]};
        double e1xe2_s1[3] = {e_1_s1[1]*e_2_s1[2]-e_1_s1[2]*e_2_s1[1], e_1_s1[2]*e_2_s1[0]-e_1_s1[0]*e_2_s1[2], e_1_s1[0]*e_2_s1[1]-e_1_s1[1]*e_2_s1[0]};
        double e1_s1[3] = {Jinv_s1 * e2xe3_s1[0], Jinv_s1 * e2xe3_s1[1], Jinv_s1 * e2xe3_s1[2]};  // e^1
        double e2_s1[3] = {Jinv_s1 * e3xe1_s1[0], Jinv_s1 * e3xe1_s1[1], Jinv_s1 * e3xe1_s1[2]};  // e^2
        double e3_s1[3] = {Jinv_s1 * e1xe2_s1[0], Jinv_s1 * e1xe2_s1[1], Jinv_s1 * e1xe2_s1[2]};  // e^3
        
        // Metric tensor components (only need g_33 for some calculations)
        double g_33_s1 = e_3_s1[0]*e_3_s1[0] + e_3_s1[1]*e_3_s1[1] + e_3_s1[2]*e_3_s1[2];
        double g_13_s1 = e_1_s1[0]*e_3_s1[0] + e_1_s1[1]*e_3_s1[1] + e_1_s1[2]*e_3_s1[2];
        double g_23_s1 = e_2_s1[0]*e_3_s1[0] + e_2_s1[1]*e_3_s1[1] + e_2_s1[2]*e_3_s1[2];
        
        double sqrt_g33_s1 = sqrt(g_33_s1);
        if (sqrt_g33_s1 < 1e-12) sqrt_g33_s1 = 1e-12;
        
        // Normal vectors: n^i = e^i / |e^i|
        double mage1_s1 = sqrt(e1_s1[0]*e1_s1[0] + e1_s1[1]*e1_s1[1] + e1_s1[2]*e1_s1[2]);
        double mage2_s1 = sqrt(e2_s1[0]*e2_s1[0] + e2_s1[1]*e2_s1[1] + e2_s1[2]*e2_s1[2]);
        double mage3_s1 = sqrt(e3_s1[0]*e3_s1[0] + e3_s1[1]*e3_s1[1] + e3_s1[2]*e3_s1[2]);
        
        if (mage1_s1 > 1e-12) {
            for (int i=0; i<3; i++) normals_data1[idx*9+0*3+i] = e1_s1[i] / mage1_s1;
        } else {
            for (int i=0; i<3; i++) normals_data1[idx*9+0*3+i] = 0.0;
        }
        if (mage2_s1 > 1e-12) {
            for (int i=0; i<3; i++) normals_data1[idx*9+1*3+i] = e2_s1[i] / mage2_s1;
        } else {
            for (int i=0; i<3; i++) normals_data1[idx*9+1*3+i] = 0.0;
        }
        if (mage3_s1 > 1e-12) {
            for (int i=0; i<3; i++) normals_data1[idx*9+2*3+i] = e3_s1[i] / mage3_s1;
        } else {
            for (int i=0; i<3; i++) normals_data1[idx*9+2*3+i] = 0.0;
        }
        
        // b_i = g_{i3} / sqrt(g_33)
        b_i_data1[idx*3 + 0] = g_13_s1 / sqrt_g33_s1;  // b_1
        b_i_data1[idx*3 + 1] = g_23_s1 / sqrt_g33_s1;  // b_2
        b_i_data1[idx*3 + 2] = sqrt_g33_s1;            // b_3
        
        // B3 = bmag / (sqrt(g_33) * |e^3|)
        B3_data1[idx] = bmag_s1 / (sqrt_g33_s1 * mage3_s1);
        
        // cmag = J * B / sqrt(g_33)
        double JB_s1 = J_s1 * bmag_s1;
        cmag_data1[idx] = JB_s1 / sqrt_g33_s1;
        
        // jacobtot_inv = 1/(J * B)
        if (JB_s1 > 1e-12) {
            jacobtot_inv_data1[idx] = 1.0 / JB_s1;
            } else {
            jacobtot_inv_data1[idx] = 1.0 / 1e-12;
        }
        
        // normcurlbhat = n^dir · curl(b̂) where dir=0 for surface 1
        normcurlbhat_data1[idx] = normals_data1[idx*9+0*3+0]*curlb_s1[0] + 
                                   normals_data1[idx*9+0*3+1]*curlb_s1[1] + 
                                   normals_data1[idx*9+0*3+2]*curlb_s1[2];
        
        // lenr = J * |e^dir| where dir=1 for surface 1
        lenr_data1[idx] = J_s1 * mage1_s1;
        
        // ========== SURFACE 2 (dir=2, normal to zeta direction) ==========
        double e_1_s2[3] = {dxdz_data2[idx*9+0], dxdz_data2[idx*9+1], dxdz_data2[idx*9+2]};
        double e_2_s2[3] = {dxdz_data2[idx*9+3], dxdz_data2[idx*9+4], dxdz_data2[idx*9+5]};
        double e_3_s2[3] = {dxdz_data2[idx*9+6], dxdz_data2[idx*9+7], dxdz_data2[idx*9+8]};
        double curlb_s2[3] = {curlbhat_data2[idx*3+0], curlbhat_data2[idx*3+1], curlbhat_data2[idx*3+2]};
        double bmag_s2 = bmag_data2[idx];
        
        // Jacobian: J = e_1 · (e_2 × e_3)
        double e2xe3_s2[3] = {e_2_s2[1]*e_3_s2[2]-e_2_s2[2]*e_3_s2[1], e_2_s2[2]*e_3_s2[0]-e_2_s2[0]*e_3_s2[2], e_2_s2[0]*e_3_s2[1]-e_2_s2[1]*e_3_s2[0]};
        double J_s2 = e_1_s2[0]*e2xe3_s2[0] + e_1_s2[1]*e2xe3_s2[1] + e_1_s2[2]*e2xe3_s2[2];
        if (fabs(J_s2) < 1e-12) J_s2 = 1e-12;
        jacobgeo_data2[idx] = J_s2;
        double Jinv_s2 = 1.0/J_s2;
        
        // Dual vectors: e^i = (1/J) * (e_j × e_k)
        double e3xe1_s2[3] = {e_3_s2[1]*e_1_s2[2]-e_3_s2[2]*e_1_s2[1], e_3_s2[2]*e_1_s2[0]-e_3_s2[0]*e_1_s2[2], e_3_s2[0]*e_1_s2[1]-e_3_s2[1]*e_1_s2[0]};
        double e1xe2_s2[3] = {e_1_s2[1]*e_2_s2[2]-e_1_s2[2]*e_2_s2[1], e_1_s2[2]*e_2_s2[0]-e_1_s2[0]*e_2_s2[2], e_1_s2[0]*e_2_s2[1]-e_1_s2[1]*e_2_s2[0]};
        double e1_s2[3] = {Jinv_s2 * e2xe3_s2[0], Jinv_s2 * e2xe3_s2[1], Jinv_s2 * e2xe3_s2[2]};  // e^1
        double e2_s2[3] = {Jinv_s2 * e3xe1_s2[0], Jinv_s2 * e3xe1_s2[1], Jinv_s2 * e3xe1_s2[2]};  // e^2
        double e3_s2[3] = {Jinv_s2 * e1xe2_s2[0], Jinv_s2 * e1xe2_s2[1], Jinv_s2 * e1xe2_s2[2]};  // e^3
        
        // Metric tensor components
        double g_33_s2 = e_3_s2[0]*e_3_s2[0] + e_3_s2[1]*e_3_s2[1] + e_3_s2[2]*e_3_s2[2];
        double g_13_s2 = e_1_s2[0]*e_3_s2[0] + e_1_s2[1]*e_3_s2[1] + e_1_s2[2]*e_3_s2[2];
        double g_23_s2 = e_2_s2[0]*e_3_s2[0] + e_2_s2[1]*e_3_s2[1] + e_2_s2[2]*e_3_s2[2];
        
        double sqrt_g33_s2 = sqrt(g_33_s2);
        if (sqrt_g33_s2 < 1e-12) sqrt_g33_s2 = 1e-12;
        
        // Normal vectors: n^i = e^i / |e^i|
        double mage1_s2 = sqrt(e1_s2[0]*e1_s2[0] + e1_s2[1]*e1_s2[1] + e1_s2[2]*e1_s2[2]);
        double mage2_s2 = sqrt(e2_s2[0]*e2_s2[0] + e2_s2[1]*e2_s2[1] + e2_s2[2]*e2_s2[2]);
        double mage3_s2 = sqrt(e3_s2[0]*e3_s2[0] + e3_s2[1]*e3_s2[1] + e3_s2[2]*e3_s2[2]);
        
        if (mage1_s2 > 1e-12) {
            for (int i=0; i<3; i++) normals_data2[idx*9+0*3+i] = e1_s2[i] / mage1_s2;
        } else {
            for (int i=0; i<3; i++) normals_data2[idx*9+0*3+i] = 0.0;
        }
        if (mage2_s2 > 1e-12) {
            for (int i=0; i<3; i++) normals_data2[idx*9+1*3+i] = e2_s2[i] / mage2_s2;
        } else {
            for (int i=0; i<3; i++) normals_data2[idx*9+1*3+i] = 0.0;
        }
        if (mage3_s2 > 1e-12) {
            for (int i=0; i<3; i++) normals_data2[idx*9+2*3+i] = e3_s2[i] / mage3_s2;
        } else {
            for (int i=0; i<3; i++) normals_data2[idx*9+2*3+i] = 0.0;
        }
        
        // b_i = g_{i3} / sqrt(g_33)
        b_i_data2[idx*3 + 0] = g_13_s2 / sqrt_g33_s2;  // b_1
        b_i_data2[idx*3 + 1] = g_23_s2 / sqrt_g33_s2;  // b_2
        b_i_data2[idx*3 + 2] = sqrt_g33_s2;            // b_3
        
        // B3 = bmag / (sqrt(g_33) * |e^3|)
        B3_data2[idx] = bmag_s2 / (sqrt_g33_s2 * mage3_s2);
        
        // cmag = J * B / sqrt(g_33)
        double JB_s2 = J_s2 * bmag_s2;
        cmag_data2[idx] = JB_s2 / sqrt_g33_s2;
        
        // jacobtot_inv = 1/(J * B)
        if (JB_s2 > 1e-12) {
            jacobtot_inv_data2[idx] = 1.0 / JB_s2;
            } else {
            jacobtot_inv_data2[idx] = 1.0 / 1e-12;
        }
        
        // normcurlbhat = n^dir · curl(b̂) where dir=1 for surface 2
        normcurlbhat_data2[idx] = normals_data2[idx*9+1*3+0]*curlb_s2[0] + 
                                   normals_data2[idx*9+1*3+1]*curlb_s2[1] + 
                                   normals_data2[idx*9+1*3+2]*curlb_s2[2];
        
        // lenr = J * |e^dir| where dir=2 for surface 2
        lenr_data2[idx] = J_s2 * mage2_s2;
        
        // ========== SURFACE 3 (dir=3, normal to alpha direction) ==========
        double e_1_s3[3] = {dxdz_data3[idx*9+0], dxdz_data3[idx*9+1], dxdz_data3[idx*9+2]};
        double e_2_s3[3] = {dxdz_data3[idx*9+3], dxdz_data3[idx*9+4], dxdz_data3[idx*9+5]};
        double e_3_s3[3] = {dxdz_data3[idx*9+6], dxdz_data3[idx*9+7], dxdz_data3[idx*9+8]};
        double curlb_s3[3] = {curlbhat_data3[idx*3+0], curlbhat_data3[idx*3+1], curlbhat_data3[idx*3+2]};
        double bmag_s3 = bmag_data3[idx];
        
        // Jacobian: J = e_1 · (e_2 × e_3)
        double e2xe3_s3[3] = {e_2_s3[1]*e_3_s3[2]-e_2_s3[2]*e_3_s3[1], e_2_s3[2]*e_3_s3[0]-e_2_s3[0]*e_3_s3[2], e_2_s3[0]*e_3_s3[1]-e_2_s3[1]*e_3_s3[0]};
        double J_s3 = e_1_s3[0]*e2xe3_s3[0] + e_1_s3[1]*e2xe3_s3[1] + e_1_s3[2]*e2xe3_s3[2];
        if (fabs(J_s3) < 1e-12) J_s3 = 1e-12;
        jacobgeo_data3[idx] = J_s3;
        double Jinv_s3 = 1.0/J_s3;
        
        // Dual vectors: e^i = (1/J) * (e_j × e_k)
        double e3xe1_s3[3] = {e_3_s3[1]*e_1_s3[2]-e_3_s3[2]*e_1_s3[1], e_3_s3[2]*e_1_s3[0]-e_3_s3[0]*e_1_s3[2], e_3_s3[0]*e_1_s3[1]-e_3_s3[1]*e_1_s3[0]};
        double e1xe2_s3[3] = {e_1_s3[1]*e_2_s3[2]-e_1_s3[2]*e_2_s3[1], e_1_s3[2]*e_2_s3[0]-e_1_s3[0]*e_2_s3[2], e_1_s3[0]*e_2_s3[1]-e_1_s3[1]*e_2_s3[0]};
        double e1_s3[3] = {Jinv_s3 * e2xe3_s3[0], Jinv_s3 * e2xe3_s3[1], Jinv_s3 * e2xe3_s3[2]};  // e^1
        double e2_s3[3] = {Jinv_s3 * e3xe1_s3[0], Jinv_s3 * e3xe1_s3[1], Jinv_s3 * e3xe1_s3[2]};  // e^2
        double e3_s3[3] = {Jinv_s3 * e1xe2_s3[0], Jinv_s3 * e1xe2_s3[1], Jinv_s3 * e1xe2_s3[2]};  // e^3
        
        // Metric tensor components
        double g_33_s3 = e_3_s3[0]*e_3_s3[0] + e_3_s3[1]*e_3_s3[1] + e_3_s3[2]*e_3_s3[2];
        double g_13_s3 = e_1_s3[0]*e_3_s3[0] + e_1_s3[1]*e_3_s3[1] + e_1_s3[2]*e_3_s3[2];
        double g_23_s3 = e_2_s3[0]*e_3_s3[0] + e_2_s3[1]*e_3_s3[1] + e_2_s3[2]*e_3_s3[2];
        
        double sqrt_g33_s3 = sqrt(g_33_s3);
        if (sqrt_g33_s3 < 1e-12) sqrt_g33_s3 = 1e-12;
        
        // Normal vectors: n^i = e^i / |e^i|
        double mage1_s3 = sqrt(e1_s3[0]*e1_s3[0] + e1_s3[1]*e1_s3[1] + e1_s3[2]*e1_s3[2]);
        double mage2_s3 = sqrt(e2_s3[0]*e2_s3[0] + e2_s3[1]*e2_s3[1] + e2_s3[2]*e2_s3[2]);
        double mage3_s3 = sqrt(e3_s3[0]*e3_s3[0] + e3_s3[1]*e3_s3[1] + e3_s3[2]*e3_s3[2]);
        
        if (mage1_s3 > 1e-12) {
            for (int i=0; i<3; i++) normals_data3[idx*9+0*3+i] = e1_s3[i] / mage1_s3;
        } else {
            for (int i=0; i<3; i++) normals_data3[idx*9+0*3+i] = 0.0;
        }
        if (mage2_s3 > 1e-12) {
            for (int i=0; i<3; i++) normals_data3[idx*9+1*3+i] = e2_s3[i] / mage2_s3;
        } else {
            for (int i=0; i<3; i++) normals_data3[idx*9+1*3+i] = 0.0;
        }
        if (mage3_s3 > 1e-12) {
            for (int i=0; i<3; i++) normals_data3[idx*9+2*3+i] = e3_s3[i] / mage3_s3;
        } else {
            for (int i=0; i<3; i++) normals_data3[idx*9+2*3+i] = 0.0;
        }
        
        // b_i = g_{i3} / sqrt(g_33)
        b_i_data3[idx*3 + 0] = g_13_s3 / sqrt_g33_s3;  // b_1
        b_i_data3[idx*3 + 1] = g_23_s3 / sqrt_g33_s3;  // b_2
        b_i_data3[idx*3 + 2] = sqrt_g33_s3;            // b_3
        
        // B3 = bmag / (sqrt(g_33) * |e^3|)
        B3_data3[idx] = bmag_s3 / (sqrt_g33_s3 * mage3_s3);
        
        // cmag = J * B / sqrt(g_33)
        double JB_s3 = J_s3 * bmag_s3;
        cmag_data3[idx] = JB_s3 / sqrt_g33_s3;
        
        // jacobtot_inv = 1/(J * B)
        if (JB_s3 > 1e-12) {
            jacobtot_inv_data3[idx] = 1.0 / JB_s3;
            } else {
            jacobtot_inv_data3[idx] = 1.0 / 1e-12;
        }
        
        // normcurlbhat = n^dir · curl(b̂) where dir=2 for surface 3
        normcurlbhat_data3[idx] = normals_data3[idx*9+2*3+0]*curlb_s3[0] + 
                                   normals_data3[idx*9+2*3+1]*curlb_s3[1] + 
                                   normals_data3[idx*9+2*3+2]*curlb_s3[2];
        
        // lenr = J * |e^dir| where dir=3 for surface 3
        lenr_data3[idx] = J_s3 * mage3_s3;
        
    }


    struct gkyl_rect_grid grid_orig;  // Original grid with proper bounds for writing
    struct gkyl_rect_grid grid;       // Grid with dummy bounds for n2m conversion
    struct gkyl_array_header_info hdr;
    // Use one of the files we already successfully loaded (e.g., bmag_corner_nodal)
    // All files should have the same grid structure
    enum gkyl_array_rio_status status = gkyl_grid_sub_array_header_read(&grid_orig, &hdr, "W7-X_field-aligned_coords_gkyl/Bmag_corner.gkyl");

    
    // Copy original grid and set dummy bounds for n2m conversion
    grid = grid_orig;
    for (int d=0; d<grid.ndim; d++) {
        grid.lower[d] = 0.0;
        grid.upper[d] = 1.0;
        grid.dx[d] = 1.0 / grid.cells[d];
    }
    grid.cellVolume = 1.0;
    
    // Create basis (3D, poly_order=1)
    struct gkyl_basis basis;
    gkyl_cart_modal_serendip(&basis, 3, 1);
    // Create surface basis: for 3D grid, surface basis is 2D with same poly_order
    struct gkyl_basis surf_basis;
    gkyl_cart_modal_serendip(&surf_basis, 2, 1);
    
    // Create computational ranges (all use same grid cells)
    struct gkyl_range comp_range;
    gkyl_range_init_from_shape1(&comp_range, grid.ndim, grid.cells);
    
    // Create n2m operators
    struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &grid, false);
    
    // Helper macro to convert nodal to modal
    #define N2M_CONVERT(nrange, comp_range, ncomp, nodal, modal) \
        modal = mkarr(false, ncomp, comp_range.volume * basis.num_basis); \
        gkyl_nodal_ops_n2m(n2m, &basis, &grid, &nrange, &comp_range, ncomp, nodal, modal, false);
    
    #define N2M_CONVERT_INTERIOR(nrange, comp_range, ncomp, nodal, modal) \
        modal = mkarr(false, ncomp, comp_range.volume * basis.num_basis); \
        gkyl_nodal_ops_n2m(n2m, &basis, &grid, &nrange, &comp_range, ncomp, nodal, modal, true);
    
    #define N2M_CONVERT_SURFACE(nrange, comp_range, ncomp, nodal, modal, dir) \
        modal = mkarr(false, ncomp, comp_range.volume * surf_basis.num_basis); \
        gkyl_nodal_ops_n2m_surface(n2m, &surf_basis, &grid, &nrange, &comp_range, ncomp, nodal, modal, dir);
    
    // Convert corner arrays - use_quad=false for corners
    struct gkyl_array *mc2nu_pos_corner_modal, *mc2p_corner_modal, *bmag_corner_modal;
    N2M_CONVERT(nrange_corner, comp_range, 3, mc2nu_pos_corner_nodal, mc2nu_pos_corner_modal);
    N2M_CONVERT(nrange_corner, comp_range, 3, mc2p_corner_nodal, mc2p_corner_modal);
    N2M_CONVERT(nrange_corner, comp_range, 1, bmag_corner_nodal, bmag_corner_modal);
    
    // Convert interior arrays - use_quad=true for interior
    struct gkyl_array *dxdz_interior_modal, *bmag_interior_modal, *curlbhat_interior_modal, *bcart_interior_modal;
    struct gkyl_array *B3_interior_modal, *b_i_interior_modal, *bmag_inv_interior_modal, *bmag_inv_sq_interior_modal;
    struct gkyl_array *bioverJB_interior_modal, *cmag_interior_modal, *dualcurlbhat_interior_modal, *dualcurlbhatoverB_interior_modal;
    struct gkyl_array *dzdx_interior_modal, *eps2_interior_modal, *g_ij_interior_modal, *gij_interior_modal;
    struct gkyl_array *gxxj_interior_modal, *gxyj_interior_modal, *gxzj_interior_modal, *gyyj_interior_modal;
    struct gkyl_array *jacobgeo_interior_modal, *jacobgeo_inv_interior_modal, *jacobtot_interior_modal, *jacobtot_inv_interior_modal;
    struct gkyl_array *normals_interior_modal, *rtg33inv_interior_modal;
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 9, dxdz_interior_nodal, dxdz_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 1, bmag_interior_nodal, bmag_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 3, curlbhat_interior_nodal, curlbhat_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 3, bcart_interior_nodal, bcart_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 1, B3_interior_nodal, B3_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 3, b_i_interior_nodal, b_i_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 1, bmag_inv_interior_nodal, bmag_inv_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 1, bmag_inv_sq_interior_nodal, bmag_inv_sq_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 3, bioverJB_interior_nodal, bioverJB_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 1, cmag_interior_nodal, cmag_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 3, dualcurlbhat_interior_nodal, dualcurlbhat_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 3, dualcurlbhatoverB_interior_nodal, dualcurlbhatoverB_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 9, dzdx_interior_nodal, dzdx_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 1, eps2_interior_nodal, eps2_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 6, g_ij_interior_nodal, g_ij_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 6, gij_interior_nodal, gij_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 1, gxxj_interior_nodal, gxxj_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 1, gxyj_interior_nodal, gxyj_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 1, gxzj_interior_nodal, gxzj_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 1, gyyj_interior_nodal, gyyj_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 1, jacobgeo_interior_nodal, jacobgeo_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 1, jacobgeo_inv_interior_nodal, jacobgeo_inv_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 1, jacobtot_interior_nodal, jacobtot_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 1, jacobtot_inv_interior_nodal, jacobtot_inv_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 9, normals_interior_nodal, normals_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 1, rtg33inv_interior_nodal, rtg33inv_interior_modal);
    
    // Convert surface arrays - use_quad=false for surfaces
    // Surface 1 arrays - use surface n2m with direction 0
    struct gkyl_array *dxdz_s1_modal, *bmag_s1_modal, *curlbhat_s1_modal;
    struct gkyl_array *B3_s1_modal, *b_i_s1_modal, *cmag_s1_modal, *normcurlbhat_s1_modal;
    struct gkyl_array *jacobgeo_s1_modal, *jacobtot_inv_s1_modal, *normals_s1_modal, *lenr_s1_modal;
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 9, dxdz_surface1_nodal, dxdz_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, bmag_surface1_nodal, bmag_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, curlbhat_surface1_nodal, curlbhat_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, B3_surface1_nodal, B3_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, b_i_surface1_nodal, b_i_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, cmag_surface1_nodal, cmag_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, normcurlbhat_surface1_nodal, normcurlbhat_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, jacobgeo_surface1_nodal, jacobgeo_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, jacobtot_inv_surface1_nodal, jacobtot_inv_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 9, normals_surface1_nodal, normals_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, lenr_surface1_nodal, lenr_s1_modal, 0);
    
    // Surface 2 arrays - use surface n2m with direction 1
    struct gkyl_array *dxdz_s2_modal, *bmag_s2_modal, *curlbhat_s2_modal;
    struct gkyl_array *B3_s2_modal, *b_i_s2_modal, *cmag_s2_modal, *normcurlbhat_s2_modal;
    struct gkyl_array *jacobgeo_s2_modal, *jacobtot_inv_s2_modal, *normals_s2_modal, *lenr_s2_modal;
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 9, dxdz_surface2_nodal, dxdz_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, bmag_surface2_nodal, bmag_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, curlbhat_surface2_nodal, curlbhat_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, B3_surface2_nodal, B3_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, b_i_surface2_nodal, b_i_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, cmag_surface2_nodal, cmag_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, normcurlbhat_surface2_nodal, normcurlbhat_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, jacobgeo_surface2_nodal, jacobgeo_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, jacobtot_inv_surface2_nodal, jacobtot_inv_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 9, normals_surface2_nodal, normals_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, lenr_surface2_nodal, lenr_s2_modal, 1);
    
    // Surface 3 arrays - use surface n2m with direction 2
    struct gkyl_array *dxdz_s3_modal, *bmag_s3_modal, *curlbhat_s3_modal;
    struct gkyl_array *B3_s3_modal, *b_i_s3_modal, *cmag_s3_modal, *normcurlbhat_s3_modal;
    struct gkyl_array *jacobgeo_s3_modal, *jacobtot_inv_s3_modal, *normals_s3_modal, *lenr_s3_modal;
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 9, dxdz_surface3_nodal, dxdz_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, bmag_surface3_nodal, bmag_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, curlbhat_surface3_nodal, curlbhat_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, B3_surface3_nodal, B3_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, b_i_surface3_nodal, b_i_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, cmag_surface3_nodal, cmag_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, normcurlbhat_surface3_nodal, normcurlbhat_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, jacobgeo_surface3_nodal, jacobgeo_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, jacobtot_inv_surface3_nodal, jacobtot_inv_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 9, normals_surface3_nodal, normals_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, lenr_surface3_nodal, lenr_s3_modal, 2);
    
    // Debug: Print to both stdout and stderr to ensure we see output

    // Write out all nodal and modal arrays 
    // Create output directories
    // Note: Test is run from gkeyll/ directory, so create directories there
    const char *output_dir = "W7X-nodal_modal_arrays";
    const char *nodal_dir = "W7X-nodal_modal_arrays/nodal";
    const char *modal_dir = "W7X-nodal_modal_arrays/modal";
    
    // Create directories (must create parent first)
    #ifdef _WIN32
        if (_mkdir(output_dir) != 0 && errno != EEXIST) {
            return; // Exit if we can't create the main directory
        }
        if (_mkdir(nodal_dir) != 0 && errno != EEXIST) {
            return;
        }
        if (_mkdir(modal_dir) != 0 && errno != EEXIST) {
            return;
        }
    #else
        if (mkdir(output_dir, 0755) != 0 && errno != EEXIST) {
            return; // Exit if we can't create the main directory
        }
        if (mkdir(nodal_dir, 0755) != 0 && errno != EEXIST) {
            return;
        }
        if (mkdir(modal_dir, 0755) != 0 && errno != EEXIST) {
            return;
        }
    #endif
    
    // Helper macro to write nodal array
    #define WRITE_NODAL(ngrid, nrange, arr, name) \
        do { \
            if (arr != NULL) { \
                char fname[512]; \
                snprintf(fname, sizeof(fname), "%s/%s.gkyl", nodal_dir, name); \
                gkyl_grid_sub_array_write(ngrid, nrange, 0, arr, fname); \
            } \
        } while(0)
    
    // Helper macro to write modal array
    #define WRITE_MODAL(grid, comp_range, arr, name) \
        do { \
            if (arr != NULL) { \
                char fname[512]; \
                snprintf(fname, sizeof(fname), "%s/%s.gkyl", modal_dir, name); \
                gkyl_grid_sub_array_write(grid, comp_range, 0, arr, fname); \
            } \
        } while(0)
    
    // Create nodal grids for writing nodal arrays
    struct gkyl_rect_grid ngrid_corner, ngrid_interior, ngrid_surface;
    gkyl_gk_geometry_init_nodal_grid(&ngrid_corner, &grid_orig, &nrange_corner);
    gkyl_gk_geometry_init_nodal_grid(&ngrid_interior, &grid_orig, &nrange_interior);
    gkyl_gk_geometry_init_nodal_grid(&ngrid_surface, &grid_orig, &nrange_surface);
    
    // Write corner nodal arrays 
    WRITE_NODAL(&ngrid_corner, &nrange_corner, mc2nu_pos_corner_nodal, "mc2nu_pos_corner_nodal");
    WRITE_NODAL(&ngrid_corner, &nrange_corner, mc2p_corner_nodal, "mc2p_corner_nodal");
    WRITE_NODAL(&ngrid_corner, &nrange_corner, bmag_corner_nodal, "bmag_corner_nodal");
    
    // Write interior nodal arrays
    WRITE_NODAL(&ngrid_interior, &nrange_interior, dxdz_interior_nodal, "dxdz_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, bmag_interior_nodal, "bmag_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, curlbhat_interior_nodal, "curlbhat_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, bcart_interior_nodal, "bcart_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, B3_interior_nodal, "B3_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, b_i_interior_nodal, "b_i_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, bmag_inv_interior_nodal, "bmag_inv_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, bmag_inv_sq_interior_nodal, "bmag_inv_sq_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, bioverJB_interior_nodal, "bioverJB_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, cmag_interior_nodal, "cmag_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, dualcurlbhat_interior_nodal, "dualcurlbhat_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, dualcurlbhatoverB_interior_nodal, "dualcurlbhatoverB_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, dzdx_interior_nodal, "dzdx_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, eps2_interior_nodal, "eps2_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, g_ij_interior_nodal, "g_ij_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, gij_interior_nodal, "gij_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, gxxj_interior_nodal, "gxxj_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, gxyj_interior_nodal, "gxyj_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, gxzj_interior_nodal, "gxzj_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, gyyj_interior_nodal, "gyyj_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, jacobgeo_interior_nodal, "jacobgeo_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, jacobgeo_inv_interior_nodal, "jacobgeo_inv_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, jacobtot_interior_nodal, "jacobtot_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, jacobtot_inv_interior_nodal, "jacobtot_inv_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, normals_interior_nodal, "normals_interior_nodal");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, rtg33inv_interior_nodal, "rtg33inv_interior_nodal");
    
    // Write surface 1 nodal arrays 
    WRITE_NODAL(&ngrid_surface, &nrange_surface, dxdz_surface1_nodal, "dxdz_surface1_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, bmag_surface1_nodal, "bmag_surface1_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, curlbhat_surface1_nodal, "curlbhat_surface1_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, B3_surface1_nodal, "B3_surface1_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, b_i_surface1_nodal, "b_i_surface1_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, cmag_surface1_nodal, "cmag_surface1_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, normcurlbhat_surface1_nodal, "normcurlbhat_surface1_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, jacobgeo_surface1_nodal, "jacobgeo_surface1_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, jacobtot_inv_surface1_nodal, "jacobtot_inv_surface1_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, normals_surface1_nodal, "normals_surface1_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, lenr_surface1_nodal, "lenr_surface1_nodal");
    
    // Write surface 2 nodal arrays 
    WRITE_NODAL(&ngrid_surface, &nrange_surface, dxdz_surface2_nodal, "dxdz_surface2_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, bmag_surface2_nodal, "bmag_surface2_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, curlbhat_surface2_nodal, "curlbhat_surface2_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, B3_surface2_nodal, "B3_surface2_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, b_i_surface2_nodal, "b_i_surface2_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, cmag_surface2_nodal, "cmag_surface2_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, normcurlbhat_surface2_nodal, "normcurlbhat_surface2_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, jacobgeo_surface2_nodal, "jacobgeo_surface2_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, jacobtot_inv_surface2_nodal, "jacobtot_inv_surface2_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, normals_surface2_nodal, "normals_surface2_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, lenr_surface2_nodal, "lenr_surface2_nodal");
    
    // Write surface 3 nodal arrays 
    WRITE_NODAL(&ngrid_surface, &nrange_surface, dxdz_surface3_nodal, "dxdz_surface3_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, bmag_surface3_nodal, "bmag_surface3_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, curlbhat_surface3_nodal, "curlbhat_surface3_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, B3_surface3_nodal, "B3_surface3_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, b_i_surface3_nodal, "b_i_surface3_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, cmag_surface3_nodal, "cmag_surface3_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, normcurlbhat_surface3_nodal, "normcurlbhat_surface3_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, jacobgeo_surface3_nodal, "jacobgeo_surface3_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, jacobtot_inv_surface3_nodal, "jacobtot_inv_surface3_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, normals_surface3_nodal, "normals_surface3_nodal");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, lenr_surface3_nodal, "lenr_surface3_nodal");
    
    // Write corner modal arrays 
    WRITE_MODAL(&grid_orig, &comp_range, mc2nu_pos_corner_modal, "mc2nu_pos_corner_modal");
    WRITE_MODAL(&grid_orig, &comp_range, mc2p_corner_modal, "mc2p_corner_modal");
    WRITE_MODAL(&grid_orig, &comp_range, bmag_corner_modal, "bmag_corner_modal");
    
    // Write interior modal arrays 
    WRITE_MODAL(&grid_orig, &comp_range, dxdz_interior_modal, "dxdz_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, bmag_interior_modal, "bmag_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, curlbhat_interior_modal, "curlbhat_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, bcart_interior_modal, "bcart_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, B3_interior_modal, "B3_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, b_i_interior_modal, "b_i_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, bmag_inv_interior_modal, "bmag_inv_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, bmag_inv_sq_interior_modal, "bmag_inv_sq_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, bioverJB_interior_modal, "bioverJB_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, cmag_interior_modal, "cmag_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualcurlbhat_interior_modal, "dualcurlbhat_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualcurlbhatoverB_interior_modal, "dualcurlbhatoverB_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dzdx_interior_modal, "dzdx_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, eps2_interior_modal, "eps2_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, g_ij_interior_modal, "g_ij_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, gij_interior_modal, "gij_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, gxxj_interior_modal, "gxxj_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, gxyj_interior_modal, "gxyj_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, gxzj_interior_modal, "gxzj_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, gyyj_interior_modal, "gyyj_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, jacobgeo_interior_modal, "jacobgeo_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, jacobgeo_inv_interior_modal, "jacobgeo_inv_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, jacobtot_interior_modal, "jacobtot_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, jacobtot_inv_interior_modal, "jacobtot_inv_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, normals_interior_modal, "normals_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, rtg33inv_interior_modal, "rtg33inv_interior_modal");
    
    // Write surface 1 modal arrays 
    WRITE_MODAL(&grid_orig, &comp_range, dxdz_s1_modal, "dxdz_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, bmag_s1_modal, "bmag_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, curlbhat_s1_modal, "curlbhat_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, B3_s1_modal, "B3_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, b_i_s1_modal, "b_i_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, cmag_s1_modal, "cmag_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, normcurlbhat_s1_modal, "normcurlbhat_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, jacobgeo_s1_modal, "jacobgeo_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, jacobtot_inv_s1_modal, "jacobtot_inv_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, normals_s1_modal, "normals_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, lenr_s1_modal, "lenr_s1_modal");
    
    // Write surface 2 modal arrays 
    WRITE_MODAL(&grid_orig, &comp_range, dxdz_s2_modal, "dxdz_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, bmag_s2_modal, "bmag_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, curlbhat_s2_modal, "curlbhat_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, B3_s2_modal, "B3_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, b_i_s2_modal, "b_i_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, cmag_s2_modal, "cmag_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, normcurlbhat_s2_modal, "normcurlbhat_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, jacobgeo_s2_modal, "jacobgeo_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, jacobtot_inv_s2_modal, "jacobtot_inv_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, normals_s2_modal, "normals_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, lenr_s2_modal, "lenr_s2_modal");
    
    // Write surface 3 modal arrays 
    WRITE_MODAL(&grid_orig, &comp_range, dxdz_s3_modal, "dxdz_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, bmag_s3_modal, "bmag_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, curlbhat_s3_modal, "curlbhat_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, B3_s3_modal, "B3_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, b_i_s3_modal, "b_i_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, cmag_s3_modal, "cmag_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, normcurlbhat_s3_modal, "normcurlbhat_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, jacobgeo_s3_modal, "jacobgeo_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, jacobtot_inv_s3_modal, "jacobtot_inv_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, normals_s3_modal, "normals_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, lenr_s3_modal, "lenr_s3_modal");

    // Clean up everything
    // Release corner loaded arrays
    gkyl_array_release(mc2nu_pos_corner_nodal);
    gkyl_array_release(mc2p_corner_nodal);
    gkyl_array_release(bmag_corner_nodal);
    
    // Release corner modal arrays
    gkyl_array_release(mc2nu_pos_corner_modal);
    gkyl_array_release(mc2p_corner_modal);
    gkyl_array_release(bmag_corner_modal);
    
    // Release interior loaded arrays
    gkyl_array_release(dxdz_interior_nodal);
    gkyl_array_release(bmag_interior_nodal);
    gkyl_array_release(curlbhat_interior_nodal);
    gkyl_array_release(bcart_interior_nodal);
    
    // Release interior calculated arrays
    gkyl_array_release(B3_interior_nodal);
    gkyl_array_release(b_i_interior_nodal);
    gkyl_array_release(bcart_interior_nodal);
    gkyl_array_release(bmag_inv_interior_nodal);
    gkyl_array_release(bmag_inv_sq_interior_nodal);
    gkyl_array_release(bioverJB_interior_nodal);
    gkyl_array_release(cmag_interior_nodal);
    gkyl_array_release(dualcurlbhat_interior_nodal);
    gkyl_array_release(dualcurlbhatoverB_interior_nodal);
    gkyl_array_release(dzdx_interior_nodal);
    gkyl_array_release(eps2_interior_nodal);
    gkyl_array_release(g_ij_interior_nodal);
    gkyl_array_release(gij_interior_nodal);
    gkyl_array_release(gxxj_interior_nodal);
    gkyl_array_release(gxyj_interior_nodal);
    gkyl_array_release(gxzj_interior_nodal);
    gkyl_array_release(gyyj_interior_nodal);
    gkyl_array_release(jacobgeo_interior_nodal);
    gkyl_array_release(jacobgeo_inv_interior_nodal);
    gkyl_array_release(jacobtot_interior_nodal);
    gkyl_array_release(jacobtot_inv_interior_nodal);
    gkyl_array_release(normals_interior_nodal);
    gkyl_array_release(rtg33inv_interior_nodal);
    
    // Release interior modal arrays
    gkyl_array_release(dxdz_interior_modal);
    gkyl_array_release(bmag_interior_modal);
    gkyl_array_release(curlbhat_interior_modal);
    gkyl_array_release(bcart_interior_modal);
    gkyl_array_release(B3_interior_modal);
    gkyl_array_release(b_i_interior_modal);
    gkyl_array_release(bmag_inv_interior_modal);
    gkyl_array_release(bmag_inv_sq_interior_modal);
    gkyl_array_release(bioverJB_interior_modal);
    gkyl_array_release(cmag_interior_modal);
    gkyl_array_release(dualcurlbhat_interior_modal);
    gkyl_array_release(dualcurlbhatoverB_interior_modal);
    gkyl_array_release(dzdx_interior_modal);
    gkyl_array_release(eps2_interior_modal);
    gkyl_array_release(g_ij_interior_modal);
    gkyl_array_release(gij_interior_modal);
    gkyl_array_release(gxxj_interior_modal);
    gkyl_array_release(gxyj_interior_modal);
    gkyl_array_release(gxzj_interior_modal);
    gkyl_array_release(gyyj_interior_modal);
    gkyl_array_release(jacobgeo_interior_modal);
    gkyl_array_release(jacobgeo_inv_interior_modal);
    gkyl_array_release(jacobtot_interior_modal);
    gkyl_array_release(jacobtot_inv_interior_modal);
    gkyl_array_release(normals_interior_modal);
    gkyl_array_release(rtg33inv_interior_modal);
    
    // Release surface 1 loaded arrays
    gkyl_array_release(dxdz_surface1_nodal);
    gkyl_array_release(bmag_surface1_nodal);
    gkyl_array_release(curlbhat_surface1_nodal);
    
    // Release surface 1 calculated arrays
    gkyl_array_release(B3_surface1_nodal);
    gkyl_array_release(b_i_surface1_nodal);
    gkyl_array_release(cmag_surface1_nodal);
    gkyl_array_release(normcurlbhat_surface1_nodal);
    gkyl_array_release(jacobgeo_surface1_nodal);
    gkyl_array_release(jacobtot_inv_surface1_nodal);
    gkyl_array_release(normals_surface1_nodal);
    gkyl_array_release(lenr_surface1_nodal);
    
    // Release surface 1 modal arrays
    gkyl_array_release(dxdz_s1_modal);
    gkyl_array_release(bmag_s1_modal);
    gkyl_array_release(curlbhat_s1_modal);
    gkyl_array_release(B3_s1_modal);
    gkyl_array_release(b_i_s1_modal);
    gkyl_array_release(cmag_s1_modal);
    gkyl_array_release(normcurlbhat_s1_modal);
    gkyl_array_release(jacobgeo_s1_modal);
    gkyl_array_release(jacobtot_inv_s1_modal);
    gkyl_array_release(normals_s1_modal);
    gkyl_array_release(lenr_s1_modal);
    
    // Release surface 2 loaded arrays
    gkyl_array_release(dxdz_surface2_nodal);
    gkyl_array_release(bmag_surface2_nodal);
    gkyl_array_release(curlbhat_surface2_nodal);
    
    // Release surface 2 calculated arrays
    gkyl_array_release(B3_surface2_nodal);
    gkyl_array_release(b_i_surface2_nodal);
    gkyl_array_release(cmag_surface2_nodal);
    gkyl_array_release(normcurlbhat_surface2_nodal);
    gkyl_array_release(jacobgeo_surface2_nodal);
    gkyl_array_release(jacobtot_inv_surface2_nodal);
    gkyl_array_release(normals_surface2_nodal);
    gkyl_array_release(lenr_surface2_nodal);
    
    // Release surface 2 modal arrays
    gkyl_array_release(dxdz_s2_modal);
    gkyl_array_release(bmag_s2_modal);
    gkyl_array_release(curlbhat_s2_modal);
    gkyl_array_release(B3_s2_modal);
    gkyl_array_release(b_i_s2_modal);
    gkyl_array_release(cmag_s2_modal);
    gkyl_array_release(normcurlbhat_s2_modal);
    gkyl_array_release(jacobgeo_s2_modal);
    gkyl_array_release(jacobtot_inv_s2_modal);
    gkyl_array_release(normals_s2_modal);
    gkyl_array_release(lenr_s2_modal);
    
    // Release surface 3 loaded arrays
    gkyl_array_release(dxdz_surface3_nodal);
    gkyl_array_release(bmag_surface3_nodal);
    gkyl_array_release(curlbhat_surface3_nodal);
    
    // Release surface 3 calculated arrays
    gkyl_array_release(B3_surface3_nodal);
    gkyl_array_release(b_i_surface3_nodal);
    gkyl_array_release(cmag_surface3_nodal);
    gkyl_array_release(normcurlbhat_surface3_nodal);
    gkyl_array_release(jacobgeo_surface3_nodal);
    gkyl_array_release(jacobtot_inv_surface3_nodal);
    gkyl_array_release(normals_surface3_nodal);
    gkyl_array_release(lenr_surface3_nodal);
    
    // Release surface 3 modal arrays
    gkyl_array_release(dxdz_s3_modal);
    gkyl_array_release(bmag_s3_modal);
    gkyl_array_release(curlbhat_s3_modal);
    gkyl_array_release(B3_s3_modal);
    gkyl_array_release(b_i_s3_modal);
    gkyl_array_release(cmag_s3_modal);
    gkyl_array_release(normcurlbhat_s3_modal);
    gkyl_array_release(jacobgeo_s3_modal);
    gkyl_array_release(jacobtot_inv_s3_modal);
    gkyl_array_release(normals_s3_modal);
    gkyl_array_release(lenr_s3_modal);
    
    // Release other resources
    gkyl_array_header_info_release(&hdr);
    gkyl_nodal_ops_release(n2m);
}

TEST_LIST = {
    { "test_3x_p1", test_3xp1},
    { NULL, NULL },
  };



