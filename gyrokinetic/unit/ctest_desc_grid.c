#include <acutest.h>
#include <stdio.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>

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
    struct gkyl_array* tangents_corner = load_gkyl_file_to_nodal("drpz_draz_corner.gkyl");
    struct gkyl_array* Bmag_corner = load_gkyl_file_to_nodal("Bmag_corner.gkyl");
    struct gkyl_array* curlBhat_corner = load_gkyl_file_to_nodal("curl_B_hat_xyz_corner.gkyl");
    struct gkyl_array* raz_corner = load_gkyl_file_to_nodal("raz_corner.gkyl");
    struct gkyl_array* rpz_corner = load_gkyl_file_to_nodal("rpz_corner.gkyl");
    
    // Check for NULL pointers (file loading failures)
    TEST_CHECK(tangents_corner != NULL);
    TEST_CHECK(Bmag_corner != NULL);
    TEST_CHECK(curlBhat_corner != NULL);
    TEST_CHECK(raz_corner != NULL);
    TEST_CHECK(rpz_corner != NULL);
    if (!tangents_corner || !Bmag_corner || !curlBhat_corner || !raz_corner || !rpz_corner) {
        return; // Exit early if files failed to load
    }
    // Calculate derived quantities for corner nodes
    // Create a simple 1D range for 100000 nodes
    int shape_corner[1] = {100000};
    struct gkyl_range nrange_corner;
    gkyl_range_init_from_shape1(&nrange_corner, 1, shape_corner);
    struct gkyl_range_iter iter_corner;
    gkyl_range_iter_init(&iter_corner, &nrange_corner);  
    // initialize jacobian and other arrays which are to be calculated (mkarr(false, number of components, nrange_corner.volume))
    struct gkyl_array* duals_corner = mkarr(false, 9, nrange_corner.volume);
    struct gkyl_array* jacobgeo_corner = mkarr(false, 1, nrange_corner.volume);
    struct gkyl_array* g_ij_corner = mkarr(false, 9, nrange_corner.volume);
    struct gkyl_array* B3_corner = mkarr(false, 1, nrange_corner.volume);
    struct gkyl_array* dualcurlBhat_corner = mkarr(false, 3, nrange_corner.volume);
    struct gkyl_array* dualcurlBhatoverB_corner = mkarr(false, 3, nrange_corner.volume);
    struct gkyl_array* rtg33inv_corner = mkarr(false, 1, nrange_corner.volume);
    struct gkyl_array* bioverJB_corner = mkarr(false, 3, nrange_corner.volume);
    struct gkyl_array* dualmags_corner = mkarr(false, 3, nrange_corner.volume);
    struct gkyl_array* dualnormals_corner = mkarr(false, 9, nrange_corner.volume);
    
    while (gkyl_range_iter_next(&iter_corner)) {
        // Get current node index
        long idx = gkyl_range_idx(&nrange_corner, iter_corner.idx);
        
        double *tang_data = (double*)tangents_corner->data;
        double *bmag_data = (double*)Bmag_corner->data;
        double *curl_data = (double*)curlBhat_corner->data;
        double *dual_data = (double*)duals_corner->data;
        double *gij_data = (double*)g_ij_corner->data;
        double *jacob_data = (double*)jacobgeo_corner->data;
        double *b3_data = (double*)B3_corner->data;
        double *dcb_data = (double*)dualcurlBhat_corner->data;
        double *dcbob_data = (double*)dualcurlBhatoverB_corner->data;
        double *rtg_data = (double*)rtg33inv_corner->data;
        double *biojb_data = (double*)bioverJB_corner->data;
        double *dmag_data = (double*)dualmags_corner->data;
        double *dnorm_data = (double*)dualnormals_corner->data;
        
        double e_1[3] = {tang_data[idx*9+0], tang_data[idx*9+1], tang_data[idx*9+2]};
        double e_2[3] = {tang_data[idx*9+3], tang_data[idx*9+4], tang_data[idx*9+5]};
        double e_3[3] = {tang_data[idx*9+6], tang_data[idx*9+7], tang_data[idx*9+8]};
        double curlB[3] = {curl_data[idx*3+0], curl_data[idx*3+1], curl_data[idx*3+2]};
        double Bmag = bmag_data[idx];
        
        // Jacobian: J = e_1 · (e_2 × e_3)
        double e2xe3[3] = {e_2[1]*e_3[2]-e_2[2]*e_3[1], e_2[2]*e_3[0]-e_2[0]*e_3[2], e_2[0]*e_3[1]-e_2[1]*e_3[0]};
        double J = e_1[0]*e2xe3[0] + e_1[1]*e2xe3[1] + e_1[2]*e2xe3[2];
        jacob_data[idx] = J;
        if (fabs(J) < 1e-12) {
            // Jacobian is too small, skip this node or set to a small value
            J = 1e-12;
            jacob_data[idx] = J;
        }
        double Jinv = 1.0/J;
        
        // Dual vectors: e^i = (1/J) * (e_j × e_k)
        double e3xe1[3] = {e_3[1]*e_1[2]-e_3[2]*e_1[1], e_3[2]*e_1[0]-e_3[0]*e_1[2], e_3[0]*e_1[1]-e_3[1]*e_1[0]};
        double e1xe2[3] = {e_1[1]*e_2[2]-e_1[2]*e_2[1], e_1[2]*e_2[0]-e_1[0]*e_2[2], e_1[0]*e_2[1]-e_1[1]*e_2[0]};
        for (int i=0; i<3; i++) {
            dual_data[idx*9+0*3+i] = Jinv * e2xe3[i];  // e^1
            dual_data[idx*9+1*3+i] = Jinv * e3xe1[i];  // e^2
            dual_data[idx*9+2*3+i] = Jinv * e1xe2[i];  // e^3
        }
        
        // Metric tensor: g_ij = e_i · e_j
        for (int i=0; i<3; i++) for (int j=0; j<3; j++) {
            double *ei = (i==0)?e_1:(i==1)?e_2:e_3;
            double *ej = (j==0)?e_1:(j==1)?e_2:e_3;
            gij_data[idx*9+i*3+j] = ei[0]*ej[0] + ei[1]*ej[1] + ei[2]*ej[2];
        }
        
        double g33 = gij_data[idx*9+2*3+2];
        double sqrt_g33 = sqrt(g33);
        double rtg33inv = 1.0/sqrt_g33;
        rtg_data[idx] = rtg33inv;
        
        // B3 = Bmag / sqrt(g_33)
        b3_data[idx] = Bmag * rtg33inv;
        
        // dualcurlbhat = e^i · curlBhat and dualcurlBhatoverB
        for (int i=0; i<3; i++) {
            double *ei = &dual_data[idx*9+i*3];
            dcb_data[idx*3+i] = ei[0]*curlB[0] + ei[1]*curlB[1] + ei[2]*curlB[2];
            dcbob_data[idx*3+i] = dcb_data[idx*3+i] / Bmag;
        }
        
        // bioverJB = g_i3 / sqrt(g_33) / J / Bmag 
        for (int i=0; i<3; i++) {
            biojb_data[idx*3+i] = gij_data[idx*9+i*3+2] / sqrt_g33 / J / Bmag;
        }
        
        // Dual vector magnitudes: |e^i|
        for (int i=0; i<3; i++) {
            double *ei = &dual_data[idx*9+i*3];
            dmag_data[idx*3+i] = sqrt(ei[0]*ei[0] + ei[1]*ei[1] + ei[2]*ei[2]);
        }
        
        // Normal vectors: e^i / |e^i|
        for (int i=0; i<3; i++) {
            double mag = dmag_data[idx*3+i];
            double *ei = &dual_data[idx*9+i*3];
            if (mag > 1e-12) {  // Avoid division by zero
                for (int j=0; j<3; j++) dnorm_data[idx*9+i*3+j] = ei[j] / mag;
            } else {
                // If magnitude is zero, set normal to zero vector
                for (int j=0; j<3; j++) dnorm_data[idx*9+i*3+j] = 0.0;
            }
        }
    }

    // Do for interior
    struct gkyl_array* tangents_interior = load_gkyl_file_to_nodal("drpz_draz_I.gkyl");
    struct gkyl_array* Bmag_interior = load_gkyl_file_to_nodal("Bmag_I.gkyl");
    struct gkyl_array* curlBhat_interior = load_gkyl_file_to_nodal("curl_B_hat_xyz_I.gkyl");
    struct gkyl_array* raz_interior = load_gkyl_file_to_nodal("raz_I.gkyl");
    struct gkyl_array* rpz_interior = load_gkyl_file_to_nodal("rpz_I.gkyl");
    
    TEST_CHECK(tangents_interior != NULL);
    TEST_CHECK(Bmag_interior != NULL);
    TEST_CHECK(curlBhat_interior != NULL);
    TEST_CHECK(raz_interior != NULL);
    TEST_CHECK(rpz_interior != NULL);
    if (!tangents_interior || !Bmag_interior || !curlBhat_interior || !raz_interior || !rpz_interior) {
        return;
    }
    // Then calculate derived quantities
    // Create a simple 1D range for 705672 nodes
    int shape_interior[1] = {705672};
    struct gkyl_range nrange_interior;
    gkyl_range_init_from_shape1(&nrange_interior, 1, shape_interior);
    struct gkyl_range_iter iter_interior;
    gkyl_range_iter_init(&iter_interior, &nrange_interior);  
    // initialize jacobian and other arrays which are to be calculated (mkarr(false, number of components, nrange_interior.volume))
    struct gkyl_array* duals_interior = mkarr(false, 9, nrange_interior.volume);
    struct gkyl_array* jacobgeo_interior = mkarr(false, 1, nrange_interior.volume);
    struct gkyl_array* g_ij_interior = mkarr(false, 9, nrange_interior.volume);
    struct gkyl_array* B3_interior = mkarr(false, 1, nrange_interior.volume);
    struct gkyl_array* dualcurlBhat_interior = mkarr(false, 3, nrange_interior.volume);
    struct gkyl_array* dualcurlBhatoverB_interior = mkarr(false, 3, nrange_interior.volume);
    struct gkyl_array* rtg33inv_interior = mkarr(false, 1, nrange_interior.volume);
    struct gkyl_array* bioverJB_interior = mkarr(false, 3, nrange_interior.volume);
    struct gkyl_array* dualmags_interior = mkarr(false, 3, nrange_interior.volume);
    struct gkyl_array* dualnormals_interior = mkarr(false, 9, nrange_interior.volume);
    
    while (gkyl_range_iter_next(&iter_interior)) {
        // Get current node index
        long idx = gkyl_range_idx(&nrange_interior, iter_interior.idx);
        
        double *tang_data = (double*)tangents_interior->data;
        double *bmag_data = (double*)Bmag_interior->data;
        double *curl_data = (double*)curlBhat_interior->data;
        double *dual_data = (double*)duals_interior->data;
        double *gij_data = (double*)g_ij_interior->data;
        double *jacob_data = (double*)jacobgeo_interior->data;
        double *b3_data = (double*)B3_interior->data;
        double *dcb_data = (double*)dualcurlBhat_interior->data;
        double *dcbob_data = (double*)dualcurlBhatoverB_interior->data;
        double *rtg_data = (double*)rtg33inv_interior->data;
        double *biojb_data = (double*)bioverJB_interior->data;
        double *dmag_data = (double*)dualmags_interior->data;
        double *dnorm_data = (double*)dualnormals_interior->data;
        
        double e_1[3] = {tang_data[idx*9+0], tang_data[idx*9+1], tang_data[idx*9+2]};
        double e_2[3] = {tang_data[idx*9+3], tang_data[idx*9+4], tang_data[idx*9+5]};
        double e_3[3] = {tang_data[idx*9+6], tang_data[idx*9+7], tang_data[idx*9+8]};
        double curlB[3] = {curl_data[idx*3+0], curl_data[idx*3+1], curl_data[idx*3+2]};
        double Bmag = bmag_data[idx];
        
        // Jacobian: J = e_1 · (e_2 × e_3)
        double e2xe3[3] = {e_2[1]*e_3[2]-e_2[2]*e_3[1], e_2[2]*e_3[0]-e_2[0]*e_3[2], e_2[0]*e_3[1]-e_2[1]*e_3[0]};
        double J = e_1[0]*e2xe3[0] + e_1[1]*e2xe3[1] + e_1[2]*e2xe3[2];
        jacob_data[idx] = J;
        if (fabs(J) < 1e-12) {
            // Jacobian is too small, skip this node or set to a small value
            J = 1e-12;
            jacob_data[idx] = J;
        }
        double Jinv = 1.0/J;
        
        // Dual vectors: e^i = (1/J) * (e_j × e_k)
        double e3xe1[3] = {e_3[1]*e_1[2]-e_3[2]*e_1[1], e_3[2]*e_1[0]-e_3[0]*e_1[2], e_3[0]*e_1[1]-e_3[1]*e_1[0]};
        double e1xe2[3] = {e_1[1]*e_2[2]-e_1[2]*e_2[1], e_1[2]*e_2[0]-e_1[0]*e_2[2], e_1[0]*e_2[1]-e_1[1]*e_2[0]};
        for (int i=0; i<3; i++) {
            dual_data[idx*9+0*3+i] = Jinv * e2xe3[i];  // e^1
            dual_data[idx*9+1*3+i] = Jinv * e3xe1[i];  // e^2
            dual_data[idx*9+2*3+i] = Jinv * e1xe2[i];  // e^3
        }
        
        // Metric tensor: g_ij = e_i · e_j
        for (int i=0; i<3; i++) for (int j=0; j<3; j++) {
            double *ei = (i==0)?e_1:(i==1)?e_2:e_3;
            double *ej = (j==0)?e_1:(j==1)?e_2:e_3;
            gij_data[idx*9+i*3+j] = ei[0]*ej[0] + ei[1]*ej[1] + ei[2]*ej[2];
        }
        
        double g33 = gij_data[idx*9+2*3+2];
        double sqrt_g33 = sqrt(g33);
        double rtg33inv = 1.0/sqrt_g33;
        rtg_data[idx] = rtg33inv;
        
        // B3 = Bmag / sqrt(g_33)
        b3_data[idx] = Bmag * rtg33inv;
        
        // dualcurlbhat = e^i · curlBhat and dualcurlBhatoverB
        for (int i=0; i<3; i++) {
            double *ei = &dual_data[idx*9+i*3];
            dcb_data[idx*3+i] = ei[0]*curlB[0] + ei[1]*curlB[1] + ei[2]*curlB[2];
            dcbob_data[idx*3+i] = dcb_data[idx*3+i] / Bmag;
        }
        
        // bioverJB = g_i3 / sqrt(g_33) / J / Bmag 
        for (int i=0; i<3; i++) {
            biojb_data[idx*3+i] = gij_data[idx*9+i*3+2] / sqrt_g33 / J / Bmag;
        }
        
        // Dual vector magnitudes: |e^i|
        for (int i=0; i<3; i++) {
            double *ei = &dual_data[idx*9+i*3];
            dmag_data[idx*3+i] = sqrt(ei[0]*ei[0] + ei[1]*ei[1] + ei[2]*ei[2]);
        }
        
        // Normal vectors: e^i / |e^i|
        for (int i=0; i<3; i++) {
            double mag = dmag_data[idx*3+i];
            double *ei = &dual_data[idx*9+i*3];
            if (mag > 1e-12) {  // Avoid division by zero
                for (int j=0; j<3; j++) dnorm_data[idx*9+i*3+j] = ei[j] / mag;
            } else {
                // If magnitude is zero, set normal to zero vector
                for (int j=0; j<3; j++) dnorm_data[idx*9+i*3+j] = 0.0;
            }
        }
    }


    // Do again for surfaces
    struct gkyl_array* tangents_surface1 = load_gkyl_file_to_nodal("drpz_draz_S1.gkyl");
    struct gkyl_array* Bmag_surface1 = load_gkyl_file_to_nodal("Bmag_S1.gkyl");
    struct gkyl_array* curlBhat_surface1 = load_gkyl_file_to_nodal("curl_B_hat_xyz_S1.gkyl");
    struct gkyl_array* raz_surface1 = load_gkyl_file_to_nodal("raz_S1.gkyl");
    struct gkyl_array* rpz_surface1 = load_gkyl_file_to_nodal("rpz_S1.gkyl");
    struct gkyl_array* tangents_surface2 = load_gkyl_file_to_nodal("drpz_draz_S2.gkyl");
    struct gkyl_array* Bmag_surface2 = load_gkyl_file_to_nodal("Bmag_S2.gkyl");
    struct gkyl_array* curlBhat_surface2 = load_gkyl_file_to_nodal("curl_B_hat_xyz_S2.gkyl");
    struct gkyl_array* raz_surface2 = load_gkyl_file_to_nodal("raz_S2.gkyl");
    struct gkyl_array* rpz_surface2 = load_gkyl_file_to_nodal("rpz_S2.gkyl");
    struct gkyl_array* tangents_surface3 = load_gkyl_file_to_nodal("drpz_draz_S3.gkyl");
    struct gkyl_array* Bmag_surface3 = load_gkyl_file_to_nodal("Bmag_S3.gkyl");
    struct gkyl_array* curlBhat_surface3 = load_gkyl_file_to_nodal("curl_B_hat_xyz_S3.gkyl");
    struct gkyl_array* raz_surface3 = load_gkyl_file_to_nodal("raz_S3.gkyl");
    struct gkyl_array* rpz_surface3 = load_gkyl_file_to_nodal("rpz_S3.gkyl");
    
    // Check all surface arrays loaded successfully
    TEST_CHECK(tangents_surface1 != NULL && Bmag_surface1 != NULL && curlBhat_surface1 != NULL && 
               raz_surface1 != NULL && rpz_surface1 != NULL);
    TEST_CHECK(tangents_surface2 != NULL && Bmag_surface2 != NULL && curlBhat_surface2 != NULL && 
               raz_surface2 != NULL && rpz_surface2 != NULL);
    TEST_CHECK(tangents_surface3 != NULL && Bmag_surface3 != NULL && curlBhat_surface3 != NULL && 
               raz_surface3 != NULL && rpz_surface3 != NULL);
    if (!tangents_surface1 || !Bmag_surface1 || !curlBhat_surface1 || !raz_surface1 || !rpz_surface1 ||
        !tangents_surface2 || !Bmag_surface2 || !curlBhat_surface2 || !raz_surface2 || !rpz_surface2 ||
        !tangents_surface3 || !Bmag_surface3 || !curlBhat_surface3 || !raz_surface3 || !rpz_surface3) {
        return;
    }
    // Derived quantities
    // Create a simple 1D range for 356400 nodes
    int shape_surface[1] = {356400};
    struct gkyl_range nrange_surface;
    gkyl_range_init_from_shape1(&nrange_surface, 1, shape_surface);
    struct gkyl_range_iter iter_surface;
    gkyl_range_iter_init(&iter_surface, &nrange_surface);  

    // calculate quantities for surface 1
    struct gkyl_array* duals_surface1 = mkarr(false, 9, nrange_surface.volume);
    struct gkyl_array* jacobgeo_surface1 = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* g_ij_surface1 = mkarr(false, 9, nrange_surface.volume);
    struct gkyl_array* B3_surface1 = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* dualcurlBhat_surface1 = mkarr(false, 3, nrange_surface.volume);
    struct gkyl_array* dualcurlBhatoverB_surface1 = mkarr(false, 3, nrange_surface.volume);
    struct gkyl_array* rtg33inv_surface1 = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* bioverJB_surface1 = mkarr(false, 3, nrange_surface.volume);
    struct gkyl_array* dualmags_surface1 = mkarr(false, 3, nrange_surface.volume);
    struct gkyl_array* dualnormals_surface1 = mkarr(false, 9, nrange_surface.volume);
    
    while (gkyl_range_iter_next(&iter_surface)) {
        // Get current node index
        long idx = gkyl_range_idx(&nrange_surface, iter_surface.idx);
        
        double *tang_data = (double*)tangents_surface1->data;
        double *bmag_data = (double*)Bmag_surface1->data;
        double *curl_data = (double*)curlBhat_surface1->data;
        double *dual_data = (double*)duals_surface1->data;
        double *gij_data = (double*)g_ij_surface1->data;
        double *jacob_data = (double*)jacobgeo_surface1->data;
        double *b3_data = (double*)B3_surface1->data;
        double *dcb_data = (double*)dualcurlBhat_surface1->data;
        double *dcbob_data = (double*)dualcurlBhatoverB_surface1->data;
        double *rtg_data = (double*)rtg33inv_surface1->data;
        double *biojb_data = (double*)bioverJB_surface1->data;
        double *dmag_data = (double*)dualmags_surface1->data;
        double *dnorm_data = (double*)dualnormals_surface1->data;
        
        double e_1[3] = {tang_data[idx*9+0], tang_data[idx*9+1], tang_data[idx*9+2]};
        double e_2[3] = {tang_data[idx*9+3], tang_data[idx*9+4], tang_data[idx*9+5]};
        double e_3[3] = {tang_data[idx*9+6], tang_data[idx*9+7], tang_data[idx*9+8]};
        double curlB[3] = {curl_data[idx*3+0], curl_data[idx*3+1], curl_data[idx*3+2]};
        double Bmag = bmag_data[idx];
        
        // Jacobian: J = e_1 · (e_2 × e_3)
        double e2xe3[3] = {e_2[1]*e_3[2]-e_2[2]*e_3[1], e_2[2]*e_3[0]-e_2[0]*e_3[2], e_2[0]*e_3[1]-e_2[1]*e_3[0]};
        double J = e_1[0]*e2xe3[0] + e_1[1]*e2xe3[1] + e_1[2]*e2xe3[2];
        jacob_data[idx] = J;
        if (fabs(J) < 1e-12) {
            // Jacobian is too small, skip this node or set to a small value
            J = 1e-12;
            jacob_data[idx] = J;
        }
        double Jinv = 1.0/J;
        
        // Dual vectors: e^i = (1/J) * (e_j × e_k)
        double e3xe1[3] = {e_3[1]*e_1[2]-e_3[2]*e_1[1], e_3[2]*e_1[0]-e_3[0]*e_1[2], e_3[0]*e_1[1]-e_3[1]*e_1[0]};
        double e1xe2[3] = {e_1[1]*e_2[2]-e_1[2]*e_2[1], e_1[2]*e_2[0]-e_1[0]*e_2[2], e_1[0]*e_2[1]-e_1[1]*e_2[0]};
        for (int i=0; i<3; i++) {
            dual_data[idx*9+0*3+i] = Jinv * e2xe3[i];  // e^1
            dual_data[idx*9+1*3+i] = Jinv * e3xe1[i];  // e^2
            dual_data[idx*9+2*3+i] = Jinv * e1xe2[i];  // e^3
        }
        
        // Metric tensor: g_ij = e_i · e_j
        for (int i=0; i<3; i++) for (int j=0; j<3; j++) {
            double *ei = (i==0)?e_1:(i==1)?e_2:e_3;
            double *ej = (j==0)?e_1:(j==1)?e_2:e_3;
            gij_data[idx*9+i*3+j] = ei[0]*ej[0] + ei[1]*ej[1] + ei[2]*ej[2];
        }
        
        double g33 = gij_data[idx*9+2*3+2];
        double sqrt_g33 = sqrt(g33);
        double rtg33inv = 1.0/sqrt_g33;
        rtg_data[idx] = rtg33inv;
        
        // B3 = Bmag / sqrt(g_33)
        b3_data[idx] = Bmag * rtg33inv;
        
        // dualcurlbhat = e^i · curlBhat and dualcurlBhatoverB
        for (int i=0; i<3; i++) {
            double *ei = &dual_data[idx*9+i*3];
            dcb_data[idx*3+i] = ei[0]*curlB[0] + ei[1]*curlB[1] + ei[2]*curlB[2];
            dcbob_data[idx*3+i] = dcb_data[idx*3+i] / Bmag;
        }
        
        // bioverJB = g_i3 / sqrt(g_33) / J / Bmag 
        for (int i=0; i<3; i++) {
            biojb_data[idx*3+i] = gij_data[idx*9+i*3+2] / sqrt_g33 / J / Bmag;
        }
        
        // Dual vector magnitudes: |e^i|
        for (int i=0; i<3; i++) {
            double *ei = &dual_data[idx*9+i*3];
            dmag_data[idx*3+i] = sqrt(ei[0]*ei[0] + ei[1]*ei[1] + ei[2]*ei[2]);
        }
        
        // Normal vectors: e^i / |e^i|
        for (int i=0; i<3; i++) {
            double mag = dmag_data[idx*3+i];
            double *ei = &dual_data[idx*9+i*3];
            if (mag > 1e-12) {  // Avoid division by zero
                for (int j=0; j<3; j++) dnorm_data[idx*9+i*3+j] = ei[j] / mag;
            } else {
                // If magnitude is zero, set normal to zero vector
                for (int j=0; j<3; j++) dnorm_data[idx*9+i*3+j] = 0.0;
            }
        }
    }

     // calculate quantities for surface 2
    struct gkyl_array* duals_surface2 = mkarr(false, 9, nrange_surface.volume);
    struct gkyl_array* jacobgeo_surface2 = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* g_ij_surface2 = mkarr(false, 9, nrange_surface.volume);
    struct gkyl_array* B3_surface2 = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* dualcurlBhat_surface2 = mkarr(false, 3, nrange_surface.volume);
    struct gkyl_array* dualcurlBhatoverB_surface2 = mkarr(false, 3, nrange_surface.volume);
    struct gkyl_array* rtg33inv_surface2 = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* bioverJB_surface2 = mkarr(false, 3, nrange_surface.volume);
    struct gkyl_array* dualmags_surface2 = mkarr(false, 3, nrange_surface.volume);
    struct gkyl_array* dualnormals_surface2 = mkarr(false, 9, nrange_surface.volume);

    while (gkyl_range_iter_next(&iter_surface)) {
        // Get current node index
        long idx = gkyl_range_idx(&nrange_surface, iter_surface.idx);
        
        double *tang_data = (double*)tangents_surface2->data;
        double *bmag_data = (double*)Bmag_surface2->data;
        double *curl_data = (double*)curlBhat_surface2->data;
        double *dual_data = (double*)duals_surface2->data;
        double *gij_data = (double*)g_ij_surface2->data;
        double *jacob_data = (double*)jacobgeo_surface2->data;
        double *b3_data = (double*)B3_surface2->data;
        double *dcb_data = (double*)dualcurlBhat_surface2->data;
        double *dcbob_data = (double*)dualcurlBhatoverB_surface2->data;
        double *rtg_data = (double*)rtg33inv_surface2->data;
        double *biojb_data = (double*)bioverJB_surface2->data;
        double *dmag_data = (double*)dualmags_surface2->data;
        double *dnorm_data = (double*)dualnormals_surface2->data;
        
        double e_1[3] = {tang_data[idx*9+0], tang_data[idx*9+1], tang_data[idx*9+2]};
        double e_2[3] = {tang_data[idx*9+3], tang_data[idx*9+4], tang_data[idx*9+5]};
        double e_3[3] = {tang_data[idx*9+6], tang_data[idx*9+7], tang_data[idx*9+8]};
        double curlB[3] = {curl_data[idx*3+0], curl_data[idx*3+1], curl_data[idx*3+2]};
        double Bmag = bmag_data[idx];
        
        // Jacobian: J = e_1 · (e_2 × e_3)
        double e2xe3[3] = {e_2[1]*e_3[2]-e_2[2]*e_3[1], e_2[2]*e_3[0]-e_2[0]*e_3[2], e_2[0]*e_3[1]-e_2[1]*e_3[0]};
        double J = e_1[0]*e2xe3[0] + e_1[1]*e2xe3[1] + e_1[2]*e2xe3[2];
        jacob_data[idx] = J;
        if (fabs(J) < 1e-12) {
            // Jacobian is too small, skip this node or set to a small value
            J = 1e-12;
            jacob_data[idx] = J;
        }
        double Jinv = 1.0/J;
        
        // Dual vectors: e^i = (1/J) * (e_j × e_k)
        double e3xe1[3] = {e_3[1]*e_1[2]-e_3[2]*e_1[1], e_3[2]*e_1[0]-e_3[0]*e_1[2], e_3[0]*e_1[1]-e_3[1]*e_1[0]};
        double e1xe2[3] = {e_1[1]*e_2[2]-e_1[2]*e_2[1], e_1[2]*e_2[0]-e_1[0]*e_2[2], e_1[0]*e_2[1]-e_1[1]*e_2[0]};
        for (int i=0; i<3; i++) {
            dual_data[idx*9+0*3+i] = Jinv * e2xe3[i];  // e^1
            dual_data[idx*9+1*3+i] = Jinv * e3xe1[i];  // e^2
            dual_data[idx*9+2*3+i] = Jinv * e1xe2[i];  // e^3
        }
        
        // Metric tensor: g_ij = e_i · e_j
        for (int i=0; i<3; i++) for (int j=0; j<3; j++) {
            double *ei = (i==0)?e_1:(i==1)?e_2:e_3;
            double *ej = (j==0)?e_1:(j==1)?e_2:e_3;
            gij_data[idx*9+i*3+j] = ei[0]*ej[0] + ei[1]*ej[1] + ei[2]*ej[2];
        }
        
        double g33 = gij_data[idx*9+2*3+2];
        double sqrt_g33 = sqrt(g33);
        double rtg33inv = 1.0/sqrt_g33;
        rtg_data[idx] = rtg33inv;
        
        // B3 = Bmag / sqrt(g_33)
        b3_data[idx] = Bmag * rtg33inv;
        
        // dualcurlbhat = e^i · curlBhat and dualcurlBhatoverB
        for (int i=0; i<3; i++) {
            double *ei = &dual_data[idx*9+i*3];
            dcb_data[idx*3+i] = ei[0]*curlB[0] + ei[1]*curlB[1] + ei[2]*curlB[2];
            dcbob_data[idx*3+i] = dcb_data[idx*3+i] / Bmag;
        }
        
        // bioverJB = g_i3 / sqrt(g_33) / J / Bmag 
        for (int i=0; i<3; i++) {
            biojb_data[idx*3+i] = gij_data[idx*9+i*3+2] / sqrt_g33 / J / Bmag;
        }
        
        // Dual vector magnitudes: |e^i|
        for (int i=0; i<3; i++) {
            double *ei = &dual_data[idx*9+i*3];
            dmag_data[idx*3+i] = sqrt(ei[0]*ei[0] + ei[1]*ei[1] + ei[2]*ei[2]);
        }
        
        // Normal vectors: e^i / |e^i|
        for (int i=0; i<3; i++) {
            double mag = dmag_data[idx*3+i];
            double *ei = &dual_data[idx*9+i*3];
            if (mag > 1e-12) {  // Avoid division by zero
                for (int j=0; j<3; j++) dnorm_data[idx*9+i*3+j] = ei[j] / mag;
            } else {
                // If magnitude is zero, set normal to zero vector
                for (int j=0; j<3; j++) dnorm_data[idx*9+i*3+j] = 0.0;
            }
        }
    }

     // calculate quantities for surface 3
    struct gkyl_array* duals_surface3 = mkarr(false, 9, nrange_surface.volume);
    struct gkyl_array* jacobgeo_surface3 = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* g_ij_surface3 = mkarr(false, 9, nrange_surface.volume);
    struct gkyl_array* B3_surface3 = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* dualcurlBhat_surface3 = mkarr(false, 3, nrange_surface.volume);
    struct gkyl_array* dualcurlBhatoverB_surface3 = mkarr(false, 3, nrange_surface.volume);
    struct gkyl_array* rtg33inv_surface3 = mkarr(false, 1, nrange_surface.volume);
    struct gkyl_array* bioverJB_surface3 = mkarr(false, 3, nrange_surface.volume);
    struct gkyl_array* dualmags_surface3 = mkarr(false, 3, nrange_surface.volume);
    struct gkyl_array* dualnormals_surface3 = mkarr(false, 9, nrange_surface.volume);

    while (gkyl_range_iter_next(&iter_surface)) {
        // Get current node index
        long idx = gkyl_range_idx(&nrange_surface, iter_surface.idx);
        
        double *tang_data = (double*)tangents_surface3->data;
        double *bmag_data = (double*)Bmag_surface3->data;
        double *curl_data = (double*)curlBhat_surface3->data;
        double *dual_data = (double*)duals_surface3->data;
        double *gij_data = (double*)g_ij_surface3->data;
        double *jacob_data = (double*)jacobgeo_surface3->data;
        double *b3_data = (double*)B3_surface3->data;
        double *dcb_data = (double*)dualcurlBhat_surface3->data;
        double *dcbob_data = (double*)dualcurlBhatoverB_surface3->data;
        double *rtg_data = (double*)rtg33inv_surface3->data;
        double *biojb_data = (double*)bioverJB_surface3->data;
        double *dmag_data = (double*)dualmags_surface3->data;
        double *dnorm_data = (double*)dualnormals_surface3->data;
        
        double e_1[3] = {tang_data[idx*9+0], tang_data[idx*9+1], tang_data[idx*9+2]};
        double e_2[3] = {tang_data[idx*9+3], tang_data[idx*9+4], tang_data[idx*9+5]};
        double e_3[3] = {tang_data[idx*9+6], tang_data[idx*9+7], tang_data[idx*9+8]};
        double curlB[3] = {curl_data[idx*3+0], curl_data[idx*3+1], curl_data[idx*3+2]};
        double Bmag = bmag_data[idx];
        
        // Jacobian: J = e_1 · (e_2 × e_3)
        double e2xe3[3] = {e_2[1]*e_3[2]-e_2[2]*e_3[1], e_2[2]*e_3[0]-e_2[0]*e_3[2], e_2[0]*e_3[1]-e_2[1]*e_3[0]};
        double J = e_1[0]*e2xe3[0] + e_1[1]*e2xe3[1] + e_1[2]*e2xe3[2];
        jacob_data[idx] = J;
        if (fabs(J) < 1e-12) {
            // Jacobian is too small, skip this node or set to a small value
            J = 1e-12;
            jacob_data[idx] = J;
        }
        double Jinv = 1.0/J;
        
        // Dual vectors: e^i = (1/J) * (e_j × e_k)
        double e3xe1[3] = {e_3[1]*e_1[2]-e_3[2]*e_1[1], e_3[2]*e_1[0]-e_3[0]*e_1[2], e_3[0]*e_1[1]-e_3[1]*e_1[0]};
        double e1xe2[3] = {e_1[1]*e_2[2]-e_1[2]*e_2[1], e_1[2]*e_2[0]-e_1[0]*e_2[2], e_1[0]*e_2[1]-e_1[1]*e_2[0]};
        for (int i=0; i<3; i++) {
            dual_data[idx*9+0*3+i] = Jinv * e2xe3[i];  // e^1
            dual_data[idx*9+1*3+i] = Jinv * e3xe1[i];  // e^2
            dual_data[idx*9+2*3+i] = Jinv * e1xe2[i];  // e^3
        }
        
        // Metric tensor: g_ij = e_i · e_j
        for (int i=0; i<3; i++) for (int j=0; j<3; j++) {
            double *ei = (i==0)?e_1:(i==1)?e_2:e_3;
            double *ej = (j==0)?e_1:(j==1)?e_2:e_3;
            gij_data[idx*9+i*3+j] = ei[0]*ej[0] + ei[1]*ej[1] + ei[2]*ej[2];
        }
        
        double g33 = gij_data[idx*9+2*3+2];
        double sqrt_g33 = sqrt(g33);
        double rtg33inv = 1.0/sqrt_g33;
        rtg_data[idx] = rtg33inv;
        
        // B3 = Bmag / sqrt(g_33)
        b3_data[idx] = Bmag * rtg33inv;
        
        // dualcurlbhat = e^i · curlBhat and dualcurlBhatoverB
        for (int i=0; i<3; i++) {
            double *ei = &dual_data[idx*9+i*3];
            dcb_data[idx*3+i] = ei[0]*curlB[0] + ei[1]*curlB[1] + ei[2]*curlB[2];
            dcbob_data[idx*3+i] = dcb_data[idx*3+i] / Bmag;
        }
        
        // bioverJB = g_i3 / sqrt(g_33) / J / Bmag 
        for (int i=0; i<3; i++) {
            biojb_data[idx*3+i] = gij_data[idx*9+i*3+2] / sqrt_g33 / J / Bmag;
        }
        
        // Dual vector magnitudes: |e^i|
        for (int i=0; i<3; i++) {
            double *ei = &dual_data[idx*9+i*3];
            dmag_data[idx*3+i] = sqrt(ei[0]*ei[0] + ei[1]*ei[1] + ei[2]*ei[2]);
        }
        
        // Normal vectors: e^i / |e^i|
        for (int i=0; i<3; i++) {
            double mag = dmag_data[idx*3+i];
            double *ei = &dual_data[idx*9+i*3];
            if (mag > 1e-12) {  // Avoid division by zero
                for (int j=0; j<3; j++) dnorm_data[idx*9+i*3+j] = ei[j] / mag;
            } else {
                // If magnitude is zero, set normal to zero vector
                for (int j=0; j<3; j++) dnorm_data[idx*9+i*3+j] = 0.0;
            }
        }
    }

    // Now we have all of the following quantities as nodal arrays (75 arrays total) 
    //for corner, interior and three surfaces:
    // raz
    // rpz
    // tangents
    // duals
    // Bmag
    // curlBhat
    // jacobgeo
    // g_ij
    // B3
    // dualcurlBhat
    // dualcurlBhatoverB
    // rtg33inv
    // bioverJB
    // dualmags
    // dualnormals
    
    // Do a n2m (use appropriate for interior,corner, and surface) on all neede quantities
    // write out expansions on the modal grid
    
    // Get cell counts from grid header (we already know array sizes, but need 3D cell structure)
    // Read grid from one file to get cells[] - bounds don't matter, only cell counts
    // NOTE: All files must have the same cells[] and ndim (spatial dimensions) for n2m conversion.
    // The computational range (based on cells[]) must be the same for corner, interior, and surface.
    // Physical bounds (lower/upper) may differ between quantities (e.g., rpz vs raz coordinates).
    // The number of components per node (ncomp) can differ (e.g., 1 for Bmag, 9 for g_ij) and is
    // handled separately as a parameter to the n2m functions - it does not affect the grid structure.
    struct gkyl_rect_grid grid_orig;  // Original grid with proper bounds for writing
    struct gkyl_rect_grid grid;       // Grid with dummy bounds for n2m conversion
    struct gkyl_array_header_info hdr;
    gkyl_grid_sub_array_header_read(&grid_orig, &hdr, "W7-X_field-aligned_coords_gkyl/drpz_draz_corner.gkyl");
    
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
        modal = mkarr(false, ncomp, comp_range.volume * basis.num_basis); \
        gkyl_nodal_ops_n2m_surface(n2m, &basis, &grid, &nrange, &comp_range, ncomp, nodal, modal, dir);
    
    // Convert corner arrays (15 total) - use_quad=false for corners
    struct gkyl_array *raz_corner_modal, *rpz_corner_modal, *tangents_corner_modal, *duals_corner_modal;
    struct gkyl_array *Bmag_corner_modal, *curlBhat_corner_modal, *jacobgeo_corner_modal, *g_ij_corner_modal;
    struct gkyl_array *B3_corner_modal, *dualcurlBhat_corner_modal, *dualcurlBhatoverB_corner_modal;
    struct gkyl_array *rtg33inv_corner_modal, *bioverJB_corner_modal, *dualmags_corner_modal, *dualnormals_corner_modal;
    N2M_CONVERT(nrange_corner, comp_range, 3, raz_corner, raz_corner_modal);
    N2M_CONVERT(nrange_corner, comp_range, 3, rpz_corner, rpz_corner_modal);
    N2M_CONVERT(nrange_corner, comp_range, 9, tangents_corner, tangents_corner_modal);
    N2M_CONVERT(nrange_corner, comp_range, 9, duals_corner, duals_corner_modal);
    N2M_CONVERT(nrange_corner, comp_range, 1, Bmag_corner, Bmag_corner_modal);
    N2M_CONVERT(nrange_corner, comp_range, 3, curlBhat_corner, curlBhat_corner_modal);
    N2M_CONVERT(nrange_corner, comp_range, 1, jacobgeo_corner, jacobgeo_corner_modal);
    N2M_CONVERT(nrange_corner, comp_range, 9, g_ij_corner, g_ij_corner_modal);
    N2M_CONVERT(nrange_corner, comp_range, 1, B3_corner, B3_corner_modal);
    N2M_CONVERT(nrange_corner, comp_range, 3, dualcurlBhat_corner, dualcurlBhat_corner_modal);
    N2M_CONVERT(nrange_corner, comp_range, 3, dualcurlBhatoverB_corner, dualcurlBhatoverB_corner_modal);
    N2M_CONVERT(nrange_corner, comp_range, 1, rtg33inv_corner, rtg33inv_corner_modal);
    N2M_CONVERT(nrange_corner, comp_range, 3, bioverJB_corner, bioverJB_corner_modal);
    N2M_CONVERT(nrange_corner, comp_range, 3, dualmags_corner, dualmags_corner_modal);
    N2M_CONVERT(nrange_corner, comp_range, 9, dualnormals_corner, dualnormals_corner_modal);
    
    // Convert interior arrays (15 total) - use_quad=true for interior
    struct gkyl_array *raz_interior_modal, *rpz_interior_modal, *tangents_interior_modal, *duals_interior_modal;
    struct gkyl_array *Bmag_interior_modal, *curlBhat_interior_modal, *jacobgeo_interior_modal, *g_ij_interior_modal;
    struct gkyl_array *B3_interior_modal, *dualcurlBhat_interior_modal, *dualcurlBhatoverB_interior_modal;
    struct gkyl_array *rtg33inv_interior_modal, *bioverJB_interior_modal, *dualmags_interior_modal, *dualnormals_interior_modal;
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 3, raz_interior, raz_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 3, rpz_interior, rpz_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 9, tangents_interior, tangents_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 9, duals_interior, duals_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 1, Bmag_interior, Bmag_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 3, curlBhat_interior, curlBhat_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 1, jacobgeo_interior, jacobgeo_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 9, g_ij_interior, g_ij_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 1, B3_interior, B3_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 3, dualcurlBhat_interior, dualcurlBhat_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 3, dualcurlBhatoverB_interior, dualcurlBhatoverB_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 1, rtg33inv_interior, rtg33inv_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 3, bioverJB_interior, bioverJB_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 3, dualmags_interior, dualmags_interior_modal);
    N2M_CONVERT_INTERIOR(nrange_interior, comp_range, 9, dualnormals_interior, dualnormals_interior_modal);
    
    // Convert surface arrays (45 total = 15 per surface × 3 surfaces) - use_quad=false for surfaces
    // Surface 1 arrays (15 total) - use surface n2m with direction 0
    struct gkyl_array *raz_s1_modal, *rpz_s1_modal, *tangents_s1_modal, *duals_s1_modal;
    struct gkyl_array *Bmag_s1_modal, *curlBhat_s1_modal, *jacobgeo_s1_modal, *g_ij_s1_modal;
    struct gkyl_array *B3_s1_modal, *dualcurlBhat_s1_modal, *dualcurlBhatoverB_s1_modal;
    struct gkyl_array *rtg33inv_s1_modal, *bioverJB_s1_modal, *dualmags_s1_modal, *dualnormals_s1_modal;
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, raz_surface1, raz_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, rpz_surface1, rpz_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 9, tangents_surface1, tangents_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 9, duals_surface1, duals_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, Bmag_surface1, Bmag_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, curlBhat_surface1, curlBhat_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, jacobgeo_surface1, jacobgeo_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 9, g_ij_surface1, g_ij_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, B3_surface1, B3_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, dualcurlBhat_surface1, dualcurlBhat_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, dualcurlBhatoverB_surface1, dualcurlBhatoverB_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, rtg33inv_surface1, rtg33inv_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, bioverJB_surface1, bioverJB_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, dualmags_surface1, dualmags_s1_modal, 0);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 9, dualnormals_surface1, dualnormals_s1_modal, 0);
    
    // Surface 2 arrays (15 total) - use surface n2m with direction 1
    struct gkyl_array *raz_s2_modal, *rpz_s2_modal, *tangents_s2_modal, *duals_s2_modal;
    struct gkyl_array *Bmag_s2_modal, *curlBhat_s2_modal, *jacobgeo_s2_modal, *g_ij_s2_modal;
    struct gkyl_array *B3_s2_modal, *dualcurlBhat_s2_modal, *dualcurlBhatoverB_s2_modal;
    struct gkyl_array *rtg33inv_s2_modal, *bioverJB_s2_modal, *dualmags_s2_modal, *dualnormals_s2_modal;
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, raz_surface2, raz_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, rpz_surface2, rpz_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 9, tangents_surface2, tangents_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 9, duals_surface2, duals_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, Bmag_surface2, Bmag_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, curlBhat_surface2, curlBhat_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, jacobgeo_surface2, jacobgeo_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 9, g_ij_surface2, g_ij_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, B3_surface2, B3_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, dualcurlBhat_surface2, dualcurlBhat_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, dualcurlBhatoverB_surface2, dualcurlBhatoverB_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, rtg33inv_surface2, rtg33inv_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, bioverJB_surface2, bioverJB_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, dualmags_surface2, dualmags_s2_modal, 1);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 9, dualnormals_surface2, dualnormals_s2_modal, 1);
    
    // Surface 3 arrays (15 total) - use surface n2m with direction 2
    struct gkyl_array *raz_s3_modal, *rpz_s3_modal, *tangents_s3_modal, *duals_s3_modal;
    struct gkyl_array *Bmag_s3_modal, *curlBhat_s3_modal, *jacobgeo_s3_modal, *g_ij_s3_modal;
    struct gkyl_array *B3_s3_modal, *dualcurlBhat_s3_modal, *dualcurlBhatoverB_s3_modal;
    struct gkyl_array *rtg33inv_s3_modal, *bioverJB_s3_modal, *dualmags_s3_modal, *dualnormals_s3_modal;
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, raz_surface3, raz_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, rpz_surface3, rpz_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 9, tangents_surface3, tangents_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 9, duals_surface3, duals_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, Bmag_surface3, Bmag_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, curlBhat_surface3, curlBhat_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, jacobgeo_surface3, jacobgeo_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 9, g_ij_surface3, g_ij_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, B3_surface3, B3_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, dualcurlBhat_surface3, dualcurlBhat_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, dualcurlBhatoverB_surface3, dualcurlBhatoverB_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 1, rtg33inv_surface3, rtg33inv_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, bioverJB_surface3, bioverJB_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 3, dualmags_surface3, dualmags_s3_modal, 2);
    N2M_CONVERT_SURFACE(nrange_surface, comp_range, 9, dualnormals_surface3, dualnormals_s3_modal, 2);
    

    // Write out all nodal and modal arrays 
    // Create output directories
    const char *output_dir = "W7X-nodal_modal_arrays";
    const char *nodal_dir = "W7X-nodal_modal_arrays/nodal";
    const char *modal_dir = "W7X-nodal_modal_arrays/modal";
    #ifdef _WIN32
        _mkdir(output_dir);
        _mkdir(nodal_dir);
        _mkdir(modal_dir);
    #else
        mkdir(output_dir, 0755);
        mkdir(nodal_dir, 0755);
        mkdir(modal_dir, 0755);
    #endif
    
    // Helper macro to write nodal array
    #define WRITE_NODAL(ngrid, nrange, arr, name) \
        do { \
            char fname[512]; \
            snprintf(fname, sizeof(fname), "%s/%s.gkyl", nodal_dir, name); \
            gkyl_grid_sub_array_write(ngrid, nrange, 0, arr, fname); \
        } while(0)
    
    // Helper macro to write modal array
    #define WRITE_MODAL(grid, comp_range, arr, name) \
        do { \
            char fname[512]; \
            snprintf(fname, sizeof(fname), "%s/%s.gkyl", modal_dir, name); \
            gkyl_grid_sub_array_write(grid, comp_range, 0, arr, fname); \
        } while(0)
    
    // Create nodal grids for writing nodal arrays
    struct gkyl_rect_grid ngrid_corner, ngrid_interior, ngrid_surface;
    gkyl_gk_geometry_init_nodal_grid(&ngrid_corner, &grid_orig, &nrange_corner);
    gkyl_gk_geometry_init_nodal_grid(&ngrid_interior, &grid_orig, &nrange_interior);
    gkyl_gk_geometry_init_nodal_grid(&ngrid_surface, &grid_orig, &nrange_surface);
    
    // Write corner nodal arrays (15 total)
    WRITE_NODAL(&ngrid_corner, &nrange_corner, raz_corner, "raz_corner");
    WRITE_NODAL(&ngrid_corner, &nrange_corner, rpz_corner, "rpz_corner");
    WRITE_NODAL(&ngrid_corner, &nrange_corner, tangents_corner, "tangents_corner");
    WRITE_NODAL(&ngrid_corner, &nrange_corner, duals_corner, "duals_corner");
    WRITE_NODAL(&ngrid_corner, &nrange_corner, Bmag_corner, "Bmag_corner");
    WRITE_NODAL(&ngrid_corner, &nrange_corner, curlBhat_corner, "curlBhat_corner");
    WRITE_NODAL(&ngrid_corner, &nrange_corner, jacobgeo_corner, "jacobgeo_corner");
    WRITE_NODAL(&ngrid_corner, &nrange_corner, g_ij_corner, "g_ij_corner");
    WRITE_NODAL(&ngrid_corner, &nrange_corner, B3_corner, "B3_corner");
    WRITE_NODAL(&ngrid_corner, &nrange_corner, dualcurlBhat_corner, "dualcurlBhat_corner");
    WRITE_NODAL(&ngrid_corner, &nrange_corner, dualcurlBhatoverB_corner, "dualcurlBhatoverB_corner");
    WRITE_NODAL(&ngrid_corner, &nrange_corner, rtg33inv_corner, "rtg33inv_corner");
    WRITE_NODAL(&ngrid_corner, &nrange_corner, bioverJB_corner, "bioverJB_corner");
    WRITE_NODAL(&ngrid_corner, &nrange_corner, dualmags_corner, "dualmags_corner");
    WRITE_NODAL(&ngrid_corner, &nrange_corner, dualnormals_corner, "dualnormals_corner");
    
    // Write interior nodal arrays (15 total)
    WRITE_NODAL(&ngrid_interior, &nrange_interior, raz_interior, "raz_interior");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, rpz_interior, "rpz_interior");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, tangents_interior, "tangents_interior");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, duals_interior, "duals_interior");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, Bmag_interior, "Bmag_interior");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, curlBhat_interior, "curlBhat_interior");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, jacobgeo_interior, "jacobgeo_interior");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, g_ij_interior, "g_ij_interior");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, B3_interior, "B3_interior");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, dualcurlBhat_interior, "dualcurlBhat_interior");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, dualcurlBhatoverB_interior, "dualcurlBhatoverB_interior");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, rtg33inv_interior, "rtg33inv_interior");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, bioverJB_interior, "bioverJB_interior");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, dualmags_interior, "dualmags_interior");
    WRITE_NODAL(&ngrid_interior, &nrange_interior, dualnormals_interior, "dualnormals_interior");
    
    // Write surface 1 nodal arrays (15 total)
    WRITE_NODAL(&ngrid_surface, &nrange_surface, raz_surface1, "raz_surface1");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, rpz_surface1, "rpz_surface1");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, tangents_surface1, "tangents_surface1");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, duals_surface1, "duals_surface1");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, Bmag_surface1, "Bmag_surface1");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, curlBhat_surface1, "curlBhat_surface1");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, jacobgeo_surface1, "jacobgeo_surface1");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, g_ij_surface1, "g_ij_surface1");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, B3_surface1, "B3_surface1");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, dualcurlBhat_surface1, "dualcurlBhat_surface1");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, dualcurlBhatoverB_surface1, "dualcurlBhatoverB_surface1");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, rtg33inv_surface1, "rtg33inv_surface1");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, bioverJB_surface1, "bioverJB_surface1");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, dualmags_surface1, "dualmags_surface1");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, dualnormals_surface1, "dualnormals_surface1");
    
    // Write surface 2 nodal arrays (15 total)
    WRITE_NODAL(&ngrid_surface, &nrange_surface, raz_surface2, "raz_surface2");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, rpz_surface2, "rpz_surface2");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, tangents_surface2, "tangents_surface2");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, duals_surface2, "duals_surface2");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, Bmag_surface2, "Bmag_surface2");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, curlBhat_surface2, "curlBhat_surface2");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, jacobgeo_surface2, "jacobgeo_surface2");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, g_ij_surface2, "g_ij_surface2");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, B3_surface2, "B3_surface2");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, dualcurlBhat_surface2, "dualcurlBhat_surface2");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, dualcurlBhatoverB_surface2, "dualcurlBhatoverB_surface2");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, rtg33inv_surface2, "rtg33inv_surface2");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, bioverJB_surface2, "bioverJB_surface2");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, dualmags_surface2, "dualmags_surface2");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, dualnormals_surface2, "dualnormals_surface2");
    
    // Write surface 3 nodal arrays (15 total)
    WRITE_NODAL(&ngrid_surface, &nrange_surface, raz_surface3, "raz_surface3");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, rpz_surface3, "rpz_surface3");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, tangents_surface3, "tangents_surface3");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, duals_surface3, "duals_surface3");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, Bmag_surface3, "Bmag_surface3");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, curlBhat_surface3, "curlBhat_surface3");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, jacobgeo_surface3, "jacobgeo_surface3");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, g_ij_surface3, "g_ij_surface3");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, B3_surface3, "B3_surface3");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, dualcurlBhat_surface3, "dualcurlBhat_surface3");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, dualcurlBhatoverB_surface3, "dualcurlBhatoverB_surface3");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, rtg33inv_surface3, "rtg33inv_surface3");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, bioverJB_surface3, "bioverJB_surface3");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, dualmags_surface3, "dualmags_surface3");
    WRITE_NODAL(&ngrid_surface, &nrange_surface, dualnormals_surface3, "dualnormals_surface3");
    
    // Write corner modal arrays (15 total)
    WRITE_MODAL(&grid_orig, &comp_range, raz_corner_modal, "raz_corner_modal");
    WRITE_MODAL(&grid_orig, &comp_range, rpz_corner_modal, "rpz_corner_modal");
    WRITE_MODAL(&grid_orig, &comp_range, tangents_corner_modal, "tangents_corner_modal");
    WRITE_MODAL(&grid_orig, &comp_range, duals_corner_modal, "duals_corner_modal");
    WRITE_MODAL(&grid_orig, &comp_range, Bmag_corner_modal, "Bmag_corner_modal");
    WRITE_MODAL(&grid_orig, &comp_range, curlBhat_corner_modal, "curlBhat_corner_modal");
    WRITE_MODAL(&grid_orig, &comp_range, jacobgeo_corner_modal, "jacobgeo_corner_modal");
    WRITE_MODAL(&grid_orig, &comp_range, g_ij_corner_modal, "g_ij_corner_modal");
    WRITE_MODAL(&grid_orig, &comp_range, B3_corner_modal, "B3_corner_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualcurlBhat_corner_modal, "dualcurlBhat_corner_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualcurlBhatoverB_corner_modal, "dualcurlBhatoverB_corner_modal");
    WRITE_MODAL(&grid_orig, &comp_range, rtg33inv_corner_modal, "rtg33inv_corner_modal");
    WRITE_MODAL(&grid_orig, &comp_range, bioverJB_corner_modal, "bioverJB_corner_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualmags_corner_modal, "dualmags_corner_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualnormals_corner_modal, "dualnormals_corner_modal");
    
    // Write interior modal arrays (15 total)
    WRITE_MODAL(&grid_orig, &comp_range, raz_interior_modal, "raz_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, rpz_interior_modal, "rpz_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, tangents_interior_modal, "tangents_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, duals_interior_modal, "duals_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, Bmag_interior_modal, "Bmag_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, curlBhat_interior_modal, "curlBhat_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, jacobgeo_interior_modal, "jacobgeo_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, g_ij_interior_modal, "g_ij_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, B3_interior_modal, "B3_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualcurlBhat_interior_modal, "dualcurlBhat_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualcurlBhatoverB_interior_modal, "dualcurlBhatoverB_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, rtg33inv_interior_modal, "rtg33inv_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, bioverJB_interior_modal, "bioverJB_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualmags_interior_modal, "dualmags_interior_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualnormals_interior_modal, "dualnormals_interior_modal");
    
    // Write surface 1 modal arrays (15 total)
    WRITE_MODAL(&grid_orig, &comp_range, raz_s1_modal, "raz_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, rpz_s1_modal, "rpz_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, tangents_s1_modal, "tangents_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, duals_s1_modal, "duals_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, Bmag_s1_modal, "Bmag_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, curlBhat_s1_modal, "curlBhat_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, jacobgeo_s1_modal, "jacobgeo_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, g_ij_s1_modal, "g_ij_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, B3_s1_modal, "B3_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualcurlBhat_s1_modal, "dualcurlBhat_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualcurlBhatoverB_s1_modal, "dualcurlBhatoverB_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, rtg33inv_s1_modal, "rtg33inv_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, bioverJB_s1_modal, "bioverJB_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualmags_s1_modal, "dualmags_s1_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualnormals_s1_modal, "dualnormals_s1_modal");
    
    // Write surface 2 modal arrays (15 total)
    WRITE_MODAL(&grid_orig, &comp_range, raz_s2_modal, "raz_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, rpz_s2_modal, "rpz_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, tangents_s2_modal, "tangents_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, duals_s2_modal, "duals_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, Bmag_s2_modal, "Bmag_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, curlBhat_s2_modal, "curlBhat_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, jacobgeo_s2_modal, "jacobgeo_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, g_ij_s2_modal, "g_ij_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, B3_s2_modal, "B3_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualcurlBhat_s2_modal, "dualcurlBhat_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualcurlBhatoverB_s2_modal, "dualcurlBhatoverB_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, rtg33inv_s2_modal, "rtg33inv_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, bioverJB_s2_modal, "bioverJB_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualmags_s2_modal, "dualmags_s2_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualnormals_s2_modal, "dualnormals_s2_modal");
    
    // Write surface 3 modal arrays (15 total)
    WRITE_MODAL(&grid_orig, &comp_range, raz_s3_modal, "raz_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, rpz_s3_modal, "rpz_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, tangents_s3_modal, "tangents_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, duals_s3_modal, "duals_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, Bmag_s3_modal, "Bmag_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, curlBhat_s3_modal, "curlBhat_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, jacobgeo_s3_modal, "jacobgeo_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, g_ij_s3_modal, "g_ij_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, B3_s3_modal, "B3_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualcurlBhat_s3_modal, "dualcurlBhat_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualcurlBhatoverB_s3_modal, "dualcurlBhatoverB_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, rtg33inv_s3_modal, "rtg33inv_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, bioverJB_s3_modal, "bioverJB_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualmags_s3_modal, "dualmags_s3_modal");
    WRITE_MODAL(&grid_orig, &comp_range, dualnormals_s3_modal, "dualnormals_s3_modal");

    // Clean up everything
    // Release corner loaded arrays (5)
    gkyl_array_release(tangents_corner);
    gkyl_array_release(Bmag_corner);
    gkyl_array_release(curlBhat_corner);
    gkyl_array_release(raz_corner);
    gkyl_array_release(rpz_corner);
    
    // Release corner derived arrays (10)
    gkyl_array_release(duals_corner);
    gkyl_array_release(jacobgeo_corner);
    gkyl_array_release(g_ij_corner);
    gkyl_array_release(B3_corner);
    gkyl_array_release(dualcurlBhat_corner);
    gkyl_array_release(dualcurlBhatoverB_corner);
    gkyl_array_release(rtg33inv_corner);
    gkyl_array_release(bioverJB_corner);
    gkyl_array_release(dualmags_corner);
    gkyl_array_release(dualnormals_corner);
    
    // Release corner modal arrays (15)
    gkyl_array_release(raz_corner_modal);
    gkyl_array_release(rpz_corner_modal);
    gkyl_array_release(tangents_corner_modal);
    gkyl_array_release(duals_corner_modal);
    gkyl_array_release(Bmag_corner_modal);
    gkyl_array_release(curlBhat_corner_modal);
    gkyl_array_release(jacobgeo_corner_modal);
    gkyl_array_release(g_ij_corner_modal);
    gkyl_array_release(B3_corner_modal);
    gkyl_array_release(dualcurlBhat_corner_modal);
    gkyl_array_release(dualcurlBhatoverB_corner_modal);
    gkyl_array_release(rtg33inv_corner_modal);
    gkyl_array_release(bioverJB_corner_modal);
    gkyl_array_release(dualmags_corner_modal);
    gkyl_array_release(dualnormals_corner_modal);
    
    // Release interior loaded arrays (5)
    gkyl_array_release(tangents_interior);
    gkyl_array_release(Bmag_interior);
    gkyl_array_release(curlBhat_interior);
    gkyl_array_release(raz_interior);
    gkyl_array_release(rpz_interior);
    
    // Release interior derived arrays (10)
    gkyl_array_release(duals_interior);
    gkyl_array_release(jacobgeo_interior);
    gkyl_array_release(g_ij_interior);
    gkyl_array_release(B3_interior);
    gkyl_array_release(dualcurlBhat_interior);
    gkyl_array_release(dualcurlBhatoverB_interior);
    gkyl_array_release(rtg33inv_interior);
    gkyl_array_release(bioverJB_interior);
    gkyl_array_release(dualmags_interior);
    gkyl_array_release(dualnormals_interior);
    
    // Release interior modal arrays (15)
    gkyl_array_release(raz_interior_modal);
    gkyl_array_release(rpz_interior_modal);
    gkyl_array_release(tangents_interior_modal);
    gkyl_array_release(duals_interior_modal);
    gkyl_array_release(Bmag_interior_modal);
    gkyl_array_release(curlBhat_interior_modal);
    gkyl_array_release(jacobgeo_interior_modal);
    gkyl_array_release(g_ij_interior_modal);
    gkyl_array_release(B3_interior_modal);
    gkyl_array_release(dualcurlBhat_interior_modal);
    gkyl_array_release(dualcurlBhatoverB_interior_modal);
    gkyl_array_release(rtg33inv_interior_modal);
    gkyl_array_release(bioverJB_interior_modal);
    gkyl_array_release(dualmags_interior_modal);
    gkyl_array_release(dualnormals_interior_modal);
    
    // Release surface 1 loaded arrays (5)
    gkyl_array_release(tangents_surface1);
    gkyl_array_release(Bmag_surface1);
    gkyl_array_release(curlBhat_surface1);
    gkyl_array_release(raz_surface1);
    gkyl_array_release(rpz_surface1);
    
    // Release surface 1 derived arrays (10)
    gkyl_array_release(duals_surface1);
    gkyl_array_release(jacobgeo_surface1);
    gkyl_array_release(g_ij_surface1);
    gkyl_array_release(B3_surface1);
    gkyl_array_release(dualcurlBhat_surface1);
    gkyl_array_release(dualcurlBhatoverB_surface1);
    gkyl_array_release(rtg33inv_surface1);
    gkyl_array_release(bioverJB_surface1);
    gkyl_array_release(dualmags_surface1);
    gkyl_array_release(dualnormals_surface1);
    
    // Release surface 1 modal arrays (15)
    gkyl_array_release(raz_s1_modal);
    gkyl_array_release(rpz_s1_modal);
    gkyl_array_release(tangents_s1_modal);
    gkyl_array_release(duals_s1_modal);
    gkyl_array_release(Bmag_s1_modal);
    gkyl_array_release(curlBhat_s1_modal);
    gkyl_array_release(jacobgeo_s1_modal);
    gkyl_array_release(g_ij_s1_modal);
    gkyl_array_release(B3_s1_modal);
    gkyl_array_release(dualcurlBhat_s1_modal);
    gkyl_array_release(dualcurlBhatoverB_s1_modal);
    gkyl_array_release(rtg33inv_s1_modal);
    gkyl_array_release(bioverJB_s1_modal);
    gkyl_array_release(dualmags_s1_modal);
    gkyl_array_release(dualnormals_s1_modal);
    
    // Release surface 2 loaded arrays (5)
    gkyl_array_release(tangents_surface2);
    gkyl_array_release(Bmag_surface2);
    gkyl_array_release(curlBhat_surface2);
    gkyl_array_release(raz_surface2);
    gkyl_array_release(rpz_surface2);
    
    // Release surface 2 derived arrays (10)
    gkyl_array_release(duals_surface2);
    gkyl_array_release(jacobgeo_surface2);
    gkyl_array_release(g_ij_surface2);
    gkyl_array_release(B3_surface2);
    gkyl_array_release(dualcurlBhat_surface2);
    gkyl_array_release(dualcurlBhatoverB_surface2);
    gkyl_array_release(rtg33inv_surface2);
    gkyl_array_release(bioverJB_surface2);
    gkyl_array_release(dualmags_surface2);
    gkyl_array_release(dualnormals_surface2);
    
    // Release surface 2 modal arrays (15)
    gkyl_array_release(raz_s2_modal);
    gkyl_array_release(rpz_s2_modal);
    gkyl_array_release(tangents_s2_modal);
    gkyl_array_release(duals_s2_modal);
    gkyl_array_release(Bmag_s2_modal);
    gkyl_array_release(curlBhat_s2_modal);
    gkyl_array_release(jacobgeo_s2_modal);
    gkyl_array_release(g_ij_s2_modal);
    gkyl_array_release(B3_s2_modal);
    gkyl_array_release(dualcurlBhat_s2_modal);
    gkyl_array_release(dualcurlBhatoverB_s2_modal);
    gkyl_array_release(rtg33inv_s2_modal);
    gkyl_array_release(bioverJB_s2_modal);
    gkyl_array_release(dualmags_s2_modal);
    gkyl_array_release(dualnormals_s2_modal);
    
    // Release surface 3 loaded arrays (5)
    gkyl_array_release(tangents_surface3);
    gkyl_array_release(Bmag_surface3);
    gkyl_array_release(curlBhat_surface3);
    gkyl_array_release(raz_surface3);
    gkyl_array_release(rpz_surface3);
    
    // Release surface 3 derived arrays (10)
    gkyl_array_release(duals_surface3);
    gkyl_array_release(jacobgeo_surface3);
    gkyl_array_release(g_ij_surface3);
    gkyl_array_release(B3_surface3);
    gkyl_array_release(dualcurlBhat_surface3);
    gkyl_array_release(dualcurlBhatoverB_surface3);
    gkyl_array_release(rtg33inv_surface3);
    gkyl_array_release(bioverJB_surface3);
    gkyl_array_release(dualmags_surface3);
    gkyl_array_release(dualnormals_surface3);
    
    // Release surface 3 modal arrays (15)
    gkyl_array_release(raz_s3_modal);
    gkyl_array_release(rpz_s3_modal);
    gkyl_array_release(tangents_s3_modal);
    gkyl_array_release(duals_s3_modal);
    gkyl_array_release(Bmag_s3_modal);
    gkyl_array_release(curlBhat_s3_modal);
    gkyl_array_release(jacobgeo_s3_modal);
    gkyl_array_release(g_ij_s3_modal);
    gkyl_array_release(B3_s3_modal);
    gkyl_array_release(dualcurlBhat_s3_modal);
    gkyl_array_release(dualcurlBhatoverB_s3_modal);
    gkyl_array_release(rtg33inv_s3_modal);
    gkyl_array_release(bioverJB_s3_modal);
    gkyl_array_release(dualmags_s3_modal);
    gkyl_array_release(dualnormals_s3_modal);
    
    // Release other resources
    gkyl_array_header_info_release(&hdr);
    gkyl_nodal_ops_release(n2m);
}

TEST_LIST = {
    { "test_3x_p1", test_3xp1},
    { NULL, NULL },
  };



