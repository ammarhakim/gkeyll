#!/usr/bin/env python3
"""
Convert all .npy (NumPy array) files in a folder to .gkyl files using postgkyl.

Usage:
    python numpy-to-gkyl.py <folder_path>
"""

import os
import sys
import glob
import argparse
import numpy as np
import postgkyl as pg


def create_grid_from_shape(shape, lower=None, upper=None):
    """Create a grid (list of coordinate arrays) from array shape.
    
    Args:
        shape: Tuple of array dimensions (excluding components)
        lower: Lower bounds for each dimension (default: zeros)
        upper: Upper bounds for each dimension (default: ones)
    
    Returns:
        List of numpy arrays, one per spatial dimension
    """
    num_dims = len(shape)
    
    if lower is None:
        lower = np.zeros(num_dims)
    if upper is None:
        upper = np.ones(num_dims)
    
    grid = []
    for d in range(num_dims):
        # Create coordinate array for this dimension
        coords = np.linspace(lower[d], upper[d], shape[d])
        grid.append(coords)
    
    return grid


def convert_np_to_gkyl(inp_dir, filename, N, dims, lower, upper, output_dir=None, verbose=False, verify=False):
    """Convert a single .np file to .gkyl format.
    
    Args:
        inp_dir: Input directory path (e.g., "PATH/TO/FOLDER/")
        filename: Name of the .npy file (e.g., "file.npy")
        N: Total number of nodes (should equal dims[0] * dims[1] * dims[2])
        dims: List of dimensions [nx, ny, nz] where:
              nx: Number of nodes in x dimension (varies slowest)
              ny: Number of nodes in y dimension (varies middle)
              nz: Number of nodes in z dimension (varies fastest)
        lower: Lower bounds for grid [nx_lower, ny_lower, nz_lower]
        upper: Upper bounds for grid [nx_upper, ny_upper, nz_upper]
        output_dir: Directory to write output .gkyl file (default: same directory as input)
        verbose: Print detailed shape information
        verify: Read back the written file to verify correctness
    
    Returns:
        Path to output file
    """
    # Combine input directory and filename to get full path
    np_file_path = os.path.join(inp_dir, filename)
    
    # Validate inputs
    if len(dims) != 3:
        raise ValueError(f"dims must have 3 elements (got {len(dims)})")
    if len(lower) != 3:
        raise ValueError(f"lower must have 3 elements (got {len(lower)})")
    if len(upper) != 3:
        raise ValueError(f"upper must have 3 elements (got {len(upper)})")
    
    nx, ny, nz = dims[0], dims[1], dims[2]
    expected_N = nx * ny * nz
    if N != expected_N:
        raise ValueError(
            f"N ({N}) doesn't match dims[0] * dims[1] * dims[2] ({nx} * {ny} * {nz} = {expected_N})"
        )
    
    # Load numpy array
    arr = np.load(np_file_path)
    original_shape = arr.shape
    original_dtype = arr.dtype
    
    # Always print numpy array shape
    print(f"\n  File: {filename}")
    print(f"  NumPy array shape: {original_shape}, dtype: {original_dtype}")
    
    if verbose:
        print(f"  Full path: {np_file_path}")
        print(f"  Expected: N={N}, nx={nx}, ny={ny}, nz={nz}")
    
    # Handle structured arrays
    if arr.dtype.names is not None:
        # Structured array - take first field
        arr = arr[arr.dtype.names[0]]
        if verbose:
            print(f"  Structured array, using field: {arr.dtype.names[0]}")
    
    # Handle both 1D (N,) and 2D (N, Z) arrays
    # If 1D, interpret as (N, 1) with 1 component
    if arr.ndim == 1:
        # 1D array: interpret as (N, 1) with 1 component
        arr_N = arr.shape[0]
        num_comps = 1
        if arr_N != N:
            raise ValueError(
                f"Array size ({arr_N}) doesn't match expected N ({N})"
            )
        # Reshape to (N, 1) for consistent handling
        arr = arr.reshape(N, 1)
        if verbose:
            print(f"  Interpreted 1D array as (N, 1) with 1 component")
    elif arr.ndim == 2:
        # 2D array: shape should be (N, Z)
        arr_N, num_comps = arr.shape
        
        if arr_N != N:
            raise ValueError(
                f"Array first dimension ({arr_N}) doesn't match expected N ({N})"
            )
        
        if num_comps not in [1, 3, 9]:
            if verbose:
                print(f"  Warning: num_comps={num_comps} is not 1, 3, or 9 (but continuing anyway)")
    else:
        raise ValueError(
            f"Expected 1D array (N,) or 2D array (N, Z), got shape {arr.shape} "
            f"(ndim={arr.ndim})"
        )
    
    if verbose:
        print(f"  Array has {num_comps} components per node")
    
    # Reshape from (N, Z) to (nx, ny, nz, Z)
    # The flattening order is: nz varies fastest, then ny, then nx varies slowest
    # This corresponds to C-order (row-major) where last dimension varies fastest
    # So reshaping (N, Z) -> (nx, ny, nz, Z) is correct
    try:
        values = arr.reshape(nx, ny, nz, num_comps)
    except ValueError as e:
        raise ValueError(
            f"Cannot reshape array from shape {arr.shape} to ({nx}, {ny}, {nz}, {num_comps}). "
            f"Expected total size: {nx * ny * nz * num_comps}, got: {arr.size}. "
            f"Original array shape: {original_shape}"
        ) from e
    
    spatial_shape = (nx, ny, nz)
    
    # Always print shape information for debugging
    print(f"  Reshaped from {arr.shape} to {values.shape}")
    print(f"  Spatial shape: {spatial_shape}, components: {num_comps}")
    print(f"  Total values: {values.size} elements")
    
    if verbose:
        print(f"  Original array shape: {original_shape}")
    
    # Create grid from spatial shape with provided bounds
    grid = create_grid_from_shape(spatial_shape, lower, upper)
    
    # Validate grid dimensions match values dimensions
    if len(grid) != len(spatial_shape):
        raise ValueError(f"Grid has {len(grid)} dimensions but values have {len(spatial_shape)} spatial dimensions")
    
    for d, (grid_size, val_size) in enumerate(zip([len(g) for g in grid], spatial_shape)):
        if grid_size != val_size:
            raise ValueError(f"Grid dimension {d} size ({grid_size}) doesn't match values dimension {d} size ({val_size})")
    
    # Create GData object
    gdata = pg.GData(load=False)
    gdata.push(grid, values)
    
    # Explicitly set cells to match values shape (gkeyll uses cells, not nodes)
    # The values shape spatial dimensions represent the number of cells
    gdata.ctx["cells"] = np.array(spatial_shape, dtype=np.int32)
    gdata.ctx["num_comps"] = num_comps
    
    # Verify GData shape matches expectations
    gdata_values = gdata.get_values()
    gdata_grid = gdata.get_grid()
    
    # Always print gkyl array shape
    print(f"  Gkyl array shape: {gdata_values.shape}")
    print(f"  Gkyl dimensions: {gdata.num_dims}D, components: {gdata.num_comps}, cells: {gdata.num_cells}")
    
    if verbose:
        print(f"  Grid dimensions: {[len(g) for g in grid]}")
        print(f"  Grid bounds: lower={[g.min() for g in grid]}, upper={[g.max() for g in grid]}")
    
    # Validate GData internal consistency
    if gdata_values.shape[:-1] != tuple(spatial_shape):
        raise ValueError(f"GData values spatial shape {gdata_values.shape[:-1]} doesn't match expected {spatial_shape}")
    if gdata_values.shape[-1] != num_comps:
        raise ValueError(f"GData values components {gdata_values.shape[-1]} doesn't match expected {num_comps}")
    if len(gdata_grid) != len(spatial_shape):
        raise ValueError(f"GData grid dimensions {len(gdata_grid)} doesn't match spatial dimensions {len(spatial_shape)}")
    
    # Determine output filename
    if output_dir is None:
        # Default: same directory as input
        base_name = os.path.splitext(np_file_path)[0]
        output_path = base_name + '.gkyl'
    else:
        # Use specified output directory
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        base_name = os.path.splitext(filename)[0]
        output_path = os.path.join(output_dir, base_name + '.gkyl')
    
    # Write to .gkyl file manually (postgkyl's write has a bug with len(values) for multi-dim arrays)
    # We'll write it correctly ourselves
    dti = np.dtype("i8")
    dtf = np.dtype("f8")
    
    fh = open(output_path, "wb")  # Use binary mode
    
    # Write header according to gkyl format version 1
    # Magic bytes "gkyl0"
    np.array([103, 107, 121, 108, 48], dtype=np.dtype("b")).tofile(fh)
    # version 1
    np.array([1], dtype=dti).tofile(fh)
    # file_type 1 (field)
    np.array([1], dtype=dti).tofile(fh)
    # meta_size (0 for no metadata)
    np.array([0], dtype=dti).tofile(fh)
    # real_type (2 = double)
    np.array([2], dtype=dti).tofile(fh)
    # num_dims
    np.array([len(spatial_shape)], dtype=dti).tofile(fh)
    # num_cells (per dimension)
    np.array(spatial_shape, dtype=dti).tofile(fh)
    # lower bounds
    np.array(lower, dtype=dtf).tofile(fh)
    # upper bounds
    np.array(upper, dtype=dtf).tofile(fh)
    # elem_sz (element size * num components, 8 bytes per double)
    np.array([num_comps * 8], dtype=dti).tofile(fh)
    # asize (total number of cells = product of spatial dimensions)
    asize = np.prod(spatial_shape)
    np.array([asize], dtype=dti).tofile(fh)
    # data (flattened array)
    np.array(values, dtype=dtf).tofile(fh)
    
    fh.close()
    
    # Verify by reading back if requested
    if verify:
        try:
            gdata_read = pg.GData(output_path)
            read_values = gdata_read.get_values()
            read_grid = gdata_read.get_grid()
            
            # Check shapes match
            if read_values.shape != values.shape:
                raise ValueError(f"Verification failed: read shape {read_values.shape} != original {values.shape}")
            
            # Check values match (within numerical precision)
            if not np.allclose(read_values, values, rtol=1e-10, atol=1e-12):
                max_diff = np.max(np.abs(read_values - values))
                raise ValueError(f"Verification failed: values differ by up to {max_diff}")
            
            # Check grid dimensions match
            if len(read_grid) != len(grid):
                raise ValueError(f"Verification failed: read grid dims {len(read_grid)} != original {len(grid)}")
            
            for d, (read_g, orig_g) in enumerate(zip(read_grid, grid)):
                if len(read_g) != len(orig_g):
                    raise ValueError(f"Verification failed: grid dim {d} length {len(read_g)} != {len(orig_g)}")
            
            if verbose:
                print(f"  ✓ Verification passed: file read back correctly")
        except Exception as e:
            print(f"  ✗ Verification failed: {e}")
            raise
    
    print(f"  Converted: {filename} -> {os.path.basename(output_path)}")
    return output_path


if __name__ == '__main__':
      parser = argparse.ArgumentParser(description="Convert .npy files from DESC to .gkyl format.")
      parser.add_argument("input_dir", help="Path to directory containing .npy files from get_field_coords.py")
      parser.add_argument("output_dir", help="Path to output directory for .gkyl files")
      args = parser.parse_args()

      inp_dir = args.input_dir
      out_dir = args.output_dir

      nfp = 5 # set to 1 or 5 depending on what NFP in DESC code was set to
      lower = [0, 0, 0]
      upper = [1, 2*np.pi, (2*np.pi)/nfp]

      num_rho = 10   #5
      num_alpha = 50 # 20
      num_zeta = 7 # 7

      dims_corner = [num_rho, num_alpha, num_zeta]
      dims_S1 = [num_rho, 2*(num_alpha-1), 2*(num_zeta-1)]
      dims_S2 = [2*(num_rho-1), num_alpha, 2*(num_zeta-1)]
      dims_S3 = [2*(num_rho-1), 2*(num_alpha-1), num_zeta]
      dims_I =  [2*(num_rho-1), 2*(num_alpha-1), 2*(num_zeta-1)]
      N_corner = num_rho*num_alpha*num_zeta
      N_s1 = num_rho * 2*(num_alpha-1) *  2*(num_zeta-1)
      N_s2 = 2*(num_rho-1) * num_alpha * 2*(num_zeta-1)
      N_s3 = 2*(num_rho-1) * 2*(num_alpha-1) * num_zeta
      N_i = 2*(num_rho-1) * 2*(num_alpha-1) * 2*(num_zeta-1)

      convert_np_to_gkyl(inp_dir, "raz_corner.npy", N_corner, dims_corner, lower, upper, out_dir)
      convert_np_to_gkyl(inp_dir, "rzp_corner.npy", N_corner, dims_corner, lower, upper, out_dir)
      convert_np_to_gkyl(inp_dir, "Bmag_corner.npy", N_corner, dims_corner, lower, upper, out_dir)

      convert_np_to_gkyl(inp_dir, "tangents_I.npy", N_i, dims_I, lower, upper, out_dir)
      convert_np_to_gkyl(inp_dir, "Bmag_I.npy", N_i, dims_I, lower, upper, out_dir)
      convert_np_to_gkyl(inp_dir, "curl_B_hat_xyz_I.npy", N_i, dims_I, lower, upper, out_dir)
      convert_np_to_gkyl(inp_dir, "b_xyz_I.npy", N_i, dims_I, lower, upper, out_dir)

      convert_np_to_gkyl(inp_dir, "tangents_S1.npy", N_s1, dims_S1, lower, upper, out_dir)
      convert_np_to_gkyl(inp_dir, "Bmag_S1.npy", N_s1, dims_S1, lower, upper, out_dir)
      convert_np_to_gkyl(inp_dir, "curl_B_hat_xyz_S1.npy", N_s1, dims_S1, lower, upper, out_dir)

      convert_np_to_gkyl(inp_dir, "tangents_S2.npy", N_s2, dims_S2, lower, upper, out_dir)
      convert_np_to_gkyl(inp_dir, "Bmag_S2.npy", N_s2, dims_S2, lower, upper, out_dir)
      convert_np_to_gkyl(inp_dir, "curl_B_hat_xyz_S2.npy", N_s2, dims_S2, lower, upper, out_dir)

      convert_np_to_gkyl(inp_dir, "tangents_S3.npy", N_s3, dims_S3, lower, upper, out_dir)
      convert_np_to_gkyl(inp_dir, "Bmag_S3.npy", N_s3, dims_S3, lower, upper, out_dir)
      convert_np_to_gkyl(inp_dir, "curl_B_hat_xyz_S3.npy", N_s3, dims_S3, lower, upper, out_dir)

