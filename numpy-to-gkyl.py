#!/usr/bin/env python3
"""
Convert all .npy (NumPy array) files in a folder to .gkyl files using postgkyl.

Usage:
    python numpy-to-gkyl.py <folder_path>
"""

import os
import sys
import glob
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


def convert_np_to_gkyl(np_file_path, output_dir=None, lower=None, upper=None, verbose=False, verify=False):
    """Convert a single .np file to .gkyl format.
    
    Args:
        np_file_path: Path to input .np file
        output_dir: Directory to write output .gkyl file (default: same directory as input)
        lower: Lower bounds for grid (default: zeros)
        upper: Upper bounds for grid (default: ones)
        verbose: Print detailed shape information
        verify: Read back the written file to verify correctness
    
    Returns:
        Path to output file
    """
    # Load numpy array
    arr = np.load(np_file_path)
    original_shape = arr.shape
    original_dtype = arr.dtype
    
    # Always print numpy array shape
    print(f"\n  File: {os.path.basename(np_file_path)}")
    print(f"  NumPy array shape: {original_shape}, dtype: {original_dtype}")
    
    if verbose:
        print(f"  Full path: {np_file_path}")
    
    # Handle different array formats
    if arr.dtype.names is not None:
        # Structured array - take first field
        arr = arr[arr.dtype.names[0]]
        if verbose:
            print(f"  Structured array, using field: {arr.dtype.names[0]}")
    
    # Determine shape: last dimension is components, rest are spatial
    # If array is 1D, treat entire array as single component
    if arr.ndim == 1:
        # 1D array: treat as single component
        values = arr.reshape(-1, 1)
        spatial_shape = (arr.shape[0],)
        num_comps = 1
    else:
        # Multi-dimensional: last dimension is components
        spatial_shape = arr.shape[:-1]
        values = arr
        num_comps = arr.shape[-1]
    
    if verbose:
        print(f"  Interpreted as: spatial_shape={spatial_shape}, num_components={num_comps}")
        print(f"  Values shape: {values.shape}")
    
    # Validate that we have at least one spatial dimension
    if len(spatial_shape) == 0:
        raise ValueError(f"Cannot convert array with shape {arr.shape}: no spatial dimensions")
    
    # Validate bounds if provided
    if lower is not None and len(lower) != len(spatial_shape):
        raise ValueError(f"Lower bounds length ({len(lower)}) doesn't match spatial dimensions ({len(spatial_shape)})")
    if upper is not None and len(upper) != len(spatial_shape):
        raise ValueError(f"Upper bounds length ({len(upper)}) doesn't match spatial dimensions ({len(spatial_shape)})")
    
    # Create grid from spatial shape
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
        file_name = os.path.basename(np_file_path)
        base_name = os.path.splitext(file_name)[0]
        output_path = os.path.join(output_dir, base_name + '.gkyl')
    
    # Write to .gkyl file
    gdata.write(output_path, extension='gkyl')
    
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
    
    print(f"  Converted: {os.path.basename(np_file_path)} -> {os.path.basename(output_path)}")
    return output_path


def convert_folder(folder_path, lower=None, upper=None, verbose=False, verify=False):
    """Convert all .npy files in a folder to .gkyl files.
    
    Args:
        folder_path: Path to folder containing .npy files
        lower: Lower bounds for grid (default: zeros)
        upper: Upper bounds for grid (default: ones)
        verbose: Print detailed shape information for each file
        verify: Read back each written file to verify correctness
    """
    if not os.path.isdir(folder_path):
        print(f"Error: {folder_path} is not a valid directory")
        return
    
    # Create output folder with _gkyl suffix in current working directory
    folder_abs_path = os.path.abspath(folder_path)
    folder_name = os.path.basename(folder_abs_path)
    # Use current working directory instead of parent directory
    current_dir = os.getcwd()
    output_folder = os.path.join(current_dir, folder_name + '_gkyl')
    
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    print(f"Output folder: {output_folder}")
    
    # Find all .npy files
    pattern = os.path.join(folder_path, '*.npy')
    np_files = glob.glob(pattern)
    
    if not np_files:
        print(f"No .npy files found in {folder_path}")
        return
    
    print(f"Found {len(np_files)} .npy file(s) in {folder_path}")
    if verify:
        print("Verification enabled: will read back each file to check correctness")
    
    # Convert each file
    success_count = 0
    error_count = 0
    
    for np_file in np_files:
        try:
            convert_np_to_gkyl(np_file, output_dir=output_folder, 
                              lower=lower, upper=upper, verbose=verbose, verify=verify)
            success_count += 1
        except Exception as e:
            error_count += 1
            print(f"  ✗ Error converting {os.path.basename(np_file)}: {e}")
            if verbose:
                import traceback
                traceback.print_exc()
    
    print(f"\nConversion complete!")
    print(f"  Successfully converted: {success_count}/{len(np_files)} files")
    if error_count > 0:
        print(f"  Errors: {error_count} files")
    print(f"  Output written to: {output_folder}")


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Convert all .npy (NumPy array) files in a folder to .gkyl files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python numpy-to-gkyl.py ./data
  python numpy-to-gkyl.py ./data --lower 0,0,0 --upper 1,1,1
  python numpy-to-gkyl.py ./data --verbose --verify
        """
    )
    parser.add_argument('folder_path', help='Path to folder containing .npy files')
    parser.add_argument('--lower', type=str, help='Lower bounds for grid (comma-separated, e.g., "0,0,0")')
    parser.add_argument('--upper', type=str, help='Upper bounds for grid (comma-separated, e.g., "1,1,1")')
    parser.add_argument('--verbose', '-v', action='store_true', 
                       help='Print detailed shape information for each file')
    parser.add_argument('--verify', action='store_true',
                       help='Read back each written file to verify correctness')
    
    args = parser.parse_args()
    
    # Parse bounds
    lower = None
    upper = None
    
    if args.lower:
        lower = np.array([float(x) for x in args.lower.split(',')])
    
    if args.upper:
        upper = np.array([float(x) for x in args.upper.split(',')])
    
    convert_folder(args.folder_path, lower=lower, upper=upper, 
                  verbose=args.verbose, verify=args.verify)

