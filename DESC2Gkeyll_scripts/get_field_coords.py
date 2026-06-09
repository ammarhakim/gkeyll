import os
import time
import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import jax
import jax.numpy as jnp

try:
    from desc.io import load
    from desc.grid import LinearGrid, Grid
    from desc.utils import rpz2xyz, rpz2xyz_vec, xyz2rpz
    from desc.compute.utils import _compute as compute_fun, get_transforms, get_profiles, get_params
    from desc.backend import root, jit, vmap
except ImportError:
    print("Failed to import DESC. Activate DESC environment (e.g. conda activate DESC)")
    exit(1)


def _rpz_to_rtz(eq, rpz, params=None, guess=None, tol=1e-6, maxiter=30):

    if params is None:
        params = eq.params_dict
    
    # Ensure rpz is a single point (3,)
    rpz = jnp.atleast_1d(rpz)
    if rpz.ndim != 1 or rpz.shape[0] != 3:
        raise ValueError(f"rpz must be shape (3,), got {rpz.shape}")
    
    # Residual function: R(rtz) - R_target, phi(rtz) - phi_target, Z(rtz) - Z_target
    # Works on a single point: rtz is (3,), rpz_target is (3,)
    @jit
    def residual(rtz, rpz_target):
        """Compute residual for root finding: rpz(rtz) - rpz_target."""
        rtz_2d = jnp.atleast_2d(rtz)  # (1, 3)
        # Provide dummy spacing for single point (needed for some computations like iota)
        spacing = jnp.array([[0.01, 0.01, 0.01]])  # (1, 3) - small dummy values
        grid = Grid(rtz_2d, spacing=spacing, sort=False, jitable=True)
        transforms = get_transforms(["R", "phi", "Z"], eq, grid, jitable=True)
        profiles = get_profiles(["R", "phi", "Z"], obj=eq, grid=grid, basis="rpz")
        data = compute_fun(
            eq,
            ["R", "phi", "Z"],
            params=params,
            transforms=transforms,
            profiles=profiles,
        )
        rpz_computed = jnp.array([data["R"][0], data["phi"][0], data["Z"][0]])
        return rpz_computed - rpz_target
    
    # Jacobian: derivatives of R, phi, Z w.r.t. rho, theta, zeta
    # Works on a single point: rtz is (3,), rpz_target is (3,) but not used
    @jit
    def jacobian(rtz, rpz_target):
        """Compute Jacobian matrix for root finding."""
        rtz_2d = jnp.atleast_2d(rtz)  # (1, 3)
        # Provide dummy spacing for single point (needed for some computations like iota)
        spacing = jnp.array([[0.01, 0.01, 0.01]])  # (1, 3) - small dummy values
        grid = Grid(rtz_2d, spacing=spacing, sort=False, jitable=True)
        transforms = get_transforms(["R_r", "R_t", "R_z", "phi_r", "phi_t", "phi_z", "Z_r", "Z_t", "Z_z"], eq, grid, jitable=True)
        profiles = get_profiles(["R_r", "R_t", "R_z", "phi_r", "phi_t", "phi_z", "Z_r", "Z_t", "Z_z"], obj=eq, grid=grid, basis="rpz")
        data = compute_fun(
            eq,
            ["R_r", "R_t", "R_z", "phi_r", "phi_t", "phi_z", "Z_r", "Z_t", "Z_z"],
            params=params,
            transforms=transforms,
            profiles=profiles,
        )
        # Reshape to (3, 3) Jacobian matrix
        J = jnp.array([
            [data["R_r"][0], data["R_t"][0], data["R_z"][0]],
            [data["phi_r"][0], data["phi_t"][0], data["phi_z"][0]],
            [data["Z_r"][0], data["Z_t"][0], data["Z_z"][0]],
        ])
        return J
    
    # Fixup function to keep rho in [0, 1] and handle periodicity
    @jit
    def fixup(rtz, *args):
        rho, theta, zeta = rtz
        rho = jnp.clip(rho, 1e-6, 1.0)  # Keep rho in valid range
        # Handle periodicity for theta and zeta
        theta = jnp.mod(theta, 2 * jnp.pi)
        zeta = jnp.mod(zeta, 2 * jnp.pi / eq.NFP)
        return jnp.array([rho, theta, zeta])
    
    # Initial guess: use rho from R (rough estimate), theta=0, zeta from phi
    if guess is None:
        R, phi, Z = rpz
        # Rough estimate: rho ~ sqrt((R - R0)^2 + Z^2) / a
        # R_lmn[0] is the (0,0,0) coefficient which is the major radius at axis
        R0 = eq.R_lmn[0] if len(eq.R_lmn) > 0 else 1.0
        rho_guess = jnp.clip(jnp.abs(R - R0) / 2.0, 0.1, 0.9)
        theta_guess = 0.0
        zeta_guess = phi  # Assuming phi = zeta for now
        guess = jnp.array([rho_guess, theta_guess, zeta_guess])
    
    # Ensure guess is a single point (3,)
    guess = jnp.atleast_1d(guess)
    if guess.ndim != 1 or guess.shape[0] != 3:
        raise ValueError(f"guess must be shape (3,), got {guess.shape}")
    
    # Root finding for single point
    rtz = root(
        residual,
        guess,
        jac=jacobian,
        args=(rpz,),
        fixup=fixup,
        tol=tol,
        maxiter=maxiter,
    )
    
    return rtz


def _compute_Bhat_xyz_from_xyz(eq, xyz, params=None, guess=None):

    if params is None:
        params = eq.params_dict
    
    # Ensure xyz is a single point (3,)
    xyz = jnp.atleast_1d(xyz)
    if xyz.ndim != 1 or xyz.shape[0] != 3:
        raise ValueError(f"xyz must be shape (3,), got {xyz.shape}")
    
    # Step 1: Convert xyz -> rpz
    rpz = xyz2rpz(xyz)  # Returns (3,)
    
    # Step 2: Convert rpz -> rtz (using root finding)
    rtz = _rpz_to_rtz(eq, rpz, params=params, guess=guess)  # Returns (3,)
    
    # Step 3: Compute B_R, B_phi, B_Z, |B| at rtz coordinates
    # Grid needs (N, 3) shape, so add dimension
    rtz_2d = jnp.atleast_2d(rtz)  # (1, 3)
    # Provide dummy spacing for single point (needed for some computations like iota)
    spacing = jnp.array([[0.01, 0.01, 0.01]])  # (1, 3) - small dummy values
    grid = Grid(rtz_2d, spacing=spacing, sort=False, jitable=True)
    transforms = get_transforms(["B_R", "B_phi", "B_Z", "|B|", "phi"], eq, grid, jitable=True)
    profiles = get_profiles(["B_R", "B_phi", "B_Z", "|B|"], obj=eq, grid=grid, basis="rpz")
    
    data = compute_fun(
        eq,
        ["B_R", "B_phi", "B_Z", "|B|", "phi"],
        params=params,
        transforms=transforms,
        profiles=profiles,
    )
    
    # Extract scalar values (first element since grid has 1 point)
    B_R = data["B_R"][0]
    B_phi = data["B_phi"][0]
    B_Z = data["B_Z"][0]
    Bmag = data["|B|"][0]
    phi = data["phi"][0]
    
    # Step 4: Compute Bhat_rpz = [B_R/|B|, B_phi/|B|, B_Z/|B|]
    eps = 1e-12
    Bhat_rpz = jnp.array([B_R, B_phi, B_Z]) / (Bmag + eps)  # Shape (3,)
    
    # Step 5: Convert Bhat_rpz -> Bhat_xyz
    # rpz2xyz_vec needs (N, 3) shape, so add dimension then remove it
    Bhat_rpz_2d = jnp.atleast_2d(Bhat_rpz)  # (1, 3)
    phi_1d = jnp.atleast_1d(phi)  # (1,)
    Bhat_xyz_2d = rpz2xyz_vec(Bhat_rpz_2d, phi=phi_1d)  # (1, 3)
    Bhat_xyz = Bhat_xyz_2d[0]  # (3,)
    
    return Bhat_xyz


def get_Bhat_x(xyz, eq, params=None, guess=None):
    Bhat_xyz = _compute_Bhat_xyz_from_xyz(eq, xyz, params=params, guess=guess)
    return Bhat_xyz[0]


def get_Bhat_y(xyz, eq, params=None, guess=None):
    Bhat_xyz = _compute_Bhat_xyz_from_xyz(eq, xyz, params=params, guess=guess)
    return Bhat_xyz[1]


def get_Bhat_z(xyz, eq, params=None, guess=None):
    Bhat_xyz = _compute_Bhat_xyz_from_xyz(eq, xyz, params=params, guess=guess)
    return Bhat_xyz[2]


def grad_Bhat_x(xyz, eq, params=None, guess=None):
    return jax.grad(get_Bhat_x, argnums=0)(xyz, eq, params=params, guess=guess)


def grad_Bhat_y(xyz, eq, params=None, guess=None):
    return jax.grad(get_Bhat_y, argnums=0)(xyz, eq, params=params, guess=guess)


def grad_Bhat_z(xyz, eq, params=None, guess=None):
    return jax.grad(get_Bhat_z, argnums=0)(xyz, eq, params=params, guess=guess)


def curl_Bhat(xyz, eq, params=None, guess=None, chunk_size=None):

    xyz = jnp.atleast_2d(xyz)  # Ensure (N, 3) shape
    N = xyz.shape[0]
    
    # If chunk_size is specified and array is large, process in chunks
    if chunk_size is not None and N > chunk_size:
        num_chunks = (N + chunk_size - 1) // chunk_size
        curl_chunks = []
        
        for i in range(num_chunks):
            start_idx = i * chunk_size
            end_idx = min((i + 1) * chunk_size, N)
            xyz_chunk = xyz[start_idx:end_idx]
            
            # Compute gradients for this chunk
            grad_x = jax.vmap(lambda x: grad_Bhat_x(x, eq, params=params, guess=guess))(xyz_chunk)
            grad_y = jax.vmap(lambda x: grad_Bhat_y(x, eq, params=params, guess=guess))(xyz_chunk)
            grad_z = jax.vmap(lambda x: grad_Bhat_z(x, eq, params=params, guess=guess))(xyz_chunk)
            
            # Compute curl for this chunk
            curl_x = grad_z[:, 1] - grad_y[:, 2]
            curl_y = grad_x[:, 2] - grad_z[:, 0]
            curl_z = grad_y[:, 0] - grad_x[:, 1]
            curl_chunk = jnp.stack([curl_x, curl_y, curl_z], axis=1)
            curl_chunks.append(curl_chunk)
        
        curl_Bhat_xyz = jnp.concatenate(curl_chunks, axis=0)
    else:
        # Process all points at once (memory intensive for large N)
        grad_x = jax.vmap(lambda x: grad_Bhat_x(x, eq, params=params, guess=guess))(xyz)
        grad_y = jax.vmap(lambda x: grad_Bhat_y(x, eq, params=params, guess=guess))(xyz)
        grad_z = jax.vmap(lambda x: grad_Bhat_z(x, eq, params=params, guess=guess))(xyz)
        
        curl_x = grad_z[:, 1] - grad_y[:, 2]
        curl_y = grad_x[:, 2] - grad_z[:, 0]
        curl_z = grad_y[:, 0] - grad_x[:, 1]
        curl_Bhat_xyz = jnp.stack([curl_x, curl_y, curl_z], axis=1)
    
    # Return with same shape as input
    if xyz.shape[0] == 1:
        return curl_Bhat_xyz[0]  # (3,)
    return curl_Bhat_xyz  # (N, 3)


def get_quad_points(raz):
    raz = np.asarray(raz)
    if raz.ndim != 2 or raz.shape[1] != 3:
        raise ValueError(f"raz must be Nx3 array, got shape {raz.shape}")
    
    # Initialize output arrays
    raz_S1_list = []
    raz_S2_list = []
    raz_S3_list = []
    raz_I_list = []
    
    # Find min and max values for each coordinate
    coord_mins = raz.min(axis=0)
    coord_maxs = raz.max(axis=0)
    r_min, a_min, z_min = coord_mins
    r_max, a_max, z_max = coord_maxs
    
    # Calculate delta values (spacing between linearly spaced points)
    # For linearly spaced grids: delta = (max - min) / (n_unique_points - 1)
    n_unique = np.array([len(np.unique(raz[:, i])) for i in range(3)])
    deltas = np.where(n_unique > 1, (coord_maxs - coord_mins) / (n_unique - 1), 0.0)
    delta_r, delta_a, delta_z = deltas
    
    # Tolerance for comparing floats
    tol = 1e-10
    # Calculate e values: e = 1/(2*sqrt(3))
    e_base = 1.0 / (2.0 * np.sqrt(3.0))
    e_r = delta_r * e_base
    e_a = delta_a * e_base
    e_z = delta_z * e_base
    
    # Process each point in raz
    for point in raz:
        r, a, z = point
        
        # Check if point is at boundaries
        at_r_min = abs(r - r_min) < tol
        at_r_max = abs(r - r_max) < tol
        at_a_min = abs(a - a_min) < tol
        at_a_max = abs(a - a_max) < tol
        at_z_min = abs(z - z_min) < tol
        at_z_max = abs(z - z_max) < tol
        
        # Create S3 points (vary rho and alpha, keep zeta same)
        if at_r_min and at_a_min:
            # Only r+e_r, a+e_a
            raz_S3_list.append([r + e_r, a + e_a, z])
        elif at_r_min and at_a_max:
            # Only r+e_r, a-e_a
            raz_S3_list.append([r + e_r, a - e_a, z])
        elif at_r_max and at_a_min:
            # Only r-e_r, a+e_a
            raz_S3_list.append([r - e_r, a + e_a, z])
        elif at_r_max and at_a_max:
            # Only r-e_r, a-e_a
            raz_S3_list.append([r - e_r, a - e_a, z])
        elif at_r_min:
            # Only r+e_r points (both a+e_a and a-e_a)
            raz_S3_list.append([r + e_r, a + e_a, z])
            raz_S3_list.append([r + e_r, a - e_a, z])
        elif at_r_max:
            # Only r-e_r points (both a+e_a and a-e_a)
            raz_S3_list.append([r - e_r, a + e_a, z])
            raz_S3_list.append([r - e_r, a - e_a, z])
        elif at_a_min:
            # Only a+e_a points (both r+e_r and r-e_r)
            raz_S3_list.append([r + e_r, a + e_a, z])
            raz_S3_list.append([r - e_r, a + e_a, z])
        elif at_a_max:
            # Only a-e_a points (both r+e_r and r-e_r)
            raz_S3_list.append([r + e_r, a - e_a, z])
            raz_S3_list.append([r - e_r, a - e_a, z])
        else:
            # Interior: all 4 combinations
            raz_S3_list.append([r + e_r, a + e_a, z])
            raz_S3_list.append([r + e_r, a - e_a, z])
            raz_S3_list.append([r - e_r, a + e_a, z])
            raz_S3_list.append([r - e_r, a - e_a, z])
        
        # Create S1 points (vary alpha and zeta, keep rho same)
        if at_z_min and at_a_min:
            # Only z+e_z, a+e_a
            raz_S1_list.append([r, a + e_a, z + e_z])
        elif at_z_min and at_a_max:
            # Only z+e_z, a-e_a
            raz_S1_list.append([r, a - e_a, z + e_z])
        elif at_z_max and at_a_min:
            # Only z-e_z, a+e_a
            raz_S1_list.append([r, a + e_a, z - e_z])
        elif at_z_max and at_a_max:
            # Only z-e_z, a-e_a
            raz_S1_list.append([r, a - e_a, z - e_z])
        elif at_z_min:
            # Only z+e_z points (both a+e_a and a-e_a)
            raz_S1_list.append([r, a + e_a, z + e_z])
            raz_S1_list.append([r, a - e_a, z + e_z])
        elif at_z_max:
            # Only z-e_z points (both a+e_a and a-e_a)
            raz_S1_list.append([r, a + e_a, z - e_z])
            raz_S1_list.append([r, a - e_a, z - e_z])
        elif at_a_min:
            # Only a+e_a points (both z+e_z and z-e_z)
            raz_S1_list.append([r, a + e_a, z + e_z])
            raz_S1_list.append([r, a + e_a, z - e_z])
        elif at_a_max:
            # Only a-e_a points (both z+e_z and z-e_z)
            raz_S1_list.append([r, a - e_a, z + e_z])
            raz_S1_list.append([r, a - e_a, z - e_z])
        else:
            # Interior: all 4 combinations
            raz_S1_list.append([r, a + e_a, z + e_z])
            raz_S1_list.append([r, a - e_a, z + e_z])
            raz_S1_list.append([r, a + e_a, z - e_z])
            raz_S1_list.append([r, a - e_a, z - e_z])
        
        # Create S2 points (vary rho and zeta, keep alpha same)
        if at_r_min and at_z_min:
            # Only r+e_r, z+e_z
            raz_S2_list.append([r + e_r, a, z + e_z])
        elif at_r_min and at_z_max:
            # Only r+e_r, z-e_z
            raz_S2_list.append([r + e_r, a, z - e_z])
        elif at_r_max and at_z_min:
            # Only r-e_r, z+e_z
            raz_S2_list.append([r - e_r, a, z + e_z])
        elif at_r_max and at_z_max:
            # Only r-e_r, z-e_z
            raz_S2_list.append([r - e_r, a, z - e_z])
        elif at_r_min:
            # Only r+e_r points (both z+e_z and z-e_z)
            raz_S2_list.append([r + e_r, a, z + e_z])
            raz_S2_list.append([r + e_r, a, z - e_z])
        elif at_r_max:
            # Only r-e_r points (both z+e_z and z-e_z)
            raz_S2_list.append([r - e_r, a, z + e_z])
            raz_S2_list.append([r - e_r, a, z - e_z])
        elif at_z_min:
            # Only z+e_z points (both r+e_r and r-e_r)
            raz_S2_list.append([r + e_r, a, z + e_z])
            raz_S2_list.append([r - e_r, a, z + e_z])
        elif at_z_max:
            # Only z-e_z points (both r+e_r and r-e_r)
            raz_S2_list.append([r + e_r, a, z - e_z])
            raz_S2_list.append([r - e_r, a, z - e_z])
        else:
            # Interior: all 4 combinations
            raz_S2_list.append([r + e_r, a, z + e_z])
            raz_S2_list.append([r + e_r, a, z - e_z])
            raz_S2_list.append([r - e_r, a, z + e_z])
            raz_S2_list.append([r - e_r, a, z - e_z])
    
    # Convert lists to arrays
    raz_S1 = np.array(raz_S1_list) if raz_S1_list else np.empty((0, 3))
    raz_S2 = np.array(raz_S2_list) if raz_S2_list else np.empty((0, 3))
    raz_S3 = np.array(raz_S3_list) if raz_S3_list else np.empty((0, 3))
    
    # Create interior points from S3 points
    # For each point in S2, create 2 points by moving z +/- e_z
    # At z_min we only do z+e_z, at z_max we only do z-e_z
    for s3_point in raz_S3:
        r_s3, a_s3, z_s3 = s3_point
        at_z_min_s3 = abs(z_s3 - z_min) < tol
        at_z_max_s3 = abs(z_s3 - z_max) < tol
        
        if at_z_min_s3:
            # Only z+e_z
            raz_I_list.append([r_s3, a_s3, z_s3 + e_z])
        elif at_z_max_s3:
            # Only z-e_z
            raz_I_list.append([r_s3, a_s3, z_s3 - e_z])
        else:
            # Both z+e_z and z-e_z
            raz_I_list.append([r_s3, a_s3, z_s3 + e_z])
            raz_I_list.append([r_s3, a_s3, z_s3 - e_z])
    
    raz_I = np.array(raz_I_list) if raz_I_list else np.empty((0, 3))
    
    return raz_S1, raz_S2, raz_S3, raz_I


def compute_gkeyll_nodes(raz, eq, node_type = "corner"):

    # Convert raz to rtz coordinates
    rtz = eq.map_coordinates(raz, inbasis=("rho", "alpha", "zeta"), outbasis=("rho", "theta", "zeta"))
    grid_rtz = Grid(rtz, coordinates="rtz", sort=False, jitable=False)    # Compute R, phi, Z and all first-order partial derivatives    
    
    names = ["R", "phi", "Z", "X_r", "X_t", "X_z", "Y_r", "Y_t", "Y_z", "Z_r", "Z_t", "Z_z", "alpha_t","|B|", "B", "B_R", "B_phi", "B_Z", "e_rho", "e_alpha", "e_zeta"]
    data = eq.compute(names, grid=grid_rtz)

    # Assemble outputs in rzp order (R, Z, phi) instead of rpz (R, phi, Z)
    rzp = np.column_stack([data["R"], data["Z"], data["phi"]])
    # Convert rzp node locations to xyz coordinates (need to convert back to rpz for rpz2xyz)
    xyz = jnp.asarray(rpz2xyz(np.column_stack([data["R"], data["phi"], data["Z"]])))
    
    Bmag = np.asarray(data["|B|"])   
    # Convert B field from rpz to xyz coordinates
    # rpz2xyz_vec expects rpz format [B_R, B_phi, B_Z] and phi coordinate
    B_R = np.asarray(data["B_R"])
    B_phi = np.asarray(data["B_phi"])
    B_Z = np.asarray(data["B_Z"])
    B_rpz = np.column_stack([B_R, B_phi, B_Z])  # rpz order: [B_R, B_phi, B_Z]
    phi_array = np.asarray(data["phi"])
    B_xyz = np.asarray(rpz2xyz_vec(B_rpz, phi=phi_array))
    b_xyz = B_xyz / Bmag[:, np.newaxis]
    
    if node_type != "corner":
        curl_B_hat_xyz = np.asarray(curl_Bhat(xyz, eq, chunk_size=5000)) 
    else: 
        curl_B_hat_xyz = None   


    dtheta_dalpha = 1.0 / (data["alpha_t"])   
    # Tangents in order: [e_1, e_2, e_3]
    # e_1: dXdrho, dYdrho, dZdrho
    # e_2: dXdalpha, dYdalpha, dZdalpha
    # e_3: dXdzeta, dYdzeta, dZdzeta
    tangents = np.column_stack([
        data["X_r"],   data["Y_r"],   data["Z_r"],
        data["X_t"] * dtheta_dalpha,   data["Y_t"] * dtheta_dalpha,   data["Z_t"] * dtheta_dalpha,
        data["X_z"], data["Y_z"], data["Z_z"]
    ])
    
    
    if node_type == "corner":
        return raz, rzp, Bmag
    elif node_type == "surface":
        return tangents, Bmag, curl_B_hat_xyz
    elif node_type == "interior":
        return tangents, Bmag, curl_B_hat_xyz, b_xyz

def main():
    parser = argparse.ArgumentParser(description="Compute field-aligned coordinates from a DESC equilibrium.")
    parser.add_argument("eq_file", help="Path to a DESC equilibrium file (H5)")
    parser.add_argument("out_dir", help="Path to output directory for .npy files")
    args = parser.parse_args()

    start_time = time.time()
    # Load equilibrium from file
    out_dir = args.out_dir
    eq = load(args.eq_file)
    # Desired grid sizes
    num_rho = 10
    num_alpha = 50
    num_zeta = 7
    # set periods = eq.NFP to plot only 1/5th of the plasma (5 field periodicity)
    periods = eq.NFP
    #Create linspace arrays for rho aloha zeta and then use them to make a grid of raz points
    rho = np.linspace(0.0, 1.0, num_rho)
    alpha = np.linspace(0.0, 2 * np.pi, num_alpha, endpoint=True) 
    zeta = np.linspace(0.0, 2 * np.pi / periods, num_zeta, endpoint=True) 
    grid_raz = LinearGrid(rho=rho, theta=alpha, zeta=zeta, NFP=periods)
    raz_corner = grid_raz.nodes
    
    #using the linearly spaced raz nodes as corner nodes, generate the quadrature nodes
    #For quad arrays num_rho, num_alpha, num_zetas
    #For S1: nrho = num_rho; nalpha = 2*(num_alpha-1); nzeta = 2*(num_zeta-1)
    #For S2: nrho = 2*(num_rho-1); nalpha = num_alpha; nzeta = 2*(num_zeta-1)
    #For S3: nrho = 2*(num_rho-1); nalpha = 2*(num_alpha-1); nzeta = num_zeta
  
    raz_S1, raz_S2, raz_S3, raz_I = get_quad_points(raz_corner)
    # To sort the points such that zeta varies fastest, then theta, then rho,
    # we need to lexsort by (rho, theta, zeta); but since lexsort sorts by last key first,
    # provide keys in order [zeta, theta, rho].
    def sort_raz_fortran_order(raz_arr):
        # raz_arr is shape (N, 3): (rho, theta, zeta)
        idx = np.lexsort((raz_arr[:, 2], raz_arr[:, 1], raz_arr[:, 0]))  # sort by (rho, theta, zeta)
        return raz_arr[idx]
    
    raz_corner = sort_raz_fortran_order(raz_corner)
    raz_S1 = sort_raz_fortran_order(raz_S1)
    raz_S2 = sort_raz_fortran_order(raz_S2)
    raz_S3 = sort_raz_fortran_order(raz_S3)
    raz_I = sort_raz_fortran_order(raz_I)
    print("Starting to compute nodal values")

    t0 = time.time()
    # Compute desired values for corner nodes
    raz_corner, rzp_corner, Bmag_corner = compute_gkeyll_nodes(raz_corner, eq, "corner")
    t1 = time.time()
    print(f"Finished computing corner nodes in {t1 - t0:.2f} seconds")
    # Compute desired values for surface quadrature nodes    
    t2 = time.time()
    tangents_S1, Bmag_S1, curl_B_hat_xyz_S1 = compute_gkeyll_nodes(raz_S1, eq, "surface")
    t3 = time.time()
    print(f"Finished computing S1 nodes in {t3 - t2:.2f} seconds")
    t4 = time.time()
    tangents_S2, Bmag_S2, curl_B_hat_xyz_S2 = compute_gkeyll_nodes(raz_S2, eq, "surface")
    t5 = time.time()
    print(f"Finished computing S2 nodes in {t5 - t4:.2f} seconds")
    t6 = time.time()
    tangents_S3, Bmag_S3, curl_B_hat_xyz_S3 = compute_gkeyll_nodes(raz_S3, eq, "surface")
    t7 = time.time()
    print(f"Finished computing S3 nodes in {t7 - t6:.2f} seconds")
    # Compute desired values for interior quadrature nodes
    t8 = time.time()
    tangents_I, Bmag_I, curl_B_hat_xyz_I, b_xyz_I = compute_gkeyll_nodes(raz_I, eq, "interior")
    t9 = time.time()
    print(f"Finished computing interior nodes in {t9 - t8:.2f} seconds")
    
    # Output directory
    out_dir = os.path.join(os.getcwd(), out_dir)
    os.makedirs(out_dir, exist_ok=True)

    # Organize all arrays: {array_name: {suffix: array}}
    # Only created above: rza_corner, rzp_corner, Bmag_corner, drzp_draz_S1, Bmag_S1, curl_B_hat_xyz_S1, drzp_draz_S2, Bmag_S2, curl_B_hat_xyz_S2, drzp_draz_S3, Bmag_S3, curl_B_hat_xyz_S3, drzp_draz_I, Bmag_I, curl_B_hat_xyz_I, B_xyz_I
    # Write raz_corner as a .txt file (keep original 2D format for text)
    np.savetxt(os.path.join(out_dir, "raz_corner.txt"), raz_corner, fmt="%.18e")
    np.savetxt(os.path.join(out_dir, "raz_S1.txt"), raz_S1, fmt="%.18e")
    np.savetxt(os.path.join(out_dir, "raz_S2.txt"), raz_S2, fmt="%.18e")
    np.savetxt(os.path.join(out_dir, "raz_S3.txt"), raz_S3, fmt="%.18e")
    np.savetxt(os.path.join(out_dir, "raz_I.txt"), raz_I, fmt="%.18e")
    # Save the computed arrays as .npy files
    np.save(os.path.join(out_dir, "raz_corner.npy"), raz_corner)
    np.save(os.path.join(out_dir, "rzp_corner.npy"), rzp_corner)
    np.save(os.path.join(out_dir, "Bmag_corner.npy"), Bmag_corner)
    
    np.save(os.path.join(out_dir, "tangents_S1.npy"), tangents_S1)
    np.save(os.path.join(out_dir, "Bmag_S1.npy"), Bmag_S1)
    np.save(os.path.join(out_dir, "curl_B_hat_xyz_S1.npy"), curl_B_hat_xyz_S1)
    
    np.save(os.path.join(out_dir, "tangents_S2.npy"), tangents_S2)
    np.save(os.path.join(out_dir, "Bmag_S2.npy"), Bmag_S2)
    np.save(os.path.join(out_dir, "curl_B_hat_xyz_S2.npy"), curl_B_hat_xyz_S2)
    
    np.save(os.path.join(out_dir, "tangents_S3.npy"), tangents_S3)
    np.save(os.path.join(out_dir, "Bmag_S3.npy"), Bmag_S3)
    np.save(os.path.join(out_dir, "curl_B_hat_xyz_S3.npy"), curl_B_hat_xyz_S3)

    np.save(os.path.join(out_dir, "tangents_I.npy"), tangents_I)
    np.save(os.path.join(out_dir, "Bmag_I.npy"), Bmag_I)
    np.save(os.path.join(out_dir, "curl_B_hat_xyz_I.npy"), curl_B_hat_xyz_I)
    np.save(os.path.join(out_dir, "b_xyz_I.npy"), b_xyz_I)


    elapsed_time = time.time() - start_time
    print(f"\nTotal execution time: {elapsed_time:.2f} seconds")


if __name__ == "__main__":
    main()


