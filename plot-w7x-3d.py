import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import postgkyl as pg

# Read the gkyl files
data = pg.GData("./W7-X_geometry/W7X-nodes.gkyl")  # (R, Z, phi)
vals = data.get_values()

# Companion file with the *same nodes* in (rho, alpha, zeta)
nu_data = pg.GData("./W7-X_geometry/W7X-nodes-computational.gkyl")     # (rho, alpha, zeta)
nu_vals = nu_data.get_values()



# Extract coordinates - assuming the data has R, Z, phi or x, y, z coordinates
# The shape depends on the data structure
if vals.ndim == 4:
    # If 4D array (nx, ny, nz, 3), extract coordinates from the last dimension and flatten spatial dimensions
    if vals.shape[3] >= 3:
        R = vals[:, :, :, 0].flatten()
        Z = vals[:, :, :, 1].flatten()
        phi = vals[:, :, :, 2].flatten()
    else:
        # If fewer components, try to flatten
        R = vals.flatten()
        Z = np.zeros_like(R)
        phi = np.zeros_like(R)
elif vals.ndim == 3:
    # If 3D array, extract coordinates from the last dimension and flatten spatial dimensions
    if vals.shape[2] >= 3:
        R = vals[:, :, 0].flatten()
        Z = vals[:, :, 1].flatten()
        phi = vals[:, :, 2].flatten()
    else:
        # If fewer dimensions, try to flatten
        R = vals.flatten()
        Z = np.zeros_like(R)
        phi = np.zeros_like(R)
elif vals.ndim == 2:
    # If 2D array, assume first column is R, second is Z, third is phi
    if vals.shape[1] >= 3:
        R = vals[:, 0]
        Z = vals[:, 1]
        phi = vals[:, 2]
    else:
        R = vals[:, 0]
        Z = vals[:, 1] if vals.shape[1] > 1 else np.zeros_like(R)
        phi = np.zeros_like(R)
else:
    # Flatten and try to reshape
    R = vals.flatten()
    Z = np.zeros_like(R)
    phi = np.zeros_like(R)

# Convert to Cartesian coordinates if needed (R, Z, phi -> x, y, z)
# If the data is already in Cartesian, use directly
# For cylindrical coordinates: x = R*cos(phi), y = R*sin(phi), z = Z
if np.any(phi != 0) or np.max(np.abs(phi)) > 0.1:
    # Assume cylindrical coordinates
    x = R * np.cos(phi)
    y = R * np.sin(phi)
    z = Z
else:
    # Assume already Cartesian or use R, Z, phi directly
    x = R
    y = Z
    z = phi if np.any(phi != 0) else np.zeros_like(x)

# Create 3D plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot the points from first gkyl file (larger markers for visibility)
ax.scatter(x, y, z, marker='.', s=6, alpha=0.7, label='W7-X nodes', color='blue')

# Set labels
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_zlabel('Z [m]')
ax.set_title('W7-X nodes 3D Plot')

# Set equal aspect ratio for better visualization
# Ensure all arrays are 1D before concatenating
all_x = np.atleast_1d(x).flatten()
all_y = np.atleast_1d(y).flatten()
all_z = np.atleast_1d(z).flatten()

max_range = np.array([all_x.max()-all_x.min(), all_y.max()-all_y.min(), all_z.max()-all_z.min()]).max() / 2.0
mid_x = (all_x.max()+all_x.min()) * 0.5
mid_y = (all_y.max()+all_y.min()) * 0.5
mid_z = (all_z.max()+all_z.min()) * 0.5
ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

# Add legend
ax.legend()

plt.tight_layout()

# ---------------------------------------------------------------------------
# 2D slices at multiple PHI / ZETA values (R–Z cross-sections) USING NODAL
# ---------------------------------------------------------------------------

# Tolerance for selecting the zeta slice (in radians)
phi_tol = 1e-2

if vals.ndim == 4 and nu_vals.ndim == 4 and vals.shape[:3] == nu_vals.shape[:3] and nu_vals.shape[3] >= 3:
    # Geometry in R-Z-phi
    R3 = vals[:, :, :, 0]
    Z3 = vals[:, :, :, 1]
    phi3 = vals[:, :, :, 2]

    # Logical coordinates in (rho, alpha, zeta)
    rho3 = nu_vals[:, :, :, 0]
    alpha3 = nu_vals[:, :, :, 1]
    zeta3 = nu_vals[:, :, :, 2]

    nx, ny, nz = zeta3.shape

    # ---------------------------------------------------------------------------
    # Rectangular grid plot of nu_data (rho, alpha, zeta)
    # ---------------------------------------------------------------------------
    
    # Create 3D plot for rectangular grid
    fig_nu = plt.figure(figsize=(10, 8))
    ax_nu = fig_nu.add_subplot(111, projection='3d')
    
    # Flatten the coordinates for scatter plot
    rho_flat = rho3.flatten()
    alpha_flat = alpha3.flatten()
    zeta_flat = zeta3.flatten()
    
    # Plot the points in rectangular (rho, alpha, zeta) space
    # Swap alpha and zeta: alpha goes upwards (z-axis), zeta goes depthwards (y-axis)
    ax_nu.scatter(rho_flat, zeta_flat, alpha_flat, marker='.', s=6, alpha=0.7, 
                  label='raz nodes', color='red')
    
    # Set labels (rho=x, zeta=y, alpha=z)
    ax_nu.set_xlabel('rho')
    ax_nu.set_ylabel('zeta [rad]')
    ax_nu.set_zlabel('alpha [rad]')
    ax_nu.set_title('W7-X nodes in rho alpha zeta space')
    
    # Set equal aspect ratio for better visualization
    all_rho = np.atleast_1d(rho_flat).flatten()
    all_alpha = np.atleast_1d(alpha_flat).flatten()
    all_zeta = np.atleast_1d(zeta_flat).flatten()
    
    max_range = np.array([all_rho.max()-all_rho.min(), 
                          all_alpha.max()-all_alpha.min(), 
                          all_zeta.max()-all_zeta.min()]).max() / 2.0
    mid_alpha = (all_alpha.max()+all_alpha.min()) * 0.5
    mid_zeta = (all_zeta.max()+all_zeta.min()) * 0.5
    ax_nu.set_xlim(0, 1)  # rho axis fixed to 0-1
    ax_nu.set_ylim(mid_zeta - max_range, mid_zeta + max_range)
    ax_nu.set_zlim(mid_alpha - max_range, mid_alpha + max_range)
    
    # Add legend
    ax_nu.legend()
    
    plt.tight_layout()

    # Choose slices: use at most the first 6 phi/zeta slices
    slice_indices = list(range(min(nz, 6)))

    # Create a single figure with multiple subplots (2 rows x 3 columns)
    n_slices = len(slice_indices)
    nrows, ncols = 2, 3
    fig2, axs2 = plt.subplots(nrows, ncols, figsize=(16, 8))
    axs2 = axs2.flatten()


    # Global R/Z limits so every subplot has identical bounds, with a small margin
    R_min, R_max = R3.min(), R3.max()
    Z_min, Z_max = Z3.min(), Z3.max()
    R_pad = 0.05 * (R_max - R_min)
    Z_pad = 0.05 * (Z_max - Z_min)
    # Global span (with padding) that will be used for all subplots
    R_span = (R_max - R_min) + 2 * R_pad
    Z_span = (Z_max - Z_min) + 2 * Z_pad

    for idx, k in enumerate(slice_indices):
        ax2 = axs2[idx]

        # Use all nodes on this (i, j, k) slice
        Rk = R3[:, :, k].ravel()
        Zk = Z3[:, :, k].ravel()

        # Scatter of the nodes at this phi/zeta
        ax2.scatter(Rk, Zk, s=10, alpha=0.9, color='blue', label='nodes')

        # Draw line segments between neighbours in rho / alpha directions
        for i in range(nx):
            for j in range(ny):
                R0 = R3[i, j, k]
                Z0 = Z3[i, j, k]

                # Neighbour in +rho direction (i+1, j, k)
                if i + 1 < nx:
                    ax2.plot(
                        [R0, R3[i + 1, j, k]],
                        [Z0, Z3[i + 1, j, k]],
                        color='black',
                        linewidth=0.5,
                    )

                # Neighbour in +alpha direction (periodic in alpha: wrap last to first)
                j_next = (j + 1) % ny
                ax2.plot(
                    [R0, R3[i, j_next, k]],
                    [Z0, Z3[i, j_next, k]],
                    color='black',
                    linewidth=0.5,
                )

        # Title with phi value for this slice (use mean value on the slice)
        phi_k = np.mean(((phi3[:, :, k] + np.pi) % (2 * np.pi) - np.pi))
        ax2.set_xlabel('R [m]')
        ax2.set_ylabel('Z [m]')
        ax2.set_title(f'Bmag for nodes at $\\phi = {phi_k:.3f}$ rad')

        # Enforce identical scale on all subplots, but center each slice
        R_center = 0.5 * (R3[:, :, k].min() + R3[:, :, k].max())
        Z_center = 0.5 * (Z3[:, :, k].min() + Z3[:, :, k].max())
        ax2.set_xlim(R_center - 0.5 * R_span, R_center + 0.5 * R_span)
        ax2.set_ylim(Z_center - 0.5 * Z_span, Z_center + 0.5 * Z_span)
        ax2.set_aspect('equal', 'box')
        ax2.legend(loc='best')

    plt.tight_layout()

plt.show()
