import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import postgkyl as pg

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <name>")
    sys.exit(1)
name = sys.argv[1]
geom_dir = f"{name}_geometry"

def basis(x,y,z):
    return 1/2**1.5, np.sqrt(3)*x/2**1.5, np.sqrt(3)*y/2**1.5, np.sqrt(3)*z/2**1.5, 3*x*y/2**1.5, 3*x*z/2**1.5, 3*y*z/2**1.5, x*y*z*(3/2)**1.5

rzpdata_corner = pg.GData(f'./{geom_dir}/{name}-nodes.gkyl')
razdata_corner = pg.GData(f'./{geom_dir}/{name}-nodes-computational.gkyl')

rzp = rzpdata_corner.get_values()
raz = razdata_corner.get_values()

print(rzp.shape)
R = rzp[:,:,:,0]
Z = rzp[:,:,:,1]
phi = rzp[:,:,:,2]

rho = raz[:,:,:,0]
alpha = raz[:,:,:,1]
zeta = raz[:,:,:,2]

bdata = pg.GData(f'./{geom_dir}/{name}-bmag_corn.gkyl')
bcoeffs = bdata.get_values()
print(bcoeffs.shape)
bvals = np.sum(bcoeffs*basis(0,0,-1), axis=-1)
print(bvals.shape)

Jdata = pg.GData(f'./{geom_dir}/{name}-jacobgeo.gkyl')
Jcoeffs = Jdata.get_values()
print(Jcoeffs.shape)
Jvals = np.sum(Jcoeffs*basis(0,0,-1), axis=-1)
print(Jvals.shape)

# Create a modified inferno colormap with brighter dark values
# Cut off the darkest part so the minimum color is dark purple instead of black
try:
    # For newer matplotlib versions
    inferno = plt.colormaps['inferno']
except (AttributeError, KeyError):
    # Fallback for older versions
    inferno = plt.cm.get_cmap('inferno')
# Sample from inferno starting at 0.1 (dark purple) instead of 0 (black)
# This cuts off the darkest 10% of the colormap
cutoff = 0.15  # Fraction of darkest colors to skip
inferno_colors = inferno(np.linspace(cutoff, 1, 256))
inferno_modified = mcolors.LinearSegmentedColormap.from_list('inferno_modified', inferno_colors)

# Global R/Z limits so every subplot has identical bounds, with a small margin
R_min, R_max = R.min(), R.max()
Z_min, Z_max = Z.min(), Z.max()
R_pad = 0.05 * (R_max - R_min)
Z_pad = 0.05 * (Z_max - Z_min)
# Global span (with padding) that will be used for all subplots
R_span = (R_max - R_min) + 2 * R_pad
Z_span = (Z_max - Z_min) + 2 * Z_pad

# Use the maximum span to ensure equal physical scales on both axes
max_span = max(R_span, Z_span)

# Get shape for determining nx, ny
nx, ny = R[:,:,0].shape
rho_tol = 1e-6  # Tolerance for checking if rho is zero

# Create figure for all |B| plots (2 rows x 3 columns = 6 subplots)
fig_bmag = plt.figure(figsize=(18, 12))
axs_bmag = fig_bmag.subplots(2, 3)

# Create figure for all Jacobian plots (2 rows x 3 columns = 6 subplots)
fig_jacob = plt.figure(figsize=(18, 12))
axs_jacob = fig_jacob.subplots(2, 3)

# Flatten the axes arrays for easier indexing
axs_bmag = axs_bmag.flatten()
axs_jacob = axs_jacob.flatten()

# Loop over zeta indices 0-5 to create 6 plots in each figure
for idx, zeta_idx in enumerate(range(6)):
    # Get the axes for this subplot
    ax_bmag = axs_bmag[idx]
    ax_jacob = axs_jacob[idx]
    
    # ========== B-field magnitude plot ==========
    pcm0 = ax_bmag.pcolor(R[:,:,zeta_idx], Z[:,:,zeta_idx], bvals[:,:,zeta_idx], cmap = inferno_modified)
    ax_bmag.scatter(R[:,:,zeta_idx], Z[:,:,zeta_idx], c='black', s=10, marker='o', edgecolors='none')

    # Draw line segments between neighbours in rho / alpha directions
    for i in range(nx):
        for j in range(ny):
            R0 = R[i, j, zeta_idx]
            Z0 = Z[i, j, zeta_idx]
            rho0 = rho[i, j, zeta_idx]
            
            # Neighbour in +rho direction (i+1, j, zeta_idx)
            if i + 1 < nx:
                rho1 = rho[i + 1, j, zeta_idx]
                # Only draw line if neither point is at rho = 0
                if abs(rho0) > rho_tol and abs(rho1) > rho_tol:
                    ax_bmag.plot(
                        [R0, R[i + 1, j, zeta_idx]],
                        [Z0, Z[i + 1, j, zeta_idx]],
                        color='black',
                        linewidth=0.5,
                    )
            
            # Neighbour in +alpha direction (i, j+1, zeta_idx)
            if j + 1 < ny:
                rho1 = rho[i, j + 1, zeta_idx]
                # Only draw line if neither point is at rho = 0
                if abs(rho0) > rho_tol and abs(rho1) > rho_tol:
                    ax_bmag.plot(
                        [R0, R[i, j + 1, zeta_idx]],
                        [Z0, Z[i, j + 1, zeta_idx]],
                        color='black',
                        linewidth=0.5,
                    )

    # Enforce identical physical scale on both axes (1 meter in R = 1 meter in Z)
    R_center = 0.5 * (R[:, :, zeta_idx].min() + R[:, :, zeta_idx].max())
    Z_center = 0.5 * (Z[:, :, zeta_idx].min() + Z[:, :, zeta_idx].max())
    ax_bmag.set_xlim(R_center - 0.5 * max_span, R_center + 0.5 * max_span)
    ax_bmag.set_ylim(Z_center - 0.5 * max_span, Z_center + 0.5 * max_span)
    ax_bmag.set_aspect('equal', 'box')
    ax_bmag.set_xlabel('R [m]')
    ax_bmag.set_ylabel('Z [m]')

    # Title with phi value for this slice
    phi_k = np.mean(((phi[:, :, zeta_idx] + np.pi) % (2 * np.pi) - np.pi))
    ax_bmag.set_title(f'{name} |B| at $\\phi = {phi_k:.3f}$ rad')

    # Add colorbar for this subplot
    cbar0 = fig_bmag.colorbar(pcm0, ax=ax_bmag)
    cbar0.set_label('|B|')

    # ========== Jacobian plot ==========
    pcm1 = ax_jacob.pcolor(R[:,:,zeta_idx], Z[:,:,zeta_idx], Jvals[:,:,zeta_idx], cmap = inferno_modified)
    ax_jacob.scatter(R[:,:,zeta_idx], Z[:,:,zeta_idx], c='black', s=10, marker='o', edgecolors='none')

    # Draw line segments between neighbours in rho / alpha directions
    for i in range(nx):
        for j in range(ny):
            R0 = R[i, j, zeta_idx]
            Z0 = Z[i, j, zeta_idx]
            rho0 = rho[i, j, zeta_idx]
            
            # Neighbour in +rho direction (i+1, j, zeta_idx)
            if i + 1 < nx:
                rho1 = rho[i + 1, j, zeta_idx]
                # Only draw line if neither point is at rho = 0
                if abs(rho0) > rho_tol and abs(rho1) > rho_tol:
                    ax_jacob.plot(
                        [R0, R[i + 1, j, zeta_idx]],
                        [Z0, Z[i + 1, j, zeta_idx]],
                        color='black',
                        linewidth=0.5,
                    )
            
            # Neighbour in +alpha direction (i, j+1, zeta_idx)
            if j + 1 < ny:
                rho1 = rho[i, j + 1, zeta_idx]
                # Only draw line if neither point is at rho = 0
                if abs(rho0) > rho_tol and abs(rho1) > rho_tol:
                    ax_jacob.plot(
                        [R0, R[i, j + 1, zeta_idx]],
                        [Z0, Z[i, j + 1, zeta_idx]],
                        color='black',
                        linewidth=0.5,
                    )

    # Enforce identical physical scale on both axes (1 meter in R = 1 meter in Z)
    ax_jacob.set_xlim(R_center - 0.5 * max_span, R_center + 0.5 * max_span)
    ax_jacob.set_ylim(Z_center - 0.5 * max_span, Z_center + 0.5 * max_span)
    ax_jacob.set_aspect('equal', 'box')
    ax_jacob.set_xlabel('R [m]')
    ax_jacob.set_ylabel('Z [m]')

    # Title with phi value for this slice
    ax_jacob.set_title(f'{name} Jacobian at $\\phi = {phi_k:.3f}$ rad')

    # Add colorbar for this subplot
    cbar1 = fig_jacob.colorbar(pcm1, ax=ax_jacob)
    cbar1.set_label('Jacobian')

# Adjust layout for both figures
fig_bmag.tight_layout()
fig_jacob.tight_layout()

plt.show()

