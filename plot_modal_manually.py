import numpy as np
import matplotlib.pyplot as plt
import postgkyl as pg

def basis(x,y,z):
    return 1/2**1.5, np.sqrt(3)*x/2**1.5, np.sqrt(3)*y/2**1.5, np.sqrt(3)*z/2**1.5, 3*x*y/2**1.5, 3*x*z/2**1.5, 3*y*z/2**1.5, x*y*z*(3/2)**1.5

rzpdata_corner = pg.GData('./W7X-nodal_modal_arrays/nodal/mc2p_corner_nodal.gkyl')
razdata_corner = pg.GData('./W7X-nodal_modal_arrays/nodal/mc2nu_pos_corner_nodal.gkyl')

rzp = rzpdata_corner.get_values()
raz = razdata_corner.get_values()

print(rzp.shape)
R = rzp[:,:,:,0]
Z = rzp[:,:,:,1]
phi = rzp[:,:,:,2]

rho = raz[:,:,:,0]
alpha = raz[:,:,:,1]
zeta = raz[:,:,:,2]

bdata = pg.GData('./W7X-nodal_modal_arrays/modal/bmag_corner_modal.gkyl')
bcoeffs = bdata.get_values()
print(bcoeffs.shape)
bvals = np.sum(bcoeffs*basis(0,0,-1), axis=-1)
print(bvals.shape)

Jdata = pg.GData('./W7X-nodal_modal_arrays/modal/jacobgeo_interior_modal.gkyl')
Jcoeffs = Jdata.get_values()
print(Jcoeffs.shape)
Jvals = np.sum(Jcoeffs*basis(0,0,-1), axis=-1)
print(Jvals.shape)

zeta_idx = 0

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

# Create a square figure
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

pcm0 = ax.pcolor(R[:,:,zeta_idx], Z[:,:,zeta_idx], bvals[:,:,zeta_idx], cmap = 'inferno')
ax.scatter(R[:,:,zeta_idx], Z[:,:,zeta_idx], c='black', s=10, marker='o', edgecolors='none')

# Draw line segments between neighbours in rho / alpha directions
nx, ny = R[:,:,zeta_idx].shape
rho_tol = 1e-6  # Tolerance for checking if rho is zero
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
                ax.plot(
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
                ax.plot(
                    [R0, R[i, j + 1, zeta_idx]],
                    [Z0, Z[i, j + 1, zeta_idx]],
                    color='black',
                    linewidth=0.5,
                )

# Enforce identical physical scale on both axes (1 meter in R = 1 meter in Z)
R_center = 0.5 * (R[:, :, zeta_idx].min() + R[:, :, zeta_idx].max())
Z_center = 0.5 * (Z[:, :, zeta_idx].min() + Z[:, :, zeta_idx].max())
ax.set_xlim(R_center - 0.5 * max_span, R_center + 0.5 * max_span)
ax.set_ylim(Z_center - 0.5 * max_span, Z_center + 0.5 * max_span)
ax.set_aspect('equal', 'box')
ax.set_xlabel('R [m]')
ax.set_ylabel('Z [m]')

# Title with phi value for this slice
phi_k = np.mean(((phi[:, :, zeta_idx] + np.pi) % (2 * np.pi) - np.pi))
ax.set_title(f'W7-X |B| at $\\phi = {phi_k:.3f}$ rad')

# Add colorbar
cbar0 = fig.colorbar(pcm0, ax=ax)
cbar0.set_label('|B|')

# Create a second square figure for Jacobian
fig2 = plt.figure(figsize=(8, 8))
ax2 = fig2.add_subplot(111)

pcm1 = ax2.pcolor(R[:,:,zeta_idx], Z[:,:,zeta_idx], Jvals[:,:,zeta_idx], cmap = 'inferno')
ax2.scatter(R[:,:,zeta_idx], Z[:,:,zeta_idx], c='black', s=10, marker='o', edgecolors='none')

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
                ax2.plot(
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
                ax2.plot(
                    [R0, R[i, j + 1, zeta_idx]],
                    [Z0, Z[i, j + 1, zeta_idx]],
                    color='black',
                    linewidth=0.5,
                )

# Enforce identical physical scale on both axes (1 meter in R = 1 meter in Z)
ax2.set_xlim(R_center - 0.5 * max_span, R_center + 0.5 * max_span)
ax2.set_ylim(Z_center - 0.5 * max_span, Z_center + 0.5 * max_span)
ax2.set_aspect('equal', 'box')
ax2.set_xlabel('R [m]')
ax2.set_ylabel('Z [m]')

# Title with phi value for this slice
ax2.set_title(f'W7-X Jacobian at $\\phi = {phi_k:.3f}$ rad')

# Add colorbar
cbar1 = fig2.colorbar(pcm1, ax=ax2)
cbar1.set_label('Jacobian')

plt.show()

