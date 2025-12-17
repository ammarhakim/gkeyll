# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import postgkyl as pg
import os
import math
import sys

import scipy.interpolate 
from scipy.interpolate import RegularGridInterpolator, interp1d
import scipy.integrate as sci


import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable


def fix_gridvals(grid):
    """Output file grids have cell-edge coordinates by default, but the values are at cell centers.
    This function will return the cell-center coordinate grid.
    Usage: cell_center_grid = fix_gridvals(cell_edge_grid)
    """
    grid = np.array(grid).squeeze()
    grid = (grid[0:-1] + grid[1:])/2
    return grid

def myinterpolate(raw_x,raw_z,x,z,data):
    myinterpolator = RegularGridInterpolator((raw_x,raw_z), data, bounds_error=False, fill_value=None)
    x0,z0 = np.meshgrid(x,z)
    return myinterpolator((x0,z0)).T


data = {}

#Load geometry

zeta_idx = 0

bcartdata = pg.GData('./W7X-nodal_modal_arrays/modal/bcart_interior_modal.gkyl')
grid,b0val = pg.data.GInterpModal(bcartdata,poly_order=1,basis_type='ms').interpolate(0)
Bx = b0val.squeeze()
grid,b1val = pg.data.GInterpModal(bcartdata,poly_order=1,basis_type='ms').interpolate(1)
By = b1val.squeeze()
grid,b2val = pg.data.GInterpModal(bcartdata,poly_order=1,basis_type='ms').interpolate(2)
Bz = b2val.squeeze()

bcartdata = pg.GData('./W7X-nodal_modal_arrays/modal/dxdz_interior_modal.gkyl')
grid,e_3_0val = pg.data.GInterpModal(bcartdata,poly_order=1,basis_type='ms').interpolate(6)
e_3x = e_3_0val.squeeze()
grid,e_3_1val = pg.data.GInterpModal(bcartdata,poly_order=1,basis_type='ms').interpolate(7)
e_3y = e_3_1val.squeeze()
grid,e_3_2val = pg.data.GInterpModal(bcartdata,poly_order=1,basis_type='ms').interpolate(8)
e_3z = e_3_2val.squeeze()



mc2pdata = pg.GData('./W7X-nodal_modal_arrays/modal/mc2p_corner_modal.gkyl')
grid,rvals = pg.data.GInterpModal(mc2pdata,poly_order=1,basis_type='ms').interpolate(0)
grid,zvals = pg.data.GInterpModal(mc2pdata,poly_order=1,basis_type='ms').interpolate(1)
grid,phivals = pg.data.GInterpModal(mc2pdata,poly_order=1,basis_type='ms').interpolate(2)


R = rvals.squeeze()
Z = zvals.squeeze()
phi = phivals.squeeze()

X = R * np.cos(phi)
Y = R * np.sin(phi)
Z = Z

from mpl_toolkits.mplot3d import Axes3D

downsample_rho = 3
downsample_alpha = 3
downsample_zeta = 2

# Downsample arrays for both B and e_3 (make sure your variables are not typos!)
X_ds = X[::downsample_rho, ::downsample_alpha, ::downsample_zeta]
Y_ds = Y[::downsample_rho, ::downsample_alpha, ::downsample_zeta]
Z_ds = Z[::downsample_rho, ::downsample_alpha, ::downsample_zeta]
Bx_ds = Bx[::downsample_rho, ::downsample_alpha, ::downsample_zeta]
By_ds = By[::downsample_rho, ::downsample_alpha, ::downsample_zeta]
Bz_ds = Bz[::downsample_rho, ::downsample_alpha, ::downsample_zeta]

# Downsample e_3 vector components
e_3x_ds = e_3x[::downsample_rho, ::downsample_alpha, ::downsample_zeta]
e_3y_ds = e_3y[::downsample_rho, ::downsample_alpha, ::downsample_zeta]
e_3z_ds = e_3z[::downsample_rho, ::downsample_alpha, ::downsample_zeta]

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Quiver for B vector field
ax.quiver(
    X_ds.flatten(), Y_ds.flatten(), Z_ds.flatten(),
    Bx_ds.flatten(), By_ds.flatten(), Bz_ds.flatten(),
    length=0.25, normalize=True, linewidth=0.5, color='tab:blue', label='$\\hat{\\mathbf{b}}$ field'
)

# Overlay quiver for e_3 vectors at the same nodes
ax.quiver(
    X_ds.flatten(), Y_ds.flatten(), Z_ds.flatten(),
    e_3x_ds.flatten(), e_3y_ds.flatten(), e_3z_ds.flatten(),
    length=0.18, normalize=True, linewidth=0.5, color='tab:orange', label='e₃ vectors'
)

# Add custom legend for quiver (since quiver does not add entry by itself)
import matplotlib.lines as mlines
quiver_b_legend = mlines.Line2D([], [], color='tab:blue', label='$\\hat{\\mathbf{b}}$ field', linewidth=2)
quiver_e3_legend = mlines.Line2D([], [], color='tab:orange', label='e₃ vectors', linewidth=2)
ax.legend(handles=[quiver_b_legend, quiver_e3_legend])


ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('3D Vector Field of $\\hat{\\mathbf{b}}$ and e₃')


plt.tight_layout()

plt.show()

