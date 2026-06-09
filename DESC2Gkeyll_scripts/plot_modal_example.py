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

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <name>")
    sys.exit(1)
name = sys.argv[1]
geom_dir = f"{name}_geometry"

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

bdata = pg.GData(f'./{geom_dir}/{name}-bmag.gkyl')
grid,val = pg.data.GInterpModal(bdata,poly_order=1,basis_type='ms').interpolate(0)
data["B"] = val.squeeze()
#data["B"] = np.column_stack((data["B"][:,:], data["B"][:,0]))

rho = fix_gridvals(grid[0])
alpha = fix_gridvals(grid[1])
print("rho shape:", rho.shape)
print("alpha shape:", alpha.shape)

bdata_corner = pg.GData(f'./{geom_dir}/{name}-bmag_corn.gkyl')
grid,val = pg.data.GInterpModal(bdata_corner,poly_order=1,basis_type='ms').interpolate(0)
data["Bcorn"] = val.squeeze()
#data["Bcorn"] = np.column_stack((data["Bcorn"][:,:], data["Bcorn"][:,0]))

jdata = pg.GData(f'./{geom_dir}/{name}-jacobgeo.gkyl')
grid,val = pg.data.GInterpModal(jdata,poly_order=1,basis_type='ms').interpolate(0)
data["J"] = val.squeeze()
#data["J"] = np.column_stack((data["J"][:,:], data["J"][:,0]))

mc2pdata = pg.GData(f'./{geom_dir}/{name}-nodes.gkyl')
grid,rvals = pg.data.GInterpModal(mc2pdata,poly_order=1,basis_type='ms').interpolate(0)
grid,zvals = pg.data.GInterpModal(mc2pdata,poly_order=1,basis_type='ms').interpolate(1)
grid,phivals = pg.data.GInterpModal(mc2pdata,poly_order=1,basis_type='ms').interpolate(2)


data["R"] = rvals.squeeze()
data["Z"] = zvals.squeeze()
data["phi"] = phivals.squeeze()

#data["R"] = np.column_stack((data["R"][:,:], data["R"][:,0]))
#data["Z"] = np.column_stack((data["Z"][:,:], data["Z"][:,0]))
#data["phi"] = np.column_stack((data["phi"][:,:], data["phi"][:,0]))



#plt.scatter(data["R"][:,:,zeta_idx].flatten(), data["Z"][:,:,zeta_idx].flatten(), c=data["J"][:,:,zeta_idx].flatten(), cmap = 'inferno')

fig, axs = plt.subplots(1, 2, figsize=(12, 5))

pcm0 = axs[0].pcolor(data["R"][:,:,zeta_idx], data["Z"][:,:,zeta_idx], data["J"][:,:,zeta_idx], cmap = 'inferno')
#pcm0 = axs[0].scatter(data["R"][:,:].flatten(), data["Z"][:,:].flatten(), c=data["B"][:,:].flatten(), cmap = 'inferno')

axs[0].set_title("Bmag from interior points")
fig.colorbar(pcm0, ax=axs[0], orientation='vertical', label='B')

pcm1 = axs[1].pcolormesh(data["R"][:,:,zeta_idx], data["Z"][:,:,zeta_idx], data["Bcorn"][:,:,zeta_idx], cmap = 'inferno')
#pcm1 = axs[1].scatter(data["R"][:,:].flatten(), data["Z"][:,:].flatten(), c=data["Bcorn"][:,:].flatten(), cmap = 'inferno')

axs[1].set_title("Bmag from corner points")
fig.colorbar(pcm1, ax=axs[1], orientation='vertical', label='Bcorn')


for ax in axs:
    ax.set_xlabel('R')
    ax.set_ylabel('Z')

"""
fig, axs = plt.subplots(1, 2, figsize=(12, 5))

pcm0 = axs[0].pcolor(rho, alpha, data["B"][:,:,zeta_idx].T, cmap='inferno')
axs[0].set_title("Bmag from interior points")
fig.colorbar(pcm0, ax=axs[0], orientation='vertical', label='B')

pcm1 = axs[1].pcolor(rho, alpha, data["Bcorn"][:,:,zeta_idx].T, cmap='inferno')
axs[1].set_title("Bmag from corner points")
fig.colorbar(pcm1, ax=axs[1], orientation='vertical', label='Bcorn')

for ax in axs:
    ax.set_xlabel('rho')
    ax.set_ylabel('alpha')
"""

#plt.colorbar()
plt.show()

