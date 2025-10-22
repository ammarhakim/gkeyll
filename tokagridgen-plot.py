import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import postgkyl as pg
import scipy.interpolate 
from scipy.interpolate import RegularGridInterpolator
import matplotlib.lines as mlines
import sys

# Get name from command line argument
if len(sys.argv) < 2:
    print("Usage: python plotmb-nodes-step.py <name>")
    print("Available names: nstxu_DN, asdex, step")
    sys.exit(1)

name = sys.argv[1]  # Get name from command line argument
if name in ["nstxu_DN", "step"]:
    DN = True
else:
    DN = False


psid = pg.GData(name+"_psi.gkyl")
interp = pg.GInterpModal(psid,2,"mt")
grid, psi = interp.interpolate()

for d in range(len(grid)):
    grid[d] = 0.5*(grid[d][:-1] + grid[d][1:])



psisep = 1.5093065418975686 # change this to read separatrix value from the GData


fig, ax = plt.subplots(figsize = (4,9))



#Plot Separatrix
ax.contour(grid[0], grid[1], psi[:,:,0].transpose(), levels=np.r_[psisep], colors="r", linestyles='dashed')



colors = ["tab:orange","tab:blue","tab:green", "tab:brown", "tab:purple"]
sim_dir = "./"
baseName = sim_dir+name #'nstxu_DN or 'asdex'
bmin = 0
bmax = 8 if DN == True else 6   # 8 for DN and 6 for SN
simNames = ['%s_b%d'%(baseName,i) for i in range(bmin,bmax)]
Rlist = []
Zlist = []
jlist = []
jmax = 0.0
jmin = np.inf
for i, simName in enumerate(simNames):
    data = pg.GData(simName+"-nodes.gkyl")
    vals = data.get_values()
    R = vals[:,:,0]
    Z = vals[:,:,1]
    phi = vals[:,:,2]

    temp_nodal_grid = data.get_grid()
    nodal_grid = []
    for d in range(0,len(temp_nodal_grid),2):
        nodal_grid.append( np.linspace(temp_nodal_grid[d][0], temp_nodal_grid[d][-1], len(temp_nodal_grid[d])-1) )



    
    
    #Construct lines
    segs1 = np.stack((R,Z), axis=2)
    segs2 = segs1.transpose(1,0,2)

    #Plot all nodes and cell boundaries
    #ax.scatter(R,Z, marker=".")
    ax.add_collection(LineCollection(segs1, linewidth=0.4, color = colors[i%len(colors)]))
    ax.add_collection(LineCollection(segs2, linewidth=0.4, color = colors[i%len(colors)]))
    ax.plot(R,Z,marker=".", color="k", linestyle="none", markersize=2.0)

    # For double null configuration, plot mirrored data with Z -> -Z
    if DN == True:
        Zm = -Z
        segs1m = np.stack((R,Zm), axis=2)
        segs2m = segs1m.transpose(1,0,2)
        ax.add_collection(LineCollection(segs1m, linewidth=0.4, color = colors[i%len(colors)]))
        ax.add_collection(LineCollection(segs2m, linewidth=0.4, color = colors[i%len(colors)]))
        ax.plot(R,Zm,marker=".", color="k", linestyle="none", markersize=2.0)

  

    Rlist.append(R)
    Zlist.append(Z)




div_handle = mlines.Line2D([],[], color = 'r', label = "Divertor", linewidth=3.0)
sep_handle = mlines.Line2D([],[], color = 'r', label = "Separatrix", linestyle="dashed")
handles= [sep_handle, div_handle]
ax.grid()
ax.set_xlabel("R [m]")
ax.set_ylabel("Z [m]")
ax.axis("tight")
ax.axis("image")
ax.legend(handles = handles, loc = "center")
#ax.set_xlim(0.4,1.9)
#ax.set_ylim(-2.5,2.5)

fig.tight_layout()
plt.show()
