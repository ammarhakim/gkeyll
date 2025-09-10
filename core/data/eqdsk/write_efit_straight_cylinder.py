import postgkyl as pg
import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.integrate as integrate
from scipy.interpolate import pchip_interpolate, interp1d
from scipy.special import ellipk, ellipe
import scipy.constants
import scipy.optimize as sco
import fortranformat as ff
from datetime import date
from typing import Any, Generator, Iterable, List, TextIO, Union

outFileName = 'straight_cylinder.geqdsk'

B0 = 0.0
R0 = .02

def psi_f(R, Z):
    Bmag = 0.5
    return Bmag/2 * R**2

#RZ box
NW = 257
NH = 257
RMIN,RMAX = .001, 1
ZMIN,ZMAX = -1, 1
RDIM = RMAX - RMIN
ZDIM = ZMAX - ZMIN
RLEFT = RMIN
ZMID = (ZMAX+ZMIN)/2.0
RMAXIS = 0.0
ZMAXIS = 0.0
SIMAG = psi_f(2.0,0) 
SIBRY = psi_f(4.0,0)
NPSI = NW #don't write
RCENTR = 0.0
BCENTR = B0
CURRENT = 0

print("simag = %g"%SIMAG)
print("sibry = %g"%SIBRY)

#Solve GS in RZ coords
Rgrid = np.linspace(RMIN,RMAX,NW)
Zgrid = np.linspace(ZMIN,ZMAX,NH)
dR = Rgrid[1] - Rgrid[0]
dZ = Zgrid[1] - Zgrid[0]
#rthetagrid = np.zeros((len(Rgrid),len(Zgrid),2))
psiRZ = np.zeros((len(Rgrid),len(Zgrid)))
BzRZ = np.zeros((len(Rgrid),len(Zgrid)))
BrRZ = np.zeros((len(Rgrid),len(Zgrid)))
BRZ = np.zeros((len(Rgrid),len(Zgrid)))
for i,Ri in enumerate(Rgrid):
    for j,Zj in enumerate(Zgrid):
        psiRZ[i,j] = psi_f(Ri,Zj)
        if Ri < 1e-17:
            Ri = Rgrid[1]
        BzRZ[i,j] = 1/Ri * (psi_f(Ri+dR,Zj) - psi_f(Ri-dR,Zj))/(2*dR)
        BrRZ[i,j] = -1/Ri * (psi_f(Ri,Zj+dZ) - psi_f(Ri,Zj-dZ))/(2*dZ)
        BRZ[i,j] = np.sqrt(BrRZ[i,j]**2 + BzRZ[i,j]**2)

plt.figure()
contour = plt.contour(Rgrid,Zgrid, psiRZ.T)
plt.xlabel('R')
plt.ylabel('Z')
plt.colorbar()
plt.title("Psi calculated in RZ coords")
plt.savefig('straight_cylinder_psi_RZ.png')
plt.close()

plt.figure()
plt.pcolormesh(Rgrid,Zgrid, BRZ.T)
plt.xlabel('R')
plt.ylabel('Z')
plt.colorbar()
plt.title("B calculated in RZ coords")
plt.savefig('straight_cylinder_B_RZ.png')
plt.close()

# Find where Rgrid = 0.04 within some tolerance
Rind = np.argmin(np.abs(Rgrid - 0.5))

plt.figure()
plt.plot(Zgrid, BRZ[Rind,:])
plt.xlabel('Z')
plt.ylabel('B')
plt.title("B calculated in RZ coords at R = "+str(Rgrid[Rind]))
plt.savefig('straight_cylinder_B_RZ_at_R.png')
plt.close()

#PSI quantities
PSIGRID = np.linspace(SIMAG, SIBRY,NPSI)
FPOL = (B0*R0/Rgrid)*Rgrid # F = RBphi
FFPRIM = np.repeat(0.0, NPSI)
PPRIME = np.repeat(0.0,NPSI)
PRES = integrate.cumtrapz(PPRIME,PSIGRID,initial=0)
PSIZR = psiRZ.T


QPSI = np.zeros(NPSI)
for i in range(NPSI):
    QPSI[i] = 0



writeList = [NW, NH,                                        #3i4
             RDIM, ZDIM, RCENTR, RLEFT, ZMID,               #5E16.9
             RMAXIS, ZMAXIS, SIMAG, SIBRY, BCENTR,          #5e16.9
             CURRENT, SIMAG, 0, RMAXIS, 0,                  #5E16.9
             ZMAXIS, 0, SIBRY, 0, 0,                        #5E16.9
             FPOL, PRES, FFPRIM, PPRIME,                    #5E16.9
             PSIZR, QPSI]                                   #5E16.9

#Header stuff
header_fmt = "(a48,3i4)"
label = 'FREEGS'
creation_date = date.today().strftime("%d/%m/%Y")
shot = int(0)
time = int(0)
shot_str = f"# {shot:d}"
time_str = f"  {time:d}ms"
comment = f"{label:11}{creation_date:10s}   {shot_str:>8s}{time_str:16s}"
def write_line(data: Iterable[Any], fh: TextIO, fmt: str) -> None:
    r"""
    Writes to a Fortran formatted ASCII data file. The file handle will be left on a
    newline.

    Parameters
    ---------
    data:
        The data to write.
    fh:
        File handle. Should be in a text write mode, i.e. ``open(filename, "w")``.
    fmt:
        A Fortran IO format string, such as ``'(6a8,3i3)'``.
    """
    fh.write(ff.FortranRecordWriter(fmt).write(data))
    fh.write("\n")
#write_line((comment, 3, NW, NH), fh, header_fmt) #3 is idum


#Now write the EFIT FILE
with open(outFileName,'w',newline='') as f:
    ##NW,NH
    ##for i in range(2):
    ##    f.write('%d '%writeList[i])
    #f.write('FREEGS     19/06/2023        # 0  0ms              3  80  90')
    #f.write('\n')
    write_line((comment, 3, NW, NH), f, header_fmt) #3 is idum
    # rdim,zdim,rcentr,rleft,zmid
    for i in range(2,7):
        f.write('%16.9E'%writeList[i])
    f.write('\n')
    # rmaxis,zmaxis,simag,sibry,bcentr
    for i in range(7,12):
        f.write('%16.9E'%writeList[i])
    f.write('\n')
    # current, simag, xdum, rmaxis, xdum
    for i in range(12,17):
        f.write('%16.9E'%writeList[i])
    f.write('\n')
    #zmaxis,xdum,sibry,xdum,xdum
    for i in range(17,22):
         f.write('%16.9E'%writeList[i])
    f.write('\n')
    #FPOL,PRES,FFPRIM,PPRIME
    for i in range(22,26):
        count=0
        for j in range(0,NW):
            f.write('%16.9E'%writeList[i][j])
            count = count+1
            if count==5:
                f.write('\n')
                count=0
    #PSIZR
    for i in range(26,27):
        count = 0
        for j in range(0,NH):
            for k in range(0,NW):
                f.write('%16.9E'%writeList[i][j][k])
                count = count+1
                if count==5:
                    f.write('\n')
                    count=0
    #QPSI
    for i in range(27,28):
        count=0
        for j in range(0,NW):
            f.write('%16.9E'%writeList[i][j])
            count = count+1
            if count==5:
                f.write('\n')
                count=0
