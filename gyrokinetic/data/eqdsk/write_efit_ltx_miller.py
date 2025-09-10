#import postgkyl as pg
import numpy as np
import matplotlib.pyplot as plt
import sys
#import pgkylUtil as pgu
import scipy.integrate as integrate
from scipy.interpolate import pchip_interpolate
import scipy.optimize as sco
import fortranformat as ff
from datetime import date
from typing import Any, Generator, Iterable, List, TextIO, Union
from scipy.interpolate import griddata


#Geometry and magnetic field.
R_axis     = 0.406051632        # [m]
B_axis     = 0.240606108        # [T]
R_LCFSmid  = 0.61183            # Major radius of the LCFS at the outboard midplane [m].
Rmid_min   = R_axis # Minimum midplane major radius of simulation box [m].
Rmid_max   = 0.7   # Maximum midplane major radius of simulation box [m].
R0         = 0.5*(Rmid_min+Rmid_max)    # Major radius of the simulation box [m].
a_mid      = R_LCFSmid-R_axis   # Minor radius at outboard midplane [m].
r0         = R0-R_axis          # Minor radius of the simulation box [m].
B0         = B_axis*(R_axis/R0) # Magnetic field magnitude in the simulation box [T].
q_LCFS     = 4.3153848          # Safety factor at the LCFS.
s_LCFS     = 2.6899871          # Magnetic shear at the LCFS. Should be ~6.6899871 but that makes qprofile change sign in the SOL.
kappa      = 1.3                # Elongation (=1 for no elongation).
delta      = 0.4                # Triangularity (=0 for no triangularity).
x_LCFS = R_LCFSmid - Rmid_min


#.Magnetic safety factor profile.
def qprofile(r):
    return q_LCFS/(1.-s_LCFS*((r-a_mid)/a_mid))

q0 = qprofile(r0)
q_min = qprofile(Rmid_min-R_axis)

#..................... NO MORE USER INPUTS BELOW (maybe) ....................#

def r_x(x):
    return Rmid_min-R_axis + x
def R_f(r, theta):
  #return R_axis + r*np.cos(theta + np.arcsin(delta)*np.sin(theta))
  return R_axis + np.outer(r,np.cos(theta + np.arcsin(delta)*np.sin(theta)))
def Z_f(r, theta):
  #return kappa*r*np.sin(theta)
  return kappa*np.outer(r,np.sin(theta))

#.Analytic derivatives.
def R_f_r(r,theta):
  return np.cos(theta + np.arcsin(delta)*np.sin(theta))
def R_f_theta(r,theta):
  return -r*(np.arcsin(delta)*np.cos(theta)+1.)*np.sin(np.arcsin(delta)*np.sin(theta)+theta)
def Z_f_r(r,theta):
  return kappa*np.sin(theta)
def Z_f_theta(r,theta):
  return kappa*r*np.cos(theta)

def Jr_f(r, theta):
  return R_f(r,theta)*( R_f_r(r,theta)*Z_f_theta(r,theta)-Z_f_r(r,theta)*R_f_theta(r,theta) )

def J_f(r, theta):
    return Jr_f(r, theta)/dPsidr_f(r, theta)

def integrand(t, r):
  return Jr_f(r,t)/np.power(R_f(r,t),2)

def dPsidr_f(r, theta):
  integral, _ = integrate.quad(integrand, 0., 2.*np.pi, args=(r), epsabs=1.e-8)
  return B0*R_axis/(2.*np.pi*qprofile(r))*integral

def Bphi_f(R):
    return B0*R0/R

def psi_fi(r,theta):
    psi_i,_ = integrate.quad(dPsidr_f, 0, r, args=(theta), epsabs=1.e-8)
    return psi_i





#RZ box
NW = 81
NH = 91
RMIN,RMAX = 0.1,0.7
ZMIN,ZMAX = -0.4,0.4
RDIM = RMAX - RMIN
ZDIM = ZMAX - ZMIN
RLEFT=RMIN
ZMID = (ZMAX+ZMIN)/2.0
RMAXIS = R_axis
ZMAXIS = 0.0
SIMAG = psi_fi(0.0,0.0)
SIBRY = psi_fi(Rmid_min - R_axis + x_LCFS, 0.0)
NPSI = NW #don't write
RCENTR = R_axis
BCENTR = B_axis
CURRENT = 0


# setup r, theta grids
xmin = 0
xmax = Rmid_max-Rmid_min
x = np.linspace(0,xmax,NW)
r = Rmid_min-R_axis+x
theta=np.linspace(-np.pi,np.pi,65)
Rm = R_f(r,theta)
Zm = Z_f(r,theta)

#Calculate psi(r,theta)
psi = np.zeros(len(r))
for i,ri in enumerate(r):
    psi[i] = psi_fi(ri,0.0)
psi_plot = np.repeat(psi,len(theta)).reshape(len(r),len(theta))
fig, ax = plt.subplots()
cax = ax.contour(Rm,Zm,psi_plot, cmap = "inferno")
#cax = ax.pcolor(Rm,Zm,psi_plot, cmap = "inferno")
plt.title("Original Psi")
plt.colorbar(cax)
plt.show()




#Interpolate to get psi in RZ coords
Rgrid = np.linspace(RMIN,RMAX,NW)
Zgrid = np.linspace(ZMIN,ZMAX,NH)
grid_r, grid_z = np.meshgrid(Rgrid,Zgrid)
flat_points = np.stack((Rm,Zm), axis = -1).reshape(-1, 2)
flat_psi = psi_plot.flatten()
psiRZ = griddata(flat_points, psi_plot.flatten(), (grid_r, grid_z), method='cubic').T

fig, ax = plt.subplots()
cax = ax.contour(Rgrid,Zgrid, psiRZ.T, cmap = "inferno")
#cax = ax.pcolor(Rgrid, Zgrid, psiRZ.T, cmap = "inferno")
plt.title("Interpolated Psi")
plt.colorbar(cax)
plt.show()




#PSI quantities
PSIGRID = np.linspace(SIMAG, SIBRY,NPSI)
FPOL = np.repeat(B_axis, NPSI)
FFPRIM = np.repeat(0.0, NPSI)
PPRIME = np.repeat(-1e-6,NPSI)
PRES = integrate.cumulative_trapezoid(PPRIME,PSIGRID,initial=0)
PSIZR = psiRZ.T



#LIMITER and boundary STUFF
r_LCFS = r_x(x_LCFS)
RLCFS = R_f(r_LCFS,theta)
ZLCFS = Z_f(r_LCFS,theta)
plt.scatter(RLCFS,ZLCFS)

#QPSI = qprofile(r)
def rlossfunc(r, psi0):
    return psi0 - psi_fi(r,0)
def rfunc(psi):
    return sco.ridder(rlossfunc,0, r_LCFS, args = (psi) )
QPSI = np.zeros(NPSI)
for i in range(NPSI):
    QPSI[i] = qprofile(rfunc(PSIGRID[i]))



NBBBS = 65
LIMITR = 10
RBBBS = Rm[-1,:]
ZBBBS = Zm[-1,:]
#RLCFS = R_f(r_LCFS,theta).squeeze()
#ZLCFS = Z_f(r_LCFS,theta).squeeze()
RLIM = np.linspace(RLEFT,R_f(r_LCFS,np.pi), LIMITR).squeeze()
ZLIM = np.repeat(0.0,LIMITR)


plt.figure()
plt.scatter(RLIM,ZLIM)
plt.scatter(RBBBS,ZBBBS)
#plt.show()

writeList = [NW, NH,                                        #3i4
             RDIM, ZDIM, RCENTR, RLEFT, ZMID,               #5E16.9
             RMAXIS, ZMAXIS, SIMAG, SIBRY, BCENTR,          #5e16.9
             CURRENT, SIMAG, 0, RMAXIS, 0,                  #5E16.9
             ZMAXIS, 0, SIBRY, 0, 0,                        #5E16.9
             FPOL, PRES, FFPRIM, PPRIME,                    #5E16.9
             PSIZR, QPSI,                                   #5E16.9
             NBBBS, LIMITR,                                 #2i5
             RBBBS, ZBBBS,                                  #5E16.9
             RLIM, ZLIM]                                    #5E16.9

#Header stuff
header_fmt = "(a48,3i4)"
label='FREEGS'
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
with open('ltx_miller.geqdsk','w',newline='') as f:
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
    #NBBBS,LIMITR
    for i in range(28,30):
        #f.write('  %d  '%writeList[i])
        f.write('%5i'%writeList[i])
    f.write('\n')
    #RBBBS,ZBBBS
    for i in range(30,32):
        count=0
        for j in range(0,NBBBS):
            f.write('%16.9E'%writeList[i][j])
            count = count+1
            if count==5:
                f.write('\n')
                count=0
    #RLIM,ZLIM
    for i in range(32,34):
        count=0
        for j in range(0,LIMITR):
            f.write('%16.9E'%writeList[i][j])
            count = count+1
            if count==5:
                f.write('\n')
                count=0

#WRITE the Wall and Limiter if desired.
#WallCoords = np.c_[RBBBS,ZBBBS]
#with open('miller_wall','w',newline='') as f:
#    for i in range(NBBBS):
#        f.write('%16.9E'%WallCoords[i,0])
#        f.write('%16.9E'%WallCoords[i,1])
#        f.write('\n')
#
#LimCoords = np.c_[RLIM,ZLIM]
#with open('miller_limiter','w',newline='') as f:
#    for i in range(LIMITR):
#        f.write('%16.9E'%LimCoords[i,0])
#        f.write('%16.9E'%LimCoords[i,1])
#        f.write('\n')



plt.show()
