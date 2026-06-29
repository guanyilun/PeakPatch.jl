import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.interpolate
import glob
import os
import sys 

halofile    = "data/halos.pksc"
abmatchfile = "data/abundance_match_table_10000Mbins_46zbins.npz"

h      = 0.68
omegam = 0.31

rho    = 2.775e11 * omegam * h**2 

ncellcut = 10.
mcellcut = (7700./6144)**3 * rho *ncellcut

Mmin = mcellcut
Mmax = 1.e16
print "Mmin = ", Mmin


#distance to redshift
zinterp=np.linspace(0,5,10000)
def hubble(z):
    return h*100*np.sqrt(omegam*(1+z)**3+1-omegam)

def drdz(z):
    return 299792.458 / hubble(z)

rinterp  = np.cumsum(drdz(zinterp)*(zinterp[1]-zinterp[0]))
rinterp -= rinterp[0]
z_to_r   = sp.interpolate.interp1d(zinterp,rinterp)
r_to_z   = sp.interpolate.interp1d(rinterp,zinterp)

#GET ppdata
def get_pp_catalogue(filename,Mmin,Mmax):

    pkfile   = open(filename,"rb")
    Non      = np.fromfile(pkfile, dtype=np.int32, count=1)[0]
    RTHMAXin = np.fromfile(pkfile, dtype=np.float32, count=1)
    redshift = np.fromfile(pkfile, dtype=np.float32, count=1)

    outnum   = 11
    npkdata  = outnum*Non

    peakdata = np.fromfile(pkfile, dtype=np.float32, count=npkdata)
    peakdata = np.reshape(peakdata,(Non,outnum))

    peakdata[:,6] = 4.0/3*np.pi * peakdata[:,6]**3 * rho

    dm  = (peakdata[:,6] > Mmin) & (peakdata[:,6] < Mmax) #cut mass
    peakdata = peakdata[dm]

    rpp   = np.sqrt(peakdata[:,0]**2+peakdata[:,1]**2+peakdata[:,2]**2)
    Zpp   = r_to_z(rpp)
    
    print "number of peak patch halos = ", len(Zpp)

    Non = np.array(len(Zpp))
    return Non, RTHMAXin, redshift, peakdata, Zpp

M_data            = np.load(abmatchfile)
M_interp          = M_data['M']
z_interp          = M_data['z']
Mvir_of_MTH_and_z = M_data['M200_of_MTH_and_z']
MTH_to_Mtinker    = scipy.interpolate.RectBivariateSpline(M_interp,z_interp,M200_of_MTH_and_z,kx=1,ky=1)

#load halos
Non, RTHMAXin, redshift, peakdata, Zpp = get_pp_catalogue(halofile,Mmin,Mmax)

# abundance match halos
M_pp           = peakdata[:,6].copy()
M_pp_correct   = MTH_to_Mtinker.ev(M_pp,Zpp)

# convert back to Rth
peakdata[:,6]  = M_pp_correct
peakdata[:,6] = (3./4/np.pi/rho * peakdata[:,6])**(1./3)

# save abundancematched catalogue
outfile = open('data/'+halofile+'_abundancematched',"wb")

Non      = Non.astype('int32')
RTHMAXin = RTHMAXin.astype('float32')
redshift = redshift.astype('float32')
peakdata = peakdata.flatten().astype('float32')

Non.tofile(outfile)
RTHMAXin.tofile(outfile)
redshift.tofile(outfile)
peakdata.tofile(outfile)
