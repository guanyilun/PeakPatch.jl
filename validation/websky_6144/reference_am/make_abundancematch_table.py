import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.optimize import curve_fit

from matplotlib import gridspec
import matplotlib.colors as colors

import hmf
from astropy import units as u

import sys
import os

h      = 0.68
omegab = 0.043
omegam = 0.31
omegal = 1-omegam
sigma8 = 0.81
ns     = 0.965

rho_mean    = 2.775e11 * omegam * h**2 

bx = 7700.

halofile  = "data/halos.pksc" #needed if Nhalofile does not exist (contains halo catalogue N(M,z) table)
Nhalofile = "data/Nhalos_websky.txt"
pkfile    = "data/planck2018_powerspectrum.dat"

# Load in power spectrum used for sim
pk_data = np.loadtxt(pkfile)
k  = pk_data[:,0]
pk = pk_data[:,1]*(2*np.pi)**3

# Set number of mass and redshift bins for lightcone abundance matching
nMbins = 10000 
z1     = 0
z2     = 4.6
nzbin  = 46 
Mmin   = 5.e11
Mmax   = 1.e16

# Redshift functions
zinterp=np.linspace(0,5,10000)
def hubble(z,h,omegam):
    return h*100*np.sqrt(omegam*(1+z)**3+1-omegam)

def drdz(z,h,omegam):
    return 299792.458 / hubble(z,h,omegam)

rinterp  = np.cumsum(drdz(zinterp,h,omegam)*(zinterp[1]-zinterp[0]))
rinterp -= rinterp[0]
z_to_r   = sp.interpolate.interp1d(zinterp,rinterp)
r_to_z   = sp.interpolate.interp1d(rinterp,zinterp)

#Massbins
Mbins   = nMbins+1
Medge   = np.logspace(np.log10(Mmin),np.log10(Mmax),Mbins)
Mcent   = 10**((np.log10(Medge)[1:]+np.log10(Medge)[:-1])/2)

#Redshift Bins
zedge   = np.linspace(z1,z2,nzbin+1)
redge   = z_to_r(zedge)
zcent   = (zedge[1:]+zedge[:-1])/2
rcent   = (redge[1:]+redge[:-1])/2

Vol     = 4./3*np.pi*(redge[1:]**3-redge[:-1]**3) #assume full sky 
# print "volume of spherical shell ", Vol

# Halo Catalogue Functions
def get_pp_catalogue(filename,Mmin,Mmax):

    pkfile   = open(filename,"rb")
    Non      = np.fromfile(pkfile, dtype=np.int32, count=1)[0]
    RTHMAXin = np.fromfile(pkfile, dtype=np.float32, count=1)
    redshift = np.fromfile(pkfile, dtype=np.float32, count=1)

    outnum   = 11

    npkdata  = outnum*Non

    peakdata = np.fromfile(pkfile, dtype=np.float32, count=npkdata)
    peakdata = np.reshape(peakdata,(Non,outnum))

    RTH   = peakdata[:,6]
    Mpp   = 4.0/3*np.pi * RTH**3 * rho
    xpp   = peakdata[:,0]
    ypp   = peakdata[:,1]
    zpp   = peakdata[:,2]

    dm  = (Mpp > Mmin) & (Mpp < Mmax) #cut mass
    Mpp = Mpp[dm]
    xpp = xpp[dm]
    ypp = ypp[dm]
    zpp = zpp[dm]

    print "number of peak patch halos = ", len(Mpp)

    rpp   = np.sqrt(xpp**2+ypp**2+zpp**2)

    dm  = (Mpp > Mmin) & (Mpp < Mmax) #cut mass

    Mpp = Mpp[dm]
    rpp = rpp[dm]

    return Mpp, rpp


def cull_pp_catalogue(M_cat,r_cat,rmin,rmax):
    """
    keep only halos in distance range
    """

    dm  = (r_cat>rmin) & (r_cat<rmax) #cut distance

    M_cat = M_cat[dm]
    r_cat = r_cat[dm]

    if len(M_cat>0.): print M_cat.max()/1e15

    return M_cat, r_cat


def get_hmf(M,bins):
    """
    get N(M) histogram for hmf
    """
    hist, binedg = np.histogram(M,bins=bins)
    hist = hist+0.0

    return hist


# Load peakpatch hmf data
# This assumes you have a file containing an array of N(M,z) from the simulation already, with the same Mbins and Zbins 
if os.path.isfile(Nhalofile):
    ppdata = np.loadtxt(Nhalofile)
else:
    M_pp, r_pp   = get_pp_catalogue(halofile,Mmin,Mmax)

    N_pp  = np.zeros((nMbins,len(zcent)))
    for i in range(len(zcent)): #set up array of [NM_z0, NM_z1, NMz2, ... ]                                                                                              
        print "working on slice", zedge[i], " < z < ",zedge[i+1]
        print "r of zedges = ", z_to_r(zedge[i]), z_to_r(zedge[i+1])

        M_ppi, r_ppi = cull_pp_catalogue(M_pp,r_pp,redge[i],redge[i+1])
        
        N_pp[:,i]    = get_hmf(M_ppi,Medges)

    np.savetxt(Nhalofile,N_pp)



# Calculate functions needed for tinker hmf myself, instead of using hmfcalc
# When given the same input powerspectrum, and using the same Tinker parameters
# As hmfcalc (which seem to be slightly different than the Tinker 2008 paper)
# These functions match the results of hmfcalc to <1% 

def mass_to_radius(m, rho_mean):
    """                                                                                                                   
    Returns the radius of a sphere of uniform denstiy rho and mass m                                                      
    """
    r = (3 * m / 4. / np.pi / rho_mean)**(1./3.)

    return r

def windowfunction(x):
    """                                                                                                                   
    Computes the window function in Fourier space (top hat in real space).                                                
    """
    W = (3. / x**3) * (np.sin(x) - x * np.cos(x))

    return W

def M_to_sigma(k, pk, Medge, Mcent, omegaM, h):
    """                                                                                                                   
    Returns mass and sigma(mass)                          
    """

    rho_mean = 2.775e11 * omegaM * h**2

    # Normalization                                                                                                       
    Anorm    = 1./(2.*np.pi**2)

    # Now compute sigma(M)                                                                                                
    sigmaedge = np.zeros(Medge.shape[0])
    sigmacent = np.zeros(Mcent.shape[0])

    for i in range(sigmaedge.shape[0]):

        radius = mass_to_radius(Medge[i], rho_mean)

        x = k * radius
        y = pk * (k * windowfunction(x))**2

        sigmaedge[i] = np.sqrt(Anorm * sp.integrate.simps(y, k, even="avg"))

    for i in range(sigmacent.shape[0]):

        radius = mass_to_radius(Mcent[i], rho_mean)

        x = k * radius
        y = pk * (k * windowfunction(x))**2

        sigmacent [i] = np.sqrt(Anorm * sp.integrate.simps(y, k, even="avg"))
 
    return sigmaedge, sigmacent

def dlnsigmainv_dM(Medge, sigmaedge):
    lnsigmainv = np.log(1/sigmaedge)
    
    diff_sig = np.diff(lnsigmainv)
    diff_M   = np.diff(Medge)

    return diff_sig/diff_M

def growth_factor(omegaM, omegaL, z):
    """                                                                                                                   
    Returns growth factor using fitting formulae from Carrol, Press & Turner (1992)                                       
    """
    # returns growth factor using fitting formulae from Carrol, Press & Turner (1992)                                     

    w = -1.
    x = 1.+z
    x2 = x**2
    x3 = x**3
    x3w = x**(3*w)

    #calculate omegaM(z) and omegaL(z)                                                                                    
    denom = omegaM*x3 + (1-omegaM-omegaL)*x2 + omegaL*x3*x3w
    omega = omegaM*x3/denom
    lamb = omegaL*x3*x3w/denom

    #fitting formulae for w=-1 universe                                                                                   
    g = 2.5*omega/(omega**(4.0/7.0) - lamb + (1+(omega/2))*(1+(lamb/70)))
    g0 = 2.5*omegaM/(omegaM**(4.0/7.0) - omegaL + (1+(omegaM/2))*(1+(omegaL/70)))

    D = (g/x)/g0

    return D

def tinker_func(x, z):
    """                                                                                                                  
    Uses fitting coefficients from Table 2 of Tinker et al. (2008) for Delta = 200 and                                   
    redshift evolution equations 5 through 8.                                                                            
    """

    A_hmfcalc = 1.858659e-01
    a_hmfcalc = 1.466904
    b_hmfcalc = 2.571104
    c_hmfcalc = 1.193958
    
    A =  0.186
    a = 1.47
    b = 2.57
    c = 1.19
    
#     A =  A_hmfcalc
#     a =  a_hmfcalc
#     b =  b_hmfcalc
#     c =  c_hmfcalc

    amp = A * (1. + z)**-0.14
    a = a * (1. + z)**-0.06
    alpha = 10**(-1. * (0.75 / np.log10(200. / 75.))**1.2)

    b = b * (1. + z)**-alpha
    c = c

    f = amp * ((x / b)**(-a) + 1.) * np.exp(-c / x**2)

    return f

def tinker_hmf(Medge, Mcent, sigmaedge, sigmacent, redshift, dlnsigmainv, omegaM, omegaL, h):
    """                                                                                                                  
    Returns dn/dm for Tinker et al. (2008) mass function using Eqs. (3, 5-8) of their paper.                             
    """

    rho_mean = 2.775e11 * omegaM * h**2


    dndm = tinker_func(sigmacent, redshift) * rho_mean/Mcent * dlnsigmainv 
    
    ngtm = np.cumsum( (dndm * np.diff(Medge))[::-1])[::-1]
    
    return dndm, ngtm

def get_tinker(k, pk, Medge, Mcent, sigmaedge, sigmacent, redshift):
   
    D = growth_factor(omegam, omegal, redshift)

    sigmaedge_i = sigmaedge*D
    sigmacent_i = sigmacent*D

    fsigma = tinker_func(sigmacent_i, redshift)

    # Get dlnsigma^-1/dM
    dlnsigmainv = dlnsigmainv_dM(Medge, sigmaedge_i)

    # Get dndM, ngtM
    dndM, ngtM = tinker_hmf(Medge, Mcent, sigmaedge_i, sigmacent_i, redshift, dlnsigmainv, omegam, omegal, h)

    return dndM, ngtM


### Get Tinker in loop over redshifts

# Get sigma(M), scale later by growth factor
redshift = 0.
sigmaedge, sigmacent = M_to_sigma(k, pk, Medge, Mcent, omegam, h)

NgtM_tinker_total = np.zeros((len(Medge)-1,len(zcent)))
dNdM_tinker_total = np.zeros((len(Mcent),len(zcent)))

nslic_integral    = 10 # number of slices each shell is split into to get integrated mass function across shell

#set up array of [Mtinker, NgtM_z1, NgtMz2, ... ]
for i in range(len(zcent)): 
    if i%10==0: print "done zslic = ", i
    # integrate N(>M|z) over redshift range of shell by performing sum
    dr = (redge[i+1]-redge[i])/nslic_integral
    NgtM_tinker_i = np.zeros(len(Medge)-1)
    dNdM_tinker_i = np.zeros(len(Mcent))
                        
    for j in range(nslic_integral):
        zcent_j = r_to_z(redge[i]+(j+1./2)*dr)


        dNdM_i, ngtM_i = get_tinker(k, pk, Medge, Mcent, sigmaedge, sigmacent,  zcent_j) # get N(>M | z) [1/Mpc^3]

        NgtM_tinker_i += ngtM_i * 4./3*np.pi * ( (redge[i]+(j+1)*dr)**3 - (redge[i]+(j)*dr)**3) #in units of total number in fullsky sim
        dNdM_tinker_i += dNdM_i * 4./3*np.pi * ( (redge[i]+(j+1)*dr)**3 - (redge[i]+(j)*dr)**3) #in units of total number in fullsky sim

    NgtM_tinker_total[:,i] = NgtM_tinker_i
    dNdM_tinker_total[:,i] = dNdM_tinker_i

 

# Load Peakpatch N(M,z) data
Medges = Medge[:-1]
dNdM_pp_total  = (np.diff(ppdata,axis=0).T/np.diff(Medges)).T
NgtM_pp_total  = np.cumsum(ppdata[::-1,:],axis=0)[::-1,:]

# MAKE ABUNDANCE MATCH TABLE
M200_of_MTH_and_z = (np.zeros((nMbins,len(zcent))).T + Medge[:-1]).T #+Medges added so if above interpolation range just keep the same mass      


for i in range(len(zcent)): #set up array of [Mtinker, NgtM_z1, NgtMz2, ... ]                                                                     

    NgtM_pp     = NgtM_pp_total[:,i]
    NgtM_tinker = NgtM_tinker_total[:,i]
    
    Medgesi     = Medge[:-1]
    
    NgtM_tinker_interp = interp1d(Medgesi, NgtM_tinker)

    MTHofNtinkergtM  = interp1d(NgtM_tinker[::-1],Medgesi[::-1]) #as x needs to be monotonically increasing                                      

    dm          = NgtM_pp > 0
    NgtM_pp     = NgtM_pp[dm]
    Medgesi     = Medgesi[dm]

    MTHofN_pp     = MTHofNtinkergtM(NgtM_pp[::-1])[::-1]

    M200ofMTH      = interp1d(Medgesi,MTHofN_pp)

    M200_of_MTH_and_z[:len(Medgesi),i] = M200ofMTH(Medgesi)


# save abundance match table
np.savez("abundance_match_table_"+str(nMbins)+"Mbins_"+str(nzbin)+"zbins", M=Medge[:-1], z=zcent, M200_of_MTH_and_z=M200_of_MTH_and_z)


M200_of_MTH_and_z_interp =  sp.interpolate.RectBivariateSpline(Medge[:-1],zcent,M200_of_MTH_and_z,kx=1,ky=1)
Mtest = 1e12
print "mass before /1e12 %f, mass after /1e12 %f" % (Mtest/1e12, M200_of_MTH_and_z_interp(Mtest,0)/1e12)


#Plot mass completion
ncellcut = 10.
mcellcut = (7700./6144)**3 * rho_mean * ncellcut

mm = np.ones(len(zcent)) * mcellcut
zz = np.linspace(0,4.6,len(zcent))
mm_abmatch = M200_of_MTH_and_z_interp.ev(mm,zz)

np.savetxt("mass_completion_ncellcut"+str(int(ncellcut))+".txt", np.c_[zz,mm_abmatch], header="redshift\t Mass_min")
plt.semilogy(zz,mm_abmatch, "k", label = "mass completion "+str(ncellcut)+" cell")
plt.xlabel("redshift")
plt.ylabel("$M_{min}$(z) for 10cell halo abundance match")
plt.savefig("completeness.pdf")

