import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate,special

# Plot the galaxy - Lyalpha forest cross correlation using the X-Shooter and UVES+HiRES datasets.
# Plot both hubble distance xi(s) and projected xi(sigma)
# Created 21/01/2014

from matplotlib import rcParams
rcParams['axes.labelsize'] = 16
rcParams['xtick.labelsize'] = 14
rcParams['ytick.labelsize'] = 14
rcParams['legend.fontsize'] = 14


modvdisp = 120.
modbeta  = 0.4
modr0    = 2.3
modgam   = 1.2

# Integration parameters
imin = -8000.
imax = -imin
h0   = 100.
def xiint(v,sigma,pi,r0,gamma,vdisp,betagal,betalya):
    r        = (sigma**2.+(pi-v/h0)**2.)**0.5
    costheta = (pi-v/h0)/r
    xir      = plaw(r,r0=r0,gamma=gamma)
    xi0      = (1.0+1.0*(betagal+betalya)/3.0+betagal*betalya/5.0)
    xi2      = (2.0*(betagal+betalya)/3.0 + 4.0*betagal*betalya/7.0)*(gamma/(gamma-3.0))
    xi4      = 8.0*betagal*betalya/35.0*(gamma*(2.0+gamma)/(3.0-gamma)/(5.0-gamma))
    xiprime  = xi0*p0(costheta) + xi2*p2(costheta) + xi4*p4(costheta)
    if vdisp > 0.:
        f = np.exp(-2.**0.5*np.abs(v)/vdisp)/(vdisp*2.**0.5)
    else:
        f = 1.0
    return xiprime*f*xir

def plaw(r,r0=0.25,gamma=1.1):
    return (r/r0)**(-np.abs(gamma))

def p0(mu):
    return 1.

def p2(mu):
    return 0.5*(3.*mu**2.-1.)

def p4(mu):
    return 0.125*(35.*mu**4. - 30.*mu**2.+3.)

def xismod(s,vdisp,betalya):
    ns      = len(s)
    ntheta  = 9
    dtheta  = np.pi/2./(ntheta)
    ximodel = np.zeros((ns,ntheta))
    xismean = np.zeros(ns)
    xismedi = np.zeros(ns)
    dsig    = 0.5*np.log10(s[1]/s[0])
    r0      = modr0
    gamma   = modgam
    betagal = 0.38
    for i in np.arange(ns):
        for j in np.arange(ntheta):
            sigma = s[i]*np.cos((dtheta*(j)))
            pi    = s[i]*np.sin((dtheta*(j)))
            args         = (sigma,pi,r0,gamma,vdisp,betagal,betalya)
            ximodel[i,j] = integrate.romberg(xiint,imin,imax,args=args,divmax=20)
            # print i,j,sigarr[i][j],piarr[i][j],np.log10(ximodel[i,j]),np.log10(((sigma[i]**2+sigma[j]**2)**0.5/modr0)**(-modgam))
        # xismean[i] = np.nanmean(ximodel[i,:])
        xismedi[i] = np.nanmedian(ximodel[i,:])
    return xismedi

s = 10.**np.arange(-1,1.6,0.1)
fig,axis = plt.subplots(1,1)
axis.loglog(s,plaw(s,modr0,modgam))
print 'Calculating xis model now!'
axis.loglog(s,xismod(s,modvdisp,modbeta))
plt.show()
