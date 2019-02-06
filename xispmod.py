import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import scipy.optimize as opt
from optparse import OptionParser
import argparse
import sys

# Calculate xi(sigma,pi) model.
# Example command line run:
# python xispmod.py -r 2.5 -g 1.2 -a 120 -b 0.4
# i.e. clustering length 2.5 Mpc, slope -1.2, velocity dispersion 120 km/s and
# Kaiser Boost beta 0.4.

from matplotlib import rcParams
rcParams['axes.labelsize'] = 16
rcParams['xtick.labelsize'] = 14
rcParams['ytick.labelsize'] = 14
rcParams['legend.fontsize'] = 12

h0   = 70.

parser = OptionParser()
parser.add_option('-r','--r0',dest='r0',
                  help='Power-law r0',type='float',
                  default=0.40)
parser.add_option('-g','--gamma',dest='gamma',
                  help='Power-law slope',type='float',
                  default=1.2)
parser.add_option('-a','--aparam',dest='aparam',
                  help='Initial guess for velocity dispersion',type='float',
                  default=120.)
parser.add_option('-b','--bparam',dest='bparam',
                  help='Initial guess for beta param',type='float',
                  default=0.4)

try:
    options,args = parser.parse_args(sys.argv[1:])
except getopt.GetoptError:
    print 'test.py -d <inputfile> -c <outputfile>'
    sys.exit(2)
print sys.argv[1:]
print parser.parse_args(sys.argv[1:])

# Integration parameters
imin = -8000.
imax = -imin

# Plot parameters
cmap = 'RdYlBu_r'
interp = "bessel"
vmin = np.log10(0.004)
vmax = np.log10(0.14)
smooth = 1.0
format = 'png'
pltmax = 100.
pltmin = 1.0

bin    = '0.167'
scale  = 'log'
ic     = 0.0048 #481

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

def plaw(r,r0=0.5,gamma=1.0):
    return (r/r0)**(-np.abs(gamma))

def p0(mu):
    return 1.

def p2(mu):
    return 0.5*(3.*mu**2.-1.)

def p4(mu):
    return 0.125*(35.*mu**4. - 30.*mu**2.+3.)


def xispmod((sigarr,piarr),vdisp,betalya):
    nsigma  = len(sigarr)

    print piarr
    ximodel = np.zeros((nsigma,nsigma))
    dsig = 0.5*np.log10(sigarr[1]/sigarr[0])
    r0      = options.r0
    gamma   = options.gamma
    betagal = 0.36
    for i in np.arange(nsigma):
        for j in np.arange(nsigma):
            print sigarr[i][j],piarr[i][j]
            args         = (sigarr[i][j],piarr[i][j],r0,gamma,vdisp,betagal,betalya)
            ximodel[i,j] = integrate.romberg(xiint,imin,imax,args=args,divmax=36)
    print 'Integrated model'
    return ximodel.ravel()

log = False
if log == True:
    nx, ny = (12, 12)
    x = np.linspace(-1, 2, nx)
    y = np.linspace(-1, 2, ny)
    xv, yv = np.meshgrid(10.**x, 10.**y)
    xplt, yplt = np.log10(xv),np.log10(yv)
else:
    nx, ny = (40, 40)
    x = np.linspace(0.5, 20, nx)
    y = np.linspace(0.5, 20, ny)
    xv, yv = np.meshgrid(x, y)
    xplt, yplt = xv,yv

print 'Computing model.'

veldisp = options.aparam
kaiser  = options.bparam

popt = (veldisp,kaiser)
ximodel = xispmod((xv,yv), *popt)
print 'Done model'

nsigma = len(xv)
model_image = np.log10(ximodel.reshape(nsigma, nsigma))
fig,axis = plt.subplots(1,1)
axis.contour(xplt, yplt, model_image, 8)
plt.show()
