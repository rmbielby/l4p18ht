import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
def plaw(r,r0,gamma):
    return  (r/r0)**gamma


npoints = 10
s = 10.**(np.arange(10)*0.2-1.)
xis = plaw(s,2.5,-1.2)
xie = (1.+xis)/(np.arange(npoints)*24.+20.)**0.5

xispe = xis+xie*np.random.randn(npoints)
for i in np.arange(npoints):
    print s[i],xis[i],xispe[i],xie[i]


out = curve_fit(plaw,s,xispe,sigma=xie)
print 'r_0   = {0:5.2f}\pm{1:4.2f}'.format(out[0][0],out[1][0,0]**0.5)
print 'gamma = {0:5.2f}\pm{1:4.2f}'.format(out[0][1],out[1][1,1]**0.5)


f,ax = plt.subplots(1,1)
ax.loglog()
ax.errorbar(s,xispe,xie,linestyle=' ',marker='o')
ax.plot(s,plaw(s,*out[0]))
plt.show()
